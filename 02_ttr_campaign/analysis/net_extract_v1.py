#!/usr/bin/env python3
"""네트워크 전수 데이터 추출 (2026-10-04) — ★읽기 전용 · 새 계산 없음 · 새 판정 없음.

동결 원본(JSON·TSV·md.log·mdp)의 수치를 그대로 TSV 로 옮긴다.
  · u_nk 추출 · MBAR · 교환 로그 재집계 · 왕복 재구성 없음.
  · 원본에 없는 값은 '미계산' / '미보유', 못 읽으면 '측정 불가'.
  · ★표시용 산술은 C2 의 레그 합(Σ복합체 − Σ용매) 하나뿐이다(지시서 C 가 요구). 열 이름에 '산술' 표기.
  · 최솟값·그 λ 쌍은 원본 목록에서 원소를 고른 것이다(평균·표준편차는 새로 내지 않는다).
  · 판정 문구는 원본 필드를 그대로 옮긴 것이고, 이 스크립트가 새로 붙인 판정은 없다.

사용:  python3 net_extract_v1.py <출력 디렉터리>      (없는 디렉터리여야 한다 · 덮어쓰기 차단)
       python3 net_extract_v1.py --selftest           (대조 함수가 변조 입력에서 실패하는지 = 음성 대조)
       NET_EXTRACT_RP=<reports 변조 사본> python3 net_extract_v1.py <임시 출력>   (값 하나를 바꾼 사본에서 종료코드 1 이어야 한다)
"""
import csv, hashlib, json, os, re, subprocess, sys
from pathlib import Path

T = Path("/home/nudge/Project/CADD/ttr")
Q = T / "qsar"
B = Q / "data/fep_m1"
RP = Path(os.environ.get("NET_EXTRACT_RP", Q / "reports"))   # ★음성 대조용: 변조 사본 디렉터리를 가리키게 할 수 있다(기본 = 원본)
FR = Q / "docs/stage1_out/final_report_20260921"
P1 = Q / "docs/stage1_out/pass1"
PREREG = Q / "models/prereg.json"
NC, NH, NM = "미계산", "미보유", "측정 불가"
LEGS = [(s, l) for s in ("cplx", "solv") for l in (1, 2, 3)]
LK = [f"{s}L{l}" for s, l in LEGS]
N16 = [f"N{i:02d}" for i in range(1, 17)]


def rel(p):
    try: return str(Path(p).relative_to(T))
    except ValueError: return str(p)


def J(p): return json.load(open(p, encoding="utf-8"))


def head(p, n=400_000):
    with open(p, "rb") as f: return f.read(n).decode("utf-8", "replace")


def tail(p, n=400_000):
    with open(p, "rb") as f:
        f.seek(0, 2); sz = f.tell(); f.seek(max(0, sz - n)); return f.read().decode("utf-8", "replace")


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""): h.update(c)
    d = h.hexdigest()
    assert len(d) == 32, "md5 길이"
    return d


# ───────────────────────── 대조(검사) — 실패를 모아 끝에 비-0 종료 ─────────────────────────
CHECKS = []


def chk(name, a, b, tol=0.0):
    """a 와 b 가 같은가(수치는 tol 이내 · 그 밖은 ==). 결과를 CHECKS 에 쌓고 bool 을 돌려준다."""
    try:
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            ok = abs(a - b) <= tol
        elif isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
            ok = len(a) == len(b) and all(
                (abs(x - y) <= tol) if isinstance(x, (int, float)) and isinstance(y, (int, float)) else x == y
                for x, y in zip(a, b))
        else:
            ok = a == b
    except Exception:
        ok = False
    CHECKS.append((name, "일치" if ok else "★불일치", repr(a)[:120], repr(b)[:120], tol))
    return ok


def selftest():
    """음성 대조 — 변조 입력에서 chk 가 실패를 내야 한다. 하나라도 통과하면 종료코드 1."""
    bad = 0
    cases = [("수치 1e-3 변조", 1.5223773659636777, 1.5233773659636777, 1e-9),
             ("목록 원소 변조", [0.1, 0.2, 0.3], [0.1, 0.2, 0.31], 1e-9),
             ("목록 길이 변조", [1, 2, 3], [1, 2], 0),
             ("문자열 변조", "통과", "★미달", 0),
             ("None 대 수치", None, 0.0, 0)]
    for nm, a, b, tol in cases:
        if chk("selftest:" + nm, a, b, tol): bad += 1; print("★음성 대조 실패(통과해 버림):", nm)
    pos = chk("selftest:양성 대조", 0.25, 0.25, 0) and chk("selftest:양성 목록", [1, 2.0], [1, 2.0], 0)
    if not pos: bad += 1; print("★양성 대조 실패")
    try:
        assert len("abc") == 32
        bad += 1; print("★해시 길이 음성 대조 실패")
    except AssertionError:
        pass
    print(f"selftest: 음성 대조 {len(cases)}건 + 해시 길이 1건 + 양성 2건 → 이상 {bad}")
    sys.exit(1 if bad else 0)


# ───────────────────────── md.log 읽기 ─────────────────────────
def log_params(p):
    """md.log 머리의 입력 파라미터 덤프에서 λ 벡터·nsteps·dt·ld-seed 등을 읽는다(★첫 기동 것)."""
    if not Path(p).exists(): return None
    h = head(p)
    d = {}
    for k, pat in (("gmx", r"GROMACS version:\s+(\S+)"), ("nsteps", r"^\s+nsteps\s+=\s+(-?\d+)"),
                   ("dt", r"^\s+dt\s+=\s+(\S+)"), ("ld_seed", r"^\s+ld-seed\s+=\s+(-?\d+)"),
                   ("n_lambdas", r"^\s+n-lambdas\s+=\s+(\d+)"), ("init_state", r"^\s+init-lambda-state\s+=\s+(-?\d+)"),
                   ("nstdhdl", r"^\s+nstdhdl\s+=\s+(\d+)")):
        m = re.search(pat, h, re.M); d[k] = m.group(1) if m else NM
    m = re.search(r"Command line:\s*\n\s*(.+)", h)
    cmd = m.group(1) if m else ""
    m2 = re.search(r"-replex\s+(\d+)", cmd); d["replex"] = m2.group(1) if m2 else NM
    m2 = re.search(r"-nex\s+(\d+)", cmd); d["nex"] = m2.group(1) if m2 else NM
    i = h.find("all-lambdas:")
    for c in ("fep", "coul", "vdw", "bonded"):
        d[c] = NM
    if i >= 0:
        blk = h[i:i + 6000]
        for c in ("fep", "coul", "vdw", "bonded"):
            m = re.search(rf"^\s*{c}-lambdas =\s+(.*)$", blk, re.M)
            if m: d[c] = " ".join(m.group(1).split())
    return d


def log_tail(p):
    """md.log 끝의 종료 통계와 GROMACS 가 스스로 찍은 교환 요약 블록(마지막 것)."""
    if not Path(p).exists(): return None
    t = tail(p)
    d = {}
    m = re.findall(r"Statistics over (\d+) steps", t); d["stat_steps"] = m[-1] if m else NM
    d["n_summary_blocks_in_tail"] = len(re.findall(r"Repl\s+average probabilities:", t))
    m = list(re.finditer(r"Repl\s+(\d+) attempts, (\d+) odd, (\d+) even", t))
    d["attempts"] = m[-1].group(1) if m else NM
    d["odd"] = m[-1].group(2) if m else NM; d["even"] = m[-1].group(3) if m else NM
    m = list(re.finditer(r"Repl\s+average probabilities:\s*\nRepl\s+([\d ]+)\nRepl\s+(.+)", t))
    d["probs"] = m[-1].group(2).split() if m else None
    m = list(re.finditer(r"Repl\s+number of exchanges:\s*\nRepl\s+([\d ]+)\nRepl\s+(.+)", t))
    d["nexch"] = m[-1].group(2).split() if m else None
    return d


def varying(d):
    """레그에서 실제로 변하는 λ 성분과 그 값 목록(문자열 그대로)."""
    for c in ("coul", "vdw", "bonded", "fep"):
        v = d.get(c, NM)
        if v != NM and len(set(v.split())) > 1: return c, v.split()
    return NM, []


def mdp_vals(p):
    out = {"gen_vel": "키 없음", "gen_seed": "키 없음", "ld_seed": "키 없음", "continuation": "키 없음"}
    if not Path(p).exists(): return None
    for ln in open(p, encoding="utf-8", errors="replace"):
        ln = ln.split(";")[0]
        if "=" not in ln: continue
        k, v = ln.split("=", 1); k = k.strip().replace("-", "_"); v = v.strip()
        if k in out: out[k] = v
    return out


def mdp_lam(p):
    d = {}
    if not Path(p).exists(): return d
    for ln in open(p, encoding="utf-8", errors="replace"):
        m = re.match(r"\s*(fep|coul|vdw|bonded)-lambdas\s*=\s*(.*)$", ln)
        if m: d[m.group(1)] = " ".join(m.group(2).split())
        m = re.match(r"\s*nsteps\s*=\s*(\d+)", ln)
        if m: d["nsteps"] = m.group(1)
    return d


def W(out, name, header, rows):
    p = out / name
    with open(p, "x", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(header); w.writerows(rows)
    return p


def f(x, nd=None):
    if x is None: return NC
    if isinstance(x, float): return repr(x) if nd is None else f"{x:.{nd}f}"
    return str(x)


def argmins(vals):
    """최솟값과 같은 값을 가진 모든 색인(동률을 숨기지 않는다)."""
    mn = min(vals); return mn, [i for i, v in enumerate(vals) if v == mn]


def pairs_str(idx, lam, n):
    if len(lam) != n + 1: return NM
    return " · ".join(f"{lam[i]}–{lam[i + 1]}" for i in idx)


def edge_id_re(e): return re.compile(rf"(?<![A-Za-z0-9]){e}(?!\d)")


def main(out):
    out.mkdir(parents=True, exist_ok=False)
    snap = subprocess.run(["date", "-Is"], capture_output=True, text=True).stdout.strip()
    boot = subprocess.run(["uptime", "-s"], capture_output=True, text=True).stdout.strip()

    # ── 엣지 목록 ──
    EDGE = {}
    for fn in ("edges_v7.txt", "edges_v12_add.txt"):
        for ln in open(B / fn):
            q = ln.split()
            if len(q) == 3: EDGE[q[0]] = (q[1], q[2], f"qsar/data/fep_m1/{fn}")
    ALL = sorted(EDGE)

    # ── 원본 적재 ──
    final_edges = J(P1 / "final_edges.json")
    cur_bar = J(FR / "current_bar.json")
    std_an = J(FR / "std_analysis.json")
    replex = J(FR / "replex.json")["per_leg"]
    rtrip = J(FR / "roundtrip.json")
    resid = J(FR / "residuals_v2.json")
    preproc = J(FR / "preproc_sensitivity.json")
    cyc_loc = J(FR / "cycle_localization.json")
    signv = J(FR / "sign_verification.json")
    headline = J(FR / "headline.json")
    c6raw = J(RP / "c6_raw_v13.json")
    c6loo = J(RP / "c6_loo_v13.json")["rows"]
    rescore = J(RP / "rescore_v13.json")
    assign = J(P1 / "assignment.json")
    lf_res = J(P1 / "local_fill_result.json"); lf_der = J(P1 / "local_fill_derived.json")
    n01n5 = J(FR / "n01_n5_standard.json")
    ext_edges = J(P1 / "ext_edges_20260921.json"); ext_state = J(P1 / "ext_state.json")
    S1 = J(RP / "S1_ddg_all.json")
    uq1 = J(RP / "uq1_prod_v14_N17.json")
    s3m = J(RP / "s3_mbar_prod_v14_N17.json")
    s3x = J(RP / "s3_exchange_prod_v14_N17.json")
    hr1t = {(r["rep"], r["leg"]): r for r in J(RP / "hr1_table_prod_v14_N17.json")}
    hr1w = {(r["rep"], r["leg"]): r for r in J(RP / "hr1_wide_prod_v14_N17.json")}
    s1r5 = {r["leg"]: r for r in J(RP / "S1_r5_diag.json")}
    prereg = J(PREREG)
    v31 = next(e for e in prereg["changelog"] if e.get("version") == 31)
    v31b, v31c = v31["(b)★최종 값"], v31["(c)★판정 (등재 문구 그대로)"]

    # ═════════ A2 — λ 실값 (md.log 의 all-lambdas · 회차별) ═════════
    LOG = {}        # (edge, rep, leg키) → dict(log 파라미터)
    a2 = []
    for e in N16:
        for s, l in LEGS:
            for rep in range(1, 6):
                tag = f"{e}_{s}_L{l}_r{rep}"; d0 = B / "pass1" / tag
                lg = d0 / "w00/prod.log"
                if not lg.exists(): continue
                d = log_params(lg); t = log_tail(lg); d.update(t)
                d["ndir"] = len([x for x in d0.glob("w[0-9][0-9]") if x.is_dir()])
                d["ok"] = (d0 / "OK").exists(); d["log"] = rel(lg)
                LOG[(e, rep, f"{s}L{l}")] = d
    for rep in range(1, 6):
        for s, l in LEGS:
            base = B / f"prod_v14_r{rep}/N17/{s}/r{rep}"
            lg = base / f"L{l}w00/prod.log"
            if not lg.exists(): continue
            d = log_params(lg); t = log_tail(lg); d.update(t)
            d["ndir"] = len([x for x in base.glob(f"L{l}w[0-9][0-9]") if x.is_dir()])
            d["ok"] = (B / f"prod_v14_r{rep}/.rep_done").exists(); d["log"] = rel(lg)
            ml = mdp_lam(base / f"L{l}w00/prod.mdp")
            d["mdp"] = ml; d["mdp_path"] = rel(base / f"L{l}w00/prod.mdp")
            LOG[("N17", rep, f"{s}L{l}")] = d
    for (e, rep, k), d in sorted(LOG.items()):
        comp, vals = varying(d)
        if k.endswith("L2"):
            chk(f"A:{e} r{rep} {k} L2 는 vdw=bonded=fep 가 변하고 coul 은 전부 0", (d["vdw"] == d["bonded"] == d["fep"], comp, set(d["coul"].split())), (True, "vdw", {"0"}))
        else:
            chk(f"A:{e} r{rep} {k} L1·L3 은 coul 만 변함", (comp, set(d["vdw"].split()), set(d["fep"].split())), ("coul", {"0"}, {"0"}))
        chk(f"A:{e} r{rep} {k} n-lambdas = 창 디렉터리 수", int(d["n_lambdas"]), d["ndir"])
        chk(f"A:{e} r{rep} {k} 종료 통계 steps = nsteps+1", int(d["stat_steps"]), int(d["nsteps"]) + 1)
        try: ns = f"{int(d['nsteps']) * float(d['dt']) / 1000:g}"
        except Exception: ns = NM
        ml = d.get("mdp", {})
        a2.append([e, f"r{rep}", k, d["gmx"], d["n_lambdas"], d["ndir"], comp, " ".join(vals) if vals else NM,
                   d["coul"], d["vdw"], d["bonded"], d["fep"], d["nsteps"], d["dt"], ns, d["stat_steps"],
                   d["replex"], d["nex"], d["ld_seed"],
                   ml.get("coul", "—") if e == "N17" else "—", ml.get("vdw", "—") if e == "N17" else "—", d["log"]])
    W(out, "A2_lambda_values.tsv",
      ["edge", "rep", "leg", "gromacs(log)", "n-lambdas(log)", "창 디렉터리 수", "변하는 성분", "λ 실값(변하는 성분)",
       "coul-lambdas(log)", "vdw-lambdas(log)", "bonded-lambdas(log)", "fep-lambdas(log)", "nsteps(log·tpr)", "dt",
       "nsteps×dt (ns)", "종료 통계 steps(log 끝)", "-replex(log 명령줄·첫 기동)", "-nex", "ld-seed(log)",
       "coul-lambdas(prod.mdp·N17만)", "vdw-lambdas(prod.mdp·N17만)", "출처 log"], a2)

    # λ 요약: 엣지×레그 — 회차 간 λ 문자열이 같은가
    a2s = []
    for e in ALL:
        for k in LK:
            reps = sorted(r for (ee, r, kk) in LOG if ee == e and kk == k)
            strs = {(LOG[(e, r, k)]["coul"], LOG[(e, r, k)]["vdw"], LOG[(e, r, k)]["bonded"], LOG[(e, r, k)]["fep"]) for r in reps}
            comp, vals = varying(LOG[(e, reps[0], k)]) if reps else (NM, [])
            nl = sorted({LOG[(e, r, k)]["n_lambdas"] for r in reps}); nd = sorted({str(LOG[(e, r, k)]["ndir"]) for r in reps})
            ns = sorted({LOG[(e, r, k)]["nsteps"] for r in reps})
            a2s.append([e, k, " ".join(f"r{r}" for r in reps), len(strs), "같다" if len(strs) == 1 else "★다르다",
                        "/".join(nl), "/".join(nd), comp, " ".join(vals) if vals else NM, "/".join(ns)])
    W(out, "A2s_lambda_by_leg.tsv",
      ["edge", "leg", "log 가 있는 회차", "서로 다른 λ 문자열 묶음 수", "회차 간 λ(문자열 일치)", "n-lambdas", "창 디렉터리 수",
       "변하는 성분", "λ 실값(첫 회차)", "nsteps(log)"], a2s)

    # ═════════ A3 — prereg 에서 엣지 이름이 등장하는 항목 ═════════
    a3 = []; net7 = prereg["fep_network_v7"]["엣지"]
    top_keys = [k for k in prereg if k != "changelog"]
    for e in ALL:
        A_, B_, _ = EDGE[e]; rx = edge_id_re(e)
        idx = [str(i) for i, x in enumerate(net7) if {x["A"], x["B"]} == {A_, B_}]
        a3.append([e, "fep_network_v7.엣지[색인]", ",".join(idx) if idx else "없음(16엣지 목록에 이 쌍 없음)", "쌍 일치"])
        for k in top_keys:
            if rx.search(json.dumps(prereg[k], ensure_ascii=False)): a3.append([e, "최상위 키", k, "엣지 이름 등장"])
        for i, c in enumerate(prereg["changelog"]):
            if rx.search(json.dumps(c, ensure_ascii=False)):
                a3.append([e, "changelog", f"[{i}] v{c.get('version')} · {c.get('at') or c.get('date')} · {c.get('키', '')}", "엣지 이름 등장"])
    for i, c in enumerate(prereg["changelog"]):
        if "prod_v14" in json.dumps(c, ensure_ascii=False):
            a3.append(["N17", "changelog", f"[{i}] v{c.get('version')} · {c.get('at') or c.get('date')} · {c.get('키', '')}", "'prod_v14' 문자열 등장"])
    W(out, "A3_prereg_mentions.tsv", ["edge", "위치", "항목", "일치 방식"], a3)
    COMMON = [k for k in top_keys if k.startswith("fep_")]

    # ═════════ A1 — 엣지 현황 ═════════
    lab = {}
    for ln in open(P1 / "final_report_20260921.txt", encoding="utf-8"):
        m = re.match(r"(N\d\d)\s+\d+\s+\S+\s+\S+.*?\s(\S*(?:완료|연장후보)\S*)\s*$", ln)
        if m: lab[m.group(1)] = m.group(2)
    a1v = {}
    for ln in open(P1 / "EDGE_TABLE_20260921.md", encoding="utf-8"):
        q = [x.strip() for x in ln.split("|")]
        if len(q) > 10 and re.fullmatch(r"N\d\d", q[1]): a1v[q[1]] = (q[9], q[10])
    flags = {}
    for ln in open(FR / "edge_flags.md", encoding="utf-8"):
        q = [x.strip() for x in ln.split("|")]
        if len(q) > 5 and re.fullmatch(r"★?\*{0,2}N\d\d\*{0,2}", q[1]): flags[re.search(r"N\d\d", q[1]).group(0)] = q[5]
    extn = {}
    for k in ext_state.get("ok", {}):
        extn[k[:3]] = extn.get(k[:3], 0) + 1
    a1 = []
    for e in ALL:
        A_, B_, src = EDGE[e]
        Ks = {k: "/".join(sorted({LOG[(ee, r, kk)]["n_lambdas"] for (ee, r, kk) in LOG if ee == e and kk == k})) for k in LK}
        nsset = sorted({a[14] for a in a2 if a[0] == e})
        ments = sorted({r[2].split(" · ")[0] for r in a3 if r[0] == e and r[1] == "changelog"})
        tk = [r[2] for r in a3 if r[0] == e and r[1] == "최상위 키"]
        if e == "N17":
            reps = [f"r{r}" for r in range(1, 6) if f"r{r}" in S1]
            done = [f"r{r}" for r in range(1, 6) if (B / f"prod_v14_r{r}/.rep_done").exists()]
            a1.append([e, f"{A_} → {B_}", src, "완료",
                       f"S1_ddg_all.json 에 ΔΔG {len(reps)}개 · .rep_done 마커 {','.join(done) or '없음'} · prereg v31 (a) 동결 선언",
                       "— (동결 채점표 final_edges.json 밖)", len(reps), " ".join(reps),
                       f"{Ks['cplxL1']}/{Ks['cplxL2']}/{Ks['cplxL3']}", f"{Ks['solvL1']}/{Ks['solvL2']}/{Ks['solvL3']}",
                       "/".join(nsset), "prod_v14", "qsar/data/fep_m1/prod_v14_r{1..5}/N17 (r2·r3 = remote_pull/r{2,3} 심볼릭 링크)",
                       "—", "; ".join(tk) or "없음", "; ".join(ments)])
        else:
            fe = final_edges[e]; tags = [f"{e}_{s}_L{l}_r{r}" for s, l in LEGS for r in range(1, 6)]
            have = [t for t in tags if (B / "pass1" / t / "w00/prod.log").exists()]
            ok = [t for t in have if (B / "pass1" / t / "OK").exists()]
            a1.append([e, f"{A_} → {B_}", src, "완료",
                       f"final_edges.json 에 n={fe['n']} 로 채점됨 · pass1 태그 {len(have)}/30 보유 · OK 마커 {len(ok)}/{len(have)}",
                       f"{lab.get(e, NH)} · A1 판정 '{a1v.get(e, (NH, NH))[0]}' · 문제란 '{a1v.get(e, (NH, NH))[1]}' · 수렴 플래그 '{flags.get(e, NH)}'",
                       fe["n"], " ".join(f"r{r}" for r in c6loo[e]["reps"]),
                       f"{Ks['cplxL1']}/{Ks['cplxL2']}/{Ks['cplxL3']}", f"{Ks['solvL1']}/{Ks['solvL2']}/{Ks['solvL3']}",
                       "/".join(nsset), "pass1 (원 트리 prod_v9c_r1~r5 · assignment.json 의 dir)",
                       "qsar/data/fep_m1/pass1/<TAG>/ (xvg·log·cpt) · pass1_full/<TAG>/ (tpr·gro·xtc·edr)",
                       f"{extn.get(e, 0)} 태그" if extn.get(e) else "0", "; ".join(tk) or "없음", "; ".join(ments) or "없음"])
    W(out, "A1_edges.tsv",
      ["edge", "리간드 쌍(A → B)", "엣지 목록 출처", "상태", "상태 근거(원본 사실)", "원본에 적힌 표기(옮김)", "채점 회차 수 n", "채점 회차",
       "n-lambdas 복합체 L1/L2/L3 (log)", "n-lambdas 용매 L1/L2/L3 (log)", "창당 ns (nsteps×dt · log)", "prod 버전", "산출 트리",
       "pass_ext 2.0 ns 연장 보유(ext_state.json ok)", "prereg 최상위 키(엣지 이름 등장)", "prereg changelog(엣지 이름 등장)"], a1)

    # ═════════ B1 — 회차별 ΔΔG ═════════
    b1 = []
    for e in N16:
        dd = c6raw[e]["dd"]; pp = preproc[e]; reps = c6loo[e]["reps"]
        chk(f"B:{e} final_edges.n = c6_raw.dd 개수", final_edges[e]["n"], len(dd))
        chk(f"B:{e} preproc.per_cur = c6_raw.dd", [float(x) for x in pp["per_cur"]], [dd[str(r)] for r in reps], 1e-9)
        for i, r in enumerate(reps):
            b1.append([e, f"r{r}", f(dd[str(r)]), "gmx bar -b 0 (동결 채점)", NC, "qsar/reports/c6_raw_v13.json", f"{e}.dd.{r}",
                       f(pp["per_dec"][i]) if len(pp["per_dec"]) == len(reps) else NC,
                       "preproc_sensitivity.json per_dec[회차 오름차순]"])
    b1.append(["N01", "r1 (로컬 보충)", f(lf_res["n2"][0]), "gmx bar (로컬 보충분 · 채택 n=5 의 r1)", NC,
               "qsar/docs/stage1_out/pass1/local_fill_result.json", "n2[0]", NC, "—"])
    b1.append(["N01", "r2 (로컬 보충)", f(lf_res["n2"][1]), "gmx bar (로컬 보충분 · 채택 n=5 의 r2)", NC,
               "qsar/docs/stage1_out/pass1/local_fill_result.json", "n2[1]", NC, "—"])
    if ext_edges.get("N04"):
        for i, v in enumerate(ext_edges["N04"]["per"]):
            b1.append(["N04", f"per[{i}] (2.0 ns 연장판 · 회차 표지 원본에 없음)", f(v), "ext_bar (2.0 ns)", NC,
                       "qsar/docs/stage1_out/pass1/ext_edges_20260921.json", f"N04.per[{i}]", NC, "—"])
    for r in range(1, 6):
        k = f"r{r}"
        chk(f"B:N17 {k} v31 등재 ΔΔG = S1_ddg_all", v31b["ΔΔG (S1_ddg_all.json 그대로)"][k], S1[k], 0)
        b1.append(["N17", k, f(S1[k]), "MBAR (alchemlyb · 8 ns 전량)", f(v31b["MBAR σ (보조 — 과소추정 가능)"][k]),
                   "qsar/reports/S1_ddg_all.json · MBAR σ = prereg v31 (b) (S1_tableD.txt 인쇄값)", k, NC, "—"])
    W(out, "B1_ddg_per_rep.tsv",
      ["edge", "rep", "ΔΔG (kcal/mol)", "추정기", "ΔΔG 의 MBAR σ (kcal/mol)", "출처", "키",
       "ΔΔG 표준 전처리판(MBAR·평형검출·상관제거) (kcal/mol)", "그 출처"], b1)

    # ═════════ B2 — 회차 통계 (등재·저장된 값만) ═════════
    b2 = []
    for e in N16:
        fe = final_edges[e]; rs = resid[e]
        chk(f"B:{e} residuals_v2.sd = final_edges.sd", rs["sd"], fe["sd"], 1e-12)
        b2.append([e, "동결 채점(패스1 0.5 ns)", fe["n"], f(fe["dd"]), f(fe["sd"]), f(rs["se"]), NC, NC, NC,
                   "final_edges.json {n,dd,sd} · SEM = final_report_20260921/residuals_v2.json 'se'"])
        pp = preproc[e]
        b2.append([e, "표준 전처리(MBAR·평형검출·상관제거)", pp["n"], f(pp["dec"]), f(pp["sd_dec"]), NC, NC, NC, NC,
                   "final_report_20260921/preproc_sensitivity.json {n,dec,sd_dec}"])
    b2.append(["N01", "채택 n=5 (동결 채점 방식 · 로컬 보충 포함)", 5, f(lf_der["dd5"]), f(lf_der["sd5"]), NC, NC, NC, NC,
               "pass1/local_fill_derived.json {dd5,sd5}"])
    b2.append(["N01", "채택 n=5 · MBAR raw", n01n5["raw"]["n"], f(n01n5["raw"]["dd"]), f(n01n5["raw"]["sd"]), NC, NC, NC, NC,
               "final_report_20260921/n01_n5_standard.json raw"])
    b2.append(["N01", "채택 n=5 · MBAR 상관제거", n01n5["dec"]["n"], f(n01n5["dec"]["dd"]), f(n01n5["dec"]["sd"]), NC, NC, NC, NC,
               "final_report_20260921/n01_n5_standard.json dec"])
    if ext_edges.get("N04"):
        x = ext_edges["N04"]
        b2.append(["N04", "2.0 ns 연장판", x["n"], f(x["dd"]), f(x["sd"]), NC, NC, NC, NC, "pass1/ext_edges_20260921.json N04"])
    for nm, key in (("주 집합 {1,4,5} (등재 σ 정의)", "주 집합 {1,4,5}"), ("보조 집합 {1,2,3,4,5} (등재 정의 밖)", "보조 집합 {1,2,3,4,5}")):
        x = v31b[key]
        b2.append(["N17", nm, x["n"], f(x["평균"]), f(x["s"]), f(x["SEM"]), f(x["t"]), f(x["95% 구간"][0]), f(x["95% 구간"][1]),
                   "prereg.json changelog v31 (b) — S1_tableD.txt 인쇄값(4자리)"])
    W(out, "B2_ddg_stats.tsv",
      ["edge", "집합/채점", "n", "평균 ΔΔG (kcal/mol)", "표준편차 s (ddof=1)", "SEM", "t", "95% 하한", "95% 상한", "출처"], b2)

    # ═════════ C1 — 회차 × 레그 ΔG ═════════
    c1 = []
    for e in N16:
        sa = {k: dict(std_an[e][k]) for k in LK}
        for k in LK:
            s, l = k[:4], k[5]
            for rep in range(1, 6):
                tag = f"{e}_{s}_L{l}_r{rep}"
                if tag not in cur_bar and tag not in sa[k]: continue
                cb = cur_bar.get(tag); x = sa[k].get(tag, {})
                if cb is not None and str(rep) in c6raw[e]["legs"].get(k, {}):
                    chk(f"C:{tag} c6_raw.legs = current_bar", c6raw[e]["legs"][k][str(rep)], cb, 0)
                c1.append([e, f"r{rep}", k, f(cb), NH + "(gmx bar 오차 저장 안 됨)", f(x.get("MBAR_raw")), f(x.get("MBAR_raw_err")),
                           f(x.get("MBAR_dec")), f(x.get("MBAR_dec_err")), f(x.get("nwin")), f(x.get("nframe_raw")), f(x.get("nframe_dec")),
                           "current_bar.json[TAG] · std_analysis.json[edge][leg][TAG]"])
    for r in range(1, 6):
        for k in LK:
            x = uq1[f"r{r}"][k]
            chk(f"C:N17 r{r} {k} s3_mbar.dG = uq1.dG", s3m[f"r{r}"][k]["dG"], x["dG"], 1e-9)
            c1.append(["N17", f"r{r}", k, "—", "—", f(x["dG"]), f(x["mbar_sig"]), "—", "—", f(x["nwin"]),
                       f(sum(s3m[f"r{r}"][k]["nrow"])), "—", f"uq1_prod_v14_N17.json r{r}.{k}.{{dG,mbar_sig}}"])
    W(out, "C1_leg_dG.tsv",
      ["edge", "rep", "leg", "ΔG gmx bar (kcal/mol·동결 채점)", "gmx bar σ", "ΔG MBAR 전량 (kcal/mol)", "MBAR σ (kcal/mol)",
       "ΔG MBAR 상관제거 (kcal/mol)", "MBAR σ 상관제거", "창 수(파일 수)", "프레임 수(raw)", "프레임 수(상관제거 뒤)", "출처·키"], c1)

    # ═════════ C2 — 회차별 레그 합(표시용 산술)과 저장된 ΔΔG ═════════
    c2 = []
    for e in N16:
        for rep in range(1, 6):
            v = {k: cur_bar.get(f"{e}_{k[:4]}_L{k[5]}_r{rep}") for k in LK}
            have = [k for k in LK if v[k] is not None]
            dd = c6raw[e]["dd"].get(str(rep))
            if not have: continue
            if len(have) == 6:
                sc = v["cplxL1"] + v["cplxL2"] + v["cplxL3"]; ss = v["solvL1"] + v["solvL2"] + v["solvL3"]
                if dd is not None: chk(f"C:{e} r{rep} 레그 합 = 저장 ΔΔG", sc - ss, dd, 1e-9)
                c2.append([e, f"r{rep}", "gmx bar", f(sc), f(ss), f(sc - ss), f(dd) if dd is not None else NH, "6/6", "current_bar.json · c6_raw_v13.json dd"])
            else:
                c2.append([e, f"r{rep}", "gmx bar", NC, NC, NC, f(dd) if dd is not None else NH,
                           f"{len(have)}/6 (없는 레그: {' '.join(k for k in LK if v[k] is None)})", "current_bar.json · c6_raw_v13.json dd"])
    for r in range(1, 6):
        v = {k: uq1[f"r{r}"][k]["dG"] for k in LK}
        sc = v["cplxL1"] + v["cplxL2"] + v["cplxL3"]; ss = v["solvL1"] + v["solvL2"] + v["solvL3"]
        chk(f"C:N17 r{r} 레그 합 = S1_ddg_all", sc - ss, S1[f"r{r}"], 1e-9)
        c2.append(["N17", f"r{r}", "MBAR", f(sc), f(ss), f(sc - ss), f(S1[f"r{r}"]), "6/6", "uq1_prod_v14_N17.json dG · S1_ddg_all.json"])
    W(out, "C2_leg_sum_vs_ddg.tsv",
      ["edge", "rep", "추정기", "Σ복합체 3레그 (산술·kcal/mol)", "Σ용매 3레그 (산술)", "Σ복합체−Σ용매 (산술)", "저장된 그 회차 ΔΔG", "보유 레그", "출처"], c2)

    # ═════════ D1 — 겹침(1차 비대각) ═════════
    d1 = []
    for e in N16:
        for k in LK:
            for tag, x in std_an[e][k]:
                rep = int(tag.split("_r")[-1]); ov = x.get("overlap_adj")
                lg = LOG.get((e, rep, k))
                if not ov: d1.append([e, f"r{rep}", k, NM, NM, NM, NM, NM, "std_analysis.json overlap_adj 없음"]); continue
                mn, ii = argmins(ov); chk(f"D:{tag} overlap_min = min(overlap_adj)", x["overlap_min"], mn, 0)
                comp, vals = varying(lg) if lg else (NM, [])
                pair = pairs_str(ii, vals, len(ov))
                d1.append([e, f"r{rep}", k, len(ov) + 1, f(mn), " · ".join(f"{i}–{i + 1}" for i in ii), comp, pair,
                           "std_analysis.json overlap_adj(=O[i,i+1] · MBAR 전량) · λ = pass1 md.log"])
    tb = list(csv.DictReader(open(RP / "s3_tableB_overlap.tsv", encoding="utf-8"), delimiter="\t"))
    for row in tb:
        x = s3m[row["rep"]][row["leg"]]
        chk(f"D:N17 {row['rep']} {row['leg']} tableB ov_sup_min = s3_mbar.ov_min", float(row["ov_sup_min"]), x["ov_min"], 5e-7)
        lam = x["lam"]; lo = x["ov_min_pair"][0]; i = lam.index(lo)
        d1.append(["N17", row["rep"], row["leg"], row["K"], f(x["ov_min"]), f"{i}–{i + 1}", "coul" if row["leg"][4:] != "L2" else "vdw",
                   f"{row['lam_lo']}–{row['lam_hi']}", "s3_tableB_overlap.tsv (ov_sup_min·lam_lo·lam_hi) · 전정밀 = s3_mbar_prod_v14_N17.json ov_min"])
    W(out, "D1_overlap.tsv",
      ["edge", "rep", "leg", "K", "겹침 1차 비대각 최솟값", "창 색인 쌍", "변하는 성분", "λ 쌍", "출처"], d1)

    # ═════════ D2 — 교환 ═════════
    d2 = []
    for e in N16:
        for k in LK:
            rx = replex[e].get(k, {})
            for rep in range(1, 6):
                lg = LOG.get((e, rep, k))
                if not lg: continue
                comp, vals = varying(lg); pr = lg.get("probs"); nx = lg.get("nexch")
                if not pr:
                    d2.append([e, f"r{rep}", k, lg["n_lambdas"], NM, NC, NM, NM, NM, NM, lg["attempts"], lg["n_summary_blocks_in_tail"], lg["log"]]); continue
                pv = [float(x) for x in pr]; mn, ii = argmins(pv)
                pair = pairs_str(ii, vals, len(pv))
                zero = [f"{vals[j]}–{vals[j + 1]}" if len(vals) == len(nx) + 1 else str(j) for j, c in enumerate(nx) if int(c) == 0] if nx else NM
                d2.append([e, f"r{rep}", k, lg["n_lambdas"], " ".join(pr), NC + "(회차별 평균 저장 안 됨)", f"{mn:.2f}", pair,
                           " ".join(nx) if nx else NM, ("없음" if zero == [] else (";".join(zero) if zero != NM else NM)),
                           f"{lg['attempts']} (odd {lg['odd']}·even {lg['even']})", lg["n_summary_blocks_in_tail"], lg["log"]])
            if rx:
                am = rx.get("acc_mean") or []
                lg0 = next((LOG[(e, r, k)] for r in range(1, 6) if (e, r, k) in LOG), None)
                comp, vals = varying(lg0) if lg0 else (NM, [])
                i = rx.get("acc_min_pair")
                pair = f"{vals[i]}–{vals[i + 1]}" if i is not None and len(vals) == len(am) + 1 else NM
                d2.append([e, f"회차 평균(nrep={rx.get('nrep')})", k, len(am) + 1, " ".join(f(x) for x in am), NC, f(rx.get("acc_min")), pair,
                           NC, NC, str(rx.get("n_attempt")), "—", "final_report_20260921/replex.json per_leg (쌍별 값을 회차에 걸쳐 평균한 저장값)"])
    c5 = v31c["C-5"]
    for r in range(1, 6):
        for k in LK:
            x = s3x[f"r{r}"][k]; rate = x["rate"]; lam = x["lam"]; mn, ii = argmins(rate)
            chk(f"D:N17 r{r} {k} s3_exchange 최소 = v31 쌍 최소", mn, c5["레그별 쌍 최소 (병기 — 판정 기준량 아님)"][f"r{r}"][k], 1e-12)
            d2.append(["N17", f"r{r}", k, len(lam), " ".join(f(v) for v in rate), f(c5["레그별 평균"][f"r{r}"][k]), f(mn), pairs_str(ii, lam, len(rate)),
                       " ".join(str(v) for v in x["acc"]) + " / 시도 " + " ".join(str(v) for v in x["att"]),
                       "없음" if not x["zero_pairs"] else json.dumps(x["zero_pairs"], ensure_ascii=False),
                       f"교환 블록 {x['nblock']} · 중복 {x['ndup']} · 기동 {x['nstart']}", "—",
                       "s3_exchange_prod_v14_N17.json {rate,acc,att,zero_pairs,lam} · 평균 = prereg v31 (c) C-5 레그별 평균"])
    W(out, "D2_exchange.tsv",
      ["edge", "rep", "leg", "K", "이웃 쌍별 수락률(λ 오름차순)", "이웃 수락률 평균", "이웃 수락률 최솟값", "그 λ 쌍",
       "쌍별 성사 수", "성사 0회 쌍", "시도 수", "log 끝 요약 블록 수", "출처"], d2)

    # ═════════ D3 — 왕복 ═════════
    d3 = []
    for e in N16:
        for k in LK:
            x = rescore[e].get(k); rt = rtrip[e].get(k, {})
            reps = [r for r in range(1, 6) if (e, r, k) in LOG]
            if x:
                chk(f"D:{e} {k} rescore.n_reps = 보유 log 수", x["n_reps"], len(reps))
                for j, r in enumerate(reps[:len(x["C_tot"])]):
                    d3.append([e, f"r{r}", k, x["K"], NH + "(워커별 값 저장 안 됨)", x["C_tot"][j], NC, x["A_tot"][j], x["nattempt"],
                               "rescore_v13.json {C_tot,A_tot} (회차 1..5 중 log 있는 순)"])
                d3.append([e, "회차 요약", k, x["K"], "—", f"C_min {x['C_min']} · C_mean {f(x['C_mean'])} · flag_C '{x['flag_C']}'", NC,
                           f"A_min {x['A_min']} · A_mean {f(x['A_mean'])} · flag_A '{x['flag_A']}'", x["nattempt"], "rescore_v13.json"])
            if rt:
                d3.append([e, "회차 요약(동결 9/21)", k, rt.get("K"), "—", "—", NC,
                           f"rt_mean {rt.get('rt_mean')} · rt_min {rt.get('rt_min')} · rt_max {rt.get('rt_max')} · flag '{rt.get('flag')}'",
                           "—", "final_report_20260921/roundtrip.json (nrep_run %s)" % rt.get("nrep_run")])
    z31 = v31c["사후 진단: 워커당 왕복 ≥1 (★사전등재 아님 — v29 ①)"]["레그별 왕복 0회 워커 수"]
    for r in range(1, 6):
        for k in LK:
            if r <= 4:
                t_ = hr1t[(f"r{r}", k)]; w_ = hr1w[(f"r{r}", k)]
                per, tot, zero, K = t_["rt_per"], t_["rt_total"], w_["rtzero"], t_["nwin"]
                src = "hr1_table_prod_v14_N17.json {rt_per,rt_total} · hr1_wide {rtzero}"
            else:
                t_ = s1r5[k]; per, tot, zero, K = t_["rt_per"], t_["rt_total"], t_["rt_zero"], t_["K"]
                src = "S1_r5_diag.json {rt_per,rt_total,rt_zero}"
            chk(f"D:N17 r{r} {k} 왕복 0회 워커 수 = v31", zero, z31[f"r{r}"][k])
            d3.append(["N17", f"r{r}", k, K, " ".join(str(v) for v in per), tot, zero, "—", "—", src])
    W(out, "D3_roundtrip.tsv",
      ["edge", "rep", "leg", "K", "워커별 왕복(워커 순)", "워커 왕복 합(C·복제축)", "왕복 0회 워커 수", "A(창위치축·옛 산정) 합/요약", "교환 시도 수", "출처"], d3)

    # ═════════ E — 네트워크 ═════════
    nodes = sorted({x for e in ALL for x in EDGE[e][:2]})
    W(out, "E1_graph.tsv", ["종류", "이름", "A", "B", "출처"],
      [["노드", n, "", "", "edges_v7.txt · edges_v12_add.txt"] for n in nodes] +
      [["엣지", e, EDGE[e][0], EDGE[e][1], EDGE[e][2]] for e in ALL])
    e2 = []
    for cid in sorted(signv["derived"], key=int):
        comp = " ".join(("+" if s > 0 else "−") + e for e, s in signv["derived"][cid].items())
        e2.append([f"고리{cid}", len(signv["derived"][cid]), comp, f(cyc_loc["closures"][cid]), f(cyc_loc["sprop"][cid]),
                   "sign_verification.json derived(구성·부호) · cycle_localization.json closures·sprop (동결 채점 · N01 n=3)"])
    e2.append(["앵커 삼각형(N01·N04·N07) · 채택 n=5 판", 3, NH + "(부호 목록이 이 파일에 없음)", f(lf_der["tri"]), f(lf_der["sprop"]),
               "pass1/local_fill_derived.json {tri,sprop} (tri_sigma %s)" % f(lf_der["tri_sigma"])])
    e2.append(["N17 이 드는 고리", NH, NH + "(구성 엣지 목록이 저장된 파일 없음 · N17_JUDGMENT_MATERIALS_20260928.md §6-6 은 개수만: 17간선 독립 사이클 8)",
               NC, NC, "prod_v14 N17 값을 넣은 닫힘 = 저장된 값 없음"])
    W(out, "E2_cycles.tsv", ["고리", "엣지 수", "구성 엣지(도는 부호)", "닫힘 값 (kcal/mol)", "전파 불확실도", "출처"], e2)
    src_true = None
    for i, ln in enumerate(open(Q / "src/edge_table_v11.py", encoding="utf-8"), 1):
        if ln.startswith("TRUE="): src_true = (i, ln.strip())
    TRUE = json.loads(src_true[1].split("=", 1)[1].replace("'", '"')) if src_true else {}
    e3 = []
    row0 = headline["rows"][0]
    for j, e in enumerate(("N01", "N04", "N07")):
        e3.append([e, f(TRUE.get(e)), f(row0["err"][j]), f"qsar/src/edge_table_v11.py:{src_true[0]} TRUE · 오차(계산−실험·동결 n=3) = headline.json rows[0].err[{j}]",
                   "Gupta 외, J. Med. Chem. 50, 5589 (2007) Table 1 — REPORT.md §6 기재"])
    for e in ALL:
        if e not in ("N01", "N04", "N07"): e3.append([e, NH, "—", "edge_table_v11.py TRUE 에 없음", "—"])
    W(out, "E3_exp.tsv", ["edge", "ΔΔG_exp (kcal/mol)", "오차(계산−실험)", "값의 출처", "문헌 출처(문서 기재)"], e3)

    # ═════════ F — 시작 구조와 난수 ═════════
    f1, f2, f2s = [], [], []
    tagdir = {}
    for pod, lst in assign["lists"].items():
        for j in lst: tagdir[j["tag"]] = (j["dir"], j.get("canon", ""), "lists")
    for j in assign["unresolved"]: tagdir[j["tag"]] = (j["dir"], j.get("canon", ""), "unresolved")

    def add_file(e, rep, s, role, p):
        p = Path(p)
        if p.exists() and p.is_file(): f1.append([e, rep, s, role, rel(p), p.stat().st_size, md5(p)])
        else: f1.append([e, rep, s, role, rel(p), NH, NM])

    def add_mdps(e, rep, s, d, canon):
        d = Path(d); kinds = {}
        files = sorted(d.glob("L[123]w[0-9][0-9]/seed.mdp")) + sorted(d.glob("L[123]w[0-9][0-9]/prod.mdp")) + sorted(d.glob("*.mdp"))
        if s == "cplx": files += sorted(d.parent.glob("*.mdp"))
        for p in files:
            if p.name in ("x.mdp", "xp.mdp", "mdout.mdp"): continue
            v = mdp_vals(p)
            if v is None: continue
            f2.append([e, rep, s, rel(p), v["gen_vel"], v["gen_seed"], v["ld_seed"], v["continuation"]])
            kinds.setdefault(p.name, []).append(v)
        for nm, vs in sorted(kinds.items()):
            f2s.append([e, rep, s, canon, nm, len(vs), "/".join(sorted({v["gen_vel"] for v in vs})),
                        "/".join(sorted({v["gen_seed"] for v in vs})), "/".join(sorted({v["ld_seed"] for v in vs}))])
        if not kinds: f2s.append([e, rep, s, canon, NH, 0, NM, NM, NM])

    for e in N16:
        for rep in range(1, 6):
            for s in ("cplx", "solv"):
                dirs = sorted({(tagdir[t][0], tagdir[t][1], tagdir[t][2]) for t in (f"{e}_{s}_L{l}_r{rep}" for l in (1, 2, 3)) if t in tagdir})
                if not dirs:
                    f1.append([e, f"r{rep}", s, "assignment.json 에 태그 없음", NH, NH, NM]); continue
                for d, canon, where in dirs:
                    d = Path(d); rr = f"r{rep}" + ("" if where == "lists" else " (assignment.unresolved)")
                    if s == "cplx":
                        add_file(e, rr, s, "c.gro (엣지 조립)", d.parent / "c.gro")
                        add_file(e, rr, s, "npt.gro (시딩 연쇄 진입 -c)", d.parent / "npt.gro")
                        add_file(e, rr, s, "npt.cpt (시딩 연쇄 진입 -t)", d.parent / "npt.cpt")
                    else:
                        add_file(e, rr, s, "lig.gro", d / "lig.gro")
                        add_file(e, rr, s, "npt.gro (시딩 연쇄 진입)", d / "npt.gro")
                    for l in (1, 2, 3): add_file(e, rr, s, f"L{l}w00/seed.gro (레그 {l} 첫 창 프로덕션 입력)", d / f"L{l}w00/seed.gro")
                    add_mdps(e, rr, s, d, canon)
    for rep in (1, 2, 3):   # N01 로컬 보충 (용매 r1·r2 · 대조 r3)
        d = B / f"local_fill_20260921/r{rep}"
        for l in (1, 2, 3): add_file("N01", f"r{rep} (local_fill_20260921)", "solv", f"L{l}w00/seed.gro", d / f"L{l}w00/seed.gro")
        add_mdps("N01", f"r{rep} (local_fill_20260921)", "solv", d, "local_fill")
    L3S = {2: B / "seeds_v14/L3entry_t1760.gro", 3: B / "seeds_v14/L3entry_t2000.gro"}
    for rep in range(1, 6):
        P = B / f"prod_v14_r{rep}/N17"; rr = f"r{rep}"
        add_file("N17", rr, "cplx", "c.gro (엣지 조립)", P / "cplx/c.gro")
        add_file("N17", rr, "cplx", "npt.gro (시딩 연쇄 진입 -c)", P / "cplx/npt.gro")
        add_file("N17", rr, "cplx", "npt.cpt (시딩 연쇄 진입 -t)", P / "cplx/npt.cpt")
        if rep in L3S: add_file("N17", rr, "cplx", "L3SEED (레그 3 진입 구조 교체 · n17_chain_v18.sh case)", L3S[rep])
        else: f1.append(["N17", rr, "cplx", "L3SEED", "미사용(연쇄 그대로 · n17_chain_v18.sh: 1|4|5) unset L3SEED)", "—", "—"])
        for l in (1, 2, 3): add_file("N17", rr, "cplx", f"L{l}w00/seed.gro (레그 {l} 첫 창 프로덕션 입력)", P / f"cplx/r{rep}/L{l}w00/seed.gro")
        add_mdps("N17", rr, "cplx", P / f"cplx/r{rep}", "prod_v14")
        for nm in ("lig.gro", "solv.gro", "em.gro", "nvt.gro", "npt.gro"):
            add_file("N17", rr, "solv", nm + (" (시딩 연쇄 진입)" if nm == "npt.gro" else ""), P / f"solv/r{rep}/{nm}")
        for l in (1, 2, 3): add_file("N17", rr, "solv", f"L{l}w00/seed.gro (레그 {l} 첫 창 프로덕션 입력)", P / f"solv/r{rep}/L{l}w00/seed.gro")
        add_mdps("N17", rr, "solv", P / f"solv/r{rep}", "prod_v14")
    W(out, "F1_input_coords.tsv", ["edge", "rep", "sys", "역할", "경로", "바이트", "md5"], f1)
    W(out, "F2_seeds_per_file.tsv", ["edge", "rep", "sys", "mdp 경로", "gen_vel", "gen_seed", "ld_seed", "continuation"], f2)
    W(out, "F2s_seeds_summary.tsv", ["edge", "rep", "sys", "정본 표기(assignment canon)", "mdp 종류", "파일 수", "gen_vel", "gen_seed", "ld_seed"], f2s)
    # 회차 간 md5 일치 요약 (같은 엣지·계·역할)
    grp = {}
    for e, rep, s, role, p, sz, h in f1:
        if h in (NM, "—"): continue
        grp.setdefault((e, s, role), []).append((rep, h))
    f1s = []
    for (e, s, role), v in sorted(grp.items()):
        hs = sorted({h for _, h in v})
        f1s.append([e, s, role, len(v), len(hs), "같다" if len(hs) == 1 and len(v) > 1 else ("다르다" if len(hs) == len(v) else "일부 같다"),
                    " ".join(f"{r}:{h[:8]}" for r, h in v)])
    W(out, "F1s_md5_across_reps.tsv", ["edge", "sys", "역할", "회차 수", "서로 다른 md5 수", "회차 간(md5)", "회차:md5 앞 8자"], f1s)

    # ═════════ A0 — 현행 엣지 목록 밖의 세대·이전 판 (원본 문구 옮김) ═════════
    a0 = []
    mdirs = {t: len(list((B / t).glob("M[0-9][0-9]"))) for t in ("stage0", "junction")}
    nprod_m = sum(1 for t in ("stage0", "junction") for _ in (B / t).rglob("prod.xvg"))
    medges = [ln.split() for ln in open(B / "edges_all.txt") if len(ln.split()) == 3]
    cl7 = next(c for c in prereg["changelog"] if c.get("version") == 7)
    a0.append(["M01~M17", "qsar/data/fep_m1/edges_all.txt", len(medges),
               f"stage0/M?? {mdirs['stage0']}개 · junction/M?? {mdirs['junction']}개 · 그 아래 prod.xvg {nprod_m}개",
               "prereg changelog v7 변경[0]: '" + cl7["변경"][0] + "' · fep_network_v5.엣지 = " + str(prereg["fep_network_v5"]["엣지"])])
    hold = [l.strip() for l in open(Q / "docs/STAGE1_HOLD_CONTAMINATED_FEP.md", encoding="utf-8")
            if ("구축된 FEP 계" in l or "**입력 오염**" in l or "dH/dλ 산출물" in l)]
    a0.append(["fep3b 세대 (C·E·P·PN 계열)", "qsar/docs/STAGE1_HOLD_CONTAMINATED_FEP.md §1 (2026-09-13)", NM + "(이번에 fep3b/ 내부를 열지 않음)",
               "fep3b/ 원본 무접촉", " / ".join(hold)])
    for v in ("prod_v12", "prod_v12b", "prod_v12c", "prod_v12d", "prod_v12e", "prod_v12f"):
        d = B / f"{v}_r1"; st = B / v / "status.txt"
        a0.append([f"N17 이전 판 {v}_r1", rel(d), 1, f".rep_done {'있음' if (d / '.rep_done').exists() else '없음'}",
                   "status.txt: " + (st.read_text(encoding="utf-8").strip() if st.exists() else NH)])
    W(out, "A0_other_generations.tsv", ["세대/판", "출처", "엣지 수", "디렉터리 사실", "원본 문구(옮김)"], a0)

    # ═════════ G — 정방향/역방향 메타 · μ2·τ2 ═════════
    fbl = list(csv.reader(open(RP / "s3_fb_legs.tsv", encoding="utf-8"), delimiter="\t"))
    fbd = list(csv.reader(open(RP / "s3_fb_ddg.tsv", encoding="utf-8"), delimiter="\t"))
    fr_l = sorted({r[2] for r in fbl[1:]}, key=float); fr_d = sorted({r[1] for r in fbd[1:]}, key=float)
    v29 = next(e for e in prereg["changelog"] if e.get("version") == 29)
    v29k = next(k for k in v29 if k.startswith("②"))
    nrow_set = sorted({n for r in s3m.values() for x in r.values() for n in x["nrow"]})
    tmax_set = sorted({x["tmax"] for r in uq1.values() for x in r.values()})
    g = [["qsar/reports/s3_fb_legs.tsv 헤더", " | ".join(fbl[0])],
         ["s3_fb_legs.tsv 자료 행 수", f"{len(fbl) - 1} (회차 {len({r[0] for r in fbl[1:]})} × 레그 {len({r[1] for r in fbl[1:]})} × 지점 {len(fr_l)})"],
         ["qsar/reports/s3_fb_ddg.tsv 헤더", " | ".join(fbd[0])],
         ["s3_fb_ddg.tsv 자료 행 수", f"{len(fbd) - 1}"],
         ["지점(data_fraction) 목록", "legs: " + " ".join(fr_l) + " · ddg: " + " ".join(fr_d)],
         ["50% 지점", "있음" if "0.5" in fr_l and "0.5" in fr_d else "없음"],
         ["생성 스크립트", "qsar/src/s3_mbar_v1.py:58-64 (alchemlyb forward_backward_convergence(us,'MBAR',num=4) → s3_mbar_prod_v14_N17.json 의 fb) → "
                      "qsar/src/s3_assemble_v1.py:38-63 (TSV 조립 · ΔΔG(f)=Σcplx 레그(f)−Σsolv 레그(f) · MBAR σ 쿼드러처)"],
         ["단위", "kcal/mol (s3_mbar_v1.py: fb 값 × kt2kcal)"],
         ["구간 정의 (prereg v29 ② 등재 문구)", v29[v29k]["★계산 전에 알고 있는 구현 사실 (alchemlyb 2.5.0 소스 직독)"]],
         ["분할 수 (prereg v29 ② 변경 2)", json.dumps(v29[v29k]["변경 2 — 정방향/역방향 분할 수"], ensure_ascii=False)],
         ["창당 행 수 · 길이 (원본 값)", f"nrow {nrow_set} (s3_mbar) · tmax {tmax_set} ps (uq1) · nstdhdl {sorted({d['nstdhdl'] for (e_, r_, k_), d in LOG.items() if e_ == 'N17'})} × dt {sorted({d['dt'] for (e_, r_, k_), d in LOG.items() if e_ == 'N17'})} ps (md.log)"],
         ["대상", "N17 · prod_v14 · r1~r5 뿐. N01~N16(패스1)에는 s3_fb 산출물이 없다 → 정방향/역방향 " + NC],
         ["μ2 · τ2 (prereg v28 (c) HA-4)", NC + " — prereg v29 ② 변경 1: " + json.dumps(v29[v29k]["변경 1 — μ2·τ2 미계산"], ensure_ascii=False)],
         ["★같은 종류가 아닌 양(참고로 위치만)", "qsar/reports/mix_v13.json — src/mix_v13.py 가 손으로 적은 쌍별 수락률(ACC)로 만든 생사 전이행렬의 t2(시도 단위). "
          "로그에서 만든 경험적 전이행렬의 μ2·τ2 가 아니다. 키: " + " / ".join(J(RP / "mix_v13.json").keys())]]
    W(out, "G1_fb_meta.tsv", ["항목", "내용"], g)

    # ═════════ REPORT.md — 위 행들에서 생성(재타이핑 없음) ═════════
    def md(header, rows):
        o = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
        o += ["| " + " | ".join(str(c).replace("|", "/") for c in r) + " |" for r in rows]
        return o + [""]

    def n4(x):
        try: return f"{float(x):+.4f}"
        except Exception: return str(x)

    def u4(x):
        try: return f"{float(x):.4f}"
        except Exception: return str(x)

    ngmx = 0
    for pd_ in Path("/proc").glob("[0-9]*"):
        try:
            if (pd_ / "comm").read_text().strip() in ("gmx", "gmx_mpi"): ngmx += 1
        except Exception: pass
    R = [f"스냅샷 {snap} · 재부팅 후 (uptime -s = {boot}) · gmx/gmx_mpi 프로세스 {ngmx}개(/proc/*/comm)", "",
         "# 네트워크 전수 데이터 추출 — 읽기 전용 · 새 계산 없음 · 새 판정 없음", "",
         f"생성 `qsar/src/net_extract_v1.py` → `{rel(out)}/` · 원본끼리의 대조 {len(CHECKS)}건(불일치 {sum(1 for c in CHECKS if c[1] != '일치')}) · 수치는 전부 아래 TSV 의 행에서 찍었다(표는 4자리 표시 · 전정밀은 TSV).",
         "표기: 미계산 = 계산된 적 없음 · 미보유 = 저장된 값/파일 없음 · 측정 불가 = 읽지 못함. 이 문서의 '통과/미달/수렴실패' 류 문구는 전부 원본 필드를 옮긴 것이다.", ""]
    # ---- A
    R += ["## A. 네트워크 현황", "",
          "엣지 목록 출처: `qsar/data/fep_m1/edges_v7.txt` (N01~N16 · prereg `fep_network_v7` '★확정' 16엣지) + `qsar/data/fep_m1/edges_v12_add.txt` (N17).",
          "상태 칸의 기준(이 추출에서 쓴 사실 기준): **완료** = 계획 회차의 프로덕션이 끝나 저장된 채점값이 있음 / 실행 중 = gmx 프로세스 있음(현재 0개) / 미시작 = 프로덕션 산출물 없음 / 폐기 = 원본에 폐기 등재 문구 있음.",
          "원본이 스스로 붙인 표기(✓완료·★연장후보·A1 판정·수렴 플래그)는 '원본 표기' 칸에 그대로 옮겼다. N01~N16 과 N17 은 prod 버전·창당 ns·λ 사다리가 다르다.", ""]
    rows = []
    for a in a1:
        e = a[0]
        cl = " ".join(sorted({"cl[" + r[2].split("]")[0].lstrip("[") + "]" for r in a3 if r[0] == e and r[1] == "changelog"},
                             key=lambda x: int(x[3:-1])))
        net = next(r[2] for r in a3 if r[0] == e and r[1].startswith("fep_network_v7"))
        rows.append([e, a[1], a[3], f"{a[6]} ({a[7]})", f"복합체 {a[8]} · 용매 {a[9]}", a[10], a[11],
                     f"fep_network_v7.엣지[{net}] · 키: {a[14]} · {cl}", a[5]])
    R += md(["엣지", "리간드 쌍", "상태", "회차 수 (채점 회차)", "레그 구성 n-lambdas L1/L2/L3", "창당 ns", "prod 버전", "prereg (이름 등장 항목)", "원본 표기(옮김)"], rows)
    R += ["prereg 항목 표기: `cl[i]` = `prereg.json` changelog 색인 i (버전·시각·키는 `A3_prereg_mentions.tsv`). 엣지 이름 없이 그물 전체에 걸리는 최상위 키: "
          + ", ".join(f"`{k}`" for k in COMMON) + ".", "",
          "상태 근거(엣지별): `A1_edges.tsv` '상태 근거' 칸. N17 의 근거: " + a1[-1][4] + ".", "",
          "### A-(e) 레그별 λ 실값 (md.log `all-lambdas` · 회차 전부 대조)", "",
          "출처: N01~N16 = `qsar/data/fep_m1/pass1/<TAG>/w00/prod.log` · N17 = `qsar/data/fep_m1/prod_v14_r<R>/N17/<계>/r<R>/L<l>w00/prod.log` (머리의 파라미터 덤프 = 첫 기동 것). "
          "레그 구성: L1 = coul-lambdas 만 변함 · L2 = vdw=bonded=fep-lambdas 가 변함 · L3 = coul-lambdas 만 변함(`A2_lambda_values.tsv` 에 네 성분 전부).", ""]
    R += md(["엣지", "레그", "log 있는 회차", "회차 간 λ 문자열", "n-lambdas", "변하는 성분", "λ 실값", "nsteps"],
            [[r[0], r[1], r[2], r[4], r[5], r[7], r[8], r[9]] for r in a2s])
    R += ["### A-0 현행 엣지 목록 밖의 세대·이전 판 (B~G 대상 아님)", ""]
    R += md(["세대/판", "출처", "엣지 수", "디렉터리 사실", "원본 문구(옮김)"], a0)
    # ---- B
    R += ["## B. 엣지별 ΔΔG (kcal/mol)", "",
          "출처: N01~N16 회차별 = `qsar/reports/c6_raw_v13.json` dd (= `preproc_sensitivity.json` per_cur · 대조 일치) · n·평균·s = `docs/stage1_out/pass1/final_edges.json` · SEM = `final_report_20260921/residuals_v2.json` se · t-95% 구간 = 미계산.",
          "N17 회차별 = `qsar/reports/S1_ddg_all.json` · 집합 통계·MBAR σ = `prereg.json` changelog v31 (b). N01~N16 의 추정기는 gmx bar(-b 0) 이고 ΔΔG 의 MBAR σ 는 미계산.", ""]
    main = {(r[0], r[1]): r for r in b1 if re.fullmatch(r"r\d", r[1])}
    rows = []
    for e in ALL:
        st = [r for r in b2 if r[0] == e][0] if e != "N17" else None
        vals = [n4(main[(e, f"r{i}")][2]) if (e, f"r{i}") in main else "—" for i in range(1, 6)]
        if e != "N17":
            rows.append([e, "gmx bar"] + vals + [st[2], n4(st[3]), u4(st[4]), u4(st[5]), NC, NC])
        else:
            for st in [r for r in b2 if r[0] == "N17"]:
                rows.append([e + " · " + st[1], "MBAR"] + vals + [st[2], st[3], st[4], st[5], st[6], f"[{st[7]}, {st[8]}]"])
    R += md(["엣지", "추정기", "r1", "r2", "r3", "r4", "r5", "n", "평균", "s", "SEM", "t", "95% 구간"], rows)
    R += ["N17 회차별 MBAR σ (보조 · prereg v31 (b) 인쇄값): " + " · ".join(f"{r[1]} {r[4]}" for r in b1 if r[0] == "N17") + ".", "",
          "저장된 다른 판(같은 엣지 · 채점/길이/회차 구성이 다름):", ""]
    R += md(["엣지", "집합/채점", "n", "평균", "s", "출처"], [[r[0], r[1], r[2], n4(r[3]), u4(r[4]), r[9]] for r in b2
                                                   if (r[0] != "N17" and not r[1].startswith("동결 채점"))])
    R += ["회차별 값 전부(표준 전처리판·로컬 보충·N04 2.0 ns 포함): `B1_ddg_per_rep.tsv` · 통계: `B2_ddg_stats.tsv`.", ""]
    # ---- C
    R += ["## C. 레그 분해 (kcal/mol)", "",
          "N17 — 출처 `qsar/reports/uq1_prod_v14_N17.json` 키 `r<R>.<레그>.{dG, mbar_sig}` (MBAR · 8 ns 전량).", ""]
    cN = {(r[1], r[2]): r for r in c1 if r[0] == "N17"}
    R += md(["레그"] + [f"r{i} ΔG ± MBAR σ" for i in range(1, 6)],
            [[k] + [f"{n4(cN[(f'r{i}', k)][5])} ± {u4(cN[(f'r{i}', k)][6])}" for i in range(1, 6)] for k in LK])
    R += ["회차별 레그 합(표시용 산술)과 저장된 그 회차 ΔΔG — 전 엣지. 출처: N01~N16 레그 = `final_report_20260921/current_bar.json` (gmx bar · 오차 미보유) · ΔΔG = `c6_raw_v13.json` dd / N17 = `uq1` · `S1_ddg_all.json`.", ""]
    R += md(["엣지", "회차", "추정기", "Σ복합체", "Σ용매", "Σ복합체−Σ용매", "저장된 ΔΔG", "보유 레그"],
            [[r[0], r[1], r[2], n4(r[3]), n4(r[4]), n4(r[5]), n4(r[6]), r[7]] for r in c2])
    R += [f"N01~N16 회차 × 레그 전부({sum(1 for r in c1 if r[0] != 'N17')}행 · gmx bar ΔG + MBAR 전량 ΔG·σ + MBAR 상관제거 ΔG·σ): `C1_leg_dG.tsv`. "
          "읽은 키: `current_bar.json[TAG]` · `std_analysis.json[edge][leg][TAG].{MBAR_raw, MBAR_raw_err, MBAR_dec, MBAR_dec_err}`.", ""]
    # ---- D
    R += ["## D. 교환 진단", "",
          "N17 (prod_v14) — 겹침 = `s3_tableB_overlap.tsv`/`s3_mbar_prod_v14_N17.json` · 수락 = `s3_exchange_prod_v14_N17.json` (평균은 prereg v31 (c) C-5 등재값) · "
          "왕복 = r1~r4 `hr1_table_prod_v14_N17.json`·`hr1_wide` / r5 `S1_r5_diag.json`.", ""]
    dO = {(r[1], r[2]): r for r in d1 if r[0] == "N17"}; dX = {(r[1], r[2]): r for r in d2 if r[0] == "N17"}
    dR = {(r[1], r[2]): r for r in d3 if r[0] == "N17"}
    rows = []
    for i in range(1, 6):
        for k in LK:
            o, x, t_ = dO[(f"r{i}", k)], dX[(f"r{i}", k)], dR[(f"r{i}", k)]
            rows.append([f"r{i}", k, o[3], u4(o[4]), o[7], u4(x[5]), u4(x[6]), x[7], x[9], t_[4], t_[5], t_[6]])
    R += md(["회차", "레그", "K", "겹침 1차 비대각 최소", "그 λ 쌍", "이웃 수락률 평균", "수락률 최소", "그 λ 쌍", "수락 0회 쌍", "워커별 왕복", "왕복 합", "왕복 0회 워커 수"], rows)
    n_d1 = sum(1 for r in d1 if r[0] != "N17"); n_d2 = sum(1 for r in d2 if r[0] != "N17" and r[1].startswith("r"))
    nz = sum(1 for r in d2 if r[0] != "N17" and r[1].startswith("r") and r[9] not in ("없음", NC))
    R += [f"N01~N16 (패스1) — 회차 × 레그 전부를 파일로: 겹침 `D1_overlap.tsv` ({n_d1}행 · `std_analysis.json` overlap_adj=O[i,i+1]) · "
          f"교환 `D2_exchange.tsv` ({n_d2}행 · GROMACS 가 `prod.log` 끝에 찍은 `Repl average probabilities`(2자리)·`number of exchanges` 를 그대로 · 회차별 평균은 미계산 · "
          f"수락 0회 쌍이 있는 행 {nz}개) + 회차 평균 저장값(`replex.json`) · 왕복 `D3_roundtrip.tsv` (`rescore_v13.json` 회차별 합 C·A · 워커별 값은 미보유 · 왕복 0회 워커 수는 미계산 · 동결 `roundtrip.json` 요약 병기).", ""]
    # ---- E
    R += ["## E. 네트워크 수준", "", f"노드 {len(nodes)}: " + " · ".join(nodes) + f". 엣지 {len(ALL)}: 위 A 표 (`E1_graph.tsv`).", ""]
    R += md(["고리", "엣지 수", "구성 엣지(부호)", "닫힘 값", "전파 불확실도", "출처"], [[r[0], r[1], r[2], n4(r[3]), u4(r[4]), r[5]] for r in e2])
    R += md(["엣지", "ΔΔG_exp", "오차(계산−실험·동결 n=3)", "값의 출처", "문헌 출처(문서 기재)"], [[r[0], n4(r[1]), n4(r[2]), r[3], r[4]] for r in e3 if r[1] != NH])
    R += ["그 밖의 엣지 " + " ".join(r[0] for r in e3 if r[1] == NH) + " 의 ΔΔG_exp = 미보유.", ""]
    # ---- F
    R += ["## F. 시작 구조와 난수", "",
          "N17 (prod_v14) — 회차 간 md5 (`F1_input_coords.tsv` 에 경로·바이트·md5 32자 전부):", ""]
    R += md(["계", "역할", "회차 수", "서로 다른 md5 수", "회차 간", "회차:md5 앞 8자"], [r[1:] for r in f1s if r[0] == "N17"])
    R += md(["회차", "계", "mdp 종류", "파일 수", "gen_vel", "gen_seed", "ld_seed"], [[r[1], r[2], r[4], r[5], r[6], r[7], r[8]] for r in f2s if r[0] == "N17"])
    cnt = {}
    for r in f2s:
        if r[0] == "N17": continue
        key = (r[2], r[4], r[6], r[7], r[8]); cnt[key] = cnt.get(key, 0) + 1
    R += ["N01~N16 (패스1 · 원 트리 `prod_v9c_r<R>/<엣지>/<계>/<정본 하위폴더>` = `assignment.json` 의 dir) — mdp 값 조합별 (엣지·회차·계) 묶음 수:", ""]
    R += md(["계", "mdp 종류", "gen_vel", "gen_seed", "ld_seed", "묶음 수"], [list(k) + [v] for k, v in sorted(cnt.items())])
    cn2 = {}
    for r in f1s:
        if r[0] == "N17": continue
        key = (r[1], r[2].split(" (")[0], r[5]); cn2[key] = cn2.get(key, 0) + 1
    R += ["N01~N16 입력 좌표의 회차 간 md5 — (계·역할·회차 간) 별 엣지 수 (엣지별 행은 `F1s_md5_across_reps.tsv` · 파일별 md5 는 `F1_input_coords.tsv`):", ""]
    R += md(["계", "역할", "회차 간(md5)", "엣지 수"], [list(k) + [v] for k, v in sorted(cn2.items())])
    t1 = []
    for i in range(1, 6):
        pth = B / f"prod_v14_r{i}/N17/solv/r{i}/lig.gro"
        t1.append(f"r{i} '" + (open(pth, encoding="utf-8", errors="replace").readline().strip() if pth.exists() else NH) + "'")
    R += ["N17 용매 `lig.gro` 의 첫 줄(제목 줄): " + " · ".join(t1) + ".", ""]
    l3 = []
    for i in range(1, 6):
        pth = B / f"prod_v14_r{i}/progress_prod_r{i}.log"
        ls_ = sorted({ln.split("]", 1)[1].strip() for ln in open(pth, encoding="utf-8", errors="replace") if "레그3 진입 구조 교체" in ln}) if pth.exists() else []
        l3.append(f"r{i}: " + (" / ".join(ls_) if ls_ else "해당 줄 없음"))
    R += ["`progress_prod_r<R>.log` 의 '레그3 진입 구조 교체' 줄: " + " · ".join(l3) + ".", ""]
    R += ["md.log 의 `ld-seed`(tpr 에 들어간 값): " + " · ".join(
        f"{lab_} {','.join(sorted({d['ld_seed'] for (e_, r_, k_), d in LOG.items() if sel(e_, r_)}))}"
        for lab_, sel in [("N01~N16 전 472 log", lambda e_, r_: e_ != "N17")] + [(f"N17 r{i}", (lambda i: (lambda e_, r_: e_ == "N17" and r_ == i))(i)) for i in range(1, 6)]) + ".", ""]
    # ---- G
    R += ["## G. 이미 계산된 것만", ""]
    R += md(["항목", "내용"], g)
    R += ["정방향/역방향 ΔΔG — `qsar/reports/s3_fb_ddg.tsv` 전 행(원본 그대로 · N17 · kcal/mol):", ""]
    R += md(fbd[0], fbd[1:])
    R += [f"엣지 × 회차 × 레그의 전 지점({len(fbl) - 1}행 · 50% 지점 포함): `qsar/reports/s3_fb_legs.tsv` (원본 · 복사하지 않았다).", ""]
    R += ["## 파일", ""] + [f"- `{rel(out)}/{nm}`" for nm in sorted(p.name for p in out.iterdir())] + [f"- `{rel(out)}/REPORT.md` · `CHECKS.tsv` · `META.json`", ""]
    with open(out / "REPORT.md", "x", encoding="utf-8") as fh:
        fh.write("\n".join(R) + "\n")

    # ═════════ 대조 결과 ═════════
    W(out, "CHECKS.tsv", ["대조", "결과", "a", "b", "허용오차"], CHECKS)
    nbad = sum(1 for c in CHECKS if c[1] != "일치")
    meta = dict(snapshot=snap, uptime_s=boot, n_checks=len(CHECKS), n_bad=nbad, outdir=rel(out),
                common_prereg_fep_keys=COMMON, n_log=len(LOG),
                script="qsar/src/net_extract_v1.py", script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with open(out / "META.json", "x", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=1)
    print(json.dumps(meta, ensure_ascii=False, indent=1))
    for nm in sorted(p.name for p in out.iterdir()): print("  ", nm)
    sys.exit(1 if nbad else 0)


if __name__ == "__main__":
    if len(sys.argv) != 2: print(__doc__); sys.exit(2)
    if sys.argv[1] == "--selftest": selftest()
    main(Path(sys.argv[1]).resolve())
