#!/usr/bin/env python3
"""조건 1(결정수 포함) 앵커 채점 — 동결 채점과 같은 방식: `gmx bar -b 0` (평형 버림 없음·상관제거 없음) · kJ→kcal ÷4.184.
  ΔΔG_r = Σ(새 복합체 레그 L1~L3, 회차 r) − Σ(기존 용매 레그 L1~L3, 회차 r)
  용매 레그는 다시 돌리지 않았다: current_bar.json 의 `<엣지>_solv_L<l>_r<r>` · N01 의 r1·r2 는 local_fill_legs.json new_solv(채택 n=5 판과 같은 값).
  통계: 회차 평균 · 표준편차(ddof=1) · 오차 = 계산 − 실험(src/edge_table_v11.py TRUE) · RMSE·MUE(앵커 3엣지).
  N01 공시험: prereg fep_null_test_v9 의 문턱 문구에 |ΔΔG| 를 대어 등급을 그대로 적는다(새 문턱 없음).
대조: (양성) 같은 bar 함수로 패스1 태그를 다시 채점해 current_bar.json 과 같아야 한다 · (음성) 창 하나를 빼면 달라져야 한다
      · 새 실행의 λ 벡터 = 패스1 md.log 의 λ · 종료 통계 250001 스텝 · 창 수 = n-lambdas.
사용: xw1_score_v1.py <뿌리(xw1)> <출력 디렉터리(없어야 한다)>
"""
import csv, glob, hashlib, json, math, os, re, statistics as st, subprocess, sys, tempfile
from pathlib import Path

T = Path("/home/nudge/Project/CADD/ttr"); Q = T / "qsar"; B = Q / "data/fep_m1"
FR = Q / "docs/stage1_out/final_report_20260921"; P1 = Q / "docs/stage1_out/pass1"
G = "/home/nudge/miniforge3/envs/md/bin/gmx"
REGKEY = os.environ.get("REGKEY", "anchor_rerun_condition1_crystal_water_v34")
EDGES = ["N01", "N04", "N07"]; NSTEP = 250000
os.environ["GMX_MAXBACKUP"] = "-1"
CHECKS = []


def chk(name, a, b, tol=0.0):
    try: ok = (abs(a - b) <= tol) if isinstance(a, (int, float)) and isinstance(b, (int, float)) else (a == b)
    except Exception: ok = False
    CHECKS.append((name, "일치" if ok else "★불일치", repr(a)[:100], repr(b)[:100], tol)); return ok


def J(p): return json.load(open(p, encoding="utf-8"))


def sha(p):
    h = hashlib.sha256(Path(p).read_bytes()).hexdigest(); assert len(h) == 64; return h


def bar(xs, W):
    """reproduce_edges_v11.py 의 bar() 와 같은 호출."""
    if len(xs) < 2: return None
    r = subprocess.run([G, "bar", "-f", *xs, "-b", "0", "-o", f"{W}/a.xvg", "-oi", f"{W}/i.xvg", "-oh", f"{W}/h.xvg"],
                       capture_output=True, text=True, timeout=1800)
    t = [l for l in r.stdout.splitlines() if l.startswith("total")]
    if not t: return None
    return float(t[-1].split("DG")[1].split("+/-")[0]) / 4.184


def log_lams(p):
    h = open(p, "rb").read(400_000).decode("utf-8", "replace"); i = h.find("all-lambdas:"); d = {}
    for c in ("coul", "vdw"):
        m = re.search(rf"^\s*{c}-lambdas =\s+(.*)$", h[i:i + 6000], re.M); d[c] = [float(x) for x in m.group(1).split()]
    return d


def main(root, out):
    pr = J(Q / "models/prereg.json")
    SELFTEST = os.environ.get("XW_SCORE_SELFTEST") == "1"   # ★자기시험: 패스1 원자료를 같은 배치로 넣어 동결값이 나오는지 본다(등재 가드 생략)
    if not SELFTEST:
        reg = [c for c in pr["changelog"] if c.get("키") == REGKEY]
        if len(reg) != 1: print(f"★prereg 에 {REGKEY} 항목이 {len(reg)}개 — 채점하지 않는다"); sys.exit(2)
        for s in reg[0]["(h)★스크립트"]:
            if sha(T / s["파일"]) != s["sha256"]: print("★등재 뒤 바뀐 스크립트:", s["파일"]); sys.exit(2)
    out.mkdir(parents=True, exist_ok=False)
    W = tempfile.mkdtemp(prefix="xw1bar_")
    cur = J(FR / "current_bar.json"); lf = J(P1 / "local_fill_legs.json"); c6 = J(Q / "reports/c6_raw_v13.json")
    fe = J(P1 / "final_edges.json"); lfr = J(P1 / "local_fill_result.json")
    src_true = [l for l in open(Q / "src/edge_table_v11.py", encoding="utf-8") if l.startswith("TRUE=")][0]
    TRUE = json.loads(src_true.split("=", 1)[1].replace("'", '"'))
    # ── 대조: 같은 bar 로 패스1 태그 재채점
    for e in EDGES:
        tag = f"{e}_cplx_L1_r3"; xs = sorted(glob.glob(str(B / f"pass1/{tag}/w*/prod.xvg")))
        v = bar(xs, W); chk(f"양성:{tag} 재채점 = current_bar.json", v, cur[tag], 1e-9)
        v2 = bar(xs[:-1], W) if len(xs) > 2 else None
        if v2 is not None: chk(f"음성:{tag} 창 하나를 빼면 달라진다", abs(v2 - cur[tag]) > 1e-6, True)
    legs, ddg = [], []
    per = {e: {} for e in EDGES}
    for e in EDGES:
        for r in range(1, 6):
            C = root / f"prod_r{r}/{e}/cplx"
            if not (C / "CANON").exists():
                ddg.append([e, f"r{r}", "미완(CANON 없음)" + (" · 재시도에서도 게이트 미달" if (C / ".chain_failed").exists() else ""), "", "", "", "", ""]); continue
            sub = (C / "CANON").read_text().strip(); v = {}
            for l in (1, 2, 3):
                ws = sorted(glob.glob(str(C / sub / f"L{l}w[0-9][0-9]")))
                xs = [w + "/prod.xvg" for w in ws]
                ref = log_lams(B / f"pass1/{e}_cplx_L{l}_r3/w00/prod.log")
                mdp = {m.group(1): [float(x) for x in m.group(2).split()] for m in
                       (re.match(r"\s*(coul|vdw)-lambdas\s*=\s*(.*)$", ln) for ln in open(ws[0] + "/prod.mdp")) if m}
                chk(f"{e} r{r} L{l} 창 수 = 패스1 n-lambdas", len(ws), len(ref["coul"]))
                for c in ("coul", "vdw"):
                    chk(f"{e} r{r} L{l} {c}-lambdas = 패스1 md.log", all(abs(a - b) < 5e-6 for a, b in zip(mdp[c], ref[c])) and len(mdp[c]) == len(ref[c]), True)
                for w in ws:
                    tl = open(w + "/prod.log", "rb").read()[-200_000:].decode("utf-8", "replace")
                    m = re.findall(r"Statistics over (\d+) steps", tl)
                    chk(f"{e} r{r} L{l} {Path(w).name} 종료 통계 = {NSTEP + 1}", int(m[-1]) if m else None, NSTEP + 1)
                v[l] = bar(xs, W)
                old = cur.get(f"{e}_cplx_L{l}_r{r}")
                legs.append([e, f"r{r}", f"cplxL{l}", sub, len(ws), repr(v[l]), repr(old) if old is not None else "미보유",
                             repr(v[l] - old) if (old is not None and v[l] is not None) else ""])
            if e == "N01" and r in (1, 2):
                sv = [lf["new_solv"][str(r)][str(l)] for l in (1, 2, 3)]; ssrc = "local_fill_legs.json new_solv"
                oldd = lfr["n2"][r - 1]
            else:
                sv = [cur.get(f"{e}_solv_L{l}_r{r}") for l in (1, 2, 3)]; ssrc = "current_bar.json"
                oldd = c6[e]["dd"].get(str(r))
            if any(x is None for x in sv) or any(v[l] is None for l in v):
                ddg.append([e, f"r{r}", "용매 레그 값 없음" if any(x is None for x in sv) else "bar 실패", "", "", "", "", ""]); continue
            sc, ss = sum(v.values()), sum(sv); d = sc - ss; per[e][r] = d
            ddg.append([e, f"r{r}", sub, repr(sc), repr(ss), repr(d), repr(oldd) if oldd is not None else "미보유", ssrc])
    summ, errs_new = [], {}
    for e in EDGES:
        vals = [per[e][r] for r in sorted(per[e])]
        n = len(vals); mean = st.mean(vals) if vals else None; sd = st.stdev(vals) if n > 1 else None
        exp = TRUE[e]; err = (mean - exp) if mean is not None else None; errs_new[e] = err
        summ.append([e, n, " ".join(f"r{r}" for r in sorted(per[e])), repr(mean) if mean is not None else "미계산", repr(sd) if sd is not None else "미계산",
                     repr(exp), repr(err) if err is not None else "미계산", fe[e]["n"], repr(fe[e]["dd"]), repr(fe[e]["sd"]), repr(fe[e]["dd"] - exp)])
    nt = pr["fep_null_test_v9"]["문턱"]; a = abs(errs_new["N01"]) if errs_new["N01"] is not None else None
    if a is None: cls = "미계산"
    elif a <= 0.5: cls = "≤0.5 → " + nt["≤0.5"]
    elif a <= 1.0: cls = "0.5~1.0 → " + nt["0.5~1.0"]
    elif a <= 2.0: cls = ">1.0 → " + nt[">1.0"]
    else: cls = ">2.0 → " + nt[">2.0"] + " (>1.0: " + nt[">1.0"] + ")"
    ok3 = all(errs_new[e] is not None for e in EDGES)
    rmse = math.sqrt(sum(errs_new[e] ** 2 for e in EDGES) / 3) if ok3 else None
    mue = sum(abs(errs_new[e]) for e in EDGES) / 3 if ok3 else None
    hl = J(FR / "headline.json")["rows"][0]

    def Wt(name, header, rows):
        with open(out / name, "x", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t"); w.writerow(header); w.writerows(rows)
    Wt("X1_legs.tsv", ["edge", "rep", "leg", "정본 하위폴더", "창 수", "ΔG 새 실행 (kcal/mol · gmx bar)", "ΔG 패스1 같은 태그", "차(새 − 패스1)"], legs)
    Wt("X2_ddg.tsv", ["edge", "rep", "정본 하위폴더/상태", "Σ복합체(새)", "Σ용매(기존)", "ΔΔG 새", "ΔΔG 패스1 같은 회차", "용매 값 출처"], ddg)
    Wt("X3_summary.tsv", ["edge", "n(새)", "회차", "평균 ΔΔG(새)", "s(새 · ddof=1)", "실험 ΔΔG", "오차(새) = 계산 − 실험", "n(동결)", "평균(동결)", "s(동결)", "오차(동결)"], summ)
    Wt("CHECKS.tsv", ["대조", "결과", "a", "b", "허용오차"], CHECKS)
    nb = sum(1 for c in CHECKS if c[1] != "일치")
    snap = subprocess.run(["date", "-Is"], capture_output=True, text=True).stdout.strip()
    boot = subprocess.run(["uptime", "-s"], capture_output=True, text=True).stdout.strip()
    f3 = lambda x: "미계산" if x is None else f"{float(x):+.3f}"; u3 = lambda x: "미계산" if x is None else f"{float(x):.3f}"
    R = [f"스냅샷 {snap} · 재부팅 후 (uptime -s = {boot})", "", f"# 조건 1(결정수 포함) — 앵커 3엣지 채점 (prereg {REGKEY})", "",
         f"생성 `qsar/src/xw1_score_v1.py` (sha256 `{sha(__file__)}`) · 대조 {len(CHECKS)}건(불일치 {nb}) · 단위 kcal/mol · 3자리 표시(전정밀 TSV).", "",
         "| 엣지 | n | 평균 ΔΔG | s | 실험 | 오차(새) | 동결 n | 동결 평균 | 동결 s | 오차(동결) |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in summ:
        R.append(f"| {r[0]} | {r[1]} ({r[2]}) | {f3(None if r[3] == '미계산' else r[3])} | {u3(None if r[4] == '미계산' else r[4])} | {f3(r[5])} | "
                 f"{f3(None if r[6] == '미계산' else r[6])} | {r[7]} | {f3(r[8])} | {u3(r[9])} | {f3(r[10])} |")
    R += ["", f"앵커 3엣지 RMSE {u3(rmse)} · MUE {u3(mue)} (동결 n=3 채점: RMSE {hl['rmse']:.3f} · MUE {hl['mue']:.3f} — headline.json rows[0]).",
          f"N01 공시험 |ΔΔG| = {u3(a)} → prereg fep_null_test_v9 문턱: {cls}.", "",
          "회차별 값 `X2_ddg.tsv` · 레그별 값과 패스1 같은 태그와의 차 `X1_legs.tsv`.", ""]
    with open(out / "REPORT.md", "x", encoding="utf-8") as fh: fh.write("\n".join(R) + "\n")
    meta = dict(at=snap, selftest=SELFTEST, n_checks=len(CHECKS), n_bad=nb, script_sha256=sha(__file__), regkey=REGKEY, rmse=rmse, mue=mue,
                n01_abs=a, n01_class=cls, per_rep=per)
    with open(out / "META.json", "x", encoding="utf-8") as fh: json.dump(meta, fh, ensure_ascii=False, indent=1)
    print("\n".join(R)); sys.exit(1 if nb else 0)


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
