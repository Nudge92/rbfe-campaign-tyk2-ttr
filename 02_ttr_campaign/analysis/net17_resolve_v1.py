#!/usr/bin/env python3
"""1단계 — N17 을 넣고 그물(노드 전위)을 다시 푼다. ★새 시뮬레이션 0 · 저장된 엣지 값 위의 산술만.

정의는 9/21 후보 표(docs/stage1_out/pass1/CANDIDATE_TABLE_20260921.md)와 같다 — 16엣지 재현 대조로 확인한다:
  · ΔΔG(A→B) = G_B − G_A · 게이지 G(CHEMBL241454) = 0
  · 무가중(OLS): G = pinv(A)·y · 전파 σ = sqrt(diag(M·diag(σ̂²)·Mᵀ)), M = pinv(A)
  · 가중(GLS):   W = diag(1/σ̂²) · G = (AᵀWA)⁻¹AᵀW·y · 전파 σ = sqrt(diag((AᵀWA)⁻¹))
  · σ̂ = 그 엣지의 회차 간 표준편차(ddof=1) · χ² = Σ(잔차/σ̂)² · 잔차 = 관측 − 적합 · dof = 엣지 수 − rank(A)
  · 잔차 스케일 σ = 전파 σ × sqrt(χ²/dof)
  · 고리 닫힘 = Σ 부호·ΔΔG · 전파 불확실도 = sqrt(Σ σ̂²)
판정 문구(구분 불가·유의·닫힘/미닫힘 등)는 붙이지 않는다 — 값만 낸다.

모드:
  --repro16 [--json <경로>]  N17 없이 16엣지만 — 9/21 저장값 재현 대조. 종료코드 0(전부 일치)/1.
  --dry <출력 디렉터리>      ★합성 자료 시험: 임의의 노드 값에서 만든 무잡음 엣지 값으로 전 경로 실행(실제 엣지 값·N17 값 미사용).
  --run <출력 디렉터리>      본 실행. prereg changelog 마지막 항목의 키가 REGKEY 이고 거기 적힌 이 스크립트 sha256 이
                            자기 것과 같아야 한다(등재 전 실행·등재 후 수정 차단). 출력 디렉터리는 없어야 한다.
환경변수 NET17_FE=<final_edges.json 변조 사본> : 음성 대조용(기본 = 원본).
"""
import csv, hashlib, json, math, os, re, subprocess, sys
from pathlib import Path
import numpy as np

T = Path("/home/nudge/Project/CADD/ttr"); Q = T / "qsar"; B = Q / "data/fep_m1"
P1 = Q / "docs/stage1_out/pass1"; FR = Q / "docs/stage1_out/final_report_20260921"
PREREG = Q / "models/prereg.json"
FE_PATH = Path(os.environ.get("NET17_FE", P1 / "final_edges.json"))
REGKEY = "network_resolve_with_n17_step1_v33"
REF = "CHEMBL241454"
ANCHORS = ["CHEMBL240808", "CHEMBL438498"]
CANDS = ["cand0023", "cand0028", "cand0030", "cand0042", "cand0118", "cand0132", "cand0271"]
CONFIGS = [("A", "17엣지 전부", []), ("B", "N14 제외 (16엣지)", ["N14"]), ("C", "N08·N11·N14 제외 (14엣지)", ["N08", "N11", "N14"])]
FITS = [("OLS", "무가중"), ("GLS", "1/σ̂² 가중")]
NEWCYC = {"8": {"N17": 1, "N16": 1, "N11": -1, "N05": -1},
          "9": {"N17": 1, "N14": -1, "N10": -1, "N08": 1, "N07": -1},
          "10": {"N17": 1, "N16": 1, "N15": 1, "N13": -1, "N10": -1, "N08": 1, "N07": -1}}

EDGE = {}
for _fn in ("edges_v7.txt", "edges_v12_add.txt"):
    for _ln in open(B / _fn):
        _q = _ln.split()
        if len(_q) == 3: EDGE[_q[0]] = (_q[1], _q[2])
N16 = [f"N{i:02d}" for i in range(1, 17)]; ALL17 = N16 + ["N17"]
assert sorted(EDGE) == ALL17, sorted(EDGE)

CHECKS = []


def chk(name, a, b, tol=0.0):
    try:
        ok = (abs(a - b) <= tol) if isinstance(a, (int, float)) and isinstance(b, (int, float)) else (a == b)
    except Exception:
        ok = False
    CHECKS.append((name, "일치" if ok else "★불일치", repr(a), repr(b), tol)); return ok


def J(p): return json.load(open(p, encoding="utf-8"))


def sha(p):
    h = hashlib.sha256(Path(p).read_bytes()).hexdigest(); assert len(h) == 64; return h


def rel(p):
    try: return str(Path(p).relative_to(T))
    except ValueError: return str(p)


def solve(names, val, sdv, wmode):
    """names: 엣지 이름 목록 · val/sdv: 이름→ΔΔG/σ̂. 돌려주는 G·σ 는 게이지(REF=0) 기준."""
    nodes = sorted({x for e in names for x in EDGE[e]} - {REF})
    A = np.zeros((len(names), len(nodes)))
    for i, e in enumerate(names):
        a, b = EDGE[e]
        if b != REF: A[i, nodes.index(b)] += 1.0
        if a != REF: A[i, nodes.index(a)] -= 1.0
    y = np.array([val[e] for e in names], float); sd = np.array([sdv[e] for e in names], float)
    rank = int(np.linalg.matrix_rank(A)); assert rank == len(nodes), f"그물이 끊겼다: rank {rank} < 노드 {len(nodes)}"
    if wmode == "OLS":
        M = np.linalg.pinv(A); G = M @ y; C = M @ np.diag(sd ** 2) @ M.T
    elif wmode == "GLS":
        W = np.diag(1.0 / sd ** 2); N = np.linalg.inv(A.T @ W @ A); G = N @ A.T @ W @ y; C = N
    else:
        raise ValueError(wmode)
    fit = A @ G; res = y - fit
    chi2 = float(((res / sd) ** 2).sum()); dof = len(names) - rank
    scale = math.sqrt(chi2 / dof) if dof > 0 else float("nan")
    sp = np.sqrt(np.diag(C))
    return dict(nodes=nodes, G=dict(zip(nodes, map(float, G))), sprop=dict(zip(nodes, map(float, sp))),
                sscaled={n: float(s * scale) for n, s in zip(nodes, sp)}, fit=dict(zip(names, map(float, fit))),
                res=dict(zip(names, map(float, res))), chi2=chi2, dof=dof, scale=scale, n_edges=len(names))


def cycle_ok(c):
    bal = {}
    for e, s in c.items():
        a, b = EDGE[e]; bal[b] = bal.get(b, 0) + s; bal[a] = bal.get(a, 0) - s
    return all(v == 0 for v in bal.values())


def closure(c, val, sdv):
    return float(sum(s * val[e] for e, s in c.items())), math.sqrt(sum(sdv[e] ** 2 for e in c))


def num(s):
    return float(s.replace("−", "-").replace("+", "").replace("*", "").replace("★", "").strip())


def repro16():
    """N17 없이 16엣지 — 9/21 저장값과 대조. 실제 N17 값은 읽지 않는다."""
    FE = J(FE_PATH); val = {e: FE[e]["dd"] for e in N16}; sdv = {e: FE[e]["sd"] for e in N16}
    RV = J(FR / "residuals_v2.json"); CL = J(FR / "cycle_localization.json"); SV = J(FR / "sign_verification.json")["derived"]
    R = {w: solve(N16, val, sdv, w) for w, _ in FITS}
    for e in N16:
        chk(f"재현:{e} OLS 적합 = residuals_v2.fit_u", R["OLS"]["fit"][e], RV[e]["fit_u"], 1e-9)
        chk(f"재현:{e} OLS 잔차 = residuals_v2.res_u", R["OLS"]["res"][e], RV[e]["res_u"], 1e-9)
        chk(f"재현:{e} GLS 적합 = residuals_v2.fit_w", R["GLS"]["fit"][e], RV[e]["fit_w"], 1e-9)
        chk(f"재현:{e} GLS 잔차 = residuals_v2.res_w", R["GLS"]["res"][e], RV[e]["res_w"], 1e-9)
    for cid, c in SV.items():
        chk(f"재현:고리{cid} 부호가 닫힌 경로", cycle_ok(c), True)
        cs, sp = closure(c, val, sdv)
        chk(f"재현:고리{cid} 닫힘 = cycle_localization", cs, CL["closures"][cid], 1e-9)
        chk(f"재현:고리{cid} 전파 = cycle_localization", sp, CL["sprop"][cid], 1e-9)
    # 후보 표 §3·§4 의 인쇄값(2자리)
    txt = open(P1 / "CANDIDATE_TABLE_20260921.md", encoding="utf-8").read()
    pm = r"\*\*([−+\-]?[\d.]+) ± ([\d.]+)\*\*\s*\|\s*\*([\d.]+)\*"
    nrow = 0
    for nd in ANCHORS + CANDS:
        m = re.search(r"^\|\s*\*{0,2}" + nd + r"\*{0,2}[^|]*\|\s*★?" + pm + r"\s*\|\s*★?" + pm, txt, re.M)
        if not chk(f"재현:후보 표 §4 에 {nd} 행이 있다", m is not None, True): continue
        nrow += 1; g = [num(x) for x in m.groups()]
        for k, w in ((0, "OLS"), (3, "GLS")):
            chk(f"재현:{nd} {w} G(2자리)", round(R[w]["G"][nd], 2), g[k], 0.0051)
            chk(f"재현:{nd} {w} 잔차 스케일 σ(2자리)", round(R[w]["sscaled"][nd], 2), g[k + 1], 0.0051)
            chk(f"재현:{nd} {w} 전파 σ(2자리)", round(R[w]["sprop"][nd], 2), g[k + 2], 0.0051)
    for lab, w in (("16간선 · 무가중 OLS", "OLS"), ("16간선 · 1/σ² 가중 GLS", "GLS")):
        m = re.search(r"^\|\s*" + re.escape(lab) + r"\s*\|\s*(\d+)\s*\|\s*([\d.]+)\s*\|\s*(\d+)\s*\|\s*★?\*\*([\d.]+)\*\*", txt, re.M)
        if not chk(f"재현:후보 표 §3 에 '{lab}' 행이 있다", m is not None, True): continue
        chk(f"재현:{w} 엣지 수", R[w]["n_edges"], int(m.group(1))); chk(f"재현:{w} χ²(2자리)", round(R[w]["chi2"], 2), float(m.group(2)), 0.0051)
        chk(f"재현:{w} dof", R[w]["dof"], int(m.group(3))); chk(f"재현:{w} χ²/dof(2자리)", round(R[w]["chi2"] / R[w]["dof"], 2), float(m.group(4)), 0.0051)
    return R, val, sdv


def W_(out, name, header, rows):
    with open(out / name, "x", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(header); w.writerows(rows)


def md(header, rows):
    o = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    return o + ["| " + " | ".join(str(c).replace("|", "/") for c in r) + " |" for r in rows] + [""]


def s3(x): return f"{x:+.3f}"
def u3(x): return f"{x:.3f}"


def compute(out, val, sdv, nv, n17sets, label, ref16=None):
    """전 조합을 풀고 TSV·REPORT.md 를 쓴다. n17sets: [(이름, 평균, s, n)]."""
    t1, t2, t3, t4, t5, t7 = [], [], [], [], [], []
    ranks = {}; res_all = {}
    for sname, m17, s17, n17 in n17sets:
        v = dict(val); s = dict(sdv); v["N17"] = m17; s["N17"] = s17
        for cid, c in NEWCYC.items():
            assert cycle_ok(c), cid
            cs, sp = closure(c, v, s)
            t4.append([sname, f"고리{cid}", len(c), " ".join(("+" if sg > 0 else "−") + e for e, sg in c.items()), repr(cs), repr(sp),
                       repr(abs(cs) / sp)])
        for ck, cname, drop in CONFIGS:
            names = [e for e in ALL17 if e not in drop]
            for w, wname in FITS:
                r = solve(names, v, s, w); key = (sname, ck, w); res_all[key] = r
                t2.append([sname, ck, cname, w, r["n_edges"], r["dof"], repr(r["chi2"]), repr(r["chi2"] / r["dof"]), repr(r["scale"])])
                order = sorted(CANDS, key=lambda nd: r["G"][nd]); ranks[key] = {nd: i + 1 for i, nd in enumerate(order)}
                for nd in ANCHORS + CANDS:
                    t1.append([sname, ck, w, nd, "앵커" if nd in ANCHORS else "후보", repr(r["G"][nd]), repr(r["sprop"][nd]), repr(r["sscaled"][nd]),
                               ranks[key].get(nd, "—")])
                for e in names:
                    t3.append([sname, ck, w, e, f"{EDGE[e][0]} → {EDGE[e][1]}", repr(v[e]), repr(r["fit"][e]), repr(r["res"][e]), repr(s[e]),
                               repr(r["res"][e] / s[e])])
                t5.append([sname, ck, w] + order)
                for i, a in enumerate(CANDS):
                    for b in CANDS[i + 1:]:
                        d = r["G"][a] - r["G"][b]; sg = math.hypot(r["sscaled"][a], r["sscaled"][b])
                        t7.append([sname, ck, w, a, b, repr(d), repr(sg), repr(abs(d) / sg)])
    W_(out, "T1_nodes.tsv", ["N17 집합", "구성", "적합", "노드", "종류", "G (kcal/mol·게이지 CHEMBL241454=0)", "전파 σ", "잔차 스케일 σ", "후보 중 순위(G 오름차순)"], t1)
    W_(out, "T2_fit_quality.tsv", ["N17 집합", "구성", "구성 설명", "적합", "엣지 수", "dof", "χ²", "χ²/dof", "스케일 인자 sqrt(χ²/dof)"], t2)
    W_(out, "T3_edge_residuals.tsv", ["N17 집합", "구성", "적합", "엣지", "노드 쌍", "관측 ΔΔG", "적합", "잔차(관측−적합)", "σ̂", "잔차/σ̂"], t3)
    W_(out, "T4_cycles_with_N17.tsv", ["N17 집합", "고리", "엣지 수", "구성 엣지(부호)", "닫힘 값", "전파 불확실도 sqrt(Σσ̂²)", "|닫힘|/전파"], t4)
    W_(out, "T5_candidate_order.tsv", ["N17 집합", "구성", "적합", "1위", "2위", "3위", "4위", "5위", "6위", "7위"], t5)
    W_(out, "T7_candidate_pairs.tsv", ["N17 집합", "구성", "적합", "후보 a", "후보 b", "G_a − G_b", "σ합 = hypot(잔차 스케일 σ)", "|차|/σ합"], t7)
    t6 = []
    if ref16 is not None:
        for w, _ in FITS:
            o16 = sorted(CANDS, key=lambda nd: ref16[w]["G"][nd])
            for nd in ANCHORS + CANDS:
                row = [w, nd, repr(ref16[w]["G"][nd]), repr(ref16[w]["sscaled"][nd]), (o16.index(nd) + 1) if nd in CANDS else "—"]
                for sname, *_ in n17sets:
                    for ck, _, _ in CONFIGS:
                        r = res_all[(sname, ck, w)]; row += [repr(r["G"][nd]), repr(r["sscaled"][nd]), ranks[(sname, ck, w)].get(nd, "—")]
                t6.append(row)
        hdr = ["적합", "노드", "9/21 16엣지 G", "9/21 잔차 스케일 σ", "9/21 순위"]
        for sname, *_ in n17sets:
            for ck, _, _ in CONFIGS: hdr += [f"{sname}·{ck} G", f"{sname}·{ck} σ", f"{sname}·{ck} 순위"]
        W_(out, "T6_compare_with_0921.tsv", hdr, t6)

    # ── REPORT.md (위 행에서 생성)
    snap = subprocess.run(["date", "-Is"], capture_output=True, text=True).stdout.strip()
    boot = subprocess.run(["uptime", "-s"], capture_output=True, text=True).stdout.strip()
    ngmx = 0
    for pd_ in Path("/proc").glob("[0-9]*"):
        try:
            if (pd_ / "comm").read_text().strip() in ("gmx", "gmx_mpi"): ngmx += 1
        except Exception: pass
    nb = sum(1 for c in CHECKS if c[1] != "일치")
    R = [f"스냅샷 {snap} · 재부팅 후 (uptime -s = {boot}) · gmx/gmx_mpi 프로세스 {ngmx}개(/proc/*/comm)", "",
         f"# {label}", "",
         f"생성 `qsar/src/net17_resolve_v1.py` (sha256 `{sha(__file__)}`) → `{rel(out)}/` · 대조 {len(CHECKS)}건(불일치 {nb}). 표는 3자리 표시 · 전정밀은 TSV.",
         "단위 kcal/mol · 게이지 G(CHEMBL241454) = 0 · G 가 낮을수록 강한 결합 · ± 는 잔차 스케일 σ(괄호 = 전파 σ) · 판정 문구 없음(값만).", "",
         "## 입력", ""]
    R += md(["엣지", "노드 쌍", "ΔΔG", "σ̂", "n"], [[e, f"{EDGE[e][0]} → {EDGE[e][1]}", s3(val[e]), u3(sdv[e]), nv[e]] for e in N16] +
            [[f"N17 · {sn}", f"{EDGE['N17'][0]} → {EDGE['N17'][1]}", s3(m), u3(s_), n_] for sn, m, s_, n_ in n17sets])
    R += ["## 후보 순서 (G 오름차순)", ""]
    R += md(["N17 집합", "구성", "적합", "1위", "2위", "3위", "4위", "5위", "6위", "7위"], t5)
    if ref16 is not None:
        R += ["9/21 16엣지(N17 없음) 순서: " + " · ".join(f"{w}: " + " < ".join(sorted(CANDS, key=lambda nd: ref16[w]["G"][nd])) for w, _ in FITS) + ".", ""]
    R += ["## 후보별 순위 (조합별)", ""]
    keys = [(sn, ck, w) for sn, *_ in n17sets for ck, _, _ in CONFIGS for w, _ in FITS]
    R += md(["후보"] + [f"{sn}·{ck}·{w}" for sn, ck, w in keys], [[nd] + [ranks[k][nd] for k in keys] for nd in CANDS])
    R += ["## 노드 값", ""]
    for sn, *_ in n17sets:
        for ck, cname, _ in CONFIGS:
            R += [f"### N17 = {sn} · 구성 {ck} ({cname})", ""]
            rows = []
            for nd in ANCHORS + CANDS:
                row = [nd]
                for w, _ in FITS:
                    r = res_all[(sn, ck, w)]
                    row += [f"{r['G'][nd]:+.3f} ± {r['sscaled'][nd]:.3f} ({r['sprop'][nd]:.3f})", ranks[(sn, ck, w)].get(nd, "—")]
                if ref16 is not None: row += [f"{ref16['OLS']['G'][nd]:+.3f} ± {ref16['OLS']['sscaled'][nd]:.3f}", f"{ref16['GLS']['G'][nd]:+.3f} ± {ref16['GLS']['sscaled'][nd]:.3f}"]
                rows.append(row)
            R += md(["노드", "OLS G ± σ (전파)", "순위", "GLS G ± σ (전파)", "순위"] + (["9/21 OLS", "9/21 GLS"] if ref16 is not None else []), rows)
    R += ["## 적합도", ""]
    R += md(["N17 집합", "구성", "적합", "엣지 수", "dof", "χ²", "χ²/dof", "스케일 인자"], [[r[0], r[1], r[3], r[4], r[5], u3(float(r[6])), u3(float(r[7])), u3(float(r[8]))] for r in t2])
    if ref16 is not None:
        R += ["9/21 16엣지: " + " · ".join(f"{w} χ² {ref16[w]['chi2']:.3f} / dof {ref16[w]['dof']} = {ref16[w]['chi2'] / ref16[w]['dof']:.3f}" for w, _ in FITS) + ".", ""]
    R += ["## N17 이 드는 고리", ""]
    R += md(["N17 집합", "고리", "구성 엣지(부호)", "닫힘 값", "전파 불확실도", "닫힘의 절댓값 ÷ 전파"], [[r[0], r[1], r[3], s3(float(r[4])), u3(float(r[5])), u3(float(r[6]))] for r in t4])
    R += ["## N17 엣지의 잔차 (관측 − 적합)", ""]
    R += md(["N17 집합", "구성", "적합", "관측", "적합", "잔차", "σ̂", "잔차/σ̂"], [[r[0], r[1], r[2], s3(float(r[5])), s3(float(r[6])), s3(float(r[7])), u3(float(r[8])), s3(float(r[9]))] for r in t3 if r[3] == "N17"])
    R += ["전 엣지 잔차 `T3_edge_residuals.tsv` · 후보 쌍 21개의 차와 σ합 `T7_candidate_pairs.tsv` · 9/21 값과 나란히 `T6_compare_with_0921.tsv`.", ""]
    with open(out / "REPORT.md", "x", encoding="utf-8") as fh: fh.write("\n".join(R) + "\n")
    return res_all, ranks


def finish(out, extra=None):
    W_(out, "CHECKS.tsv", ["대조", "결과", "a", "b", "허용오차"], CHECKS)
    nb = sum(1 for c in CHECKS if c[1] != "일치")
    meta = dict(at=subprocess.run(["date", "-Is"], capture_output=True, text=True).stdout.strip(), n_checks=len(CHECKS), n_bad=nb,
                script="qsar/src/net17_resolve_v1.py", script_sha256=sha(__file__), final_edges=rel(FE_PATH), final_edges_sha256=sha(FE_PATH),
                numpy=np.__version__, **(extra or {}))
    with open(out / "META.json", "x", encoding="utf-8") as fh: json.dump(meta, fh, ensure_ascii=False, indent=1)
    print(json.dumps(meta, ensure_ascii=False, indent=1))
    return nb


def main():
    a = sys.argv[1:]
    if not a or a[0] not in ("--repro16", "--dry", "--run"): print(__doc__); sys.exit(2)
    if a[0] == "--repro16":
        repro16(); nb = sum(1 for c in CHECKS if c[1] != "일치")
        for c in CHECKS:
            if c[1] != "일치": print("★불일치", c[0], c[2], c[3])
        res = dict(n_checks=len(CHECKS), n_bad=nb, final_edges=rel(FE_PATH), final_edges_sha256=sha(FE_PATH), script_sha256=sha(__file__))
        if "--json" in a: Path(a[a.index("--json") + 1]).write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
        print(json.dumps(res, ensure_ascii=False)); sys.exit(1 if nb else 0)
    out = Path(a[1]).resolve()
    if a[0] == "--dry":
        # 합성 자료: 임의의 노드 값 → 무잡음 엣지 값. 실제 final_edges·N17 값을 읽지 않는다.
        out.mkdir(parents=True, exist_ok=False)
        rng = np.random.default_rng(20261004)
        nodes = sorted({x for e in ALL17 for x in EDGE[e]}); truth = {n: (0.0 if n == REF else float(rng.normal(0, 3))) for n in nodes}
        val = {e: truth[EDGE[e][1]] - truth[EDGE[e][0]] for e in N16}; sdv = {e: float(rng.uniform(0.3, 4.0)) for e in N16}
        nv = {e: 5 for e in N16}; t17 = truth[EDGE["N17"][1]] - truth[EDGE["N17"][0]]
        res_all, _ = compute(out, val, sdv, nv, [("합성-가", t17, 0.9, 3), ("합성-나", t17, 0.8, 5)], "★합성 자료 시험 (실제 값 아님)")
        mx_g = max(abs(r["G"][n] - truth[n]) for r in res_all.values() for n in r["nodes"]); mx_c2 = max(r["chi2"] for r in res_all.values())
        v = dict(val); v["N17"] = t17; s = dict(sdv); s["N17"] = 1.0
        SV = J(FR / "sign_verification.json")["derived"]
        mx_cl = max(abs(closure(c, v, s)[0]) for c in list(SV.values()) + list(NEWCYC.values()))
        chk("합성:노드 값 복원 오차 최대 ≤ 1e-9", mx_g <= 1e-9, True); chk("합성:χ² 최대 ≤ 1e-12", mx_c2 <= 1e-12, True)
        chk("합성:고리 10개 닫힘 최대 ≤ 1e-9", mx_cl <= 1e-9, True)
        # 음성 대조(합성): 엣지 하나의 부호를 뒤집으면 복원이 깨져야 한다
        bad = dict(val); bad["N16"] = -bad["N16"] + 1.0
        rb = solve(ALL17, {**bad, "N17": t17}, {**sdv, "N17": 1.0}, "OLS")
        chk("합성 음성 대조:부호를 뒤집은 엣지가 있으면 χ² > 1e-6", rb["chi2"] > 1e-6, True)
        nb = finish(out, dict(mode="dry", max_node_err=mx_g, max_chi2=mx_c2, max_closure=mx_cl, negctl_chi2=rb["chi2"])); sys.exit(1 if nb else 0)
    # --run : 등재 확인
    pr = J(PREREG); last = pr["changelog"][-1]
    if last.get("키") != REGKEY:
        print(f"★거부: prereg 마지막 항목의 키가 {REGKEY} 가 아니다 ({last.get('키')}) — 등재 전에는 본 실행을 하지 않는다"); sys.exit(2)
    reg_sha = last["(h)★스크립트"]["sha256"]
    if reg_sha != sha(__file__):
        print(f"★거부: 등재된 스크립트 sha256 과 다르다 (등재 {reg_sha[:16]}… · 현재 {sha(__file__)[:16]}…)"); sys.exit(2)
    ref16, val, sdv = repro16()
    if any(c[1] != "일치" for c in CHECKS):
        print("★재현 대조 불일치 — 본 실행을 하지 않는다"); [print("  ", c[0], c[2], c[3]) for c in CHECKS if c[1] != "일치"]; sys.exit(1)
    FE = J(FE_PATH); nv = {e: FE[e]["n"] for e in N16}
    v31 = next(e for e in pr["changelog"] if e.get("version") == 31)["(b)★최종 값"]
    sets = []
    for nm, key in (("주{1,4,5}", "주 집합 {1,4,5}"), ("보조{1..5}", "보조 집합 {1,2,3,4,5}")):
        x = v31[key]; sets.append((nm, float(x["평균"]), float(x["s"]), int(x["n"])))
    reg_in = last["(a)★입력"]["N17 (prereg v31 (b) 등재값 그대로)"]
    for nm, m, s_, n_ in sets:
        chk(f"등재:{nm} 평균 = v33 (a)", m, reg_in[nm]["평균"]); chk(f"등재:{nm} s = v33 (a)", s_, reg_in[nm]["s"]); chk(f"등재:{nm} n = v33 (a)", n_, reg_in[nm]["n"])
    chk("등재:final_edges.json sha256 = v33 (a)", sha(FE_PATH), last["(a)★입력"]["N01~N16"]["sha256"])
    out.mkdir(parents=True, exist_ok=False)
    compute(out, val, sdv, nv, sets, "1단계 — N17 을 넣고 그물 다시 풀기 (prereg v33)", ref16=ref16)
    nb = finish(out, dict(mode="run", regkey=REGKEY)); sys.exit(1 if nb else 0)


if __name__ == "__main__":
    main()
