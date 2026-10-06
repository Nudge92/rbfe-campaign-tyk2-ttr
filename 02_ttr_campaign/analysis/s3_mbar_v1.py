#!/usr/bin/env python3
"""S-3 (a)(b)(c) — 일관성 대조 · 겹침 행렬 · 정방향/역방향. ★새 MD 0 · 기존 prod.xvg 재적합.
입력은 n17_uq_v1.py 와 같다: 창별 extract_u_nk(310 K) 를 ★절단 없이·솎음 없이 이어 붙인다.
(a) 전체 적합 ΔG 를 보고된 레그 ΔG(qsar/reports/uq1_<pre>_<edge>.json)와 대조.
(b) 겹침: 등재 계산(c_criteria_v13.py:27-28)과 같은 정의 — adj = O[i,i+1], C-1 = min(adj) ≥ 0.03.
    삼중대각 확인용으로 아래 대각 O[i+1,i] 도 함께 낸다.
(c) alchemlyb.convergence.forward_backward_convergence(us, 'MBAR', num=4) — kT → kcal/mol 은
    같은 적합의 to_kcalmol 비로 환산(같은 상수).
사용: s3_mbar_v1.py <prod접두> <edge> <rep...> [--negctl]
"""
import sys, glob, json, warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd

T = Path("/home/nudge/Project/CADD/ttr"); TEMP = 310.0; NUM = 4
LEGS = [(s, l) for s in ("cplx", "solv") for l in (1, 2, 3)]
COORD = {1: "coul-lambdas", 2: "vdw-lambdas", 3: "coul-lambdas"}


def lam(mdp, L):
    for ln in Path(mdp).read_text().splitlines():
        if ln.lower().startswith(COORD[L]):
            return [float(x) for x in ln.split("=", 1)[1].split()]


def fit_leg(pre, edge, rep, s, l, drop=None, conv=True):
    from alchemlyb.parsing.gmx import extract_u_nk
    from alchemlyb.estimators import MBAR
    from alchemlyb.postprocessors.units import to_kcalmol
    D = T/f"qsar/data/fep_m1/{pre}_r{rep}/{edge}/{s}/r{rep}"
    fs = sorted(glob.glob(str(D/f"L{l}w*/prod.xvg")))
    dropped = None
    if drop is not None:
        dropped = extract_u_nk(fs[drop], TEMP)             # 빠질 창의 상태 라벨을 얻는다
        fs = [f for i, f in enumerate(fs) if i != drop]
    us = [extract_u_nk(f, TEMP) for f in fs]
    if dropped is not None:                                 # ★표본 0 인 상태 열도 함께 뺀다(안 빼면 alchemlyb 가 예외)
        key = dropped.index[0][1:]; key = key[0] if len(key) == 1 else key
        us = [x.drop(columns=[key]) for x in us]
    u = pd.concat(us)
    m = MBAR().fit(u)
    dkt = float(m.delta_f_.iloc[0, -1]); dG = float(to_kcalmol(m.delta_f_).iloc[0, -1])
    fac = dG/dkt                                             # kT → kcal/mol (alchemlyb 상수 그대로)
    out = dict(nwin=len(fs), nrow=[len(x) for x in us], dG=dG, kt2kcal=fac)
    if drop is None:
        O = np.asarray(m.overlap_matrix, float); K = O.shape[0]
        lv = lam(D/f"L{l}w00"/"prod.mdp", l)
        sup = [float(O[i, i+1]) for i in range(K-1)]; sub = [float(O[i+1, i]) for i in range(K-1)]
        i = int(np.argmin(sup))
        out.update(K=K, lam=lv, ov_sup=sup, ov_sub=sub, ov_diag=[float(O[i, i]) for i in range(K)],
                   ov_min=float(min(sup)), ov_min_pair=[lv[i], lv[i+1]],
                   ov_min_both=float(min(sup+sub)),
                   n_adj_exact_zero=int(sum(1 for x in sup+sub if x == 0.0)),
                   n_adj_lt_0p005=int(sum(1 for x in sup+sub if x < 0.005)),
                   row_sum_maxdev=float(np.abs(O.sum(1)-1).max()),
                   c1="통과" if min(sup) >= 0.03 else "★미달")
    if conv and drop is None:
        from alchemlyb.convergence import forward_backward_convergence
        fb = forward_backward_convergence(us, "MBAR", num=NUM)
        out["fb"] = dict(frac=[float(x) for x in fb["data_fraction"]],
                         fwd=[float(x)*fac for x in fb["Forward"]],
                         fwd_err=[float(x)*fac for x in fb["Forward_Error"]],
                         bwd=[float(x)*fac for x in fb["Backward"]],
                         bwd_err=[float(x)*fac for x in fb["Backward_Error"]])
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    pre, edge, reps = args[0], args[1], args[2:]
    REF = json.load(open(T/f"qsar/reports/uq1_{pre}_{edge}.json"))
    if "--negctl" in sys.argv:                               # ★음성 대조: 가운데 창 하나를 뺀다
        rep, s, l = reps[0], "cplx", 3
        full = REF[f"r{rep}"][f"{s}L{l}"]["dG"]
        r = fit_leg(pre, edge, rep, s, l, drop=4, conv=False)
        d = abs(r["dG"]-full)
        print(f"  [음성 대조] r{rep} {s}L{l} 에서 창 1개(정렬 순서 5번째)를 뺀 입력: ΔG {r['dG']:+.6f} vs 보고 {full:+.6f} "
              f"· |차| {d:.6f} → {'✓불일치 검출' if d > 1e-6 else '✗검출 실패'} (창 {r['nwin']})")
        return
    out = {}; worst = 0.0
    for rep in reps:
        out[f"r{rep}"] = {}
        for s, l in LEGS:
            key = f"{s}L{l}"
            r = fit_leg(pre, edge, rep, s, l)
            ref = REF[f"r{rep}"][key]["dG"]; d = abs(r["dG"]-ref); worst = max(worst, d)
            r["dG_reported"] = ref; r["dG_absdiff"] = d
            out[f"r{rep}"][key] = r
            fb = r["fb"]
            print(f"  r{rep} {key}: ΔG {r['dG']:+.6f} (보고 {ref:+.6f} · |차| {d:.1e}) · 겹침최소 {r['ov_min']:.4f} "
                  f"@λ {r['ov_min_pair'][0]:.4f}→{r['ov_min_pair'][1]:.4f} · {r['c1']} · "
                  f"정 {' '.join(f'{x:+.3f}' for x in fb['fwd'])} · 역 {' '.join(f'{x:+.3f}' for x in fb['bwd'])}", flush=True)
            if d > 1e-6:
                print(f"★★중단: r{rep} {key} 전체 적합 ΔG 가 보고값과 다르다 (|차| {d:.3e}) — 입력이 다르다"); 
                json.dump(out, open(T/f"qsar/reports/s3_mbar_{pre}_{edge}.PARTIAL.json", "w"), indent=1, default=float)
                sys.exit(3)
    p = T/f"qsar/reports/s3_mbar_{pre}_{edge}.json"
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float), encoding="utf-8")
    print(f"\n  일관성 대조: 30 레그 최대 |차| {worst:.3e} kcal/mol\n저장 {p}")


if __name__ == "__main__":
    main()
