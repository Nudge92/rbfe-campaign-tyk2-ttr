#!/usr/bin/env python
"""순/역방향 수렴 — 1 k_BT 게이트 (K-3 / L-3).

기준 (Mey et al. 2020, LiveCoMS 2, 18378 · §12 체크리스트):
  "disagreements larger than ★1 k_BT in the final part of the forward and reverse
   trajectories can be useful to detect unconverged results"

alchemlyb.convergence.forward_backward_convergence 를 쓴다(Klimovich 2015 절차).
반환 단위는 ★k_BT 다 — 변환하지 않는다. 그래야 문턱과 같은 단위다.

★최종 구간 = data_fraction 1.0 지점. 순방향은 전체 앞에서부터 100 %,
역방향은 뒤에서부터 100 % 이므로 ★둘 다 전체 자료를 쓴다. 따라서 그 지점의
차이는 원리상 0 에 가까워야 하고, 의미 있는 비교는 ★마지막에서 두 번째 지점
(앞/뒤 각각 90 %)이다. 두 값을 모두 보고한다.

동결 산출물 읽기 전용 · 새 MD 0 · GPU 불필요. 기존 코드 무수정.
"""
import glob
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
TEMP = 310.0
KCAL_PER_KT = 0.0019872041 * TEMP     # 310 K → 0.61603 kcal/mol


def leg_conv(leg_glob, num=10):
    from alchemlyb.parsing.gmx import extract_u_nk
    from alchemlyb.convergence import forward_backward_convergence
    fs = sorted(glob.glob(leg_glob + "/prod.xvg"))
    if not fs:
        return None, 0
    dfs = [extract_u_nk(f, TEMP) for f in fs]
    return forward_backward_convergence(dfs, "MBAR", num=num), len(fs)


def report(label, root, legs=("L1", "L2", "L3"), num=10):
    print(f"\n{'='*74}\n## {label}\n   {root}")
    for L in legs:
        try:
            c, n = leg_conv(f"{root}/{L}w*", num)
        except Exception as e:
            print(f"\n### {L}: ★실패 {type(e).__name__}: {e}")
            continue
        if c is None:
            print(f"\n### {L}: xvg 없음")
            continue
        f = np.asarray(c["Forward"], float)
        b = np.asarray(c["Backward"], float)
        fe = np.asarray(c["Forward_Error"], float)
        be = np.asarray(c["Backward_Error"], float)
        frac = np.asarray(c["data_fraction"], float)
        d = np.abs(f - b)
        print(f"\n### {L}  창 {n}")
        print(f"{'분율':>6} {'순방향':>11} {'역방향':>11} {'|차|':>9} {'순오차':>8} {'역오차':>8}")
        for i in range(len(f)):
            mark = "  ★>1kT" if d[i] > 1.0 else ""
            print(f"{frac[i]:>6.1f} {f[i]:>11.4f} {b[i]:>11.4f} {d[i]:>9.4f} "
                  f"{fe[i]:>8.4f} {be[i]:>8.4f}{mark}")
        last, prev = d[-1], d[-2] if len(d) > 1 else np.nan
        print(f"  ★최종 구간(분율 1.0) |차| = {last:.4f} k_BT = {last*KCAL_PER_KT:.4f} kcal/mol"
              f"   → {'★1 k_BT 초과' if last > 1.0 else '1 k_BT 이내'}")
        print(f"   직전 구간(분율 {frac[-2]:.1f}) |차| = {prev:.4f} k_BT = {prev*KCAL_PER_KT:.4f} kcal/mol"
              f"   → {'★1 k_BT 초과' if prev > 1.0 else '1 k_BT 이내'}")
        print(f"   |차| 최댓값(전 구간) = {d.max():.4f} k_BT at 분율 {frac[int(np.argmax(d))]:.1f}")


if __name__ == "__main__":
    B = "/home/nudge/Project/CADD/ttr/qsar/data/fep_m1"
    targets = [
        ("v12d 회차1 복합체 (L3 7창 · ★2 ns)",  f"{B}/prod_v12d_r1/N17/cplx/r1"),
        ("v12d 회차1 용매   (L3 7창 · ★2 ns)",  f"{B}/prod_v12d_r1/N17/solv/r1"),
        ("v12f 회차1 복합체 (L3 10창 · ★2 ns)", f"{B}/prod_v12f_r1/N17/cplx/r1"),
        ("v12f 회차1 용매   (L3 10창 · ★2 ns)", f"{B}/prod_v12f_r1/N17/solv/r1"),
        ("v14  회차1 복합체 (L3 10창 · ★8 ns)", f"{B}/prod_v14_r1/N17/cplx/r1"),
    ]
    only = sys.argv[1] if len(sys.argv) > 1 else None
    print(f"# 1 k_BT = {KCAL_PER_KT:.5f} kcal/mol  (T = {TEMP} K)")
    for lab, root in targets:
        if only and only not in root:
            continue
        report(lab, root)
