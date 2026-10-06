#!/usr/bin/env python
"""
3단계 edge1 분석 (b)~(f).

★창별 dH/dλ 의 '활성 성분'만 쓴다 — CDK2 감사에서 배운 것:
  vdw 블록(0~12)에서는 λ_coul 이 고정이라 dH/dλ_coul 은 적분에 기여하지 않는다.
  성분을 합쳐 SD 를 재면 기여하지 않는 잡음까지 포함되어 불확도를 부풀린다.
  vdw 블록 활성 = vdw + bonded · coul 블록 활성 = coul
★TI 가중치는 성분별 사다리꼴이며 블록당 Σw = 1 이다.
"""
import json, os, sys
import numpy as np
import pandas as pd
from alchemlyb.parsing.gmx import extract_u_nk, extract_dHdl
from alchemlyb.estimators import MBAR, BAR, TI
from pymbar import timeseries

S = os.path.dirname(os.path.abspath(__file__))
T = 300.0
KT_KJ = 8.314462618e-3*T
KCAL = 4.184
NW = 24
JUNC = 12                      # 접합부 상태 (vdw=1, coul=0)

def lam_vectors():
    fb = open(f'{S}/mdp/fep_block.txt').read()
    g = lambda k: [float(x) for x in
                   [l for l in fb.splitlines() if l.startswith(k)][0].split('=')[1].split()]
    return g('coul-lambdas'), g('vdw-lambdas')

def trap_w(l):
    n = len(l); w = []
    for i in range(n):
        lo = l[i]-l[i-1] if i > 0 else 0.0
        hi = l[i+1]-l[i] if i < n-1 else 0.0
        w.append((lo+hi)/2)
    return w

def analyze(leg):
    files = [f'{S}/{leg}/w{ i:02d}/prod.xvg'.replace(' ', '') for i in range(NW)]
    files = [f for f in files if os.path.exists(f)]
    if len(files) < NW:
        print(f"  ★{leg}: {len(files)}/{NW} 창만 존재 — 부분 분석"); 
    uk = [extract_u_nk(f, T=T) for f in files]
    dh = [extract_dHdl(f, T=T) for f in files]
    coul, vdw = lam_vectors()
    wc, wv = trap_w(coul), trap_w(vdw)

    # ---------- (d) 창별 t0 / g / Neff (활성 성분 기준) ----------
    rows = []
    for i, d in enumerate(dh):
        act = (d['vdw']+d['bonded']) if i < JUNC else d['coul']
        act = np.asarray(act, float)
        t = d.index.get_level_values('time').values
        dt = t[1]-t[0]
        t0, g, Neff = timeseries.detect_equilibration(act)
        w = wv[i] if i < JUNC else wc[i]
        sd = act[int(t0):].std(ddof=1)
        se = sd/np.sqrt(max(Neff, 1))
        rows.append(dict(state=i, block='vdw' if i < JUNC else 'coul',
                         coul=coul[i], vdw=vdw[i], N=len(act), dt=float(dt),
                         t0_ps=float(t0*dt), g_ps=float(g*dt), Neff=float(Neff),
                         mean=float(act.mean()), sd=float(sd), se=float(se),
                         w=float(w), u_kcal=float(se*w/KCAL)))
    # ---------- 평형 구간 제거 후 estimator ----------
    uk_eq = []
    for i, u in enumerate(uk):
        t0 = int(rows[i]['t0_ps']/rows[i]['dt'])
        uk_eq.append(u.iloc[t0:])
    U = pd.concat(uk_eq)
    dh_eq = [d.iloc[int(rows[i]['t0_ps']/rows[i]['dt'])] if False else
             d.iloc[int(rows[i]['t0_ps']/rows[i]['dt']):] for i, d in enumerate(dh)]
    D = pd.concat(dh_eq)

    out = {}
    mb = MBAR().fit(U)
    out['MBAR'] = (float(mb.delta_f_.iloc[0, -1])*KT_KJ/KCAL,
                   float(mb.d_delta_f_.iloc[0, -1])*KT_KJ/KCAL)
    out['overlap'] = np.asarray(mb.overlap_matrix)
    out['mbar_perwin'] = [(float(mb.delta_f_.iloc[i, i+1])*KT_KJ/KCAL,
                           float(mb.d_delta_f_.iloc[i, i+1])*KT_KJ/KCAL)
                          for i in range(len(files)-1)]
    br = BAR().fit(U)
    out['BAR'] = (float(br.delta_f_.iloc[0, -1])*KT_KJ/KCAL,
                  float(br.d_delta_f_.iloc[0, -1])*KT_KJ/KCAL)
    ti = TI().fit(D)
    out['TI'] = (float(ti.delta_f_.iloc[0, -1])*KT_KJ/KCAL,
                 float(ti.d_delta_f_.iloc[0, -1])*KT_KJ/KCAL)
    out['rows'] = rows
    return out

res = {}
for leg in ('cplx', 'solv'):
    print(f"\n{'='*80}\n### {leg} ###")
    r = analyze(leg); res[leg] = r
    print(f"  MBAR {r['MBAR'][0]:+.3f} ± {r['MBAR'][1]:.3f} kcal/mol")
    print(f"  BAR  {r['BAR'][0]:+.3f} ± {r['BAR'][1]:.3f}")
    print(f"  TI   {r['TI'][0]:+.3f} ± {r['TI'][1]:.3f}")

# ---------- ΔΔG ----------
print(f"\n{'='*80}\n### ΔΔG = complex − solvent ###")
for k in ('MBAR', 'BAR', 'TI'):
    v = res['cplx'][k][0]-res['solv'][k][0]
    e = np.hypot(res['cplx'][k][1], res['solv'][k][1])
    print(f"  {k:<5} {v:+.3f} ± {e:.3f} kcal/mol")
vals = [res['cplx'][k][0]-res['solv'][k][0] for k in ('MBAR', 'BAR', 'TI')]
print(f"  세 추정량 폭 = {max(vals)-min(vals):.3f} kcal/mol")
print(f"  실험값 ΔΔG(ejm_31→ejm_42) = -0.24 kcal/mol")

json.dump({leg: {k: (v if not isinstance(v, np.ndarray) else v.tolist())
                 for k, v in r.items()} for leg, r in res.items()},
          open(f'{S}/analysis_fep.json', 'w'), indent=1, default=float)
np.save(f'{S}/overlap_cplx.npy', res['cplx']['overlap'])
np.save(f'{S}/overlap_solv.npy', res['solv']['overlap'])
print(f"\n저장: analysis_fep.json · overlap_*.npy")
