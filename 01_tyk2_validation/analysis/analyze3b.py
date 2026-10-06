#!/usr/bin/env python
"""3블록 FEP 분석 — 레그별 ΔG 합산 + (a)~(f)"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0,'/home/nudge/Project/TYK2-PILOT/07_step3_lambda')
from block_g import g_eff_curve
from alchemlyb.parsing.gmx import extract_u_nk, extract_dHdl
from alchemlyb.estimators import MBAR, BAR, TI
from pymbar import timeseries
S='/home/nudge/Project/TYK2-PILOT/07_step3_lambda'
KT=8.314462618e-3*300; KC=4.184
NW={1:5,2:9,3:5}
CUT_FRAC=0.2      # 앞 20% 일괄 절단 (창별 t0 는 인공 구멍을 만든다 — edge1 에서 확인)

def leg_estimates(tag,leg,env):
    fs=[f'{S}/3b_{tag}_{env}/L{leg}w{i:02d}/prod.xvg' for i in range(NW[leg])]
    if not all(os.path.exists(f) for f in fs): return None
    uk=[];dh=[]
    for f in fs:
        u=extract_u_nk(f,T=300); d=extract_dHdl(f,T=300)
        c=int(len(u)*CUT_FRAC); uk.append(u.iloc[c:]); dh.append(d.iloc[c:])
    U=pd.concat(uk); D=pd.concat(dh)
    mb=MBAR().fit(U); br=BAR().fit(U); ti=TI().fit(D)
    bse=float(np.sqrt(sum(br.d_delta_f_.iloc[i,i+1]**2 for i in range(NW[leg]-1))))*KT/KC
    O=np.asarray(mb.overlap_matrix)
    return dict(MBAR=(float(mb.delta_f_.iloc[0,-1])*KT/KC,float(mb.d_delta_f_.iloc[0,-1])*KT/KC),
                BAR=(float(br.delta_f_.iloc[0,-1])*KT/KC,bse),
                TI=(float(ti.delta_f_.iloc[0,-1])*KT/KC,float(ti.d_delta_f_.iloc[0,-1])*KT/KC),
                overlap=[float(O[i,i+1]) for i in range(NW[leg]-1)],
                files=fs)

def window_stats(tag,leg,env):
    out=[]
    for i in range(NW[leg]):
        f=f'{S}/3b_{tag}_{env}/L{leg}w{i:02d}/prod.xvg'
        if not os.path.exists(f): continue
        d=extract_dHdl(f,T=300)
        a=np.asarray(d['coul'] if leg in (1,3) else (d['vdw']+d['bonded']),float)
        _,gr,_=timeseries.detect_equilibration(a); gr*=0.2
        c=g_eff_curve(a); ge=c[-1][1]
        out.append(dict(leg=leg,w=i,sd=float(a.std(ddof=1)),g_raw=float(gr),g_eff=float(ge),
                        mean=float(a.mean()),N=len(a)))
    return out

EXP={'edge1':-0.24,'edge2':None,'edge3':None,'edge4':None}
res={}
for tag in ('edge1','edge2','edge3','edge4'):
    got={}
    for env in ('cplx','solv'):
        legs={}
        for L in (1,2,3):
            e=leg_estimates(tag,L,env)
            if e is None: break
            legs[L]=e
        if len(legs)!=3: break
        got[env]=legs
    if len(got)!=2:
        print(f"{tag}: 자료 불완전 — 건너뜀"); continue
    print(f"\n{'='*84}\n### {tag} ###")
    tot={}
    for env in ('cplx','solv'):
        print(f"  --- {env} ---")
        for k in ('MBAR','BAR','TI'):
            v=sum(got[env][L][k][0] for L in (1,2,3))
            e=float(np.sqrt(sum(got[env][L][k][1]**2 for L in (1,2,3))))
            tot[(env,k)]=(v,e)
            per=" + ".join(f"{got[env][L][k][0]:+.3f}" for L in (1,2,3))
            print(f"    {k:<5} {per} = {v:+.4f} ± {e:.4f}")
    print(f"  --- ΔΔG ---")
    dd={}
    for k in ('MBAR','BAR','TI'):
        v=tot[('cplx',k)][0]-tot[('solv',k)][0]
        e=float(np.hypot(tot[('cplx',k)][1],tot[('solv',k)][1]))
        dd[k]=(v,e); print(f"    {k:<5} {v:+.4f} ± {e:.4f} kcal/mol")
    sp=max(dd[k][0] for k in dd)-min(dd[k][0] for k in dd)
    print(f"    세 추정량 폭 {sp:.4f}" + (f" · 실험 {EXP[tag]:+.2f} · MBAR 오차 {dd['MBAR'][0]-EXP[tag]:+.3f}" if EXP[tag] else ""))
    res[tag]=dict(total={f"{e}_{k}":tot[(e,k)] for e in ('cplx','solv') for k in ('MBAR','BAR','TI')},
                  ddg={k:dd[k] for k in dd}, spread=sp,
                  overlap={f"{e}_L{L}":got[e][L]['overlap'] for e in ('cplx','solv') for L in (1,2,3)},
                  win={f"{e}_L{L}":window_stats(tag,L,e) for e in ('cplx','solv') for L in (1,2,3)})
json.dump(res,open(f'{S}/analysis3b.json','w'),indent=1,default=float)
print(f"\n저장: analysis3b.json")
