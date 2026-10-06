#!/usr/bin/env python
"""4단계 분석 — 17회. 역방향은 부호 반전. 실효 g·병목 창 포함."""
import os, sys, json, itertools
import numpy as np, pandas as pd
sys.path.insert(0,'/home/nudge/Project/TYK2-PILOT/07_step3_lambda')
from block_g import g_eff_curve
from alchemlyb.parsing.gmx import extract_u_nk, extract_dHdl
from alchemlyb.estimators import MBAR, BAR, TI
S='/home/nudge/Project/TYK2-PILOT/07_step3_lambda'
KT=8.314462618e-3*300; KC=4.184; NW={1:5,2:9,3:5}; CUT=0.2
RUNS={'edge1':[1,2,3],'edge1R':[1,2,3],'edge2':[1,2,3],'edge2R':[1,2,3],
      'edge3':[1],'edge3R':[1],'edge4':[1],'edge4R':[1],'edge5':[1]}
REV={'edge1R','edge2R','edge3R','edge4R'}

def leg(tag,rep,L,env):
    fs=[f'{S}/s4_{tag}_r{rep}_{env}/L{L}w{i:02d}/prod.xvg' for i in range(NW[L])]
    if not all(map(os.path.exists,fs)): return None
    uk=[];dh=[]
    for f in fs:
        u=extract_u_nk(f,T=300); d=extract_dHdl(f,T=300)
        c=int(len(u)*CUT); uk.append(u.iloc[c:]); dh.append(d.iloc[c:])
    U=pd.concat(uk); D=pd.concat(dh)
    mb=MBAR().fit(U); br=BAR().fit(U); ti=TI().fit(D)
    bse=float(np.sqrt(sum(br.d_delta_f_.iloc[i,i+1]**2 for i in range(NW[L]-1))))*KT/KC
    return dict(MBAR=(float(mb.delta_f_.iloc[0,-1])*KT/KC,float(mb.d_delta_f_.iloc[0,-1])*KT/KC),
                BAR=(float(br.delta_f_.iloc[0,-1])*KT/KC,bse),
                TI=(float(ti.delta_f_.iloc[0,-1])*KT/KC,float(ti.d_delta_f_.iloc[0,-1])*KT/KC))

def gstat(tag,rep):
    """복합체 창별 실효 g → 병목 창"""
    out=[]
    for L in (1,2,3):
        for i in range(NW[L]):
            f=f'{S}/s4_{tag}_r{rep}_cplx/L{L}w{i:02d}/prod.xvg'
            if not os.path.exists(f): continue
            d=extract_dHdl(f,T=300)
            a=np.asarray(d['coul'] if L in (1,3) else (d['vdw']+d['bonded']),float)
            c=g_eff_curve(a)
            out.append((f'L{L}w{i}', float(c[-1][1]) if c else np.nan))
    return out

res={}
for tag,reps in RUNS.items():
    res[tag]={}
    for rep in reps:
        g={}
        ok=True
        for env in ('cplx','solv'):
            legs={}
            for L in (1,2,3):
                e=leg(tag,rep,L,env)
                if e is None: ok=False; break
                legs[L]=e
            if not ok: break
            g[env]=legs
        if not ok: print(f"★{tag} r{rep} 자료 불완전"); continue
        d={}
        for k in ('MBAR','BAR','TI'):
            c=sum(g['cplx'][L][k][0] for L in (1,2,3)); s=sum(g['solv'][L][k][0] for L in (1,2,3))
            ce=np.sqrt(sum(g['cplx'][L][k][1]**2 for L in (1,2,3)))
            se=np.sqrt(sum(g['solv'][L][k][1]**2 for L in (1,2,3)))
            v=c-s; err=float(np.hypot(ce,se))
            if tag in REV: v=-v          # ★역방향: ΔΔG(역) = −MBAR(A→B)
            d[k]=(v,err)
        gs=gstat(tag,rep)
        bn=max(gs,key=lambda x:x[1]) if gs else ('-',np.nan)
        res[tag][rep]=dict(ddg=d, legs={e:{L:g[e][L]['MBAR'][0] for L in (1,2,3)} for e in ('cplx','solv')},
                           g_max=bn[1], g_win=bn[0], g_all=gs)
        print(f"{tag:<7} r{rep}  MBAR {d['MBAR'][0]:+.4f}±{d['MBAR'][1]:.4f}  "
              f"BAR {d['BAR'][0]:+.4f}  TI {d['TI'][0]:+.4f}  g_max {bn[1]:.2f} @{bn[0]}")
json.dump({t:{r:{'ddg':{k:list(v) for k,v in d['ddg'].items()},'g_max':d['g_max'],'g_win':d['g_win'],
                'legs':d['legs'],'g_all':d['g_all']}
              for r,d in rr.items()} for t,rr in res.items()}, open('step4_raw.json','w'), indent=1)
print("\n저장: step4_raw.json")
