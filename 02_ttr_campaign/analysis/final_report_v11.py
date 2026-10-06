#!/usr/bin/env python3
"""★패스 1 최종 정리 — 16엣지 ΔΔG · leg별 σ · 7고리 닫힘 · 연장 분류 (가)(나)(다). ★새 MD 0."""
import subprocess, glob, os, sys, math, json
import statistics as st
os.environ["GMX_MAXBACKUP"]="-1"
G="/home/nudge/miniforge3/envs/md/bin/gmx"; CAN="data/fep_m1/pass1"
W="/tmp/claude-1000/-home-nudge/886b5a77-803b-4e9b-a0ab-abf11f1e6535/scratchpad/bartmp"; os.makedirs(W,exist_ok=True)
RT=0.0019872*310
EDGE={"N01":("CHEMBL241454","CHEMBL240808"),"N02":("CHEMBL240808","cand0023"),
 "N03":("cand0030","CHEMBL240808"),"N04":("CHEMBL240808","CHEMBL438498"),
 "N05":("CHEMBL241454","cand0023"),"N06":("CHEMBL241454","cand0030"),
 "N07":("CHEMBL241454","CHEMBL438498"),"N08":("cand0028","CHEMBL438498"),
 "N09":("cand0030","CHEMBL438498"),"N10":("cand0028","cand0118"),
 "N11":("cand0023","cand0042"),"N12":("cand0030","cand0023"),
 "N13":("cand0118","cand0132"),"N14":("cand0118","cand0271"),
 "N15":("cand0042","cand0132"),"N16":("cand0271","cand0042")}
CYCLES=[("1","N13 N14 N15 N16","cand0118 cand0132 cand0042 cand0271"),
 ("2","N02 N03 N12","CHEMBL240808 cand0030 cand0023"),
 ("3","N05 N06 N12","CHEMBL241454 cand0030 cand0023"),
 ("4","N03 N04 N09","CHEMBL240808 CHEMBL438498 cand0030"),
 ("5","N06 N07 N09","CHEMBL241454 CHEMBL438498 cand0030"),
 ("6","N08 N09 N10 N11 N12 N14 N16","cand0118 cand0028 CHEMBL438498 cand0030 cand0023 cand0042 cand0271"),
 ("7","N01 N02 N05","CHEMBL240808 CHEMBL241454 cand0023")]
TRUE={"N01":(0.000,0.206),"N04":(-3.022,0.307),"N07":(-3.022,0.264)}
def bar(tag):
    xs=sorted(glob.glob(f"{CAN}/{tag}/w*/prod.xvg"))
    if len(xs)<2: return None
    r=subprocess.run([G,"bar","-f",*xs,"-b","0","-o",f"{W}/a.xvg","-oi",f"{W}/i.xvg","-oh",f"{W}/h.xvg"],
                     capture_output=True,text=True,timeout=1800,env={**os.environ})
    t=[l for l in r.stdout.splitlines() if l.startswith("total")]
    return float(t[-1].split("DG")[1].split("+/-")[0])/4.184 if t else None
res={}
for e in EDGE:
    legs={}; dd=[]
    for r in range(1,6):
        v={}
        for s in ("cplx","solv"):
            for l in (1,2,3):
                x=bar(f"{e}_{s}_L{l}_r{r}")
                if x is not None: legs.setdefault((s,l),[]).append(x)
                v[(s,l)]=x
        if all(v[k] is not None for k in v):
            dd.append(sum(v[("cplx",l)] for l in (1,2,3))-sum(v[("solv",l)] for l in (1,2,3)))
    res[e]=dict(n=len(dd), dd=st.mean(dd) if dd else None,
                sd=st.stdev(dd) if len(dd)>1 else None,
                legsd={f"{s}L{l}":(st.stdev(legs[(s,l)]) if len(legs.get((s,l),[]))>1 else None)
                       for s in ("cplx","solv") for l in (1,2,3)})
inc={e:[c for c,eg,_ in CYCLES if e in eg.split()] for e in EDGE}
print("★★패스 1 최종 — 엣지 16개\n")
print(f"{'엣지':5}{'n':>3}{'ΔΔG':>9}{'★σ':>8}{'cL1':>7}{'cL2':>7}{'cL3':>7}{'sL1':>7}{'sL2':>7}{'sL3':>7}  고리")
for e in EDGE:
    R=res[e]; f=lambda x: f"{x:7.3f}" if x is not None else f"{'—':>7}"
    sd=R['sd']; mark="✓완료" if (sd is not None and sd<=0.5) else ("★연장후보" if sd is not None else "")
    print(f"{e:5}{R['n']:>3}{(f'{R[chr(100)+chr(100)]:9.3f}' if R['dd'] is not None else f'{chr(8212):>9}')}"
          f"{(f'{sd:8.3f}' if sd is not None else f'{chr(8212):>8}')}"
          +"".join(f(R['legsd'][k]) for k in ("cplxL1","cplxL2","cplxL3","solvL1","solvL2","solvL3"))
          +f"  {','.join(inc[e]) or '-':10} {mark}")
print("\n★★고리 7개 닫힘 오차")
for c,eg,path in CYCLES:
    nodes=path.split(); tot=0.0; ok=True; parts=[]
    for i in range(len(nodes)):
        a,b=nodes[i],nodes[(i+1)%len(nodes)]
        hit=None
        for e,(x,y) in EDGE.items():
            if (x,y)==(a,b): hit=(e,+1); break
            if (y,x)==(a,b): hit=(e,-1); break
        if hit is None: ok=False; parts.append(f"★없음 {a}→{b}"); continue
        if res[hit[0]]['dd'] is None: ok=False; parts.append(f"{hit[0]}:미완"); continue
        tot+=hit[1]*res[hit[0]]['dd']; parts.append(f"{hit[1]:+d}{hit[0]}")
    biggest=max((eg.split()), key=lambda e: (res[e]['sd'] or -1))
    print(f"  고리{c} ({len(eg.split())}엣지): 닫힘 {tot:+7.3f}  {'' if ok else '★미완'}  "
          f"경로 {' '.join(parts)}  · σ최대 {biggest}({res[biggest]['sd'] if res[biggest]['sd'] is not None else '—'})")
print("\n★★연장 분류 (동결 규칙 17:37)")
for e in EDGE:
    R=res[e]; sd=R['sd']
    if e in TRUE:
        cat="(가) 앵커"; note=("N01=완료·연장없음" if e=="N01" else ("승인 $6" if e=="N04" else ("σ>0.5 → 연장" if (sd or 0)>0.5 else "σ≤0.5 → 안함")))
    elif sd is not None and sd>0.5 and inc[e]:
        cat="(나) 고리"; note=f"고리 {','.join(inc[e])}"
    else:
        cat="(다) 제외"; note="실험값 없음·σ≤0.5 또는 고리 기여 작음"
    if cat!="(다) 제외" or (sd or 0)>0.5:
        print(f"  {e}: {cat}  σ={sd if sd is not None else '—'}  {note}")
json.dump({e:{k:v for k,v in res[e].items()} for e in res}, open("docs/stage1_out/pass1/final_edges.json","w"), indent=1, ensure_ascii=False)
