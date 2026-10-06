#!/usr/bin/env python3
"""(A) 고리 부호 독립 재유도 — ★새 계산 0. 원장(hyb.log · merged.itp)만 읽는다.
★경로 1: merged.itp 의 상태 A/B 원자형(할로겐·더미)
★경로 2: pmx hyb.log 의 'Reading ligand 1/2 from' (= 상태 A/B)
★두 경로가 엇갈리면 ★"확인 불가".
"""
import re, json
from pathlib import Path
Q=Path("/home/nudge/Project/CADD/ttr/qsar"); OUT=Q/"docs/stage1_out/final_report_20260921"
EDGE={"N01":("CHEMBL241454","CHEMBL240808"),"N02":("CHEMBL240808","cand0023"),
 "N03":("cand0030","CHEMBL240808"),"N04":("CHEMBL240808","CHEMBL438498"),
 "N05":("CHEMBL241454","cand0023"),"N06":("CHEMBL241454","cand0030"),
 "N07":("CHEMBL241454","CHEMBL438498"),"N08":("cand0028","CHEMBL438498"),
 "N09":("cand0030","CHEMBL438498"),"N10":("cand0028","cand0118"),
 "N11":("cand0023","cand0042"),"N12":("cand0030","cand0023"),
 "N13":("cand0118","cand0132"),"N14":("cand0118","cand0271"),
 "N15":("cand0042","cand0132"),"N16":("cand0271","cand0042")}
CYCLES=[("1","cand0118 cand0132 cand0042 cand0271"),
 ("2","CHEMBL240808 cand0030 cand0023"),("3","CHEMBL241454 cand0030 cand0023"),
 ("4","CHEMBL240808 CHEMBL438498 cand0030"),("5","CHEMBL241454 CHEMBL438498 cand0030"),
 ("6","cand0118 cand0028 CHEMBL438498 cand0030 cand0023 cand0042 cand0271"),
 ("7","CHEMBL240808 CHEMBL241454 cand0023")]
# ★할로겐 수(구조에서 확정 · TARGET_DATA_STATUS.md)
HAL={"CHEMBL241454":0,"CHEMBL240808":2,"CHEMBL438498":2,"cand0023":0,"cand0028":1,
     "cand0030":0,"cand0042":0,"cand0118":1,"cand0132":0,"cand0271":0}

def path2(N):
    """pmx hyb.log 의 ligand 1/2 = 상태 A/B"""
    f=Q/f"data/fep_m1/junction_v9/{N}/hyb.log"
    if not f.exists(): return None
    t=f.read_text(errors="replace")
    m1=re.search(r'Reading ligand 1 from: "([^"]+)"', t)
    m2=re.search(r'Reading ligand 2 from: "([^"]+)"', t)
    if not (m1 and m2): return None
    g=lambda p: Path(p).parent.name
    return (g(m1.group(1)), g(m2.group(1)))

def path1(N):
    """merged.itp 의 할로겐 상태로 A/B 판별 (할로겐 수가 다른 쌍에서만 결정적)"""
    f=Q/f"data/fep_m1/junction_v9/{N}/merged.itp"
    if not f.exists(): return None
    inb=False; nA=nB=0
    for ln in f.read_text().splitlines():
        s=ln.strip()
        if s.startswith("["): inb=s.replace(" ","").startswith("[atoms]"); continue
        if not inb or not s or s.startswith(";"): continue
        fl=s.split()
        if len(fl)<8: continue
        tA=fl[1]; tB=fl[8] if len(fl)>=11 else fl[1]
        if tA in ("cl","br","i","f"): nA+=1
        if tB in ("cl","br","i","f"): nB+=1
    a,b=EDGE[N]
    if HAL[a]==HAL[b]: return ("동수-판별불가", nA, nB)
    # 관측 (nA,nB) 가 (HAL[a],HAL[b]) 와 맞으면 정방향
    if (nA,nB)==(HAL[a],HAL[b]): return ("정방향", nA, nB)
    if (nA,nB)==(HAL[b],HAL[a]): return ("★역방향", nA, nB)
    return ("★불일치", nA, nB)

print("════ ③ 간선 방향 — 2경로 독립 확인 ════")
print(f"| 간선 | 원장 정의 (from→to) | 경로1 할로겐(A,B) | 경로2 pmx ligand1→2 | 판정 |")
print("|---|---|---|---|---|")
verdict={}
for N in sorted(EDGE):
    a,b=EDGE[N]; p1=path1(N); p2=path2(N)
    ok2 = (p2==(a,b)) if p2 else None
    s1 = f"{p1[0]} ({p1[1]},{p1[2]})" if p1 else "—"
    s2 = f"{p2[0]}→{p2[1]}" if p2 else "—"
    if p2 is None: v="★확인 불가(로그 없음)"
    elif not ok2: v="★★역방향"
    elif p1 and p1[0]=="★역방향": v="★★두 경로 엇갈림 — 확인 불가"
    elif p1 and p1[0]=="정방향": v="**일치(2경로)**"
    else: v="일치(경로2) · 경로1 판별불가"
    verdict[N]=v
    print(f"| **{N}** | {a} → {b} | {s1} | {s2} | {v} |")

print("\n════ ④ 7고리 부호 재유도 ════")
print("| 고리 | 노드 순회 | 재유도 부호 | 현재 코드 부호 | 일치 |")
print("|---|---|---|---|---|")
CUR={'1':{'N13':1,'N15':-1,'N16':-1,'N14':-1},'2':{'N03':-1,'N12':1,'N02':-1},'3':{'N06':1,'N12':1,'N05':-1},
     '4':{'N04':1,'N09':-1,'N03':1},'5':{'N07':1,'N09':-1,'N06':-1},
     '6':{'N10':-1,'N08':1,'N09':-1,'N12':1,'N11':1,'N16':-1,'N14':-1},'7':{'N01':-1,'N05':1,'N02':-1}}
allok=True
derived={}
for c,path in CYCLES:
    nodes=path.split(); sg={}; parts=[]; bad=False
    for i in range(len(nodes)):
        x,y=nodes[i], nodes[(i+1)%len(nodes)]
        hit=None
        for e,(p,q) in EDGE.items():
            if (p,q)==(x,y): hit=(e,+1); break
            if (q,p)==(x,y): hit=(e,-1); break
        if hit is None: parts.append(f"★없음 {x}→{y}"); bad=True; continue
        sg[hit[0]]=hit[1]; parts.append(f"{'+' if hit[1]>0 else '−'}{hit[0]}")
    derived[c]=sg
    same = (sg==CUR[c])
    allok &= same
    print(f"| 고리{c} | {' → '.join(n[-7:] for n in nodes)} → {nodes[0][-7:]} | {' '.join(parts)} | "
          f"{' '.join(('+' if v>0 else '−')+k for k,v in CUR[c].items())} | {'★일치' if same else '★★불일치'} |")
print(f"\n★★7고리 부호 재유도 = {'전부 일치' if allok else '★불일치 있음'}")
json.dump(dict(edge_verdict=verdict, derived={c:{k:int(v) for k,v in s.items()} for c,s in derived.items()},
               all_match=bool(allok)), open(OUT/"sign_verification.json","w"), indent=1, ensure_ascii=False)
