#!/usr/bin/env python3
"""★T 절 보고 표 생성 — ★읽기 전용. ★final_edges.json 만 근거로 한다(재계산 0).
   규격은 docs/TTR_FEP_CONTROL.md 의 ★T 절. ★열 추가·이름 변경 금지."""
import json, math, time, os, sys, glob
from pathlib import Path

R  = Path("/home/nudge/Project/CADD/ttr/qsar")
FE = R/"docs/stage1_out/pass1/final_edges.json"
AS = R/"docs/stage1_out/pass1/assignment.json"

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
TRUE={"N01":0.000,"N04":-3.022,"N07":-3.022}
# ★열 10 — 2026-09-20 파일 확인분 (docs/TTR_FEP_CONTROL.md T-1)
NOTE={"N01":"n=3 · ★solv r1·r2 가 ★어떤 로컬 트리에도 없음(pass1·pass1_full·pass_ext·pass1_2ns 전수 0건)",
      "N11":"★09-21 재실행 회수 완료(5/5 · 창 21 · OK 마커) → n=5 · ★σ̂ 4.632 는 그물 최대",
      "N08":"n=4 · cplx L3 r2 없음",
      "N14":"n=4 · solv L3 r5 없음"}
CLOSE_TOL = 1.0   # ★닫힘 기준 — F-6 널 문턱과 같은 값을 쓴다(기준 먼저 적고 판정)

def wins():
    """★엣지별 창 수 — assignment.json 의 win (lists + unresolved 둘 다)"""
    A=json.load(open(AS)); w={}
    recs=[j for pl in A["lists"].values() for j in pl]+A.get("unresolved",[])
    for j in recs:
        e=j["tag"].split("_")[0]; s="cplx" if "_cplx_" in j["tag"] else "solv"; l=int(j["leg"]) if "leg" in j else int(j["tag"].split("_L")[1][0])
        w.setdefault(e,{}).setdefault(s,{})[l]=j["win"]
    return {e:{s:sum(v.get(s,{}).values()) for s in ("cplx","solv")} for e,v in w.items()}

def verdict(n, sd):
    if n is None or n < 3: return "판정불가(표본부족)"
    if sd is None:         return "판정불가(σ없음)"
    return "마감가능" if sd <= 0.5 else "연장필요"

def main():
    if not FE.exists(): print("★final_edges.json 없음 — S1 미완"); return 2
    d=json.load(open(FE)); W=wins()
    mt=time.strftime("%F %H:%M:%S", time.localtime(FE.stat().st_mtime))
    print(f"근거 {FE.relative_to(R)} · 생성 {mt} · 기준 길이 0.5 ns (패스 1)\n")
    print("## T-1 · 엣지 표")
    print("| 엣지 | 노드 쌍 | 창 수 cplx/solv(합) | n | ΔΔG | σ (n) | 실험값 | 오차 | A1 판정 | 문제 |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    cnt={"마감가능":0,"연장필요":0,"판정불가":0}
    for e in [f"N{i:02d}" for i in range(1,17)]:
        r=d.get(e,{}); n=r.get("n",0); dd=r.get("dd"); sd=r.get("sd")
        c=W.get(e,{}).get("cplx","?"); s=W.get(e,{}).get("solv","?")
        tot=(c+s) if isinstance(c,int) and isinstance(s,int) else "?"
        v=verdict(n,sd); cnt["판정불가" if v.startswith("판정불가") else v]+=1
        exp=f"{TRUE[e]:+.3f}" if e in TRUE else "—"
        err=f"{dd-TRUE[e]:+.3f}" if (e in TRUE and dd is not None) else "—"
        print(f"| {e} | {EDGE[e][0]} → {EDGE[e][1]} | {c}/{s} ({tot}) | {n} | "
              f"{f'{dd:+.3f}' if dd is not None else '—'} | "
              f"{(f'{sd:.3f} (n={n})' if sd is not None else f'— (n={n})')} | {exp} | {err} | {v} | {NOTE.get(e,'—')} |")
    print(f"\n## T-2 · 고리 표  (닫힘 기준 |합| ≤ {CLOSE_TOL:.1f} kcal/mol — F-6 널 문턱과 동일)")
    print("| 고리 | 구성 엣지 (도는 순서·부호) | 계산 가능? | 합 | 전파 불확실도 | 닫힘 |")
    print("|---|---|---|---|---|---|")
    okc=0; blocked=[]
    for c,eg,path in CYCLES:
        nodes=path.split(); parts=[]; tot=0.0; var=0.0; can=True; why=[]
        for i in range(len(nodes)):
            a,b=nodes[i],nodes[(i+1)%len(nodes)]
            hit=None
            for e,(x,y) in EDGE.items():
                if (x,y)==(a,b): hit=(e,+1); break
                if (y,x)==(a,b): hit=(e,-1); break
            if hit is None: can=False; why.append(f"엣지없음 {a}→{b}"); parts.append(f"?{a}→{b}"); continue
            e,sgn=hit; parts.append(f"{'+' if sgn>0 else '−'}{e}")
            rr=d.get(e,{})
            if rr.get("n",0)==0 or rr.get("dd") is None:
                can=False; why.append(f"{e} n=0"); continue
            tot+=sgn*rr["dd"]
            if rr.get("sd") is not None: var+=rr["sd"]**2
        prop=math.sqrt(var) if var>0 else None
        if can:
            okc+=1
            close=("닫힘" if abs(tot)<=CLOSE_TOL else "★미닫힘")+f" (|{tot:+.3f}| vs {CLOSE_TOL:.1f})"
            tots=f"{tot:+.3f}"; props=f"{prop:.3f}" if prop else "—"
        else:
            blocked.append(c); close="—"; tots="—"; props="—"
        print(f"| 고리{c} | {' '.join(parts)} | {'예' if can else '★불가 ('+', '.join(why)+')'} | {tots} | {props} | {close} |")
    print(f"\n★N11 이 낀 고리 = 계산 불가 {len(blocked)}개 (고리 {', '.join(blocked) if blocked else '없음'})")
    # ★H-1 다리 확인 — N11 을 뺀 그래프 연결성
    es={e:v for e,v in EDGE.items() if d.get(e,{}).get("n",0)>0}
    nodes={n for p in es.values() for n in p}
    adj={n:set() for n in nodes}
    for a,b in es.values(): adj[a].add(b); adj[b].add(a)
    seen=set(); st=[next(iter(nodes))] if nodes else []
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x); st+=list(adj[x]-seen)
    allnodes={n for p in EDGE.values() for n in p}
    print(f"★H-1 다리 — n>0 엣지만으로 이룬 그래프: 노드 {len(seen)}/{len(allnodes)} 연결"
          f"{' · ★연결됨' if len(seen)==len(nodes) and nodes else ''}"
          f"{' · ★고립 노드 '+', '.join(sorted(allnodes-seen)) if allnodes-seen else ''}")
    print(f"\n★R-13 — 사이클 닫힘은 ★노드 전위형 편향 P_i 를 잡지 못한다(사이클에서 상쇄). 앵커 실험값 대조(F-14)로 본다.")
    print("\n## T-3 · 요약")
    print(f"1. 마감가능 {cnt['마감가능']} · 연장필요 {cnt['연장필요']} · 판정불가 {cnt['판정불가']}  (합 16)")
    line=[]
    for e in ("N01","N04","N07"):
        r=d.get(e,{}); dd=r.get("dd"); n=r.get("n",0)
        line.append(f"{e} 오차 {f'{dd-TRUE[e]:+.3f}' if dd is not None else '—'} (n={n})")
    print(f"2. 검문소 — {' · '.join(line)}")
    print(f"3. 계산 가능한 고리 {okc}/7")
    return 0

if __name__=="__main__": sys.exit(main())
