#!/usr/bin/env python3
"""노드 전위 최소제곱 + ★두 양의 분리 보고. ★재계산 0 — final_edges.json 만 읽는다.

★★설계 원칙(감독 지시 2026-09-21):
   ★'내부 불일치 잔차' 와 '실험 대비 오차' 를 ★같은 표에 넣지 않는다. ★다른 양이다.
   - 내부 불일치 잔차 : 계산이 ★자기 자신과 안 맞는 정도. 실험을 전혀 안 쓴다.
   - 실험 대비 오차   : 계산과 ★실험값의 차. ★실험값 자체의 불확실도를 따로 붙인다.
★F-10 / H24 — 상대 자유에너지에 ★상관계수(ρ·τ)를 쓰지 않는다.
"""
import json, sys, math
import numpy as np

EDGE={"N01":("CHEMBL241454","CHEMBL240808"),"N02":("CHEMBL240808","cand0023"),
 "N03":("cand0030","CHEMBL240808"),"N04":("CHEMBL240808","CHEMBL438498"),
 "N05":("CHEMBL241454","cand0023"),"N06":("CHEMBL241454","cand0030"),
 "N07":("CHEMBL241454","CHEMBL438498"),"N08":("cand0028","CHEMBL438498"),
 "N09":("cand0030","CHEMBL438498"),"N10":("cand0028","cand0118"),
 "N11":("cand0023","cand0042"),"N12":("cand0030","cand0023"),
 "N13":("cand0118","cand0132"),"N14":("cand0118","cand0271"),
 "N15":("cand0042","cand0132"),"N16":("cand0271","cand0042")}
REF="CHEMBL241454"

# ★실험 ΔΔG — ΔΔG(A→B) = RT·ln(Kd_B/Kd_A) @ 310 K (우리 MD 온도 · `ANCHOR_TRIANGLE_PREREG.md:44-52`)
EXP={"N01":0.000,"N04":-3.022,"N07":-3.022}
# ★실험값 불확실도 ★성분별 — ★전부 출처가 있는 것만
EXP_BUDGET={
 # ★2026-09-21 갱신 — ★원전 SE 로 교체. 반올림 추정치(±0.002/±0.157)는 ★폐기한다.
 # ★출처: /home/nudge/Project/CADD/ttr/docs/gupta2007.pdf Table 1 (3쪽)
 #        S.N. 6 = 0.2±0.08 · S.N. 7 = 27.0±8.03 · S.N. 10 = 27±4.10  (nM)
 # ★전파: σ(ΔΔG) = RT·sqrt((SE_A/Kd_A)² + (SE_B/Kd_B)²)
 # ★반올림 성분은 SE 에 이미 포함되므로 ★이중 계상하지 않는다.
 "N01":[("★원전 SE 전파 (Gupta 2007 Table 1)",0.2057,"27±4.10 / 27.0±8.03 nM · gupta2007.pdf 3쪽"),
        ("온도 미기재 298.15↔310 K",0.000,"ANCHOR_TRIANGLE_PREREG.md:45-48 · ★ΔΔG=0 이라 온도에 무관")],
 "N04":[("★원전 SE 전파 (Gupta 2007 Table 1)",0.3071,"27.0±8.03 / 0.2±0.08 nM · gupta2007.pdf 3쪽"),
        ("온도 미기재 298.15↔310 K",0.115,"ANCHOR_TRIANGLE_PREREG.md:48 원문")],
 "N07":[("★원전 SE 전파 (Gupta 2007 Table 1)",0.2636,"27±4.10 / 0.2±0.08 nM · gupta2007.pdf 3쪽"),
        ("온도 미기재 298.15↔310 K",0.115,"ANCHOR_TRIANGLE_PREREG.md:48 원문")],
}
EXP_UNQUANT=[
 ("Kd1/Kd2 자리 귀속","240808 은 Kd1 27 / ★Kd2 3.0 nM 로 역전 · Adair 2자리 적합 → 두 27.0 이 같은 자리인지 ★미확인",
  "ERRATA.md:24(A-6) · PREREG_FROZEN.md:59"),
 ("양성자화 pH 4.4 vs 계산 pH 7.4","실험은 ~61 % 음이온 혼합물 · 계산은 100 % 음이온. ★상한 ~1.0 kcal/mol · 추정 보정 ~0",
  "PROTONATION_0B_20260921.md · ef/protonate.py:1"),
 ("단백질 V30M vs 실험 WT","계산 6E72=V30M · 실험값은 변이체 제외(WT 통일). ★간선 양쪽에서 1차 상쇄 · 2차 효과 ★미측정",
  "TARGET_DATA_STATUS.md:12-15 · QSAR_WORKFLOW.md:150"),
]

def solve(edges, dd):
    nodes=sorted({n for e in edges for n in EDGE[e]})
    idx={n:i for i,n in enumerate(nodes)}
    A=np.zeros((len(edges),len(nodes))); b=np.array([dd[e] for e in edges])
    for k,e in enumerate(edges):
        a,c=EDGE[e]; A[k,idx[a]]-=1; A[k,idx[c]]+=1
    # 기준 고정
    A=np.hstack([A, np.zeros((len(edges),0))])
    r=np.zeros((1,len(nodes))); r[0,idx[REF]]=1
    A2=np.vstack([A, r*1e6]); b2=np.concatenate([b,[0.0]])
    g,*_=np.linalg.lstsq(A2,b2,rcond=None)
    pred={e: g[idx[EDGE[e][1]]]-g[idx[EDGE[e][0]]] for e in edges}
    return nodes, idx, g, pred

def main(path):
    d=json.load(open(path))
    used=[e for e in EDGE if d.get(e,{}).get("dd") is not None]
    skip=[e for e in EDGE if e not in used]
    dd={e:d[e]["dd"] for e in used}
    sd={e:d[e].get("sd") for e in used}
    n ={e:d[e].get("n") for e in used}
    # ★게이트 — 기준 노드가 사용 간선에 없으면 ★조용히 통과시키지 않는다
    if len(used) < 2:
        print(f"★★중단: 사용 가능한 간선이 {len(used)} 개다(dd 가 None 이 아닌 것). 최소 2 필요.", file=sys.stderr); return 2
    covered={n for e in used for n in EDGE[e]}
    if REF not in covered:
        print(f"★★중단: 기준 노드 {REF} 가 사용 간선에 없다.", file=sys.stderr); return 2
    nodes, idx, g, pred = solve(used, dd)
    dof=len(used)-(len(nodes)-1)
    print(f"★입력 {path}")
    print(f"★사용 간선 {len(used)} 개" + (f" · ★제외 {skip} (dd=None)" if skip else " · ★제외 없음"))
    print(f"★노드 {len(nodes)} 개 · 기준 {REF} = 0 · ★자유도 = {len(used)} − ({len(nodes)}−1) = {dof}\n")

    print("## 표 1 — 노드 전위 (kcal/mol · 기준 0)  ★실험 안 씀")
    print(f"| 노드 | G |\n|---|---:|")
    for nd in sorted(nodes, key=lambda x:g[idx[x]]):
        print(f"| {nd} | {g[idx[nd]]:+.3f} |")

    print("\n## 표 2 — ★내부 불일치 잔차  ★실험 안 씀")
    print("★잔차 = 계산 ΔΔG − 노드전위 모형. ★'계산이 자기와 안 맞는 정도'다.")
    print(f"\n| 순위 | 간선 | 노드 쌍 | n | 계산 ΔΔG | 모형 | ★내부 잔차 | σ̂(회차간) | 잔차/σ̂ |")
    print("|---|---|---|---:|---:|---:|---:|---:|---:|")
    rows=sorted(used, key=lambda e:-abs(dd[e]-pred[e]))
    for i,e in enumerate(rows,1):
        r=dd[e]-pred[e]; s=sd[e]
        rs=f"{r/s:.2f}" if s else "—"
        print(f"| {i} | **{e}** | {EDGE[e][0]} → {EDGE[e][1]} | {n[e]} | {dd[e]:+.3f} | {pred[e]:+.3f} | **{r:+.3f}** | {s:.3f} | {rs} |" if s else
              f"| {i} | **{e}** | {EDGE[e][0]} → {EDGE[e][1]} | {n[e]} | {dd[e]:+.3f} | {pred[e]:+.3f} | **{r:+.3f}** | — | — |")
    res=np.array([dd[e]-pred[e] for e in used])
    # ★자체검사 — 최소제곱 해는 잔차가 설계행렬 열공간에 ★직교해야 한다
    Ac=np.zeros((len(used),len(nodes)))
    for k,e in enumerate(used):
        a,c=EDGE[e]; Ac[k,idx[a]]-=1; Ac[k,idx[c]]+=1
    orth=np.abs(Ac.T@res); free=[i for i,nd in enumerate(nodes) if nd!=REF]
    mx=orth[free].max()
    print(f"\n★자체검사(최소제곱 직교성): max|AᵀR| over 자유노드 = {mx:.2e}  " + ("OK" if mx<1e-6 else "★★실패 — 해가 틀렸다"))
    assert mx<1e-6, "최소제곱 해가 정규방정식을 만족하지 않는다"
    print(f"\n★잔차 RMS = **{math.sqrt((res**2).sum()/max(dof,1)):.3f}** kcal/mol (자유도 {dof} 로 나눔) · 잔차합 = {res.sum():+.3e}")

    print("\n---\n")
    print("## 표 3 — ★실험 대비 오차  ★★표 2 와 ★다른 양이다. 섞지 않는다.")
    print("★실험값이 있는 간선은 ★3 개뿐이다(앵커 삼각형).")
    print(f"\n| 간선 | 계산 ΔΔG | σ̂(계산) | 실험 ΔΔG | ★실험값 불확실도 | ★실험 대비 오차 | 오차가 실험 불확실도의 몇 배 |")
    print("|---|---:|---:|---:|---:|---:|---:|")
    for e in ("N01","N04","N07"):
        if e not in dd: 
            print(f"| **{e}** | — | — | {EXP[e]:+.3f} | — | ★계산값 없음 | — |"); continue
        comps=EXP_BUDGET[e]; u=math.sqrt(sum(c[1]**2 for c in comps))
        err=dd[e]-EXP[e]
        rat = f"{abs(err)/u:.1f}×" if u >= 0.05 else "★비율 무의미 (불확실도 ≈ 0)"
        print(f"| **{e}** | {dd[e]:+.3f} | {sd[e]:.3f} | {EXP[e]:+.3f} | ±{u:.3f} | **{err:+.3f}** | {rat} |")
    print("\n### 표 3-1 — 실험값 불확실도의 ★정량 성분 (출처 있는 것만)")
    print("| 간선 | 성분 | ±kcal/mol | 출처 |\n|---|---|---:|---|")
    for e,comps in EXP_BUDGET.items():
        for nm,v,src in comps: print(f"| {e} | {nm} | {v:.3f} | {src} |")
    print("\n### 표 3-2 — ★정량화 못 한 성분 (★숫자를 지어내지 않는다)")
    print("| 성분 | 내용 | 출처 |\n|---|---|---|")
    for nm,txt,src in EXP_UNQUANT: print(f"| ★{nm} | {txt} | {src} |")
    print("\n★★`src/final_report_v11.py:24` 의 σ(0.206 · 0.307 · 0.264)는 ★원전 SE 전파값임이 확인됐다(RISK-6 정정) — ★위 표 1행과 같은 값이다.")
    print("★★표 2 의 잔차와 표 3 의 오차를 ★더하거나 비교하지 말 것 — ★전자는 실험을 안 쓰고 후자는 실험이 전부다.")
    return 0

if __name__=="__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv)>1 else "docs/stage1_out/pass1/final_edges.json"))
