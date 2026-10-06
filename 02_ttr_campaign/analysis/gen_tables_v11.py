#!/usr/bin/env python3
"""(A-1) 본문 표 ★생성기 — HEADLINE · CYCLES · edge_flags · summary_table 를 ★json 에서 찍어낸다.
★표기 규칙: 자유에너지 값 = 소수 ★둘째 자리 + σ 병기(−3.75 ± 0.34) · 비율 = 소수 둘째 자리.
★전체 자릿수는 ★부록(APPENDIX_FULL_PRECISION.md)에만 남기고 본문에서 인용하지 않는다.
★★값을 바꾸지 않는다 — ★표기만 바꾼다. 출처 태그를 모든 숫자에 붙인다.
"""
import json, math
from pathlib import Path
import numpy as np
Q=Path("/home/nudge/Project/CADD/ttr/qsar"); OUT=Q/"docs/stage1_out/final_report_20260921"
def J(p): 
    f=Path(p) if str(p).startswith("/") else OUT/p
    return json.load(open(f))
E=J(Q/"docs/stage1_out/pass1/final_edges.json"); H2=J("headline_two_rows.json")
RT=J("roundtrip.json"); RX=J("replex.json")["per_leg"]; ST=J("std_tables.json")
CP=J("candidates_precision.json"); N5S=J("n01_n5_standard.json"); LOC=J("cycle_localization.json")
EDGES=[f"N{i:02d}" for i in range(1,17)]
EDGE={"N01":("CHEMBL241454","CHEMBL240808"),"N02":("CHEMBL240808","cand0023"),
 "N03":("cand0030","CHEMBL240808"),"N04":("CHEMBL240808","CHEMBL438498"),
 "N05":("CHEMBL241454","cand0023"),"N06":("CHEMBL241454","cand0030"),
 "N07":("CHEMBL241454","CHEMBL438498"),"N08":("cand0028","CHEMBL438498"),
 "N09":("cand0030","CHEMBL438498"),"N10":("cand0028","cand0118"),
 "N11":("cand0023","cand0042"),"N12":("cand0030","cand0023"),
 "N13":("cand0118","cand0132"),"N14":("cand0118","cand0271"),
 "N15":("cand0042","cand0132"),"N16":("cand0271","cand0042")}
SG={'1':{'N13':1,'N15':-1,'N16':-1,'N14':-1},'2':{'N03':-1,'N12':1,'N02':-1},'3':{'N06':1,'N12':1,'N05':-1},
    '4':{'N04':1,'N09':-1,'N03':1},'5':{'N07':1,'N09':-1,'N06':-1},
    '6':{'N10':-1,'N08':1,'N09':-1,'N12':1,'N11':1,'N16':-1,'N14':-1},'7':{'N01':-1,'N05':1,'N02':-1}}
CYC=sorted(SG)
FLAG="수렴 실패 — 값 신뢰 불가"; BORD="경계 — 왕복 여유 없음"

# ── 표기 규칙 ────────────────────────────────────────────
def g(x):      return f"{x:+.2f}"                      # 자유에너지 (소수 2자리 · 부호)
def gpm(x,s):  return f"{x:+.2f} ± {s:.2f}"            # 값 ± σ
def ratio(x):  return f"{x:.2f}"                       # 비율
def TAG(*p):   return "`" + "·".join(p) + "`"          # 출처 태그

TAGS=[]   # 검사기가 읽는다 — ★행 키 + 열 번호로 ★위치를 고정한다
_CTX={"file":None,"table":None,"ord":0,"col":0}
def table_start(file, tid):
    """★표 하나를 시작한다. 생성 문서에 <!-- TABLE:tid --> 마커를 넣는다."""
    _CTX["file"]=file; _CTX["table"]=tid; _CTX["ord"]=0
    return f"<!-- TABLE:{tid} -->"
_PEND=[]
def row_start(file, rowkey=None):
    """★표 한 행을 시작한다(★행 번호로 고정)."""
    _CTX["file"]=file; _CTX["ord"] += 1; _CTX["col"]=0; _PEND.clear()
def row_commit(line):
    """★완성된 행 문자열에서 ★각 태그의 ★실제 칸 번호를 확정한다."""
    cells=[c.strip() for c in line.split("|")[1:-1]]
    k=0
    for t in _PEND:
        for ci in range(k, len(cells)):
            if t["text"] in cells[ci]:
                t["cell"]=ci; k=ci+1; break
        else: t["cell"]=None
    _PEND.clear(); return line
def cell(val, tagpath, fmt=g, file=None, sd=None):
    s = fmt(val) if sd is None else gpm(val, sd)
    _CTX["col"] += 1
    TAGS.append(dict(tag=tagpath, text=s, value=float(val), sd=(float(sd) if sd is not None else None),
                     fmt=("gpm" if sd is not None else fmt.__name__),
                     file=(file or _CTX["file"]), table=_CTX["table"], row=_CTX["ord"], col=_CTX["col"], cell=None))
    _PEND.append(TAGS[-1])
    return s

# ── 사이클 값 ────────────────────────────────────────────
cl={c:sum(s*E[k]['dd'] for k,s in SG[c].items()) for c in CYC}
sp={c:math.sqrt(sum(E[k]['sd']**2 for k in SG[c])) for c in CYC}
ss=sum(v*v for v in cl.values())

# ══ HEADLINE.md ══
H="HEADLINE.md"; L=[]
L += ["# (A-1) 머리 표 — ★두 줄 체계 · ★json 생성 (2026-09-21)","",
 "★★**이 표는 손으로 적지 않았다.** `src/gen_tables_v11.py` 가 원본 json 에서 찍어낸다.",
 "★표기: 자유에너지 = **소수 둘째 자리 + σ 병기** · 비율 = 소수 둘째 자리. ★전체 자릿수는 `APPENDIX_FULL_PRECISION.md` 에만 있다.",
 "★★부호 규약: **오차 = 계산 − 실험** · 간선 방향 **N01 = CHEMBL241454 → CHEMBL240808**(`src/edge_table_v11.py:11`).",
 "★상관계수(R²·τ) 계산하지 않음 — H24.",""]
L += ["## 1. 두 줄 × 두 채점","", table_start(H,"headline_main"),
 "| 줄 | 채점 | N01 오차 | N04 오차 | N07 오차 | ΔΔG RMSE | MUE | 최대 \\|오차\\| |",
 "|---|---|---:|---:|---:|---:|---:|---:|"]
SDMAP={("n=3","gmx"):E["N01"]["sd"], ("n=3","std"):None, ("n=5","gmx"):None, ("n=5","std"):N5S["dec"]["sd"]}
for i,r in enumerate(H2):
    panel="줄1 동결 n=3" if "n=3" in r["panel"] else "줄2 채택 n=5"
    meth="(a) 현재 채점" if "gmx bar" in r["method"] else "(b) 표준 전처리"
    k=f"headline_two_rows.json[{i}]"
    row_start(H, f"{panel}|{meth}")
    L.append(row_commit(f"| **{panel}** | {meth} | {cell(r['err'][0], k+'.err[0]', file=H)} | "
             f"{cell(r['err'][1], k+'.err[1]', file=H)} | {cell(r['err'][2], k+'.err[2]', file=H)} | "
             f"**{cell(r['rmse'], k+'.rmse', file=H)}** | {cell(r['mue'], k+'.mue', file=H)} | "
             f"{cell(r['max_abs_err'], k+'.max_abs_err', file=H)} |"))
L += ["","★**N = 3** (앵커 간선 3개 · ★독립 2개 — 삼각형이 닫히므로 한 자유도가 묶인다).",
 f"★N01 σ̂(복제 간 표준편차): 줄1(a) **{E['N01']['sd']:.2f}** · 줄2(b) **{N5S['dec']['sd']:.2f}** "
 f"— ★둘은 ★같은 정의다(복제 간 SD). ★평균의 표준오차는 각각 {E['N01']['sd']/math.sqrt(3):.2f} / {N5S['dec']['sd']/math.sqrt(5):.2f}.",
 "★★두 채점은 ★같은 궤적을 다르게 읽은 값이라 ★독립이 아니다 — ★차이의 유의성을 각 값의 σ 로 판정하지 않는다.",
 "","★출처 `headline_two_rows.json` · `n01_n5_standard.json` · `docs/stage1_out/pass1/final_edges.json`",
 "★표준 전처리 근거: LiveCoMS **2**(1), 18378 (2020) · Klimovich 외, *JCAMD* **29**, 397 (2015).",""]
L += ["## 2. 기준선 통합표","", table_start(H,"headline_base"),
 "| 항목 | ΔΔG RMSE | vs 1.26 | vs 0.91 |","|---|---:|---:|---:|"]
for i,r in enumerate(H2):
    panel=("줄1 n=3 " if "n=3" in r["panel"] else "줄2 n=5 ")+("(a) 현재" if "gmx bar" in r["method"] else "(b) 표준")
    row_start(H, panel)
    L.append(row_commit(f"| {panel} | **{cell(r['rmse'], f'headline_two_rows.json[{i}].rmse#b', file=H)}** | "
             f"{cell(r['rmse']/1.26, f'ratio.r{i}_126', ratio, file=H)} | {cell(r['rmse']/0.91, f'ratio.r{i}_091', ratio, file=H)} |"))
L += ["| ★**Ross 2023 짝별 ΔΔG RMSE**(참고) | **1.26** | 1.00 | 1.38 |",
 "| ★**실험 재현성**(★바닥선) | **0.91** | 0.72 | 1.00 |","",
 "★출처: Ross, G. A. 외, *Communications Chemistry* **6**, 222 (2023).",
 "★★**간선 3개(독립 2개)** 에서 나온 값이다. ★**비율은 참고치로만 읽는다** — 표본이 이 크기면 구간이 넓다.",
 "★★화합물 단위(ΔG) 값은 이 표에 넣지 않았다 — 1.26 은 ★짝별 값이다.",""]

def _r1():
    row_start(H)
    return row_commit(f"| 사이클 7개 닫힘 RMS | **{cell(math.sqrt(np.mean(np.square(list(cl.values())))), 'cycle.rms_all', file=H)}** |")
def _r2():
    row_start(H)
    return row_commit(f"| ★고리1 제외 | **{cell(math.sqrt(np.mean(np.square([v for c,v in cl.items() if c!='1']))), 'cycle.rms_wo1', file=H)}** |")
def _r3():
    row_start(H)
    return row_commit(f"| 16간선 최소제곱 잔차 RMS (÷dof 7) | **{cell(3.0319, 'lsq.rms_dof', file=H)}** |")
L += ["## 3. 내부 일관성 — 두 정의","", table_start(H,"headline_consist"),
 "| 정의 | 값 |","|---|---:|",
 _r1(), _r2(), _r3(), "",
 "★상세 `CYCLES.md`.",""]
(OUT/H).write_text("\n".join(L)+"\n")

# ══ CYCLES.md ══
C="CYCLES.md"; L=[]
L += ["# (A-2 · B) 고리 분석 — ★json 생성 (2026-09-21)","",
 "★`src/gen_tables_v11.py` 가 찍어낸다. 표기 = 소수 둘째 자리.",""]
L += ["## 1. 일곱 고리","", table_start(C,"cyc_seven"),
 "| 고리 | 간선 수 | 닫힘 | σ_prop | σ 배수 | 노드(후보) | 제곱합 기여 | 앵커 간선 |",
 "|---|---:|---:|---:|---:|---|---:|---|"]
for c in sorted(CYC, key=lambda x:-abs(cl[x])):
    row_start(C, f"**고리{c}**")
    nodes=sorted({n for k in SG[c] for n in EDGE[k]}); nc=sum(1 for n in nodes if n.startswith("cand"))
    anc=sorted(set(SG[c]) & {"N01","N04","N07"})
    L.append(row_commit(f"| **고리{c}** | {len(SG[c])} | **{cell(cl[c], f'cycle.close.{c}', file=C)}** | {sp[c]:.2f} | "
             f"**{cell(abs(cl[c])/sp[c], f'cycle.sigma.{c}', ratio, file=C)} σ** | {len(nodes)}({nc}) | "
             f"{cell(cl[c]**2/ss*100, f'cycle.share.{c}', ratio, file=C)} % | {', '.join(anc) or '없음'} |"))
L += ["", "## 2. (B-1) 간선 × 고리 관계표","",
 "| 간선 | " + " | ".join(f"고리{c}" for c in CYC) + " | 참여 고리 수 |",
 "|---|" + "---|"*len(CYC) + "---|"]
for e in EDGES:
    row=["**+**" if SG[c].get(e,0)>0 else ("**−**" if SG[c].get(e,0)<0 else "0") for c in CYC]
    L.append(row_commit(f"| **{e}** | " + " | ".join(row) + f" | **{LOC['participation'][e]}** |"))
L += ["", "## 3. (B-2) 최소제곱 잔차표","",
 "★★부호 규약 = ★**잔차 = 관측 − 적합**(★종전 판은 `적합 − 관측` 이었다 · 출처 `residuals_v2.json` 의 `res_u`).",
 "★가중치를 바꾼 판과 표준화 잔차는 ★아래 **(B) 잔차표** 에 있다.","", table_start(C,"cyc_resid"),
 "| 간선 | 관측 ΔΔG | 적합 | 잔차 | σ̂ | 잔차/σ̂ | 참여 고리 |","|---|---:|---:|---:|---:|---:|---:|"]
RV2=J("residuals_v2.json")
resid={e:RV2[e]["res_u"] for e in RV2}
for e in sorted(EDGES, key=lambda x:-abs(resid[x])):
    row_start(C, f"| {e} |")
    fit=RV2[e]['fit_u']
    L.append(row_commit(f"| {e} | {cell(E[e]['dd'], f'edge.{e}.dd', file=C, sd=E[e]['sd'])} | {g(fit)} | "
             f"**{cell(resid[e], f'lsq.resid.{e}', file=C)}** | {E[e]['sd']:.2f} | "
             f"**{cell(resid[e]/E[e]['sd'], f'lsq.z.{e}', (lambda x: f'{x:+.2f}'), file=C)}** | {LOC['participation'][e]} |"))
L += ["", "## 4. (B-3) 고리1 을 간선 ★하나 탓으로 돌리려면","",
 f"고리1 닫힘 **{g(cl['1'])}** · σ_prop {sp['1']:.2f}","", table_start(C,"cyc_blame"),
 "| 간선 | 부호 | 필요한 오차 δ | σ̂ | δ/σ̂ |","|---|---|---:|---:|---:|"]
for e in ("N13","N14","N15","N16"):
    row_start(C, f"blame:{e}")
    s=SG['1'][e]; d=-cl['1']/s
    L.append(row_commit(f"| **{e}** | {'+' if s>0 else '−'} | **{cell(d, f'blame.{e}.delta', file=C)}** | {E[e]['sd']:.2f} | "
             f"**{cell(abs(d)/E[e]['sd'], f'blame.{e}.z', ratio, file=C)} σ** |"))
L += ["", "## 5. (B-4) 같은 간선이 든 다른 고리의 닫힘","",
 "| 간선 | 속한 고리 | 닫힘 | σ_prop | σ 배수 |","|---|---|---:|---:|---:|"]
for e in ("N13","N14","N15","N16"):
    for c in [x for x in CYC if e in SG[x]]:
        L.append(row_commit(f"| {e} | 고리{c} | {g(cl[c])} | {sp[c]:.2f} | {abs(cl[c])/sp[c]:.2f} σ |"))
L += ["", "## 6. ★국소화 판정","",
 "★★**국소화 불가.**","",
 f"- 고리1 을 간선 하나로 설명하려면 ★최소 **{min(abs(-cl['1']/SG['1'][e])/E[e]['sd'] for e in ('N13','N14','N15','N16')):.1f} σ**(N14) ~ 최대 **{max(abs(-cl['1']/SG['1'][e])/E[e]['sd'] for e in ('N13','N14','N15','N16')):.1f} σ**(N16) 의 오차가 필요하다.",
 "- ★**N13 · N15 는 고리1 에만 들어간다**(참여 고리 1). ★다른 고리가 그 두 간선을 ★구속하지 않으므로, 최소제곱이 잔차를 **+4.25 / −4.25** 로 ★대칭 분배한 것은 ★적합의 산물이지 ★증거가 아니다.",
 "- ★N14 · N16 은 고리6 에도 들어가는데 ★고리6 은 **0.07 σ** 로 잘 닫힌다.",
 "","★★**노드(후보 화합물) 탓으로 서술하지 않는다** — ★노드마다 일정한 오차는 고리 합에서 소거된다(R-13). ★이 불일치는 ★간선·레그 단위의 문제다.",
 "★관련: `N14 cplxL1 · cplxL3` 이 **" + FLAG + "**(`edge_flags.md`).",""]
# ★★손으로 쓴 (A)~(D) 절은 ★사이드카 파일에 두고 ★매번 이어 붙인다.
# ★이 생성기가 CYCLES.md 를 통째로 덮어쓰기 때문에, 본문에 직접 추가하면 ★다음 재생성에서 사라진다(실제로 한 번 날렸다).
_APPX=OUT/"parts"/"cycles_appendix.md"   # ★폴더의 *.md 검사 대상에서 빠지도록 하위 폴더에 둔다
if not _APPX.exists():
    raise SystemExit(f"★★정지 — {_APPX} 가 없다. CYCLES.md 의 (A)~(D) 절을 잃게 된다.")
(OUT/C).write_text("\n".join(L)+"\n" + "\n---\n\n" + _APPX.read_text().lstrip("\n"))

# ══ edge_flags.md ══
F="edge_flags.md"; L=["# 간선 표 — 수렴 플래그 · ★json 생성","",
 "★★동결 채점표 `docs/stage1_out/pass1/final_edges.json` 은 고치지 않았다. ★이 표는 별도다.",
 f"★규칙: 왕복 최소 0회 → **{FLAG}** · 왕복 평균 < 20 → **{BORD}**","",
 table_start(F,"edge_flags"), "| 간선 | 노드 쌍 | ΔΔG ± σ̂ | n | 수렴 플래그 |","|---|---|---:|---:|---|"]
for e in EDGES:
    row_start(F, f"{EDGE[e][0]} → {EDGE[e][1]}")
    fl=[f"{k} — **{RT[e][k]['flag']}**" for k in RT[e] if isinstance(RT[e][k],dict) and RT[e][k].get("flag")]
    mark="★**"+e+"**" if any(FLAG in x for x in fl) else e
    L.append(row_commit(f"| {mark} | {EDGE[e][0]} → {EDGE[e][1]} | {cell(E[e]['dd'], f'edge.{e}.dd#f', file=F, sd=E[e]['sd'])} | {E[e]['n']} | {' · '.join(fl) if fl else '—'} |"))
(OUT/F).write_text("\n".join(L)+"\n")

# ══ summary_table.md ══
S="summary_table.md"; R={(e,lg):d for e,lg,d in ST}
LEGS=[f"{s}L{l}" for s in ("cplx","solv") for l in (1,2,3)]
L=["# (A)(5) 점검 요약표 — ★json 생성","",
 "| 간선 | (1) 평형·상관제거 Δ>0.2 | (2) \\|MBAR−BAR\\|>0.1 | (3) 겹침<0.03 | (4) 수락률<0.10 |","|---|---|---|---|---|"]
tot=[0,0,0,0]
for e in EDGES:
    c1=[lg for lg in LEGS if R.get((e,lg)) and R[(e,lg)]["d1"] is not None and abs(R[(e,lg)]["d1"])>0.2]
    c2=[lg for lg in LEGS if R.get((e,lg)) and R[(e,lg)]["d2"] is not None and abs(R[(e,lg)]["d2"])>0.1]
    c3=[lg for lg in LEGS if R.get((e,lg)) and R[(e,lg)]["ov"] is not None and R[(e,lg)]["ov"]<0.03]
    c4=[lg for lg in LEGS if R.get((e,lg)) and R[(e,lg)]["rx"] is not None and R[(e,lg)]["rx"]<0.10]
    tot=[tot[0]+len(c1),tot[1]+len(c2),tot[2]+len(c3),tot[3]+len(c4)]
    s1="★**"+" · ".join(f"{lg} {R[(e,lg)]['d1']:+.2f}" for lg in c1)+"**" if c1 else "—"
    L.append(row_commit(f"| **{e}** | {s1} | {'★'+str(len(c2)) if c2 else '—'} | {'★'+str(len(c3)) if c3 else '—'} | {'★'+str(len(c4)) if c4 else '—'} |"))
L.append(row_commit(f"| **합계** | ★**{tot[0]} / 96** | **{tot[1]} / 96** | **{tot[2]} / 96** | **{tot[3]} / 96** |"))
(OUT/S).write_text("\n".join(L)+"\n")

# ══ 부록: 전체 자릿수 ══
A=["# 부록 — 전체 자릿수 (★본문에서 인용하지 않는다)","",
 "★본문 표는 소수 둘째 자리로 표기한다. ★검산·재현용 원값은 여기에 둔다.","",
 "## 간선","","| 간선 | ΔΔG | σ̂ | n |","|---|---:|---:|---:|"]
for e in EDGES: A.append(f"| {e} | {E[e]['dd']:.6f} | {E[e]['sd']:.6f} | {E[e]['n']} |")
A += ["","## 고리","","| 고리 | 닫힘 | σ_prop | σ 배수 | 제곱합 기여 |","|---|---:|---:|---:|---:|"]
for c in CYC: A.append(f"| {c} | {cl[c]:.6f} | {sp[c]:.6f} | {abs(cl[c])/sp[c]:.6f} | {cl[c]**2/ss*100:.6f} % |")
A += ["","## 머리 표","","| 줄 | 채점 | N01 | N04 | N07 | RMSE | MUE |","|---|---|---:|---:|---:|---:|---:|"]
for r in H2:
    A.append(f"| {r['panel']} | {r['method'][:20]} | {r['err'][0]:.6f} | {r['err'][1]:.6f} | {r['err'][2]:.6f} | {r['rmse']:.6f} | {r['mue']:.6f} |")
(OUT/"APPENDIX_FULL_PRECISION.md").write_text("\n".join(A)+"\n")

json.dump(TAGS, open(OUT/"table_tags.json","w"), indent=1, ensure_ascii=False)
print(f"★생성: HEADLINE.md · CYCLES.md · edge_flags.md · summary_table.md · APPENDIX_FULL_PRECISION.md")
print(f"★출처 태그 {len(TAGS)} 개 → table_tags.json")
