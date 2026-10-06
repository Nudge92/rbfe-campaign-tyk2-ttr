#!/usr/bin/env python
"""RESP 하이브리드 1단계 검증 — 원자 수·dummy 방향·bonded A/B·B 상태 열 존재."""
import os
HD='/home/nudge/Project/TYK2-PILOT/04_hybrid_resp'
OLD='/home/nudge/Project/TYK2-PILOT/04_hybrid'
E={'edge1':'1_ejm31_ejm42','edge2':'2_ejm46_ejm42','edge3':'3_ejm46_ejm50',
   'edge4':'4_ejm46_ejm47','edge5':'5_ejm31_ejm50'}
NIDX={'bonds':2,'pairs':2,'angles':3,'dihedrals':4}
def scan(f):
    sec=None; at=[]; nb={'bonds':0,'angles':0,'dihedrals':0}; bstate={'bonds':0,'angles':0,'dihedrals':0}
    for ln in open(f):
        t=ln.split(';')[0].strip()
        if t.startswith('['): sec=t.strip('[] '); continue
        if not t: continue
        p=t.split()
        if sec=='atoms' and p[0].isdigit():
            at.append((p[1],p[8],float(p[6]),float(p[9])))
        if sec in nb and p and p[0].isdigit():
            nb[sec]+=1
            rest=p[NIDX[sec]+1:]
            if rest and len(rest)%2==0: bstate[sec]+=1
    return at,nb,bstate
print(f"{'edge':<8}{'원자':>5}{'A더미':>6}{'B더미':>6}{'qA':>10}{'qB':>10}"
      f"{'bonds':>7}{'angles':>7}{'dih':>7}{'B열':>6}  판정")
ok_all=True
for tag,d in E.items():
    f=f'{HD}/{d}/merged.itp'
    if not os.path.exists(f): print(f"{tag:<8} ★파일 없음"); ok_all=False; continue
    at,nb,bs=scan(f)
    da=sum(1 for a in at if a[0].startswith('DUM')); db=sum(1 for a in at if a[1].startswith('DUM'))
    qa=sum(a[2] for a in at); qb=sum(a[3] for a in at)
    # 원본과 원자 수·더미 구성 비교
    oat,onb,_=scan(f'{OLD}/{d}/merged.itp')
    oda=sum(1 for a in oat if a[0].startswith('DUM')); odb=sum(1 for a in oat if a[1].startswith('DUM'))
    prob=[]
    if len(at)!=len(oat): prob.append(f"원자수 {len(at)}≠{len(oat)}")
    if (da,db)!=(oda,odb): prob.append(f"더미 {da}/{db}≠{oda}/{odb}")
    if abs(qa)>1e-4 or abs(qb)>1e-4: prob.append(f"순전하 {qa:+.5f}/{qb:+.5f}")
    for s in nb:
        if nb[s]!=bs[s]: prob.append(f"{s} B열 {bs[s]}/{nb[s]}")
        if nb[s]!=onb[s]: prob.append(f"{s} 수 {nb[s]}≠{onb[s]}")
    if prob: ok_all=False
    print(f"{tag:<8}{len(at):>5}{da:>6}{db:>6}{qa:>10.6f}{qb:>10.6f}"
          f"{nb['bonds']:>7}{nb['angles']:>7}{nb['dihedrals']:>7}"
          f"{'전부':>6}  {'통과' if not prob else '★'+' '.join(prob)}")
print("검증 " + ("전부 통과" if ok_all else "★실패 있음"))
