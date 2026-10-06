#!/usr/bin/env python
"""edge5(ejm_31+ejm_50) RESP 하이브리드 1단계 검증.
   원자 수 · dummy 방향 · bonded A/B(dihedral 은 합산) · B 상태 열 존재."""
import numpy as np
import sys
DIR=sys.argv[1] if len(sys.argv)>1 else '5_ejm31_ejm50'
NEW=f'/home/nudge/Project/TYK2-PILOT/04_hybrid_resp/{DIR}/merged.itp'
OLD=f'/home/nudge/Project/TYK2-PILOT/04_hybrid/{DIR}/merged.itp'
NIDX={'bonds':2,'pairs':2,'angles':3,'dihedrals':4}
def scan(f):
    sec=None; at=[]; rows={k:[] for k in NIDX}
    for ln in open(f):
        t=ln.split(';')[0].strip()
        if t.startswith('['): sec=t.strip('[] '); continue
        if not t: continue
        p=t.split()
        if sec=='atoms' and p[0].isdigit(): at.append((p[1],p[4],p[8],float(p[6]),float(p[9])))
        if sec in NIDX and p and p[0].isdigit(): rows[sec].append(p)
    return at,rows
na,ra=scan(NEW); no,ro=scan(OLD)
bad=[]
print(f"  원자 수  new {len(na)} · old {len(no)}")
if len(na)!=len(no): bad.append("원자 수 불일치")
# dummy 방향
dn=[(a[1],a[0],a[2]) for a in na if a[0].startswith('DUM') or a[2].startswith('DUM')]
do=[(a[1],a[0],a[2]) for a in no if a[0].startswith('DUM') or a[2].startswith('DUM')]
print(f"  dummy new {len(dn)} · old {len(do)}")
for x,y in zip(dn,do):
    if x!=y: bad.append(f"dummy 방향 {x} vs {y}")
if len(dn)!=len(do): bad.append("dummy 개수 불일치")
# 순전하
qa=sum(a[3] for a in na); qb=sum(a[4] for a in na)
print(f"  순전하 A {qa:+.6f} · B {qb:+.6f}")
if abs(qa)>1e-4 or abs(qb)>1e-4: bad.append(f"순전하 {qa:+.5f}/{qb:+.5f}")
# bonded: 개수 + B 상태 열 존재
for s in NIDX:
    nb=sum(1 for p in ra[s] if len(p[NIDX[s]+1:])>0 and len(p[NIDX[s]+1:])%2==0)
    print(f"  {s}: new {len(ra[s])} (B열 {nb}) · old {len(ro[s])}")
    if len(ra[s])!=len(ro[s]): bad.append(f"{s} 개수 {len(ra[s])}≠{len(ro[s])}")
    if s!='pairs' and nb!=len(ra[s]): bad.append(f"{s} B열 {nb}/{len(ra[s])}")
# dihedral 은 합산 비교 (pmx 가 같은 쿼드를 여러 줄로 쪼갠다)
import collections
def dsum(rows):
    d=collections.defaultdict(lambda:[0.0,0.0])
    for p in rows['dihedrals']:
        rest=p[5:]
        if len(rest)!=6: continue
        k=(tuple(p[:4]),p[4],rest[2],rest[5])
        d[k][0]+=float(rest[1]); d[k][1]+=float(rest[4])
    return d
da,do_=dsum(ra),dsum(ro)
print(f"  dihedral 합산 그룹 new {len(da)} · old {len(do_)}")
if set(da)!=set(do_): bad.append("dihedral 쿼드 집합 불일치")
print("★검증 실패: "+" · ".join(bad) if bad else "  ✅1단계 검증 전부 통과")
