#!/usr/bin/env python
"""edge2R 검증: 3블록 방향이 뒤집혔을 때 소멸/생성 원자의 역할과 레그 순서."""
import sys
B='/home/nudge/Project/TYK2-PILOT'
def atoms(p):
    out=[]; on=False
    for ln in open(p):
        t=ln.split(';')[0].rstrip()
        if t.strip().startswith('['): on='atoms' in t; continue
        f=t.split()
        if on and len(f)>=8 and f[0].isdigit():
            out.append(dict(i=int(f[0]),tA=f[1],nm=f[4],qA=float(f[6]),
                            tB=f[8] if len(f)>8 else f[1], qB=float(f[9]) if len(f)>9 else float(f[6])))
    return out
M=atoms(f'{B}/04_hybrid/2_ejm46_ejm42/merged.itp')
L={n:atoms(f'{B}/07_step3_lambda/3block/edge2/edge2_leg{n}.itp') for n in (1,2,3)}
dumA=[a['i'] for a in M if a['tA'].startswith('DUM')]     # A(ejm_46)에서 dummy = 46→42 로 갈 때 생성
dumB=[a['i'] for a in M if a['tB'].startswith('DUM')]     # B(ejm_42)에서 dummy = 46→42 로 갈 때 소멸
common=[a['i'] for a in M if not a['tA'].startswith('DUM') and not a['tB'].startswith('DUM')]
print(f"merged.itp (A=ejm_46 · B=ejm_42)  원자 {len(M)}")
print(f"  DUM in A (46→42 에서 **생성**, 42→46 에서 **소멸**): {len(dumA)}  {dumA}")
print(f"  DUM in B (46→42 에서 **소멸**, 42→46 에서 **생성**): {len(dumB)}  {dumB}")
print(f"  공통: {len(common)}")
def dq(leg,idx):  # 레그 내 전하 변화량 합
    d={a['i']:(a['qA'],a['qB']) for a in L[leg]}
    return sum(abs(d[i][1]-d[i][0]) for i in idx)
print(f"\n{'레그':<6}{'설명':<34}{'Σ|Δq| dumA':>12}{'Σ|Δq| dumB':>12}{'Σ|Δq| 공통':>12}")
desc={1:'typeA=typeB=origA · q: origA→mid',2:'type A→B · q 고정(mid)',3:'typeA=typeB=origB · q: mid→origB'}
for n in (1,2,3):
    print(f"L{n:<5}{desc[n]:<34}{dq(n,dumA):12.4f}{dq(n,dumB):12.4f}{dq(n,common):12.4f}")
# 타입 변화 위치
for n in (1,2,3):
    ch=[a['i'] for a in L[n] if a['tA']!=a['tB']]
    print(f"  L{n} typeA≠typeB: {len(ch)}개 {ch[:12]}")
print(f"""
=== 정방향(46→42) 경로 ===
  L1 (dumB 방전) → L2 (LJ) → L3 (dumA 충전)
=== 역방향(42→46) 경로 — 같은 토폴로지를 거꾸로 훑는다 ===
  L3 역주행 (dumA 방전) → L2 역주행 (LJ) → L1 역주행 (dumB 충전)
  ★레그 순서 3→2→1, 각 레그 안에서 init-lambda-state 를 N-1 → 0
  ★ΔΔG(42→46) = −MBAR(공유 토폴로지의 A→B)""")
# 맨전하 부하 (접합부)
def bare(leg, at_end):
    """레그 경계에서 LJ 가 0 인데 전하가 남아있는 원자의 Σ|q|"""
    d={a['i']:a for a in L[leg]}
    s=0.0; n=0
    for i in dumA+dumB:
        a=d[i]; t = a['tB'] if at_end else a['tA']; q = a['qB'] if at_end else a['qA']
        if t.startswith('DUM') and abs(q)>1e-6: s+=abs(q); n+=1
    return s,n
print("\n=== 접합부 맨전하 부하 (LJ 0 인데 전하 남은 원자) ===")
for n in (1,2,3):
    s0,n0=bare(n,False); s1,n1=bare(n,True)
    print(f"  L{n}  시작 Σ|q|={s0:.4f}({n0}개)   끝 Σ|q|={s1:.4f}({n1}개)")
