#!/usr/bin/env python3
"""
TYK2 3단계 λ 스케줄 — coul/vdw 순차 분리.

  1구간: coul 0 -> 1 (vdw 는 A 상태에 고정)   ← 전하를 먼저 바꾼다
  2구간: coul 1 고정, vdw 0 -> 1              ← 그다음 LJ
  bonded·fep = vdw 와 동일 (결합길이 변화는 원자 정체 변화와 함께 가야 안전)
  접합부(coul=1, vdw=0)는 두 구간이 공유 → 총 창 = N_coul + N_vdw - 1

CDK2 26창은 네 벡터가 전부 동일한 동시 소거였고 그것이 L00 문제의 배경이었다.
여기서는 처음부터 순차로 간다. sc-coul = no — coul 을 vdw 가 A 상태인 동안
바꾸므로 전하 소프트코어가 필요 없다.

★창 개수는 잠정치다. 짧은 벤치마크에서 인접 ΔH 를 재고 확정한다.
"""
import argparse

def uniform(n):
    return [round(i/(n-1), 4) for i in range(n)]

def endpoint_dense(n):
    """양 끝 조밀 — vdW 소프트코어 구간용"""
    if n < 4: return uniform(n)
    inner = [i/(n-3) for i in range(n-2)]
    pts = [0.0, inner[1]/2] + inner[1:-1] + [1-(1-inner[-2])/2, 1.0]
    return [round(x, 4) for x in pts]

def build(ncoul, nvdw, spacing='endpoint', order='vdw_first'):
    cpts = uniform(ncoul)
    vpts = endpoint_dense(nvdw) if spacing == 'endpoint' else uniform(nvdw)
    if order == 'coul_first':
        coul = list(cpts) + [1.0]*(nvdw-1)
        vdw  = [0.0]*ncoul + list(vpts[1:])
    else:
        # ★vdw 먼저. coul-first 는 접합부(coul=1,vdw=0)에서 생성 원자가
        #   LJ 없이 전하만 갖게 되어 붕괴한다(실측: 5000 step 사망).
        #   vdw-first 접합부(coul=0,vdw=1)는 1 ns 완주 확인.
        vdw  = list(vpts) + [1.0]*(ncoul-1)
        coul = [0.0]*nvdw + list(cpts[1:])
    assert len(coul) == len(vdw) == ncoul+nvdw-1
    return coul, vdw

def fmt(name, v):
    return f"{name:<24} = " + " ".join(f"{x:g}" for x in v)

def emit(coul, vdw, ncoul, nvdw, state='__STATE__'):
    L = [
     "; ===== TYK2 3단계 RBFE · λ 순차 분리 (vdw -> coul) =====",
     f"; gen_lambdas.py --ncoul {ncoul} --nvdw {nvdw}"
     f" · 총 {len(coul)} 창 (vdw {nvdw} + coul {ncoul} - 접합부 1)",
     "free-energy              = yes",
     f"init-lambda-state        = {state}",
     fmt("fep-lambdas", vdw),
     fmt("coul-lambdas", coul),
     fmt("vdw-lambdas", vdw),
     fmt("bonded-lambdas", vdw),
     "sc-alpha                 = 0.5",
     "sc-power                 = 1",
     "sc-sigma                 = 0.3",
     "sc-coul                  = no",
     "nstdhdl                  = 100",
     "calc-lambda-neighbors    = -1",      # MBAR: 전 상태 ΔH 저장
     "dhdl-derivatives         = yes",
     "separate-dhdl-file       = yes",
    ]
    return "\n".join(L)+"\n"

def trap_w(l):
    n=len(l); w=[]
    for i in range(n):
        lo = l[i]-l[i-1] if i > 0 else 0.0
        hi = l[i+1]-l[i] if i < n-1 else 0.0
        w.append((lo+hi)/2)
    return w

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--ncoul', type=int, required=True)
    p.add_argument('--nvdw', type=int, required=True)
    p.add_argument('--spacing', choices=('endpoint','uniform'), default='endpoint')
    p.add_argument('--order', choices=('vdw_first','coul_first'), default='vdw_first')
    p.add_argument('--weights', action='store_true')
    a = p.parse_args()
    coul, vdw = build(a.ncoul, a.nvdw, a.spacing, a.order)
    if a.weights:
        wc, wv = trap_w(coul), trap_w(vdw)
        print(f"{'state':>5}{'coul':>9}{'vdw':>9}{'w_coul':>9}{'w_vdw':>9}")
        for i,(c,v) in enumerate(zip(coul,vdw)):
            print(f"{i:>5}{c:>9.4f}{v:>9.4f}{wc[i]:>9.4f}{wv[i]:>9.4f}")
        print(f"\nΣw_coul={sum(wc):.6f}  Σw_vdw={sum(wv):.6f}  총 {len(coul)} 창")
    else:
        print(emit(coul, vdw, a.ncoul, a.nvdw), end='')
