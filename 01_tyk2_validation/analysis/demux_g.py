#!/usr/bin/env python
"""HREMD 상관시간 — 워커 연속 궤적 기준 g 대 상태 고정 기준 g.

★GROMACS 복제교환은 **좌표를 교환**하고 각 디렉터리는 **고정 λ** 를 유지한다.
  → MBAR 입력(디렉터리별 dhdl)은 그대로 유효.
  → 상관시간만 워커를 따라 다시 계산해야 한다.
"""
import numpy as np, re, os, sys
from block_g import g_eff_curve

def parse_exchange(logpath):
    """복제교환 로그에서 시간별 순열을 읽는다.
       'Replica exchange at step N time T' 다음 줄들의 'Repl ex' 패턴."""
    times=[]; perms=[]
    cur=None
    for ln in open(logpath, errors='ignore'):
        m=re.match(r'Replica exchange at step (\d+) time ([\d.]+)', ln)
        if m:
            cur=float(m.group(2)); continue
        if cur is not None and ln.startswith('Repl ex'):
            # 'Repl ex  0 x  1    2 x  3 ...'  x=교환 성사, 빈칸=실패
            toks=ln[7:].rstrip('\n')
            times.append(cur); perms.append(toks); cur=None
    return times, perms

def build_walkers(nrep, times, perms):
    """워커 w 가 시각 t 에 어느 디렉터리(=상태)에 있는지."""
    state=list(range(nrep))          # state[w] = 디렉터리 인덱스
    traj=[list(state)]
    for tk in perms:
        # 'i x j' 쌍 추출
        for m in re.finditer(r'(\d+)\s+x\s+(\d+)', tk):
            a,b=int(m.group(1)),int(m.group(2))
            wa=state.index(a); wb=state.index(b)
            state[wa],state[wb]=state[wb],state[wa]
        traj.append(list(state))
    return np.array(traj)             # (n_exchange+1, nrep)

def acceptance(perms, nrep):
    att=0; acc=0
    for tk in perms:
        att += len(re.findall(r'\d+\s+[x ]\s+\d+', tk)) or (nrep-1)
        acc += len(re.findall(r'\d+\s+x\s+\d+', tk))
    return acc, att, (acc/att if att else float('nan'))

def analyze(rd, nrep, dt=0.2):
    """rd/wXX/prod.xvg + rd/w00/prod.log"""
    from alchemlyb.parsing.gmx import extract_dHdl
    D=[]
    for i in range(nrep):
        d=extract_dHdl(f'{rd}/w{i:02d}/prod.xvg', T=300)
        c=[x for x in d.columns]
        a=np.asarray(d[c].sum(axis=1), float)
        D.append(a)
    n=min(len(x) for x in D); D=np.array([x[:n] for x in D])   # (nrep, nframe)
    # (a) 상태 고정 기준 — 기존 방식
    g_state=[g_eff_curve(D[i])[-1][1] for i in range(nrep)]
    # (b) 워커 기준
    times,perms=parse_exchange(f'{rd}/w00/prod.log')
    acc,att,rate=acceptance(perms,nrep)
    if not perms: return dict(g_state=g_state, rate=float('nan'), g_walk=None, note='교환 기록 없음')
    tr=build_walkers(nrep,times,perms)          # (nex+1, nrep)
    # 각 프레임을 가장 가까운 교환 시각에 매핑
    fr_t=np.arange(n)*dt
    idx=np.searchsorted(np.array([0]+times), fr_t, side='right')-1
    idx=np.clip(idx,0,len(tr)-1)
    g_walk=[]
    for w in range(nrep):
        st=tr[idx,w]                            # 프레임별 이 워커가 있는 디렉터리
        series=D[st, np.arange(n)]              # 워커를 따라간 dH/dl
        g_walk.append(g_eff_curve(series)[-1][1])
    # 상태 지수 혼합시간
    g_mix=[g_eff_curve(tr[idx,w].astype(float))[-1][1] for w in range(nrep)]
    return dict(g_state=g_state, g_walk=g_walk, g_mix=g_mix, rate=rate, acc=acc, att=att,
                nex=len(perms))

if __name__=='__main__':
    rd=sys.argv[1]; nrep=int(sys.argv[2])
    r=analyze(rd,nrep)
    print(f"교환 시도 {r.get('att','?')} · 성사 {r.get('acc','?')} · **수락률 {r['rate']*100:.1f}%**")
    if r.get('g_walk'):
        gs=np.array(r['g_state']); gw=np.array(r['g_walk']); gm=np.array(r['g_mix'])
        print(f"\n{'창':<6}{'g_state(기존)':>14}{'g_walker(demux)':>17}{'배율':>8}{'g_mix(상태지수)':>16}")
        for i in range(nrep):
            print(f"w{i:<5}{gs[i]:14.2f}{gw[i]:17.2f}{gw[i]/gs[i]:8.2f}{gm[i]:16.2f}")
        print(f"\n  평균 g_state {gs.mean():.2f} · g_walker {gw.mean():.2f} → "
              f"★기존 방식이 상관시간을 **{gw.mean()/gs.mean():.2f}배 과소평가**")
        print(f"  = Neff 를 {gw.mean()/gs.mean():.2f}배 과대평가")
