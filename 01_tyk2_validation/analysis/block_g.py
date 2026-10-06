#!/usr/bin/env python
"""
블록 평균법으로 '평균의 실효 상관시간' g_eff 를 구한다.

원시 시계열의 자기상관 g 는 신호가 (큰 빠른 성분 + 작은 느린 성분) 으로
되어 있을 때 큰 쪽만 본다. 평균의 불확도를 정하는 것은 느린 쪽이다.

  g_eff(b) = b · Var(블록평균) / Var(원시)      [b = 블록 길이, 프레임]

b 를 키우며 g_eff(b) 가 평평해지는 구간(plateau)의 값이 실효 g 다.
★블록 수가 적으면 Var(블록평균) 추정이 불안정하다 —
  상대오차 ≈ sqrt(2/(n_b−1)). n_b >= 20 을 유지한다.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
★★ HREMD 를 쓸 때 — 이 파일의 핵심 주의 (2026-09-08 추가)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
GROMACS 복제교환은 **좌표를 교환**하고 각 디렉터리는 **고정 λ** 를 유지한다.
디렉터리별 dH/dλ 시계열은 MBAR 입력으로는 옳지만, **상관시간 계산에는 쓰면 안 된다.**
교환이 일어날 때마다 그 디렉터리의 물리적 계가 통째로 바뀌므로 시계열이
인위적으로 탈상관되어 **g 를 크게 과소평가**하고, 그만큼 Neff·오차막대가 낙관적이 된다.

  올바른 g  = 워커(연속 궤적)를 따라간 시계열에서 계산한 것
  틀린 g    = 디렉터리 고정(state-fixed) 시계열에서 계산한 것

→ HREMD 데이터에는 `g_eff_hremd(...)` 를 쓴다. `g_eff(...)` 를 그대로 쓰면
   `hremd_log=` 를 넘겨 경고를 받거나, 아래 가드에 걸려 멈춘다.
"""
import re
import numpy as np


# ─────────────────────────── 기본 블록 평균 ───────────────────────────

def g_eff_curve(a, dt=0.2, min_blocks=20):
    a = np.asarray(a, float); N = len(a); v = a.var(ddof=1)
    out = []
    b = 2
    while N//b >= min_blocks:
        nb = N//b
        bm = a[:nb*b].reshape(nb, b).mean(axis=1)
        g = b*bm.var(ddof=1)/v
        rel = np.sqrt(2.0/(nb-1))
        out.append((b*dt, g*dt, nb, rel))
        b = int(np.ceil(b*1.4))
    return out


def plateau(curve, ntail=3):
    """plateau 값 = 꼬리 ntail 점의 중앙값. 곡선이 짧으면 nan.

    ★ 마지막 한 점(`curve[-1][1]`)을 쓰면 블록 수가 가장 적은 지점이라
      상대오차가 가장 큰 값을 그대로 채택하게 된다. 항상 이 함수를 쓸 것."""
    if len(curve) < ntail:
        return (float('nan'), float('nan'))
    tail = [x[1] for x in curve[-ntail:]]
    return (float(np.median(tail)), float(np.std(tail)))


def g_eff(a, dt=0.2, min_blocks=20, hremd_log=None):
    """plateau 값 = 마지막 3점의 중앙값 (블록 수 조건 만족 구간에서).

    ★ hremd_log 를 주면 '이 시계열은 state-fixed 이므로 g 가 과소평가된다'는
      뜻으로 해석해 예외를 던진다. HREMD 는 g_eff_hremd 를 쓸 것."""
    if hremd_log is not None:
        raise ValueError(
            "HREMD 데이터에 g_eff() 를 쓰면 상관시간을 과소평가한다. "
            "g_eff_hremd(D, logpath, dt) 를 사용할 것.")
    c = g_eff_curve(a, dt, min_blocks)
    m, s = plateau(c)
    return (m, s, c)


# ─────────────────────────── HREMD demux ───────────────────────────

def parse_exchange(logpath):
    """복제교환 로그 → (교환 시각[ps], 순열 문자열) 목록."""
    times = []; perms = []; cur = None
    for ln in open(logpath, errors='ignore'):
        m = re.match(r'Replica exchange at step (\d+) time ([\d.]+)', ln)
        if m:
            cur = float(m.group(2)); continue
        if cur is not None and ln.startswith('Repl ex'):
            times.append(cur); perms.append(ln[7:].rstrip('\n')); cur = None
    return times, perms


def build_walkers(nrep, perms):
    """워커 w 가 각 교환 시점에 어느 디렉터리(=고정 λ 상태)에 있는지.
    반환 (n_exchange+1, nrep) 정수 배열."""
    state = list(range(nrep))          # state[w] = 디렉터리 인덱스
    traj = [list(state)]
    for tk in perms:
        for m in re.finditer(r'(\d+)\s+x\s+(\d+)', tk):
            a, b = int(m.group(1)), int(m.group(2))
            wa = state.index(a); wb = state.index(b)
            state[wa], state[wb] = state[wb], state[wa]
        traj.append(list(state))
    return np.array(traj)


def acceptance(perms, nrep):
    att = 0; acc = 0
    for tk in perms:
        att += len(re.findall(r'\d+\s+[x ]\s+\d+', tk)) or (nrep-1)
        acc += len(re.findall(r'\d+\s+x\s+\d+', tk))
    return acc, att, (acc/att if att else float('nan'))


def demux_series(D, logpath, dt=0.2):
    """state-fixed 시계열 D (nrep, nframe) → 워커를 따라간 시계열 (nrep, nframe).

    ★ 이것이 상관시간 계산의 올바른 입력이다."""
    D = np.asarray(D, float); nrep, n = D.shape
    times, perms = parse_exchange(logpath)
    if not perms:
        raise ValueError(f"{logpath}: 교환 기록 없음 — HREMD 가 아니거나 로그가 잘렸다")
    tr = build_walkers(nrep, perms)
    fr_t = np.arange(n)*dt
    idx = np.clip(np.searchsorted(np.array([0.0]+times), fr_t, side='right')-1, 0, len(tr)-1)
    return np.array([D[tr[idx, w], np.arange(n)] for w in range(nrep)]), tr, idx, (times, perms)


def g_eff_hremd(D, logpath, dt=0.2, min_blocks=20):
    """HREMD 의 올바른 g. state-fixed g 도 함께 돌려주어 배율을 볼 수 있게 한다.

    반환 dict: g_state / g_walk / g_mix / ratio / rate / nrep / nframe"""
    D = np.asarray(D, float); nrep, n = D.shape
    W, tr, idx, (times, perms) = demux_series(D, logpath, dt)
    acc, att, rate = acceptance(perms, nrep)
    g_state = [plateau(g_eff_curve(D[i], dt, min_blocks))[0] for i in range(nrep)]
    g_walk  = [plateau(g_eff_curve(W[i], dt, min_blocks))[0] for i in range(nrep)]
    g_mix   = [plateau(g_eff_curve(tr[idx, w].astype(float), dt, min_blocks))[0] for w in range(nrep)]
    gs = np.array(g_state, float); gw = np.array(g_walk, float)
    return dict(g_state=g_state, g_walk=g_walk, g_mix=g_mix,
                ratio=float(np.nanmean(gw)/np.nanmean(gs)),
                ratio_per_window=(gw/gs).tolist(),
                rate=rate, acc=acc, att=att, nex=len(perms), nrep=nrep, nframe=n)


def neff(n, g):
    """실효 표본 수. g 는 시간 단위(ps)이고 n 은 프레임 수이므로 dt 로 나눈 값을 쓴다."""
    return n/max(g, 1e-12)
