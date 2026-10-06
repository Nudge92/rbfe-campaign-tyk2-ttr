#!/usr/bin/env python3
"""
3블록 FEP 토폴로지 생성기 — merged.itp 한 벌에서 레그 3벌을 만든다.

논문 처방 (LiveCoMS Best Practices):
  1. 사라지는 원자의 정전기 끄기
  2. Lennard-Jones 변경 (삽입·제거 동시)
  3. 도입되는 원자의 정전기 켜기

GROMACS 는 원자별 λ 그룹을 지원하지 않으므로 단일 사다리로는 불가능하다.
토폴로지 3벌을 만들어 레그 3개로 나눈다.

  q_mid = (공통 원자: 원본 qB · 소멸/생성 원자: 0)

  레그1 방전: typeA=typeB=원본A · qA=원본qA → qB=q_mid · bonded B←A
  레그2 LJ  : typeA=원본A typeB=원본B · qA=qB=q_mid · bonded 원본 그대로
  레그3 충전: typeA=typeB=원본B · qA=q_mid → qB=원본qB · bonded A←B

이렇게 하면 어느 상태에서도 'LJ 없이 전하만 가진 원자'가 없다.
"""
import sys, os, re

NIDX = {'bonds': 2, 'pairs': 2, 'angles': 3, 'dihedrals': 4}


def parse(path):
    """섹션별 (원본줄, 파싱된 필드) 보존"""
    out = []          # [(sec, raw, fields|None)]
    sec = None
    for ln in open(path):
        s = ln.rstrip('\n')
        body = s.split(';')[0]
        m = re.match(r'\s*\[\s*(\S+)\s*\]', body)
        if m:
            sec = m.group(1); out.append((sec, s, None)); continue
        if not body.strip():
            out.append((sec, s, None)); continue
        f = body.split()
        out.append((sec, s, f if f[0].isdigit() else None))
    return out


def build(src, outdir, tag):
    rows = parse(src)
    # ---- 1) 원본 atoms 읽기 ----
    atoms = {}
    for sec, raw, f in rows:
        if sec == 'atoms' and f and len(f) >= 11:
            atoms[int(f[0])] = dict(tA=f[1], qA=float(f[6]), mA=float(f[7]),
                                    tB=f[8], qB=float(f[9]), mB=float(f[10]))
    # ---- 2) q_mid ----
    qmid = {}
    for i, a in atoms.items():
        uniq = a['tA'].startswith('DUM') or a['tB'].startswith('DUM')
        qmid[i] = 0.0 if uniq else a['qB']

    def emit(leg):
        lines = []
        for sec, raw, f in rows:
            if sec == 'atoms' and f and len(f) >= 11:
                i = int(f[0]); a = atoms[i]
                if leg == 1:
                    tA, tB = a['tA'], a['tA']; qa, qb = a['qA'], qmid[i]; mA = mB = a['mA']
                elif leg == 2:
                    tA, tB = a['tA'], a['tB']; qa, qb = qmid[i], qmid[i]; mA, mB = a['mA'], a['mB']
                else:
                    tA, tB = a['tB'], a['tB']; qa, qb = qmid[i], a['qB']; mA = mB = a['mB']
                lines.append(f"{i:6d} {tA:>10s} {f[2]:>6s} {f[3]:>6s} {f[4]:>6s} {f[5]:>6s}"
                             f" {qa:12.6f} {mA:10.4f} {tB:>10s} {qb:12.6f} {mB:10.4f}")
            elif sec in NIDX and f and leg != 2:
                n = NIDX[sec]; rest = f[n+1:]
                if not rest or len(rest) % 2:
                    lines.append(raw); continue
                h = len(rest)//2
                # 레그1: B←A · 레그3: A←B  (결합 섭동 없음)
                keep = rest[:h] if leg == 1 else rest[h:]
                lines.append(" ".join(f[:n+1]) + " " + " ".join(keep) + " " + " ".join(keep))
            else:
                lines.append(raw)
        return "\n".join(lines) + "\n"

    os.makedirs(outdir, exist_ok=True)
    for leg in (1, 2, 3):
        open(f"{outdir}/{tag}_leg{leg}.itp", 'w').write(emit(leg))
    return atoms, qmid


def verify(src, outdir, tag):
    """레그 경계 연속성 + 끝점이 원본 A/B 와 일치하는지"""
    orig = {}
    for sec, raw, f in parse(src):
        if sec == 'atoms' and f and len(f) >= 11:
            orig[int(f[0])] = (f[1], float(f[6]), f[8], float(f[9]))
    L = {}
    for leg in (1, 2, 3):
        d = {}
        for sec, raw, f in parse(f"{outdir}/{tag}_leg{leg}.itp"):
            if sec == 'atoms' and f and len(f) >= 11:
                d[int(f[0])] = (f[1], float(f[6]), f[8], float(f[9]))
        L[leg] = d
    ok = True; msg = []
    n = len(orig)
    # 시작점 = 원본 A
    bad = [i for i in orig if L[1][i][0] != orig[i][0] or abs(L[1][i][1]-orig[i][1]) > 1e-6]
    if bad: ok = False; msg.append(f"레그1 시작 ≠ 원본A ({len(bad)})")
    # 끝점 = 원본 B
    bad = [i for i in orig if L[3][i][2] != orig[i][2] or abs(L[3][i][3]-orig[i][3]) > 1e-6]
    if bad: ok = False; msg.append(f"레그3 끝 ≠ 원본B ({len(bad)})")
    # 경계 연속 (레그1끝 == 레그2시작, 레그2끝 == 레그3시작)
    for a, b in ((1, 2), (2, 3)):
        bad = [i for i in orig if L[a][i][2] != L[b][i][0] or abs(L[a][i][3]-L[b][i][1]) > 1e-6]
        if bad: ok = False; msg.append(f"레그{a}끝 ≠ 레그{b}시작 ({len(bad)})")
    # 순전하
    nets = [sum(L[1][i][1] for i in orig), sum(L[1][i][3] for i in orig),
            sum(L[3][i][3] for i in orig)]
    return ok, msg, nets, n


if __name__ == '__main__':
    B = '/home/nudge/Project/TYK2-PILOT'
    import sys
    E = {'edge1': '1_ejm31_ejm42', 'edge2': '2_ejm46_ejm42',
         'edge3': '3_ejm46_ejm50', 'edge4': '4_ejm46_ejm47',
         'edge5': '5_ejm31_ejm50', 'edge6': '6_ejm42_ejm50'}
    if len(sys.argv) > 1: E = {t: E[t] for t in sys.argv[1:]}
    print(f"{'edge':<8}{'원자':>5}{'검증':>8}{'net(A)':>10}{'net(중간)':>11}{'net(B)':>10}  비고")
    for tag, d in E.items():
        src = f'{B}/04_hybrid/{d}/merged.itp'
        out = f'{B}/07_step3_lambda/3block/{tag}'
        build(src, out, tag)
        ok, msg, nets, n = verify(src, out, tag)
        print(f"{tag:<8}{n:>5}{'통과' if ok else '★실패':>8}"
              f"{nets[0]:>10.4f}{nets[1]:>11.4f}{nets[2]:>10.4f}  {' '.join(msg)}")
