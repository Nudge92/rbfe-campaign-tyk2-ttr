#!/bin/bash
# =============================================================================
#  역방향 edge (ejm_42 → ejm_46) 준비
#
#  ★핵심: 하이브리드 토폴로지는 (2) 를 그대로 재사용한다. 새로 만들지 않는다.
#    pmx merged.itp 는 A(ejm_46) / B(ejm_42) 양쪽 파라미터를 모두 담고 있고,
#    검증에서 atoms 43/43 · bonds 45/45 · angles 78/78 · dihedrals 274/274 가
#    전부 B 상태 열을 가진 것을 확인했다.
#    (CDK2 self-RBFE 는 B 열 자체가 없어 무효였다 — 그 전례를 피하는 검사)
#
#  ★★시작 좌표는 반드시 mergedB.pdb 를 쓴다★★
#    mergedA.pdb = A 상태(ejm_46) 기하 — 정방향용
#    mergedB.pdb = B 상태(ejm_42) 기하 — 역방향용
#    검증: mergedA 앞 36원자가 ejm_46_ord.pdb 와 좌표차 0.000000 Å,
#          mergedB 의 매핑 28쌍이 ejm_42_ord.pdb 와 좌표차 0.000000 Å.
#    정방향 좌표로 역방향을 시작하면 λ=0 끝점이 B 상태인데 기하는 A 라
#    첫 창부터 큰 왜곡이 걸린다.
#
#  ★λ 벡터는 뒤집는다
#    정방향: fep/coul/vdw-lambdas = 0.0 … 1.0   (init-lambda-state 0 → N)
#    역방향: 같은 벡터를 역순으로 쓰거나, 같은 벡터에 대해
#            init-lambda-state 를 N → 0 으로 훑는다.
#    둘 중 어느 쪽이든 토폴로지는 동일하다.
#
#  이 스크립트는 파일 배치만 한다. 계 구축과 MD 는 하지 않는다.
# =============================================================================
set -euo pipefail

B=/home/nudge/Project/TYK2-PILOT
SRC=$B/04_hybrid/2_ejm46_ejm42
DST=$B/04_hybrid/2R_ejm42_ejm46

[ -f "$SRC/merged.itp" ] || { echo "없음: $SRC/merged.itp" >&2; exit 1; }
mkdir -p "$DST"

# 토폴로지는 심볼릭 링크로 공유 — 복사본이 갈라지는 것을 막는다
for f in merged.itp ffmerged.itp hybrid_atomtypes.itp; do
    ln -sf "$SRC/$f" "$DST/$f"
done
# 시작 좌표만 B 상태
cp "$SRC/mergedB.pdb" "$DST/start.pdb"
cp "$SRC/mergedA.pdb" "$DST/end_reference.pdb"

cat > "$DST/README.txt" <<'T'
역방향 edge: ejm_42 → ejm_46

토폴로지  : ../2_ejm46_ejm42/merged.itp (심볼릭 링크, 공유)
시작 좌표 : start.pdb  = 2_ejm46_ejm42/mergedB.pdb  ← ★반드시 이것
참고      : end_reference.pdb = mergedA.pdb (정방향 시작점)

주의
  - merged.itp 의 A 상태는 ejm_46, B 상태는 ejm_42 다.
    역방향은 "B 에서 출발해 A 로 간다" 이므로 λ 진행을 뒤집는다.
  - 토폴로지를 새로 만들지 않는다. 새로 만들면 매핑이 달라져
    정·역 비교(hysteresis)가 성립하지 않는다.
  - 정·역이 같은 매핑·같은 파라미터를 쓰는 것이 이 대조의 전제다.
T

echo "역방향 준비 완료: $DST"
ls -la "$DST" | sed 's/^/  /'
echo
echo "== B 상태 완전성 재확인 =="
"${PY:-$HOME/miniforge3/envs/md/bin/python}" - "$SRC/merged.itp" <<'PY'
import sys, re
p=sys.argv[1]
def sect(name, nidx):
    rows=[]; on=False
    for l in open(p):
        s=l.strip()
        if s.startswith('['): on=re.match(rf'\[\s*{name}\s*\]',s) is not None; continue
        if not on or not s or s.startswith(';'): continue
        f=re.sub(r';.*','',s).split()
        if f: rows.append(f)
    withB=sum(1 for f in rows if (len(f)-nidx-1)>0 and (len(f)-nidx-1)%2==0)
    return withB, len(rows)
for n,i in (('atoms',None),('bonds',2),('angles',3),('dihedrals',4)):
    if n=='atoms':
        rows=[];on=False
        for l in open(p):
            s=l.strip()
            if s.startswith('['): on='atoms' in s; continue
            if on and s and not s.startswith(';'):
                f=s.split()
                if len(f)>=5 and f[0].isdigit(): rows.append(f)
        w=sum(1 for f in rows if len(f)>=11)
        print(f"  atoms      {w}/{len(rows)}  (typeB·chargeB·massB 존재)")
    else:
        w,t=sect(n,i); print(f"  {n:<10} {w}/{t}")
PY
