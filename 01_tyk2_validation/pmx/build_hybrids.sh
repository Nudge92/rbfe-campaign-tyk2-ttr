#!/bin/bash
# =============================================================================
#  TYK2-PILOT 1단계 — pmx 하이브리드 생성 (4벌)
#
#  CDK2-FEP/03_hybrid/build.sh 에서 복원한 명령·옵션을 그대로 쓴다:
#    atomMapping  : --H2Hpolar · alignment+MCS(기본) · chirality check(기본)
#                   --RingsOnly 미사용 · --swap 미사용 · --d 0.05 · --timeout 10
#    ligandHybrid : -pairs 명시 · --fit 미사용(CDK2 에서 바이트 일치로 확정)
#
#  순서대로 진행하고 앞의 것이 통과해야 다음으로 간다.
#    (1) ejm_31 + ejm_42   LOMAP 0.905  대조군
#    (2) ejm_46 + ejm_42   LOMAP 0.905  ★정·역 양방향 공유
#    (3) ejm_46 + ejm_50   LOMAP 0.333
#    (4) ejm_46 + ejm_47   LOMAP 0.050  ★3원환→4원환, 마지막
#
#  ★(4) 주의: 사이클로프로필(cx 3개) → 사이클로부틸(cy 4개) 는 삼각형과
#    사각형의 원자 대응에 유일한 답이 없다. 매핑이 화학적으로 말이 안 되면
#    거기서 멈춘다. (1)~(3) 만으로도 2단계에 갈 수 있다.
# =============================================================================
set -uo pipefail

B=/home/nudge/Project/TYK2-PILOT
PMX=${PMX:-$HOME/miniforge3/envs/md/bin/pmx}
GMX=${GMX:-$HOME/miniforge3/envs/md_cuda/bin/gmx}
PY=${PY:-$HOME/miniforge3/envs/md/bin/python}
P=$B/03_params
OUT=$B/04_hybrid

for exe in "$PMX" "$GMX"; do [ -x "$exe" ] || { echo "없음: $exe" >&2; exit 1; }; done
mkdir -p "$OUT"

# ---------------------------------------------------------------- 0단계
# GROMACS 원자순서 PDB. pmx 는 원자 인덱스로 pairs 를 쓰므로
# lig.itp 와 인덱스가 1:1 이어야 한다 → lig_gmx.gro 를 그대로 변환.
echo "== 0단계: 원자순서 정렬 PDB"
for L in ejm_31 ejm_42 ejm_46 ejm_47 ejm_50; do
    [ -f "$OUT/${L}_ord.pdb" ] && continue
    "$GMX" editconf -f "$P/$L/lig_gmx.gro" -o "$OUT/${L}_ord.pdb" >/dev/null 2>&1
    echo "   ${L}_ord.pdb  $(grep -c '^ATOM\|^HETATM' "$OUT/${L}_ord.pdb") atoms"
done

# ---------------------------------------------------------------- 하이브리드
build_one () {
    local n=$1 A=$2 Bl=$3 note=$4
    local d=$OUT/$n
    echo
    echo "════════ ($n) $A + $Bl   $note ════════"
    mkdir -p "$d"; cd "$d"

    # --- 매핑 ---
    # --H2Hpolar : 극성 수소끼리 매핑 허용. 빼면 아마이드 N-H 가 dummy 가 되어
    #              매핑이 줄고 하이브리드가 불필요하게 커진다 (CDK2 v1 전례).
    "$PMX" atomMapping \
        -i1 "$OUT/${A}_ord.pdb" -i2 "$OUT/${Bl}_ord.pdb" \
        -o1 pairs1.dat -o2 pairs2.dat \
        -opdbm1 m1.pdb -opdbm2 m2.pdb \
        -score score.dat -log map.log \
        --H2Hpolar > pmx_map.out 2>&1 \
        || { echo "   ★atomMapping 실패"; tail -15 pmx_map.out; return 1; }

    local np; np=$(wc -l < pairs1.dat)
    echo "   매핑 $np 쌍 · $(cat score.dat)"

    # --- 하이브리드 ---
    "$PMX" ligandHybrid \
        -i1 "$OUT/${A}_ord.pdb" -i2 "$OUT/${Bl}_ord.pdb" \
        -itp1 "$P/$A/lig.itp" -itp2 "$P/$Bl/lig.itp" \
        -pairs pairs1.dat \
        -oA mergedA.pdb -oB mergedB.pdb \
        -oitp merged.itp -offitp ffmerged.itp \
        -log hyb.log > pmx_hyb.out 2>&1 \
        || { echo "   ★ligandHybrid 실패"; tail -15 pmx_hyb.out; return 1; }

    # --- atomtypes 합집합 (pmx 산출물 아님) ---
    # ★정렬하지 않는다 — 등장 순서 보존 (CDK2 에서 sort 때문에 diff 가 났다)
    {
        echo "[ atomtypes ]"
        cat "$P/$A/lig_atomtypes.itp" "$P/$Bl/lig_atomtypes.itp" \
            | grep -v '^\[' | grep -v '^[[:space:]]*;' | grep -v '^[[:space:]]*$' \
            | awk '!seen[$1]++'
    } > hybrid_atomtypes.itp

    echo "   → merged.itp $(awk '/^\[ *atoms *\]/{a=1;next} /^\[/{a=0} a&&$0!~/^;/&&NF>3{c++} END{print c}' merged.itp) atoms"
    cd "$OUT"; return 0
}

FAIL=0
build_one 1_ejm31_ejm42 ejm_31 ejm_42 "LOMAP 0.905 대조군"      || FAIL=1
[ $FAIL -eq 0 ] && { build_one 2_ejm46_ejm42 ejm_46 ejm_42 "LOMAP 0.905 정·역 공유" || FAIL=2; }
[ $FAIL -eq 0 ] && { build_one 3_ejm46_ejm50 ejm_46 ejm_50 "LOMAP 0.333"            || FAIL=3; }
[ $FAIL -eq 0 ] && { build_one 4_ejm46_ejm47 ejm_46 ejm_47 "LOMAP 0.050 3원환→4원환" || FAIL=4; }

echo
if [ $FAIL -eq 0 ]; then echo "4벌 전부 생성 완료"
else echo "★($FAIL) 에서 중단 — 앞 단계 산출물은 유효"; fi
echo "다음: verify_hybrids.py"
exit 0
