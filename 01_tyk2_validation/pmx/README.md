# pmx — 하이브리드 토폴로지 생성과 검증

`build_hybrids.sh` 가 `pmx atomMapping --H2Hpolar --d 0.05` 와 `pmx ligandHybrid -pairs` 를
부르고, 옵션은 파일 머리말에 전부 적혀 있다. `verify_*.py` 는 매핑·방향·역방향 토폴로지를
자동으로 대조한다(더미 중원자 수 대조 포함).

경로는 실행 당시 그대로이며 재현하려면 수정이 필요하다.
