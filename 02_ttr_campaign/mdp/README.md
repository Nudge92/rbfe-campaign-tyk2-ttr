# mdp — 두 프로토콜의 실행 설정

`net16_*` 는 1차 그물 16엣지(GROMACS 2025.2 · 창당 0.5 ns · `-replex 200`),
`N17_*` 는 추가 엣지(2026.3 · 창당 8 ns · 창 8/17/10)의 복합체 레그 1·2·3 프로덕션 설정이다.
soft-core 는 두 프로토콜과 TYK2 가 모두 같다(`sc-alpha 0.5` · `sc-power 1` · `sc-sigma 0.3`).
`net16_{em,nvt}_restrained_k500.mdp` 는 단백질 구속(k=500) 완화 단계다.
