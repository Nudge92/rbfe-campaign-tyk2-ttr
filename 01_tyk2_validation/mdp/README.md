# mdp — λ 스케줄 · soft-core · 교환 설정

`fep_block.txt` 가 24창 3블록 λ 사다리(전하 끄기 → vdW·결합 → 전하 켜기)와 soft-core
(`sc-alpha 0.5` · `sc-power 1` · `sc-sigma 0.3` · `sc-coul no`)를 정의한다.
`win_prod.mdp` 는 창당 1 ns 프로덕션, 교환은 실행 시 `-replex 1000`(2 ps)으로 준다.
`hremd_legs/` 는 레그별 실제 창 수(6/11/7)를 담은 세 블록이다.
