# figures

`network_final_17edges.png` — 돌린 17엣지 그물(앵커 3 파랑 · 후보 7 주황 · 점선 = 수렴 플래그
N08·N11·N14)과 재설계 제안 2엣지. `../02_ttr_campaign/analysis/netdesign_v1.py` 산출물.
`prevalidation_roc_auc.png` · `tyk2_calc_vs_exp.png` — `make_figures.py` 가 저장소 안 CSV 에서
직접 읽어 그린다. 하드코딩한 값이 없고, 라벨은 한글 폰트가 없는 환경을 고려해 영문이다.

`overview.svg` — 캠페인 전체 흐름도. 이 저장소에서 **유일하게 손으로 작성한 그림**이고
`make_figures.py` 가 만들지 않는다. 그림 안의 수치는 최상위 `README.md` 본문에서 옮겨 적은
것이므로, **본문 수치가 바뀌면 이 파일도 같이 고쳐야 한다.** 라이트/다크 양쪽에 대응하지만
(`prefers-color-scheme`) OS 테마를 따르므로, GitHub 만 다크로 쓰는 환경에서는 밝게 보인다.
