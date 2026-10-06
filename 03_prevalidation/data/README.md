# data — 벤치마크 원표

`docking_scores.csv`(actives 53 + decoy 2,640 · 3D 생성 실패 10종 제외)와
`similarity_all_fp_v2.csv`(53 + 2,650)가 ROC-AUC 의 근거다. 두 파일의 decoy 수가 다른 이유는
그것뿐이다.
쌍 정확도는 시험이 셋으로 나뉜다 — `benzophenone_pairs.csv`(8쌍) ·
`diphenylether_pairs.csv`(21쌍) · `scaffold_pairs.csv`(43쌍, micro/macro 별도).
`prevalidation_summary.csv` 는 `../analysis/recompute_metrics.py` 가 만든 재계산 결과다.
