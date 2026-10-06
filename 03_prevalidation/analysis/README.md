# analysis — 재계산

`recompute_metrics.py` 는 `../data/` 의 CSV 만 읽어 ROC-AUC(Mann-Whitney U · 동점 0.5)와
쌍 정확도를 다시 계산한다. 하드코딩한 수치가 없으므로, README 의 숫자가 원표에서 나오는지
이 스크립트 하나로 확인할 수 있다. 표준 라이브러리만 쓴다.
