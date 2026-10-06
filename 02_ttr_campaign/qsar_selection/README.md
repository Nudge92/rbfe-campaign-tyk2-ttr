# qsar_selection — FEP 후보를 고른 축과 그 축의 전이 실패

`model_loso.py` 가 모델 9종(A1~C3)을 사전 선언한다 — 기술자 3안 × 계열이지만 완전교차가
아니다(2열이 A는 GP, B·C는 kNN). `dock_targets.py:12,41-42` 가 대표 3종(A1 Ridge · B1 KRR ·
A3 RF)의 Borda 순위로 `cand####` 이름을 붙인다 — A1·A3 이 같은 기술자(물성20)를 쓰므로
기술자는 2종이다. `D_consensus_rep3.md` 는 다른 분석의 대표 3종(A1·B1·C3)이라 혼동 주의.

`tox24_external.py` + `D_external.json` 이 외부 검증 0/9 의 근거다(평가 64종 · InChIKey 겹침 0).
`A3_summary_endpoint_switch.md` 는 엔드포인트만 바꿨을 때의 0/9 → 2/9 비교다.

경로는 실행 당시 그대로이며 재현하려면 수정이 필요하다.
