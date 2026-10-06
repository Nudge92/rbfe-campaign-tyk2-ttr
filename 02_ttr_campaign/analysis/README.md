# analysis — 추출 · 고리 분석 · 재설계

`net_extract_v1.py` 가 원 로그·JSON 에서 엣지별 값을 전수 추출하고(`--selftest` 내장),
`verify_cycle_signs_v11.py` 는 고리 부호를 그래프 순회로 독립 재유도한다.
`node_potential_v11.py` 는 최소제곱으로 노드 전위를, `netdesign_v1.py` 는 재설계안을 낸다.
`xw1_score_v1.py` 는 결정수 포함 조건의 재채점이다.

경로는 실행 당시 그대로이며 재현하려면 수정이 필요하다.
