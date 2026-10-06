# analysis — 자유에너지 추정과 신뢰구간

`analyze*.py` 가 alchemlyb/pymbar 로 MBAR·BAR·TI 3종을 병행 계산하고, `demux_g.py` 는
HREMD 궤적을 워커 단위로 되살려 상관시간을 다시 잰다(디렉터리 단위 g 는 7~46배 과소평가).
`build_edges_summary.py` 는 원 JSON 에서 6엣지 표를, `bootstrap_rmse_ci.py` 는 RMSE 95% 구간을
낸다 — 시드 고정 10,000회와 6⁶ 전수 열거가 모두 0.434–1.053 을 준다.

경로는 실행 당시 그대로이며 재현하려면 수정이 필요하다.
