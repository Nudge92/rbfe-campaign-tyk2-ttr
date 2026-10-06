# results — 엣지별·회차별 값과 수렴 진단

`edges_summary.csv` 가 README 의 RMSE·MUE·부호 일치의 근거다(생성 = `../analysis/build_edges_summary.py`).
회차별 원값은 `step4_raw.json`(독립 창) · `hrA/hrB/hr2_ddg.json`(HREMD), 수렴 진단은
`overlap24.txt`(겹침) · `hremd_accept.txt`·`accept24.txt`(교환 수락률) · `geff.json`(상관시간)에 있다.
`ligands_ki.yml` 은 실험 Ki 원본(doi 10.1016/j.ejmech.2013.03.070 Table 4)이다.
