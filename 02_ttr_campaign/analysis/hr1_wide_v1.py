#!/usr/bin/env python3
"""hr1_wide — 24 레그 와이드 표의 파생 계산. 즉석 실행분(2026-10-04 09:44:45 KST · 세션 전사 줄 7117)을 그대로 옮겨 저장했다.
입력: qsar/reports/hr1_table_prod_v14_N17.json (hr1_table_v1.py 의 출력). ★교환 로그·xvg 를 다시 읽지 않는다 — 그 파일의 값에서 파생량만 낸다.
파생량: ratio·sq·cor·res·drift·c5·rtmin·rtmed·rtmax·rtzero·rtt·c7·gd_rtt·blind. 나머지 필드는 입력 값을 이름만 바꿔 옮긴다.
★ratio·sq·cor·res·gd_rtt 는 근거가 철회된 양이다(HR_RETRACTION_20261004.md). 동결 파일을 그대로 다시 만들기 위해서만 남긴다.
★필드 이름 c7 은 사전등재 기준이 아니다 — "사후 진단: 워커당 왕복 ≥1" 으로 읽는다(같은 문서 §5).
사용: hr1_wide_v1.py <출력 디렉터리> [--table <대체 hr1_table json>]  > <출력 디렉터리>/hr1_wide.txt
★<출력 디렉터리>/hr1_wide_prod_v14_N17.json 이 이미 있으면 죽는다(덮어쓰지 않는다).
"""
import sys
from pathlib import Path
_T = Path("/home/nudge/Project/CADD/ttr")
OUT = Path(sys.argv[1]).resolve()
assert OUT.is_dir(), f"출력 디렉터리가 없다: {OUT}"
assert OUT != _T and _T not in OUT.parents, "출력은 저장소 밖 임시 디렉터리여야 한다 — 동결 산출물 경로에는 쓰지 않는다"
TABLE = sys.argv[sys.argv.index("--table")+1] if "--table" in sys.argv else str(_T/"qsar/reports/hr1_table_prod_v14_N17.json")
# ──── 아래는 즉석 코드 그대로 (바뀐 줄 3: 입력 경로 · 출력 경로[쓰기 모드 w→x] · 마지막 저장 print) ────
import json, math, numpy as np
R=json.load(open(TABLE))
T_PS=8000.0
OILED={('r1','cplxL1'),('r1','cplxL2'),('r1','cplxL3'),('r1','solvL1')}   # v27 §③ 오염분
rows=[]
for r in R:
    rt=np.array(r['rt_per']) if r['rt_per'] else None
    sq=math.sqrt(r['g_demux']/r['g_naive'])
    cor=r['mbar_sig']*sq
    res=r['blk_se']/cor
    b=r['blocks']; drift=b[3]-b[1]
    rtt=T_PS/np.mean(rt) if rt is not None and np.mean(rt)>0 else float('inf')
    rows.append(dict(rep=r['rep'], leg=r['leg'], K=r['nwin'],
        dlmin=r['dl_min'], dlmax=r['dl_max'], dlmean=r['dl_mean'],
        sig=r['mbar_sig'], gn=r['g_naive'], gd=r['g_demux'], ratio=r['g_demux']/r['g_naive'],
        sq=sq, cor=cor, bse=r['blk_se'], res=res,
        b=b, drift=drift,
        amean=r['acc_mean'], amin=r['acc_min'], apair=r['acc_min_pair'],
        c5='통과' if r['acc_mean']>=0.20 else '★미달',
        rtmin=int(rt.min()), rtmed=float(np.median(rt)), rtmax=int(rt.max()),
        rtzero=int((rt==0).sum()), rttot=r['rt_total'], rtt=rtt,
        c7='통과' if rt.min()>=1 else '★미달',
        gd_rtt=r['g_demux']/rtt if np.isfinite(rtt) else float('nan'),
        blind='오염' if (r['rep'],r['leg']) in OILED else '맹검'))

print("═"*204)
print("HR-1 · 24레그 와이드 표   (prod_v14 · N17 · 8 ns 전량 · 40,001 프레임/창)")
print("═"*204)
h=(f"{'회차':<4}{'다리':<8}{'창':>3}{'Δλmin':>8}{'Δλmax':>8}{'Δλ평균':>8}│"
   f"{'MBARσ':>8}{'g_naive':>9}{'g_demux':>9}{'비':>8}{'√비':>7}{'보정σ':>8}{'블록SE':>8}{'★잔차비':>9}│"
   f"{'드리프트':>9}│{'수락평균':>9}{'★수락최소':>10}{'쌍':>4}{'C-5':>6}│"
   f"{'왕복계':>7}{'최소':>5}{'중앙':>6}{'최대':>5}{'★0회':>6}{'C-7':>6}│{'RTT(ps)':>9}{'gd/RTT':>8}│{'구분':>5}")
print(h); print("─"*204)
for x in rows:
    print(f"{x['rep']:<4}{x['leg']:<8}{x['K']:>3}{x['dlmin']:>8.4f}{x['dlmax']:>8.4f}{x['dlmean']:>8.4f}│"
          f"{x['sig']:>8.4f}{x['gn']:>9.3f}{x['gd']:>9.2f}{x['ratio']:>8.2f}{x['sq']:>7.2f}{x['cor']:>8.4f}{x['bse']:>8.4f}{x['res']:>9.2f}│"
          f"{x['drift']:>+9.3f}│{x['amean']:>9.4f}{x['amin']:>10.4f}{x['apair']:>4}{x['c5']:>6}│"
          f"{x['rttot']:>7}{x['rtmin']:>5}{x['rtmed']:>6.1f}{x['rtmax']:>5}{x['rtzero']:>6}{x['c7']:>6}│"
          f"{x['rtt']:>9.1f}{x['gd_rtt']:>8.3f}│{x['blind']:>5}")
print("═"*204)
json.dump(rows, open(OUT/'hr1_wide_prod_v14_N17.json','x'), ensure_ascii=False, indent=1, default=float)
print(f"저장 {OUT/'hr1_wide_prod_v14_N17.json'}")
