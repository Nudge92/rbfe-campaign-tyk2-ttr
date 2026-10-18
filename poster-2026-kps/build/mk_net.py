# -*- coding: utf-8 -*-
"""FEP 섭동 그물 — 노드 10 · 엣지 17"""
import matplotlib; matplotlib.use('Agg')
from matplotlib import font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import math, csv
fm.fontManager.addfont('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
KR='Noto Sans CJK KR'; plt.rcParams['font.family']=KR
INK,INK2,MUTE,LINE='#15202B','#3E4C57','#6E7A83','#D5DCD6'
TEAL,TEAL_S='#1A6460','#E4EFEE'; AMB,AMB_S='#B4652A','#FBF0E3'; SUNK='#F3F5F2'

import os
F=os.path.join(os.path.dirname(__file__),'..','data','slide_edges.csv')
rows=list(csv.DictReader(open(F,encoding='utf-8-sig')))
EDGES=[(r['A_실행방향'],r['B_실행방향'],r['실험ΔΔG'] not in ('미보유','')) for r in rows]
ORDER=['CHEMBL241454','CHEMBL240808','CHEMBL438498','cand0023','cand0030',
       'cand0028','cand0118','cand0132','cand0042','cand0271']
SHORT={n:(n.replace('CHEMBL','') if n.startswith('CHEMBL') else n.replace('cand','')) for n in ORDER}

W,H=2400,1120
fig=plt.figure(figsize=(W/100,H/100),dpi=100)
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,W); ax.set_ylim(H,0); ax.axis('off')
fig.patch.set_facecolor('white')

CX,CY,RX,RY=1200,580,880,300
POS={}
for i,n in enumerate(ORDER):
    th=-math.pi/2+2*math.pi*i/len(ORDER)
    POS[n]=(CX+RX*math.cos(th), CY+RY*math.sin(th))

ax.text(1200,62,'선 하나 = 두 분자의 결합력 차이를 계산한 1건',
        fontsize=40,fontweight='bold',ha='center',va='center',color=INK,zorder=6,fontfamily=KR)
ax.text(1200,128,'결합력은 직접 못 구한다. 닮은 쌍끼리만 이을 수 있어 사슬처럼 엮어 전체를 비교한다',
        fontsize=28,ha='center',va='center',color=MUTE,zorder=6,fontfamily=KR)

for a,b,hot in EDGES:
    (x1,y1),(x2,y2)=POS[a],POS[b]
    ax.plot([x1,x2],[y1,y2], color=(AMB if hot else LINE),
            lw=(6.0 if hot else 2.2), zorder=(3 if hot else 1), solid_capstyle='round')

NW,NH=258,122
for n in ORDER:
    x,y=POS[n]; anchor=n.startswith('CHEMBL')
    ax.add_patch(FancyBboxPatch((x-NW/2,y+NH/2),NW,-NH,
        boxstyle="round,pad=0,rounding_size=9",
        facecolor=(SUNK if anchor else TEAL), edgecolor=(MUTE if anchor else 'none'),
        linewidth=(2.4 if anchor else 0), zorder=5))
    ax.text(x,y-26,'앵커' if anchor else '후보',fontsize=23,fontweight='bold',
            ha='center',va='center',color=(MUTE if anchor else '#9DC6C2'),
            zorder=6,fontfamily=KR)
    ax.text(x,y+22,SHORT[n],fontsize=33,fontweight='bold',ha='center',va='center',
            color=(INK2 if anchor else 'white'),zorder=6,fontfamily=KR)

# 범례
ax.plot([690,790],[1022,1022],color=AMB,lw=7,solid_capstyle='round',zorder=5)
ax.text(815,1022,'실험값으로 대조할 수 있었던 엣지 3개 — 모두 앵커끼리다',fontsize=27,
        fontweight='bold',ha='left',va='center',color=AMB,zorder=6,fontfamily=KR)
ax.text(1200,1076,'타파미디스는 이 그물에 들어 있지 않다',fontsize=26,
        ha='center',va='center',color=MUTE,zorder=6,fontfamily=KR)

fig.savefig(os.path.join(os.path.dirname(__file__),'..','fig','network.png'),dpi=100,facecolor='white'); plt.close(fig)
print('network.png  엣지', len(EDGES), '· 대조 엣지', sum(h for _,_,h in EDGES))
