#!/usr/bin/env python3
"""out/TTR_포스터_95x150.pdf → readme/*.jpg

README 에 스크롤로 읽히는 블록별 확대 이미지를 만든다.
좌표는 Playwright 로 poster.html 의 각 블록 경계를 실측한 값(mm)이다.
포스터 레이아웃을 바꾸면 좌표도 다시 재야 한다.
"""
import os, subprocess
from PIL import Image, ImageChops

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDF  = os.path.join(ROOT, 'out', 'TTR_포스터_95x150.pdf')
OUT  = os.path.join(ROOT, 'readme')
DPI  = 110
K    = DPI / 25.4                      # mm → px

# (파일명, x, y, w, h) 단위 mm — 페이지 950 × 1500
CROPS = [
    ('02_표적',      20,  163, 342, 303),
    ('03_착안',     356,  163, 574, 303),
    ('04_깔때기',    20,  460, 628, 304),
    ('05_생성모델',  650,  552, 276, 142),
    # 후보 7종 줄은 가로 898mm 라 통짜로는 글씨가 안 읽힌다. 분자 4 | 3 으로 쪼갠다.
    ('06a_후보7종_좌',150,  784, 378,  89),
    ('06b_후보7종_우',513,  784, 290,  89),
    ('07_구조확인',  20,  927, 405, 278),
    ('08_결합력계산',418,  927, 512, 278),
    ('09_검증과한계',20, 1199, 503, 233),
    ('10_추후과제', 516, 1199, 414, 233),
]

def trim(path, margin=14):
    """가장자리의 순백 여백만 깎는다. 카드 테두리·틴트 배경은 내용이라 남는다."""
    im = Image.open(path).convert('RGB')
    bg = Image.new('RGB', im.size, (255, 255, 255))
    mask = ImageChops.difference(im, bg).convert('L').point(lambda v: 255 if v > 6 else 0)
    bb = mask.getbbox()
    if not bb:
        return
    box = (max(0, bb[0]-margin), max(0, bb[1]-margin),
           min(im.width, bb[2]+margin), min(im.height, bb[3]+margin))
    if box != (0, 0, im.width, im.height):
        im.crop(box).save(path, quality=90)

def main():
    os.makedirs(OUT, exist_ok=True)
    for name, x, y, w, h in CROPS:
        px = [str(round(v * K)) for v in (x, y, w, h)]
        subprocess.run(['pdftoppm', '-r', str(DPI), '-jpeg', '-jpegopt', 'quality=90',
                        '-x', px[0], '-y', px[1], '-W', px[2], '-H', px[3],
                        '-singlefile', PDF, os.path.join(OUT, name)], check=True)
        f = os.path.join(OUT, name + '.jpg')
        trim(f)
        print(f'{name:14} {Image.open(f).size}  {os.path.getsize(f)//1024} KB')

if __name__ == '__main__':
    main()
