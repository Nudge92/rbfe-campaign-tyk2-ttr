const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch({ executablePath: '/opt/pw-browsers/chromium' });
  const p = await b.newPage();
  const errs = [];
  p.on('pageerror', e => errs.push(e.message));
  await p.goto('file://' + process.cwd() + '/poster.html');
  await p.waitForTimeout(1500);
  const h = await p.evaluate(() => {
    const pg = document.querySelector('.page');
    return { scroll: document.body.scrollHeight, page: pg.scrollHeight,
             over: pg.scrollHeight - pg.clientHeight };
  });
  console.log('본문 높이', h.page, 'px · 넘침', h.over, 'px  (1500mm ≈', Math.round(1500/25.4*96), 'px)');
  const imgs = await p.$$eval('img', a => a.map(i => i.naturalWidth > 0));
  console.log('이미지', imgs.filter(Boolean).length + '/' + imgs.length, errs.length ? 'JS오류:'+errs : '· JS오류 없음');
  await p.pdf({ path: 'out/TTR_포스터_95x150.pdf', width: '950mm', height: '1500mm',
                printBackground: true, pageRanges: '1' });
  await b.close();
})();
