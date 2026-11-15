"""check.py <html>: every piece of text must lie inside its box (SVG rectangles; HTML cards, boxes, cells and code blocks).
Example: python check.py "$PWD/tutorial_en.html". Prints one line per problem and the count; "tight" SVG texts (margin under 3 px) are listed apart."""
import asyncio, sys, json
from playwright.async_api import async_playwright
JS = r"""
() => {
  const out = [];
  const inside = (a, b, tol) => a.left >= b.left - tol && a.right <= b.right + tol && a.top >= b.top - tol && a.bottom <= b.bottom + tol;
  // SVG texts vs the smallest rect that contains their centre
  document.querySelectorAll('svg').forEach((svg, si) => {
    const rects = [...svg.querySelectorAll('rect')].map(r => r.getBoundingClientRect());
    const sb = svg.getBoundingClientRect();
    svg.querySelectorAll('text').forEach(t => {
      const tb = t.getBoundingClientRect(); if (!tb.width) return;
      const cx = (tb.left + tb.right) / 2, cy = (tb.top + tb.bottom) / 2;
      const cont = rects.filter(r => cx >= r.left && cx <= r.right && cy >= r.top && cy <= r.bottom).sort((a, b) => a.width * a.height - b.width * b.height);
      const box = cont.length ? cont[0] : sb;
      if (!inside(tb, box, 0.5)) out.push({kind: 'svg' + si, text: t.textContent.slice(0, 60), over: [Math.round(box.left - tb.left), Math.round(tb.right - box.right), Math.round(box.top - tb.top), Math.round(tb.bottom - box.bottom)]});
      else { const m = Math.min(tb.left - box.left, box.right - tb.right); if (cont.length && m < 3) out.push({kind: 'svg-tight' + si, text: t.textContent.slice(0, 60), margin: +m.toFixed(1)}); }
    });
  });
  // HTML text nodes vs their nearest box-like ancestor
  const SEL = '.card,.nbcard,.box,.trail,.kpi>div,.legend>div,pre,td,th,figcaption,.cover,.toc a,.step>div,dd,dt,.pill,.chapnum';
  const w = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
  let n;
  while ((n = w.nextNode())) {
    if (!n.textContent.trim() || n.parentElement.closest('svg')) continue;
    const c = n.parentElement.closest(SEL); if (!c) continue;
    const cb = c.getBoundingClientRect();
    const rg = document.createRange(); rg.selectNodeContents(n);
    for (const r of rg.getClientRects()) {
      if (!r.width) continue;
      if (!inside(r, cb, 0.6)) { out.push({kind: 'html:' + c.className + c.tagName, text: n.textContent.trim().slice(0, 60), over: [Math.round(cb.left - r.left), Math.round(r.right - cb.right), Math.round(cb.top - r.top), Math.round(r.bottom - cb.bottom)]}); break; }
    }
  }
  document.querySelectorAll('pre,td,.card,.nbcard,.box').forEach(e => { if (e.scrollWidth > e.clientWidth + 1) out.push({kind: 'scroll:' + e.tagName + '.' + e.className, text: e.textContent.trim().slice(0, 60), over: [e.scrollWidth - e.clientWidth]}); });
  return out;
}
"""
async def main(path):
    async with async_playwright() as p:
        b = await p.chromium.launch(); pg = await b.new_page(viewport={'width': 665, 'height': 1000})
        await pg.emulate_media(media='print')
        await pg.goto('file://' + path, wait_until='networkidle')
        res = await pg.evaluate(JS)
        await b.close()
    bad = [r for r in res if 'tight' not in r['kind']]
    for r in res: print(json.dumps(r, ensure_ascii=False))
    print('problems:', len(bad), '| tight svg texts:', len(res) - len(bad))
asyncio.run(main(sys.argv[1]))
