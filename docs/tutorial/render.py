"""Render both tutorial PDFs from this folder (docs/tutorial/) into the repository root.
Usage: python render.py [--anonymous]    (needs the Python package playwright with its Chromium: pip install playwright;
playwright install chromium). --anonymous writes "anonymous" in place of the GitHub owner of the repository (read from
repository-code in CITATION.cff) in the PDFs only, for double-blind review; the HTML files are not changed.
check.py <file.html> verifies that every piece of text lies inside its box (SVG rectangles and HTML cards, boxes and cells)."""
import asyncio, os, re, sys
from playwright.async_api import async_playwright
HERE = os.path.dirname(os.path.abspath(__file__))
JOBS = [("tutorial.html", "../../TUTORIAL_MEGAN_SC_pt.pdf", "MEGAN_SC · tutorial ilustrado"),
        ("tutorial_en.html", "../../TUTORIAL_MEGAN_SC_en.pdf", "MEGAN_SC · illustrated tutorial")]
FOOT = ('<div style="width:100%;font-family:Inter,sans-serif;font-size:8px;color:#8a94a6;padding:0 17mm;display:flex;'
        'justify-content:space-between"><span>{}</span><span><span class="pageNumber"></span> / <span class="totalPages"></span></span></div>')
ANON = '--anonymous' in sys.argv
OWNER = re.search(r'github\.com/([^/"\s]+)/', open(os.path.join(HERE, '..', '..', 'CITATION.cff'), encoding='utf-8').read()).group(1)


def source(html):
    """The page to render: the HTML file itself, or a temporary copy without the owner when --anonymous."""
    if not ANON:
        return os.path.join(HERE, html)
    tmp = os.path.join(HERE, '_anon_' + html)
    open(tmp, 'w', encoding='utf-8').write(open(os.path.join(HERE, html), encoding='utf-8').read().replace(OWNER, 'anonymous'))
    return tmp


async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch()
        for html, pdf, label in JOBS:
            if not os.path.exists(os.path.join(HERE, html)):
                print("skipped (not found):", html); continue
            pg = await b.new_page()
            src = source(html)
            await pg.goto("file://" + src, wait_until="networkidle")
            await pg.pdf(path=os.path.join(HERE, pdf), format="A4", print_background=True, display_header_footer=True,
                         header_template="<div></div>", footer_template=FOOT.format(label),
                         margin={"top": "14mm", "bottom": "16mm", "left": "0mm", "right": "0mm"})
            await pg.close()
            if ANON: os.remove(src)
            print("written", os.path.normpath(os.path.join(HERE, pdf)))
        await b.close()
asyncio.run(main())
