#!/usr/bin/env python3
"""Convert percent-format sources (# %% / # %% [markdown]) to .ipynb.  Placeholders:
   __FIXES_PY__  -> contents of src/megan_sc/fixes.py ;  __FILE:<relpath>__ -> contents of that file (path relative to the
   repository root).
usage: python build_ipynb.py <src.py> [...]   (writes ../<name>.ipynb)"""
import sys, os, json, re
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
def expand(text):
    def rep(m):
        p = os.path.join(ROOT, m.group(1)); s = open(p, encoding='utf-8').read()
        assert "'''" not in s, f"{p} contains triple single quotes"
        return s
    text = text.replace('__FIXES_PY__', '__FILE:src/megan_sc/fixes.py__')
    return re.sub(r'__FILE:([^_][^ ]*?\.(?:py|sh|csv))__', rep, text)
def convert(src):
    lines = open(src, encoding='utf-8').read().split('\n')
    cells, cur, kind = [], [], None
    def flush():
        if kind is None: return
        body = '\n'.join(cur).strip('\n')
        if kind == 'markdown':
            body = '\n'.join(l[2:] if l.startswith('# ') else ('' if l == '#' else l) for l in body.split('\n'))
            cells.append({'cell_type': 'markdown', 'metadata': {}, 'source': body.splitlines(True)})
        else:
            body = expand(body)
            cells.append({'cell_type': 'code', 'metadata': {}, 'execution_count': None, 'outputs': [], 'source': body.splitlines(True)})
    for l in lines:
        if l.startswith('# %%'):
            flush(); cur = []; kind = 'markdown' if '[markdown]' in l else 'code'; continue
        cur.append(l)
    flush()
    nb = {'cells': cells, 'metadata': {'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
                                       'language_info': {'name': 'python'}, 'colab': {'provenance': []}},
          'nbformat': 4, 'nbformat_minor': 5}
    out = os.path.join(os.path.dirname(HERE), os.path.splitext(os.path.basename(src))[0] + '.ipynb')
    json.dump(nb, open(out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    # syntax check of all code cells (magics stripped)
    import ast
    for i, c in enumerate(cells):
        if c['cell_type'] != 'code': continue
        code = ''.join(c['source'])
        code = '\n'.join('pass' if l.lstrip().startswith(('%', '!')) else l for l in code.split('\n'))
        try: ast.parse(code)
        except SyntaxError as e: raise SystemExit(f'{out}: cell {i}: {e}')
    print('wrote', out, len(cells), 'cells')
for s in sys.argv[1:]: convert(s)
