"""Build presentation/prelim-slides/outline.html from presentation/prelim-slides/OUTLINE.md.

Script: build_outline_html.py
Purpose: One visual page listing every planned slide, section by section, with its
         figure inline, so the outline can be reviewed before any slide is built.
Author: Andrew Walther
Created: 2026-10-25 (adapted for Part III 2026-11-03 from SpatialCRT IncidenceDesign)
Dependencies: Python 3 standard library only

Run from the project root:  python3 presentation/prelim-slides/tools/build_outline_html.py

OUTLINE.md is the single source. Its structure:
  ## <section>                      a talk section (or "Backup")
  ### <number>. <title>             one slide
  - <Field>: <value>                Figure, Slide text, Notes, Time, Carry-over, Prelim note
Anything before the first "## 1." heading (open questions, timing table) is shown
as a preamble. Figure values that name a file in presentation/prelim-slides/figures/ are shown
inline; any other value is printed as text. A slide without a figure, or whose
figure file is missing, is reported on stdout so nothing is dropped silently.
Figures are embedded as base64 data URIs, so outline.html is a single portable file.
"""
import base64
import html
import re
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
SRC = HERE / 'OUTLINE.md'
OUT = HERE / 'outline.html'
FIG_DIR = HERE / 'figures'

TAG_CLASS = {'prelim': 'tag-prelim', 'technical': 'tag-tech', 'talk-only': 'tag-talk'}


def inline(text):
    """Escape text and render **bold** and `code` spans."""
    t = html.escape(text)
    t = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', t)
    return re.sub(r'`([^`]+)`', r'<code>\1</code>', t)


def parse(lines):
    """Split OUTLINE.md into a preamble and a list of (section, [slides])."""
    preamble, sections, slide = [], [], None
    for line in lines:
        if line.startswith('## ') and (sections or re.match(r'## (\d|Backup)', line)):
            sections.append((line[3:].strip(), []))
            slide = None
        elif line.startswith('### ') and sections:
            slide = {'title': line[4:].strip(), 'fields': {}}
            sections[-1][1].append(slide)
        elif slide is not None and line.startswith('- ') and ':' in line:
            key, value = line[2:].split(':', 1)
            slide['fields'][key.strip()] = value.strip()
        elif not sections:
            preamble.append(line)
    return preamble, sections


def render_preamble(lines):
    """Render the preamble's headings, bullets, tables and paragraphs."""
    out, table, para = [], [], []

    def flush():
        if para:
            out.append('<p>' + inline(' '.join(para)) + '</p>')
            para.clear()
        if table:
            rows = [r for r in table if not re.match(r'^\|[-| ]+\|$', r)]
            cells = [[c.strip() for c in r.strip('|').split('|')] for r in rows]
            out.append('<table><thead><tr>' + ''.join(f'<th>{inline(c)}</th>' for c in cells[0]) +
                       '</tr></thead><tbody>' + ''.join(
                           '<tr>' + ''.join(f'<td>{inline(c)}</td>' for c in r) + '</tr>' for r in cells[1:]) +
                       '</tbody></table>')
            table.clear()

    bullets = []
    for line in lines:
        if line.startswith('|'):
            table.append(line)
            continue
        if line.startswith('- '):
            if para:
                flush()
            bullets.append(line[2:])
            continue
        if line.startswith('  ') and bullets:
            bullets[-1] += ' ' + line.strip()
            continue
        if bullets:
            out.append('<ul>' + ''.join(f'<li>{inline(b)}</li>' for b in bullets) + '</ul>')
            bullets = []
        if line.startswith('# '):
            flush(); out.append(f'<h1>{inline(line[2:])}</h1>')
        elif line.startswith('## '):
            flush(); out.append(f'<h2>{inline(line[3:])}</h2>')
        elif not line.strip():
            flush()
        else:
            para.append(line.strip())
    if bullets:
        out.append('<ul>' + ''.join(f'<li>{inline(b)}</li>' for b in bullets) + '</ul>')
    flush()
    return '\n'.join(out)


def render_slide(s, problems):
    """Render one slide card; record missing figures in problems."""
    f = s['fields']
    fig = f.get('Figure', '')
    fig_file = fig.split()[0] if fig else ''
    if fig_file and (FIG_DIR / fig_file).exists():
        # Embedded as a data URI so the page is one self-contained file
        data = base64.b64encode((FIG_DIR / fig_file).read_bytes()).decode()
        figure = f'<img src="data:image/png;base64,{data}" alt="{html.escape(fig_file)}" loading="lazy">'
        extra = fig[len(fig_file):].strip()
        if extra:
            figure += f'<p class="figalt">{inline(extra)}</p>'
    else:
        if fig_file.endswith('.png'):
            problems.append(f'missing figure file: {fig_file} ({s["title"]})')
        figure = f'<div class="nofig">{inline(fig or "no figure listed")}</div>'
    tag = f.get('Carry-over', '').strip()
    tag_html = f'<span class="tag {TAG_CLASS.get(tag, "")}">{inline(tag)}</span>' if tag else ''
    time = f'<span class="time">{inline(f["Time"])} min</span>' if 'Time' in f else ''
    rows = ''.join(f'<dt>{inline(k)}</dt><dd>{inline(f[k])}</dd>'
                   for k in ('Slide text', 'Notes', 'Prelim note') if k in f)
    return (f'<article class="slide"><header><h3>{inline(s["title"])}</h3>'
            f'<div class="meta">{time}{tag_html}</div></header>'
            f'<div class="body"><div class="fig">{figure}</div><dl>{rows}</dl></div></article>')


CSS = """
:root { --bg:#fbfbfa; --card:#ffffff; --ink:#1d1f21; --muted:#5b6168; --line:#dcdfe3;
  --accent:#13294b; --prelim:#1f6f4a; --tech:#8a4b00; --talk:#5b4a8a; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg:#16181b; --card:#1f2226; --ink:#e7e9ec; --muted:#a3a9b0; --line:#33383e;
  --accent:#9cc3ff; --prelim:#6fd3a0; --tech:#f0b469; --talk:#c2b2f2; } }
:root[data-theme="dark"] { --bg:#16181b; --card:#1f2226; --ink:#e7e9ec; --muted:#a3a9b0;
  --line:#33383e; --accent:#9cc3ff; --prelim:#6fd3a0; --tech:#f0b469; --talk:#c2b2f2; }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
  font:16px/1.5 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; }
main { max-width:1180px; margin:0 auto; padding:24px 16px 64px; }
h1 { font-size:1.6rem; margin:0 0 8px; color:var(--accent); }
h2 { font-size:1.25rem; margin:36px 0 12px; padding-bottom:6px; border-bottom:2px solid var(--line); }
table { border-collapse:collapse; margin:8px 0 16px; }
th, td { border:1px solid var(--line); padding:4px 10px; text-align:left; }
nav { position:sticky; top:0; background:var(--bg); padding:8px 0; border-bottom:1px solid var(--line);
  display:flex; flex-wrap:wrap; gap:8px 16px; font-size:.9rem; z-index:1; }
nav a { color:var(--accent); text-decoration:none; }
.slide { background:var(--card); border:1px solid var(--line); border-radius:8px; margin:14px 0; padding:14px 16px; }
.slide header { display:flex; justify-content:space-between; gap:12px; align-items:baseline; flex-wrap:wrap; }
.slide h3 { margin:0; font-size:1.05rem; }
.meta { display:flex; gap:8px; align-items:center; font-size:.8rem; color:var(--muted); }
.tag { border:1px solid currentColor; border-radius:10px; padding:0 8px; }
.tag-prelim { color:var(--prelim); } .tag-tech { color:var(--tech); } .tag-talk { color:var(--talk); }
.body { display:grid; grid-template-columns:minmax(0,1.1fr) minmax(0,1fr); gap:16px; margin-top:10px; }
@media (max-width:760px) { .body { grid-template-columns:1fr; } }
.fig img { width:100%; height:auto; border:1px solid var(--line); background:#fff; }
.nofig { border:1px dashed var(--line); color:var(--muted); padding:18px; font-size:.9rem; }
.figalt { color:var(--muted); font-size:.85rem; margin:4px 0 0; }
dl { margin:0; } dt { font-weight:600; font-size:.85rem; color:var(--muted); margin-top:6px; }
dd { margin:0 0 4px; }
code { font-size:.9em; }
"""


def main():
    """Read OUTLINE.md, write outline.html, and report anything not rendered."""
    lines = SRC.read_text().splitlines()
    preamble, sections = parse(lines)
    problems = []
    nav = ''.join(f'<a href="#s{i}">{inline(name)}</a>' for i, (name, _) in enumerate(sections))
    body = []
    for i, (name, slides) in enumerate(sections):
        body.append(f'<h2 id="s{i}">{inline(name)}</h2>')
        body.extend(render_slide(s, problems) for s in slides)
    n_slides = sum(len(s) for _, s in sections)
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>Part III Slide Outline</title><style>{CSS}</style></head><body><main>'
            f'{render_preamble(preamble)}<nav>{nav}</nav>{"".join(body)}'
            f'<p class="figalt">Built from OUTLINE.md: {n_slides} slides in {len(sections)} sections.</p>'
            f'</main></body></html>')
    OUT.write_text(page)
    print(f'Wrote {OUT} ({n_slides} slides, {len(sections)} sections)')
    for p in problems:
        print('WARNING:', p)


if __name__ == '__main__':
    main()
