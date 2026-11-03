"""Rebuild the slide entries of presentation/prelim-slides/OUTLINE.md from the deck source.

Script: outline_from_deck.py
Purpose: Keep OUTLINE.md's slide-by-slide section in step with
         ssbmf_talk.qmd (prelim Part III), so the outline describes the deck as built.
Author: Andrew Walther
Created: 2026-10-28 (adapted for Part III 2026-11-03 from SpatialCRT IncidenceDesign)
Dependencies: Python 3 standard library only

Run from the project root:  python3 presentation/prelim-slides/tools/outline_from_deck.py
Then rebuild the visual page:  python3 presentation/prelim-slides/tools/build_outline_html.py

Everything in OUTLINE.md above the line "<!-- SLIDES (generated) -->" is kept as
written (preamble, timing, decisions). Below it, one "### n. title" entry per slide is
written with Figure, Slide text (the slide's visible text, condensed), Notes (the full
script, or the backup note), Time and Carry-over. Carry-over tags come
from the CARRY dict below; a slide missing from it is reported, not silently tagged.
"""
import html
import re
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
QMD = HERE / 'ssbmf_talk.qmd'
OUT = HERE / 'OUTLINE.md'
MARK = '<!-- SLIDES (generated) -->'

# Keys are each main slide's data-menu-title; backup slides are tagged technical.
CARRY = {
    'Title': 'talk-only',
    'Outline': 'talk-only',
    'PDAC': 'prelim',
    'Two-step analysis': 'prelim',
    'Existing methods': 'prelim',
    'Aims': 'prelim',
    'Notation': 'prelim',
    'Joint model': 'prelim',
    'Priors': 'prelim',
    'Fitting': 'prelim',
    'Prediction': 'prelim',
    'Fitted vs finalized': 'prelim',
    'Several data types': 'prelim',
    'Simulation design': 'prelim',
    'Simulation concordance': 'prelim',
    'Starting rank': 'prelim',
    'Planned simulation test': 'prelim',
    'Application pipeline': 'prelim',
    'External concordance': 'prelim',
    'Risk groups': 'prelim',
    'Subtype programs': 'prelim',
    'Tissue programs': 'prelim',
    'Planned results': 'prelim',
    'Summary': 'prelim',
    'Limitations': 'prelim',
    'Remaining work': 'prelim',
}


def visible_text(block):
    """Strip notes, HTML tags, fenced-div markers and markdown emphasis from a slide body."""
    block = re.sub(r'::: \{\.notes\}.*?\n:::', '', block, flags=re.S)
    block = re.sub(r'<svg.*?</svg>', ' ', block, flags=re.S)
    block = re.sub(r'<!--.*?-->', ' ', block, flags=re.S)
    block = re.sub(r'<[^>]+>', ' ', block)
    block = re.sub(r'^:::.*$', ' ', block, flags=re.M)
    block = re.sub(r'^#+\s*', '', block, flags=re.M)
    block = block.replace('**', '').replace('&middot;', '·').replace('&nbsp;', ' ')
    block = html.unescape(re.sub(r'^\s*- ', '', block, flags=re.M))
    return re.sub(r'\s+', ' ', block).strip()


def main():
    """Parse the deck into slides and rewrite the generated part of OUTLINE.md."""
    src = QMD.read_text()
    slides = re.split(r'(?m)^(?=## )', src)[1:]
    entries, section, n_main, n_back, missing = [], None, 0, 0, []
    for sl in slides:
        head = sl.split('\n', 1)[0]
        title = re.sub(r'\s*\{.*\}\s*$', '', head[3:]).replace('Backup &middot; ', '').strip() or 'Title'
        menu = re.search(r'data-menu-title="([^"]+)"', head).group(1)
        part = re.search(r'data-footer-part="([^"]*)"', head).group(1)
        backup = menu.startswith('B:')
        part = part or ('Introduction' if menu in ('Title', 'Outline') else 'Discussion')
        if part != section:
            section = part
            if not backup:
                n_sec = len([e for e in entries if e.startswith('\n## ')]) + 1
            entries.append('\n## ' + ('Backup' if backup else f'{n_sec}. {part}') + '\n')
        figs = re.findall(r'src="figures/([^"]+)"', sl)
        notes = re.search(r'::: \{\.notes\}\n(.*?)\n:::', sl, flags=re.S)
        notes = notes.group(1) if notes else ''
        tm = re.match(r'\[([\d.]+) min', notes.strip())
        body = re.sub(r'^\[[^\]]*\]\s*', '', notes.strip())
        sents = [x[2:] if x.startswith('- ') else x for x in body.split('\n') if x.strip()]
        note = ' '.join(sents)  # full script, so outline.html shows what is said
        text = visible_text(sl.split('\n', 1)[1] if '\n' in sl else '')
        if len(text) > 420:
            text = text[:417].rsplit(' ', 1)[0] + ' …'
        if backup:
            n_back += 1
            label = f'B{n_back}'
        else:
            n_main += 1
            label = str(n_main)
        tag = 'technical' if backup else CARRY.get(menu)
        if tag is None:
            missing.append(menu)
        lines = [f'### {label}. {title}',
                 f'- Figure: {figs[0] if figs else "none"}',
                 f'- Slide text: {text or "(title or image only)"}',
                 f'- Notes: {note}']
        if tm:
            lines.append(f'- Time: {tm.group(1)}')
        lines.append(f'- Carry-over: {tag or "UNTAGGED"}')
        entries.append('\n'.join(lines) + '\n')
    pre = OUT.read_text().split(MARK)[0].rstrip() + '\n\n' + MARK + '\n'
    OUT.write_text(pre + '\n'.join(entries))
    print(f'Wrote {OUT}: {n_main} main, {n_back} backup slides')
    for m in missing:
        print('WARNING: no carry-over tag for', m)


if __name__ == '__main__':
    main()
