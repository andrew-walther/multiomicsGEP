#!/usr/bin/env python3
"""
prelim_transform.py -- generate the bios-dissertation Chapter 4 draft from the
standalone multiomicsGEP chapter.

Called by sync_to_prelim.sh; not normally run by hand. Standard library only,
so it runs under any python3 (including inside the git post-commit hook).

Steps, in order (nothing is written unless step 1 passes):

  1. Bib check (fail loud). Every @key cited in the chapter must exist in the
     prelim master bib, and for each cited key the core fields (author, title,
     year, journal, volume, number, pages, doi) must agree between the
     chapter bib and the master bib. Printing-only fields and `annote` are
     ignored on purpose (the multiomicsGEP bib drops some and adds `annote`).
     Any problem -> message on stderr, exit 2, no files touched.
  2. Transform the chapter text (see transform() for the numbered T-steps).
  3. Write <out-dir>/project3-ssbmf-draft.qmd if its content changed.
  4. Copy every figure the draft \\includegraphics-references into
     <out-dir>/figures/ (same subpath), and delete files there that are no
     longer referenced.

Prints one machine-readable line on success:
  RESULT qmd=<changed|unchanged> figs_copied=N figs_removed=M keys=K
"""

import argparse
import filecmp
import os
import re
import shutil
import sys

CORE_FIELDS = ("author", "title", "year", "journal", "volume", "number", "pages", "doi")


def die(msg, code=2):
    sys.stderr.write("prelim_transform: ERROR: " + msg + "\n")
    sys.exit(code)


# ---------------------------------------------------------------------------
# Step 1: BibTeX parsing and the citation-key check
# ---------------------------------------------------------------------------

def parse_bib(path):
    """Minimal brace-aware BibTeX parser.

    Returns {key: {field_lower: raw_value}}. Values keep their inner text
    with the outermost {...} or "..." delimiters removed. @comment/@string/
    @preamble blocks are skipped. Duplicate keys are an error (pandoc would
    silently keep one of them).
    """
    text = open(path, encoding="utf-8").read()
    entries = {}
    i, n = 0, len(text)
    while True:
        at = text.find("@", i)
        if at < 0:
            break
        m = re.match(r"@(\w+)\s*([{(])", text[at:])
        if not m:
            i = at + 1
            continue
        etype = m.group(1).lower()
        open_ch = m.group(2)
        close_ch = "}" if open_ch == "{" else ")"
        j = at + m.end()
        # find the matching close of the whole entry
        depth, k = 1, j
        while k < n and depth:
            if text[k] == open_ch:
                depth += 1
            elif text[k] == close_ch:
                depth -= 1
            k += 1
        body = text[j:k - 1]
        i = k
        if etype in ("comment", "string", "preamble"):
            continue
        key, _, rest = body.partition(",")
        key = key.strip()
        if key in entries:
            die("duplicate citekey '%s' in %s" % (key, path))
        entries[key] = parse_fields(rest)
    return entries


def parse_fields(s):
    fields = {}
    i, n = 0, len(s)
    while i < n:
        m = re.compile(r"\s*,?\s*([A-Za-z][\w\-]*)\s*=\s*").match(s, i)
        if not m:
            break
        name = m.group(1).lower()
        i = m.end()
        if i < n and s[i] == "{":
            depth, k = 1, i + 1
            while k < n and depth:
                if s[k] == "{":
                    depth += 1
                elif s[k] == "}":
                    depth -= 1
                k += 1
            val = s[i + 1:k - 1]
            i = k
        elif i < n and s[i] == '"':
            k = s.index('"', i + 1)
            val = s[i + 1:k]
            i = k + 1
        else:
            m2 = re.compile(r"[^,\s]+").match(s, i)
            val = m2.group(0) if m2 else ""
            i = m2.end() if m2 else n
        fields[name] = val
    return fields


def norm(v):
    """Whitespace-insensitive comparison only; no other normalization."""
    return " ".join(v.split())


def cited_keys(qmd_text):
    """@keys cited in the chapter body (YAML, HTML comments, code excluded)."""
    t = strip_yaml(qmd_text)[1]
    t = re.sub(r"<!--.*?-->", "", t, flags=re.S)
    t = re.sub(r"^```.*?^```", "", t, flags=re.S | re.M)
    keys = set()
    for m in re.finditer(r"(?<![\w@])@([A-Za-z0-9_][A-Za-z0-9_:.#$%&\-+?<>~/]*)", t):
        k = m.group(1).rstrip(".:;,?")
        # Quarto cross-references (@fig-x, @tbl-x, @sec-x, @eq-x) are not citations
        if re.match(r"(fig|tbl|sec|eq)-", k):
            continue
        keys.add(k)
    return keys


def check_bib(qmd_text, chapter_bib, master_bib):
    ch = parse_bib(chapter_bib)
    ms = parse_bib(master_bib)
    keys = sorted(cited_keys(qmd_text))
    problems = []
    for k in keys:
        if k not in ms:
            problems.append("cited key '%s' is missing from the master bib %s" % (k, master_bib))
            continue
        if k not in ch:
            problems.append("cited key '%s' is missing from the chapter bib %s" % (k, chapter_bib))
            continue
        for f in CORE_FIELDS:
            a, b = ch[k].get(f), ms[k].get(f)
            if (a is None) != (b is None) or (a is not None and norm(a) != norm(b)):
                problems.append("key '%s' field '%s' differs:\n      chapter bib: %r\n      master bib:  %r"
                                % (k, f, a, b))
    if problems:
        die("bib check failed (%d problem(s)); nothing was written:\n  - %s"
            % (len(problems), "\n  - ".join(problems)))
    return keys


# ---------------------------------------------------------------------------
# Step 2: the text transforms
# ---------------------------------------------------------------------------

def strip_yaml(text):
    """Split off the leading YAML block. Returns (yaml_text, body)."""
    if not text.startswith("---\n"):
        die("chapter does not start with a YAML block")
    end = text.index("\n---\n", 4)
    return text[4:end], text[end + 5:]


def yaml_title(yaml_text):
    m = re.search(r'^title:\s*"(.*)"\s*$', yaml_text, flags=re.M)
    if not m:
        die("could not find a double-quoted title: in the chapter YAML")
    return m.group(1)


def demote_headings(text):
    """T4: add one '#' to every ATX heading outside code fences and comments."""
    out, in_fence, in_comment = [], False, False
    for line in text.split("\n"):
        s = line.lstrip()
        if not in_comment and s.startswith("```"):
            in_fence = not in_fence
        elif not in_fence:
            if "<!--" in line and "-->" not in line.split("<!--", 1)[1]:
                in_comment = True
            elif in_comment and "-->" in line:
                in_comment = False
            elif not in_comment and re.match(r"^#{1,5} ", line):
                line = "#" + line
        out.append(line)
    return "\n".join(out)


REFS_BLOCK_STANDALONE = "# References {.unnumbered}\n\n::: {#refs}\n:::\n"

APPENDIX_BLOCK = """```{=latex}
% Appendix: \\appendix resets the chapter counter and switches the TOC prefix
% to APPENDIX (bios-prelim.cls). \\setcounter{chapter}{1} makes this chapter's
% appendix APPENDIX B, after Chapter 2's APPENDIX A (plan step 7).
\\appendix
\\addtocontents{toc}{\\protect\\renewcommand{\\protect\\uncTocChapPrefix}{APPENDIX}}
\\setcounter{chapter}{1}
```

# Supplementary Material for Chapter 4
"""

REFS_BLOCK_PRELIM = """```{=latex}
\\chapter*{References}
\\addcontentsline{toc}{chapter}{References}
```

::: {#refs}
:::
"""

HARNESS = """```{=latex}
% Standalone preview harness (as in Chapter 2): delete this block when the
% chapters are assembled into one prelim document. \\setcounter must follow
% \\mainmatter, or the heading silently reads CHAPTER 1.
\\frontmatter
\\renewcommand{\\contentsname}{TABLE OF CONTENTS}
\\tableofcontents
\\mainmatter
\\setcounter{chapter}{2}
```
"""


def transform(chapter_text, header_text, source_label):
    # T1: drop the standalone YAML (author block, \linenumbers, natbib, a4
    #     geometry); keep only its title.
    yaml_text, body = strip_yaml(chapter_text)
    title = yaml_title(yaml_text)

    # T2: remove the standalone References heading + refs div + \clearpage.
    #     The prelim puts References at the very end, after the appendix,
    #     as an unnumbered chapter (Chapter 2's convention; see T6).
    if body.count(REFS_BLOCK_STANDALONE) != 1:
        die("expected exactly one standalone References block:\n" + REFS_BLOCK_STANDALONE)
    body = body.replace(REFS_BLOCK_STANDALONE, "")

    # T3: split at the single '# Appendix' heading.
    #     The appendix is optional: the chapter currently has none.
    parts = re.split(r"^# Appendix\s*$", body, flags=re.M)
    if len(parts) > 2:
        die("expected at most one '# Appendix' heading")
    main, appendix = parts if len(parts) == 2 else (body, None)

    # T4: the chapter title becomes the level-1 (\chapter) heading, so every
    #     body heading moves down one level: '#' -> '##' (3.1), '##' -> '###'.
    main = demote_headings(main)

    # T5: the appendix becomes APPENDIX B with its own chapter-level heading;
    #     '## A1. Title' -> '## Title' (the class numbers it B.1). The
    #     \applabel anchors are redefined as plain \label in the header, so
    #     'Appendix \ref{...}' prints B.1 here and A1 in the standalone.
    if appendix is not None:
        appendix, n_app = re.subn(r"^## A\d+\. ", "## ", appendix, flags=re.M)
        if n_app == 0:
            die("found no '## A<n>. ' appendix subsection headings")

    # T6: assemble. Header template (title substituted), generated-file notice,
    #     preview harness, chapter heading, body, appendix, References.
    header = header_text.replace("{{TITLE}}", title.replace('"', '\\"'))
    notice = ("<!--\nGENERATED FILE -- DO NOT EDIT. Regenerated from\n" + source_label +
              "\nby tools/sync_to_prelim.sh in multiomicsGEP. Edit the multiomicsGEP chapter and\n"
              "commit; the post-commit hook regenerates and commits this copy.\n-->\n")
    out = (header.rstrip("\n") + "\n" + notice + "\n" + HARNESS + "\n# " + title + "\n" +
           main.rstrip("\n") + "\n\n" +
           (APPENDIX_BLOCK + appendix.rstrip("\n") + "\n\n" if appendix is not None else "") +
           REFS_BLOCK_PRELIM)
    return out


# ---------------------------------------------------------------------------
# Step 4: figures
# ---------------------------------------------------------------------------

def sync_figures(draft_text, src_dir, out_dir):
    # Both raw-LaTeX \includegraphics and Markdown ![caption](path) images.
    refs = sorted(set(re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", draft_text)) |
                  set(re.findall(r"!\[[^\]]*\]\((figures/[^)\s]+)\)", draft_text)))
    for r in refs:
        if not r.startswith("figures/"):
            die("figure path outside figures/: " + r)
        if not os.path.isfile(os.path.join(src_dir, r)):
            die("referenced figure not found in the chapter directory: " + r)
    copied = 0
    for r in refs:
        src, dst = os.path.join(src_dir, r), os.path.join(out_dir, r)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if not (os.path.exists(dst) and filecmp.cmp(src, dst, shallow=False)):
            shutil.copy2(src, dst)
            copied += 1
    removed = 0
    fig_root = os.path.join(out_dir, "figures")
    keep = {os.path.normpath(os.path.join(out_dir, r)) for r in refs}
    for dirpath, _, files in os.walk(fig_root, topdown=False):
        for f in files:
            p = os.path.normpath(os.path.join(dirpath, f))
            if p not in keep:
                os.remove(p)
                removed += 1
        if dirpath != fig_root and not os.listdir(dirpath):
            os.rmdir(dirpath)
    return copied, removed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chapter", required=True)
    ap.add_argument("--chapter-bib", required=True)
    ap.add_argument("--master-bib", required=True)
    ap.add_argument("--header", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--source-label", required=True)
    ap.add_argument("--check-only", action="store_true")
    a = ap.parse_args()

    chapter_text = open(a.chapter, encoding="utf-8").read()
    keys = check_bib(chapter_text, a.chapter_bib, a.master_bib)
    if a.check_only:
        print("RESULT bib-check-ok keys=%d" % len(keys))
        return

    draft = transform(chapter_text, open(a.header, encoding="utf-8").read(), a.source_label)
    os.makedirs(a.out_dir, exist_ok=True)
    out_qmd = os.path.join(a.out_dir, "project3-ssbmf-draft.qmd")
    old = open(out_qmd, encoding="utf-8").read() if os.path.exists(out_qmd) else None
    if old != draft:
        with open(out_qmd, "w", encoding="utf-8") as fh:
            fh.write(draft)
    copied, removed = sync_figures(draft, os.path.dirname(os.path.abspath(a.chapter)), a.out_dir)
    print("RESULT qmd=%s figs_copied=%d figs_removed=%d keys=%d"
          % ("unchanged" if old == draft else "changed", copied, removed, len(keys)))


if __name__ == "__main__":
    main()
