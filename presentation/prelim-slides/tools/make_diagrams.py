"""Draw the Part III deck's schematic diagrams as PNG files.

Script: make_diagrams.py
Purpose: Plain-label versions of the progress-book diagrams
         (docs/progress_book/figures/build_diagrams.py) for the prelim deck:
         no in-figure titles and no project shorthand. Each diagram is written
         as SVG and rendered to PNG with headless Google Chrome.
Author: Andrew Walther
Created: 2026-11-03
Dependencies: Python 3 standard library; Google Chrome (headless) for rendering

Run from the multiomicsGEP root:
    python3 presentation/prelim-slides/tools/make_diagrams.py
Writes presentation/prelim-slides/figures/diagram_*.png. The SVG sources are
kept in a temporary directory only. No data are read.
"""
import subprocess
import sys
import tempfile
from html import escape
from pathlib import Path

FIG = Path(__file__).resolve().parents[1] / "figures"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
INK, ARROW, SURV = "#102a43", "#395775", "#9c631b"
FILLS = ("#e7f1fc", "#e7f5ed", "#fff2d9", "#eeeafa", "#fde8e8")


class Canvas:
    """Collect SVG elements with the shared typography and arrow style."""

    def __init__(self, width, height):
        self.width, self.height = width, height
        self.svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
                    f'viewBox="0 0 {width} {height}">',
                    '<defs><marker id="tip" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
                    'markerHeight="7" orient="auto"><path d="M0 0L10 5L0 10Z" fill="context-stroke"/></marker></defs>',
                    f'<rect width="{width}" height="{height}" fill="white"/>']

    def text(self, x, y, value, size=20, weight=400, anchor="middle", color=INK):
        """Add one line of escaped text."""
        self.svg.append(f'<text x="{x}" y="{y}" font-family="Helvetica, Arial, sans-serif" font-size="{size}" '
                        f'font-weight="{weight}" text-anchor="{anchor}" fill="{color}">{escape(value)}</text>')

    def box(self, x, y, w, h, fill, stroke="#b4c7da"):
        """Add a rounded node."""
        self.svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{fill}" stroke="{stroke}"/>')

    def node(self, x, y, w, h, fill, lines):
        """Add a node whose first line is bold and the rest regular, vertically centred."""
        self.box(x, y, w, h, fill)
        top = y + h / 2 - 13 * (len(lines) - 1) + 7
        for i, line in enumerate(lines):
            self.text(x + w / 2, top + 26 * i, line, 21 if i == 0 else 18, 700 if i == 0 else 400)

    def arrow(self, path, color=ARROW, dashed=False):
        """Add a connector path ending in an arrowhead."""
        dash = ' stroke-dasharray="7 5"' if dashed else ''
        self.svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2.5"{dash} marker-end="url(#tip)"/>')

    def tiles(self, x, y, rows, cols, color, size=14):
        """Draw a schematic matrix (no data)."""
        for i in range(rows):
            for j in range(cols):
                op = 0.45 + 0.17 * ((i + 2 * j) % 4)
                self.svg.append(f'<rect x="{x + (size + 1) * j}" y="{y + (size + 1) * i}" width="{size}" '
                                f'height="{size}" fill="{color}" opacity="{op:.2f}"/>')

    def write(self, name, tmp):
        """Save the SVG and render it to a PNG at twice the nominal size."""
        src = Path(tmp) / f"{name}.svg"
        src.write_text("\n".join(self.svg + ["</svg>"]), encoding="utf-8")
        out = FIG / f"{name}.png"
        cmd = [CHROME, "--headless=new", "--disable-gpu", "--hide-scrollbars", "--force-device-scale-factor=2",
               f"--window-size={self.width},{self.height}", f"--screenshot={out}", src.as_uri()]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0 or not out.exists():
            sys.exit(f"Chrome failed to render {name}: {res.stderr[-400:]}")
        print("Wrote", out)


def flow(name, tmp, steps, loop=None, note=None):
    """A left-to-right row of nodes, with an optional dashed loop from the last back to the first."""
    w, gap, h = 330, 60, 150
    width = 40 + len(steps) * w + (len(steps) - 1) * gap
    c = Canvas(width, 300 if (loop or note) else 230)
    for i, lines in enumerate(steps):
        x = 20 + i * (w + gap)
        c.node(x, 40, w, h, FILLS[i % len(FILLS)], lines)
        if i < len(steps) - 1:
            c.arrow(f"M{x + w + 4} {40 + h / 2}h{gap - 10}")
    if loop:
        last = 20 + (len(steps) - 1) * (w + gap) + w / 2
        c.arrow(f"M{last} {40 + h + 4}v38H{20 + w / 2}v-34", dashed=True)
        c.text(width / 2, 40 + h + 34, loop, 19, 600)
    elif note:
        c.text(width / 2, 40 + h + 60, note, 19, 600)
    c.write(name, tmp)


def two_step_vs_joint(tmp):
    """Two-step analysis (top) and joint fit (bottom), with survival paths in brown."""
    c = Canvas(1240, 560)
    c.box(10, 10, 1220, 240, "#f7f9fc")
    c.text(30, 44, "Two-step analysis", 24, 700, "start")
    c.tiles(40, 90, 5, 7, "#668db8")
    c.text(95, 190, "Expression Y", 19)
    c.node(200, 80, 220, 90, FILLS[0], ["Factorize", "Y ≈ LFᵀ"])
    c.tiles(470, 95, 4, 3, "#328456")
    c.tiles(525, 102, 3, 6, "#668db8")
    c.text(550, 190, "Fixed L and F", 19, 700)
    c.node(680, 80, 220, 90, FILLS[2], ["Cox model", "on fixed scores"])
    c.tiles(1010, 95, 4, 1, "#ba7a23")
    c.text(1017, 190, "Coefficients β", 19)
    for p in ("M160 125H195", "M424 125H462", "M630 125H675", "M904 125H1000"):
        c.arrow(p)
    c.text(790, 232, "Survival (t, δ) enters only here", 18, 600, color=SURV)
    c.arrow("M790 212V176", SURV)
    c.box(10, 270, 1220, 280, "#f7f9fc")
    c.text(30, 304, "Joint model", 24, 700, "start")
    c.tiles(40, 330, 4, 7, "#668db8")
    c.text(95, 412, "Expression Y", 19)
    c.tiles(80, 440, 3, 2, "#ba7a23")
    c.text(95, 505, "Survival (t, δ)", 19)
    c.node(220, 325, 240, 80, FILLS[0], ["Expression term", "Y ≈ LFᵀ"])
    c.node(220, 430, 240, 80, FILLS[2], ["Survival term", "risk from Y, F and β"])
    c.arrow("M160 365H215")
    c.arrow("M135 470H215", SURV)
    c.node(560, 360, 280, 120, FILLS[1], ["Fit together", "update L, F and β", "F appears in both terms"])
    c.arrow("M464 365H510V400H555")
    c.arrow("M464 470H510V440H555", SURV)
    c.arrow("M844 420H905")
    c.tiles(925, 390, 4, 3, "#328456")
    c.tiles(980, 397, 3, 6, "#668db8")
    c.tiles(1085, 390, 4, 1, "#ba7a23")
    c.text(1020, 480, "L, F and β estimated jointly", 19)
    c.write("diagram_two_step_vs_joint", tmp)


def notation(tmp):
    """Matrix sketch: Y ≈ L × Fᵀ, paired with each patient's survival outcome."""
    c = Canvas(1300, 330)
    items = [(30, 8, 12, "#668db8", "Y: n patients × p genes"),
             (360, 8, 3, "#328456", "L: n × K activities"),
             (600, 3, 12, "#9dbfe6", "Fᵀ: K × p loadings"),
             (1130, 8, 2, "#ba7a23", "(t, δ): time, event")]
    for x, rows, cols, color, cap in items:
        c.tiles(x, 40 if rows == 8 else 85, rows, cols, color, 18)
        c.text(x, 225, cap, 20, 700, "start")
    c.text(310, 120, "≈", 44)
    c.text(530, 120, "×", 40)
    c.text(1030, 125, "paired with", 20)
    c.text(650, 285, "One row per patient. A column of F is a program: its weight on every gene.", 20)
    c.write("diagram_notation", tmp)


def multi_data(tmp):
    """Several data types share patient activities L; each has its own loadings and prior."""
    c = Canvas(1240, 420)
    c.node(20, 40, 300, 110, FILLS[0], ["Expression Y⁽¹⁾", "n × p₁ genes"])
    c.node(20, 230, 300, 110, FILLS[3], ["Methylation Y⁽²⁾", "n × p₂ CpG sites"])
    c.node(470, 135, 300, 110, FILLS[1], ["Shared activities L", "n × K, one row per patient"])
    c.node(920, 20, 300, 110, FILLS[0], ["Loadings F⁽¹⁾", "prior family for genes"])
    c.node(920, 150, 300, 110, FILLS[3], ["Loadings F⁽²⁾", "sparser prior for CpG sites"])
    c.node(920, 280, 300, 120, FILLS[2], ["One risk score", "η = Σₘ Y⁽ᵐ⁾F⁽ᵐ⁾D⁽ᵐ⁾ β"])
    c.arrow("M324 95H400V175H465")
    c.arrow("M324 285H400V205H465")
    c.arrow("M774 170H850V75H915")
    c.arrow("M774 205H915")
    c.arrow("M774 225H850V340H915", SURV)
    c.write("diagram_multi_data", tmp)


def main():
    """Draw every diagram."""
    if not Path(CHROME).exists():
        sys.exit(f"Google Chrome not found at {CHROME}")
    FIG.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        two_step_vs_joint(tmp)
        notation(tmp)
        multi_data(tmp)
        flow("diagram_fitting", tmp, [
            ["Approximate the Cox term", "weighted least squares", "around the current risk"],
            ["Update each program k", "βₖ, then Lₖ, then Fₖ", "(empirical Bayes)"],
            ["Update noise", "one precision per gene"],
            ["Check convergence", "stop or repeat"]],
            loop="Repeat each sweep. Start with extra programs and let the priors remove the unneeded ones.")
        flow("diagram_prediction", tmp, [
            ["New cohort", "keep training genes,", "standardize within cohort"],
            ["Project", "Z = Y F̂ D̂"],
            ["Score", "η = Z β̂", "training direction"],
            ["Evaluate", "C-index against (t, δ)"]],
            note="Nothing about the new patient is estimated.")


if __name__ == "__main__":
    main()
