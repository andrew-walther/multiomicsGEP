"""Build standalone SVG diagrams for the project tracker from the repository root."""

from html import escape
from pathlib import Path

OUTPUT = Path(__file__).resolve().parent


def diagram(filename, title, panels):
    """Write one SVG with vertically stacked panels of connected model steps."""
    matrix_intro = filename == "supervision-states.svg"
    top = 205 if matrix_intro else 60
    width, height = 1140, top + 8 + 178 * len(panels)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img">',
           f'<title>{escape(title)}</title>',
           '<rect width="100%" height="100%" rx="12" fill="#f5f8fc"/>',
           '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0L10 5L0 10Z" fill="#395775"/></marker></defs>']

    def label(x, y, value, size=17, weight=400, anchor="start"):
        """Append escaped text using the diagram’s shared typography."""
        svg.append(f'<text x="{x}" y="{y}" font-family="system-ui,sans-serif" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="#102a43">{escape(value)}</text>')

    label(24, 37, title, 25, 700)
    if matrix_intro:
        for x, rows, cols, color, caption in [
            (64, 5, 10, "#c6d7e8", "Y: patients × genes"),
            (277, 5, 3, "#91caa9", "L: patient scores"),
            (458, 3, 10, "#9dbfe6", "Fᵀ: gene loadings"),
            (721, 5, 2, "#f3d084", "(t, δ): same patients"),
        ]:
            for row in range(rows):
                for col in range(cols):
                    opacity = 0.48 + 0.17 * ((row + 2 * col) % 4)
                    svg.append(f'<rect x="{x+col*11}" y="{72+row*11}" width="10" height="10" fill="{color}" opacity="{opacity}"/>')
            label(x, 153, caption, 15, 700)
        label(220, 107, "≈", 28)
        label(392, 107, "×", 26)
        label(620, 104, "paired with", 16)
        label(64, 181, "Rows represent patients. Each patient contributes a molecular profile and an observed time and event indicator.", 16)
    for i, (heading, boxes, footer) in enumerate(panels):
        y = top + i * 178
        svg.append(f'<rect x="18" y="{y}" width="1104" height="162" rx="12" fill="white" stroke="#b4c7da"/>')
        label(34, y + 27, heading, 19, 700)
        for j, (name, line1, line2) in enumerate(boxes):
            x = 34 + j * 271
            svg.append(f'<rect x="{x}" y="{y+42}" width="244" height="80" rx="8" fill="{("#e7f1fc", "#e7f5ed", "#fff2d9", "#eeeafa")[j]}"/>')
            label(x + 122, y + 64, name, 17, 700, "middle")
            label(x + 122, y + 87, line1, 15, anchor="middle")
            label(x + 122, y + 108, line2, 15, anchor="middle")
            if j == 2 and heading.startswith("A "):
                svg.append(f'<path d="M{x+268} {y+82}h-21" stroke="#395775" stroke-width="2" marker-end="url(#arrow)"/>')
            elif j < len(boxes) - 1:
                svg.append(f'<path d="M{x+246} {y+82}h21" stroke="#395775" stroke-width="2" marker-end="url(#arrow)"/>')
        if matrix_intro and heading.startswith("B "):
            svg.append(f'<path d="M970 {y+123}v7H430v-7" fill="none" stroke="#395775" stroke-width="2" stroke-dasharray="5 4" marker-end="url(#arrow)"/>')
        label(34, y + 147, footer, 15, 600)
    (OUTPUT / filename).write_text("\n".join(svg + ["</svg>"]) + "\n", encoding="utf-8")


diagram("supervision-states.svg", "From a predictive benchmark to survival-guided programs", [
    ("A  Historical single-modality YFB benchmark", [
        ("Expression", "Y · n patients × p genes", "Training molecular profiles"),
        ("EBMF reconstruction", "Y ≈ LFᵀ", "L: n × K   F: p × K"),
        ("Cox predictor", "η = Zβ, Z from Y and F", "β: K survival coefficients"),
        ("Survival observations", "(tᵢ, δᵢ) · n paired outcomes", "Cox risk sets inform β"),
    ], "Feedback into F is off (α_F = 0). The same patients contribute molecular profiles and survival observations."),
    ("B  Target single-modality method", [
        ("Molecular term", "Fit reconstruction LFᵀ", "Gene-wise noise precision τ"),
        ("Prior and posterior", "Estimate moments of L, F, β", "Carry posterior variances"),
        ("Cox term", "YFB: η = Z(F)β", "LB alternative: η = Lβ"),
        ("Supervised updates", "YFB: feedback to F", "LB: direct feedback to L"),
    ], "Goal: the survival term changes the molecular representation while preserving recovery, calibration and convergence."),
    ("C  Multimodal extension with modality-specific priors", [
        ("Paired modalities", "Y¹: expression, n × p₁", "Y²: methylation, n × p₂"),
        ("Shared patient factors", "Yᵐ ≈ L(Fᵐ)ᵀ", "One L across modalities"),
        ("Modality loadings", "Fᵐ: pₘ × K", "Prior and noise by modality"),
        ("Joint survival model", "Combine modality projections", "One risk score per patient"),
    ], "M = 1 is the single-modality case. Evaluate the same supervision gate before adopting a unified fitter."),
])
diagram("variational-updates.svg", "How the objective becomes a fitting algorithm", [
    ("Objective components", [
        ("Molecular likelihood", "Observed Y given L, F, τ", "Residuals drive reconstruction"),
        ("Cox partial likelihood", "Time, event and risk sets", "η links molecular and survival"),
        ("Empirical Bayes priors", "Separate priors for L, F, β", "Shrink noisy coordinates"),
        ("Variational objective", "Expected log terms and KL", "Working Cox approximation"),
    ], "The molecular and Cox weights depend on the model variant. S-9 is a planned change to the genomic weighting."),
    ("One sweep of the experimental YFB fitter", [
        ("Refresh Cox terms", "Compute η and risk-set terms", "Form local working quantities"),
        ("Update L and F", "Use residuals and priors", "Include survival terms in F"),
        ("Update β and τ", "Use projected-score moments", "Update noise or hold it fixed"),
        ("Check the fit", "Track objective and movement", "Record explicit failures"),
    ], "Repeat until the stopping criteria are met or the fit reaches a recorded cap or failure."),
])
diagram("held-out-prediction.svg", "Evaluate frozen YFB predictions on new patients", [
    ("For a fixed fit, held-out outcomes evaluate its predictions", [
        ("Align molecular inputs", "Apply benchmark preprocessing", "Match held-out genes to F"),
        ("Project new patients", "Z* = Y*F̃", "F̃ uses training norms"),
        ("Score risk", "η* = Z*β̂", "Use frozen coefficients"),
        ("Evaluate discrimination", "Compare η* with (t*, δ*)", "Compute cohort C-index"),
    ], "Held-out outcomes enter evaluation. Use the training-fixed risk orientation for concordance."),
])

def comparison():
    """Draw separate data paths for the sequential and joint YFB fits."""
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1140 506" role="img">',
           '<title>Two-step versus joint EBMF and Cox fitting</title>',
           '<defs><marker id="tip" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0L10 5L0 10Z" fill="#395775"/></marker></defs>',
           '<rect width="1140" height="506" rx="14" fill="#f5f8fc"/>']

    def label(x, y, text, size=17, weight=400):
        """Add a centered diagram label."""
        svg.append(f'<text x="{x}" y="{y}" text-anchor="middle" font-family="system-ui,sans-serif" font-size="{size}" font-weight="{weight}" fill="#102a43">{escape(text)}</text>')

    def box(x, y, width, height, fill):
        """Draw a node or panel background."""
        svg.append(f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="10" fill="{fill}" stroke="#b4c7da"/>')

    def tiles(x, y, rows, columns, color):
        """Draw a schematic matrix with no underlying patient data."""
        for i in range(rows):
            for j in range(columns):
                opacity = 0.50 + 0.15 * ((i + 2 * j) % 4)
                svg.append(f'<rect x="{x+12*j}" y="{y+12*i}" width="11" height="11" fill="{color}" opacity="{opacity}"/>')

    def arrow(path, color="#395775"):
        """Connect the output of one step to the input of another."""
        svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2.5" marker-end="url(#tip)"/>')

    box(18, 16, 1104, 221, "white")
    label(570, 46, "A  Two-step fit", 23, 700)
    tiles(59, 83, 5, 7, "#668db8")
    label(100, 164, "Expression Y")
    box(197, 77, 168, 81, "#e7f1fc")
    label(281, 110, "Fit EBMF", 19, 700)
    label(281, 137, "Y ≈ LFᵀ")
    tiles(425, 91, 4, 3, "#328456")
    tiles(477, 91, 3, 6, "#668db8")
    label(484, 164, "Freeze L and F", 17, 700)
    box(634, 77, 173, 81, "#fff2d9")
    label(720, 109, "Fit Cox", 19, 700)
    label(720, 137, "Use fixed scores")
    tiles(933, 97, 4, 1, "#ba7a23")
    label(939, 165, "Coefficients β")
    for path in ("M148 115H192", "M369 115H412", "M557 115H629", "M812 115H916"):
        arrow(path)
    label(720, 207, "Time and event (t, δ)", 16, 600)
    arrow("M720 185V162", "#9c631b")
    label(306, 207, "Survival outcomes enter only at the Cox step.", 16)

    box(18, 251, 1104, 238, "white")
    label(570, 281, "B  Proposed joint YFB fit", 23, 700)
    tiles(59, 306, 4, 7, "#668db8")
    label(100, 374, "Expression Y")
    tiles(85, 398, 3, 2, "#ba7a23")
    label(100, 454, "Time and event (t, δ)", 16)
    box(221, 305, 204, 67, "#e7f1fc")
    label(323, 331, "EBMF reconstruction", 17, 700)
    label(323, 355, "Y ≈ LFᵀ")
    box(221, 389, 204, 67, "#fff2d9")
    label(323, 415, "Cox survival term", 17, 700)
    label(323, 440, "Risk from Y, F and β")
    arrow("M148 337H216")
    arrow("M169 414H216", "#9c631b")
    box(546, 325, 248, 111, "#e7f5ed")
    label(670, 358, "Fit both terms together", 19, 700)
    label(670, 384, "Update L, F and β")
    label(670, 408, "Cox information enters the F update", 13)
    arrow("M430 337H482V355H541")
    arrow("M430 422H482V405H541", "#9c631b")
    arrow("M799 380H867")
    tiles(883, 352, 4, 3, "#328456")
    tiles(934, 352, 3, 6, "#668db8")
    tiles(1023, 352, 4, 1, "#ba7a23")
    label(955, 428, "Joint estimates of L, F and β")
    svg.append('<path d="M950 438V449H670V441" fill="none" stroke="#395775" stroke-width="2" stroke-dasharray="5 4" marker-end="url(#tip)"/>')
    label(658, 473, "The molecular and survival terms share F throughout fitting.", 16)
    (OUTPUT / "two-step-vs-joint.svg").write_text("\n".join(svg + ["</svg>"]) + "\n", encoding="utf-8")


comparison()
