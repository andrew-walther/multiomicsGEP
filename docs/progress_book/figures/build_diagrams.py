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
    ("Training ends before held-out survival is used", [
        ("Align molecular inputs", "Apply benchmark preprocessing", "Match held-out genes to F"),
        ("Project new patients", "Z* = Y*F̃", "F̃ uses training norms"),
        ("Score risk", "η* = Z*β̂", "Use frozen coefficients"),
        ("Evaluate discrimination", "Compare η* with (t*, δ*)", "Compute cohort C-index"),
    ], "Held-out outcomes enter evaluation. Use the training-fixed risk orientation for concordance."),
])
