#!/usr/bin/env python3
"""Compose multi-panel main figures from existing single-panel PNGs."""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

ROOT = Path("/mnt/sda1/data/TNT/analysis")
FIG = ROOT / "figures"
OUT = FIG / "panels"
OUT.mkdir(exist_ok=True)

def panel(name, grid, images, figsize, title=None):
    rows, cols = grid
    fig, axes = plt.subplots(rows, cols, figsize=figsize, dpi=200)
    axes = axes.flatten() if rows*cols > 1 else [axes]
    for ax, (label, path) in zip(axes, images):
        ax.axis("off")
        if path and Path(path).exists():
            ax.imshow(mpimg.imread(path))
        ax.set_title(label, fontsize=14, fontweight="bold", loc="left", x=-0.02, y=1.0)
    for ax in axes[len(images):]:
        ax.axis("off")
    if title:
        fig.suptitle(title, fontsize=15, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {name}")

# Fig 2: WES landscape
panel("Fig2_WES_landscape", (1,2),
      [("A. TMB by response", FIG/"Fig_TMB_by_response.png"),
       ("B. Driver oncoprint (pre)", FIG/"Fig_driver_oncoprint_pre.png")],
      figsize=(14,6),
      title="Figure 2. Somatic landscape of MSS LARC TNT cohort")

# Fig 3: RNA immune / signature
panel("Fig3_RNA_signatures", (1,2),
      [("A. Pre-treatment signature boxplots", FIG/"Fig_pre_signatures_boxplot.png"),
       ("B. Pre-treatment signature heatmap", FIG/"Fig_pre_signature_heatmap.png")],
      figsize=(16,7),
      title="Figure 3. Pre-treatment immune and pathway signatures")

# Fig 4: Main result — GSEA + integrated features
panel("Fig4_GSEA_integration", (1,2),
      [("A. Hallmark GSEA (pre, good vs bad)", FIG/"Fig_GSEA_Hallmark.png"),
       ("B. Feature-response associations", FIG/"Fig_response_features_barh.png")],
      figsize=(16,7),
      title="Figure 4. DNA repair / cell-cycle programs predict TNT response")

# Fig 5: Integration + paired delta
panel("Fig5_integration_delta", (2,2),
      [("A. Feature correlation", FIG/"Fig_feature_correlation_heatmap.png"),
       ("B. Paired Δ ssGSEA", FIG/"Fig_paired_delta_ssgsea.png"),
       ("C. Paired Δ TCR/BCR (TRUST4)", FIG/"Fig_paired_delta_TRUST4.png"),
       ("D. placeholder", None)],
      figsize=(16,13),
      title="Figure 5. Feature integration and treatment-induced changes")

print("done")
