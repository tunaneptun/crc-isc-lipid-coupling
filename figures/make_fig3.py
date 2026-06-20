#!/usr/bin/env python3
"""Generate Figure 3: composite decomposition -- retention bar chart."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "core_validation_v2/results/tables/gate3_variant_effects.csv"
OUT_DIR = Path(__file__).resolve().parent

df = pd.read_csv(DATA)

# ── Step 1: verify ───────────────────────────────────────────────────
key_variants = ["baseline", "hmgcs2_removed", "fabp1_removed", "core_beta_oxidation_5"]
for v in key_variants:
    row = df[df["variant"] == v].iloc[0]
    print(f"{v}: retention = {row['retention_ratio']:.4f}")

# ── Step 2: plot ─────────────────────────────────────────────────────
plot_order = ["baseline", "hmgcs2_removed", "fabp1_removed", "core_beta_oxidation_5"]
labels = [
    "Baseline\n(full composite)",
    "HMGCS2\nremoved",
    "FABP1\nremoved",
    "Core beta-oxidation\n(5 genes only)",
]

retentions = []
for v in plot_order:
    retentions.append(df.loc[df["variant"] == v, "retention_ratio"].values[0])

threshold = 0.50
clr_above = "#5B7B94"   # slate
clr_below = "#C47A3A"   # muted orange
colors = [clr_above if r >= threshold else clr_below for r in retentions]

fig, ax = plt.subplots(figsize=(6, 4.5), dpi=300)

x = range(len(plot_order))
bars = ax.bar(x, retentions, color=colors, width=0.6, edgecolor="white", linewidth=0.5)

# Value labels on bars
for xi, r in zip(x, retentions):
    ax.text(xi, r + 0.02, f"{r:.2f}", ha="center", va="bottom", fontsize=10, fontweight="bold",
            color="#333333")

# Threshold line
ax.axhline(threshold, color="#999999", ls="--", lw=1.2, zorder=0)
ax.text(len(plot_order) - 0.55, threshold + 0.02, "50% retention threshold",
        ha="right", va="bottom", fontsize=8.5, color="#777777")

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=9.5)
ax.set_ylabel("Effect retained (fraction of baseline shift)", fontsize=11)
ax.set_ylim(0, 1.15)
ax.tick_params(axis="y", labelsize=10)

ax.set_title("What carries the shift: composite decomposition", fontsize=13, pad=12)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
fig.savefig(OUT_DIR / "fig3_decomposition.png", dpi=300,
            bbox_inches="tight", pad_inches=0.15, facecolor="white")
fig.savefig(OUT_DIR / "fig3_decomposition.svg",
            bbox_inches="tight", pad_inches=0.15, facecolor="white")
plt.close()

print(f"\nSaved: {OUT_DIR / 'fig3_decomposition.png'}")
print(f"Saved: {OUT_DIR / 'fig3_decomposition.svg'}")
