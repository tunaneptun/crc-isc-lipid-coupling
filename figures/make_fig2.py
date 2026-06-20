#!/usr/bin/env python3
"""Generate Figure 2: ISC vs lipid/FAO hexbin scatter, normal vs tumor."""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "core_validation_v2/results/gate2b_clean_background_scores_cE01.parquet"
OUT_DIR = Path(__file__).resolve().parent

df = pd.read_parquet(DATA)

# ── Step 1: verify ───────────────────────────────────────────────────
print("Columns:", list(df.columns))
print(f"Shape:   {df.shape}")
print()

panels = [
    ("N", "Adjacent-normal"),
    ("T", "Tumor"),
]

for code, label in panels:
    mask = df["SPECIMEN_TYPE"] == code
    n = mask.sum()
    isc = df.loc[mask, "ISC_scanpy_cleanbg"].values
    fao = df.loc[mask, "FAO_scanpy_cleanbg"].values
    r, _ = stats.pearsonr(isc, fao)
    print(f"{label} ({code}): n={n:,}, Pearson r={r:.6f}")

# ── Step 2: plot ─────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 5), dpi=300, sharey=True)

# Shared limits across panels
isc_all = df["ISC_scanpy_cleanbg"].values
fao_all = df["FAO_scanpy_cleanbg"].values
pad = 0.05
xlim = (np.percentile(isc_all, 0.5) - pad, np.percentile(isc_all, 99.5) + pad)
ylim = (np.percentile(fao_all, 0.5) - pad, np.percentile(fao_all, 99.5) + pad)

hb_collections = []

for ax, (code, label) in zip(axes, panels):
    mask = df["SPECIMEN_TYPE"] == code
    isc = df.loc[mask, "ISC_scanpy_cleanbg"].values
    fao = df.loc[mask, "FAO_scanpy_cleanbg"].values
    n = len(isc)
    r, _ = stats.pearsonr(isc, fao)

    hb = ax.hexbin(isc, fao, gridsize=50, cmap="viridis",
                   norm=mcolors.LogNorm(), mincnt=1,
                   extent=[xlim[0], xlim[1], ylim[0], ylim[1]])
    hb_collections.append(hb)

    # Linear fit
    m, b = np.polyfit(isc, fao, 1)
    x_fit = np.linspace(xlim[0], xlim[1], 100)
    ax.plot(x_fit, m * x_fit + b, color="crimson", lw=1.8, zorder=5)

    ax.set_title(f"{label}  (r = {r:.2f})", fontsize=12)
    ax.set_xlabel("ISC score", fontsize=11)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.tick_params(labelsize=10)

    ax.text(0.97, 0.03, f"n = {n:,}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9.5, color="#444444")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes[0].set_ylabel("lipid/FAO score", fontsize=11)

# Shared colorbar
fig.subplots_adjust(right=0.88)
cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
cb = fig.colorbar(hb_collections[0], cax=cbar_ax)
cb.set_label("cells", fontsize=10)
cb.ax.tick_params(labelsize=9)

fig.subplots_adjust(left=0.08, right=0.87, wspace=0.08)
fig.savefig(OUT_DIR / "fig2_isc_fao_scatter.png", dpi=300,
            bbox_inches="tight", pad_inches=0.15, facecolor="white")
fig.savefig(OUT_DIR / "fig2_isc_fao_scatter.svg",
            bbox_inches="tight", pad_inches=0.15, facecolor="white")
plt.close()

print(f"\nSaved: {OUT_DIR / 'fig2_isc_fao_scatter.png'}")
print(f"Saved: {OUT_DIR / 'fig2_isc_fao_scatter.svg'}")
