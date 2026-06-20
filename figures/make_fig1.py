#!/usr/bin/env python3
"""Generate Figure 1: paired slopegraph of per-patient ISC-lipid/FAO coupling."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "core_validation_v2/results/tables/gate4_per_patient.csv"
OUT_DIR = Path(__file__).resolve().parent

df = pd.read_csv(DATA)
df = df[(df["method"] == "scanpy") & (df["retained_for_paired_analysis"] == True)]

r_normal = df["r_normal"].values
r_tumor = df["r_tumor"].values
n = len(df)

mean_n = r_normal.mean()
mean_t = r_tumor.mean()
n_increase = (r_tumor > r_normal).sum()   # shift toward zero (majority)
n_decrease = (r_tumor < r_normal).sum()   # exception

# ── Verification ─────────────────────────────────────────────────────
print(f"n_patients:     {n}")
print(f"mean r_normal:  {mean_n:.6f}")
print(f"mean r_tumor:   {mean_t:.6f}")
print(f"tumor > normal: {n_increase}  (shift toward zero)")
print(f"tumor < normal: {n_decrease}  (exception)")

# ── Figure ───────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 6), dpi=300)

x_left, x_right = 0, 1
clr_majority = "#888888"
clr_exception = "#C47A3A"

for rn, rt in zip(r_normal, r_tumor):
    if rt >= rn:  # majority: coupling moves toward zero
        ax.plot([x_left, x_right], [rn, rt],
                color=clr_majority, lw=0.8, alpha=0.4, zorder=1)
    else:
        ax.plot([x_left, x_right], [rn, rt],
                color=clr_exception, lw=0.8, alpha=0.7, zorder=2)

# Mean trajectory
ax.plot([x_left, x_right], [mean_n, mean_t],
        color="#1a1a2e", lw=2.5, zorder=5, solid_capstyle="round")
ax.plot(x_left, mean_n, "o", color="#1a1a2e", ms=7, zorder=6)
ax.plot(x_right, mean_t, "o", color="#1a1a2e", ms=7, zorder=6)

ax.annotate(f"mean {mean_n:.2f}", xy=(x_left, mean_n),
            xytext=(-0.12, mean_n), textcoords="data",
            fontsize=9.5, color="#1a1a2e", fontweight="bold",
            ha="right", va="center")
ax.annotate(f"mean {mean_t:.2f}", xy=(x_right, mean_t),
            xytext=(1.12, mean_t), textcoords="data",
            fontsize=9.5, color="#1a1a2e", fontweight="bold",
            ha="left", va="center")

# Reference line at y = 0
ax.axhline(0, color="#cccccc", ls=":", lw=1, zorder=0)

# Axes
ax.set_xlim(-0.6, 1.6)
ax.set_xticks([x_left, x_right])
ax.set_xticklabels(["Adjacent-normal", "Tumor"], fontsize=11)
ax.set_ylabel("Per-patient correlation r  (ISC score vs lipid/FAO score)",
              fontsize=11)
ax.tick_params(axis="y", labelsize=10)

# Spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Title / subtitle
fig.suptitle("ISC\u2013lipid/FAO coupling per donor: tumor vs adjacent-normal",
             fontsize=13, y=0.98)
fig.text(0.5, 0.925, "Pelka 2021 cE01 Stem/TA-like cells  \u00b7  n = 34 paired donors",
         ha="center", fontsize=10, color="grey")

# Legend
from matplotlib.lines import Line2D
handles = [
    Line2D([0], [0], color=clr_majority, lw=1, alpha=0.6,
           label=f"Coupling less negative in tumor ({n_increase} patients)"),
    Line2D([0], [0], color=clr_exception, lw=1, alpha=0.8,
           label=f"Coupling more negative in tumor ({n_decrease} patient)"),
]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14),
          frameon=False, ncol=1, fontsize=9, handlelength=1.8)

plt.tight_layout()
fig.savefig(OUT_DIR / "fig1_per_patient_coupling.png", dpi=300,
            bbox_inches="tight", pad_inches=0.15, facecolor="white")
fig.savefig(OUT_DIR / "fig1_per_patient_coupling.svg",
            bbox_inches="tight", pad_inches=0.15, facecolor="white")
plt.close()

print(f"\nSaved: {OUT_DIR / 'fig1_per_patient_coupling.png'}")
print(f"Saved: {OUT_DIR / 'fig1_per_patient_coupling.svg'}")
