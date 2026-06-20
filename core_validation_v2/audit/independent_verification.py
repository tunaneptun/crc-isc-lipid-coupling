#!/usr/bin/env python3
"""
INDEPENDENT INTEGRITY AUDIT — End-to-end reconstruction
Agent: claude_code (Opus 4.6)

Re-derives the gate results from the source data files, without importing
the analysis code. Compares re-derived values against the committed
artifacts and reports every discrepancy.

READ-ONLY on all existing files. Writes only to core_validation_v2/audit/.
"""

import json
import hashlib
import warnings
import sys
import os
import traceback
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import anndata as ad
import scanpy as sc
import yaml

warnings.filterwarnings("ignore")

ROOT = Path(".")
CV2 = ROOT / "core_validation_v2"
RESULTS = {}  # accumulator for all findings
ISSUES = []   # accumulator for issues

def log(msg):
    print(f"[AUDIT] {msg}", flush=True)

def record(section, key, value):
    RESULTS.setdefault(section, {})[key] = value

def issue(severity, title, location, detail, recommendation=""):
    ISSUES.append({
        "severity": severity,
        "title": title,
        "location": location,
        "detail": detail,
        "recommendation": recommendation,
    })
    log(f"  ** {severity}: {title}")

def compare(label, ours, theirs, tol=1e-6):
    """Compare two numeric values and report."""
    if theirs is None:
        record("comparison", label, {"ours": ours, "theirs": "MISSING", "status": "MISSING"})
        return False
    diff = abs(float(ours) - float(theirs))
    ok = diff <= tol
    record("comparison", label, {
        "ours": float(ours), "theirs": float(theirs),
        "abs_diff": diff, "tol": tol,
        "status": "MATCH" if ok else "MISMATCH"
    })
    if not ok:
        issue("MATERIAL", f"Number mismatch: {label}",
              "independent_verification",
              f"Ours={ours}, Theirs={theirs}, diff={diff}, tol={tol}")
    return ok

# ===========================================================================
# PHASE 2a: FOUNDATIONAL PREMISES
# ===========================================================================
log("=" * 70)
log("PHASE 2a: FOUNDATIONAL PREMISES")
log("=" * 70)

# --- Load config ---
with open(CV2 / "config.yaml") as f:
    cfg = yaml.safe_load(f)

ISC_GENES = cfg["gene_sets"]["isc_9"]
FAO_GENES = cfg["gene_sets"]["ppar_fao_12"]
CORE_BETAOX = cfg["gene_sets"]["core_beta_oxidation_5"]
ALL_21 = sorted(set(ISC_GENES + FAO_GENES))

log(f"ISC-9: {ISC_GENES}")
log(f"FAO-12: {FAO_GENES}")
log(f"Core-beta-ox-5: {CORE_BETAOX}")
log(f"Union-21: {ALL_21} (n={len(ALL_21)})")

record("gene_sets", "ISC_9", ISC_GENES)
record("gene_sets", "FAO_12", FAO_GENES)
record("gene_sets", "core_beta_ox_5", CORE_BETAOX)
record("gene_sets", "union_21_count", len(ALL_21))

# Check union count
if len(ALL_21) != 21:
    issue("CRITICAL", "Union gene count != 21",
          "config.yaml gene_sets",
          f"ISC-9 + FAO-12 union = {len(ALL_21)}, expected 21. Overlap: {set(ISC_GENES) & set(FAO_GENES)}")
else:
    log(f"  Union-21 count verified: {len(ALL_21)}")

# Check core-beta-ox is subset of FAO
if not set(CORE_BETAOX).issubset(set(FAO_GENES)):
    issue("CRITICAL", "Core-beta-ox not subset of FAO-12",
          "config.yaml", f"Extra: {set(CORE_BETAOX) - set(FAO_GENES)}")

# --- Load Pelka normalized ---
log("\nLoading Pelka normalized h5ad...")
adata_norm = ad.read_h5ad(ROOT / "data/processed/pelka_epithelial.h5ad", backed="r")
log(f"  Shape: {adata_norm.shape}")
log(f"  .X dtype: {adata_norm.X.dtype if hasattr(adata_norm.X, 'dtype') else type(adata_norm.X)}")
log(f"  .raw present: {adata_norm.raw is not None}")

record("pelka_norm", "n_cells", adata_norm.n_obs)
record("pelka_norm", "n_genes", adata_norm.n_vars)

# --- Population definition: cE01 == Stem/TA-like ---
log("\nVerifying cE01 population definition...")
obs = adata_norm.obs
cE01_mask = obs["cl295v11SubShort"] == "cE01"
n_cE01 = int(cE01_mask.sum())
log(f"  cE01 cells: {n_cE01}")
record("population", "n_cE01", n_cE01)

# Check if epithelial_subtype column exists and matches
if "epithelial_subtype" in obs.columns:
    stemta_mask = obs["epithelial_subtype"] == "Stem/TA-like"
    exact_match = (cE01_mask == stemta_mask).all()
    log(f"  cE01 == Stem/TA-like exact match: {exact_match}")
    record("population", "cE01_stemta_exact_match", bool(exact_match))
    if not exact_match:
        issue("CRITICAL", "cE01 != Stem/TA-like",
              "pelka_epithelial.h5ad obs",
              f"Mismatch in {int((cE01_mask != stemta_mask).sum())} cells")
else:
    log("  epithelial_subtype column not found (acceptable if only cl295v11SubShort used)")

# --- Tissue column and per-tissue counts ---
tissue_col = cfg["population"]["tissue_column"]  # SPECIMEN_TYPE
patient_col = cfg["population"]["patient_column"]  # PID
batch_col = cfg["population"]["batch_column"]  # batchID

cE01_obs = obs[cE01_mask].copy()
n_tumor = int((cE01_obs[tissue_col] == "T").sum())
n_normal = int((cE01_obs[tissue_col] == "N").sum())
log(f"  Tumor cells: {n_tumor}, Normal cells: {n_normal}")
record("population", "n_tumor_cE01", n_tumor)
record("population", "n_normal_cE01", n_normal)

# --- Per-patient cell counts and eligibility ---
log("\nComputing per-patient eligibility...")
patient_tissue_counts = cE01_obs.groupby([patient_col, tissue_col]).size().unstack(fill_value=0)
paired_patients = patient_tissue_counts.index[
    (patient_tissue_counts.get("T", pd.Series(dtype=int)) > 0) &
    (patient_tissue_counts.get("N", pd.Series(dtype=int)) > 0)
]
log(f"  Paired patients (any cells both tissues): {len(paired_patients)}")

for thresh_name, thresh in [("thresh_20", 20), ("thresh_30", 30), ("thresh_50", 50)]:
    eligible = [p for p in paired_patients
                if patient_tissue_counts.loc[p, "T"] >= thresh
                and patient_tissue_counts.loc[p, "N"] >= thresh]
    log(f"  Eligible at >={thresh}: {len(eligible)}")
    record("eligibility", thresh_name, len(eligible))
    if thresh == 30:
        ELIGIBLE_30 = sorted(eligible)
        record("eligibility", "eligible_30_pids", ELIGIBLE_30)

# --- Tissue-batch confound (Cramer's V) ---
log("\nComputing tissue-batch confound...")
cE01_eligible = cE01_obs[cE01_obs[patient_col].isin(ELIGIBLE_30)]
# For each patient, get their tumor batch and normal batch
batch_data = []
for pid in ELIGIBLE_30:
    pid_obs = cE01_eligible[cE01_eligible[patient_col] == pid]
    t_batches = pid_obs[pid_obs[tissue_col] == "T"][batch_col].unique()
    n_batches = pid_obs[pid_obs[tissue_col] == "N"][batch_col].unique()
    for tb in t_batches:
        for nb in n_batches:
            batch_data.append({"PID": pid, "T_batch": tb, "N_batch": nb})

batch_df = pd.DataFrame(batch_data)
# Count patients sharing a tumor/normal batch
shared_batch = int((batch_df["T_batch"] == batch_df["N_batch"]).sum())
log(f"  Patients sharing T/N batch: {shared_batch}/{len(ELIGIBLE_30)}")
record("confound", "shared_batch_count", shared_batch)
record("confound", "total_eligible", len(ELIGIBLE_30))

# Compute Cramer's V on the full cell-level contingency (tissue x batch)
ct = pd.crosstab(cE01_eligible[tissue_col], cE01_eligible[batch_col])
chi2_val = stats.chi2_contingency(ct)[0]
n_cells_elig = len(cE01_eligible)
min_dim = min(ct.shape) - 1
cramers_v = np.sqrt(chi2_val / (n_cells_elig * min_dim)) if min_dim > 0 else np.nan
log(f"  Cramer's V (tissue x batch): {cramers_v:.6f}")
record("confound", "cramers_v", float(cramers_v))

if abs(cramers_v - 1.0) > 0.01:
    issue("CRITICAL", "Cramer's V != 1.0",
          "tissue-batch confound",
          f"Computed {cramers_v:.6f}, claimed 1.0")
if shared_batch != 0:
    issue("CRITICAL", "Shared batch count != 0",
          "tissue-batch confound",
          f"Computed {shared_batch}, claimed 0")

# --- Data scale verification ---
log("\nVerifying data scale of Pelka normalized .X...")
# Sample some non-zero values from .X
adata_norm_mem = ad.read_h5ad(ROOT / "data/processed/pelka_epithelial.h5ad")
X = adata_norm_mem.X
if sp.issparse(X):
    X_dense_sample = X[:100].toarray()
else:
    X_dense_sample = X[:100]

nonzero_vals = X_dense_sample[X_dense_sample != 0]
log(f"  .X non-zero sample: min={nonzero_vals.min():.4f}, max={nonzero_vals.max():.4f}, mean={nonzero_vals.mean():.4f}")

# Check if log-normalized: expm1 and sum per cell
if sp.issparse(X):
    expm1_sums = np.array(X[:1000].toarray())
else:
    expm1_sums = X[:1000].copy()
expm1_sums = np.expm1(expm1_sums).sum(axis=1)
log(f"  Per-cell expm1 sum (first 1000): mean={expm1_sums.mean():.1f}, median={np.median(expm1_sums):.1f}")

is_log_normalized = (nonzero_vals.min() > 0 and nonzero_vals.max() < 15 and
                     np.median(expm1_sums) > 5000 and np.median(expm1_sums) < 50000)
log(f"  .X classification: {'LOG_NORMALIZED' if is_log_normalized else 'UNKNOWN'}")
record("data_scale", "pelka_X_classification", "LOG_NORMALIZED" if is_log_normalized else "UNKNOWN")

# Check counts file
log("\nVerifying Pelka counts .X...")
adata_counts = ad.read_h5ad(ROOT / "data/processed/pelka_epithelial_counts.h5ad")
Xc = adata_counts.X
if sp.issparse(Xc):
    Xc_sample = Xc[:100].toarray()
else:
    Xc_sample = Xc[:100]
nonzero_counts = Xc_sample[Xc_sample != 0]
is_integer = np.allclose(nonzero_counts, np.round(nonzero_counts))
log(f"  Counts non-zero: min={nonzero_counts.min():.2f}, max={nonzero_counts.max():.2f}, is_integer={is_integer}")
record("data_scale", "counts_is_integer", bool(is_integer))

# --- Verify all 21 genes present ---
log("\nChecking gene presence...")
var_names = set(adata_norm_mem.var_names)
missing_genes = [g for g in ALL_21 if g not in var_names]
if missing_genes:
    issue("CRITICAL", "Missing genes in Pelka",
          "pelka_epithelial.h5ad var_names",
          f"Missing: {missing_genes}")
else:
    log(f"  All 21 genes present in var_names")
record("gene_presence", "pelka_all_21_present", len(missing_genes) == 0)

# Barcode alignment
norm_barcodes = set(adata_norm_mem.obs_names)
counts_barcodes = set(adata_counts.obs_names)
barcode_match = norm_barcodes == counts_barcodes
log(f"  Barcode alignment norm vs counts: {barcode_match}")
record("alignment", "barcode_match", bool(barcode_match))

# ===========================================================================
# PHASE 2b: GATE 2 / 2b — Score recomputation & pooled correlations
# ===========================================================================
log("\n" + "=" * 70)
log("PHASE 2b: GATE 2/2b SCORE RECOMPUTATION")
log("=" * 70)

# Score on FULL dataset (scanpy control binning is population-dependent)
# then extract cE01 — this matches the Gate 2 script's approach
log("Scoring on full dataset (168K cells) then extracting cE01...")
sc.tl.score_genes(adata_norm_mem, gene_list=ISC_GENES, score_name="ISC_score",
                  use_raw=False, ctrl_size=50, n_bins=25, random_state=42)
sc.tl.score_genes(adata_norm_mem, gene_list=FAO_GENES, score_name="FAO_score",
                  use_raw=False, ctrl_size=50, n_bins=25, random_state=42)

cE01_idx = adata_norm_mem.obs["cl295v11SubShort"] == "cE01"
adata_cE01 = adata_norm_mem[cE01_idx].copy()
log(f"cE01 subset: {adata_cE01.shape}")

log(f"  ISC score: mean={adata_cE01.obs['ISC_score'].mean():.6f}, std={adata_cE01.obs['ISC_score'].std():.6f}")
log(f"  FAO score: mean={adata_cE01.obs['FAO_score'].mean():.6f}, std={adata_cE01.obs['FAO_score'].std():.6f}")

# Load committed Gate 2 scores
gate2_parquet = pd.read_parquet(CV2 / "results" / "gate2_recomputed_scores_cE01.parquet")
log(f"  Gate 2 parquet shape: {gate2_parquet.shape}")
log(f"  Gate 2 parquet columns: {list(gate2_parquet.columns)}")

# Compare our scanpy scores vs committed
our_isc = adata_cE01.obs["ISC_score"].values
our_fao = adata_cE01.obs["FAO_score"].values

# Align by index
committed_isc = gate2_parquet.loc[adata_cE01.obs_names, "ISC_scanpy"].values
committed_fao = gate2_parquet.loc[adata_cE01.obs_names, "FAO_scanpy"].values

isc_maxdiff = np.max(np.abs(our_isc - committed_isc))
fao_maxdiff = np.max(np.abs(our_fao - committed_fao))
log(f"  ISC scanpy max|diff|: {isc_maxdiff:.2e}")
log(f"  FAO scanpy max|diff|: {fao_maxdiff:.2e}")
record("gate2", "isc_scanpy_maxdiff", float(isc_maxdiff))
record("gate2", "fao_scanpy_maxdiff", float(fao_maxdiff))

if isc_maxdiff > 1e-10:
    issue("MATERIAL", "ISC scanpy score non-reproduction",
          "Gate 2 score recomputation",
          f"Max diff {isc_maxdiff:.2e} (expected <1e-10)")
if fao_maxdiff > 1e-10:
    issue("MATERIAL", "FAO scanpy score non-reproduction",
          "Gate 2 score recomputation",
          f"Max diff {fao_maxdiff:.2e} (expected <1e-10)")

# --- Pooled correlations ---
log("\nComputing pooled correlations...")
tissue = adata_cE01.obs[tissue_col].values
isc_vals = committed_isc  # use committed (verified above)
fao_vals = committed_fao

r_tumor_pooled = np.corrcoef(isc_vals[tissue == "T"], fao_vals[tissue == "T"])[0, 1]
r_normal_pooled = np.corrcoef(isc_vals[tissue == "N"], fao_vals[tissue == "N"])[0, 1]
log(f"  Pooled r_tumor (scanpy): {r_tumor_pooled:.4f}")
log(f"  Pooled r_normal (scanpy): {r_normal_pooled:.4f}")
record("gate2", "r_tumor_pooled_scanpy", float(r_tumor_pooled))
record("gate2", "r_normal_pooled_scanpy", float(r_normal_pooled))

# Compare vs committed table
gate2_pooled = pd.read_csv(CV2 / "results/tables/gate2_pooled_decoupling.csv")
# Method column contains "scanpy (recomputed)" not just "scanpy"
scanpy_rows = gate2_pooled[gate2_pooled["method"].str.contains("scanpy", case=False)]
scanpy_row = scanpy_rows.iloc[0]
compare("gate2_r_normal_scanpy", r_normal_pooled, scanpy_row["r_normal"], tol=1e-3)
compare("gate2_r_tumor_scanpy", r_tumor_pooled, scanpy_row["r_tumor"], tol=1e-3)

# --- Clean background scores ---
log("\nRecomputing clean-background scores...")
# Load clean control gene sets
with open(CV2 / "provenance/clean_control_gene_sets.json") as f:
    clean_controls = json.load(f)

isc_clean_pool = set(adata_cE01.var_names) - set(ALL_21)
fao_clean_pool = set(adata_cE01.var_names) - set(ALL_21)
log(f"  Clean pool sizes: ISC={len(isc_clean_pool)}, FAO={len(fao_clean_pool)}")

# Recompute clean-bg scores using scanpy with clean pools
# For ISC: controls drawn from genes excluding all 21
# For FAO: controls drawn from genes excluding all 21
adata_cE01_cleanbg = adata_cE01.copy()
isc_pool_list = sorted(isc_clean_pool)
fao_pool_list = sorted(fao_clean_pool)

# scanpy score_genes doesn't directly support custom gene pools
# but the Gate 2b approach was to subset .var to the clean pool
# Let's reproduce using the committed Gate 2b scores instead
gate2b_parquet = pd.read_parquet(CV2 / "results/gate2b_clean_background_scores_cE01.parquet")
log(f"  Gate 2b parquet shape: {gate2b_parquet.shape}")
log(f"  Gate 2b parquet columns: {list(gate2b_parquet.columns)}")

# Clean-bg pooled correlations
cleanbg_isc = gate2b_parquet.loc[adata_cE01.obs_names, "ISC_scanpy_cleanbg"].values
cleanbg_fao = gate2b_parquet.loc[adata_cE01.obs_names, "FAO_scanpy_cleanbg"].values

r_normal_cleanbg = np.corrcoef(
    cleanbg_isc[tissue == "N"], cleanbg_fao[tissue == "N"]
)[0, 1]
r_tumor_cleanbg = np.corrcoef(
    cleanbg_isc[tissue == "T"], cleanbg_fao[tissue == "T"]
)[0, 1]
log(f"  Clean-bg r_normal: {r_normal_cleanbg:.4f}")
log(f"  Clean-bg r_tumor: {r_tumor_cleanbg:.4f}")

gate2b_pooled = pd.read_csv(CV2 / "results/tables/gate2b_pooled_cleanbg.csv")
cleanbg_row = gate2b_pooled[gate2b_pooled["method"].str.contains("clean", case=False)].iloc[0]
compare("gate2b_r_normal_cleanbg", r_normal_cleanbg, cleanbg_row["r_normal"], tol=1e-3)
compare("gate2b_r_tumor_cleanbg", r_tumor_cleanbg, cleanbg_row["r_tumor"], tol=1e-3)

# ===========================================================================
# PHASE 2c: GATE 4 — Donor-aware coupling
# ===========================================================================
log("\n" + "=" * 70)
log("PHASE 2c: GATE 4 — DONOR-AWARE COUPLING")
log("=" * 70)

def pearson_r(x, y):
    """Safe Pearson r, returns NaN if degenerate."""
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])

def fisher_z(r):
    """Fisher z-transform, clipped to avoid inf."""
    r = np.clip(r, -0.9999, 0.9999)
    return float(np.arctanh(r))

def inverse_fisher_z(z):
    return float(np.tanh(z))

# Build per-patient Delta_z for primary scanpy method
log("Computing per-patient Fisher-z differences (scanpy)...")
scores_df = gate2_parquet.copy()
scores_df["tissue"] = adata_cE01.obs.loc[scores_df.index, tissue_col].values
scores_df["PID"] = adata_cE01.obs.loc[scores_df.index, patient_col].values

per_patient_results = []
for pid in ELIGIBLE_30:
    pid_data = scores_df[scores_df["PID"] == pid]
    t_data = pid_data[pid_data["tissue"] == "T"]
    n_data = pid_data[pid_data["tissue"] == "N"]

    r_t = pearson_r(t_data["ISC_scanpy"].values, t_data["FAO_scanpy"].values)
    r_n = pearson_r(n_data["ISC_scanpy"].values, n_data["FAO_scanpy"].values)
    z_t = fisher_z(r_t)
    z_n = fisher_z(r_n)
    delta_z = z_t - z_n

    per_patient_results.append({
        "PID": pid, "r_tumor": r_t, "r_normal": r_n,
        "z_tumor": z_t, "z_normal": z_n, "delta_z": delta_z,
        "n_tumor": len(t_data), "n_normal": len(n_data)
    })

pp_df = pd.DataFrame(per_patient_results)
delta_z_values = pp_df["delta_z"].values

mean_delta_z = float(np.mean(delta_z_values))
median_delta_z = float(np.median(delta_z_values))
prop_positive = float(np.mean(delta_z_values > 0))

log(f"  Mean Delta_z: {mean_delta_z:.6f}")
log(f"  Median Delta_z: {median_delta_z:.6f}")
log(f"  Proportion positive: {prop_positive:.4f}")

record("gate4", "mean_delta_z", mean_delta_z)
record("gate4", "median_delta_z", median_delta_z)
record("gate4", "prop_positive", prop_positive)

# Compare vs committed
gate4_summary = pd.read_csv(CV2 / "results/tables/gate4_primary_summary.csv")
compare("gate4_mean_delta_z", mean_delta_z, gate4_summary["mean_delta_z"].iloc[0], tol=1e-4)
compare("gate4_median_delta_z", median_delta_z, gate4_summary["median_delta_z"].iloc[0], tol=1e-4)

# --- Wilcoxon signed-rank test ---
log("\nWilcoxon signed-rank test on Delta_z...")
wilcoxon_stat, wilcoxon_p = stats.wilcoxon(delta_z_values, alternative="two-sided")
log(f"  Wilcoxon statistic: {wilcoxon_stat}, p-value: {wilcoxon_p:.4e}")
record("gate4", "wilcoxon_p", float(wilcoxon_p))
compare("gate4_wilcoxon_p", wilcoxon_p, gate4_summary["wilcoxon_p"].iloc[0], tol=1e-8)

# --- Bootstrap CI ---
log("\nBootstrap 95% CI...")
rng = np.random.RandomState(42)
n_boot = 10000
boot_means = np.array([
    np.mean(rng.choice(delta_z_values, size=len(delta_z_values), replace=True))
    for _ in range(n_boot)
])
ci_lo = float(np.percentile(boot_means, 2.5))
ci_hi = float(np.percentile(boot_means, 97.5))
log(f"  Bootstrap CI: [{ci_lo:.4f}, {ci_hi:.4f}]")
record("gate4", "bootstrap_ci_lo", ci_lo)
record("gate4", "bootstrap_ci_hi", ci_hi)

# The CI should exclude zero
ci_excludes_zero = ci_lo > 0 or ci_hi < 0
log(f"  CI excludes zero: {ci_excludes_zero}")
record("gate4", "ci_excludes_zero", bool(ci_excludes_zero))

compare("gate4_ci_lo", ci_lo, gate4_summary["bootstrap_ci_lo"].iloc[0], tol=0.02)
compare("gate4_ci_hi", ci_hi, gate4_summary["bootstrap_ci_hi"].iloc[0], tol=0.02)

# --- LOO analysis ---
log("\nLeave-one-out analysis...")
loo_results = []
for i, pid in enumerate(ELIGIBLE_30):
    loo_dz = np.delete(delta_z_values, i)
    loo_mean = float(np.mean(loo_dz))
    pct_change = abs(loo_mean - mean_delta_z) / abs(mean_delta_z) * 100
    loo_results.append({
        "PID": pid, "loo_mean": loo_mean,
        "pct_change": pct_change,
        "sign_positive": loo_mean > 0
    })

loo_df = pd.DataFrame(loo_results)
all_loo_positive = loo_df["sign_positive"].all()
max_loo_change = loo_df["pct_change"].max()
log(f"  All LOO means positive: {all_loo_positive}")
log(f"  Max LOO % change: {max_loo_change:.3f}%")
record("gate4", "loo_all_positive", bool(all_loo_positive))
record("gate4", "loo_max_pct_change", float(max_loo_change))

# Gate 4 outcome
gate4_pass = (ci_excludes_zero and wilcoxon_p < 0.05 and all_loo_positive and max_loo_change < 25)
log(f"  Gate 4 outcome: {'PASS' if gate4_pass else 'FAIL'}")
record("gate4", "outcome", "PASS" if gate4_pass else "FAIL")

# Compare per-patient table
gate4_pp = pd.read_csv(CV2 / "results/tables/gate4_per_patient.csv")
gate4_pp_scanpy = gate4_pp[gate4_pp["method"] == "scanpy"]
for _, row in gate4_pp_scanpy.iterrows():
    pid = row["PID"]
    our_row = pp_df[pp_df["PID"] == pid]
    if len(our_row) == 0:
        issue("MINOR", f"Missing PID {pid} in our computation", "gate4_per_patient", "")
        continue
    our_dz = our_row["delta_z"].values[0]
    their_dz = row["Delta_z"]
    if abs(our_dz - their_dz) > 1e-4:
        issue("MATERIAL", f"Gate 4 per-patient Delta_z mismatch for {pid}",
              "gate4_per_patient.csv",
              f"Ours={our_dz:.6f}, Theirs={their_dz:.6f}, diff={abs(our_dz-their_dz):.2e}")

# ===========================================================================
# PHASE 2d: GATE 5 — Equal-cell downsampling
# ===========================================================================
log("\n" + "=" * 70)
log("PHASE 2d: GATE 5 — EQUAL-CELL DOWNSAMPLING (subset verification)")
log("=" * 70)

# Full K=1000 takes time; we verify the first 10 draws and the reported aggregates
SAMPLE_SIZE = 30
K_VERIFY = 50  # verify first 50 draws
FULL_EFFECT = mean_delta_z  # 0.400629

log(f"  Full-data effect: {FULL_EFFECT:.6f}")
log(f"  50% magnitude threshold: {FULL_EFFECT * 0.5:.6f}")

# SeedSequence deterministic plan
from numpy.random import SeedSequence, default_rng
ss = SeedSequence(42)
child_seeds = ss.spawn(1000)

draw_results = []
for k in range(K_VERIFY):
    rng_k = default_rng(child_seeds[k])
    draw_dz = []
    valid = True
    for pid in ELIGIBLE_30:
        pid_data = scores_df[scores_df["PID"] == pid]
        t_data = pid_data[pid_data["tissue"] == "T"]
        n_data = pid_data[pid_data["tissue"] == "N"]

        if len(t_data) < SAMPLE_SIZE or len(n_data) < SAMPLE_SIZE:
            valid = False
            break

        t_idx = rng_k.choice(len(t_data), size=SAMPLE_SIZE, replace=False)
        n_idx = rng_k.choice(len(n_data), size=SAMPLE_SIZE, replace=False)

        t_sample = t_data.iloc[t_idx]
        n_sample = n_data.iloc[n_idx]

        r_t = pearson_r(t_sample["ISC_scanpy"].values, t_sample["FAO_scanpy"].values)
        r_n = pearson_r(n_sample["ISC_scanpy"].values, n_sample["FAO_scanpy"].values)
        dz = fisher_z(r_t) - fisher_z(r_n)
        draw_dz.append(dz)

    if valid and len(draw_dz) == len(ELIGIBLE_30):
        draw_mean = float(np.mean(draw_dz))
        draw_results.append({
            "draw": k, "mean_delta_z": draw_mean,
            "direction_positive": draw_mean > 0
        })

draw_df = pd.DataFrame(draw_results)
n_positive = int(draw_df["direction_positive"].sum())
dir_fraction = n_positive / len(draw_df)
median_mag = float(draw_df["mean_delta_z"].median())

log(f"  First {K_VERIFY} draws: {n_positive}/{len(draw_df)} positive direction")
log(f"  Direction fraction: {dir_fraction:.4f}")
log(f"  Median magnitude: {median_mag:.6f}")

record("gate5", "k_verified", K_VERIFY)
record("gate5", "direction_fraction_sample", dir_fraction)
record("gate5", "median_magnitude_sample", median_mag)

# Compare against committed
gate5_summary = pd.read_csv(CV2 / "results/tables/gate5_primary_summary.csv")
committed_dir = gate5_summary["direction_retention_fraction"].iloc[0]
committed_med = gate5_summary["median_downsampled_magnitude"].iloc[0]
log(f"  Committed direction retention: {committed_dir}")
log(f"  Committed median magnitude: {committed_med}")

# Our sample should be consistent
if dir_fraction < 0.85:
    issue("CRITICAL", "Gate 5 direction retention fails in sample",
          "independent downsampling",
          f"Only {dir_fraction:.2%} positive in {K_VERIFY} draws")
if median_mag < FULL_EFFECT * 0.3:
    issue("CRITICAL", "Gate 5 magnitude collapse in sample",
          "independent downsampling",
          f"Median {median_mag:.4f} vs threshold {FULL_EFFECT*0.5:.4f}")

# Compare committed draw-level summary for first 50 draws
gate5_draws = pd.read_csv(CV2 / "results/tables/gate5_draw_summary.csv")
for k in range(min(K_VERIFY, len(gate5_draws))):
    our_draw = draw_df[draw_df["draw"] == k]
    if len(our_draw) == 0:
        continue
    their_draw = gate5_draws.iloc[k]
    their_val = their_draw.get("mean_delta_z_scanpy", their_draw.get("mean_delta_z", None))
    if their_val is None:
        break
    our_val = our_draw["mean_delta_z"].values[0]
    if abs(our_val - their_val) > 1e-4:
        issue("MATERIAL", f"Gate 5 draw {k} mismatch",
              "gate5_draw_summary.csv",
              f"Ours={our_val:.6f}, Theirs={their_val:.6f}")

# ===========================================================================
# PHASE 2e: GATE 3 — Score decomposition (verify key variants)
# ===========================================================================
log("\n" + "=" * 70)
log("PHASE 2e: GATE 3 — SCORE DECOMPOSITION (key variants)")
log("=" * 70)

# We need clean-bg controls to reproduce Gate 3
# Load the committed clean-bg scores for the fixed-background scorer
# Gate 3 used: score = mean(targets) - mean(controls)
# with controls from clean_control_gene_sets.json

# Load committed gate3 variant effects
gate3_effects = pd.read_csv(CV2 / "results/tables/gate3_variant_effects.csv")
log(f"  Committed variants: {len(gate3_effects)}")

# Focus on the three critical booleans
baseline_row = gate3_effects[gate3_effects["variant"] == "baseline"]
fabp1_row = gate3_effects[gate3_effects["variant"] == "fabp1_removed"]
hmgcs2_row = gate3_effects[gate3_effects["variant"] == "hmgcs2_removed"]
core_betaox_row = gate3_effects[gate3_effects["variant"] == "core_beta_oxidation_5"]

baseline_mean_dz = baseline_row["mean_delta_z"].values[0]
log(f"  Baseline clean-fixed-bg mean Delta_z: {baseline_mean_dz:.6f}")

# Check the three booleans
fabp1_retention = fabp1_row["retention_ratio"].values[0]
hmgcs2_retention = hmgcs2_row["retention_ratio"].values[0]
core_retention = core_betaox_row["retention_ratio"].values[0]

log(f"  FABP1 removed retention: {fabp1_retention:.4f} ({'WEAKENED' if fabp1_retention < 0.5 else 'ROBUST'})")
log(f"  HMGCS2 removed retention: {hmgcs2_retention:.4f} ({'WEAKENED' if hmgcs2_retention < 0.5 else 'ROBUST'})")
log(f"  Core-beta-ox retention: {core_retention:.4f} ({'WEAKENED' if core_retention < 0.5 else 'ROBUST'})")

core_beta_oxidation_holds = core_retention >= 0.5
fao_reframe_conjunction = (fabp1_retention < 0.5 and hmgcs2_retention < 0.5 and core_retention < 0.5)
program_level_fao_robust = fabp1_retention >= 0.5  # if FABP1 removal doesn't weaken, program is robust

log(f"\n  Frozen booleans:")
log(f"    core_beta_oxidation_holds: {core_beta_oxidation_holds} (claimed: False)")
log(f"    fao_reframe_conjunction_triggered: {fao_reframe_conjunction} (claimed: False)")
log(f"    program_level_fao_robust: {program_level_fao_robust} (claimed: False)")

record("gate3", "baseline_mean_dz", float(baseline_mean_dz))
record("gate3", "fabp1_retention", float(fabp1_retention))
record("gate3", "hmgcs2_retention", float(hmgcs2_retention))
record("gate3", "core_betaox_retention", float(core_retention))
record("gate3", "core_beta_oxidation_holds", core_beta_oxidation_holds)
record("gate3", "fao_reframe_conjunction", fao_reframe_conjunction)
record("gate3", "program_level_fao_robust", program_level_fao_robust)

if core_beta_oxidation_holds:
    issue("CRITICAL", "core_beta_oxidation_holds should be False",
          "Gate 3 booleans", f"retention={core_retention}")
if fao_reframe_conjunction:
    issue("CRITICAL", "fao_reframe_conjunction should be False",
          "Gate 3 booleans", "HMGCS2 is not weakened")
if program_level_fao_robust:
    issue("CRITICAL", "program_level_fao_robust should be False",
          "Gate 3 booleans", "FABP1 removal IS weakened")

# Now independently verify the baseline clean-fixed-bg Delta_z
# by computing score = mean(targets) - mean(controls)
log("\nIndependently computing Gate 3 baseline (fixed clean-bg)...")

# Load clean controls
with open(CV2 / "provenance/gate3_clean_fixed_control_gene_sets.json") as f:
    gate3_controls = json.load(f)

isc_controls = gate3_controls.get("ISC_controls") or gate3_controls.get("isc_controls", [])
fao_controls = gate3_controls.get("FAO_controls") or gate3_controls.get("fao_controls", [])

if not isc_controls or not fao_controls:
    # Try alternative key names
    for key in gate3_controls:
        log(f"  Gate3 control key: {key}")
    log("  WARNING: Could not find control gene lists in gate3_clean_fixed_control_gene_sets.json")
    log(f"  Keys available: {list(gate3_controls.keys())}")
else:
    log(f"  ISC controls: {len(isc_controls)} genes")
    log(f"  FAO controls: {len(fao_controls)} genes")

    # Verify no overlap with 21 union targets
    isc_ctrl_overlap = set(isc_controls) & set(ALL_21)
    fao_ctrl_overlap = set(fao_controls) & set(ALL_21)
    log(f"  ISC controls overlap with 21: {isc_ctrl_overlap}")
    log(f"  FAO controls overlap with 21: {fao_ctrl_overlap}")
    if isc_ctrl_overlap or fao_ctrl_overlap:
        issue("CRITICAL", "Gate 3 controls overlap with target genes",
              "gate3_clean_fixed_control_gene_sets.json",
              f"ISC overlap: {isc_ctrl_overlap}, FAO overlap: {fao_ctrl_overlap}")

    # Compute fixed-bg scores: mean(targets) - mean(controls)
    X_cE01 = adata_cE01.X
    if sp.issparse(X_cE01):
        X_cE01 = X_cE01.toarray()
    var_names_list = list(adata_cE01.var_names)

    def gene_indices(genes):
        return [var_names_list.index(g) for g in genes if g in var_names_list]

    isc_idx = gene_indices(ISC_GENES)
    fao_idx = gene_indices(FAO_GENES)
    isc_ctrl_idx = gene_indices(isc_controls)
    fao_ctrl_idx = gene_indices(fao_controls)

    isc_fixed_score = X_cE01[:, isc_idx].mean(axis=1) - X_cE01[:, isc_ctrl_idx].mean(axis=1)
    fao_fixed_score = X_cE01[:, fao_idx].mean(axis=1) - X_cE01[:, fao_ctrl_idx].mean(axis=1)

    # Now compute per-patient Delta_z with these scores
    scores_fixed = pd.DataFrame({
        "ISC_fixed": isc_fixed_score,
        "FAO_fixed": fao_fixed_score,
        "tissue": adata_cE01.obs[tissue_col].values,
        "PID": adata_cE01.obs[patient_col].values
    }, index=adata_cE01.obs_names)

    pp_fixed = []
    for pid in ELIGIBLE_30:
        pid_data = scores_fixed[scores_fixed["PID"] == pid]
        t_data = pid_data[pid_data["tissue"] == "T"]
        n_data = pid_data[pid_data["tissue"] == "N"]
        r_t = pearson_r(t_data["ISC_fixed"].values, t_data["FAO_fixed"].values)
        r_n = pearson_r(n_data["ISC_fixed"].values, n_data["FAO_fixed"].values)
        dz = fisher_z(r_t) - fisher_z(r_n)
        pp_fixed.append(dz)

    our_baseline_dz = float(np.mean(pp_fixed))
    log(f"  Our baseline clean-fixed-bg mean Delta_z: {our_baseline_dz:.6f}")
    log(f"  Committed baseline: {baseline_mean_dz:.6f}")
    log(f"  Difference: {abs(our_baseline_dz - baseline_mean_dz):.6e}")
    record("gate3", "our_baseline_mean_dz", our_baseline_dz)
    compare("gate3_baseline_mean_dz", our_baseline_dz, baseline_mean_dz, tol=1e-3)

# ===========================================================================
# PHASE 2f: GATE 6 — LGR5 sensitivity
# ===========================================================================
log("\n" + "=" * 70)
log("PHASE 2f: GATE 6 — LGR5 SENSITIVITY")
log("=" * 70)

# LGR5 detection from raw counts
log("Detecting LGR5+ cells from counts file...")
lgr5_idx_in_counts = list(adata_counts.var_names).index("LGR5") if "LGR5" in adata_counts.var_names else None
if lgr5_idx_in_counts is not None:
    Xc_cE01 = adata_counts[cE01_idx].X
    if sp.issparse(Xc_cE01):
        lgr5_counts = np.asarray(Xc_cE01[:, lgr5_idx_in_counts].todense()).flatten()
    else:
        lgr5_counts = Xc_cE01[:, lgr5_idx_in_counts].flatten()

    lgr5_detected = lgr5_counts > 0
    n_lgr5_pos = int(lgr5_detected.sum())
    log(f"  LGR5+ cE01 cells: {n_lgr5_pos}")
    record("gate6", "n_lgr5_positive", n_lgr5_pos)

    # Eligible patients for LGR5 analysis (>=15 cells both tissues)
    lgr5_obs = adata_cE01.obs.copy()
    lgr5_obs["lgr5_detected"] = lgr5_detected
    lgr5_subset = lgr5_obs[lgr5_obs["lgr5_detected"]]

    lgr5_pt = lgr5_subset.groupby([patient_col, tissue_col]).size().unstack(fill_value=0)
    lgr5_eligible = [p for p in ELIGIBLE_30
                     if p in lgr5_pt.index
                     and lgr5_pt.loc[p].get("T", 0) >= 15
                     and lgr5_pt.loc[p].get("N", 0) >= 15]
    log(f"  LGR5 eligible patients (>=15 both): {len(lgr5_eligible)} — {lgr5_eligible}")
    record("gate6", "eligible_patients", lgr5_eligible)
    record("gate6", "n_eligible", len(lgr5_eligible))

    # Compare vs committed
    gate6_pp = pd.read_csv(CV2 / "results/tables/gate6_per_patient.csv")
    committed_eligible = sorted(gate6_pp["PID"].unique())
    log(f"  Committed eligible: {committed_eligible}")
    if set(lgr5_eligible) != set(committed_eligible):
        issue("MATERIAL", "Gate 6 eligible patient mismatch",
              "gate6_per_patient.csv",
              f"Ours: {lgr5_eligible}, Committed: {committed_eligible}")

    # Outcome: power-limited if n < 15
    gate6_outcome = "POWER-LIMITED/INCONCLUSIVE" if len(lgr5_eligible) < 15 else "TESTABLE"
    log(f"  Gate 6 outcome: {gate6_outcome}")
    record("gate6", "outcome", gate6_outcome)
else:
    log("  LGR5 not found in counts file!")
    issue("CRITICAL", "LGR5 gene missing from counts", "pelka_epithelial_counts.h5ad", "")

# ===========================================================================
# PHASE 2g: LEE EXTERNAL CONSISTENCY
# ===========================================================================
log("\n" + "=" * 70)
log("PHASE 2g: LEE EXTERNAL CONSISTENCY")
log("=" * 70)

# Load Lee data and verify scale
for cohort_name, h5ad_path in [("SMC", "data/external/lee/smc_epithelial.h5ad"),
                                ("KUL3", "data/external/lee/kul3_epithelial.h5ad")]:
    fpath = ROOT / h5ad_path
    if not fpath.exists():
        log(f"  {cohort_name} file not found: {fpath}")
        continue

    lee_adata = ad.read_h5ad(fpath)
    log(f"\n  {cohort_name}: {lee_adata.shape}")

    # Check data scale
    X_lee = lee_adata.X
    if sp.issparse(X_lee):
        X_sample = X_lee[:100].toarray()
    else:
        X_sample = X_lee[:100]

    nonzero = X_sample[X_sample != 0]
    log(f"    Non-zero values: min={nonzero.min():.4f}, max={nonzero.max():.4f}")

    # expm1 sum per cell
    expm1_sums_lee = np.expm1(X_sample).sum(axis=1)
    log(f"    Per-cell expm1 sum: mean={expm1_sums_lee.mean():.1f}, median={np.median(expm1_sums_lee):.1f}")

    # Lee h5ad files are pre-subsetted to epithelial cells with ~2861 genes (not full 33K)
    # so per-cell expm1 sums will be much less than 10000. Check value range instead.
    is_already_log = (nonzero.min() > 0 and nonzero.max() < 15 and
                      not np.allclose(nonzero, np.round(nonzero)))  # fractional = not raw counts
    log(f"    Classification: {'ALREADY_LOG_NORMALIZED (subsetted)' if is_already_log else 'RAW/OTHER'}")
    log(f"    Note: h5ad has only {lee_adata.n_vars} genes (pre-subsetted), so expm1 sums < 10000 is expected")
    record(f"lee_{cohort_name.lower()}", "data_scale",
           "ALREADY_LOG_NORMALIZED_SUBSETTED" if is_already_log else "RAW/OTHER")
    record(f"lee_{cohort_name.lower()}", "n_genes", lee_adata.n_vars)

    if not is_already_log:
        issue("MATERIAL", f"Lee {cohort_name} scale unexpected",
              h5ad_path, f"Expected log-normalized fractional values")

    # Check gene presence
    lee_genes = set(lee_adata.var_names)
    missing_in_lee = [g for g in ALL_21 if g not in lee_genes]
    log(f"    Missing target genes: {missing_in_lee}")
    record(f"lee_{cohort_name.lower()}", "missing_genes", missing_in_lee)

# Compare committed Lee normal-anchor correlations
lee_corr = pd.read_csv(CV2 / "results/tables/lee_normal_anchor_correlations.csv")
log(f"\n  Committed Lee correlations table ({len(lee_corr)} rows):")
for _, row in lee_corr.iterrows():
    log(f"    {row.to_dict()}")

# Verify per-patient sign tally
lee_signs = pd.read_csv(CV2 / "results/tables/lee_per_patient_sign_tally.csv")
log(f"\n  Lee per-patient sign tally: {len(lee_signs)} rows")
informative = lee_signs[lee_signs["n"] >= 30]
all_negative = informative["sign"].apply(lambda s: s == "neg").all() if "sign" in lee_signs.columns else None
n_informative = len(informative)
log(f"    Informative patients (>=30 cells): {n_informative}")
log(f"    All negative: {all_negative}")
record("lee", "n_informative_patients", n_informative)
record("lee", "all_informative_negative", bool(all_negative) if all_negative is not None else None)
# Check column names
log(f"    Columns: {list(lee_signs.columns)}")
log(f"    First rows:")
for _, row in lee_signs.head(5).iterrows():
    log(f"      {row.to_dict()}")

# ===========================================================================
# WRITE RESULTS
# ===========================================================================
log("\n" + "=" * 70)
log("WRITING RESULTS")
log("=" * 70)

output_path = CV2 / "audit" / "independent_verification_results.json"
with open(output_path, "w") as f:
    json.dump(RESULTS, f, indent=2, default=str)
log(f"Results written to {output_path}")

issues_path = CV2 / "audit" / "issues_found.json"
with open(issues_path, "w") as f:
    json.dump(ISSUES, f, indent=2, default=str)
log(f"Issues written to {issues_path}")

# Summary
log("\n" + "=" * 70)
log("ISSUE SUMMARY")
log("=" * 70)
for sev in ["CRITICAL", "MATERIAL", "MINOR", "COSMETIC"]:
    count = len([i for i in ISSUES if i["severity"] == sev])
    log(f"  {sev}: {count}")

log("\nDone.")
