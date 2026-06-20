"""
Gate 6 -- LGR5 Sensitivity Without Circularity
================================================
Protocol: protocol_frozen_v9 (A008)
Execution: after Gate 4 PASS + Gate 5 PASS + Gate 3

Gate 6 is NOT a core hard kill. INCONCLUSIVE or non-supported does not
re-open A9 (Gate 4 AND Gate 5).

Tests whether the donor-aware ISC-FAO coupling shift persists in the
LGR5-transcript-detected cE01 subset with LGR5 removed from the ISC
score (8-gene). Two circularity controls: subset defined by raw count,
not score; LGR5 removed from ISC.
"""

import hashlib
import json
import warnings
from pathlib import Path
from datetime import date

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import yaml

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CV2 = Path(__file__).resolve().parent
RESULTS_DIR = CV2 / "results"
TABLES_DIR = RESULTS_DIR / "tables"
PROV_DIR = CV2 / "provenance"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
N_BOOT = 10000

# ═══════════════════════════════════════════════════════════════════════
# Copied helpers (source scripts not import-safe)
# From 04_gate4_donor_aware.py + 02_gate2_recompute_scores.py
# ═══════════════════════════════════════════════════════════════════════

def get_versions():
    import importlib.metadata as im
    return {k: im.version(k) for k in ["numpy", "scipy", "pandas", "scanpy"]}

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def sha256_str(s):
    return hashlib.sha256(s.encode()).hexdigest()

def is_near_zero_var(arr, atol=1e-15):
    return np.isclose(np.var(arr, ddof=0), 0.0, rtol=0.0, atol=atol)

def safe_pearsonr(x, y):
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < 4:
        return np.nan, "uninformative", f"n={len(x)}<4"
    if is_near_zero_var(x):
        return np.nan, "uninformative", "ISC near-zero variance"
    if is_near_zero_var(y):
        return np.nan, "uninformative", "FAO near-zero variance"
    r, _ = stats.pearsonr(x, y)
    if not np.isfinite(r):
        return np.nan, "uninformative", "non-finite r"
    return r, "ok", ""

def clip_r(r):
    if abs(r) >= 1.0:
        return np.sign(r) * 0.999999, True
    return r, False

def fisher_z(r):
    return np.arctanh(r)

def fixed_bg_score(X_sparse, var_names, gene_set, control_genes):
    var_idx = {g: i for i, g in enumerate(var_names)}
    target_idx = [var_idx[g] for g in gene_set if g in var_idx]
    ctrl_idx = [var_idx[g] for g in control_genes if g in var_idx]
    if len(target_idx) == 0:
        return np.zeros(X_sparse.shape[0])
    target_mat = X_sparse[:, target_idx]
    if sp.issparse(target_mat):
        target_mean = np.asarray(target_mat.mean(axis=1)).ravel()
    else:
        target_mean = np.asarray(target_mat).mean(axis=1)
    if len(ctrl_idx) == 0:
        return target_mean
    ctrl_mat = X_sparse[:, ctrl_idx]
    if sp.issparse(ctrl_mat):
        ctrl_mean = np.asarray(ctrl_mat.mean(axis=1)).ravel()
    else:
        ctrl_mean = np.asarray(ctrl_mat).mean(axis=1)
    return target_mean - ctrl_mean

def compute_per_patient_corr(score_df, isc_col, fao_col, pid_col, tissue_col,
                             tumor_code, normal_code):
    """Per-patient r, z, Delta_z. Returns list of dicts."""
    rows = []
    for pid in sorted(score_df[pid_col].unique()):
        pdata = score_df[score_df[pid_col] == pid]
        rec = {"PID": pid}
        for tissue_name, code, prefix in [("normal", normal_code, "normal"),
                                          ("tumor", tumor_code, "tumor")]:
            tdata = pdata[pdata[tissue_col] == code]
            n = len(tdata)
            rec[f"n_{tissue_name}"] = n
            if n == 0:
                rec[f"{prefix}_status"] = "uninformative"
                rec[f"{prefix}_reason"] = "no cells"
                rec[f"r_{prefix}"] = np.nan
                continue
            isc = tdata[isc_col].values.astype(np.float64)
            fao = tdata[fao_col].values.astype(np.float64)
            r, status, reason = safe_pearsonr(isc, fao)
            rec[f"{prefix}_status"] = status
            rec[f"{prefix}_reason"] = reason
            rec[f"r_{prefix}"] = r
        if rec.get("normal_status") == "ok" and rec.get("tumor_status") == "ok":
            rec["informative"] = True
            rn, _ = clip_r(rec["r_normal"])
            rt, _ = clip_r(rec["r_tumor"])
            rec["z_normal"] = fisher_z(rn)
            rec["z_tumor"] = fisher_z(rt)
            rec["delta_z"] = rec["z_tumor"] - rec["z_normal"]
        else:
            rec["informative"] = False
            rec["z_normal"] = np.nan
            rec["z_tumor"] = np.nan
            rec["delta_z"] = np.nan
            reasons = []
            if rec.get("normal_status") == "uninformative":
                reasons.append(f"normal:{rec.get('normal_reason','')}")
            if rec.get("tumor_status") == "uninformative":
                reasons.append(f"tumor:{rec.get('tumor_reason','')}")
            rec["exclusion_reason"] = "; ".join(reasons)
        rows.append(rec)
    return rows

def donor_aware_stats(dz_arr, seed=SEED, n_boot=N_BOOT):
    """Wilcoxon, bootstrap CI, LOO on Δz array."""
    dz = np.array(dz_arr, dtype=np.float64)
    n = len(dz)
    mean_dz = float(np.mean(dz))
    median_dz = float(np.median(dz))
    try:
        wstat, wp = stats.wilcoxon(dz, alternative="two-sided",
                                   zero_method="wilcox", method="auto")
    except ValueError:
        wstat, wp = np.nan, np.nan
    rng = np.random.RandomState(seed)
    boot_means = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        boot_means[i] = np.mean(dz[idx])
    ci_lo = float(np.percentile(boot_means, 2.5))
    ci_hi = float(np.percentile(boot_means, 97.5))
    ci_excludes_0 = (ci_lo > 0) or (ci_hi < 0)
    loo_signs = []
    loo_max_change = 0.0
    for i in range(n):
        loo_mean = float(np.mean(np.delete(dz, i)))
        loo_signs.append(loo_mean > 0)
        if abs(mean_dz) > 0:
            rel = abs(loo_mean - mean_dz) / abs(mean_dz) * 100
            loo_max_change = max(loo_max_change, rel)
    loo_all_positive = all(loo_signs)
    return {
        "n": n, "mean_delta_z": mean_dz, "median_delta_z": median_dz,
        "wilcoxon_p": float(wp) if np.isfinite(wp) else None,
        "bootstrap_ci_lo": ci_lo, "bootstrap_ci_hi": ci_hi,
        "ci_excludes_0": ci_excludes_0,
        "loo_all_positive": loo_all_positive,
        "loo_max_rel_change_pct": loo_max_change,
        "proportion_delta_z_positive": float(np.mean(dz > 0)),
    }


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("Gate 6 -- LGR5 Sensitivity Without Circularity")
    print("=" * 70)
    versions = get_versions()
    print(f"Versions: {versions}")

    with open(CV2 / "config.yaml") as f:
        cfg = yaml.safe_load(f)
    isc_9 = cfg["gene_sets"]["isc_9"]
    fao_12 = cfg["gene_sets"]["ppar_fao_12"]
    isc_8 = [g for g in isc_9 if g != "LGR5"]
    assert len(isc_8) == 8
    targets21 = set(isc_9) | set(fao_12)
    pid_col = cfg["population"]["patient_column"]
    tissue_col = cfg["population"]["tissue_column"]
    tumor_code = cfg["population"]["tumor_code"]
    normal_code = cfg["population"]["normal_code"]

    # Load controls
    clean_ctrl_path = PROV_DIR / "clean_control_gene_sets.json"
    clean_ctrl_hash = sha256_file(clean_ctrl_path)
    with open(clean_ctrl_path) as f:
        ccd = json.load(f)
    isc_controls = ccd["recovered_selected_control_genes"]["isc_controls"]
    fao_controls = ccd["recovered_selected_control_genes"]["fao_controls"]

    # Load Gate 4 PIDs
    with open(PROV_DIR / "04_gate4_donor_aware.json") as f:
        g4p = json.load(f)
    primary_34 = set(g4p["eligible_pids"])

    # Load Gate 3 reference
    g3_ve = pd.read_csv(TABLES_DIR / "gate3_variant_effects.csv")
    g3_lgr5 = g3_ve[g3_ve["variant"] == "isc_loo_LGR5"].iloc[0]
    g3_ref_mean_dz = g3_lgr5["mean_delta_z"]
    g3_ref_retention = g3_lgr5["retention_ratio"]
    g3_pp = pd.read_csv(TABLES_DIR / "gate3_per_patient_variant.csv")
    g3_pp_lgr5 = g3_pp[g3_pp["variant"] == "isc_loo_LGR5"].set_index("PID")

    # ══════════════════════════════════════════════════════════════════
    # GUARD 1: Controls
    # ══════════════════════════════════════════════════════════════════
    print("\n--- Guard 1: Controls ---")
    assert len(isc_controls) == 300, f"ISC controls: {len(isc_controls)}"
    assert len(fao_controls) == 400, f"FAO controls: {len(fao_controls)}"
    assert len(set(isc_controls)) == 300, "ISC duplicates"
    assert len(set(fao_controls)) == 400, "FAO duplicates"
    assert set(isc_controls) & targets21 == set(), f"ISC overlap: {set(isc_controls) & targets21}"
    assert set(fao_controls) & targets21 == set(), f"FAO overlap: {set(fao_controls) & targets21}"
    assert "LGR5" not in set(isc_controls), "LGR5 in ISC controls"
    assert "LGR5" not in set(fao_controls), "LGR5 in FAO controls"
    isc_ctrl_hash = sha256_str(",".join(sorted(isc_controls)))
    fao_ctrl_hash = sha256_str(",".join(sorted(fao_controls)))
    print(f"  ISC=300, FAO=400, no dups, no target overlap, LGR5 absent. PASS")

    # ══════════════════════════════════════════════════════════════════
    # GUARD 2: Alignment
    # ══════════════════════════════════════════════════════════════════
    print("\n--- Guard 2: Alignment ---")
    h5ad_path = ROOT / cfg["inputs"]["epithelial_normalized"]
    counts_path = ROOT / cfg["inputs"]["epithelial_counts"]
    adata_norm = ad.read_h5ad(h5ad_path, backed="r")
    adata_counts = ad.read_h5ad(counts_path, backed="r")

    var_names = list(adata_norm.var_names)
    # All controls in var_names
    for g in isc_controls:
        assert g in set(var_names), f"ISC control {g} missing"
    for g in fao_controls:
        assert g in set(var_names), f"FAO control {g} missing"

    assert list(adata_norm.obs_names) == list(adata_counts.obs_names), "Barcode order mismatch"
    norm_cE01 = adata_norm.obs[cfg["population"]["subtype_column_short"]] == "cE01"
    counts_cE01 = adata_counts.obs[cfg["population"]["subtype_column_short"]] == "cE01"
    assert norm_cE01.sum() == 61953
    assert counts_cE01.sum() == 61953
    assert (norm_cE01.values == counts_cE01.values).all(), "cE01 mask differs"
    # PID/SPECIMEN_TYPE
    assert (adata_norm.obs["PID"].values == adata_counts.obs["PID"].values).all()
    assert (adata_norm.obs["SPECIMEN_TYPE"].values == adata_counts.obs["SPECIMEN_TYPE"].values).all()
    # Counts integer-like (spot check)
    sample = adata_counts.X[:500, :200]
    if sp.issparse(sample):
        sample = sample.toarray()
    nz = sample[sample != 0]
    assert np.all(nz >= 0) and np.all(nz == np.round(nz)), "Counts not integer-like"
    # LGR5 unique in counts
    counts_vn = list(adata_counts.var_names)
    assert counts_vn.count("LGR5") == 1, f"LGR5 count: {counts_vn.count('LGR5')}"
    print("  Barcodes aligned, cE01=61953 both, PID/tissue match, counts integer, LGR5 unique. PASS")

    # ── Extract cE01 data ─────────────────────────────────────────────
    cE01_barcodes = adata_norm.obs_names[norm_cE01]
    obs_cE01 = adata_norm.obs.loc[cE01_barcodes].copy()
    print(f"\nExtracting cE01 expression matrix...")
    X_cE01 = adata_norm[cE01_barcodes].X
    if sp.issparse(X_cE01):
        X_cE01 = X_cE01.tocsc()

    # ── Compute per-cell scores (full cE01) ───────────────────────────
    print("Computing ISC-8 and FAO-12 per-cell scores on full cE01...")
    isc8_scores = fixed_bg_score(X_cE01, var_names, isc_8, isc_controls)
    fao12_scores = fixed_bg_score(X_cE01, var_names, fao_12, fao_controls)

    obs_cE01["ISC8"] = isc8_scores
    obs_cE01["FAO12"] = fao12_scores

    # ══════════════════════════════════════════════════════════════════
    # GUARD 3: Gate 3 reproduction
    # ══════════════════════════════════════════════════════════════════
    print("\n--- Guard 3: Gate 3 reproduction ---")
    # Donor-aware on full cE01, 34 primary PIDs
    full34 = obs_cE01[obs_cE01[pid_col].isin(primary_34)].copy()
    pp_full34 = compute_per_patient_corr(full34, "ISC8", "FAO12", pid_col, tissue_col,
                                         tumor_code, normal_code)
    paired34 = [r for r in pp_full34 if r["informative"]]
    assert len(paired34) == 34, f"Expected 34 informative, got {len(paired34)}"
    dz34 = np.array([r["delta_z"] for r in paired34])
    repro_mean_dz = float(np.mean(dz34))

    # Per-patient comparison to Gate 3
    max_pp_diff = 0.0
    for r in paired34:
        pid = r["PID"]
        if pid in g3_pp_lgr5.index:
            ref_dz = g3_pp_lgr5.loc[pid, "delta_z"]
            diff = abs(r["delta_z"] - ref_dz)
            max_pp_diff = max(max_pp_diff, diff)

    mean_dz_diff = abs(repro_mean_dz - g3_ref_mean_dz)
    print(f"  Reproduced mean dz: {repro_mean_dz:.6f} (Gate 3 ref: {g3_ref_mean_dz})")
    print(f"  Mean dz diff: {mean_dz_diff:.2e}")
    print(f"  Max per-patient dz diff: {max_pp_diff:.2e}")

    if mean_dz_diff > 1e-4 or max_pp_diff > 1e-9:
        print(f"  FAIL: reproduction tolerance exceeded")
        print("STOP: Guard 3 failed.")
        return
    print("  Guard 3: PASS")

    # ══════════════════════════════════════════════════════════════════
    # LGR5 detection and eligibility (A008 G8-2)
    # ══════════════════════════════════════════════════════════════════
    print("\n--- LGR5 Detection & Eligibility ---")
    lgr5_idx = counts_vn.index("LGR5")
    cE01_indices = np.where(counts_cE01.values)[0]
    lgr5_raw = adata_counts.X[cE01_indices, lgr5_idx]
    if sp.issparse(lgr5_raw):
        lgr5_raw = np.asarray(lgr5_raw.todense()).ravel()
    else:
        lgr5_raw = np.asarray(lgr5_raw).ravel()

    obs_cE01["LGR5_raw"] = lgr5_raw
    obs_cE01["LGR5_detected"] = lgr5_raw > 0
    n_lgr5_detected = obs_cE01["LGR5_detected"].sum()
    print(f"  LGR5-detected cE01 cells: {n_lgr5_detected} / {len(obs_cE01)}")

    # Per-PID per-tissue LGR5-detected counts
    lgr5_counts = obs_cE01[obs_cE01["LGR5_detected"]].groupby([pid_col, tissue_col]).size().unstack(fill_value=0)
    if tumor_code not in lgr5_counts.columns:
        lgr5_counts[tumor_code] = 0
    if normal_code not in lgr5_counts.columns:
        lgr5_counts[normal_code] = 0

    # Eligible: >=15 in BOTH tissues
    eligible_mask = (lgr5_counts[tumor_code] >= 15) & (lgr5_counts[normal_code] >= 15)
    eligible_pids = sorted(lgr5_counts.index[eligible_mask])
    n_eligible = len(eligible_pids)
    overlap_34 = set(eligible_pids) & primary_34
    print(f"  Eligible PIDs (>=15 LGR5+ both tissues): {n_eligible}")
    print(f"  Overlap with primary 34: {len(overlap_34)}")

    # ══════════════════════════════════════════════════════════════════
    # Main analysis: LGR5+ subset (A008 G8-4)
    # ══════════════════════════════════════════════════════════════════
    print("\n--- Main Analysis: LGR5+ Subset ---")
    lgr5_subset = obs_cE01[(obs_cE01[pid_col].isin(eligible_pids)) & obs_cE01["LGR5_detected"]].copy()
    pp_subset = compute_per_patient_corr(lgr5_subset, "ISC8", "FAO12", pid_col, tissue_col,
                                         tumor_code, normal_code)

    # Also compute full-cE01 comparator for same eligible PIDs
    full_elig = obs_cE01[obs_cE01[pid_col].isin(eligible_pids)].copy()
    pp_full = compute_per_patient_corr(full_elig, "ISC8", "FAO12", pid_col, tissue_col,
                                       tumor_code, normal_code)
    pp_full_dict = {r["PID"]: r for r in pp_full}

    # Build output table
    per_patient_rows = []
    for r in pp_subset:
        pid = r["PID"]
        full_r = pp_full_dict.get(pid, {})
        per_patient_rows.append({
            "PID": pid,
            "in_primary_34": pid in primary_34,
            "n_lgr5_tumor": int(lgr5_counts.loc[pid, tumor_code]) if pid in lgr5_counts.index else 0,
            "n_lgr5_normal": int(lgr5_counts.loc[pid, normal_code]) if pid in lgr5_counts.index else 0,
            "r_normal_subset": r.get("r_normal", np.nan),
            "r_tumor_subset": r.get("r_tumor", np.nan),
            "z_normal_subset": r.get("z_normal", np.nan),
            "z_tumor_subset": r.get("z_tumor", np.nan),
            "delta_z_subset": r.get("delta_z", np.nan),
            "informative": r.get("informative", False),
            "exclusion_reason": r.get("exclusion_reason", ""),
            "r_normal_full": full_r.get("r_normal", np.nan),
            "r_tumor_full": full_r.get("r_tumor", np.nan),
            "delta_z_full": full_r.get("delta_z", np.nan),
        })

    pp_out = pd.DataFrame(per_patient_rows)
    pp_out.to_csv(TABLES_DIR / "gate6_per_patient.csv", index=False)
    print(f"Saved: tables/gate6_per_patient.csv ({len(pp_out)} rows)")

    # Paired informative patients
    retained = pp_out[pp_out["informative"] == True]
    n_retained = len(retained)
    print(f"  Retained (both tissues informative): {n_retained}")

    # ══════════════════════════════════════════════════════════════════
    # Outcome classification (A008 G8-5)
    # ══════════════════════════════════════════════════════════════════
    if n_retained < 15:
        outcome = "POWER-LIMITED / INCONCLUSIVE"
        subset_stats = {"n": n_retained, "mean_delta_z": None, "median_delta_z": None,
                        "wilcoxon_p": None, "bootstrap_ci_lo": None, "bootstrap_ci_hi": None,
                        "ci_excludes_0": None, "loo_all_positive": None}
        if n_retained > 0:
            dz_sub = retained["delta_z_subset"].values
            subset_stats["mean_delta_z"] = float(np.mean(dz_sub))
            subset_stats["median_delta_z"] = float(np.median(dz_sub))
    else:
        dz_sub = retained["delta_z_subset"].values
        subset_stats = donor_aware_stats(dz_sub, seed=SEED, n_boot=N_BOOT)
        mean_dz = subset_stats["mean_delta_z"]

        if mean_dz <= 0:
            outcome = "DIRECTION NOT PRESERVED"
        elif (subset_stats["wilcoxon_p"] is not None and subset_stats["wilcoxon_p"] < 0.05
              and subset_stats["ci_excludes_0"]
              and subset_stats["loo_all_positive"]):
            outcome = "SUPPORTED"
        else:
            outcome = "NOT ESTABLISHED"

    print(f"\n  OUTCOME: {outcome}")
    if subset_stats.get("mean_delta_z") is not None:
        print(f"  Mean dz (subset): {subset_stats['mean_delta_z']:.6f}")
        print(f"  Median dz (subset): {subset_stats['median_delta_z']:.6f}")
    if subset_stats.get("wilcoxon_p") is not None:
        print(f"  Wilcoxon p: {subset_stats['wilcoxon_p']:.2e}")
        print(f"  Bootstrap CI: [{subset_stats['bootstrap_ci_lo']:.4f}, {subset_stats['bootstrap_ci_hi']:.4f}]")
        print(f"  CI excludes 0: {subset_stats['ci_excludes_0']}")
        print(f"  LOO all positive: {subset_stats['loo_all_positive']}")
        print(f"  Prop dz>0: {subset_stats.get('proportion_delta_z_positive', 'N/A')}")

    # ── Descriptive comparator ────────────────────────────────────────
    print("\n--- Descriptive Comparator ---")
    retained_pids = set(retained["PID"].values)
    full_retained = [r for r in pp_full if r["PID"] in retained_pids and r["informative"]]
    if len(full_retained) > 0:
        dz_full_comp = np.array([r["delta_z"] for r in full_retained])
        mean_dz_full = float(np.mean(dz_full_comp))
        if subset_stats.get("mean_delta_z") is not None and abs(mean_dz_full) > 0:
            retention_ratio = abs(subset_stats["mean_delta_z"]) / abs(mean_dz_full)
        else:
            retention_ratio = np.nan
    else:
        mean_dz_full = np.nan
        retention_ratio = np.nan

    print(f"  Mean dz (full cE01, same PIDs): {mean_dz_full:.6f}" if np.isfinite(mean_dz_full) else "  N/A")
    print(f"  Retention ratio |subset|/|full|: {retention_ratio:.4f}" if np.isfinite(retention_ratio) else "  N/A")
    print(f"  Gate 3 isc_loo_LGR5 context (34 PIDs): {g3_ref_mean_dz}")

    # ══════════════════════════════════════════════════════════════════
    # Terminal summary
    # ══════════════════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("TERMINAL SUMMARY")
    print(f"{'='*70}")
    print(f"Guard 1 (controls): PASS (300/400, no dups, no overlap, LGR5 absent)")
    print(f"Guard 2 (alignment): PASS (61953 cE01 both files, barcodes/PID/tissue match)")
    print(f"Guard 3 (Gate 3 repro): PASS (mean dz={repro_mean_dz:.6f}, ref={g3_ref_mean_dz}, "
          f"diff={mean_dz_diff:.2e}, max pp diff={max_pp_diff:.2e})")
    print(f"Eligible PIDs: {n_eligible} (overlap with 34: {len(overlap_34)})")
    print(f"Retained (both informative): {n_retained}")
    print(f"OUTCOME: {outcome}")
    if subset_stats.get("mean_delta_z") is not None:
        print(f"  mean dz (subset) = {subset_stats['mean_delta_z']:.6f}")
        print(f"  median dz (subset) = {subset_stats['median_delta_z']:.6f}")
    if subset_stats.get("wilcoxon_p") is not None:
        print(f"  Wilcoxon p = {subset_stats['wilcoxon_p']:.2e}")
        print(f"  Bootstrap CI = [{subset_stats['bootstrap_ci_lo']:.4f}, {subset_stats['bootstrap_ci_hi']:.4f}]")
        print(f"  CI excludes 0 = {subset_stats['ci_excludes_0']}")
        print(f"  LOO all positive = {subset_stats['loo_all_positive']}")
    if np.isfinite(retention_ratio):
        print(f"  Comparator retention = {retention_ratio:.4f}")
    print(f"  Gate 3 context (34 PIDs) = {g3_ref_mean_dz}")
    print(f"\nFraming: LGR5-transcript-detected cE01 subset; dropout-affected proxy;")
    print(f"  no mechanistic/causal LGR5+ stem-compartment claim.")
    print(f"Ceiling: correlational; tissue/batch Cramer's V=1.0; external replication required.")
    print(f"A9 (Gate 4+5) is unaffected by this gate.")

    # ══════════════════════════════════════════════════════════════════
    # Write report
    # ══════════════════════════════════════════════════════════════════
    lines = [
        "# Gate 6 -- LGR5 Sensitivity Report", "",
        f"**Protocol:** protocol_frozen_v9 (A008)",
        f"**Versions:** {versions}",
        f"**Seed:** {SEED}, **Bootstrap:** {N_BOOT}", "",
        "## Pre-Analysis Guards", "",
        f"1. **Controls:** ISC=300, FAO=400, no dups, zero target overlap, LGR5 absent. **PASS**",
        f"2. **Alignment:** 61953 cE01 both files, barcodes match, PID/tissue match, counts integer. **PASS**",
        f"3. **Gate 3 reproduction:** mean dz={repro_mean_dz:.6f} (ref={g3_ref_mean_dz}, "
        f"diff={mean_dz_diff:.2e}), max per-patient diff={max_pp_diff:.2e}. **PASS**", "",
        "## LGR5 Detection & Eligibility", "",
        f"- LGR5-detected cE01 cells: {n_lgr5_detected} / 61953",
        f"- Eligible PIDs (>=15 LGR5+ both tissues): **{n_eligible}**",
        f"- Overlap with primary 34: {len(overlap_34)}", "",
        "## Main Result (LGR5+ Subset)", "",
        f"- Retained paired PIDs: {n_retained}",
        f"- **Outcome: {outcome}**",
    ]
    if subset_stats.get("mean_delta_z") is not None:
        lines.extend([
            f"- Mean dz: {subset_stats['mean_delta_z']:.6f}",
            f"- Median dz: {subset_stats['median_delta_z']:.6f}",
        ])
    if subset_stats.get("wilcoxon_p") is not None:
        lines.extend([
            f"- Wilcoxon p: {subset_stats['wilcoxon_p']:.2e}",
            f"- Bootstrap 95% CI: [{subset_stats['bootstrap_ci_lo']:.4f}, {subset_stats['bootstrap_ci_hi']:.4f}]",
            f"- CI excludes 0: {subset_stats['ci_excludes_0']}",
            f"- LOO all positive: {subset_stats['loo_all_positive']}",
            f"- Proportion dz>0: {subset_stats.get('proportion_delta_z_positive', 'N/A')}",
        ])
    lines.extend(["", "## Descriptive Comparator (NOT a gate)", ""])
    if np.isfinite(mean_dz_full):
        lines.append(f"- Mean dz (full cE01, same {n_retained} PIDs): {mean_dz_full:.6f}")
    if np.isfinite(retention_ratio):
        lines.append(f"- Retention ratio |subset|/|full|: {retention_ratio:.4f}")
    lines.extend([
        f"- Gate 3 isc_loo_LGR5 (34 PIDs, full cE01): {g3_ref_mean_dz}", "",
        "## Framing & Caveats", "",
        "- LGR5-transcript-detected = raw count > 0; a dropout-affected proxy.",
        "- No mechanistic/causal LGR5+ stem-compartment claim is licensed.",
        "- **Batch ceiling (E11):** tissue and batch perfectly confounded (Cramer's V=1.0).",
        "- External replication (Lee et al.) mandatory for any strong biological claim.",
        "- **A9 (Gate 4 + Gate 5) is unaffected** by this gate's outcome.",
        f"- Gate 6 does NOT re-open the core survival determination.",
    ])
    with open(RESULTS_DIR / "gate6_lgr5_sensitivity_report.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nSaved: results/gate6_lgr5_sensitivity_report.md")

    # ══════════════════════════════════════════════════════════════════
    # Save provenance
    # ══════════════════════════════════════════════════════════════════
    prov = {
        "script": "06_gate6_lgr5_sensitivity.py",
        "gate": "Gate 6",
        "versions": versions,
        "seed": SEED,
        "n_bootstrap": N_BOOT,
        "n_ce01": 61953,
        "n_lgr5_detected_ce01": int(n_lgr5_detected),
        "n_eligible_pids": n_eligible,
        "eligible_pids": eligible_pids,
        "overlap_with_primary_34": len(overlap_34),
        "n_retained_paired": n_retained,
        "outcome": outcome,
        "guard_results": {
            "guard1_controls": "PASS",
            "guard2_alignment": "PASS",
            "guard3_gate3_reproduction": "PASS",
            "guard3_mean_dz_diff": mean_dz_diff,
            "guard3_max_pp_diff": max_pp_diff,
        },
        "subset_stats": {k: v for k, v in subset_stats.items() if k != "n"},
        "comparator_mean_dz_full": mean_dz_full if np.isfinite(mean_dz_full) else None,
        "comparator_retention_ratio": retention_ratio if np.isfinite(retention_ratio) else None,
        "gate3_context_mean_dz": g3_ref_mean_dz,
        "inputs": {
            "epithelial_normalized": str(h5ad_path),
            "epithelial_counts": str(counts_path),
            "clean_control_genes": str(clean_ctrl_path),
            "config": str(CV2 / "config.yaml"),
        },
        "input_hashes": {
            "clean_control_gene_sets_json": clean_ctrl_hash,
            "config_yaml": sha256_file(CV2 / "config.yaml"),
        },
        "status": "COMPLETE",
    }
    with open(PROV_DIR / "06_gate6_lgr5_sensitivity.json", "w") as f:
        json.dump(prov, f, indent=2)
    print(f"Saved: provenance/06_gate6_lgr5_sensitivity.json")

    print("\nDone.")


if __name__ == "__main__":
    main()
