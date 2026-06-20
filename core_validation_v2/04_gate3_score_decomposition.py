"""
Gate 3 — Donor-Aware Score Decomposition
==========================================
Protocol: protocol_frozen_v8 (v1 + A001-A007)
Execution: after Gate 4 PASS + Gate 5 PASS

Gate 3 is interpretation/decomposition — NOT a core hard kill.
The A9 survival (Gate 4 + Gate 5) is already met.
Gate 3 determines whether:
  - The FAO-specific framing holds (or should reframe to differentiation/lipid)
  - Either score is single-gene-driven (program-level robustness)

Variants (22 + baseline):
  - baseline: full ISC-9 x full FAO-12 (clean fixed-bg)
  - 9 ISC-LOO: each ISC gene removed in turn
  - 12 FAO-LOO: each FAO gene removed in turn (includes fabp1/hmgcs2/angptl4)
  - 1 core beta-oxidation sub-score: ISC-9 x core-beta-ox-5

Material weakening = direction lost OR |mean dz| < 50% of baseline.
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

# ── Paths ──────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[1]
CV2 = Path(__file__).resolve().parent
RESULTS_DIR = CV2 / "results"
TABLES_DIR = RESULTS_DIR / "tables"
PROV_DIR = CV2 / "provenance"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
N_BOOT = 10000

# ═══════════════════════════════════════════════════════════════════════
# Copied helpers (NOT imported — source scripts are not import-safe)
# Source: 04_gate4_donor_aware.py (lines 43-198) SHA-256 of that file
#         recorded in provenance.
# Source: 02_gate2_recompute_scores.py (lines 223-249) SHA-256 recorded.
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


def sha256_str(s: str) -> str:
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


def compute_per_patient(df, isc_col, fao_col, method_name, pid_col, tissue_col,
                        tumor_code, normal_code):
    rows = []
    for pid in sorted(df[pid_col].unique()):
        pdata = df[df[pid_col] == pid]
        rec = {"method": method_name, "PID": pid}
        for tissue, code, prefix in [("normal", normal_code, "normal"),
                                     ("tumor", tumor_code, "tumor")]:
            tdata = pdata[pdata[tissue_col] == code]
            n = len(tdata)
            rec[f"n_{tissue}"] = n
            if n == 0:
                rec[f"{prefix}_status"] = "uninformative"
                rec[f"{prefix}_exclusion_reason"] = "no cells"
                rec[f"r_{prefix}"] = np.nan
                continue
            isc = tdata[isc_col].values.astype(np.float64)
            fao = tdata[fao_col].values.astype(np.float64)
            r, status, reason = safe_pearsonr(isc, fao)
            rec[f"{prefix}_status"] = status
            rec[f"{prefix}_exclusion_reason"] = reason
            rec[f"r_{prefix}"] = r
        if (rec.get("normal_status") == "ok" and rec.get("tumor_status") == "ok"):
            rec["retained_for_paired_analysis"] = True
            rn, _ = clip_r(rec["r_normal"])
            rt, _ = clip_r(rec["r_tumor"])
            rec["z_normal"] = fisher_z(rn)
            rec["z_tumor"] = fisher_z(rt)
            rec["Delta_z"] = rec["z_tumor"] - rec["z_normal"]
        else:
            rec["retained_for_paired_analysis"] = False
            rec["z_normal"] = np.nan
            rec["z_tumor"] = np.nan
            rec["Delta_z"] = np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def donor_aware_test(paired_df, seed=SEED, n_boot=N_BOOT):
    dz = paired_df["Delta_z"].values.astype(np.float64)
    n = len(dz)
    mean_dz = float(np.mean(dz))
    median_dz = float(np.median(dz))
    n_zeros = int(np.sum(dz == 0))
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
    loo_rows = []
    for i in range(n):
        loo_dz = np.delete(dz, i)
        loo_mean = float(np.mean(loo_dz))
        try:
            _, loo_p = stats.wilcoxon(loo_dz, alternative="two-sided",
                                      zero_method="wilcox", method="auto")
        except ValueError:
            loo_p = np.nan
        sign_pos = loo_mean > 0
        rel_change = abs(loo_mean - mean_dz) / abs(mean_dz) * 100 if abs(mean_dz) > 0 else 0.0
        loo_rows.append({
            "excluded_PID": paired_df.iloc[i]["PID"],
            "n_remaining": n - 1,
            "mean_delta_z": round(loo_mean, 6),
            "wilcoxon_p": round(loo_p, 6) if np.isfinite(loo_p) else None,
            "sign_positive": sign_pos,
            "relative_effect_change_pct": round(rel_change, 2),
        })
    all_same_sign = all(r["sign_positive"] == (mean_dz > 0) for r in loo_rows)
    max_rel_change = max(r["relative_effect_change_pct"] for r in loo_rows)
    return {
        "n": n,
        "mean_delta_z": mean_dz,
        "median_delta_z": median_dz,
        "wilcoxon_stat": float(wstat) if np.isfinite(wstat) else None,
        "wilcoxon_p": float(wp) if np.isfinite(wp) else None,
        "n_exact_zeros": n_zeros,
        "bootstrap_ci_lo": ci_lo,
        "bootstrap_ci_hi": ci_hi,
        "ci_excludes_0": ci_excludes_0,
        "loo_all_same_sign": all_same_sign,
        "loo_max_rel_change_pct": max_rel_change,
        "proportion_delta_z_positive": float(np.mean(dz > 0)),
    }


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
        target_mean = target_mat.mean(axis=1)
    if len(ctrl_idx) == 0:
        return target_mean
    ctrl_mat = X_sparse[:, ctrl_idx]
    if sp.issparse(ctrl_mat):
        ctrl_mean = np.asarray(ctrl_mat.mean(axis=1)).ravel()
    else:
        ctrl_mean = ctrl_mat.mean(axis=1)
    return target_mean - ctrl_mean


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("Gate 3 -- Score Decomposition")
    print("=" * 70)

    versions = get_versions()
    print(f"Versions: {versions}")

    # ── Load config ────────────────────────────────────────────────────
    with open(CV2 / "config.yaml") as f:
        cfg = yaml.safe_load(f)

    isc_9 = cfg["gene_sets"]["isc_9"]
    fao_12 = cfg["gene_sets"]["ppar_fao_12"]
    core_betaox_5 = cfg["gene_sets"]["core_beta_oxidation_5"]
    targets21 = set(isc_9) | set(fao_12)

    pid_col = cfg["population"]["patient_column"]
    tissue_col = cfg["population"]["tissue_column"]
    tumor_code = cfg["population"]["tumor_code"]
    normal_code = cfg["population"]["normal_code"]
    min_cells = cfg["eligibility"]["primary_min_cells"]

    # ── Load clean control lists (A007 G7-2) ──────────────────────────
    clean_ctrl_path = PROV_DIR / "clean_control_gene_sets.json"
    clean_ctrl_hash = sha256_file(clean_ctrl_path)
    with open(clean_ctrl_path) as f:
        clean_ctrl_data = json.load(f)

    isc_controls = clean_ctrl_data["recovered_selected_control_genes"]["isc_controls"]
    fao_controls = clean_ctrl_data["recovered_selected_control_genes"]["fao_controls"]

    # ══════════════════════════════════════════════════════════════════
    # GUARD 1: Counts / uniqueness / presence / hashes
    # ══════════════════════════════════════════════════════════════════
    print("\n--- Guard 1: Counts / uniqueness / hashes ---")
    g1_pass = True
    g1_issues = []

    isc_ctrl_count = len(isc_controls)
    fao_ctrl_count = len(fao_controls)
    isc_ctrl_unique = len(set(isc_controls))
    fao_ctrl_unique = len(set(fao_controls))

    if isc_ctrl_count != 300:
        g1_issues.append(f"ISC controls: expected 300, got {isc_ctrl_count}")
        g1_pass = False
    if fao_ctrl_count != 400:
        g1_issues.append(f"FAO controls: expected 400, got {fao_ctrl_count}")
        g1_pass = False
    if isc_ctrl_unique != isc_ctrl_count:
        g1_issues.append(f"ISC controls have {isc_ctrl_count - isc_ctrl_unique} duplicates")
        g1_pass = False
    if fao_ctrl_unique != fao_ctrl_count:
        g1_issues.append(f"FAO controls have {fao_ctrl_count - fao_ctrl_unique} duplicates")
        g1_pass = False

    isc_ctrl_hash = sha256_str(",".join(sorted(isc_controls)))
    fao_ctrl_hash = sha256_str(",".join(sorted(fao_controls)))

    print(f"  ISC controls: {isc_ctrl_count} (unique: {isc_ctrl_unique}), hash: {isc_ctrl_hash[:16]}...")
    print(f"  FAO controls: {fao_ctrl_count} (unique: {fao_ctrl_unique}), hash: {fao_ctrl_hash[:16]}...")

    if not g1_pass:
        for issue in g1_issues:
            print(f"  FAIL: {issue}")
        print("STOP: Guard 1 failed.")
        return
    print("  Guard 1: PASS")

    # ── Load AnnData (read-only) ──────────────────────────────────────
    h5ad_path = ROOT / cfg["inputs"]["epithelial_normalized"]
    print(f"\nLoading {h5ad_path.name}...")
    adata = ad.read_h5ad(h5ad_path, backed="r")
    var_names = list(adata.var_names)
    var_set = set(var_names)

    # Check all control genes exist in var_names
    isc_ctrl_missing = [g for g in isc_controls if g not in var_set]
    fao_ctrl_missing = [g for g in fao_controls if g not in var_set]
    if isc_ctrl_missing:
        print(f"STOP: Guard 1 addendum -- ISC controls missing from var_names: {isc_ctrl_missing}")
        return
    if fao_ctrl_missing:
        print(f"STOP: Guard 1 addendum -- FAO controls missing from var_names: {fao_ctrl_missing}")
        return

    # ══════════════════════════════════════════════════════════════════
    # GUARD 2: Union-target exclusion
    # ══════════════════════════════════════════════════════════════════
    print("\n--- Guard 2: Union-target exclusion ---")
    isc_ctrl_target_overlap = set(isc_controls) & targets21
    fao_ctrl_target_overlap = set(fao_controls) & targets21
    if isc_ctrl_target_overlap or fao_ctrl_target_overlap:
        print(f"  FAIL: ISC controls intersection targets21 = {isc_ctrl_target_overlap}")
        print(f"  FAIL: FAO controls intersection targets21 = {fao_ctrl_target_overlap}")
        print("STOP: Guard 2 failed.")
        return
    print("  ISC controls intersection targets21 = empty set")
    print("  FAO controls intersection targets21 = empty set")
    print("  Guard 2: PASS")

    # ── Subset to cE01 ────────────────────────────────────────────────
    obs = adata.obs
    cE01_mask = obs[cfg["population"]["subtype_column_short"]] == "cE01"
    stl_mask = obs[cfg["population"]["subtype_column_derived"]] == cfg["population"]["verify_against_label"]
    assert (cE01_mask == stl_mask).all(), "cE01 != Stem/TA-like mismatch"

    cE01_barcodes = obs.index[cE01_mask]
    n_cE01 = len(cE01_barcodes)
    print(f"\ncE01 cells: {n_cE01}")
    assert n_cE01 == 61953, f"Expected 61953 cE01 cells, got {n_cE01}"

    print("Extracting cE01 expression matrix...")
    X_cE01 = adata[cE01_barcodes].X
    if sp.issparse(X_cE01):
        X_cE01 = X_cE01.tocsc()
    obs_cE01 = obs.loc[cE01_barcodes].copy()

    # ══════════════════════════════════════════════════════════════════
    # GUARD 3: Eligibility drift
    # ══════════════════════════════════════════════════════════════════
    print("\n--- Guard 3: Eligibility drift ---")
    gate4_prov_path = PROV_DIR / "04_gate4_donor_aware.json"
    with open(gate4_prov_path) as f:
        gate4_prov = json.load(f)
    gate4_pids = set(gate4_prov["eligible_pids"])

    pid_tissue_counts = obs_cE01.groupby([pid_col, tissue_col]).size().unstack(fill_value=0)
    if tumor_code in pid_tissue_counts.columns and normal_code in pid_tissue_counts.columns:
        eligible_mask = (pid_tissue_counts[tumor_code] >= min_cells) & (pid_tissue_counts[normal_code] >= min_cells)
        recomputed_pids = set(pid_tissue_counts.index[eligible_mask])
    else:
        recomputed_pids = set()

    if recomputed_pids != gate4_pids:
        extra = recomputed_pids - gate4_pids
        missing = gate4_pids - recomputed_pids
        print(f"  FAIL: PID set mismatch. Extra: {extra}, Missing: {missing}")
        print("STOP: NEEDS_HUMAN_CONFIRMATION -- eligibility drift")
        return

    eligible_pids = sorted(gate4_pids)
    print(f"  Recomputed: {len(recomputed_pids)}, Gate 4: {len(gate4_pids)} -- exact match")
    print("  Guard 3: PASS")

    eligible_mask_obs = obs_cE01[pid_col].isin(gate4_pids)
    obs_eligible = obs_cE01[eligible_mask_obs].copy()

    # ── Compute baseline scores ───────────────────────────────────────
    print("\nComputing baseline clean fixed-bg scores on all cE01 cells...")
    isc_baseline_all = fixed_bg_score(X_cE01, var_names, isc_9, isc_controls)
    fao_baseline_all = fixed_bg_score(X_cE01, var_names, fao_12, fao_controls)

    # ══════════════════════════════════════════════════════════════════
    # GUARD 4: Baseline reproduction
    # ══════════════════════════════════════════════════════════════════
    print("\n--- Guard 4: Baseline reproduction ---")
    gate2b_path = RESULTS_DIR / "gate2b_clean_background_scores_cE01.parquet"
    gate2b_df = pd.read_parquet(gate2b_path)
    if gate2b_df.index.name != "barcode":
        if "barcode" in gate2b_df.columns:
            gate2b_df = gate2b_df.set_index("barcode")

    common_barcodes = cE01_barcodes.intersection(gate2b_df.index)
    if len(common_barcodes) != n_cE01:
        print(f"  FAIL: barcode mismatch -- common={len(common_barcodes)}, cE01={n_cE01}")
        print("STOP: Guard 4 failed.")
        return

    gate2b_aligned = gate2b_df.loc[cE01_barcodes]
    isc_ref = gate2b_aligned["ISC_scanpy_cleanbg"].values
    fao_ref = gate2b_aligned["FAO_scanpy_cleanbg"].values

    isc_max_diff = float(np.max(np.abs(isc_baseline_all - isc_ref)))
    fao_max_diff = float(np.max(np.abs(fao_baseline_all - fao_ref)))
    isc_pearson = float(np.corrcoef(isc_baseline_all, isc_ref)[0, 1])
    fao_pearson = float(np.corrcoef(fao_baseline_all, fao_ref)[0, 1])
    isc_spearman = float(stats.spearmanr(isc_baseline_all, isc_ref).statistic)
    fao_spearman = float(stats.spearmanr(fao_baseline_all, fao_ref).statistic)

    print(f"  ISC max|diff|={isc_max_diff:.2e}, Pearson={isc_pearson:.10f}, Spearman={isc_spearman:.10f}")
    print(f"  FAO max|diff|={fao_max_diff:.2e}, Pearson={fao_pearson:.10f}, Spearman={fao_spearman:.10f}")

    if isc_max_diff >= 1e-9 or fao_max_diff >= 1e-9:
        print(f"  FAIL: max diff >= 1e-9 (ISC: {isc_max_diff}, FAO: {fao_max_diff})")
        print("STOP: Guard 4 failed.")
        return
    print("  Guard 4: PASS")

    # Save guard provenance
    guard_prov = {
        "source_file": str(clean_ctrl_path),
        "source_hash": clean_ctrl_hash,
        "isc_controls_count": isc_ctrl_count,
        "fao_controls_count": fao_ctrl_count,
        "isc_controls_hash": isc_ctrl_hash,
        "fao_controls_hash": fao_ctrl_hash,
        "isc_duplicates": 0,
        "fao_duplicates": 0,
        "union_target_exclusion_isc": True,
        "union_target_exclusion_fao": True,
        "scoring_formula": "score = mean(target_genes) - mean(fixed_clean_control_genes) per cell, on log-normalized .X",
        "no_regeneration": "Controls read from saved JSON; no build_fixed_control_genes() or scanpy control selection called",
        "baseline_reproduction": {
            "isc_max_abs_diff": isc_max_diff,
            "fao_max_abs_diff": fao_max_diff,
            "isc_pearson": isc_pearson,
            "fao_pearson": fao_pearson,
            "isc_spearman": isc_spearman,
            "fao_spearman": fao_spearman,
            "gating_condition_met": True,
        },
        "isc_controls": isc_controls,
        "fao_controls": fao_controls,
    }
    with open(PROV_DIR / "gate3_clean_fixed_control_gene_sets.json", "w") as f:
        json.dump(guard_prov, f, indent=2)
    print("\nSaved: provenance/gate3_clean_fixed_control_gene_sets.json")

    # ══════════════════════════════════════════════════════════════════
    # Build variant definitions (22 + baseline = 23)
    # ══════════════════════════════════════════════════════════════════
    variants = []
    variants.append({
        "variant": "baseline", "variant_family": "baseline",
        "removed_gene": None, "isc_genes": list(isc_9), "fao_genes": list(fao_12),
    })
    for gene in isc_9:
        variants.append({
            "variant": f"isc_loo_{gene}", "variant_family": "isc_loo",
            "removed_gene": gene,
            "isc_genes": [g for g in isc_9 if g != gene], "fao_genes": list(fao_12),
        })
    named = {"FABP1": "fabp1_removed", "HMGCS2": "hmgcs2_removed", "ANGPTL4": "angptl4_removed"}
    for gene in fao_12:
        vname = named.get(gene, f"fao_loo_{gene}")
        variants.append({
            "variant": vname, "variant_family": "fao_loo",
            "removed_gene": gene,
            "isc_genes": list(isc_9), "fao_genes": [g for g in fao_12 if g != gene],
        })
    variants.append({
        "variant": "core_beta_oxidation_5", "variant_family": "sub_score",
        "removed_gene": None, "isc_genes": list(isc_9), "fao_genes": list(core_betaox_5),
    })
    print(f"\nTotal variants (incl. baseline): {len(variants)}")
    assert len(variants) == 23, f"Expected 23, got {len(variants)}"

    # ══════════════════════════════════════════════════════════════════
    # Compute per-variant donor-aware effects
    # ══════════════════════════════════════════════════════════════════
    per_patient_rows = []
    variant_effects = []
    baseline_mean_dz = None

    for vi, v in enumerate(variants):
        vname = v["variant"]
        isc_genes = v["isc_genes"]
        fao_genes = v["fao_genes"]

        # Compute variant scores on ALL cE01 cells (then subset to eligible)
        isc_score_all = fixed_bg_score(X_cE01, var_names, isc_genes, isc_controls)
        fao_score_all = fixed_bg_score(X_cE01, var_names, fao_genes, fao_controls)

        score_df = obs_eligible.copy()
        score_df["ISC_variant"] = isc_score_all[eligible_mask_obs.values]
        score_df["FAO_variant"] = fao_score_all[eligible_mask_obs.values]

        pp_df = compute_per_patient(
            score_df, "ISC_variant", "FAO_variant", vname,
            pid_col, tissue_col, tumor_code, normal_code
        )

        for _, row in pp_df.iterrows():
            informative = row.get("retained_for_paired_analysis", False)
            excl_reasons = []
            if row.get("normal_status") == "uninformative":
                excl_reasons.append(f"normal:{row.get('normal_exclusion_reason', '')}")
            if row.get("tumor_status") == "uninformative":
                excl_reasons.append(f"tumor:{row.get('tumor_exclusion_reason', '')}")
            per_patient_rows.append({
                "variant": vname,
                "variant_family": v["variant_family"],
                "removed_gene": v["removed_gene"] or "",
                "PID": row["PID"],
                "n_tumor": int(row.get("n_tumor", 0)),
                "n_normal": int(row.get("n_normal", 0)),
                "r_normal": row.get("r_normal", np.nan),
                "r_tumor": row.get("r_tumor", np.nan),
                "z_normal": row.get("z_normal", np.nan),
                "z_tumor": row.get("z_tumor", np.nan),
                "delta_z": row.get("Delta_z", np.nan),
                "informative": informative,
                "exclusion_reason": "; ".join(excl_reasons) if excl_reasons else "",
            })

        paired = pp_df[pp_df["retained_for_paired_analysis"] == True].copy()
        n_informative = len(paired)
        informative_pids = set(paired["PID"].values)
        excluded_pids = gate4_pids - informative_pids
        pid_set_matches = informative_pids == gate4_pids

        if not pid_set_matches:
            print(f"\nSTOP: NEEDS_HUMAN_CONFIRMATION -- variant '{vname}' informative PID set "
                  f"differs from baseline. Excluded: {excluded_pids}")
            return

        dat = donor_aware_test(paired, seed=SEED, n_boot=N_BOOT)

        if vname == "baseline":
            baseline_mean_dz = dat["mean_delta_z"]
            print(f"\nBaseline clean fixed-bg mean dz: {baseline_mean_dz:.6f}")

        if baseline_mean_dz is not None and baseline_mean_dz != 0:
            retention_ratio = abs(dat["mean_delta_z"]) / abs(baseline_mean_dz)
        else:
            retention_ratio = np.nan
        direction_preserved = (dat["mean_delta_z"] > 0) == (baseline_mean_dz > 0) if baseline_mean_dz is not None else True
        material_weakening = (not direction_preserved) or (retention_ratio < 0.50)

        if vname == "baseline":
            interp_flag = "BASELINE"
        elif material_weakening:
            interp_flag = "MATERIALLY_WEAKENED"
        else:
            interp_flag = "ROBUST"

        variant_effects.append({
            "variant": vname,
            "variant_family": v["variant_family"],
            "removed_gene": v["removed_gene"] or "",
            "isc_target_genes": ",".join(isc_genes),
            "fao_target_genes": ",".join(fao_genes),
            "n_patients": n_informative,
            "informative_pid_set_matches_baseline": pid_set_matches,
            "excluded_pid_count": len(excluded_pids),
            "mean_delta_z": round(dat["mean_delta_z"], 6),
            "median_delta_z": round(dat["median_delta_z"], 6),
            "baseline_mean_delta_z": round(baseline_mean_dz, 6) if baseline_mean_dz is not None else None,
            "retention_ratio": round(retention_ratio, 4) if np.isfinite(retention_ratio) else None,
            "direction_preserved": direction_preserved,
            "material_weakening": material_weakening,
            "wilcoxon_p_descriptive": dat["wilcoxon_p"],
            "bootstrap_ci_low_descriptive": round(dat["bootstrap_ci_lo"], 6),
            "bootstrap_ci_high_descriptive": round(dat["bootstrap_ci_hi"], 6),
            "loo_sign_stable_descriptive": dat["loo_all_same_sign"],
            "interpretation_flag": interp_flag,
        })

        if vi > 0 and vi % 5 == 0:
            print(f"  Computed {vi}/{len(variants)-1} variants...")

    print(f"  All {len(variants)} variants computed.")

    # ══════════════════════════════════════════════════════════════════
    # Save tables
    # ══════════════════════════════════════════════════════════════════
    pp_out = pd.DataFrame(per_patient_rows)
    pp_out.to_csv(TABLES_DIR / "gate3_per_patient_variant.csv", index=False)
    print(f"\nSaved: tables/gate3_per_patient_variant.csv ({len(pp_out)} rows)")

    ve_df = pd.DataFrame(variant_effects)
    ve_df.to_csv(TABLES_DIR / "gate3_variant_effects.csv", index=False)
    print(f"Saved: tables/gate3_variant_effects.csv ({len(ve_df)} rows)")

    # ══════════════════════════════════════════════════════════════════
    # Interpretation
    # ══════════════════════════════════════════════════════════════════
    scanpy_primary_dz = 0.400629
    clean_bg_dz = 0.386688

    fabp1_row = ve_df[ve_df["variant"] == "fabp1_removed"].iloc[0]
    hmgcs2_row = ve_df[ve_df["variant"] == "hmgcs2_removed"].iloc[0]
    betaox_row = ve_df[ve_df["variant"] == "core_beta_oxidation_5"].iloc[0]
    angptl4_row = ve_df[ve_df["variant"] == "angptl4_removed"].iloc[0]

    fao_reframe_needed = (
        bool(fabp1_row["material_weakening"]) and
        bool(hmgcs2_row["material_weakening"]) and
        bool(betaox_row["material_weakening"])
    )

    isc_loo_weak = ve_df[(ve_df["variant_family"] == "isc_loo") & (ve_df["material_weakening"] == True)]
    fao_loo_weak = ve_df[(ve_df["variant_family"] == "fao_loo") & (ve_df["material_weakening"] == True)]
    isc_single_gene = len(isc_loo_weak) > 0
    fao_single_gene = len(fao_loo_weak) > 0
    isc_problem = list(isc_loo_weak["removed_gene"].values) if isc_single_gene else []
    fao_problem = list(fao_loo_weak["removed_gene"].values) if fao_single_gene else []

    if fao_reframe_needed:
        interpretation = ("REFRAME: FAO-specific framing not supported; reframe as "
                          "differentiation/lipid-handling-state relationship (Outcome 2)")
    elif isc_single_gene or fao_single_gene:
        parts = []
        if isc_single_gene:
            parts.append(f"ISC program-level not robust (driven by: {', '.join(isc_problem)})")
        if fao_single_gene:
            parts.append(f"FAO program-level not robust (driven by: {', '.join(fao_problem)})")
        interpretation = "SINGLE-GENE-DRIVEN: " + "; ".join(parts)
    else:
        interpretation = ("FAO framing retained; program-level framing robust for both scores. "
                          "Outcome 1 not declared -- Gate 6/7 pending.")

    # ── Terminal summary ──────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("INTERPRETATION")
    print(f"{'='*70}")
    print(f"Baseline mean dz: {baseline_mean_dz:.6f}")
    print(f"Context: official scanpy = {scanpy_primary_dz}, clean-bg Gate 4 = {clean_bg_dz}")
    print(f"50% threshold: {abs(baseline_mean_dz) * 0.5:.6f}")
    print(f"\nFAO reframe condition (A007 G7-5):")
    print(f"  fabp1_removed:  mw={fabp1_row['material_weakening']}, ret={fabp1_row['retention_ratio']}")
    print(f"  hmgcs2_removed: mw={hmgcs2_row['material_weakening']}, ret={hmgcs2_row['retention_ratio']}")
    print(f"  core_beta_ox_5: mw={betaox_row['material_weakening']}, ret={betaox_row['retention_ratio']}")
    print(f"  angptl4 (info): mw={angptl4_row['material_weakening']}, ret={angptl4_row['retention_ratio']}")
    print(f"  => FAO reframe needed: {fao_reframe_needed}")
    print(f"\nProgram-level:")
    print(f"  ISC single-gene-driven: {isc_single_gene} {isc_problem}")
    print(f"  FAO single-gene-driven: {fao_single_gene} {fao_problem}")
    print(f"\n=> {interpretation}")

    print(f"\n{'='*70}")
    print("VARIANT SUMMARY")
    print(f"{'='*70}")
    print(f"{'Variant':<30} {'mean_dz':>10} {'retention':>10} {'dir':>5} {'mat_wk':>8}")
    print("-" * 70)
    for _, row in ve_df.iterrows():
        vn = row["variant"][:29]
        mdz = f"{row['mean_delta_z']:.4f}"
        ret = f"{row['retention_ratio']:.4f}" if row["retention_ratio"] is not None else "N/A"
        d = "+" if row["direction_preserved"] else "-"
        mw = "YES" if row["material_weakening"] else "no"
        print(f"{vn:<30} {mdz:>10} {ret:>10} {d:>5} {mw:>8}")

    print(f"\nBatch ceiling (E11): tissue/batch perfectly confounded. External replication mandatory.")
    print(f"Outcome 1 not declared; next = Gate 6.")

    # ══════════════════════════════════════════════════════════════════
    # Write Markdown report
    # ══════════════════════════════════════════════════════════════════
    lines = [
        "# Gate 3 -- Score Decomposition Report", "",
        f"**Protocol:** protocol_frozen_v8 (A007)",
        f"**Versions:** {versions}",
        f"**Seed:** {SEED}, **Bootstrap:** {N_BOOT}",
        f"**Eligible patients:** {len(eligible_pids)}, drift check: PASSED", "",
        "## Pre-Analysis Integrity Guards", "",
        f"1. **Counts/hashes:** ISC controls={isc_ctrl_count}, FAO controls={fao_ctrl_count}. "
        f"No duplicates. All present in var_names. **PASS**",
        f"2. **Union-target exclusion:** ISC and FAO controls have zero intersection with targets21. **PASS**",
        f"3. **Eligibility drift:** 34 PIDs, exact match to Gate 4. **PASS**",
        f"4. **Baseline reproduction:** ISC max|diff|={isc_max_diff:.2e}, FAO max|diff|={fao_max_diff:.2e} "
        f"(both < 1e-9). Pearson ISC={isc_pearson:.10f}, FAO={fao_pearson:.10f}. **PASS**", "",
        "## Baseline", "",
        f"Clean fixed-bg baseline mean dz: **{baseline_mean_dz:.6f}**",
        f"Context: official scanpy primary (0.400629), clean-bg Gate 4 (0.386688).",
        f"50% material-weakening threshold: {abs(baseline_mean_dz) * 0.5:.6f}", "",
        "## Variant Effects", "",
        "| Variant | mean_dz | retention | dir | mat_weak | Wilcoxon p | CI | LOO stable |",
        "|---------|---------|-----------|-----|----------|------------|-----|------------|",
    ]
    for _, row in ve_df.iterrows():
        vn = row["variant"]
        mdz = f"{row['mean_delta_z']:.4f}"
        ret = f"{row['retention_ratio']:.4f}" if row["retention_ratio"] is not None else "N/A"
        d = "+" if row["direction_preserved"] else "-"
        mw = "YES" if row["material_weakening"] else "no"
        wp = f"{row['wilcoxon_p_descriptive']:.2e}" if row["wilcoxon_p_descriptive"] is not None else "N/A"
        ci = f"[{row['bootstrap_ci_low_descriptive']:.4f}, {row['bootstrap_ci_high_descriptive']:.4f}]"
        loo = "yes" if row["loo_sign_stable_descriptive"] else "no"
        lines.append(f"| {vn} | {mdz} | {ret} | {d} | {mw} | {wp} | {ci} | {loo} |")

    lines.extend(["", "## FAO-Interpretation Evaluation (A007 G7-5)", "",
        f"- fabp1_removed: material_weakening={fabp1_row['material_weakening']} (retention={fabp1_row['retention_ratio']})",
        f"- hmgcs2_removed: material_weakening={hmgcs2_row['material_weakening']} (retention={hmgcs2_row['retention_ratio']})",
        f"- core_beta_oxidation_5: material_weakening={betaox_row['material_weakening']} (retention={betaox_row['retention_ratio']})",
        f"- angptl4_removed (reported, not in condition): material_weakening={angptl4_row['material_weakening']} (retention={angptl4_row['retention_ratio']})", "",
        f"**Conjunction (ALL required for reframe): {fao_reframe_needed}**", "",
        "## Program-Level Evaluation", "",
        f"- ISC single-gene-driven: {isc_single_gene}" + (f" (genes: {', '.join(isc_problem)})" if isc_single_gene else ""),
        f"- FAO single-gene-driven: {fao_single_gene}" + (f" (genes: {', '.join(fao_problem)})" if fao_single_gene else ""), "",
        "## Interpretation", "",
        f"**{interpretation}**", "",
        "## Mandatory Statements", "",
        "**Batch claim ceiling (E11):** Tissue and batch are perfectly confounded",
        "(Cramer's V=1.0; 0/34 patients share T/N batch). External replication mandatory.", "",
        "**Pseudoreplication caveat (E0):** Each cell is not an independent replicate.", "",
        "**Scope:** Gate 3 is interpretation/decomposition. The A9 core survival",
        "(Gate 4 PASS + Gate 5 PASS) is not re-opened. Outcome 1 is not declared",
        "complete; Gate 6 and Gate 7 remain.",
    ])
    with open(RESULTS_DIR / "gate3_score_decomposition_report.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nSaved: results/gate3_score_decomposition_report.md")

    # ══════════════════════════════════════════════════════════════════
    # Save provenance
    # ══════════════════════════════════════════════════════════════════
    gate4_script_hash = sha256_file(CV2 / "04_gate4_donor_aware.py")
    gate2_script_hash = sha256_file(CV2 / "02_gate2_recompute_scores.py")

    prov = {
        "script": "04_gate3_score_decomposition.py",
        "gate": "Gate 3",
        "versions": versions,
        "seed": SEED,
        "n_bootstrap": N_BOOT,
        "n_eligible": len(eligible_pids),
        "eligible_pids": eligible_pids,
        "n_ce01": n_cE01,
        "n_variants": len(variants),
        "copied_helpers": {
            "donor_aware_from": "04_gate4_donor_aware.py",
            "donor_aware_source_hash": gate4_script_hash,
            "fixed_bg_score_from": "02_gate2_recompute_scores.py",
            "fixed_bg_score_source_hash": gate2_script_hash,
        },
        "guard_results": {
            "guard1_counts_hashes": "PASS",
            "guard2_union_target_exclusion": "PASS",
            "guard3_eligibility_drift": "PASS",
            "guard4_baseline_reproduction": "PASS",
            "guard4_isc_max_diff": isc_max_diff,
            "guard4_fao_max_diff": fao_max_diff,
        },
        "baseline_mean_delta_z": baseline_mean_dz,
        "fao_reframe_needed": fao_reframe_needed,
        "isc_single_gene_driven": isc_single_gene,
        "fao_single_gene_driven": fao_single_gene,
        "isc_problem_genes": isc_problem,
        "fao_problem_genes": fao_problem,
        "interpretation": interpretation,
        "inputs": {
            "epithelial_normalized": str(h5ad_path),
            "gate2b_scores": str(gate2b_path),
            "clean_control_genes": str(clean_ctrl_path),
            "gate4_provenance": str(gate4_prov_path),
            "config": str(CV2 / "config.yaml"),
        },
        "input_hashes": {
            "clean_control_gene_sets_json": clean_ctrl_hash,
            "config_yaml": sha256_file(CV2 / "config.yaml"),
        },
        "status": "COMPLETE",
    }
    with open(PROV_DIR / "04_gate3_score_decomposition.json", "w") as f:
        json.dump(prov, f, indent=2)
    print(f"Saved: provenance/04_gate3_score_decomposition.json")

    print("\nDone.")


if __name__ == "__main__":
    main()
