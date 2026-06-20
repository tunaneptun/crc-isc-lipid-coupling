"""
Gate 1 — Input Validation & Matrix Resolution
===============================================
Protocol: CORE_VALIDATION_PROTOCOL.md (tag protocol_frozen_v1) +
          CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md (tag protocol_frozen_v2)
Execution order: 1 (runs first)

Inputs
------
- data/processed/pelka_epithelial.h5ad (read-only)
- data/processed/pelka_epithelial_counts.h5ad (read-only)
- core_validation_v2/config.yaml

Outputs
-------
- core_validation_v2/results/gate1_input_validation.json
- core_validation_v2/provenance/gate1_provenance.json

What this gate does (v1 + Amendment A4, A8)
-------------------------------------------
1. Discover the real per-patient obs column; verify unique-patient count;
   check missing patient IDs; verify no barcode maps to >1 patient;
   verify each patient's T and N cells are correctly paired under one PID;
   verify MMR metadata consistent within patient.
2. Verify cl295v11SubShort == "cE01" matches epithelial_subtype ==
   "Stem/TA-like" exactly.
3. Inspect .X and .raw.X empirically: report what each contains
   (raw counts vs log-normalized) by checking value ranges, integer-likeness.
4. Check cell/gene/patient counts, barcode alignment between normalized
   and counts files, duplicate barcodes/gene names.
5. Report eligible-patient counts at >=20, >=30, >=50 (both tissues), cE01-only.
6. Apply eligibility fallback (Amendment A1): determine primary threshold.
7. Discover batch column(s); produce cross-tabs (patient x tissue, patient x batch,
   tissue x batch); assess tissue-batch confounding (A4).
8. Discover QC columns (total_counts, n_genes_by_counts, pct_counts_mt,
   S_score, G2M_score) for use in A8 sensitivity (A8).

Kill criteria
-------------
- KILL: inputs inconsistent and unreconcilable -> stop; no downstream runs.
"""

import hashlib
import json
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG_PATH = Path(__file__).resolve().parent / "config.yaml"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
PROV_DIR = Path(__file__).resolve().parent / "provenance"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def inspect_matrix(X, label, n_sample=500):
    """Inspect a matrix and report whether it looks like raw counts or log-normalized."""
    info = {"label": label, "shape": list(X.shape), "sparse": sp.issparse(X), "dtype": str(X.dtype)}
    samp = X[:min(n_sample, X.shape[0]), :min(n_sample, X.shape[1])]
    if sp.issparse(samp):
        samp = samp.toarray()
    samp = samp.astype(np.float64)
    nz = samp[samp > 0]
    if len(nz) == 0:
        info["verdict"] = "ALL_ZERO"
        return info
    info["min_nz"] = float(nz.min())
    info["max"] = float(nz.max())
    info["mean_nz"] = float(nz.mean())
    info["median_nz"] = float(np.median(nz))
    frac_int = float((np.abs(nz - np.round(nz)) < 1e-6).mean())
    info["frac_integer"] = frac_int
    if nz.max() > 50 and frac_int > 0.99:
        info["verdict"] = "RAW_COUNTS"
    elif nz.max() < 20 and frac_int < 0.5:
        info["verdict"] = "LOG_NORMALIZED"
    else:
        info["verdict"] = "AMBIGUOUS"
    return info


def main():
    print("=" * 72)
    print("Gate 1 — Input Validation & Matrix Resolution")
    print("=" * 72)

    with open(CFG_PATH) as f:
        cfg = yaml.safe_load(f)

    results = {"gate": "Gate 1", "status": "RUNNING", "checks": {}, "kill": False, "kill_reason": None}

    # ── Load both h5ad files ─────────────────────────────────────────────
    norm_path = ROOT / cfg["inputs"]["epithelial_normalized"]
    counts_path = ROOT / cfg["inputs"]["epithelial_counts"]

    print(f"\nLoading {norm_path.name} ...")
    adata_norm = ad.read_h5ad(str(norm_path))
    print(f"  Shape: {adata_norm.n_obs:,} x {adata_norm.n_vars:,}")

    print(f"Loading {counts_path.name} ...")
    adata_counts = ad.read_h5ad(str(counts_path))
    print(f"  Shape: {adata_counts.n_obs:,} x {adata_counts.n_vars:,}")

    obs = adata_norm.obs.copy()

    # ── 1. Donor identity resolution (Amendment A4) ──────────────────────
    print("\n--- 1. Donor identity resolution ---")

    # Discover patient column
    candidate_pid_cols = ["PID", "patient", "Patient", "donor", "patient_id", "orig_ident"]
    pid_col = None
    for c in candidate_pid_cols:
        if c in obs.columns:
            pid_col = c
            break
    if pid_col is None:
        for c in obs.columns:
            if "patient" in c.lower() or "pid" in c.lower() or "donor" in c.lower():
                pid_col = c
                break

    results["checks"]["patient_column"] = pid_col
    print(f"  Patient column discovered: {pid_col}")

    if pid_col is None:
        results["kill"] = True
        results["kill_reason"] = "No patient column found in obs"
        _save_and_exit(results, norm_path, counts_path, cfg)
        return

    n_patients = obs[pid_col].nunique()
    n_missing_pid = int(obs[pid_col].isna().sum())
    print(f"  Unique patients: {n_patients}")
    print(f"  Missing patient IDs: {n_missing_pid}")
    results["checks"]["n_patients"] = n_patients
    results["checks"]["n_missing_pid"] = n_missing_pid

    # Check no barcode maps to >1 patient (barcodes should be unique)
    dup_barcodes = int(obs.index.duplicated().sum())
    print(f"  Duplicate barcodes: {dup_barcodes}")
    results["checks"]["duplicate_barcodes"] = dup_barcodes

    # Verify each patient's T and N cells paired under one PID
    tissue_col = cfg["population"]["tissue_column"]
    if tissue_col not in obs.columns:
        # Try derived
        for tc in ["tissue_type", "SPECIMEN_TYPE"]:
            if tc in obs.columns:
                tissue_col = tc
                break
    results["checks"]["tissue_column_used"] = tissue_col

    patients_with_both = set()
    patients_tumor_only = set()
    patients_normal_only = set()
    for pid, grp in obs.groupby(pid_col):
        tissues = set(grp[tissue_col].dropna().unique())
        t_code = cfg["population"]["tumor_code"]
        n_code = cfg["population"]["normal_code"]
        has_t = t_code in tissues
        has_n = n_code in tissues
        if has_t and has_n:
            patients_with_both.add(pid)
        elif has_t:
            patients_tumor_only.add(pid)
        elif has_n:
            patients_normal_only.add(pid)

    print(f"  Patients with both T+N: {len(patients_with_both)}")
    print(f"  Patients tumor-only: {len(patients_tumor_only)}")
    print(f"  Patients normal-only: {len(patients_normal_only)}")
    results["checks"]["patients_with_both_tissues"] = len(patients_with_both)
    results["checks"]["patients_tumor_only"] = len(patients_tumor_only)
    results["checks"]["patients_normal_only"] = len(patients_normal_only)

    # MMR consistency within patient
    mmr_col = None
    for c in ["MMRStatus", "mmr_status", "MMR"]:
        if c in obs.columns:
            mmr_col = c
            break
    results["checks"]["mmr_column"] = mmr_col
    mmr_inconsistent = []
    if mmr_col:
        for pid, grp in obs.groupby(pid_col):
            tumor_cells = grp[grp[tissue_col] == cfg["population"]["tumor_code"]]
            if len(tumor_cells) > 0:
                vals = tumor_cells[mmr_col].dropna().unique()
                if len(vals) > 1:
                    mmr_inconsistent.append(str(pid))
        print(f"  MMR inconsistencies within patient: {len(mmr_inconsistent)}")
        if mmr_inconsistent:
            print(f"    Patients: {mmr_inconsistent[:5]}")
    results["checks"]["mmr_inconsistent_patients"] = mmr_inconsistent

    # ── 2. cE01 == Stem/TA-like verification ─────────────────────────────
    print("\n--- 2. cE01 filter verification ---")

    short_col = cfg["population"]["subtype_column_short"]
    derived_col = cfg["population"]["subtype_column_derived"]

    if short_col not in obs.columns:
        results["kill"] = True
        results["kill_reason"] = f"Column {short_col} not found in obs"
        _save_and_exit(results, norm_path, counts_path, cfg)
        return

    ce01_mask = obs[short_col].astype(str) == "cE01"
    n_ce01 = int(ce01_mask.sum())
    print(f"  cl295v11SubShort == 'cE01': {n_ce01:,} cells")

    if derived_col in obs.columns:
        stem_mask = obs[derived_col].astype(str) == cfg["population"]["verify_against_label"]
        n_stem = int(stem_mask.sum())
        print(f"  epithelial_subtype == 'Stem/TA-like': {n_stem:,} cells")

        match_exact = (ce01_mask == stem_mask).all()
        ce01_not_stem = int((ce01_mask & ~stem_mask).sum())
        stem_not_ce01 = int((stem_mask & ~ce01_mask).sum())
        print(f"  Exact match: {match_exact}")
        print(f"  cE01 but not Stem/TA-like: {ce01_not_stem}")
        print(f"  Stem/TA-like but not cE01: {stem_not_ce01}")

        results["checks"]["ce01_count"] = n_ce01
        results["checks"]["stem_ta_count"] = n_stem
        results["checks"]["ce01_stem_exact_match"] = bool(match_exact)
        results["checks"]["ce01_not_stem"] = ce01_not_stem
        results["checks"]["stem_not_ce01"] = stem_not_ce01

        if not match_exact:
            results["kill"] = True
            results["kill_reason"] = (
                f"cE01 filter does not match Stem/TA-like exactly: "
                f"{ce01_not_stem} cE01-not-Stem, {stem_not_ce01} Stem-not-cE01"
            )
            _save_and_exit(results, norm_path, counts_path, cfg)
            return
    else:
        print(f"  WARNING: derived column '{derived_col}' not found; cannot cross-verify")
        results["checks"]["ce01_count"] = n_ce01
        results["checks"]["ce01_stem_exact_match"] = "UNVERIFIABLE"

    # ── 3. Matrix inspection ─────────────────────────────────────────────
    print("\n--- 3. Matrix inspection ---")

    x_info = inspect_matrix(adata_norm.X, ".X (normalized file)")
    print(f"  .X: verdict={x_info['verdict']}, max={x_info.get('max','?')}, "
          f"frac_int={x_info.get('frac_integer','?')}")
    results["checks"]["norm_X"] = x_info

    if adata_norm.raw is not None:
        raw_info = inspect_matrix(adata_norm.raw.X, ".raw.X (normalized file)")
        print(f"  .raw.X: verdict={raw_info['verdict']}, max={raw_info.get('max','?')}, "
              f"frac_int={raw_info.get('frac_integer','?')}")
        results["checks"]["norm_raw_X"] = raw_info
    else:
        print("  .raw: None")
        results["checks"]["norm_raw_X"] = {"verdict": "NOT_PRESENT"}

    counts_x_info = inspect_matrix(adata_counts.X, ".X (counts file)")
    print(f"  counts .X: verdict={counts_x_info['verdict']}, max={counts_x_info.get('max','?')}, "
          f"frac_int={counts_x_info.get('frac_integer','?')}")
    results["checks"]["counts_X"] = counts_x_info

    # ── 4. Barcode alignment, gene names, dimensions ─────────────────────
    print("\n--- 4. Barcode alignment & gene names ---")

    norm_barcodes = set(adata_norm.obs_names)
    counts_barcodes = set(adata_counts.obs_names)
    common = norm_barcodes & counts_barcodes
    norm_only = norm_barcodes - counts_barcodes
    counts_only = counts_barcodes - norm_barcodes
    print(f"  Normalized barcodes: {len(norm_barcodes):,}")
    print(f"  Counts barcodes: {len(counts_barcodes):,}")
    print(f"  Common: {len(common):,}")
    print(f"  Norm-only: {len(norm_only):,}")
    print(f"  Counts-only: {len(counts_only):,}")
    results["checks"]["barcode_alignment"] = {
        "norm": len(norm_barcodes), "counts": len(counts_barcodes),
        "common": len(common), "norm_only": len(norm_only), "counts_only": len(counts_only),
    }

    if len(common) < len(norm_barcodes) * 0.99:
        results["kill"] = True
        results["kill_reason"] = (
            f"Barcode mismatch: only {len(common):,}/{len(norm_barcodes):,} common"
        )
        _save_and_exit(results, norm_path, counts_path, cfg)
        return

    # Gene name duplicates
    norm_gene_dups = int(adata_norm.var_names.duplicated().sum())
    counts_gene_dups = int(adata_counts.var_names.duplicated().sum())
    print(f"  Gene name duplicates (norm): {norm_gene_dups}")
    print(f"  Gene name duplicates (counts): {counts_gene_dups}")
    results["checks"]["gene_duplicates_norm"] = norm_gene_dups
    results["checks"]["gene_duplicates_counts"] = counts_gene_dups

    # Check all 21 target genes present
    all_genes = list(cfg["gene_sets"]["isc_9"]) + list(cfg["gene_sets"]["ppar_fao_12"])
    genes_present = [g for g in all_genes if g in adata_norm.var_names]
    genes_missing = [g for g in all_genes if g not in adata_norm.var_names]
    print(f"  Target genes present: {len(genes_present)}/21")
    if genes_missing:
        print(f"  MISSING: {genes_missing}")
    results["checks"]["target_genes_present"] = len(genes_present)
    results["checks"]["target_genes_missing"] = genes_missing

    # Check which target genes are among duplicated symbols
    if norm_gene_dups > 0:
        dup_names = set(adata_norm.var_names[adata_norm.var_names.duplicated(keep=False)])
        target_in_dups = [g for g in all_genes if g in dup_names]
        print(f"  Target genes in duplicated symbols: {target_in_dups}")
        results["checks"]["target_genes_in_duplicates"] = target_in_dups
    else:
        results["checks"]["target_genes_in_duplicates"] = []

    # Dimensions
    results["checks"]["dimensions"] = {
        "norm_cells": adata_norm.n_obs, "norm_genes": adata_norm.n_vars,
        "counts_cells": adata_counts.n_obs, "counts_genes": adata_counts.n_vars,
    }
    print(f"  Norm: {adata_norm.n_obs:,} x {adata_norm.n_vars:,}")
    print(f"  Counts: {adata_counts.n_obs:,} x {adata_counts.n_vars:,}")

    # ── 5. Eligible patient counts at thresholds ─────────────────────────
    print("\n--- 5. Eligible patient counts (cE01 only) ---")

    ce01_obs = obs[ce01_mask].copy()
    thresholds = [20, 30, 50]
    eligibility = {}

    for thresh in thresholds:
        eligible = []
        for pid, grp in ce01_obs.groupby(pid_col):
            n_t = int((grp[tissue_col] == cfg["population"]["tumor_code"]).sum())
            n_n = int((grp[tissue_col] == cfg["population"]["normal_code"]).sum())
            if n_t >= thresh and n_n >= thresh:
                eligible.append({"patient": str(pid), "n_tumor": n_t, "n_normal": n_n})
        eligibility[f"thresh_{thresh}"] = {
            "threshold": thresh,
            "n_eligible": len(eligible),
            "patients": eligible,
        }
        print(f"  >= {thresh} cells both T+N: {len(eligible)} patients eligible")

    results["checks"]["eligibility"] = eligibility

    # Amendment A1: determine primary threshold
    primary_thresh = cfg["eligibility"]["primary_min_cells"]
    n_at_primary = eligibility[f"thresh_{primary_thresh}"]["n_eligible"]
    power_floor = cfg["eligibility"]["power_contingency_floor"]

    if n_at_primary >= power_floor:
        power_status = "ADEQUATE"
        effective_threshold = primary_thresh
    elif eligibility["thresh_20"]["n_eligible"] >= power_floor:
        power_status = "FALLBACK_20"
        effective_threshold = 20
    else:
        power_status = "UNDERPOWERED"
        effective_threshold = 20

    results["checks"]["power_status"] = power_status
    results["checks"]["effective_threshold"] = effective_threshold
    print(f"  Power status: {power_status} (effective threshold: {effective_threshold})")

    # ── 6. Batch / tissue confound (Amendment A4) ────────────────────────
    print("\n--- 6. Batch/tissue confound assessment ---")

    batch_col = None
    for c in ["batchID", "batch", "Batch", "batch_id", "library_id"]:
        if c in obs.columns:
            batch_col = c
            break
    results["checks"]["batch_column"] = batch_col

    if batch_col:
        n_batches = obs[batch_col].nunique()
        print(f"  Batch column: {batch_col} ({n_batches} unique)")

        # Patient x tissue cross-tab (how many patients have T in one batch, N in another)
        ce01_batch = ce01_obs[[pid_col, tissue_col, batch_col]].copy()

        # Per-patient: count batches per tissue
        confound_rows = []
        for pid, grp in ce01_batch.groupby(pid_col):
            t_batches = set(grp[grp[tissue_col] == cfg["population"]["tumor_code"]][batch_col].dropna().unique())
            n_batches_set = set(grp[grp[tissue_col] == cfg["population"]["normal_code"]][batch_col].dropna().unique())
            shared = t_batches & n_batches_set
            confound_rows.append({
                "patient": str(pid),
                "n_tumor_batches": len(t_batches),
                "n_normal_batches": len(n_batches_set),
                "shared_batches": len(shared),
                "t_n_same_batch": len(shared) > 0,
            })
        df_confound = pd.DataFrame(confound_rows)
        n_same_batch = int(df_confound["t_n_same_batch"].sum())
        n_diff_batch = int((~df_confound["t_n_same_batch"]).sum())
        # Only count patients that have both tissues
        both_mask = (df_confound["n_tumor_batches"] > 0) & (df_confound["n_normal_batches"] > 0)
        n_paired = int(both_mask.sum())
        n_same_paired = int(df_confound.loc[both_mask, "t_n_same_batch"].sum())

        print(f"  Patients with both T+N in cE01: {n_paired}")
        print(f"    T and N share >= 1 batch: {n_same_paired}")
        print(f"    T and N fully in different batches: {n_paired - n_same_paired}")

        tissue_batch_confound = "STRONG" if n_same_paired < n_paired * 0.1 else \
                                "MODERATE" if n_same_paired < n_paired * 0.5 else "WEAK"
        print(f"  Tissue-batch confound: {tissue_batch_confound}")
        results["checks"]["batch_tissue_confound"] = {
            "n_batches": n_batches,
            "n_paired_patients": n_paired,
            "n_same_batch": n_same_paired,
            "n_diff_batch": n_paired - n_same_paired,
            "confound_level": tissue_batch_confound,
        }
    else:
        print("  No batch column found; batch confound cannot be assessed")
        results["checks"]["batch_tissue_confound"] = {"confound_level": "UNKNOWN"}

    # ── 7. QC column discovery (Amendment A8) ─────────────────────────────
    print("\n--- 7. QC column discovery ---")

    qc_candidates = {
        "total_counts": ["total_counts", "nCount_RNA", "n_counts"],
        "n_genes": ["n_genes_by_counts", "nFeature_RNA", "n_genes"],
        "pct_mt": ["pct_counts_mt", "percent.mt", "percent_mito"],
        "S_score": ["S_score", "S.Score", "phase_S"],
        "G2M_score": ["G2M_score", "G2M.Score", "phase_G2M"],
    }
    qc_columns_found = {}
    for qc_name, candidates in qc_candidates.items():
        found = None
        for c in candidates:
            if c in obs.columns:
                found = c
                break
        qc_columns_found[qc_name] = found
        status = f"found as '{found}'" if found else "NOT FOUND"
        print(f"  {qc_name}: {status}")

    results["checks"]["qc_columns"] = qc_columns_found

    # ── 8. Obs columns inventory ─────────────────────────────────────────
    print("\n--- 8. Full obs columns inventory ---")
    results["checks"]["obs_columns"] = list(obs.columns)
    print(f"  {len(obs.columns)} columns: {list(obs.columns)}")

    # ── Summary ──────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("GATE 1 SUMMARY")
    print("=" * 72)

    if results["kill"]:
        results["status"] = "KILL"
        print(f"  VERDICT: KILL — {results['kill_reason']}")
    else:
        results["status"] = "PASS"
        print(f"  VERDICT: PASS")

    print(f"  Patient column: {pid_col}")
    print(f"  Unique patients: {n_patients}")
    print(f"  cE01 == Stem/TA-like: {results['checks'].get('ce01_stem_exact_match', '?')}")
    print(f"  .X: {x_info['verdict']}")
    print(f"  .raw.X: {results['checks'].get('norm_raw_X', {}).get('verdict', 'N/A')}")
    print(f"  Counts .X: {counts_x_info['verdict']}")
    print(f"  Target genes: {len(genes_present)}/21 present, {len(results['checks'].get('target_genes_in_duplicates', []))} in duplicates")
    print(f"  Eligible at >=30: {eligibility['thresh_30']['n_eligible']}")
    print(f"  Eligible at >=20: {eligibility['thresh_20']['n_eligible']}")
    print(f"  Eligible at >=50: {eligibility['thresh_50']['n_eligible']}")
    print(f"  Power status: {power_status} (effective threshold: {effective_threshold})")
    print(f"  Batch column: {batch_col}")

    _save_and_exit(results, norm_path, counts_path, cfg)


def _save_and_exit(results, norm_path, counts_path, cfg):
    # Convert numpy types for JSON
    def _convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, set):
            return list(obj)
        if isinstance(obj, dict):
            return {k: _convert(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_convert(v) for v in obj]
        return obj

    results = _convert(results)

    # Save results
    results_path = RESULTS_DIR / "gate1_input_validation.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  Saved: {results_path}")

    # Save provenance
    prov = {
        "script": "01_gate1_input_validation.py",
        "gate": "Gate 1",
        "inputs": {
            "epithelial_normalized": {"path": str(norm_path)},
            "epithelial_counts": {"path": str(counts_path)},
            "config": {"path": str(CFG_PATH)},
        },
        "outputs": {
            "results": str(results_path),
        },
        "status": results["status"],
    }
    prov_path = PROV_DIR / "gate1_provenance.json"
    with open(prov_path, "w") as f:
        json.dump(prov, f, indent=2, default=str)
    print(f"  Saved: {prov_path}")


if __name__ == "__main__":
    main()
