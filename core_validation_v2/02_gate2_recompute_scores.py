"""
Gate 2 — Score Recomputation
==============================
Protocol: CORE_VALIDATION_PROTOCOL.md + AMENDMENT_001 + AMENDMENT_002
         (tag protocol_frozen_v3)
Execution order: 2 (after Gate 1)

Inputs
------
- data/processed/pelka_epithelial.h5ad (read-only; log-normalized .X)
- core_validation_v2/config.yaml

Outputs
-------
- core_validation_v2/results/gate2_recomputed_scores_cE01.parquet
- core_validation_v2/results/gate2_score_recomputation_report.md
- core_validation_v2/results/tables/gate2_pooled_decoupling.csv
- core_validation_v2/provenance/fixed_control_gene_sets.json
- core_validation_v2/provenance/02_gate2_recompute_scores.json

What this gate does (B1-B5)
---------------------------
1. Recompute ISC & FAO scores with three methods on log-normalized .X.
2. Build and save a fixed-background scorer for Gate 3 decomposition.
3. Descriptive pooled cell-level decoupling (reporting only, not a kill).
4. Write canonical cE01 score table for Gates 4-6.

Kill criteria
-------------
NOT a project-kill gate. If pooled direction is not consistent across
scoring methods, note the discrepancy but PROCEED to Gate 4.
"""

import hashlib
import json
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scipy.stats as stats
import yaml

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG_PATH = Path(__file__).resolve().parent / "config.yaml"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
TABLES_DIR = RESULTS_DIR / "tables"
PROV_DIR = Path(__file__).resolve().parent / "provenance"

TABLES_DIR.mkdir(parents=True, exist_ok=True)


def sha256_str(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


# ═══════════════════════════════════════════════════════════════════════
# AUCell — deterministic, chunked implementation (B3)
# ═══════════════════════════════════════════════════════════════════════
def aucell_score(X_sparse, var_names, gene_set, auc_threshold=0.05, seed=42):
    """
    Deterministic AUCell scoring, vectorized within chunks.
    X_sparse: sparse CSR matrix (cells x genes), log-normalized.
    var_names: array of gene names aligned with columns.
    gene_set: list of gene names.
    auc_threshold: fraction of ranked genes used for AUC.
    Returns: 1-D array of AUC scores, one per cell.
    """
    var_idx = {g: i for i, g in enumerate(var_names)}
    target_indices = np.array([var_idx[g] for g in gene_set if g in var_idx])
    n_targets = len(target_indices)
    n_genes = X_sparse.shape[1]
    n_top = max(1, int(np.ceil(n_genes * auc_threshold)))
    max_auc = n_targets * n_top  # maximum possible AUC

    if max_auc == 0:
        return np.zeros(X_sparse.shape[0])

    # Create a boolean mask for target genes (for fast lookup)
    is_target = np.zeros(n_genes, dtype=bool)
    is_target[target_indices] = True

    rng = np.random.RandomState(seed)
    scores = np.empty(X_sparse.shape[0], dtype=np.float64)

    chunk_size = 500
    n_cells = X_sparse.shape[0]
    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        chunk = X_sparse[start:end]
        if sp.issparse(chunk):
            chunk = chunk.toarray()
        n_chunk = chunk.shape[0]

        # Generate tie-breaking noise for the whole chunk
        tie_break = rng.random((n_chunk, n_genes))

        # Descending sort: negate expression, use tie_break as secondary
        # np.lexsort sorts by last key first, so we want (-expression, tie_break)
        # But lexsort works column-wise. We need row-wise sorting.
        # Use argsort on a combined key: -expression * large_factor + tie_break
        # This ensures descending expression with random tie-breaking.
        sort_key = -chunk.astype(np.float64) * 1e10 + tie_break
        # argsort along axis=1 gives the ranking (ascending of sort_key = descending expression)
        order = np.argsort(sort_key, axis=1, kind='stable')

        # For each cell, check which of the top n_top ranked genes are targets
        top_genes = order[:, :n_top]  # shape (n_chunk, n_top)
        # Check target membership
        target_hits = is_target[top_genes]  # shape (n_chunk, n_top), bool
        # AUC weight: position n_top - rank_pos for each hit
        weights = np.arange(n_top, 0, -1)  # [n_top, n_top-1, ..., 1]
        auc_values = (target_hits * weights).sum(axis=1) / max_auc
        scores[start:end] = auc_values

        if (start // chunk_size) % 50 == 0:
            print(f"    AUCell progress: {end:,}/{n_cells:,} cells", flush=True)

    return scores


# ═══════════════════════════════════════════════════════════════════════
# Mean z-score (B4)
# ═══════════════════════════════════════════════════════════════════════
def zscore_score(X_sparse, var_names, gene_set, ddof=0):
    """
    Mean z-score: standardize each target gene across all cells (ddof=0),
    average per cell.
    Returns: 1-D array, list of genes used, list of skipped genes.
    """
    var_idx = {g: i for i, g in enumerate(var_names)}
    z_cols = []
    used_genes = []
    skipped_genes = []
    for g in gene_set:
        if g not in var_idx:
            continue
        col = X_sparse[:, var_idx[g]]
        if sp.issparse(col):
            col = col.toarray().ravel()
        else:
            col = np.asarray(col).ravel()
        col = col.astype(np.float64)
        mu = col.mean()
        sigma = col.std(ddof=ddof)
        if sigma == 0:
            skipped_genes.append(g)
            continue
        z_cols.append((col - mu) / sigma)
        used_genes.append(g)

    if len(z_cols) == 0:
        return np.zeros(X_sparse.shape[0]), used_genes, skipped_genes

    z_matrix = np.column_stack(z_cols)
    return z_matrix.mean(axis=1), used_genes, skipped_genes


# ═══════════════════════════════════════════════════════════════════════
# Fixed-background scorer (B2)
# ═══════════════════════════════════════════════════════════════════════
def build_fixed_control_genes(X_sparse, var_names, gene_set, n_bins=25,
                              ctrl_size=50, seed=42):
    """
    Deterministically select expression-matched control genes for a gene set.
    Mimics scanpy's binning logic: bin genes by mean expression, then for each
    target gene, sample ctrl_size controls from the same bin (excluding targets).
    Returns the unique sorted control gene list and per-gene controls.
    """
    rng = np.random.RandomState(seed)
    var_idx = {g: i for i, g in enumerate(var_names)}

    # Compute mean expression per gene (sparse-safe)
    if sp.issparse(X_sparse):
        gene_means = np.asarray(X_sparse.mean(axis=0)).ravel()
    else:
        gene_means = X_sparse.mean(axis=0)

    # Bin genes by mean expression
    n_genes = len(var_names)
    sorted_indices = np.argsort(gene_means)
    bin_assignments = np.zeros(n_genes, dtype=int)
    bin_size = int(np.ceil(n_genes / n_bins))
    for i, idx in enumerate(sorted_indices):
        bin_assignments[idx] = i // bin_size

    target_set = set()
    for g in gene_set:
        if g in var_idx:
            target_set.add(var_idx[g])

    # For each target gene, sample ctrl_size controls from the same bin
    all_controls = set()
    per_gene_controls = {}
    for g in gene_set:
        if g not in var_idx:
            continue
        g_idx = var_idx[g]
        g_bin = bin_assignments[g_idx]
        # All genes in the same bin, excluding targets
        candidates = [j for j in range(n_genes)
                      if bin_assignments[j] == g_bin and j not in target_set]
        if len(candidates) == 0:
            continue
        if len(candidates) <= ctrl_size:
            chosen = candidates
        else:
            chosen = rng.choice(candidates, size=ctrl_size, replace=False).tolist()
        per_gene_controls[g] = [var_names[c] for c in sorted(chosen)]
        all_controls.update(chosen)

    control_genes_unique = sorted([var_names[c] for c in all_controls])
    return control_genes_unique, per_gene_controls


def fixed_bg_score(X_sparse, var_names, gene_set, control_genes):
    """
    Score = mean(target genes) - mean(control genes) per cell.
    """
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
    print("=" * 72)
    print("Gate 2 — Score Recomputation")
    print("=" * 72)

    with open(CFG_PATH) as f:
        cfg = yaml.safe_load(f)

    # ── Load data ──────────────────────────────────────────────────────
    norm_path = ROOT / cfg["inputs"]["epithelial_normalized"]
    print(f"\nLoading {norm_path.name} (read-only) ...")
    adata = ad.read_h5ad(str(norm_path))
    print(f"  Shape: {adata.n_obs:,} x {adata.n_vars:,}")

    scanpy_version = None
    try:
        import importlib.metadata
        scanpy_version = importlib.metadata.version("scanpy")
    except Exception:
        scanpy_version = getattr(sc, "__version__", "unknown")

    var_names = np.array(adata.var_names)
    gene_pool_hash = sha256_str(",".join(var_names))

    isc_genes = list(cfg["gene_sets"]["isc_9"])
    fao_genes = list(cfg["gene_sets"]["ppar_fao_12"])

    pid_col = cfg["population"]["patient_column"]
    tissue_col = cfg["population"]["tissue_column"]
    short_col = cfg["population"]["subtype_column_short"]

    report_lines = []

    def rpt(line=""):
        report_lines.append(line)
        print(line)

    rpt(f"# Gate 2 — Score Recomputation Report\n")
    rpt(f"**Protocol:** CORE_VALIDATION_PROTOCOL.md + AMENDMENT_001 + AMENDMENT_002 (tag protocol_frozen_v3)")
    rpt(f"**scanpy version:** {scanpy_version}")
    rpt(f"**Cells:** {adata.n_obs:,}")
    rpt(f"**Genes:** {adata.n_vars:,}")
    rpt(f"**gene_pool hash (sha256 of comma-joined var_names):** `{gene_pool_hash[:16]}...`")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # B1. Official scanpy score_genes
    # ══════════════════════════════════════════════════════════════════
    rpt("## B1. Official scanpy score_genes (PRIMARY)")
    rpt()

    resolved_call = (
        "sc.tl.score_genes(adata, gene_list=GENES, "
        "ctrl_as_ref=True, ctrl_size=50, gene_pool=None, "
        "n_bins=25, score_name=NAME, random_state=42, "
        "copy=False, use_raw=False, layer=None)"
    )
    rpt(f"Resolved call: `{resolved_call}`")
    rpt(f"Note: gene_pool=None means full var_names; ctrl_as_ref=True is the installed default (scanpy {scanpy_version}).")
    rpt()

    sc.tl.score_genes(adata, gene_list=isc_genes,
                      ctrl_as_ref=True, ctrl_size=50, gene_pool=None,
                      n_bins=25, score_name="ISC_scanpy", random_state=42,
                      copy=False, use_raw=False, layer=None)

    sc.tl.score_genes(adata, gene_list=fao_genes,
                      ctrl_as_ref=True, ctrl_size=50, gene_pool=None,
                      n_bins=25, score_name="FAO_scanpy", random_state=42,
                      copy=False, use_raw=False, layer=None)

    rpt(f"ISC_scanpy: mean={adata.obs['ISC_scanpy'].mean():.6f}, std={adata.obs['ISC_scanpy'].std():.6f}")
    rpt(f"FAO_scanpy: mean={adata.obs['FAO_scanpy'].mean():.6f}, std={adata.obs['FAO_scanpy'].std():.6f}")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # B2. Fixed-background scorer
    # ══════════════════════════════════════════════════════════════════
    rpt("## B2. Fixed-background scorer")
    rpt()

    X = adata.X  # sparse, log-normalized

    isc_ctrl, isc_per_gene = build_fixed_control_genes(
        X, var_names, isc_genes, n_bins=25, ctrl_size=50, seed=42)
    fao_ctrl, fao_per_gene = build_fixed_control_genes(
        X, var_names, fao_genes, n_bins=25, ctrl_size=50, seed=42)

    rpt(f"ISC fixed controls: {len(isc_ctrl)} unique genes")
    rpt(f"FAO fixed controls: {len(fao_ctrl)} unique genes")

    # Compute fixed-bg scores
    isc_fixedbg = fixed_bg_score(X, var_names, isc_genes, isc_ctrl)
    fao_fixedbg = fixed_bg_score(X, var_names, fao_genes, fao_ctrl)
    adata.obs["ISC_scanpy_fixedbg"] = isc_fixedbg
    adata.obs["FAO_scanpy_fixedbg"] = fao_fixedbg

    # Validate fixed-bg vs official
    r_isc_fb = stats.pearsonr(adata.obs["ISC_scanpy"].values, isc_fixedbg)[0]
    r_fao_fb = stats.pearsonr(adata.obs["FAO_scanpy"].values, fao_fixedbg)[0]
    rpt(f"Fixed-bg vs official agreement (all epithelial):")
    rpt(f"  ISC Pearson r = {r_isc_fb:.6f}")
    rpt(f"  FAO Pearson r = {r_fao_fb:.6f}")
    rpt()

    # Save control gene lists
    fixed_bg_prov = {
        "description": "Fixed control-gene background for Gate 3 decomposition (Amendment B2)",
        "gene_universe_n": len(var_names),
        "gene_universe_hash": gene_pool_hash,
        "duplicate_handling": "none (0 duplicates in normalized var_names)",
        "n_bins": 25,
        "ctrl_size": 50,
        "seed": 42,
        "isc_target_genes": isc_genes,
        "isc_control_genes_unique": isc_ctrl,
        "isc_per_gene_controls": isc_per_gene,
        "fao_target_genes": fao_genes,
        "fao_control_genes_unique": fao_ctrl,
        "fao_per_gene_controls": fao_per_gene,
    }
    with open(PROV_DIR / "fixed_control_gene_sets.json", "w") as f:
        json.dump(fixed_bg_prov, f, indent=2)
    rpt(f"Saved: provenance/fixed_control_gene_sets.json")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # B3. AUCell
    # ══════════════════════════════════════════════════════════════════
    rpt("## B3. AUCell (deterministic, chunked)")
    rpt()

    X_csr = X.tocsr() if sp.issparse(X) else X

    isc_auc = aucell_score(X_csr, var_names, isc_genes, auc_threshold=0.05, seed=42)
    fao_auc = aucell_score(X_csr, var_names, fao_genes, auc_threshold=0.05, seed=42)
    adata.obs["ISC_aucell"] = isc_auc
    adata.obs["FAO_aucell"] = fao_auc

    rpt(f"ISC_aucell: mean={isc_auc.mean():.6f}, std={isc_auc.std():.6f}")
    rpt(f"FAO_aucell: mean={fao_auc.mean():.6f}, std={fao_auc.std():.6f}")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # B4. Mean z-score
    # ══════════════════════════════════════════════════════════════════
    rpt("## B4. Mean z-score (ddof=0)")
    rpt()

    isc_z, isc_z_used, isc_z_skip = zscore_score(X_csr, var_names, isc_genes, ddof=0)
    fao_z, fao_z_used, fao_z_skip = zscore_score(X_csr, var_names, fao_genes, ddof=0)
    adata.obs["ISC_zscore"] = isc_z
    adata.obs["FAO_zscore"] = fao_z

    rpt(f"ISC_zscore: mean={isc_z.mean():.6f}, std={isc_z.std():.6f}, genes_used={len(isc_z_used)}/9")
    if isc_z_skip:
        rpt(f"  SKIPPED (zero variance): {isc_z_skip}")
    rpt(f"FAO_zscore: mean={fao_z.mean():.6f}, std={fao_z.std():.6f}, genes_used={len(fao_z_used)}/12")
    if fao_z_skip:
        rpt(f"  SKIPPED (zero variance): {fao_z_skip}")
    rpt()

    obs = adata.obs
    ce01_mask = obs[short_col].astype(str) == "cE01"

    # ══════════════════════════════════════════════════════════════════
    # Descriptive pooled decoupling
    # ══════════════════════════════════════════════════════════════════
    rpt("## Descriptive Pooled Decoupling (cE01, REPORTING ONLY)")
    rpt()
    rpt("Population: cE01 only (Stem/TA-like_prolif excluded per protocol).")
    rpt()

    ce01_obs = obs[ce01_mask].copy()
    tumor_mask = ce01_obs[tissue_col] == cfg["population"]["tumor_code"]
    normal_mask = ce01_obs[tissue_col] == cfg["population"]["normal_code"]

    method_pairs = [
        ("ISC_scanpy", "FAO_scanpy", "scanpy (recomputed)"),
        ("ISC_aucell", "FAO_aucell", "AUCell (recomputed)"),
        ("ISC_zscore", "FAO_zscore", "z-score (recomputed)"),
    ]

    decoupling_rows = []
    for isc_col, fao_col, method_label in method_pairs:
        if isc_col not in ce01_obs.columns or fao_col not in ce01_obs.columns:
            continue

        isc_t = ce01_obs.loc[tumor_mask, isc_col].values.astype(np.float64)
        fao_t = ce01_obs.loc[tumor_mask, fao_col].values.astype(np.float64)
        isc_n = ce01_obs.loc[normal_mask, isc_col].values.astype(np.float64)
        fao_n = ce01_obs.loc[normal_mask, fao_col].values.astype(np.float64)

        valid_t = np.isfinite(isc_t) & np.isfinite(fao_t)
        valid_n = np.isfinite(isc_n) & np.isfinite(fao_n)

        r_tumor = stats.pearsonr(isc_t[valid_t], fao_t[valid_t])[0]
        r_normal = stats.pearsonr(isc_n[valid_n], fao_n[valid_n])[0]
        diff = r_tumor - r_normal

        row = {
            "method": method_label,
            "n_tumor": int(valid_t.sum()),
            "n_normal": int(valid_n.sum()),
            "r_tumor": round(r_tumor, 6),
            "r_normal": round(r_normal, 6),
            "r_diff": round(diff, 6),
        }
        decoupling_rows.append(row)
        rpt(f"  {method_label}: r_normal={r_normal:.4f}, r_tumor={r_tumor:.4f}, diff={diff:.4f}")

    rpt()

    dec_df = pd.DataFrame(decoupling_rows)
    dec_df.to_csv(TABLES_DIR / "gate2_pooled_decoupling.csv", index=False)
    rpt(f"Saved: tables/gate2_pooled_decoupling.csv")
    rpt()

    # Check direction consistency
    recomputed_diffs = [r for r in decoupling_rows if "recomputed" in r["method"]]
    directions = [np.sign(r["r_diff"]) for r in recomputed_diffs]
    direction_consistent = len(set(directions)) == 1
    rpt(f"Direction consistency across recomputed methods: {'YES' if direction_consistent else 'NO'}")
    if recomputed_diffs:
        common_dir = ("tumor more negative" if directions[0] < 0
                      else "tumor more positive" if directions[0] > 0
                      else "zero")
        rpt(f"Common direction: r_diff = r_tumor - r_normal is {common_dir}")

    rpt()

    # ══════════════════════════════════════════════════════════════════
    # Write canonical cE01 score table (parquet)
    # ══════════════════════════════════════════════════════════════════
    rpt("## Canonical cE01 Score Table")
    rpt()

    score_cols = [
        "ISC_scanpy", "ISC_aucell", "ISC_zscore",
        "FAO_scanpy", "FAO_aucell", "FAO_zscore",
        "ISC_scanpy_fixedbg", "FAO_scanpy_fixedbg",
    ]
    meta_cols = [pid_col, tissue_col]
    qc_cols_cfg = cfg.get("qc_columns", {})
    qc_cols = [v for v in [
        qc_cols_cfg.get("total_counts"),
        qc_cols_cfg.get("n_genes_by_counts"),
        qc_cols_cfg.get("pct_counts_mt"),
    ] if v is not None and v in obs.columns]

    out_cols = meta_cols + score_cols + qc_cols
    ce01_df = ce01_obs[out_cols].copy()
    ce01_df.index.name = "barcode"

    parquet_path = RESULTS_DIR / "gate2_recomputed_scores_cE01.parquet"
    ce01_df.to_parquet(parquet_path, index=True)
    rpt(f"Saved: gate2_recomputed_scores_cE01.parquet ({len(ce01_df):,} cells, {len(out_cols)} columns)")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # Write report
    # ══════════════════════════════════════════════════════════════════
    report_path = RESULTS_DIR / "gate2_score_recomputation_report.md"
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))
    rpt(f"\nSaved: {report_path.name}")

    # ══════════════════════════════════════════════════════════════════
    # Provenance
    # ══════════════════════════════════════════════════════════════════
    prov = {
        "script": "02_gate2_recompute_scores.py",
        "gate": "Gate 2",
        "scanpy_version": scanpy_version,
        "resolved_scanpy_call": resolved_call,
        "gene_pool_hash": gene_pool_hash,
        "ctrl_as_ref": True,
        "aucell_implementation": "custom deterministic chunked (auc_threshold=0.05, seed=42)",
        "zscore_ddof": 0,
        "isc_genes": isc_genes,
        "fao_genes": fao_genes,
        "isc_zscore_genes_used": isc_z_used,
        "isc_zscore_genes_skipped": isc_z_skip,
        "fao_zscore_genes_used": fao_z_used,
        "fao_zscore_genes_skipped": fao_z_skip,
        "fixed_bg_isc_n_controls": len(isc_ctrl),
        "fixed_bg_fao_n_controls": len(fao_ctrl),
        "fixed_bg_vs_official_isc_r": round(r_isc_fb, 6),
        "fixed_bg_vs_official_fao_r": round(r_fao_fb, 6),
        "n_epithelial_cells": int(adata.n_obs),
        "n_ce01_cells": int(ce01_mask.sum()),
        "inputs": {
            "epithelial_normalized": str(norm_path),
            "config": str(CFG_PATH),
        },
        "outputs": {
            "scores_parquet": str(parquet_path),
            "report": str(report_path),
            "pooled_decoupling": str(TABLES_DIR / "gate2_pooled_decoupling.csv"),
            "fixed_control_genes": str(PROV_DIR / "fixed_control_gene_sets.json"),
        },
        "status": "COMPLETE",
    }
    prov_path = PROV_DIR / "02_gate2_recompute_scores.json"
    with open(prov_path, "w") as f:
        json.dump(prov, f, indent=2, default=str)
    rpt(f"Saved: {prov_path.name}")

    print("\n" + "=" * 72)
    print("GATE 2 COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()
