"""
Lee 2020 External Normal-Anchor Consistency Check
===================================================
Protocol: protocol_frozen_v10 (A009)
NOT a gate. Descriptive sign-only cross-cohort consistency of the
normal ISC-FAO anticorrelation.
"""

import gzip, hashlib, json, warnings
from pathlib import Path
from datetime import date

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scipy.stats as stats

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CV2  = Path(__file__).resolve().parent
RESULTS = CV2 / "results"
TABLES  = RESULTS / "tables"
PROV    = CV2 / "provenance"
TABLES.mkdir(parents=True, exist_ok=True)

# ── Gene sets (LOCKED) ───────────────────────────────────────────────
ISC_9     = ["LGR5","OLFM4","ASCL2","SOX9","CD44","SMOC2","RGMB","EPHB2","MSI1"]
ISC_8     = [g for g in ISC_9 if g != "LGR5"]
FAO_12    = ["CPT1A","HMGCS2","FABP1","ACOX1","ACADL","ACADM","PDK4","PPARD","PPARA","ANGPTL4","HADH","ACAA2"]
FAO_11    = [g for g in FAO_12 if g != "FABP1"]
CORE_BO_5 = ["CPT1A","ACADL","ACADM","HADH","ACAA2"]
UNION_21  = set(ISC_9) | set(FAO_12)

# ── Copied verbatim from 03_gate2b_clean_background.py lines 62-111 ─
def recover_scanpy_controls(X, var_names, gene_list, ctrl_as_ref=True,
                            ctrl_size=50, n_bins=25, random_state=42):
    gene_list_idx = pd.Index(gene_list)
    gene_pool = pd.Index(var_names, dtype="string")
    if sp.issparse(X):
        obs_avg_vals = np.asarray(X.mean(axis=0)).ravel()
    else:
        obs_avg_vals = X.mean(axis=0)
    obs_avg = pd.Series(obs_avg_vals, index=gene_pool)
    obs_avg = obs_avg[np.isfinite(obs_avg)]
    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    keep_ctrl_in_obs_cut = False if ctrl_as_ref else obs_cut.index.isin(gene_list_idx)
    np.random.seed(random_state)
    all_controls = pd.Index([], dtype="string")
    for cut in np.unique(obs_cut.loc[gene_list_idx.intersection(obs_avg.index)]):
        r_genes = obs_cut[(obs_cut == cut) & ~keep_ctrl_in_obs_cut].index
        if ctrl_size < len(r_genes):
            r_genes = r_genes.to_series().sample(ctrl_size).index
        if ctrl_as_ref:
            r_genes = r_genes.difference(gene_list_idx)
        all_controls = all_controls.union(r_genes)
    return sorted(all_controls.tolist())

# ── Copied verbatim from 02_gate2_recompute_scores.py lines 223-249 ─
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

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def sha256_str(s):
    return hashlib.sha256(s.encode()).hexdigest()

def get_versions():
    import importlib.metadata as im
    return {k: im.version(k) for k in ["scanpy","numpy","pandas","scipy"]}

# ═════════════════════════════════════════════════════════════════════
def load_lee_anchor(ann_path, mat_path, cohort_name, expected_n):
    """STEP 1: load annotation, stream matrix for anchor cells only."""
    print(f"\n{'='*60}")
    print(f"  {cohort_name}: Loading")
    print(f"{'='*60}")

    # 1a. Annotation
    ann = pd.read_csv(ann_path, sep="\t", compression="gzip")
    anchor_mask = (ann["Class"] == "Normal") & (ann["Cell_subtype"] == "Stem-like/TA")
    anchor_ann = ann[anchor_mask].copy()
    n_anchor = len(anchor_ann)
    print(f"  Anchor cells (Normal Stem-like/TA): {n_anchor}")
    if n_anchor != expected_n:
        raise RuntimeError(f"STOP: expected {expected_n} anchor cells, got {n_anchor}")
    anchor_barcodes = anchor_ann["Index"].tolist()
    anchor_set = set(anchor_barcodes)

    # 1b. Stream matrix — keep only anchor columns
    print(f"  Streaming {Path(mat_path).name}...")
    with gzip.open(mat_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        mat_barcodes = header[1:]
        # Verify order matches annotation
        assert mat_barcodes == ann["Index"].tolist(), "STOP: barcode order mismatch"
        # Column indices for anchor cells
        anchor_col_idx = [i for i, bc in enumerate(mat_barcodes) if bc in anchor_set]
        assert len(anchor_col_idx) == n_anchor

        genes = []
        data_rows = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            gene = parts[0]
            genes.append(gene)
            row_vals = np.array([float(parts[1 + i]) for i in anchor_col_idx], dtype=np.float32)
            data_rows.append(row_vals)

    # genes x cells -> cells x genes
    raw_mat = np.stack(data_rows, axis=0).T  # (n_anchor, n_genes)
    ordered_anchor_barcodes = [mat_barcodes[i] for i in anchor_col_idx]
    print(f"  Matrix: {raw_mat.shape[0]} cells x {raw_mat.shape[1]} genes")

    # 1c. Verify uniqueness + target presence
    gene_counts = pd.Series(genes).value_counts()
    dups = gene_counts[gene_counts > 1]
    if len(dups) > 0:
        raise RuntimeError(f"STOP: duplicate gene symbols: {dups.index.tolist()[:10]}")
    gene_set = set(genes)
    for g in UNION_21:
        assert g in gene_set, f"STOP: target gene {g} missing"
    for g in ISC_8 + FAO_11 + CORE_BO_5:
        assert g in gene_set, f"STOP: variant gene {g} missing"
    assert len(set(anchor_barcodes)) == n_anchor, "STOP: duplicate anchor barcodes"
    assert len(set(ordered_anchor_barcodes)) == n_anchor
    print(f"  Uniqueness + targets: OK")

    # 1d. Build AnnData
    # NOTE: A009 I9-4 specifies normalize_total + log1p assuming raw integer counts.
    # However, the GEO matrices (GSE132465/GSE144735) are ALREADY natural-log-TPM
    # (confirmed: nonzero values are fractional ~0.19-5.6, not integers).
    # Applying normalize_total + log1p would DOUBLE-normalize.
    # The data is ALREADY in a log-normalized space comparable to Pelka's log1p(CPM/1e4).
    # We therefore load the data AS-IS and document this discrepancy.
    adata = ad.AnnData(
        X=raw_mat,
        obs=pd.DataFrame({"barcode": ordered_anchor_barcodes}).set_index("barcode"),
        var=pd.DataFrame(index=genes),
    )
    bc_to_pid = dict(zip(anchor_ann["Index"], anchor_ann["Patient"]))
    adata.obs["PID"] = [bc_to_pid[bc] for bc in ordered_anchor_barcodes]

    # Verify data is already log-normalized (not raw counts)
    nz = raw_mat[raw_mat != 0]
    is_integer = np.all(nz == np.round(nz))
    print(f"  Matrix nonzero: min={nz.min():.4f}, max={nz.max():.4f}, integer-like={is_integer}")
    assert not is_integer, "UNEXPECTED: data looks like raw counts, not log-TPM"
    print(f"  Data is already natural-log-TPM; normalize_total+log1p SKIPPED (would double-normalize)")
    print(f"  DISCREPANCY: A009 I9-4 assumed raw counts; actual GEO files are log-TPM.")

    return adata, genes, anchor_ann


def select_controls_and_score(adata, genes, cohort_name):
    """STEP 2+3: select controls, compute 5 scores."""
    var_names = list(adata.var_names)
    X = adata.X

    # STEP 2: Control selection on anchor subset
    isc_pool = [g for g in var_names if g not in set(FAO_12)]
    fao_pool = [g for g in var_names if g not in set(ISC_9)]

    # Subset X to the pools for control selection
    isc_pool_idx = [var_names.index(g) for g in isc_pool]
    fao_pool_idx = [var_names.index(g) for g in fao_pool]
    X_isc_pool = X[:, isc_pool_idx]
    X_fao_pool = X[:, fao_pool_idx]

    isc_controls = recover_scanpy_controls(X_isc_pool, isc_pool, ISC_9)
    fao_controls = recover_scanpy_controls(X_fao_pool, fao_pool, FAO_12)

    print(f"  {cohort_name} ISC controls: {len(isc_controls)}, FAO controls: {len(fao_controls)}")

    # GATING GUARD
    isc_overlap = set(isc_controls) & UNION_21
    fao_overlap = set(fao_controls) & UNION_21
    if isc_overlap or fao_overlap:
        raise RuntimeError(f"STOP: control-target overlap ISC={isc_overlap} FAO={fao_overlap}")
    print(f"  Guard: 0 target overlap in both. PASS")

    # STEP 3: Compute 5 scores
    scores = {}
    scores["ISC_9"]     = fixed_bg_score(X, var_names, ISC_9, isc_controls)
    scores["ISC_8"]     = fixed_bg_score(X, var_names, ISC_8, isc_controls)
    scores["FAO_12"]    = fixed_bg_score(X, var_names, FAO_12, fao_controls)
    scores["FAO_11"]    = fixed_bg_score(X, var_names, FAO_11, fao_controls)
    scores["core_bo_5"] = fixed_bg_score(X, var_names, CORE_BO_5, fao_controls)

    ctrl_info = {
        "isc_n": len(isc_controls), "fao_n": len(fao_controls),
        "isc_hash": sha256_str(",".join(sorted(isc_controls))),
        "fao_hash": sha256_str(",".join(sorted(fao_controls))),
    }
    return scores, ctrl_info


def compute_correlations(scores):
    """STEP 4: pooled Pearson r for 4 pairs."""
    pairs = [
        ("ISC_9",  "FAO_12",    "primary"),
        ("ISC_8",  "FAO_12",    "ISC-8 (LGR5 dropped)"),
        ("ISC_9",  "FAO_11",    "FABP1 dropped"),
        ("ISC_9",  "core_bo_5", "core beta-ox-5"),
    ]
    rows = []
    for isc_key, fao_key, label in pairs:
        r, _ = stats.pearsonr(scores[isc_key], scores[fao_key])
        rows.append({"pair": label, "isc": isc_key, "fao": fao_key,
                      "r": round(r, 6), "n": len(scores[isc_key])})
    return rows


def per_patient_sign_tally(adata, scores, min_cells=30):
    """STEP 5: per-patient r(ISC-9, FAO-12) for patients with >= min_cells."""
    pids = adata.obs["PID"].values
    unique_pids = sorted(set(pids))
    rows = []
    for pid in unique_pids:
        mask = pids == pid
        n = int(mask.sum())
        if n >= min_cells:
            r, _ = stats.pearsonr(scores["ISC_9"][mask], scores["FAO_12"][mask])
            rows.append({"PID": pid, "n": n, "r": round(r, 6), "sign": "neg" if r < 0 else "pos"})
        else:
            rows.append({"PID": pid, "n": n, "r": None, "sign": f"<{min_cells} cells"})
    return rows


# ═════════════════════════════════════════════════════════════════════
def main():
    print("=" * 70)
    print("Lee 2020 External Normal-Anchor Consistency (A009)")
    print("=" * 70)
    versions = get_versions()
    print(f"Versions: {versions}")

    lee_dir = ROOT / "data" / "external" / "lee"
    cohorts = [
        ("SMC",  lee_dir / "GSE132465_annotation.txt.gz", lee_dir / "GSE132465_matrix.txt.gz", 406),
        ("KUL3", lee_dir / "GSE144735_annotation.txt.gz", lee_dir / "GSE144735_matrix.txt.gz", 495),
    ]

    all_corrs = []
    all_sign_tally = []
    all_percell = []
    ctrl_infos = {}
    feasibility = {}

    for cohort_name, ann_path, mat_path, expected_n in cohorts:
        # STEP 1
        adata, genes, anchor_ann = load_lee_anchor(ann_path, mat_path, cohort_name, expected_n)

        # STEP 2+3
        scores, ctrl_info = select_controls_and_score(adata, genes, cohort_name)
        ctrl_infos[cohort_name] = ctrl_info

        # Save per-cell scores
        pc_df = adata.obs[["PID"]].copy()
        pc_df["cohort"] = cohort_name
        for k, v in scores.items():
            pc_df[k] = v
        all_percell.append(pc_df)

        # STEP 4
        corr_rows = compute_correlations(scores)
        for row in corr_rows:
            row["cohort"] = cohort_name
        all_corrs.extend(corr_rows)

        # STEP 5
        sign_rows = per_patient_sign_tally(adata, scores)
        for row in sign_rows:
            row["cohort"] = cohort_name
        all_sign_tally.extend(sign_rows)

        # STEP 7 feasibility
        full_ann = pd.read_csv(ann_path, sep="\t", compression="gzip")
        epi = full_ann[full_ann["Cell_type"] == "Epithelial cells"]
        xt = pd.crosstab(epi["Cell_subtype"], epi["Class"], margins=True)
        feasibility[cohort_name] = xt
        stem_tumor = len(epi[(epi["Cell_subtype"] == "Stem-like/TA") & (epi["Class"] == "Tumor")])
        cms_normal = len(epi[epi["Cell_subtype"].str.startswith("CMS") & (epi["Class"] == "Normal")])
        print(f"\n  {cohort_name} feasibility: Stem-like/TA in Tumor={stem_tumor}, CMS in Normal={cms_normal}")

    # ── STEP 6: Pelka reference anchors ──────────────────────────────
    print(f"\n{'='*60}")
    print("  Pelka Reference Anchors")
    print(f"{'='*60}")

    # 6a. Committed Gate 2b clean-bg
    g2b = pd.read_parquet(RESULTS / "gate2b_clean_background_scores_cE01.parquet")
    pelka_norm = g2b[g2b["SPECIMEN_TYPE"] == "N"]
    n_pelka_norm = len(pelka_norm)
    print(f"  Pelka normal cE01: {n_pelka_norm}")
    assert n_pelka_norm == 11786, f"Expected 11786, got {n_pelka_norm}"

    r_primary_pelka, _ = stats.pearsonr(
        pelka_norm["ISC_scanpy_cleanbg"].values,
        pelka_norm["FAO_scanpy_cleanbg"].values)
    print(f"  Committed clean-bg r_normal: {r_primary_pelka:.4f}")
    assert abs(r_primary_pelka - (-0.3822)) < 5e-4, \
        f"STOP: reproduction guard failed: {r_primary_pelka:.6f} vs -0.3822"
    print(f"  Reproduction guard 6a: PASS ({r_primary_pelka:.6f})")

    # 6b. Recompute companions from matrix + committed controls
    pelka_h5ad = ROOT / "data" / "processed" / "pelka_epithelial.h5ad"
    adata_pelka = ad.read_h5ad(pelka_h5ad, backed="r")
    cE01_barcodes = g2b.index.tolist()
    X_cE01 = adata_pelka[cE01_barcodes].X
    if sp.issparse(X_cE01):
        X_cE01 = X_cE01.tocsc()
    pelka_var = list(adata_pelka.var_names)

    # Load committed Pelka controls
    with open(PROV / "clean_control_gene_sets.json") as f:
        pelka_ctrl = json.load(f)
    pelka_isc_ctrl = pelka_ctrl["recovered_selected_control_genes"]["isc_controls"]
    pelka_fao_ctrl = pelka_ctrl["recovered_selected_control_genes"]["fao_controls"]
    assert len(pelka_isc_ctrl) == 300
    assert len(pelka_fao_ctrl) == 400

    # Reproduction guard 6b
    isc9_pelka = fixed_bg_score(X_cE01, pelka_var, ISC_9, pelka_isc_ctrl)
    fao12_pelka = fixed_bg_score(X_cE01, pelka_var, FAO_12, pelka_fao_ctrl)
    ref_isc = g2b.loc[cE01_barcodes, "ISC_scanpy_cleanbg"].values
    ref_fao = g2b.loc[cE01_barcodes, "FAO_scanpy_cleanbg"].values
    isc_diff = float(np.max(np.abs(isc9_pelka - ref_isc)))
    fao_diff = float(np.max(np.abs(fao12_pelka - ref_fao)))
    print(f"  Guard 6b: ISC max|diff|={isc_diff:.2e}, FAO max|diff|={fao_diff:.2e}")
    assert isc_diff < 1e-9 and fao_diff < 1e-9, \
        f"STOP: guard 6b failed (ISC={isc_diff}, FAO={fao_diff})"
    print(f"  Reproduction guard 6b: PASS")

    # Pelka normal subset
    norm_mask = np.array(g2b["SPECIMEN_TYPE"].values == "N")
    pelka_scores = {
        "ISC_9":     isc9_pelka[norm_mask],
        "ISC_8":     fixed_bg_score(X_cE01, pelka_var, ISC_8, pelka_isc_ctrl)[norm_mask],
        "FAO_12":    fao12_pelka[norm_mask],
        "FAO_11":    fixed_bg_score(X_cE01, pelka_var, FAO_11, pelka_fao_ctrl)[norm_mask],
        "core_bo_5": fixed_bg_score(X_cE01, pelka_var, CORE_BO_5, pelka_fao_ctrl)[norm_mask],
    }
    pelka_corrs = compute_correlations(pelka_scores)
    for row in pelka_corrs:
        row["cohort"] = "Pelka"
    all_corrs.extend(pelka_corrs)

    # ── Save outputs ──────────────────────────────────────────────────
    corr_df = pd.DataFrame(all_corrs)
    corr_df.to_csv(TABLES / "lee_normal_anchor_correlations.csv", index=False)
    print(f"\nSaved: tables/lee_normal_anchor_correlations.csv")

    sign_df = pd.DataFrame(all_sign_tally)
    sign_df.to_csv(TABLES / "lee_per_patient_sign_tally.csv", index=False)
    print(f"Saved: tables/lee_per_patient_sign_tally.csv")

    percell_df = pd.concat(all_percell, ignore_index=False)
    percell_df.to_parquet(TABLES / "lee_anchor_percell_scores.parquet")
    print(f"Saved: tables/lee_anchor_percell_scores.parquet ({len(percell_df)} cells)")

    # ── Terminal summary ──────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("3-COHORT NORMAL-ANCHOR CORRELATION TABLE")
    print(f"{'='*70}")
    pivot = corr_df.pivot_table(index="pair", columns="cohort", values="r")
    for col in ["Pelka", "SMC", "KUL3"]:
        if col not in pivot.columns:
            pivot[col] = np.nan
    pivot = pivot[["Pelka", "SMC", "KUL3"]]
    print(pivot.to_string())

    primary_row = corr_df[corr_df["pair"] == "primary"]
    all_neg = all(primary_row["r"].values < 0)
    print(f"\nPrimary r(ISC-9, FAO-12) negative in all 3 cohorts: {all_neg}")

    # Per-patient sign tally summary
    eligible_signs = sign_df[sign_df["r"].notna()]
    n_neg = (eligible_signs["sign"] == "neg").sum()
    n_total = len(eligible_signs)
    print(f"Per-patient sign tally (>=30 cells): {n_neg}/{n_total} negative")

    # ── Write report ──────────────────────────────────────────────────
    lines = [
        "# Lee 2020 External Normal-Anchor Consistency Report", "",
        f"**Protocol:** protocol_frozen_v10 (A009)",
        f"**Versions:** {versions}",
        f"**random_state:** 42", "",
        "## Reproduction Guards", "",
        f"- **6a (committed Pelka clean-bg r_normal):** {r_primary_pelka:.6f} (expected -0.3822). **PASS**",
        f"- **6b (fixed_bg_score reproduces parquet):** ISC max|diff|={isc_diff:.2e}, FAO max|diff|={fao_diff:.2e}. **PASS**",
    ]
    for cn, ci in ctrl_infos.items():
        lines.append(f"- **{cn} controls:** ISC={ci['isc_n']}, FAO={ci['fao_n']}; "
                      f"target overlap=0. **PASS**")

    lines.extend(["", "## 3-Cohort Normal-Anchor Correlations", "",
        "| Pair | Pelka (n=11786) | SMC (n=406) | KUL3 (n=495) |",
        "|------|----------------|-------------|--------------|"])
    for pair_label in ["primary", "ISC-8 (LGR5 dropped)", "FABP1 dropped", "core beta-ox-5"]:
        vals = {}
        for _, row in corr_df[corr_df["pair"] == pair_label].iterrows():
            vals[row["cohort"]] = f"{row['r']:.4f}"
        lines.append(f"| {pair_label} | {vals.get('Pelka','N/A')} | "
                      f"{vals.get('SMC','N/A')} | {vals.get('KUL3','N/A')} |")

    lines.extend(["",
        f"**Primary r(ISC-9, FAO-12) negative in all 3 cohorts: {all_neg}**", "",
        "## Descriptive Consistency Verdict", "",
        "The inverse ISC-lipid/FAO coupling observed in Pelka normal cE01 is "
        + ("reproduced in sign" if all_neg else "NOT reproduced in sign")
        + " in both Lee normal Stem-like/TA cohorts (SMC and KUL3). "
        "This is a normal-baseline consistency anchor and does NOT externally validate "
        "the tumor-side decoupling, which remains single-cohort and tissue-batch-confounded.", "",
        "## Companion Scores (descriptive)", "",
        "- r(ISC-9, FAO-11) tests whether FABP1 drives the NORMAL anticorrelation. "
        "This is a DIFFERENT quantity from Gate 3's tumor-vs-normal shift.",
        "- r(ISC-9, core-bo-5) tests the core beta-oxidation sub-score in normal tissue.", ""])

    # Per-patient sign tally
    lines.extend(["## Per-Patient Sign Tally (>=30 cells)", "",
        "| Cohort | PID | n | r | sign |", "|--------|-----|---|---|------|"])
    for _, row in sign_df.iterrows():
        r_str = f"{row['r']:.4f}" if row["r"] is not None else "—"
        lines.append(f"| {row['cohort']} | {row['PID']} | {row['n']} | {r_str} | {row['sign']} |")

    # Feasibility
    lines.extend(["", "## Feasibility Cross-Tabs (I9-3)", ""])
    for cn, xt in feasibility.items():
        lines.append(f"### {cn}")
        lines.append("```")
        lines.append(xt.to_string())
        lines.append("```")
        epi_ann = pd.read_csv(
            lee_dir / f"{'GSE132465' if cn=='SMC' else 'GSE144735'}_annotation.txt.gz",
            sep="\t", compression="gzip")
        epi_ann = epi_ann[epi_ann["Cell_type"] == "Epithelial cells"]
        st = len(epi_ann[(epi_ann["Cell_subtype"]=="Stem-like/TA") & (epi_ann["Class"]=="Tumor")])
        cn_n = len(epi_ann[epi_ann["Cell_subtype"].str.startswith("CMS") & (epi_ann["Class"]=="Normal")])
        lines.append(f"Stem-like/TA in Tumor: **{st}**. CMS in Normal: **{cn_n}**.")
        lines.append("")

    # Claim ceiling
    lines.extend(["## Claim Ceiling (A009 I9-8)", "",
        "- The Lee normal anchor does NOT satisfy the E11 external-replication requirement "
        "for the central finding.",
        "- Even if the normal anchor is fully sign-consistent, the donor-aware TUMOR-vs-NORMAL "
        "shift remains a single-cohort (Pelka) result with NO external replication.",
        "- Lee carries no independent sequencing-batch column beyond Sample, so tissue is "
        "plausibly sample-confounded in Lee as well.",
        "- Strongest permitted statement: 'the inverse ISC-lipid/FAO coupling observed in "
        "Pelka normal cE01 is [/is not] reproduced in sign in two independent normal colon "
        "Stem-like/TA datasets; this is a normal-baseline consistency anchor and does not "
        "externally validate the tumor-side decoupling, which remains single-cohort and "
        "tissue-batch-confounded.'"])

    with open(RESULTS / "lee_external_consistency_report.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nSaved: results/lee_external_consistency_report.md")

    # ── Provenance ────────────────────────────────────────────────────
    prov = {
        "script": "08_lee_external_consistency.py",
        "protocol": "protocol_frozen_v10 (A009)",
        "versions": versions,
        "random_state": 42,
        "inputs": {
            "smc_ann": str(cohorts[0][1]), "smc_mat": str(cohorts[0][2]),
            "kul3_ann": str(cohorts[1][1]), "kul3_mat": str(cohorts[1][2]),
            "pelka_h5ad": str(pelka_h5ad),
            "gate2b_parquet": str(RESULTS / "gate2b_clean_background_scores_cE01.parquet"),
            "clean_controls": str(PROV / "clean_control_gene_sets.json"),
        },
        "input_hashes": {
            "smc_ann": sha256_file(cohorts[0][1]),
            "smc_mat": sha256_file(cohorts[0][2]),
            "kul3_ann": sha256_file(cohorts[1][1]),
            "kul3_mat": sha256_file(cohorts[1][2]),
            "clean_controls": sha256_file(PROV / "clean_control_gene_sets.json"),
        },
        "anchor_counts": {"SMC": 406, "KUL3": 495, "Pelka_normal_cE01": 11786},
        "lee_controls": ctrl_infos,
        "pelka_controls": {"isc_n": 300, "fao_n": 400},
        "reproduction_guards": {
            "guard_6a_r_normal": round(r_primary_pelka, 6),
            "guard_6b_isc_max_diff": isc_diff,
            "guard_6b_fao_max_diff": fao_diff,
        },
        "primary_all_negative": all_neg,
        "status": "COMPLETE",
    }
    with open(PROV / "lee_external_consistency.json", "w") as f:
        json.dump(prov, f, indent=2)
    print(f"Saved: provenance/lee_external_consistency.json")

    print("\nDone.")


if __name__ == "__main__":
    main()
