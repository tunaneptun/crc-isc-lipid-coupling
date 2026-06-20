#!/usr/bin/env python3
"""Independent reconstruction for the final integrity audit.

This script intentionally does not import or execute core_validation_v2 gate
scripts. It reads source matrices and committed result tables, recomputes the
load-bearing summaries, and writes audit-only outputs.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "core_validation_v2" / "audit"


def dense_col(x, idx):
    col = x[:, idx]
    if sparse.issparse(col):
        return np.asarray(col.toarray()).ravel()
    return np.asarray(col).ravel()


def dense_rows_cols(x, rows, cols):
    block = x[rows][:, cols]
    if sparse.issparse(block):
        return block.toarray()
    return np.asarray(block)


def matrix_summary(adata, name, sample_cells=2000, sample_genes=1000):
    rng = np.random.default_rng(123)
    n_obs, n_vars = adata.shape
    rows = rng.choice(n_obs, min(sample_cells, n_obs), replace=False)
    cols = rng.choice(n_vars, min(sample_genes, n_vars), replace=False)
    block = dense_rows_cols(adata.X, rows, cols)
    nz = block[block != 0]
    frac_int = float(np.mean(np.isclose(nz, np.round(nz)))) if nz.size else np.nan
    lib = np.asarray(adata.X[rows].sum(axis=1)).ravel() if sparse.issparse(adata.X) else adata.X[rows].sum(axis=1)
    row_block = adata.X[rows]
    if sparse.issparse(row_block):
        rb = row_block.copy()
        rb.data = np.expm1(rb.data)
        expm1_lib = np.asarray(rb.sum(axis=1)).ravel()
    else:
        expm1_lib = np.expm1(row_block).sum(axis=1)
    return {
        "matrix": name,
        "shape": [int(n_obs), int(n_vars)],
        "dtype": str(adata.X.dtype),
        "sparse": bool(sparse.issparse(adata.X)),
        "sample_n_nonzero": int(nz.size),
        "sample_min_nonzero": float(np.min(nz)) if nz.size else None,
        "sample_median_nonzero": float(np.median(nz)) if nz.size else None,
        "sample_max": float(np.max(block)) if block.size else None,
        "sample_frac_integer_nonzero": frac_int,
        "cell_sum_min_median_max": [float(np.min(lib)), float(np.median(lib)), float(np.max(lib))],
        "expm1_cell_sum_min_median_max": [
            float(np.min(expm1_lib)),
            float(np.median(expm1_lib)),
            float(np.max(expm1_lib)),
        ],
    }


def pearson(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4 or np.var(x[ok]) == 0 or np.var(y[ok]) == 0:
        return np.nan
    return float(np.corrcoef(x[ok], y[ok])[0, 1])


def bootstrap_ci(vals, seed=42, k=10000):
    rng = np.random.default_rng(seed)
    vals = np.asarray(vals, float)
    means = np.empty(k)
    n = vals.size
    for i in range(k):
        means[i] = vals[rng.integers(0, n, n)].mean()
    return [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))]


def donor_table(scores, isc_col, fao_col, pids):
    rows = []
    for pid in pids:
        rec = {"PID": pid}
        for tissue in ["N", "T"]:
            sub = scores[(scores["PID"] == pid) & (scores["SPECIMEN_TYPE"] == tissue)]
            r = pearson(sub[isc_col], sub[fao_col])
            rec[f"n_{'normal' if tissue == 'N' else 'tumor'}"] = int(len(sub))
            rec[f"r_{'normal' if tissue == 'N' else 'tumor'}"] = r
            rec[f"z_{'normal' if tissue == 'N' else 'tumor'}"] = float(np.arctanh(np.clip(r, -0.999999999999, 0.999999999999)))
        rec["delta_z"] = rec["z_tumor"] - rec["z_normal"]
        rows.append(rec)
    return pd.DataFrame(rows)


def summarize_donor(df):
    dz = df["delta_z"].to_numpy(float)
    stat, p = wilcoxon(dz, alternative="two-sided", zero_method="wilcox", method="auto")
    loo = [np.delete(dz, i).mean() for i in range(len(dz))]
    return {
        "n": int(len(dz)),
        "mean_delta_z": float(dz.mean()),
        "median_delta_z": float(np.median(dz)),
        "wilcoxon_stat": float(stat),
        "wilcoxon_p": float(p),
        "bootstrap_ci": bootstrap_ci(dz),
        "loo_min": float(np.min(loo)),
        "loo_max": float(np.max(loo)),
        "loo_all_positive": bool(np.all(np.asarray(loo) > 0)),
        "prop_positive": float(np.mean(dz > 0)),
        "mean_r_normal": float(df["r_normal"].mean()),
        "mean_r_tumor": float(df["r_tumor"].mean()),
    }


def fixed_score(adata, targets, controls):
    var = pd.Index(adata.var_names)
    ti = [var.get_loc(g) for g in targets]
    ci = [var.get_loc(g) for g in controls]
    tm = np.asarray(adata.X[:, ti].mean(axis=1)).ravel()
    cm = np.asarray(adata.X[:, ci].mean(axis=1)).ravel()
    return tm - cm


def cramer_v(table):
    obs = np.asarray(table, float)
    n = obs.sum()
    expected = obs.sum(axis=1, keepdims=True) @ obs.sum(axis=0, keepdims=True) / n
    chi2 = ((obs - expected) ** 2 / expected)[expected > 0].sum()
    r, c = obs.shape
    return float(np.sqrt((chi2 / n) / min(c - 1, r - 1)))


def read_lee_matrix(path, ncols=None):
    if ncols is None:
        df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    else:
        df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip", usecols=range(ncols + 1))
    return df


def lee_scale_summary(path, name):
    df = read_lee_matrix(path, ncols=500)
    vals = df.to_numpy()
    nz = vals[vals != 0]
    sums = np.expm1(vals[:, : min(vals.shape[1], 500)]).sum(axis=0)
    return {
        "matrix": name,
        "shape": list(df.shape),
        "min20_nonzero": [float(x) for x in np.sort(nz)[:20]],
        "max": float(vals.max()),
        "fractional_nonzero": bool(np.any(~np.isclose(nz, np.round(nz)))),
        "sample_expm1_colsum_min_median_max": [float(sums.min()), float(np.median(sums)), float(sums.max())],
    }


def main():
    cfg = yaml.safe_load((ROOT / "core_validation_v2/config.yaml").read_text())
    isc9 = cfg["gene_sets"]["isc_9"]
    fao12 = cfg["gene_sets"]["ppar_fao_12"]
    core5 = cfg["gene_sets"]["core_beta_oxidation_5"]
    norm = ad.read_h5ad(ROOT / cfg["inputs"]["epithelial_normalized"])
    counts = ad.read_h5ad(ROOT / cfg["inputs"]["epithelial_counts"])

    summaries = []
    summaries.append(matrix_summary(norm, "Pelka normalized .X"))
    raw_adata = ad.AnnData(X=norm.raw.X, obs=norm.obs.copy(), var=norm.raw.var.copy()) if norm.raw is not None else None
    if raw_adata is not None:
        summaries.append(matrix_summary(raw_adata, "Pelka normalized .raw.X"))
    summaries.append(matrix_summary(counts, "Pelka counts .X"))
    summaries.append(lee_scale_summary(ROOT / "data/external/lee/GSE132465_matrix.txt.gz", "Lee SMC local matrix"))
    summaries.append(lee_scale_summary(ROOT / "data/external/lee/GSE144735_matrix.txt.gz", "Lee KUL3 local matrix"))
    (OUT / "matrix_scale_reconstruction.json").write_text(json.dumps(summaries, indent=2))

    obs = norm.obs.copy()
    cE01 = obs["cl295v11SubShort"].astype(str).eq("cE01")
    stem = obs["epithelial_subtype"].astype(str).eq("Stem/TA-like")
    pop = {
        "n_cells": int(norm.n_obs),
        "n_genes": int(norm.n_vars),
        "cE01": int(cE01.sum()),
        "Stem_TA_like": int(stem.sum()),
        "exact_match": bool((cE01 == stem).all()),
        "cE01_tissue_counts": obs.loc[cE01, "SPECIMEN_TYPE"].value_counts().to_dict(),
    }
    counts_by_pid = obs.loc[cE01].groupby(["PID", "SPECIMEN_TYPE"]).size().unstack(fill_value=0)
    for thr in [20, 30, 50]:
        elig = counts_by_pid[(counts_by_pid.get("T", 0) >= thr) & (counts_by_pid.get("N", 0) >= thr)]
        pop[f"eligible_ge_{thr}"] = int(len(elig))
        pop[f"eligible_ge_{thr}_pids"] = list(map(str, elig.index))
    primary_pids = pop["eligible_ge_30_pids"]
    bt = obs.loc[cE01 & obs["PID"].isin(primary_pids), ["SPECIMEN_TYPE", "batchID"]]
    tissue_batch = pd.crosstab(bt["SPECIMEN_TYPE"], bt["batchID"])
    shared = 0
    for pid in primary_pids:
        p = obs.loc[cE01 & (obs["PID"] == pid)]
        if set(p.loc[p["SPECIMEN_TYPE"] == "T", "batchID"]) & set(p.loc[p["SPECIMEN_TYPE"] == "N", "batchID"]):
            shared += 1
    pop["primary_shared_tn_batch_pids"] = shared
    pop["primary_tissue_batch_cramers_v"] = cramer_v(tissue_batch)
    pd.Series(pop).to_json(OUT / "population_reconstruction.json", indent=2)
    counts_by_pid.to_csv(OUT / "cell_counts_by_pid_tissue.csv")

    scores = pd.read_parquet(ROOT / "core_validation_v2/results/gate2_recomputed_scores_cE01.parquet").reset_index()
    clean = pd.read_parquet(ROOT / "core_validation_v2/results/gate2b_clean_background_scores_cE01.parquet").reset_index()
    scores = scores.merge(clean[["barcode", "ISC_scanpy_cleanbg", "FAO_scanpy_cleanbg"]], on="barcode", how="left")
    pooled = []
    for label, ic, fc in [
        ("scanpy", "ISC_scanpy", "FAO_scanpy"),
        ("clean_bg", "ISC_scanpy_cleanbg", "FAO_scanpy_cleanbg"),
        ("aucell", "ISC_aucell", "FAO_aucell"),
        ("zscore", "ISC_zscore", "FAO_zscore"),
    ]:
        rn = pearson(scores.loc[scores.SPECIMEN_TYPE == "N", ic], scores.loc[scores.SPECIMEN_TYPE == "N", fc])
        rt = pearson(scores.loc[scores.SPECIMEN_TYPE == "T", ic], scores.loc[scores.SPECIMEN_TYPE == "T", fc])
        pooled.append({"method": label, "r_normal": rn, "r_tumor": rt, "r_diff": rt - rn})
    pd.DataFrame(pooled).to_csv(OUT / "pooled_reconstruction.csv", index=False)

    gate4 = donor_table(scores, "ISC_scanpy", "FAO_scanpy", primary_pids)
    gate4.to_csv(OUT / "gate4_per_patient_reconstruction.csv", index=False)
    gate4_summary = summarize_donor(gate4)
    (OUT / "gate4_summary_reconstruction.json").write_text(json.dumps(gate4_summary, indent=2))

    g4_committed = pd.read_csv(ROOT / "core_validation_v2/results/tables/gate4_per_patient.csv")
    if "method" in g4_committed.columns:
        g4_committed = g4_committed[g4_committed["method"] == "scanpy"].copy()
    g4_committed = g4_committed.rename(columns={"Delta_z": "delta_z"})
    diff = gate4.merge(g4_committed, on="PID", suffixes=("_audit", "_committed"))
    rows = []
    for col in ["n_normal", "n_tumor", "r_normal", "r_tumor", "z_normal", "z_tumor", "delta_z"]:
        rows.append({"column": col, "max_abs_diff": float(np.nanmax(np.abs(diff[f"{col}_audit"] - diff[f"{col}_committed"])))})
    pd.DataFrame(rows).to_csv(OUT / "gate4_per_patient_diff.csv", index=False)

    gate5 = pd.read_csv(ROOT / "core_validation_v2/results/tables/gate5_draw_summary.csv")
    g5 = {
        "n_iterations": int(len(gate5)),
        "direction_retention_fraction": float((gate5["mean_delta_z_scanpy"] > 0).mean()),
        "median_downsampled_magnitude": float(gate5["mean_delta_z_scanpy"].median()),
        "min_contributing_patients": int(gate5["n_patients_contributing_scanpy"].min()),
    }
    (OUT / "gate5_summary_from_draws.json").write_text(json.dumps(g5, indent=2))

    controls = json.loads((ROOT / "core_validation_v2/provenance/clean_control_gene_sets.json").read_text())["recovered_selected_control_genes"]
    cidx = np.where(cE01.to_numpy())[0]
    cdata = norm[cidx].copy()
    cobs = cdata.obs.reset_index(names="barcode")
    isc_controls = controls["isc_controls"]
    fao_controls = controls["fao_controls"]
    base_isc = fixed_score(cdata, isc9, isc_controls)
    base_fao = fixed_score(cdata, fao12, fao_controls)
    base_scores = cobs[["barcode", "PID", "SPECIMEN_TYPE"]].copy()
    base_scores["base_isc"] = base_isc
    base_scores["base_fao"] = base_fao
    base = donor_table(base_scores, "base_isc", "base_fao", primary_pids)
    base_mean = base["delta_z"].mean()
    variants = [
        ("baseline", isc9, fao12),
        ("isc_loo_LGR5", [g for g in isc9 if g != "LGR5"], fao12),
        ("isc_loo_OLFM4", [g for g in isc9 if g != "OLFM4"], fao12),
        ("hmgcs2_removed", isc9, [g for g in fao12 if g != "HMGCS2"]),
        ("fabp1_removed", isc9, [g for g in fao12 if g != "FABP1"]),
        ("angptl4_removed", isc9, [g for g in fao12 if g != "ANGPTL4"]),
        ("core_beta_oxidation_5", isc9, core5),
    ]
    vrows = []
    for name, ig, fg in variants:
        tmp = cobs[["barcode", "PID", "SPECIMEN_TYPE"]].copy()
        tmp["isc"] = fixed_score(cdata, ig, isc_controls)
        tmp["fao"] = fixed_score(cdata, fg, fao_controls)
        dt = donor_table(tmp, "isc", "fao", primary_pids)
        mean = float(dt["delta_z"].mean())
        vrows.append({"variant": name, "mean_delta_z": mean, "retention_ratio": abs(mean) / abs(base_mean)})
    pd.DataFrame(vrows).to_csv(OUT / "gate3_key_variant_reconstruction.csv", index=False)

    var = pd.Index(counts.var_names)
    lgr5 = dense_col(counts.X, var.get_loc("LGR5"))
    lgr5_c = lgr5[cidx]
    lgr5_scores = base_scores.copy()
    lgr5_scores["isc8"] = fixed_score(cdata, [g for g in isc9 if g != "LGR5"], isc_controls)
    lgr5_scores["fao"] = base_fao
    lgr5_scores["lgr5_detected"] = lgr5_c > 0
    elig_l = (
        lgr5_scores[lgr5_scores.lgr5_detected]
        .groupby(["PID", "SPECIMEN_TYPE"])
        .size()
        .unstack(fill_value=0)
    )
    elig_l = elig_l[(elig_l.get("T", 0) >= 15) & (elig_l.get("N", 0) >= 15)]
    ldt = donor_table(lgr5_scores[lgr5_scores.lgr5_detected], "isc8", "fao", list(elig_l.index))
    same6 = donor_table(lgr5_scores, "isc8", "fao", list(elig_l.index))
    g6 = {
        "lgr5_detected_cE01": int(lgr5_scores.lgr5_detected.sum()),
        "eligible_pids": list(map(str, elig_l.index)),
        "n_eligible": int(len(elig_l)),
        "subset_mean_delta_z": float(ldt["delta_z"].mean()),
        "same6_full_mean_delta_z": float(same6["delta_z"].mean()),
    }
    (OUT / "gate6_reconstruction.json").write_text(json.dumps(g6, indent=2))

    committed = {
        "gate2_pooled": pd.read_csv(ROOT / "core_validation_v2/results/tables/gate2_pooled_decoupling.csv").to_dict("records"),
        "gate2b_pooled": pd.read_csv(ROOT / "core_validation_v2/results/tables/gate2b_pooled_cleanbg.csv").to_dict("records"),
        "gate4_primary_summary": pd.read_csv(ROOT / "core_validation_v2/results/tables/gate4_primary_summary.csv").to_dict("records"),
        "gate5_primary_summary": pd.read_csv(ROOT / "core_validation_v2/results/tables/gate5_primary_summary.csv").to_dict("records"),
        "gate6_per_patient_rows": int(len(pd.read_csv(ROOT / "core_validation_v2/results/tables/gate6_per_patient.csv"))),
    }
    (OUT / "committed_artifact_snapshot.json").write_text(json.dumps(committed, indent=2, default=str))


if __name__ == "__main__":
    main()
