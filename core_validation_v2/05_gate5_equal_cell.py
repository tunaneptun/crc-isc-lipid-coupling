"""
Gate 5 — Equal-Cell Downsampling (SECOND HARD KILL)
=====================================================
Protocol: protocol_frozen_v7 (A006 pins implementation of frozen criterion)
Decision: direction_retention >= 0.90 AND median_magnitude >= 0.2003145
"""

import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats
import yaml

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG_PATH = Path(__file__).resolve().parent / "config.yaml"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
TABLES_DIR = RESULTS_DIR / "tables"
PROV_DIR = Path(__file__).resolve().parent / "provenance"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

K = 1000
SAMPLE_SIZE = 30
FULL_DATA_EFFECT = 0.400629
MAGNITUDE_THRESHOLD = 0.50 * FULL_DATA_EFFECT  # 0.2003145
DIRECTION_THRESHOLD = 0.90
BOOT_N = 10000
BOOT_SEED = 42


def get_versions():
    import importlib.metadata as im
    return {k: im.version(k) for k in ["numpy", "scipy", "pandas"]}


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


def donor_aware_test_on_dz(dz_arr, seed=BOOT_SEED, n_boot=BOOT_N):
    n = len(dz_arr)
    mean_dz = float(np.mean(dz_arr))
    median_dz = float(np.median(dz_arr))
    n_zeros = int(np.sum(dz_arr == 0))
    try:
        wstat, wp = stats.wilcoxon(dz_arr, alternative="two-sided",
                                   zero_method="wilcox", method="auto")
        wstat, wp = float(wstat), float(wp)
    except ValueError:
        wstat, wp = np.nan, np.nan
    rng = np.random.RandomState(seed)
    boot_means = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        boot_means[i] = np.mean(dz_arr[idx])
    ci_lo = float(np.percentile(boot_means, 2.5))
    ci_hi = float(np.percentile(boot_means, 97.5))
    ci_excludes_0 = (ci_lo > 0) or (ci_hi < 0)
    # LOO
    loo_rows = []
    for i in range(n):
        loo_dz = np.delete(dz_arr, i)
        loo_mean = float(np.mean(loo_dz))
        try:
            _, loo_p = stats.wilcoxon(loo_dz, alternative="two-sided",
                                      zero_method="wilcox", method="auto")
            loo_p = float(loo_p)
        except ValueError:
            loo_p = np.nan
        rel_change = abs(loo_mean - mean_dz) / abs(mean_dz) * 100 if abs(mean_dz) > 0 else 0.0
        loo_rows.append({
            "n_remaining": n - 1,
            "mean_delta_z": round(loo_mean, 6),
            "wilcoxon_p": round(loo_p, 6) if np.isfinite(loo_p) else None,
            "sign_positive": loo_mean > 0,
            "relative_effect_change_pct": round(rel_change, 2),
        })
    all_same_sign = all(r["sign_positive"] == (mean_dz > 0) for r in loo_rows)
    max_rel = max(r["relative_effect_change_pct"] for r in loo_rows)
    return {
        "n": n, "mean_delta_z": mean_dz, "median_delta_z": median_dz,
        "wilcoxon_stat": wstat, "wilcoxon_p": wp,
        "n_exact_zeros": n_zeros, "effective_n": n - n_zeros,
        "bootstrap_ci_lo": ci_lo, "bootstrap_ci_hi": ci_hi,
        "ci_excludes_0": ci_excludes_0,
        "loo_all_same_sign": all_same_sign,
        "loo_max_rel_change_pct": max_rel,
        "loo_effect_sensitive": max_rel > 25.0,
        "prop_positive": float(np.mean(dz_arr > 0)),
        "loo_rows": loo_rows,
    }


def main():
    print("=" * 72, flush=True)
    print("Gate 5 — Equal-Cell Downsampling (SECOND HARD KILL)", flush=True)
    print("=" * 72, flush=True)

    versions = get_versions()
    for k, v in versions.items():
        print(f"  {k}: {v}", flush=True)

    with open(CFG_PATH) as f:
        cfg = yaml.safe_load(f)

    pid_col = cfg["population"]["patient_column"]
    tissue_col = cfg["population"]["tissue_column"]
    tumor_code = cfg["population"]["tumor_code"]
    normal_code = cfg["population"]["normal_code"]

    # ── Load scores ────────────────────────────────────────────────────
    g2 = pd.read_parquet(RESULTS_DIR / "gate2_recomputed_scores_cE01.parquet")
    g2b = pd.read_parquet(RESULTS_DIR / "gate2b_clean_background_scores_cE01.parquet")
    assert g2.index.equals(g2b.index), "Barcode index mismatch"
    df = g2.join(g2b[["ISC_scanpy_cleanbg", "FAO_scanpy_cleanbg"]])
    assert len(df) == 61953, f"Expected 61953, got {len(df)}"
    print(f"\nJoined: {len(df):,} cE01 cells", flush=True)

    # ── Eligibility + drift check ─────────────────────────────────────
    thresh = cfg["eligibility"]["primary_min_cells"]
    pid_counts = df.groupby([pid_col, tissue_col]).size().unstack(fill_value=0)
    eligible_mask = (pid_counts.get(tumor_code, pd.Series(dtype=int)) >= thresh) & \
                    (pid_counts.get(normal_code, pd.Series(dtype=int)) >= thresh)
    eligible_pids = sorted(eligible_mask[eligible_mask].index.tolist())

    with open(RESULTS_DIR / "gate1_input_validation.json") as f:
        g1 = json.load(f)
    g1_pids = sorted([p["patient"] for p in g1["checks"]["eligibility"]["thresh_30"]["patients"]])
    drift_match = eligible_pids == g1_pids
    print(f"Eligible: {len(eligible_pids)}, Gate 1: {len(g1_pids)}, Match: {drift_match}", flush=True)

    if not drift_match:
        print("STOP: NEEDS_HUMAN_CONFIRMATION — eligibility drift", flush=True)
        return

    elig_df = df[df[pid_col].isin(eligible_pids)]

    report = []
    def rpt(line=""):
        report.append(line)
        print(line, flush=True)

    rpt("# Gate 5 — Equal-Cell Downsampling Report\n")
    rpt(f"**Protocol:** protocol_frozen_v7 (A006)")
    rpt(f"**Versions:** {versions}")
    rpt(f"**K:** {K}, **Sample size:** {SAMPLE_SIZE}/PID/tissue")
    rpt(f"**Seed plan:** SeedSequence(42).spawn({K})")
    rpt(f"**Full-data effect:** {FULL_DATA_EFFECT}")
    rpt(f"**Magnitude threshold (50%):** {MAGNITUDE_THRESHOLD}")
    rpt(f"**Direction threshold:** {DIRECTION_THRESHOLD}")
    rpt(f"**Eligible patients:** {len(eligible_pids)}, drift check: PASSED")
    rpt()

    # ── Seed plan ─────────────────────────────────────────────────────
    parent_seq = np.random.SeedSequence(42)
    child_seqs = parent_seq.spawn(K)

    # ── Pre-index cells by PID x tissue ───────────────────────────────
    methods = [
        ("scanpy", "ISC_scanpy", "FAO_scanpy"),
        ("cleanbg", "ISC_scanpy_cleanbg", "FAO_scanpy_cleanbg"),
        ("aucell", "ISC_aucell", "FAO_aucell"),
        ("zscore", "ISC_zscore", "FAO_zscore"),
    ]

    # Frozen order: PIDs lexicographic, tissues ["N", "T"]
    tissue_order = [normal_code, tumor_code]
    strata_order = [(pid, t) for pid in eligible_pids for t in tissue_order]

    # Pre-sort barcode pools per stratum
    barcode_pools = {}
    for pid, tissue in strata_order:
        mask = (elig_df[pid_col] == pid) & (elig_df[tissue_col] == tissue)
        barcodes = sorted(elig_df.index[mask].tolist())
        barcode_pools[(pid, tissue)] = barcodes

    # ══════════════════════════════════════════════════════════════════
    # TASK 1 — K iterations
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 1: Downsampling Iterations\n")

    draw_rows = []
    # Accumulators for consolidated per-patient
    # per_patient_dz[method][pid] = list of Delta_z across iterations
    per_patient_dz = {m: {pid: [] for pid in eligible_pids} for m, _, _ in methods}
    per_patient_zn = {pid: [] for pid in eligible_pids}
    per_patient_zt = {pid: [] for pid in eligible_pids}

    valid_iteration_guard_fail = False

    for draw_id in range(K):
        child = child_seqs[draw_id]
        rng = np.random.default_rng(child)

        # Sample cells in frozen deterministic order
        sampled_barcodes_all = []
        sampled_per_stratum = {}
        for pid, tissue in strata_order:
            pool = barcode_pools[(pid, tissue)]
            chosen_idx = rng.choice(len(pool), size=SAMPLE_SIZE, replace=False)
            chosen = [pool[i] for i in sorted(chosen_idx)]
            sampled_per_stratum[(pid, tissue)] = chosen
            for bc in chosen:
                sampled_barcodes_all.append(f"{pid}\t{tissue}\t{bc}")

        # Barcode hash
        sampled_barcodes_all.sort()
        hash_input = "\n".join(sampled_barcodes_all)
        bc_hash = hashlib.sha256(hash_input.encode()).hexdigest()

        # Child state for provenance (first/last only to keep file manageable)
        child_spawn_key = str(child.spawn_key)
        child_state = json.dumps({"entropy": child.entropy, "spawn_key": list(child.spawn_key),
                                  "n_children_spawned": child.n_children_spawned})

        # Per-method donor-aware effect
        draw_rec = {"draw_id": draw_id, "child_spawn_key": child_spawn_key,
                    "child_state_json": child_state, "barcode_hash": bc_hash}

        for method_name, isc_col, fao_col in methods:
            patient_dz = []
            n_contributing = 0

            for pid in eligible_pids:
                results_per_tissue = {}
                for tissue in tissue_order:
                    barcodes = sampled_per_stratum[(pid, tissue)]
                    sub = elig_df.loc[barcodes]
                    isc = sub[isc_col].values.astype(np.float64)
                    fao = sub[fao_col].values.astype(np.float64)
                    r, status, _ = safe_pearsonr(isc, fao)
                    results_per_tissue[tissue] = (r, status)

                r_n, st_n = results_per_tissue[normal_code]
                r_t, st_t = results_per_tissue[tumor_code]

                if st_n == "ok" and st_t == "ok":
                    rn_c, _ = clip_r(r_n)
                    rt_c, _ = clip_r(r_t)
                    zn = fisher_z(rn_c)
                    zt = fisher_z(rt_c)
                    dz = zt - zn
                    patient_dz.append(dz)
                    n_contributing += 1
                    per_patient_dz[method_name][pid].append(dz)
                    if method_name == "scanpy":
                        per_patient_zn[pid].append(zn)
                        per_patient_zt[pid].append(zt)

            draw_rec[f"n_patients_contributing_{method_name}"] = n_contributing

            if n_contributing > 0:
                arr = np.array(patient_dz)
                mean_dz = float(np.mean(arr))
                med_dz = float(np.median(arr))
                try:
                    _, wp = stats.wilcoxon(arr, alternative="two-sided",
                                           zero_method="wilcox", method="auto")
                    wp = float(wp)
                except ValueError:
                    wp = np.nan

                draw_rec[f"mean_delta_z_{method_name}"] = round(mean_dz, 6)
                draw_rec[f"median_delta_z_{method_name}"] = round(med_dz, 6)
                draw_rec[f"wilcoxon_p_{method_name}"] = round(wp, 6) if np.isfinite(wp) else None
                draw_rec[f"direction_positive_{method_name}"] = mean_dz > 0
                draw_rec[f"significant_positive_{method_name}"] = (wp < 0.05 and mean_dz > 0) if np.isfinite(wp) else False
            else:
                draw_rec[f"mean_delta_z_{method_name}"] = None
                draw_rec[f"direction_positive_{method_name}"] = None
                draw_rec[f"significant_positive_{method_name}"] = None

            # Valid-iteration guard for scanpy primary
            if method_name == "scanpy" and n_contributing < 15:
                valid_iteration_guard_fail = True

        draw_rows.append(draw_rec)

        if (draw_id + 1) % 200 == 0:
            print(f"  Iteration {draw_id + 1}/{K} complete", flush=True)

    draw_df = pd.DataFrame(draw_rows)
    # Trim child_state_json for CSV manageability (keep in provenance)
    draw_df_csv = draw_df.copy()
    draw_df_csv.to_csv(TABLES_DIR / "gate5_draw_summary.csv", index=False)
    rpt(f"Saved: tables/gate5_draw_summary.csv ({K} rows)")
    rpt()

    if valid_iteration_guard_fail:
        rpt("**STOP: NEEDS_HUMAN_CONFIRMATION — at least one scanpy iteration has <15 contributing patients**")
        _save_report(report, RESULTS_DIR / "gate5_equal_cell_report.md")
        return

    rpt("Valid-iteration guard: **PASSED** (all 1000 scanpy iterations have >=15 contributing patients)")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 2 — Frozen gate decision (scanpy primary only)
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 2: Frozen Gate Decision (A006 F1)\n")

    scanpy_means = draw_df["mean_delta_z_scanpy"].values.astype(np.float64)
    direction_pos = scanpy_means > 0
    direction_retention = float(np.mean(direction_pos))
    median_magnitude = float(np.median(scanpy_means))

    direction_pass = direction_retention >= DIRECTION_THRESHOLD
    magnitude_pass = median_magnitude >= MAGNITUDE_THRESHOLD

    rpt(f"direction_retention_fraction: {direction_retention:.4f} (threshold: >= {DIRECTION_THRESHOLD})")
    rpt(f"median_downsampled_magnitude: {median_magnitude:.6f} (threshold: >= {MAGNITUDE_THRESHOLD})")
    rpt(f"direction_pass: {direction_pass}")
    rpt(f"magnitude_pass: {magnitude_pass}")
    rpt()

    if direction_pass and magnitude_pass:
        outcome = "PASS"
    elif not direction_pass:
        outcome = "FAIL-DIRECTION-UNSTABLE"
    else:
        outcome = "FAIL-MAGNITUDE-COLLAPSE"

    rpt(f"**GATE 5 OUTCOME: {outcome}**")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 3 — Supplementary (NON-GATING)
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 3: Supplementary Analyses (NON-GATING)\n")
    rpt("*None of the following changes the Task 2 decision.*\n")

    # (a) Consolidated per-patient Delta_z
    consol_rows = []
    for pid in eligible_pids:
        rec = {"PID": pid}
        for method_name, _, _ in methods:
            dz_list = per_patient_dz[method_name][pid]
            rec[f"n_draws_contributing_{method_name}"] = len(dz_list)
            rec[f"consolidated_delta_z_{method_name}"] = round(float(np.mean(dz_list)), 6) if dz_list else None

        zn_list = per_patient_zn[pid]
        zt_list = per_patient_zt[pid]
        rec["consolidated_z_normal_scanpy"] = round(float(np.mean(zn_list)), 6) if zn_list else None
        rec["consolidated_z_tumor_scanpy"] = round(float(np.mean(zt_list)), 6) if zt_list else None
        consol_rows.append(rec)

    consol_df = pd.DataFrame(consol_rows)
    consol_df.to_csv(TABLES_DIR / "gate5_consolidated_per_patient.csv", index=False)

    # Gate-4-style summary on consolidated scanpy Delta_z
    consol_dz = consol_df["consolidated_delta_z_scanpy"].dropna().values.astype(np.float64)
    consol_test = donor_aware_test_on_dz(consol_dz)

    rpt("### (a) Consolidated Donor-Aware Summary (scanpy)\n")
    rpt(f"Consolidated mean Delta_z: {consol_test['mean_delta_z']:.6f}")
    rpt(f"Consolidated median Delta_z: {consol_test['median_delta_z']:.6f}")
    rpt(f"Wilcoxon stat: {consol_test['wilcoxon_stat']}, p: {consol_test['wilcoxon_p']:.6f}")
    rpt(f"  Exact zeros: {consol_test['n_exact_zeros']}, effective n: {consol_test['effective_n']}")
    rpt(f"Bootstrap 95% CI: [{consol_test['bootstrap_ci_lo']:.6f}, {consol_test['bootstrap_ci_hi']:.6f}]")
    rpt(f"  CI excludes 0: {consol_test['ci_excludes_0']}")
    rpt(f"LOO sign stable: {consol_test['loo_all_same_sign']}")
    rpt(f"LOO max % change: {consol_test['loo_max_rel_change_pct']:.2f}%")
    rpt(f"Proportion Delta_z > 0: {consol_test['prop_positive']:.4f}")
    rpt()

    # Correlation scale
    mean_z_n = float(consol_df["consolidated_z_normal_scanpy"].mean())
    mean_z_t = float(consol_df["consolidated_z_tumor_scanpy"].mean())
    rpt(f"tanh(consolidated mean z_normal): {np.tanh(mean_z_n):.6f}")
    rpt(f"tanh(consolidated mean z_tumor): {np.tanh(mean_z_t):.6f}")
    rpt(f"tanh difference: {np.tanh(mean_z_t) - np.tanh(mean_z_n):.6f}")
    rpt()

    # Save LOO
    loo_rows_out = []
    for i, pid in enumerate(eligible_pids):
        lr = consol_test["loo_rows"][i]
        lr["excluded_PID"] = pid
        loo_rows_out.append(lr)
    pd.DataFrame(loo_rows_out).to_csv(TABLES_DIR / "gate5_loo_primary.csv", index=False)
    rpt("Saved: tables/gate5_loo_primary.csv, gate5_consolidated_per_patient.csv")
    rpt()

    # (b) Descriptive draw fractions
    rpt("### (b) Descriptive Draw Fractions\n")
    sig_frac = float(draw_df["wilcoxon_p_scanpy"].dropna().apply(lambda p: p < 0.05).mean())
    sig_pos_frac = float(draw_df["significant_positive_scanpy"].dropna().mean())
    rpt(f"Significant draw fraction (Wilcoxon p<0.05): {sig_frac:.4f}")
    rpt(f"Significant-positive draw fraction: {sig_pos_frac:.4f}")
    rpt()

    # (c) Method sensitivity
    rpt("### (c) Method Sensitivity\n")
    method_sens_rows = []
    for method_name, _, _ in methods:
        col = f"mean_delta_z_{method_name}"
        vals = draw_df[col].dropna().values.astype(np.float64)
        dir_ret = float(np.mean(vals > 0))
        med_mag = float(np.median(vals))
        consol_col = f"consolidated_delta_z_{method_name}"
        consol_mean = float(consol_df[consol_col].dropna().mean())
        method_sens_rows.append({
            "method": method_name,
            "direction_retention_fraction": round(dir_ret, 4),
            "median_downsampled_magnitude": round(med_mag, 6),
            "consolidated_mean_delta_z": round(consol_mean, 6),
            "direction": "+" if consol_mean > 0 else "-" if consol_mean < 0 else "0",
        })
        rpt(f"  {method_name}: dir_ret={dir_ret:.4f}, med_mag={med_mag:.6f}, "
            f"consol_mean_Δz={consol_mean:.6f}")

    rpt()
    msens_df = pd.DataFrame(method_sens_rows)
    msens_df.to_csv(TABLES_DIR / "gate5_method_sensitivity.csv", index=False)

    all_positive = all(r["direction"] == "+" for r in method_sens_rows)
    rpt(f"All methods positive direction under equal-cell: {'YES' if all_positive else 'NO'}")
    rpt(f"Matches Gate 4 / E10 cross-method pattern: {'YES' if all_positive else 'NO'}")
    rpt()
    rpt("Saved: tables/gate5_method_sensitivity.csv")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 4 — Primary summary CSV
    # ══════════════════════════════════════════════════════════════════
    summary = {
        "K": K, "sample_size_per_PID_tissue": SAMPLE_SIZE,
        "seed_plan": "SeedSequence(42).spawn(1000)",
        "full_data_effect": FULL_DATA_EFFECT,
        "magnitude_threshold": MAGNITUDE_THRESHOLD,
        "direction_retention_fraction": round(direction_retention, 4),
        "median_downsampled_magnitude": round(median_magnitude, 6),
        "magnitude_pass": magnitude_pass,
        "direction_pass": direction_pass,
        "outcome_class": outcome,
        "consolidated_mean_delta_z": round(consol_test["mean_delta_z"], 6),
        "consolidated_wilcoxon_p": round(consol_test["wilcoxon_p"], 6),
        "consolidated_bootstrap_ci_low": round(consol_test["bootstrap_ci_lo"], 6),
        "consolidated_bootstrap_ci_high": round(consol_test["bootstrap_ci_hi"], 6),
        "consolidated_loo_sign_stable": consol_test["loo_all_same_sign"],
        "significant_draw_fraction": round(sig_frac, 4),
        "significant_positive_draw_fraction": round(sig_pos_frac, 4),
    }
    pd.DataFrame([summary]).to_csv(TABLES_DIR / "gate5_primary_summary.csv", index=False)
    rpt("Saved: tables/gate5_primary_summary.csv")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 5 — Mandatory statements
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 5: Mandatory Statements\n")
    rpt("**Batch claim ceiling (E11):** Gate 5 tests only within-patient cell-count")
    rpt("imbalance robustness. It does NOT address the perfect tissue-batch confound")
    rpt("(Cramer's V = 1.0; tissue and batch perfectly confounded, never independent).")
    rpt("Lee et al. external replication is mandatory for strong biological interpretation.")
    rpt()
    rpt("**Pseudoreplication caveat (E0):** Each cell is not an independent replicate.")
    rpt("The K-iteration distribution is a Monte Carlo subsampling stability distribution,")
    rpt("not a classical confidence interval.")
    rpt()
    rpt("**Power limitation (A11):** All analyses use the same 34-patient eligible set.")
    rpt()
    rpt("**Cell-cycle limitation (C6):** S/G2M scores absent; residual cell-cycle")
    rpt("variation within cE01 not fully assessable.")
    rpt()
    rpt("**Scope:** A PASS shows robustness to equal-cell subsampling; it does NOT prove")
    rpt("that cell-count imbalance was absent as a contributor, and it does NOT address")
    rpt("the tissue-batch confound.")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 6 — Holistic verdict
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 6: Holistic Verdict\n")
    rpt(f"**Gate 5 outcome: {outcome}**")
    rpt()

    g4_pass = True  # Gate 4 PASS established
    a9_survival = g4_pass and (outcome == "PASS")

    rpt(f"**A9 survival (Gate 4 PASS + Gate 5 {outcome}):** "
        f"{'MET — core survives' if a9_survival else 'NOT MET — core does not survive'}")
    rpt()

    if a9_survival:
        rpt("### Cross-method claim language (E10)")
        if all_positive:
            rpt("All four methods (scanpy, clean-bg, AUCell, z-score) show positive direction")
            rpt("under equal-cell subsampling, consistent with the Gate 4 pattern.")
            rpt("→ 'A shift toward a more positive ISC-FAO relationship in tumor-derived")
            rpt("  cE01 cells across scoring methods.'")
        rpt()
        rpt("### Batch ceiling (E11)")
        rpt("All claims remain hedged by the perfect tissue-batch confound.")
        rpt("Tissue and batch are perfectly confounded and not separable.")
        rpt("External replication (Lee et al.) is mandatory for any strong biological")
        rpt("interpretation.")
        rpt()
        rpt("### Project outcome pointer")
        rpt("With both Gate 4 and Gate 5 PASS, the core observation survives.")
        # NOTE: the output below is SUPERSEDED by Gate 3 and the post-audit addendum; see FINAL_REPORT.md
        rpt("Evidence points toward **Outcome 1** (proceed to secondaries) or")
        rpt("**Outcome 2** (differentiation/lipid-handling reframe pending Gate 3),")
        rpt("depending on Gate 3 decomposition results. Gate 3 is next.")
    else:
        rpt("The core does not survive Gate 5. Evidence points toward **Outcome 3** —")
        rpt("archive as methods/exploratory portfolio project.")

    rpt()

    # ── Save ──────────────────────────────────────────────────────────
    _save_report(report, RESULTS_DIR / "gate5_equal_cell_report.md")

    prov = {
        "script": "05_gate5_equal_cell.py",
        "gate": "Gate 5",
        "versions": versions,
        "seed_plan": "SeedSequence(42).spawn(1000)",
        "K": K,
        "sample_size": SAMPLE_SIZE,
        "full_data_effect": FULL_DATA_EFFECT,
        "magnitude_threshold": MAGNITUDE_THRESHOLD,
        "direction_threshold": DIRECTION_THRESHOLD,
        "sampling_order": "PIDs lexicographic, tissues [N, T], barcodes lexicographic",
        "hash_format": "SHA-256 of sorted PID<TAB>SPECIMEN_TYPE<TAB>barcode lines joined by newline",
        "drift_check": "PASSED",
        "eligible_pids": eligible_pids,
        "valid_iteration_guard": "PASSED",
        "direction_retention_fraction": round(direction_retention, 4),
        "median_downsampled_magnitude": round(median_magnitude, 6),
        "outcome_class": outcome,
        "inputs": {
            "gate2_scores": str(RESULTS_DIR / "gate2_recomputed_scores_cE01.parquet"),
            "gate2b_scores": str(RESULTS_DIR / "gate2b_clean_background_scores_cE01.parquet"),
            "gate1_json": str(RESULTS_DIR / "gate1_input_validation.json"),
            "config": str(CFG_PATH),
        },
        "status": "COMPLETE",
    }
    with open(PROV_DIR / "05_gate5_equal_cell.json", "w") as f:
        json.dump(prov, f, indent=2, default=str)
    print(f"Saved: provenance/05_gate5_equal_cell.json", flush=True)

    print(f"\n{'=' * 72}", flush=True)
    print(f"GATE 5 COMPLETE — OUTCOME: {outcome}", flush=True)
    print(f"{'=' * 72}", flush=True)


def _save_report(lines, path):
    with open(path, "w") as f:
        f.write("\n".join(lines))
    print(f"Saved: {path.name}", flush=True)


if __name__ == "__main__":
    main()
