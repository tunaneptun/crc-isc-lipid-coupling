"""
Gate 4 — Donor-Aware Coupling (THE DECISION GATE)
====================================================
Protocol: protocol_frozen_v6 (v1 + A001-A005 + DEVIATION_001)
Execution order: after Gate 2/2b, before Gate 3/5

Hard kill decided ONLY on official scanpy primary.
Other methods/models are sensitivities that bind claim language.
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import yaml

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG_PATH = Path(__file__).resolve().parent / "config.yaml"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
TABLES_DIR = RESULTS_DIR / "tables"
PROV_DIR = Path(__file__).resolve().parent / "provenance"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

SEED = 42
N_BOOT = 10000

# ═══════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════

def get_versions():
    import importlib.metadata as im
    return {k: im.version(k) for k in ["numpy", "scipy", "statsmodels", "pandas"]}


def is_near_zero_var(arr, atol=1e-15):
    """Population variance near zero, per frozen rule."""
    return np.isclose(np.var(arr, ddof=0), 0.0, rtol=0.0, atol=atol)


def safe_pearsonr(x, y):
    """Pearson r with undefined-correlation handling.
    Returns (r, status, reason).
    status: 'ok' or 'uninformative'.
    """
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
    """Clip r to +/-0.999999 for arctanh, return (clipped_r, was_clipped)."""
    if abs(r) >= 1.0:
        return np.sign(r) * 0.999999, True
    return r, False


def fisher_z(r):
    return np.arctanh(r)


def compute_per_patient(df, isc_col, fao_col, method_name, pid_col, tissue_col,
                        tumor_code, normal_code):
    """Compute per-patient r, z, Delta_z for one method."""
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

        # Determine retention
        if (rec.get("normal_status") == "ok" and rec.get("tumor_status") == "ok"):
            rec["retained_for_paired_analysis"] = True
            rn, clip_n = clip_r(rec["r_normal"])
            rt, clip_t = clip_r(rec["r_tumor"])
            rec["z_normal"] = fisher_z(rn)
            rec["z_tumor"] = fisher_z(rt)
            rec["Delta_z"] = rec["z_tumor"] - rec["z_normal"]
            rec["clip_flag"] = clip_n or clip_t
        else:
            rec["retained_for_paired_analysis"] = False
            rec["z_normal"] = np.nan
            rec["z_tumor"] = np.nan
            rec["Delta_z"] = np.nan
            rec["clip_flag"] = False

        rows.append(rec)

    return pd.DataFrame(rows)


def donor_aware_test(paired_df, seed=SEED, n_boot=N_BOOT):
    """Run Wilcoxon, bootstrap CI, LOO on paired Delta_z.
    Returns dict with results.
    """
    dz = paired_df["Delta_z"].values.astype(np.float64)
    n = len(dz)
    mean_dz = float(np.mean(dz))
    median_dz = float(np.median(dz))

    # Wilcoxon
    # Handle exact zeros
    n_zeros = int(np.sum(dz == 0))
    try:
        wstat, wp = stats.wilcoxon(dz, alternative="two-sided",
                                   zero_method="wilcox", method="auto")
    except ValueError as e:
        wstat, wp = np.nan, np.nan

    # Bootstrap
    rng = np.random.RandomState(seed)
    boot_means = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        boot_means[i] = np.mean(dz[idx])
    ci_lo = float(np.percentile(boot_means, 2.5))
    ci_hi = float(np.percentile(boot_means, 97.5))
    ci_excludes_0 = (ci_lo > 0) or (ci_hi < 0)

    # LOO
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
        if abs(mean_dz) > 0:
            rel_change = abs(loo_mean - mean_dz) / abs(mean_dz) * 100
        else:
            rel_change = 0.0
        loo_rows.append({
            "excluded_PID": paired_df.iloc[i]["PID"],
            "n_remaining": n - 1,
            "mean_delta_z": round(loo_mean, 6),
            "wilcoxon_p": round(loo_p, 6) if np.isfinite(loo_p) else None,
            "sign_positive": sign_pos,
            "relative_effect_change_pct": round(rel_change, 2),
        })

    loo_df = pd.DataFrame(loo_rows)
    all_same_sign = all(r["sign_positive"] == (mean_dz > 0) for r in loo_rows)
    max_rel_change = max(r["relative_effect_change_pct"] for r in loo_rows)
    loo_effect_sensitive = max_rel_change > 25.0

    return {
        "n": n,
        "mean_delta_z": mean_dz,
        "median_delta_z": median_dz,
        "wilcoxon_stat": float(wstat) if np.isfinite(wstat) else None,
        "wilcoxon_p": float(wp) if np.isfinite(wp) else None,
        "n_exact_zeros": n_zeros,
        "effective_n_after_zeros": n - n_zeros,
        "bootstrap_ci_lo": ci_lo,
        "bootstrap_ci_hi": ci_hi,
        "ci_excludes_0": ci_excludes_0,
        "loo_all_same_sign": all_same_sign,
        "loo_max_rel_change_pct": max_rel_change,
        "loo_effect_sensitive": loo_effect_sensitive,
        "proportion_delta_z_positive": float(np.mean(dz > 0)),
        "loo_df": loo_df,
    }


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 72, flush=True)
    print("Gate 4 — Donor-Aware Coupling (DECISION GATE)", flush=True)
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

    assert g2.index.equals(g2b.index), "Barcode index mismatch between Gate 2 and Gate 2b"
    df = g2.join(g2b[["ISC_scanpy_cleanbg", "FAO_scanpy_cleanbg"]])
    print(f"\nJoined cE01 score table: {len(df):,} cells", flush=True)
    assert len(df) == 61953, f"Expected 61953 cells, got {len(df)}"

    report = []
    def rpt(line=""):
        report.append(line)
        print(line, flush=True)

    rpt("# Gate 4 — Donor-Aware Coupling Report\n")
    rpt(f"**Protocol:** protocol_frozen_v6")
    rpt(f"**Versions:** {versions}")
    rpt(f"**Cells:** {len(df):,}")
    rpt(f"**Seed:** {SEED}, **Bootstrap iterations:** {N_BOOT}")
    rpt()

    # ── Eligibility & drift check ─────────────────────────────────────
    rpt("## Eligibility & Drift Check\n")

    thresh = cfg["eligibility"]["primary_min_cells"]
    pid_counts = df.groupby([pid_col, tissue_col]).size().unstack(fill_value=0)
    eligible_mask = (pid_counts.get(tumor_code, pd.Series(dtype=int)) >= thresh) & \
                    (pid_counts.get(normal_code, pd.Series(dtype=int)) >= thresh)
    eligible_pids = sorted(eligible_mask[eligible_mask].index.tolist())
    n_eligible = len(eligible_pids)

    rpt(f"Primary threshold: >= {thresh} cells in both T and N")
    rpt(f"Eligible patients: {n_eligible}")

    # Drift check against Gate 1
    with open(RESULTS_DIR / "gate1_input_validation.json") as f:
        g1 = json.load(f)
    g1_pids = sorted([p["patient"] for p in g1["checks"]["eligibility"]["thresh_30"]["patients"]])

    drift_match = eligible_pids == g1_pids
    rpt(f"Gate 1 eligible at >=30: {len(g1_pids)}")
    rpt(f"Gate 4 eligible at >=30: {n_eligible}")
    rpt(f"PID membership match: {'YES' if drift_match else 'NO — STOP'}")

    if not drift_match:
        rpt("\n**VERDICT: NEEDS_HUMAN_CONFIRMATION — eligibility drift detected**")
        _save_report(report, RESULTS_DIR / "gate4_donor_aware_report.md")
        return

    if n_eligible < 15:
        rpt(f"\n**VERDICT: POWER-LIMITED (n={n_eligible} < 15)**")
        _save_report(report, RESULTS_DIR / "gate4_donor_aware_report.md")
        return

    rpt(f"Drift check: **PASSED** (34 patients, identical to Gate 1)")
    rpt()

    # Report per-patient cell counts
    elig_df = df[df[pid_col].isin(eligible_pids)]
    for pid in eligible_pids:
        nt = int((elig_df[elig_df[pid_col] == pid][tissue_col] == tumor_code).sum())
        nn = int((elig_df[elig_df[pid_col] == pid][tissue_col] == normal_code).sum())

    # ══════════════════════════════════════════════════════════════════
    # TASK 1 — Per-patient Fisher-z (all four methods)
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 1: Per-Patient Fisher-z\n")

    methods = [
        ("scanpy", "ISC_scanpy", "FAO_scanpy"),
        ("clean_bg", "ISC_scanpy_cleanbg", "FAO_scanpy_cleanbg"),
        ("aucell", "ISC_aucell", "FAO_aucell"),
        ("zscore", "ISC_zscore", "FAO_zscore"),
    ]

    all_pp = []
    method_paired = {}  # method -> paired DataFrame

    for method_name, isc_col, fao_col in methods:
        pp = compute_per_patient(elig_df, isc_col, fao_col, method_name,
                                 pid_col, tissue_col, tumor_code, normal_code)
        all_pp.append(pp)
        paired = pp[pp["retained_for_paired_analysis"] == True].copy()
        excluded = pp[pp["retained_for_paired_analysis"] == False]
        method_paired[method_name] = paired

        rpt(f"**{method_name}:** {len(paired)} retained, {len(excluded)} excluded")
        if len(excluded) > 0:
            for _, row in excluded.iterrows():
                reasons = []
                if row.get("normal_status") == "uninformative":
                    reasons.append(f"normal: {row.get('normal_exclusion_reason','?')}")
                if row.get("tumor_status") == "uninformative":
                    reasons.append(f"tumor: {row.get('tumor_exclusion_reason','?')}")
                rpt(f"  Excluded {row['PID']}: {'; '.join(reasons)}")

    rpt()

    pp_all = pd.concat(all_pp, ignore_index=True)
    pp_all.to_csv(TABLES_DIR / "gate4_per_patient.csv", index=False)
    rpt("Saved: tables/gate4_per_patient.csv")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 2 — Primary donor-aware test (scanpy only)
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 2: Primary Donor-Aware Test (Official Scanpy)\n")
    rpt("Wilcoxon and the bootstrap CI evaluate complementary summaries,")
    rpt("not the same estimand: Wilcoxon tests for a consistent directional")
    rpt("shift (median-like); the bootstrap CI estimates the mean effect size.")
    rpt()

    scanpy_paired = method_paired["scanpy"]
    if len(scanpy_paired) < 15:
        rpt(f"**VERDICT: POWER-LIMITED (n={len(scanpy_paired)} < 15 retained)**")
        _save_report(report, RESULTS_DIR / "gate4_donor_aware_report.md")
        return

    primary = donor_aware_test(scanpy_paired)

    rpt(f"Retained patients: {primary['n']}")
    rpt(f"Mean Delta_z: {primary['mean_delta_z']:.6f}")
    rpt(f"Median Delta_z: {primary['median_delta_z']:.6f}")
    rpt(f"Proportion Delta_z > 0: {primary['proportion_delta_z_positive']:.4f}")
    rpt(f"Wilcoxon stat: {primary['wilcoxon_stat']}, p = {primary['wilcoxon_p']:.6f}")
    rpt(f"  Exact zeros: {primary['n_exact_zeros']}, effective n: {primary['effective_n_after_zeros']}")
    rpt(f"  scipy version: {versions['scipy']}")
    rpt(f"Bootstrap 95% CI: [{primary['bootstrap_ci_lo']:.6f}, {primary['bootstrap_ci_hi']:.6f}]")
    rpt(f"  CI excludes 0: {primary['ci_excludes_0']}")
    rpt(f"LOO sign stability: all same sign = {primary['loo_all_same_sign']}")
    rpt(f"LOO max % change: {primary['loo_max_rel_change_pct']:.2f}%")
    rpt(f"LOO effect-size-sensitive flag: {primary['loo_effect_sensitive']}")
    rpt()

    # Save LOO
    primary["loo_df"].to_csv(TABLES_DIR / "gate4_loo_primary.csv", index=False)
    rpt("Saved: tables/gate4_loo_primary.csv")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 3 — Effect size, correlation scale
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 3: Effect Size on Correlation Scale\n")

    sp = scanpy_paired
    mean_r_n = float(sp["r_normal"].mean())
    mean_r_t = float(sp["r_tumor"].mean())
    med_r_n = float(sp["r_normal"].median())
    med_r_t = float(sp["r_tumor"].median())
    mean_z_n = float(sp["z_normal"].mean())
    mean_z_t = float(sp["z_tumor"].mean())
    tanh_z_n = float(np.tanh(mean_z_n))
    tanh_z_t = float(np.tanh(mean_z_t))

    rpt(f"Mean raw r_normal: {mean_r_n:.6f}")
    rpt(f"Mean raw r_tumor: {mean_r_t:.6f}")
    rpt(f"Median raw r_normal: {med_r_n:.6f}")
    rpt(f"Median raw r_tumor: {med_r_t:.6f}")
    rpt(f"Mean z_normal: {mean_z_n:.6f}")
    rpt(f"Mean z_tumor: {mean_z_t:.6f}")
    rpt(f"Mean Delta_z: {primary['mean_delta_z']:.6f}")
    rpt(f"Median Delta_z: {primary['median_delta_z']:.6f}")
    rpt(f"tanh(mean_z_normal): {tanh_z_n:.6f}")
    rpt(f"tanh(mean_z_tumor): {tanh_z_t:.6f}")
    rpt(f"tanh(mean_z_tumor) - tanh(mean_z_normal): {tanh_z_t - tanh_z_n:.6f}")
    rpt(f"Proportion Delta_z > 0: {primary['proportion_delta_z_positive']:.4f}")
    rpt()

    summary_row = {
        "mean_r_normal": round(mean_r_n, 6), "mean_r_tumor": round(mean_r_t, 6),
        "median_r_normal": round(med_r_n, 6), "median_r_tumor": round(med_r_t, 6),
        "mean_z_normal": round(mean_z_n, 6), "mean_z_tumor": round(mean_z_t, 6),
        "mean_delta_z": round(primary["mean_delta_z"], 6),
        "median_delta_z": round(primary["median_delta_z"], 6),
        "tanh_mean_z_normal": round(tanh_z_n, 6), "tanh_mean_z_tumor": round(tanh_z_t, 6),
        "tanh_diff": round(tanh_z_t - tanh_z_n, 6),
        "prop_delta_z_positive": round(primary["proportion_delta_z_positive"], 4),
        "wilcoxon_p": primary["wilcoxon_p"],
        "bootstrap_ci_lo": primary["bootstrap_ci_lo"],
        "bootstrap_ci_hi": primary["bootstrap_ci_hi"],
        "n_patients": primary["n"],
    }
    pd.DataFrame([summary_row]).to_csv(TABLES_DIR / "gate4_primary_summary.csv", index=False)
    rpt("Saved: tables/gate4_primary_summary.csv")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 4 — Outcome classification
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 4: Primary Outcome Classification\n")
    rpt("Decision tree (frozen rules):")
    rpt()

    n_retained = primary["n"]
    mean_dz = primary["mean_delta_z"]
    wp = primary["wilcoxon_p"]
    ci_excl = primary["ci_excludes_0"]
    loo_sign = primary["loo_all_same_sign"]
    loo_sens = primary["loo_effect_sensitive"]

    rpt(f"1. n_retained = {n_retained} >= 15? {'YES' if n_retained >= 15 else 'NO -> POWER-LIMITED'}")

    if n_retained < 15:
        outcome = "POWER-LIMITED"
    elif mean_dz <= 0:
        outcome = "DIRECTION CONTRADICTED"
        rpt(f"2. mean_delta_z = {mean_dz:.6f} > 0? NO -> DIRECTION CONTRADICTED")
        # Check if the negative is itself significant
        neg_sig = wp < 0.05 and ci_excl
        rpt(f"   (Negative mean Delta_z significant? Wilcoxon p={wp:.6f}, CI excl 0={ci_excl} -> {'YES' if neg_sig else 'NO'})")
        rpt(f"   Note: DIRECTION CONTRADICTED means the point estimate is not in the hypothesized direction,")
        rpt(f"   NOT that a reverse effect is established.")
    else:
        rpt(f"2. mean_delta_z = {mean_dz:.6f} > 0? YES")
        rpt(f"3. Wilcoxon p = {wp:.6f} < 0.05? {'YES' if wp < 0.05 else 'NO'}")
        rpt(f"4. Bootstrap CI excludes 0? {ci_excl} [{primary['bootstrap_ci_lo']:.6f}, {primary['bootstrap_ci_hi']:.6f}]")
        rpt(f"5. LOO all same sign? {loo_sign}")
        rpt(f"6. LOO max % change = {primary['loo_max_rel_change_pct']:.2f}% <= 25%? {not loo_sens}")

        if wp < 0.05 and ci_excl and loo_sign and not loo_sens:
            outcome = "PASS"
        elif not loo_sign:
            outcome = "LOO-UNSTABLE"
        else:
            outcome = "NOT ESTABLISHED"

    rpt(f"\n**PRIMARY OUTCOME: {outcome}**")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 5 — Method sensitivity + clean-bg D2 classification
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 5: Method Sensitivity & Clean-BG Classification\n")

    sens_rows = []
    for method_name in ["scanpy", "clean_bg", "aucell", "zscore"]:
        paired = method_paired[method_name]
        if len(paired) < 3:
            sens_rows.append({"method": method_name, "n_retained": len(paired),
                              "mean_delta_z": np.nan, "status": "insufficient"})
            rpt(f"**{method_name}:** insufficient retained patients ({len(paired)})")
            continue
        res = donor_aware_test(paired)
        sens_rows.append({
            "method": method_name,
            "n_retained": res["n"],
            "n_excluded": n_eligible - res["n"],
            "mean_delta_z": round(res["mean_delta_z"], 6),
            "median_delta_z": round(res["median_delta_z"], 6),
            "wilcoxon_p": round(res["wilcoxon_p"], 6) if res["wilcoxon_p"] is not None else None,
            "bootstrap_ci_lo": round(res["bootstrap_ci_lo"], 6),
            "bootstrap_ci_hi": round(res["bootstrap_ci_hi"], 6),
            "ci_excludes_0": res["ci_excludes_0"],
            "proportion_delta_z_positive": round(res["proportion_delta_z_positive"], 4),
            "loo_all_same_sign": res["loo_all_same_sign"],
            "loo_max_rel_change_pct": round(res["loo_max_rel_change_pct"], 2),
        })
        rpt(f"**{method_name}:** n={res['n']}, mean_Δz={res['mean_delta_z']:.4f}, "
            f"med_Δz={res['median_delta_z']:.4f}, Wilcoxon p={res['wilcoxon_p']:.4f}, "
            f"CI=[{res['bootstrap_ci_lo']:.4f},{res['bootstrap_ci_hi']:.4f}], "
            f"prop>0={res['proportion_delta_z_positive']:.3f}")

    rpt()
    pd.DataFrame(sens_rows).to_csv(TABLES_DIR / "gate4_method_sensitivity.csv", index=False)
    rpt("Saved: tables/gate4_method_sensitivity.csv")
    rpt()

    # Clean-BG D2 classification (on COMMON retained PID set)
    rpt("### Clean-BG vs Official-Scanpy D2 Classification\n")

    scanpy_pids = set(method_paired["scanpy"]["PID"].tolist())
    cleanbg_pids = set(method_paired["clean_bg"]["PID"].tolist())
    common_pids = sorted(scanpy_pids & cleanbg_pids)
    rpt(f"Official scanpy retained: {len(scanpy_pids)}")
    rpt(f"Clean-bg retained: {len(cleanbg_pids)}")
    rpt(f"Common retained: {len(common_pids)}")
    pid_diff = scanpy_pids.symmetric_difference(cleanbg_pids)
    if pid_diff:
        rpt(f"Membership differences: {sorted(pid_diff)}")
    rpt()

    if len(common_pids) < 15:
        d2_class = "POWER-LIMITED / UNINFORMATIVE"
        rpt(f"Common retained < 15 -> {d2_class}")
    else:
        scanpy_common = method_paired["scanpy"][method_paired["scanpy"]["PID"].isin(common_pids)]
        cleanbg_common = method_paired["clean_bg"][method_paired["clean_bg"]["PID"].isin(common_pids)]
        # Align by PID
        sc_dz = scanpy_common.set_index("PID")["Delta_z"]
        cb_dz = cleanbg_common.set_index("PID")["Delta_z"]
        sc_mean = float(sc_dz.mean())
        cb_mean = float(cb_dz.mean())
        ratio = abs(cb_mean) / abs(sc_mean) if abs(sc_mean) > 0 else float("inf")

        rpt(f"Scanpy mean Delta_z (common): {sc_mean:.6f}")
        rpt(f"Clean-bg mean Delta_z (common): {cb_mean:.6f}")
        rpt(f"Magnitude ratio |clean|/|primary|: {ratio:.4f}")

        if np.sign(cb_mean) != np.sign(sc_mean):
            d2_class = "REVERSED"
        elif ratio >= 0.50:
            d2_class = "DIRECTION PRESERVED"
        else:
            d2_class = "MATERIALLY ATTENUATED"

        rpt(f"**D2 Classification: {d2_class}**")
    rpt()

    # Cross-method direction pattern (E10)
    rpt("### Cross-Method Direction Pattern\n")
    all_positive = all(
        s.get("mean_delta_z", 0) > 0
        for s in sens_rows if s.get("n_retained", 0) >= 3
        and not np.isnan(s.get("mean_delta_z", np.nan))
    )
    any_reversed = any(
        s.get("mean_delta_z", 0) < 0
        for s in sens_rows if s.get("n_retained", 0) >= 3
        and not np.isnan(s.get("mean_delta_z", np.nan))
    )

    for s in sens_rows:
        if s.get("n_retained", 0) >= 3:
            dir_str = "+" if s.get("mean_delta_z", 0) > 0 else "-" if s.get("mean_delta_z", 0) < 0 else "0"
            rpt(f"  {s['method']}: {dir_str} (mean_Δz = {s.get('mean_delta_z', 'N/A')})")

    if all_positive:
        claim_lang = "All methods agree in direction (positive Delta_z = more positive coupling in tumor). Strong cross-method consistency supports the program-level claim language."
    elif any_reversed:
        claim_lang = "At least one method reverses direction. Claim limited to method-dependent positive-shift statement."
    else:
        claim_lang = "Mixed or zero effects. No strong program-level claim."

    rpt(f"\n{claim_lang}")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 6 — IVW sensitivity
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 6: IVW Sensitivity (Descriptive)\n")

    sp = scanpy_paired.copy()
    sp["se_dz"] = np.sqrt(1.0 / (sp["n_tumor"] - 3) + 1.0 / (sp["n_normal"] - 3))
    sp["w"] = 1.0 / sp["se_dz"] ** 2

    ivw_mean = float(np.sum(sp["w"] * sp["Delta_z"]) / np.sum(sp["w"]))

    # Patient-level bootstrap for IVW
    rng = np.random.RandomState(SEED)
    dz_arr = sp["Delta_z"].values
    w_arr = sp["w"].values
    n_ivw = len(sp)
    boot_ivw = np.empty(N_BOOT)
    for i in range(N_BOOT):
        idx = rng.choice(n_ivw, size=n_ivw, replace=True)
        boot_ivw[i] = np.sum(w_arr[idx] * dz_arr[idx]) / np.sum(w_arr[idx])
    ivw_ci_lo = float(np.percentile(boot_ivw, 2.5))
    ivw_ci_hi = float(np.percentile(boot_ivw, 97.5))

    rpt(f"IVW mean Delta_z: {ivw_mean:.6f}")
    rpt(f"IVW bootstrap 95% CI: [{ivw_ci_lo:.6f}, {ivw_ci_hi:.6f}]")
    rpt(f"Note: weights approximate (cells not independent replicates)")
    rpt()

    ivw_row = {"method": "scanpy", "ivw_mean_delta_z": round(ivw_mean, 6),
               "bootstrap_ci_low": round(ivw_ci_lo, 6),
               "bootstrap_ci_high": round(ivw_ci_hi, 6),
               "n_patients": n_ivw, "seed": SEED,
               "note": "weights_approximate_cells_not_independent"}
    pd.DataFrame([ivw_row]).to_csv(TABLES_DIR / "gate4_ivw_summary.csv", index=False)
    rpt("Saved: tables/gate4_ivw_summary.csv")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 7 — Mixed-effects model
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 7: Mixed-Effects Model\n")

    me_df = elig_df[[pid_col, tissue_col, "ISC_scanpy", "FAO_scanpy"]].copy()
    me_df = me_df.dropna()
    me_df[tissue_col] = me_df[tissue_col].astype(str)

    fit_kwargs = {"reml": False, "method": "lbfgs", "maxiter": 500, "disp": False}
    me_results = []

    # Preferred: random slope
    for model_type, re_formula in [("random_slope", "~ISC_scanpy"), ("random_intercept", "~1")]:
        formula = f"FAO_scanpy ~ ISC_scanpy * C({tissue_col}, Treatment(reference='N'))"
        rpt(f"### {model_type}")
        rpt(f"Formula: {formula}")
        rpt(f"re_formula: {re_formula}")
        rpt(f"fit kwargs: {fit_kwargs}")

        converged = False
        singular = False
        warnings_list = []
        interaction_coef = interaction_se = interaction_p = None
        interaction_ci_lo = interaction_ci_hi = None
        interp_level = "UNINFORMATIVE"

        try:
            import io, contextlib
            warn_buf = io.StringIO()
            with warnings.catch_warnings(record=True) as caught_warnings:
                warnings.simplefilter("always")
                model = smf.mixedlm(formula, data=me_df, groups=me_df[pid_col],
                                    re_formula=re_formula)
                result = model.fit(**fit_kwargs)

            warnings_list = [str(w.message) for w in caught_warnings]
            converged = result.converged

            # Check singularity
            re_cov = result.cov_re
            if hasattr(re_cov, 'values'):
                re_cov_arr = re_cov.values
            else:
                re_cov_arr = np.atleast_2d(re_cov)
            if np.any(np.abs(np.diag(re_cov_arr)) < 1e-10):
                singular = True

            # Extract interaction term
            interaction_term = f"ISC_scanpy:C({tissue_col}, Treatment(reference='N'))[T.T]"
            if interaction_term in result.params:
                interaction_coef = float(result.params[interaction_term])
                interaction_se = float(result.bse[interaction_term])
                interaction_p = float(result.pvalues[interaction_term])
                ci = result.conf_int().loc[interaction_term]
                interaction_ci_lo = float(ci[0])
                interaction_ci_hi = float(ci[1])

            if converged and not singular:
                if model_type == "random_slope":
                    interp_level = "CORROBORATIVE"
                else:
                    interp_level = "SUGGESTIVE ONLY"
            elif converged and singular:
                interp_level = "UNINFORMATIVE (singular)"
            else:
                interp_level = "UNINFORMATIVE (non-converged)"

        except Exception as e:
            warnings_list.append(f"Exception: {str(e)}")
            interp_level = "UNINFORMATIVE (exception)"

        rpt(f"Converged: {converged}")
        rpt(f"Singular: {singular}")
        if warnings_list:
            rpt(f"Warnings: {warnings_list[:3]}")
        if interaction_coef is not None:
            rpt(f"Interaction ISC_scanpy:tissue[T]: coef={interaction_coef:.6f}, "
                f"SE={interaction_se:.6f}, p={interaction_p:.6f}, "
                f"CI=[{interaction_ci_lo:.6f}, {interaction_ci_hi:.6f}]")
            rpt(f"Interpretation: positive coef = more positive ISC-FAO slope in tumor")
        rpt(f"Interpretation level: {interp_level}")
        rpt()

        me_results.append({
            "model_type": model_type,
            "formula": formula,
            "groups": pid_col,
            "re_formula": re_formula,
            "fit_kwargs": str(fit_kwargs),
            "n_cells": len(me_df),
            "n_patients": me_df[pid_col].nunique(),
            "converged": converged,
            "singular_fit": singular,
            "optimizer": "lbfgs",
            "warnings": str(warnings_list[:5]),
            "interaction_term": interaction_term if interaction_coef is not None else None,
            "interaction_coef": interaction_coef,
            "interaction_se": interaction_se,
            "interaction_ci_low": interaction_ci_lo,
            "interaction_ci_high": interaction_ci_hi,
            "interaction_p": interaction_p,
            "interpretation_level": interp_level,
        })

        # If random slope succeeded, skip fallback
        if converged and not singular and model_type == "random_slope":
            rpt("Random-slope model converged and non-singular; fallback not needed.")
            break

    pd.DataFrame(me_results).to_csv(TABLES_DIR / "gate4_mixed_models.csv", index=False)
    rpt("Saved: tables/gate4_mixed_models.csv")
    rpt()
    rpt("**Note:** This does not resolve tissue-batch confounding (perfectly confounded,")
    rpt("Cramer's V = 1.0) and does not override the per-patient primary.")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 8 — QC residualization (stratum-level)
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 8: QC Residualization (Stratum-Level)\n")

    qc_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    qc_strata_rows = []
    qc_resid_per_patient = []

    for pid in eligible_pids:
        for tissue_code in [normal_code, tumor_code]:
            stratum = elig_df[(elig_df[pid_col] == pid) & (elig_df[tissue_col] == tissue_code)].copy()
            n_before = len(stratum)
            tissue_label = "normal" if tissue_code == normal_code else "tumor"
            strata_rec = {"PID": pid, "tissue": tissue_label, "n_before": n_before}

            if n_before < 10:
                strata_rec.update({"n_after": n_before, "rank": None, "cond": None,
                                   "status": "QC-UNINFORMATIVE", "reason": f"n={n_before}<10"})
                qc_strata_rows.append(strata_rec)
                continue

            # Prepare QC matrix
            stratum = stratum.dropna(subset=qc_cols + ["ISC_scanpy", "FAO_scanpy"])
            n_after = len(stratum)

            X_qc = np.column_stack([
                np.log1p(stratum["total_counts"].values),
                stratum["n_genes_by_counts"].values,
                stratum["pct_counts_mt"].values,
            ])
            X_qc = sm.add_constant(X_qc)
            rank = np.linalg.matrix_rank(X_qc)
            cond = float(np.linalg.cond(X_qc))

            if rank < X_qc.shape[1]:
                strata_rec.update({"n_after": n_after, "rank": rank,
                                   "cond": round(cond, 2),
                                   "status": "QC-UNINFORMATIVE",
                                   "reason": "rank-deficient"})
                qc_strata_rows.append(strata_rec)
                continue

            residual_df = n_after - X_qc.shape[1]
            if residual_df < 10:
                strata_rec.update({"n_after": n_after, "rank": rank,
                                   "cond": round(cond, 2),
                                   "status": "QC-UNINFORMATIVE",
                                   "reason": f"residual_df={residual_df}<10"})
                qc_strata_rows.append(strata_rec)
                continue

            # Residualize ISC and FAO
            try:
                isc_resid = sm.OLS(stratum["ISC_scanpy"].values, X_qc).fit().resid
                fao_resid = sm.OLS(stratum["FAO_scanpy"].values, X_qc).fit().resid
            except Exception as e:
                strata_rec.update({"n_after": n_after, "rank": rank,
                                   "cond": round(cond, 2),
                                   "status": "QC-UNINFORMATIVE",
                                   "reason": f"OLS error: {e}"})
                qc_strata_rows.append(strata_rec)
                continue

            # Pearson on residuals
            r_resid, status_resid, reason_resid = safe_pearsonr(isc_resid, fao_resid)

            strata_rec.update({"n_after": n_after, "rank": rank,
                               "cond": round(cond, 2),
                               "status": status_resid,
                               "reason": reason_resid,
                               "r_resid": round(r_resid, 6) if np.isfinite(r_resid) else None})
            qc_strata_rows.append(strata_rec)

            if status_resid == "ok":
                qc_resid_per_patient.append({
                    "PID": pid, "tissue": tissue_label,
                    "r_resid": r_resid, "n": n_after,
                })

    qc_strata_df = pd.DataFrame(qc_strata_rows)
    qc_strata_df.to_csv(TABLES_DIR / "gate4_qc_residualization.csv", index=False)

    # Paired analysis on QC-residualized r
    qc_rp = pd.DataFrame(qc_resid_per_patient)
    qc_paired_pids = []
    if len(qc_rp) > 0:
        for pid in eligible_pids:
            pid_qc = qc_rp[qc_rp["PID"] == pid]
            has_n = "normal" in pid_qc["tissue"].values
            has_t = "tumor" in pid_qc["tissue"].values
            if has_n and has_t:
                qc_paired_pids.append(pid)

    n_qc_paired = len(qc_paired_pids)
    n_strata_total = len(qc_strata_rows)
    n_strata_info = sum(1 for r in qc_strata_rows if r.get("status") == "ok")
    n_strata_rd = sum(1 for r in qc_strata_rows if r.get("reason", "").startswith("rank"))
    n_strata_low = sum(1 for r in qc_strata_rows if "residual_df" in r.get("reason", "") or "n=" in r.get("reason", ""))

    rpt(f"Total strata: {n_strata_total}")
    rpt(f"Informative strata: {n_strata_info}")
    rpt(f"Rank-deficient: {n_strata_rd}")
    rpt(f"Low df / small n: {n_strata_low}")
    rpt(f"Paired patients (both tissues informative): {n_qc_paired}")
    rpt()

    qc_summary = {"n_strata_total": n_strata_total, "n_strata_informative": n_strata_info,
                   "n_strata_rank_deficient": n_strata_rd, "n_strata_low_df": n_strata_low,
                   "n_paired_patients_retained": n_qc_paired}

    if n_qc_paired < 15:
        qc_summary.update({"mean_delta_z": None, "median_delta_z": None,
                           "wilcoxon_p": None, "bootstrap_ci_low": None,
                           "bootstrap_ci_high": None,
                           "qc_sensitivity_status": "POWER-LIMITED/UNINFORMATIVE"})
        rpt(f"QC sensitivity: POWER-LIMITED/UNINFORMATIVE (n={n_qc_paired} < 15)")
    else:
        # Build paired Delta_z
        qc_paired_rows = []
        for pid in qc_paired_pids:
            pid_qc = qc_rp[qc_rp["PID"] == pid]
            r_n = float(pid_qc[pid_qc["tissue"] == "normal"]["r_resid"].values[0])
            r_t = float(pid_qc[pid_qc["tissue"] == "tumor"]["r_resid"].values[0])
            rn_c, _ = clip_r(r_n)
            rt_c, _ = clip_r(r_t)
            qc_paired_rows.append({"PID": pid, "r_normal": r_n, "r_tumor": r_t,
                                   "z_normal": fisher_z(rn_c), "z_tumor": fisher_z(rt_c),
                                   "Delta_z": fisher_z(rt_c) - fisher_z(rn_c)})

        qc_pp = pd.DataFrame(qc_paired_rows)
        qc_test = donor_aware_test(qc_pp)

        qc_summary.update({
            "mean_delta_z": round(qc_test["mean_delta_z"], 6),
            "median_delta_z": round(qc_test["median_delta_z"], 6),
            "wilcoxon_p": round(qc_test["wilcoxon_p"], 6) if qc_test["wilcoxon_p"] is not None else None,
            "bootstrap_ci_low": round(qc_test["bootstrap_ci_lo"], 6),
            "bootstrap_ci_high": round(qc_test["bootstrap_ci_hi"], 6),
            "qc_sensitivity_status": "COMPUTED",
        })
        rpt(f"QC-residualized donor-aware: mean_Δz={qc_test['mean_delta_z']:.6f}, "
            f"median_Δz={qc_test['median_delta_z']:.6f}")
        rpt(f"Wilcoxon p={qc_test['wilcoxon_p']:.6f}, "
            f"CI=[{qc_test['bootstrap_ci_lo']:.6f}, {qc_test['bootstrap_ci_hi']:.6f}]")
        rpt(f"Weakening here does not by itself prove an artifact.")

    rpt()
    pd.DataFrame([qc_summary]).to_csv(TABLES_DIR / "gate4_qc_primary_summary.csv", index=False)
    rpt("Saved: tables/gate4_qc_primary_summary.csv, gate4_qc_residualization.csv")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 9 — Mandatory statements
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 9: Mandatory Statements\n")

    rpt("**Pseudoreplication caveat (E0):** Each cell is not an independent replicate.")
    rpt("The per-patient paired Fisher-z framework treats the patient as the unit of")
    rpt("inference, but cell-level correlations within a patient remain vulnerable")
    rpt("to technical and biological confounds at the cell level.")
    rpt()
    rpt("**Batch claim ceiling (E11):** Tissue type and sequencing batch (batchID)")
    rpt("are perfectly confounded within the paired patients (Cramer's V = 1.0;")
    rpt("0/34 primary patients share a tumor-normal batch). Tissue and batch are")
    rpt("NOT separable — they are perfectly confounded. Batch adjustment is")
    rpt("non-identifiable. The strongest achievable causal claim remains hedged:")
    rpt("robust to measured QC and to control-construction choices, but")
    rpt("batch-confounded by design and requiring external replication (Lee et al.)")
    rpt("for any strong biological interpretation.")
    rpt()
    rpt("**Cell-cycle limitation (C6):** S/G2M cell-cycle scores are absent from")
    rpt("the metadata; residual cell-cycle variation within cE01 cannot be fully")
    rpt("assessed. Primary analysis excludes Stem/TA-like_prolif.")
    rpt()
    rpt(f"**Power limitation (A11):** All donor-aware analyses use the same")
    rpt(f"{n_eligible}-patient eligible set. If any branch is power-limited,")
    rpt(f"all are reported as power-limited together.")
    rpt()

    # ══════════════════════════════════════════════════════════════════
    # TASK 10 — Holistic verdict
    # ══════════════════════════════════════════════════════════════════
    rpt("## Task 10: Holistic Verdict\n")

    rpt(f"**Hard-kill outcome (Task 4): {outcome}**")
    rpt()

    if outcome == "PASS":
        rpt("The primary donor-aware per-patient Fisher-z test PASSES all three criteria.")
        rpt("Survival still requires Gate 5 (cell-count-imbalance robustness).")
        rpt()
        rpt("### Sensitivity synthesis:")
        rpt(f"- Clean-bg D2 classification: {d2_class}")
        rpt(f"- Cross-method direction: {'all positive' if all_positive else 'mixed/reversed'}")
        rpt(f"- QC residualization: {qc_summary.get('qc_sensitivity_status', 'N/A')}")

        if qc_summary.get("mean_delta_z") is not None:
            rpt(f"  (mean_Δz = {qc_summary['mean_delta_z']}, p = {qc_summary.get('wilcoxon_p')})")

        me_interp = me_results[0]["interpretation_level"] if me_results else "N/A"
        me_coef = me_results[0].get("interaction_coef") if me_results else None
        rpt(f"- Mixed-model: {me_interp}")
        if me_coef is not None:
            rpt(f"  (interaction coef = {me_coef:.6f})")
        rpt(f"- IVW: mean_Δz = {ivw_mean:.6f}, CI = [{ivw_ci_lo:.6f}, {ivw_ci_hi:.6f}]")
        rpt()

        # Map to claim language
        if d2_class == "DIRECTION PRESERVED" and all_positive:
            rpt("**Claim language (E10):** Program-level claim supported — the tumor-vs-normal")
            rpt("attenuation of ISC-FAO cell-level coupling is consistent across scoring methods")
            rpt("and robust to control-pool construction.")
        elif d2_class in ("MATERIALLY ATTENUATED", "REVERSED") or any_reversed:
            rpt("**Claim language (E10):** Claim downgrades to method-dependent positive-shift")
            rpt("statement. No strong program-level claim that normal negative coupling collapses")
            rpt("in tumor.")
        else:
            rpt("**Claim language (E10):** Intermediate — some methods agree, clean-bg preserved.")
            rpt("Claim strength is moderate.")

        rpt()
        rpt("**Batch ceiling (E11):** All claims remain hedged by the perfect tissue-batch")
        rpt("confound. External replication (Lee et al.) is mandatory for any strong")
        rpt("biological interpretation. Tissue and batch are not separable.")
        rpt()
        rpt("**Project outcome pointer:** Pending Gate 5 PASS, evidence points toward")
        # NOTE: the output below is SUPERSEDED by Gate 3 and the post-audit addendum; see FINAL_REPORT.md
        rpt("Outcome 1 (proceed to secondaries) or Outcome 2 (reframe pending Gate 3),")
        rpt("depending on Gate 3 decomposition results.")

    elif outcome == "DIRECTION CONTRADICTED":
        rpt("The primary effect is not in the hypothesized direction.")
        rpt("**Project outcome:** Outcome 3 — abandon central narrative; archive as")
        rpt("methods/exploratory portfolio project.")

    elif outcome in ("NOT ESTABLISHED", "LOO-UNSTABLE"):
        rpt(f"The primary test does not meet all criteria ({outcome}).")
        rpt("**Project outcome:** Outcome 3 — central narrative not established;")
        rpt("archive as methods/exploratory portfolio project.")

    elif outcome == "POWER-LIMITED":
        rpt("Insufficient retained patients for a powered test.")
        rpt("**Project outcome:** UNDERPOWERED / INCONCLUSIVE (distinct from 'core dies').")

    rpt()

    # ══════════════════════════════════════════════════════════════════
    # Save report + provenance
    # ══════════════════════════════════════════════════════════════════
    _save_report(report, RESULTS_DIR / "gate4_donor_aware_report.md")

    prov = {
        "script": "04_gate4_donor_aware.py",
        "gate": "Gate 4",
        "versions": versions,
        "seed": SEED,
        "n_bootstrap": N_BOOT,
        "n_eligible": n_eligible,
        "eligible_pids": eligible_pids,
        "drift_check": "PASSED" if drift_match else "FAILED",
        "primary_outcome": outcome,
        "primary_mean_delta_z": primary["mean_delta_z"],
        "primary_wilcoxon_p": primary["wilcoxon_p"],
        "primary_bootstrap_ci": [primary["bootstrap_ci_lo"], primary["bootstrap_ci_hi"]],
        "primary_loo_all_same_sign": primary["loo_all_same_sign"],
        "primary_loo_max_rel_change": primary["loo_max_rel_change_pct"],
        "clean_bg_d2": d2_class if 'd2_class' in dir() else "N/A",
        "fit_kwargs_mixed": str(fit_kwargs),
        "inputs": {
            "gate2_scores": str(RESULTS_DIR / "gate2_recomputed_scores_cE01.parquet"),
            "gate2b_scores": str(RESULTS_DIR / "gate2b_clean_background_scores_cE01.parquet"),
            "config": str(CFG_PATH),
        },
        "status": "COMPLETE",
    }
    with open(PROV_DIR / "04_gate4_donor_aware.json", "w") as f:
        json.dump(prov, f, indent=2, default=str)

    print(f"\nSaved: provenance/04_gate4_donor_aware.json", flush=True)
    print("\n" + "=" * 72, flush=True)
    print(f"GATE 4 COMPLETE — OUTCOME: {outcome}", flush=True)
    print("=" * 72, flush=True)


def _save_report(lines, path):
    with open(path, "w") as f:
        f.write("\n".join(lines))
    print(f"Saved: {path.name}", flush=True)


if __name__ == "__main__":
    main()
