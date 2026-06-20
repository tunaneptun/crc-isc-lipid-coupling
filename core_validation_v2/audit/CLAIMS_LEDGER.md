# CLAIMS LEDGER

**Date:** 2026-06-19
**Method:** Read-only lookup over committed artifacts. No values re-derived or computed.

---

## Token Table

| Token | Intended meaning | Value (verbatim) | Source file | Field/line or cell | Commit hash | Notes |
|-------|-----------------|-------------------|-------------|-------------------|-------------|-------|
| **DATA / §3** | | | | | | |
| `epi_n` | Total epithelial cells in the Pelka scoring population | 168295 | `core_validation_v2/results/gate1_input_validation.json` | `.checks.dimensions.norm_cells` | `cdd0ef22798654505ee82efcd04acbee4fc4e0c4` | Also in gate1 report line 24. This is the full `pelka_epithelial.h5ad`, not just cE01. |
| `cE01_n` | Number of cE01 (Stem/TA-like) cells | 61953 | `core_validation_v2/results/gate1_input_validation.json` | `.checks.ce01.ce01_count` | `cdd0ef22798654505ee82efcd04acbee4fc4e0c4` | Agrees with gate1 report line 107. |
| `cE01_tumor_n` | cE01 cells in tumor tissue | 50167 | `core_validation_v2/results/tables/gate2_pooled_decoupling.csv` | `n_tumor` column, row "scanpy (recomputed)" | `04adaaed5f76358781d066af5692e11a98c8cd64` | Also in gate1 report line 109. Not a named field in gate1 JSON. |
| `cE01_normal_n` | cE01 cells in adjacent-normal tissue | 11786 | `core_validation_v2/results/tables/gate2_pooled_decoupling.csv` | `n_normal` column, row "scanpy (recomputed)" | `04adaaed5f76358781d066af5692e11a98c8cd64` | Also in gate1 report line 110. Not a named field in gate1 JSON. |
| `n_patients_ge30` | Eligible paired donors at >=30 threshold | 34 | `core_validation_v2/results/gate1_input_validation.json` | `.checks.eligibility.thresh_30.n_eligible` | `cdd0ef22798654505ee82efcd04acbee4fc4e0c4` | Primary set for all gates. |
| `elig_thresholds` | Per-tissue thresholds and patient counts | >=20: 35; >=30: 34; >=50: 33 | `core_validation_v2/results/gate1_input_validation.json` | `.checks.eligibility.thresh_20/30/50.n_eligible` | `cdd0ef22798654505ee82efcd04acbee4fc4e0c4` | |
| `cramers_v` | Cramer's V (tissue x batch) over eligible set | 1.000000 | `core_validation_v2/redteam/RED_TEAM_GATE1.md` | Line 63: "cE01 tissue-vs-batch Cramer's V: 1.000000" | `b8339cb95330aa8e25e9df3da3314f5fbc5df3e3` | Gate1 JSON records `confound_level: STRONG`, `n_same_batch: 0` but not the V value itself. Redteam report is the primary committed source of the numeric V. |
| `shared_batch` | Patients sharing a T/N batch | 0 / 34 | `core_validation_v2/results/gate1_input_validation.json` + `core_validation_v2/redteam/RED_TEAM_GATE1.md` | JSON: `.checks.batch_tissue_confound.n_same_batch` = 0, `.n_paired_patients` = 35 (>=20 set); Redteam line 27: "same-batch T/N=0" over 34 eligible | `cdd0ef2...` / `b8339cb...` | Gate1 JSON denominator is 35 (>=20 set); redteam and all reports use 0/34 (>=30 primary set). The correct report-level token is **0 / 34** per the primary eligible set. |
| **GATE 2 / 2b / §5** | | | | | | |
| `r_normal_pooled` | Gate 2 pooled Pearson r, normal, scanpy | -0.380809 | `core_validation_v2/results/tables/gate2_pooled_decoupling.csv` | `r_normal`, row "scanpy (recomputed)" | `04adaaed5f76358781d066af5692e11a98c8cd64` | |
| `r_tumor_pooled` | Gate 2 pooled Pearson r, tumor, scanpy | -0.018209 | `core_validation_v2/results/tables/gate2_pooled_decoupling.csv` | `r_tumor`, row "scanpy (recomputed)" | `04adaaed5f76358781d066af5692e11a98c8cd64` | |
| `r_normal_cleanbg` | Gate 2b pooled Pearson r, normal, clean-bg | -0.382179 | `core_validation_v2/results/tables/gate2b_pooled_cleanbg.csv` | `r_normal`, row "clean_bg_scanpy" | `16678de8ac85408a32ed3c876f111c81aadf9460` | This is the committed Pelka normal anchor. Agrees with Lee report guard 6a. |
| `r_tumor_cleanbg` | Gate 2b pooled Pearson r, tumor, clean-bg | -0.032731 | `core_validation_v2/results/tables/gate2b_pooled_cleanbg.csv` | `r_tumor`, row "clean_bg_scanpy" | `16678de8ac85408a32ed3c876f111c81aadf9460` | |
| **GATE 4** | | | | | | |
| `gate4_mean_dz` | Mean per-patient Delta_z | 0.400629 | `core_validation_v2/results/tables/gate4_primary_summary.csv` | `mean_delta_z` | `71b3e9eb1d09d76a213027cafaa96b39ad7a63f0` | |
| `gate4_median_dz` | Median per-patient Delta_z | 0.422797 | `core_validation_v2/results/tables/gate4_primary_summary.csv` | `median_delta_z` | `71b3e9eb1d09d76a213027cafaa96b39ad7a63f0` | |
| `gate4_wilcoxon_p` | Wilcoxon signed-rank p (two-sided) | 2.3283064365386963e-10 | `core_validation_v2/results/tables/gate4_primary_summary.csv` | `wilcoxon_p` | `71b3e9eb1d09d76a213027cafaa96b39ad7a63f0` | Gate4 report displays "0.000000" (format truncation). Full-precision value here is the authoritative source. |
| `gate4_ci_lo` | Bootstrap 95% CI lower bound | 0.3371415699522461 | `core_validation_v2/results/tables/gate4_primary_summary.csv` | `bootstrap_ci_lo` | `71b3e9eb1d09d76a213027cafaa96b39ad7a63f0` | |
| `gate4_ci_hi` | Bootstrap 95% CI upper bound | 0.46414002142135113 | `core_validation_v2/results/tables/gate4_primary_summary.csv` | `bootstrap_ci_hi` | `71b3e9eb1d09d76a213027cafaa96b39ad7a63f0` | |
| `gate4_prop_pos` | Proportion and count of patients with positive Delta_z | 0.9706 (33/34) | `core_validation_v2/results/tables/gate4_primary_summary.csv` + `gate4_per_patient.csv` | `prop_delta_z_positive` = 0.9706; per-patient CSV shows 33 positive, 1 negative | `71b3e9eb...` | 33/34 derived from per-patient CSV. |
| `gate4_neg_patient` | Single patient with negative Delta_z | C140 | `core_validation_v2/results/tables/gate4_per_patient.csv` | method=scanpy, PID=C140, Delta_z=-0.076247 | `71b3e9eb1d09d76a213027cafaa96b39ad7a63f0` | Only negative Delta_z in the scanpy rows. |
| `gate4_loo_outcome` | LOO stability: all same sign + max % change | All 34 LOO means positive; max relative change = 3.61% (C140) | `core_validation_v2/results/tables/gate4_loo_primary.csv` | All `sign_positive` = True; max `relative_effect_change_pct` = 3.61 at PID C140 | `71b3e9eb1d09d76a213027cafaa96b39ad7a63f0` | |
| **GATE 5** | | | | | | |
| `gate5_dir_retention` | Equal-cell direction retention | 1.0 | `core_validation_v2/results/tables/gate5_primary_summary.csv` | `direction_retention_fraction` | `9b293500cc5367077388bd205d23a469f663faee` | 1000/1000 draws positive. |
| `gate5_median_mag` | Median magnitude across draws | 0.414199 | `core_validation_v2/results/tables/gate5_primary_summary.csv` | `median_downsampled_magnitude` | `9b293500cc5367077388bd205d23a469f663faee` | |
| `gate5_n_draws` | Number of equal-cell draws | 1000 | `core_validation_v2/results/tables/gate5_primary_summary.csv` | `K` | `9b293500cc5367077388bd205d23a469f663faee` | |
| **GATE 3** | | | | | | |
| `gate3_baseline_dz` | Clean fixed-bg baseline mean Delta_z | 0.386688 | `core_validation_v2/results/tables/gate3_variant_effects.csv` | row variant=baseline, col `mean_delta_z` | `c4c5ce32b9ae901b49f37a8f213abf05e39bd5f2` | |
| `gate3_fabp1_retention` | Retention with FABP1 removed | 0.4303 | `core_validation_v2/results/tables/gate3_variant_effects.csv` | row variant=fabp1_removed, col `retention_ratio` | `c4c5ce32b9ae901b49f37a8f213abf05e39bd5f2` | Materially weakened (<0.50). |
| `gate3_corebox_retention` | Core beta-ox-5 retention | 0.3499 | `core_validation_v2/results/tables/gate3_variant_effects.csv` | row variant=core_beta_oxidation_5, col `retention_ratio` | `c4c5ce32b9ae901b49f37a8f213abf05e39bd5f2` | Materially weakened (<0.50). |
| `gate3_hmgcs2_retention` | Retention with HMGCS2 removed | 1.0069 | `core_validation_v2/results/tables/gate3_variant_effects.csv` | row variant=hmgcs2_removed, col `retention_ratio` | `c4c5ce32b9ae901b49f37a8f213abf05e39bd5f2` | Robust (>0.50). |
| `gate3_olfm4_retention` | ISC-LOO retention with OLFM4 removed | 0.5018 | `core_validation_v2/results/tables/gate3_variant_effects.csv` | row variant=isc_loo_OLFM4, col `retention_ratio` | `c4c5ce32b9ae901b49f37a8f213abf05e39bd5f2` | 0.0018 above 0.50 threshold; marginal. |
| `gate3_booleans` | Three frozen booleans | `core_beta_oxidation_holds` = FALSE; `fao_reframe_conjunction_triggered` = FALSE; `program_level_fao_robust` = FALSE | `INTERPRETATION_NOTE_GATE3.md` | Lines 29–31 | `c4c5ce32b9ae901b49f37a8f213abf05e39bd5f2` | **See notes below.** The provenance JSON (`04_gate3_score_decomposition.json`) stores *equivalent* booleans under different keys: `fao_reframe_needed`=false, `isc_single_gene_driven`=false, `fao_single_gene_driven`=true. The canonical three-boolean naming appears only in INTERPRETATION_NOTE_GATE3.md. |
| **GATE 6** | | | | | | |
| `gate6_lgr5_n` | LGR5-detected cE01 cells | 9804 | `core_validation_v2/provenance/06_gate6_lgr5_sensitivity.json` | `n_lgr5_detected_ce01` | `462c0f1b20113bbaab2f43ef494213d550fde08d` | |
| `gate6_n_patients` | Eligible patients (>=15 both tissues) | 6 | `core_validation_v2/provenance/06_gate6_lgr5_sensitivity.json` | `n_eligible_pids` | `462c0f1b20113bbaab2f43ef494213d550fde08d` | PIDs: C122, C124, C125, C126, C139, C143. |
| `gate6_dz` | Descriptive mean Delta_z (LGR5 subset) | 0.315498124073651 | `core_validation_v2/provenance/06_gate6_lgr5_sensitivity.json` | `subset_stats.mean_delta_z` | `462c0f1b20113bbaab2f43ef494213d550fde08d` | Agrees with gate6 report "Mean Delta z (LGR5 subset): 0.315498". |
| **LEE** | | | | | | |
| `lee_r_pelka_normal` | Normal-anchor pooled r, Pelka | -0.382179 | `core_validation_v2/results/tables/lee_normal_anchor_correlations.csv` | row pair=primary, cohort=Pelka, col `r` | `c9caebcf9bfe5a3fd9d0fdad64af224b7d46a63b` | Matches `r_normal_cleanbg` (Gate 2b). |
| `lee_r_smc_normal` | Normal-anchor pooled r, SMC | -0.25692299008369446 | `core_validation_v2/results/tables/lee_normal_anchor_correlations.csv` | row pair=primary, cohort=SMC, col `r` | `c9caebcf9bfe5a3fd9d0fdad64af224b7d46a63b` | |
| `lee_r_kul3_normal` | Normal-anchor pooled r, KUL3 | -0.4066520035266876 | `core_validation_v2/results/tables/lee_normal_anchor_correlations.csv` | row pair=primary, cohort=KUL3, col `r` | `c9caebcf9bfe5a3fd9d0fdad64af224b7d46a63b` | |
| `lee_tally` | Per-patient sign tally (neg / informative) | 12 / 12 (all informative patients negative) | `core_validation_v2/results/tables/lee_per_patient_sign_tally.csv` | 16 total rows; 12 with n>=30, all sign="neg"; 4 with n<30 excluded | `c9caebcf9bfe5a3fd9d0fdad64af224b7d46a63b` | |
| `lee_control_counts` | Per-cohort recovered controls (ISC / FAO) | SMC: 250 / 300; KUL3: 350 / 300; Pelka: 300 / 400 | `core_validation_v2/provenance/lee_external_consistency.json` | `lee_controls.SMC.isc_n`/`fao_n`, `.KUL3.isc_n`/`fao_n`, `pelka_controls.isc_n`/`fao_n` | `c9caebcf9bfe5a3fd9d0fdad64af224b7d46a63b` | |
| **ENVIRONMENT / §9** | | | | | | |
| `env_versions` | Package versions used | scanpy 1.11.5, numpy 2.4.2, pandas 2.3.3, scipy 1.17.1, anndata 0.12.10 | `core_validation_v2/provenance/lee_external_consistency.json` (scanpy/numpy/pandas/scipy) + `BASELINE_ENVIRONMENT.txt` (anndata) | `versions` dict in JSON; pip freeze in env file | `c9caebcf...` / `fcfd5e99...` | scanpy absent from gate4 JSON; present in Lee JSON. statsmodels 0.14.6 also used (gate4 JSON) but not a core token. |
| **GENE SETS / appendix** | | | | | | |
| `geneset_isc9` | 9 ISC genes | LGR5, OLFM4, ASCL2, SOX9, CD44, SMOC2, RGMB, EPHB2, MSI1 | `core_validation_v2/config.yaml` | `gene_sets.isc_9` | `aa9ada3ac3888a1379223568975555d98e902aa2` | |
| `geneset_fao12` | 12 lipid/FAO-composite genes | CPT1A, HMGCS2, FABP1, ACOX1, ACADL, ACADM, PDK4, PPARD, PPARA, ANGPTL4, HADH, ACAA2 | `core_validation_v2/config.yaml` | `gene_sets.ppar_fao_12` | `aa9ada3ac3888a1379223568975555d98e902aa2` | |
| `geneset_corebox5` | 5 core beta-oxidation genes | CPT1A, ACADL, ACADM, HADH, ACAA2 | `core_validation_v2/config.yaml` | `gene_sets.core_beta_oxidation_5` | `aa9ada3ac3888a1379223568975555d98e902aa2` | |
| `geneset_fao11` | FAO-12 minus FABP1 | CPT1A, HMGCS2, ACOX1, ACADL, ACADM, PDK4, PPARD, PPARA, ANGPTL4, HADH, ACAA2 | `core_validation_v2/config.yaml` | Derived from `gene_sets.ppar_fao_12` minus FABP1 | `aa9ada3ac3888a1379223568975555d98e902aa2` | Also explicitly constructed in `08_lee_external_consistency.py` line 33. |
| `geneset_isc8` | ISC-9 minus LGR5 | OLFM4, ASCL2, SOX9, CD44, SMOC2, RGMB, EPHB2, MSI1 | `core_validation_v2/config.yaml` | Derived from `gene_sets.isc_9` minus LGR5 | `aa9ada3ac3888a1379223568975555d98e902aa2` | Also explicitly constructed in `08_lee_external_consistency.py` line 31. |

---

## SPECIAL TOKEN: `r_normal_repr` / `r_tumor_repr`

### The question

The narrative framing uses "r_normal ~ −0.44 → r_tumor ~ −0.07" (INTERPRETATION_NOTE_GATE3.md line 53). These do NOT match the pooled Gate 2 values (-0.3808 / -0.0182) or the pooled Gate 2b values (-0.3822 / -0.0327).

### (a) Where do -0.44 / -0.07 appear in committed artifacts?

| File | Line | Exact text |
|------|------|-----------|
| `INTERPRETATION_NOTE_GATE3.md` | 53 | "r_normal ≈ −0.44 → r_tumor ≈ −0.07" |

This is the only place (the value traces to `gate4_primary_summary.csv`: mean_r_normal = −0.436952, mean_r_tumor = −0.072628).

### (b) Mean and median of per-patient r_normal and r_tumor (committed source)

From `core_validation_v2/results/tables/gate4_primary_summary.csv` (commit `71b3e9eb...`):

| Statistic | r_normal | r_tumor |
|-----------|----------|---------|
| `mean_r_normal` / `mean_r_tumor` | **-0.436952** | **-0.072628** |
| `median_r_normal` / `median_r_tumor` | **-0.434499** | **-0.086323** |

Also stated in `gate4_donor_aware_report.md` lines 49–52:
> "Mean raw r_normal: -0.436952 / Mean raw r_tumor: -0.072628 / Median raw r_normal: -0.434499 / Median raw r_tumor: -0.086323"

### (c) Assessment

| Narrative value | Closest committed quantity | Committed value | Match quality |
|----------------|---------------------------|-----------------|---------------|
| r_normal ≈ −0.44 | `mean_r_normal` | -0.436952 | Rounds to -0.44 ✓ |
| r_tumor ≈ −0.07 | `mean_r_tumor` | -0.072628 | Rounds to -0.07 ✓ |

**The narrative "~ −0.44 / ~ −0.07" corresponds to the MEAN of per-patient Pearson r values** (`mean_r_normal` = -0.436952 ≈ -0.44; `mean_r_tumor` = -0.072628 ≈ -0.07), NOT to the pooled cell-level correlations.

This is a meaningful distinction:
- **Pooled cell-level r** (Gate 2, standard scanpy): r_normal = -0.3808, r_tumor = -0.0182
- **Mean of per-patient r** (Gate 4 primary summary): mean_r_normal = -0.4370, mean_r_tumor = -0.0726

The per-patient mean r is ~15% more negative in normal than the pooled value because pooled correlation is dominated by high-cell-count patients while the mean of per-patient values weights each patient equally.

**Recommendation:** The write-up should use the committed full-precision values (`mean_r_normal` = -0.437, `mean_r_tumor` = -0.073) and explicitly label them as "mean of per-patient Pearson r" to distinguish from the pooled cell-level values. The current "≈ −0.44" is a legitimate rounding of -0.437 but omits the distinction.

**Status: RESOLVED** — the values correspond to committed quantities, but the report should cite the precise source.

---

## Counts

| Category | Count |
|----------|-------|
| **RESOLVED** | 38 |
| **UNRESOLVED** | 0 |
| **DISCREPANT** | 1 |

### DISCREPANT tokens

**`gate3_booleans`:** The three canonical boolean names (`core_beta_oxidation_holds`, `fao_reframe_conjunction_triggered`, `program_level_fao_robust`) with values (FALSE, FALSE, FALSE) appear in `INTERPRETATION_NOTE_GATE3.md` lines 29–31 (committed at `c4c5ce32...`). The provenance JSON (`04_gate3_score_decomposition.json`, same commit) stores *functionally equivalent* booleans under **different key names**: `fao_reframe_needed` = false, `isc_single_gene_driven` = false, `fao_single_gene_driven` = true. The mapping is:

| Canonical name (INTERPRETATION_NOTE) | JSON key | Canonical value | JSON value | Equivalent? |
|--------------------------------------|----------|-----------------|------------|-------------|
| `core_beta_oxidation_holds` | *(not present)* | FALSE | *(absent)* | No direct JSON counterpart |
| `fao_reframe_conjunction_triggered` | `fao_reframe_needed` | FALSE | false | Yes (same meaning) |
| `program_level_fao_robust` | `fao_single_gene_driven` | FALSE | true | Yes (inverted: single-gene-driven=true ↔ robust=false) |

The discrepancy is in **naming and key structure**, not in substance. The JSON does not store `core_beta_oxidation_holds` at all — it must be inferred from `gate3_variant_effects.csv` (variant=core_beta_oxidation_5, retention_ratio=0.3499 < 0.50 → false). The write-up should use the INTERPRETATION_NOTE canonical names and cite both sources.

### UNRESOLVED tokens

None.

---

## Robustness tokens (appended from ROBUSTNESS_LEDGER.md)

| Token | Intended meaning | Value (verbatim) | Source file | Field/line or cell | Commit hash | Notes |
|-------|-----------------|-------------------|-------------|-------------------|-------------|-------|
| `aucell_dz` | AUCell mean Delta_z | 0.151983 | `core_validation_v2/results/tables/gate4_method_sensitivity.csv` | row method=aucell, col `mean_delta_z` | `1e4ca1005ef0008cf3e0ef1b061d91ff01223c7e` | Smaller magnitude than scanpy but not reversed. |
| `zscore_dz` | z-score mean Delta_z | 0.206823 | `core_validation_v2/results/tables/gate4_method_sensitivity.csv` | row method=zscore, col `mean_delta_z` | `1e4ca1005ef0008cf3e0ef1b061d91ff01223c7e` | Smaller magnitude than scanpy but not reversed. |
| `ivw_dz` | IVW mean Delta_z | 0.386054 | `core_validation_v2/results/tables/gate4_ivw_summary.csv` | `ivw_mean_delta_z` | `1e4ca1005ef0008cf3e0ef1b061d91ff01223c7e` | Descriptive only; weights approximate. |
| `mixed_interaction` | Mixed-effects random-slope interaction coefficient | 0.361484 | `core_validation_v2/results/tables/gate4_mixed_models.csv` | `interaction_coef` (random_slope row) | `1e4ca1005ef0008cf3e0ef1b061d91ff01223c7e` | CORROBORATIVE; MLE boundary warning; does not resolve batch confound. Full precision in CSV: 0.36148427087348006. |
| `qc_dz` | QC-residualized mean Delta_z | 0.396429 | `core_validation_v2/results/tables/gate4_qc_primary_summary.csv` | `mean_delta_z` | `1e4ca1005ef0008cf3e0ef1b061d91ff01223c7e` | 34 paired patients retained; cell-cycle covariates absent. |
