# ROBUSTNESS_LEDGER

Date: 2026-06-20  
Auditor: independent re-derivation  
Scope: Gate 4 robustness/sensitivity battery only.  
Mode: read-only extraction from committed artifacts; no analysis scripts re-run.

Current commit hash for tracked source artifacts: `1e4ca1005ef0008cf3e0ef1b061d91ff01223c7e`

Tracked-source check:

> `git ls-files --stage`: `100644 2093cb64d401d49903ba33e657430a52e6419985 0	core_validation_v2/results/tables/gate4_method_sensitivity.csv`

> `git ls-files --stage`: `100644 4d06b59d9882d44c01e49642e386ae1ca6e4a1af 0	core_validation_v2/results/tables/gate4_ivw_summary.csv`

> `git ls-files --stage`: `100644 593ac7d398f85bdcb59d2d66d975c90875d672b6 0	core_validation_v2/results/tables/gate4_mixed_models.csv`

> `git ls-files --stage`: `100644 f474ce9885cbd8281217b2c97f39caa145e8372d 0	core_validation_v2/results/tables/gate4_qc_primary_summary.csv`

> `git status --short` on the listed Gate 4 robustness artifacts returned no output, indicating no working-tree modifications to those tracked files.

## Extracted Robustness Battery

| Name | What it tested | Pre-committed criterion / interpretation rule | Result, verbatim with key numbers | Directionally consistent with primary Gate 4 effect? | Source file + field/line + commit hash |
|---|---|---|---|---|---|
| Method sensitivity across scoring methods | Donor-aware Fisher-z effect across official scanpy, clean-bg, AUCell, and z-score score definitions. | A001 A6 required scanpy as primary and AUCell/z-score as mandatory sensitivities; direction consistency required for sensitivities. A005 E10 required same direction across scanpy, clean-bg, AUCell, and z-score for "shift toward a more positive ISC-FAO relationship..." language. | CSV rows: `scanpy,34,0,0.400629,0.422797,0.0,0.337142,0.46414,True,0.9706,True,3.61`; `clean_bg,34,0,0.386688,0.393259,0.0,0.321472,0.450786,True,0.9706,True,3.72`; `aucell,34,0,0.151983,0.150491,3e-06,0.09984,0.202556,True,0.8529,True,8.11`; `zscore,34,0,0.206823,0.220008,0.0,0.15404,0.258007,True,0.9118,True,6.27`. Gate 4 report summary: `All methods agree in direction (positive Delta_z = more positive coupling in tumor).` | Yes. All mean Delta_z values are positive: scanpy +0.400629, clean-bg +0.386688, AUCell +0.151983, z-score +0.206823. AUCell/z-score magnitudes are smaller than scanpy but not reversed. | `core_validation_v2/results/tables/gate4_method_sensitivity.csv:1-5`; `core_validation_v2/results/gate4_donor_aware_report.md:77-104`; commit `1e4ca1005ef0008cf3e0ef1b061d91ff01223c7e`. |
| Inverse-variance weighting (IVW) | Descriptive inverse-variance-weighted Delta_z summary using approximate cell-count-derived SE weights. | A005 E7: IVW is descriptive and non-primary; report IVW point estimate and patient-level bootstrap CI; do not use IVW normal-theory p-value as a gate. | CSV row: `scanpy,0.386054,0.312504,0.460438,34,42,weights_approximate_cells_not_independent`. Gate 4 report: `IVW mean Delta_z: 0.386054`; `IVW bootstrap 95% CI: [0.312504, 0.460438]`; `Note: weights approximate (cells not independent replicates)`. | Yes. IVW mean Delta_z is positive at +0.386054 and the bootstrap CI is fully positive [0.312504, 0.460438]. | `core_validation_v2/results/tables/gate4_ivw_summary.csv:1-2`; `core_validation_v2/results/gate4_donor_aware_report.md:106-112`; commit `1e4ca1005ef0008cf3e0ef1b061d91ff01223c7e`. |
| Mixed-effects model | Corroborating cell-level mixed-effects interaction model, preferred random slope where converged, with patient grouping and ISC-by-tissue interaction. | A001 A7: Fisher-z remains primary; interaction agreeing in direction makes core claim more credible; reversal/inconsistency makes result method-sensitive; convergence failure is uninformative and does not override primary. A005 E8: random-slope preferred; batch not identifiable; positive ISC:tissue[T] interaction means more positive ISC-FAO slope in tumor; model does not resolve tissue-batch confounding or override per-patient primary. | CSV row: `random_slope,"FAO_scanpy ~ ISC_scanpy * C(SPECIMEN_TYPE, Treatment(reference='N'))",PID,~ISC_scanpy,"{'reml': False, 'method': 'lbfgs', 'maxiter': 500, 'disp': False}",41244,34,True,False,lbfgs,['The MLE may be on the boundary of the parameter space.'],"ISC_scanpy:C(SPECIMEN_TYPE, Treatment(reference='N'))[T.T]",0.36148427087348006,0.010286799956082967,0.34132251344338926,0.38164602830357086,1.6181942597910495e-270,CORROBORATIVE`. Gate 4 report: `Interaction ISC_scanpy:tissue[T]: coef=0.361484, SE=0.010287, p=0.000000, CI=[0.341323, 0.381646]`; `Interpretation level: CORROBORATIVE`; `Random-slope model converged and non-singular; fallback not needed.` | Yes. Interaction coefficient is positive at +0.361484, with CI [0.341323, 0.381646], matching the primary direction: more positive ISC-FAO slope in tumor. Caveat: model warning says MLE may be on boundary, and report states this does not resolve batch confounding. | `core_validation_v2/results/tables/gate4_mixed_models.csv:1-2`; `core_validation_v2/results/gate4_donor_aware_report.md:114-131`; commit `1e4ca1005ef0008cf3e0ef1b061d91ff01223c7e`. |
| QC residualization | Donor-aware effect after residualizing both scores within PID x tissue strata against available QC/technical covariates. | A001 A8: direction and reasonable magnitude retained weakens simple QC confound explanation; disappearance must be reported as QC/proliferation-sensitive; absent cell-cycle scores require QC-only residualization and the gap stated. A005 E9: complete-case stratum rules; rank-deficient or residual df < 10 => QC-UNINFORMATIVE; paired QC sensitivity requires both tissues informative; if fewer than 15 paired patients remain, POWER-LIMITED / UNINFORMATIVE. | CSV row: `68,68,0,0,34,0.396429,0.404918,0.0,0.333568,0.459614,COMPUTED`. Gate 4 report: `Total strata: 68`; `Informative strata: 68`; `Rank-deficient: 0`; `Low df / small n: 0`; `Paired patients (both tissues informative): 34`; `QC-residualized donor-aware: mean_Δz=0.396429, median_Δz=0.404918`; `Wilcoxon p=0.000000, CI=[0.333568, 0.459614]`; `Weakening here does not by itself prove an artifact.` | Yes. QC-residualized mean Delta_z is positive at +0.396429, close to the primary +0.400629, with fully positive CI [0.333568, 0.459614]. It was not inconclusive by the frozen stratum rules. | `core_validation_v2/results/tables/gate4_qc_primary_summary.csv:1-2`; `core_validation_v2/results/gate4_donor_aware_report.md:133-145`; commit `1e4ca1005ef0008cf3e0ef1b061d91ff01223c7e`. |

## Pre-Registration Source Quotes

### Method Sensitivity

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:41`: `### A6. Mandatory method-sensitivity across scoring methods (EXTENDS Gate 4)`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:42`: `Gate 2 already recomputes scanpy, AUCell, and mean z-score on the same log-normalized matrix. Gate 4's donor-aware analysis is run with all three: scanpy score_genes = primary; AUCell and mean z-score = mandatory sensitivities. The hard pass-criteria (A2, A3) are required only for the primary; for the two sensitivities, DIRECTION consistency is required (magnitude agreement is not, as methods are not numerically identical). If a sensitivity reverses direction, the result is labeled method-sensitive and no strong program-level claim is made.`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_005.md:42`: `## E10. Score-method claim-language matrix`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_005.md:44`: `- scanpy, clean-bg, AUCell, z-score all same direction -> "a shift toward a more positive ISC-FAO relationship in tumor-derived cE01 cells across scoring methods".`

### IVW

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_005.md:33`: `## E7. IVW is descriptive and non-primary`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_005.md:34`: `The inverse-variance-weighted Delta_z summary (weights w = 1/SE^2, SE = sqrt(1/(n_tumor-3)+1/(n_normal-3))) is descriptive only. Weights are approximate (E0). Report the IVW point estimate and a patient-level bootstrap CI. Each IVW bootstrap iteration resamples whole paired patient records with replacement, carrying each record Delta_z and its precomputed approximate weight together, and recomputes the weighted mean; never resample cells or tissues. Do NOT use an IVW normal-theory p-value as a gate.`

### Mixed-Effects Model

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:44`: `### A7. Gate 4 corroboration: mixed-effects interaction model (EXTENDS Gate 4)`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:45`: `The per-patient Fisher-z test remains the PRIMARY donor-aware gate (most conservative against pseudoreplication; patient is the unit). As mandatory corroboration on the same eligible patients, fit:`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:50`: `Interpretation (does NOT create a separate story):`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:51`: `- Fisher-z gate and the interaction agree in direction -> core claim more credible.`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:52`: `- Fisher-z gate passes but the interaction reverses or is clearly inconsistent -> result labeled method-sensitive; no strong program-level claim.`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_005.md:36`: `## E8. Mixed-effects corroboration (random-slope preferred; batch not identifiable)`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_005.md:37`: `Conceptually a patient random intercept plus a patient random slope on ISC, with an ISC-by-tissue fixed interaction. Implement in statsmodels (NOT lme4 formula syntax): smf.mixedlm("FAO_scanpy ~ ISC_scanpy * C(SPECIMEN_TYPE, Treatment(reference='N'))", data=eligible_cells, groups=eligible_cells["PID"], re_formula="~ISC_scanpy") as preferred if it converges; otherwise re_formula="1" labeled SUGGESTIVE ONLY. N is the explicit tissue reference, so a positive ISC:tissue[T] interaction means a more positive ISC-FAO slope in tumor. Do NOT add batchID (tissue and batch are perfectly confounded and not separable). Report formula, groups, re_formula, optimizer, n obs, n PIDs, convergence, singular-fit, warnings, and the interaction coef/SE/CI/p. Neither model resolves tissue-batch confounding nor overrides the per-patient primary.`

### QC Residualization

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:56`: `### A8. QC / technical confound sensitivity (EXTENDS Gate 4)`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:57`: `Gate 1 reports presence and real names of total_counts, n_genes_by_counts, pct_counts_mt, and cell-cycle scores (S_score, G2M_score), discovering actual names if different. The core donor-aware analysis is then repeated after residualizing BOTH scores per cell against available pre-specified QC covariates:`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:60`: `- Direction and a reasonable magnitude retained after residualization -> a simple technical/QC confound explanation is weakened.`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:61`: `- Effect largely disappears -> reported explicitly as QC/proliferation-sensitive.`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_001.md:63`: `- If cell-cycle scores are absent, a QC-only residualization is reported and the gap is stated.`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_005.md:39`: `## E9. QC residualization failed-stratum rules (frozen)`

> `CORE_VALIDATION_PROTOCOL_AMENDMENT_005.md:40`: `Per PID x tissue stratum: complete cases only, intercept, regress each score on log1p(total_counts), n_genes_by_counts, pct_counts_mt. Report n before/after, design-matrix rank, condition number. If rank-deficient OR residual df < 10, mark QC-UNINFORMATIVE. No silent covariate drop. Paired QC sensitivity includes only patients with BOTH tissues informative. If fewer than 15 paired patients remain, label POWER-LIMITED / UNINFORMATIVE. Weakening under QC residualization does not by itself prove an artifact.`

## Verbatim Result Quotes

### Method Sensitivity

> `core_validation_v2/results/tables/gate4_method_sensitivity.csv:1`: `method,n_retained,n_excluded,mean_delta_z,median_delta_z,wilcoxon_p,bootstrap_ci_lo,bootstrap_ci_hi,ci_excludes_0,proportion_delta_z_positive,loo_all_same_sign,loo_max_rel_change_pct`

> `core_validation_v2/results/tables/gate4_method_sensitivity.csv:2`: `scanpy,34,0,0.400629,0.422797,0.0,0.337142,0.46414,True,0.9706,True,3.61`

> `core_validation_v2/results/tables/gate4_method_sensitivity.csv:3`: `clean_bg,34,0,0.386688,0.393259,0.0,0.321472,0.450786,True,0.9706,True,3.72`

> `core_validation_v2/results/tables/gate4_method_sensitivity.csv:4`: `aucell,34,0,0.151983,0.150491,3e-06,0.09984,0.202556,True,0.8529,True,8.11`

> `core_validation_v2/results/tables/gate4_method_sensitivity.csv:5`: `zscore,34,0,0.206823,0.220008,0.0,0.15404,0.258007,True,0.9118,True,6.27`

> `core_validation_v2/results/gate4_donor_aware_report.md:99-104`: `scanpy: + (mean_Δz = 0.400629)` / `clean_bg: + (mean_Δz = 0.386688)` / `aucell: + (mean_Δz = 0.151983)` / `zscore: + (mean_Δz = 0.206823)` / `All methods agree in direction (positive Delta_z = more positive coupling in tumor). Strong cross-method consistency supports the program-level claim language.`

Important post-Gate-3 caveat: the final clause "supports the program-level claim language" is superseded by later Gate 3 interpretation. The robustness result remains directionally consistent, but it does not support a coherent FAO/beta-oxidation program claim.

> `core_validation_v2/results/gate4_donor_aware_report.md:203`: `The score-level donor-aware shift survives (Gate 4 + Gate 5 PASS), but the coherent FAO/beta-oxidation program-level claim does NOT.`

### IVW

> `core_validation_v2/results/tables/gate4_ivw_summary.csv:1`: `method,ivw_mean_delta_z,bootstrap_ci_low,bootstrap_ci_high,n_patients,seed,note`

> `core_validation_v2/results/tables/gate4_ivw_summary.csv:2`: `scanpy,0.386054,0.312504,0.460438,34,42,weights_approximate_cells_not_independent`

> `core_validation_v2/results/gate4_donor_aware_report.md:108`: `IVW mean Delta_z: 0.386054`

> `core_validation_v2/results/gate4_donor_aware_report.md:109`: `IVW bootstrap 95% CI: [0.312504, 0.460438]`

> `core_validation_v2/results/gate4_donor_aware_report.md:110`: `Note: weights approximate (cells not independent replicates)`

### Mixed-Effects Model

> `core_validation_v2/results/tables/gate4_mixed_models.csv:1`: `model_type,formula,groups,re_formula,fit_kwargs,n_cells,n_patients,converged,singular_fit,optimizer,warnings,interaction_term,interaction_coef,interaction_se,interaction_ci_low,interaction_ci_high,interaction_p,interpretation_level`

> `core_validation_v2/results/tables/gate4_mixed_models.csv:2`: `random_slope,"FAO_scanpy ~ ISC_scanpy * C(SPECIMEN_TYPE, Treatment(reference='N'))",PID,~ISC_scanpy,"{'reml': False, 'method': 'lbfgs', 'maxiter': 500, 'disp': False}",41244,34,True,False,lbfgs,['The MLE may be on the boundary of the parameter space.'],"ISC_scanpy:C(SPECIMEN_TYPE, Treatment(reference='N'))[T.T]",0.36148427087348006,0.010286799956082967,0.34132251344338926,0.38164602830357086,1.6181942597910495e-270,CORROBORATIVE`

> `core_validation_v2/results/gate4_donor_aware_report.md:120-125`: `Converged: True` / `Singular: False` / `Warnings: ['The MLE may be on the boundary of the parameter space.']` / `Interaction ISC_scanpy:tissue[T]: coef=0.361484, SE=0.010287, p=0.000000, CI=[0.341323, 0.381646]` / `Interpretation: positive coef = more positive ISC-FAO slope in tumor` / `Interpretation level: CORROBORATIVE`

> `core_validation_v2/results/gate4_donor_aware_report.md:130-131`: `**Note:** This does not resolve tissue-batch confounding (perfectly confounded,` / `Cramer's V = 1.0) and does not override the per-patient primary.`

### QC Residualization

> `core_validation_v2/results/tables/gate4_qc_primary_summary.csv:1`: `n_strata_total,n_strata_informative,n_strata_rank_deficient,n_strata_low_df,n_paired_patients_retained,mean_delta_z,median_delta_z,wilcoxon_p,bootstrap_ci_low,bootstrap_ci_high,qc_sensitivity_status`

> `core_validation_v2/results/tables/gate4_qc_primary_summary.csv:2`: `68,68,0,0,34,0.396429,0.404918,0.0,0.333568,0.459614,COMPUTED`

> `core_validation_v2/results/gate4_donor_aware_report.md:135-143`: `Total strata: 68` / `Informative strata: 68` / `Rank-deficient: 0` / `Low df / small n: 0` / `Paired patients (both tissues informative): 34` / `QC-residualized donor-aware: mean_Δz=0.396429, median_Δz=0.404918` / `Wilcoxon p=0.000000, CI=[0.333568, 0.459614]` / `Weakening here does not by itself prove an artifact.`

## Cell-Cycle Caveat Confirmation

The requested cell-cycle caveat is present in the Gate 4 report.

> `core_validation_v2/results/gate4_donor_aware_report.md:163`: `**Cell-cycle limitation (C6):** S/G2M cell-cycle scores are absent from`

> `core_validation_v2/results/gate4_donor_aware_report.md:164`: `the metadata; residual cell-cycle variation within cE01 cannot be fully`

> `core_validation_v2/results/gate4_donor_aware_report.md:165`: `assessed. Primary analysis excludes Stem/TA-like_prolif.`

Plain-language confirmation: QC residualization was performed without S/G2M cell-cycle covariates because they were absent. The effect remained directionally positive after available-QC residualization, but residual proliferation/cell-cycle variation within cE01 cannot be excluded.

## Consistency Summary

All four executed Gate 4 robustness checks were directionally consistent with the primary Gate 4 effect:

- Method sensitivity: all methods positive; mean Delta_z scanpy +0.400629, clean-bg +0.386688, AUCell +0.151983, z-score +0.206823.
- IVW: positive IVW mean Delta_z +0.386054, CI [0.312504, 0.460438].
- Mixed-effects: positive random-slope interaction +0.361484, CI [0.341323, 0.381646], corroborative, with boundary warning and batch-confound caveat.
- QC residualization: positive mean Delta_z +0.396429, CI [0.333568, 0.459614], computed/informative across 34 paired patients.

None of the four was directionally inconsistent. None was classified as power-limited or inconclusive. The mixed model carries a boundary warning and cannot resolve the perfect tissue-batch confound; IVW is descriptive because weights are approximate; QC residualization cannot address absent cell-cycle covariates.
