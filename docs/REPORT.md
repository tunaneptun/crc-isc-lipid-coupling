# Stemness and lipid/FAO coupling in human colorectal-cancer Stem/TA-like cells

A donor-aware single-cell analysis in the Pelka 2021 colorectal-cancer atlas.

## 1. Summary

In mouse models, a high-fat diet expands intestinal stem and progenitor cells and raises tumorigenicity through a PPARδ/PPARα-associated fatty-acid-oxidation (FAO) and lipid-metabolism program. It is not established whether a comparable link between stemness and lipid metabolism is detectable in human colorectal cancer. This report asks one specific, testable version of that question in the Pelka 2021 colorectal-cancer single-cell atlas: within cE01 "Stem/TA-like" epithelial cells, does the cell-level relationship between an intestinal-stem-cell (ISC) score and a lipid/FAO-associated composite score differ between tumor and patient-matched adjacent-normal tissue, analysed with the patient as the unit of inference?

The answer is a bounded one. The relationship does shift: an inverse association in normal tissue becomes much weaker in tumor, and the shift survives donor-aware testing, equal-cell downsampling, and several alternative scoring choices. But it is not a coherent FAO or β-oxidation program (it rests largely on one gene, FABP1, which is not a β-oxidation enzyme), it is correlational only (tissue and sequencing batch are perfectly confounded), and it is not replicated on the tumor side in a second dataset. The result is best read as a hypothesis-generating transcriptional pattern, not a mechanism.

## 2. Background

High-fat diet, PPARδ/PPARα signalling, and fatty-acid oxidation have been linked to intestinal stem-cell function and tumorigenicity in mouse models, and reviewed as part of a broader relationship between diet, lipid metabolism, and intestinal stem cells [1, 2, 3]. Whether an analogous coupling between stemness and lipid-associated transcriptional state is visible in human colorectal-cancer epithelium is an open question. Most single-cell analyses of colorectal cancer focus on which cell types or pathways change in abundance or expression. This analysis instead asks whether the relationship between two transcriptional axes, stemness and lipid/FAO-associated expression, changes within a single epithelial state between tumor and matched normal tissue.

## 3. The question

Within cE01 Stem/TA-like epithelial cells of the Pelka atlas, does the per-patient correlation between an ISC score and a lipid/FAO-associated composite score differ between tumor-derived and adjacent-normal cells? The scope is deliberately narrow: one population, one pair of scores, one tumor-versus-normal contrast, examined through several independent checks rather than broadened into exploratory directions.

## 4. Data and the batch ceiling

The Pelka et al. 2021 atlas epithelial compartment contains 168,295 cells. The analysis population is cE01 "Stem/TA-like" cells (cl295v11SubShort == "cE01", matching epithelial_subtype == "Stem/TA-like"); the proliferating subset is excluded. cE01 totals 61,953 cells (50,167 tumor, 11,786 adjacent-normal). Patients require at least 30 cE01 cells in both tissues for the primary analysis, giving 34 paired donors (35 at a threshold of 20, 33 at 50).

One feature of this cohort bounds every result. Tissue type and sequencing batch are perfectly confounded: Cramér's V = 1.000000, and none of the 34 patients had tumor and normal tissue processed in the same batch. No statistical method can separate a biological tumor-versus-normal difference from a batch difference under this design. Every result below is therefore correlational, and no causal or mechanistic claim follows from it. This is stated before the results because it limits all of them.

## 5. Methods

Scores are computed on the log-normalized expression matrix with scanpy score_genes (use_raw=False, ctrl_size=50, n_bins=25, random_state=42). The gene sets are fixed: an ISC score of 9 genes (LGR5, OLFM4, ASCL2, SOX9, CD44, SMOC2, RGMB, EPHB2, MSI1), a lipid/FAO-associated composite of 12 genes (CPT1A, HMGCS2, FABP1, ACOX1, ACADL, ACADM, PDK4, PPARD, PPARA, ANGPTL4, HADH, ACAA2), and a 5-gene core β-oxidation sub-score (CPT1A, ACADL, ACADM, HADH, ACAA2).

Three correlation quantities are kept strictly separate and are not interchangeable: the pooled cell-level correlation (all cells together, which is pseudoreplicated and used only descriptively), the mean of per-patient correlations, and the donor-aware paired difference in Fisher-z-transformed correlations. Only the donor-aware quantity is used for inference, with the patient as the unit. Significance is assessed with a Wilcoxon signed-rank test over patients, and uncertainty with a patient-level bootstrap confidence interval.

Every number reported here was re-derived independently from the source data by a separate reconstruction that does not import the analysis code or trust its outputs, and the two derivations agree. Each analysis step writes a provenance record (input file hashes, parameters, seed, outputs, counts).

## 6. Results

**Pooled, descriptive.** Across all cells, the ISC-composite correlation is -0.382179 in normal and -0.032731 in tumor. Because cells within a patient are not independent, this pooled view is descriptive only and is not used for inference.

**Donor-aware (primary).** Treating each patient as one observation, the mean per-patient correlation is -0.437 in normal and -0.073 in tumor. The paired per-patient Fisher-z difference has a mean of 0.400629 (Wilcoxon p = 2.33 x 10⁻¹⁰; 95% bootstrap CI 0.337 to 0.464, excluding zero), and 33 of 34 patients show a positive shift. Leave-one-patient-out re-estimation leaves the mean essentially unchanged (largest single-patient change about 3.6%).

**Equal-cell downsampling.** Drawing an equal number of cells per patient and tissue (30 cells, 1,000 draws) retains the direction of the effect in 1,000 of 1,000 draws, with a median magnitude near 0.414. Differences in cell count between tumor and normal do not explain the shift.

**Alternative scoring and models.** The direction is consistent across four further analyses: AUCell (0.152), a z-score composite (0.207), inverse-variance-weighted aggregation (0.386), a mixed-model tissue interaction (0.361, fit at a boundary), and QC-residualized scores (0.396). The magnitude varies with method; the direction does not.

**Decomposition.** This is where the strong reading fails. The baseline composite shift (Δz about 0.387) is sensitive to a single gene. Removing FABP1 drops the retention to 0.430, below the 0.50 threshold used here, while removing HMGCS2 does not weaken it (retention about 1.0). A 5-gene core β-oxidation sub-score does not carry the shift (retention 0.350). FABP1 is a fatty-acid-binding and transport protein associated with enterocyte differentiation, not a β-oxidation enzyme such as CPT1A. The effect is therefore better described as a FABP1-linked lipid-handling signal than as a coherent FAO or β-oxidation program. On the ISC side, the score is only marginally robust: leaving out OLFM4 lowers retention to 0.502, just above threshold.

**LGR5 sensitivity.** Restricting to LGR5-detected cells to avoid circularity (LGR5 is both a stem marker and a score component) leaves only 6 patients with enough cells in both tissues. In that underpowered subset the descriptive shift is in the same direction (about 0.316), but the check is inconclusive.

**External consistency (Lee 2020).** In the Lee colorectal-cancer data, the normal-tissue correlation has the same negative sign as Pelka across both cohorts (Pelka -0.382, SMC -0.257, KUL3 -0.407; 12 of 12 informative patients negative). The tumor-versus-normal comparison cannot be made in Lee: its annotation assigns no Stem/TA-like cells to tumor tissue in either cohort. The normal anchor is consistent; the tumor-side shift is not externally replicated.

## 7. Limitations

The tissue/batch confound (Cramér's V = 1.000000) is permanent and is the dominant limitation: the biological "loss of inverse coupling" reading and a purely technical batch-to-batch reading are equally consistent with the data. The tumor-side effect rests on a single cohort (Pelka, 34 patients); Lee provides only a normal-tissue anchor. The composite signal is single-gene-sensitive (FABP1), and the ISC score is only marginally robust (OLFM4). The effect's direction is stable across scoring methods but its magnitude is not (for example scanpy 0.400629 versus AUCell 0.152), so the magnitude should be read as method-dependent. Cell-cycle scores were not available in the metadata and were not added, so residual proliferation structure within cE01 cannot be excluded. Per-patient correlations are descriptive summaries, and cells within a patient are not independent.

Going further would require a batch-decoupled design (tumor and normal of the same patient sequenced together) and a tumor-side external cohort with paired tissue.

## 8. Reproducibility

Environment: scanpy 1.11.5, numpy 2.4.2, pandas 2.3.3, scipy 1.17.1, anndata 0.12.10. The analysis lives in core_validation_v2/; each step script writes a provenance record, and every reported number traces to a committed ledger entry (CLAIMS_LEDGER.md; robustness values in ROBUSTNESS_LEDGER.md). The analysis steps run in the order 1, 2, 2b, 4, 5, 3, 6, Lee.

## References

1. Beyaz S, Mana MD, Roper J, et al. High-fat diet enhances stemness and tumorigenicity of intestinal progenitors. Nature 531:53–58 (2016).
2. Mana MD, Hussey AM, Tzouanas CN, et al. High-fat diet-activated fatty acid oxidation mediates intestinal stemness and tumorigenicity. Cell Reports 35(10):109212 (2021).
3. Shay JES, Yilmaz ÖH. Dietary and metabolic effects on intestinal stem cells in health and disease. Nature Reviews Gastroenterology & Hepatology 22(1):23–38 (2025).
4. Pelka K, Hofree M, Chen JH, et al. Spatially organized multicellular immune hubs in human colorectal cancer. Cell 184(18):4734–4752 (2021).
5. Lee HO, Hong Y, Etlioglu HE, et al. Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer. Nature Genetics 52(6):594–603 (2020). Cohorts SMC (GSE132465) and KUL3 (GSE144735).
