# Immune and stromal content does not explain the effect sizes of published inflammation- and metabolism-related prognostic signatures in non-small cell lung cancer: a calibration study of 49 signatures with positive and negative controls

**Draft v1 · 2026-10-21 · target journal: Cancer Medicine**

> Citation note: the text uses [n] numbering; the References list 9 methodological and data-source references; the 49 signature references are listed individually in Table 1 (with PMID/DOI).
> **All 58 references passed dual-source verification (Europe PMC + CrossRef) on 2026-10-21** (49 + 9, 0 issues);
> verification records in `references_verified.md` and `logs/refs_verify_49.log`; volume/issue/pages taken from EPMC core records, not from memory.
> One mis-quoted DOI corrected: the DOI for Liu 2018 is 10.1016/j.cell.2018.02.052 (the earlier record 10.1016/j.cell.2018.03.046 resolved to a different paper in the same journal).
> All values come from `results/*.tsv` (see the quick-reference table in §9 of `正文起草计划_B路线.md`); no value was transcribed by hand.

---

## Abstract

**Background:** Prognostic transcriptomic signatures for non-small cell lung cancer (NSCLC) are published continuously, and tumour microenvironment (TME) content is known to shape bulk expression. Whether immune and stromal content explains the effect sizes of published signatures has not been tested systematically with controls.

**Methods:** We assembled 49 published inflammation-related (n = 11), metabolism-related (n = 32) and immune-related (n = 6) prognostic signatures under rules fixed before analysis and locked by hash. In TCGA-LUAD (482 patients, 177 deaths) we fitted Cox models with the signature score, age, sex and stage (M1) and additionally ImmuneScore and StromalScore (M3), and quantified attenuation as A = 1 − β<sub>M3</sub>/β<sub>M1</sub>. Uncertainty was obtained by cluster bootstrap (5,000 resamples). We used three control axes and per-signature random gene panels (200 draws) as positive and negative controls, and replicated directionality in two independent cohorts (GSE68465, n = 439; GSE31210, n = 226) and across endpoints (progression-free interval, disease-specific survival).

**Results:** Median attenuation across the 49 signatures was −1.4% (95% CI −7.4% to +1.8%); 21 signatures shrank, 25 grew and 3 reversed direction, and among the 32 nominally significant signatures the split was 17:15. The lymphocyte control axis attenuated by 68.9% (P changed from 0.053 to 0.775), confirming power. Attenuation percentiles were compressed toward zero relative to random panels (median 0.440; Kolmogorov–Smirnov P = 0.001). Displacement magnitude tracked signature–ImmuneScore correlation (partial ρ = 0.373, P = 0.009) whereas its sign did not (ρ = −0.125, P = 0.391). The three-stratum ordering (immune > metabolism ≈ inflammation) reproduced in both external cohorts and across endpoints.

**Conclusions:** Immune and stromal content is not an explanatory variable for the effect sizes of published NSCLC signatures; published signatures attenuate even less than random gene panels. Signature studies should report TME-adjusted effect sizes rather than treat TME content as a confounder to be removed.

**Keywords:** non-small cell lung cancer; prognostic signature; tumour microenvironment; ESTIMATE; calibration; reproducibility

---

## 1. Introduction

Bulk transcriptomic prognostic signatures for NSCLC continue to appear at a high rate, and a substantial fraction of them are anchored to inflammation or to metabolic reprogramming. In parallel, the composition of the tumour microenvironment (TME; the fraction of immune and stromal cells within a specimen) is itself associated with outcome and is a major determinant of bulk expression profiles [1]. Because signature scores and TME estimates are often strongly correlated, any association between a published signature and survival can be decomposed into a component carried by the signature itself and a component shared with TME content. Which component dominates has direct consequences: if TME content explains most of the association, a signature is largely a surrogate for an already measurable quantity; if it explains little, the signature carries information that TME estimates do not.

Most signature publications do not report this decomposition. The omission is not incidental. Removing TME content from a signature's effect is a strong test, and the result is unfavourable in a way that is difficult to publish: reviewers reasonably ask whether the correction has any power at all, and without a positive control the answer cannot be given. Studies that do adjust for TME content tend to pool immune-anchored and metabolism-anchored signatures into a single estimate, which is itself problematic: an immune-anchored signature is expected to attenuate strongly after adjustment because it is, by construction, a proxy for immune content. Mixing such signatures into a pooled estimate inflates attenuation for the whole set, a form of predictive confounding rather than a circularity, but one that produces the same practical error.

Here we report a calibration study rather than a new signature. We assembled 49 published inflammation-, metabolism- and immune-related prognostic signatures under rules fixed before analysis, quantified how much each signature's effect size changed after adjustment for tumour purity and immune/stromal content, and surrounded that estimate with positive and negative controls, a random-panel null distribution, two independent external cohorts and two additional survival endpoints. Our question is narrow and deliberately so: not whether published signatures are valid, but whether immune and stromal content explains their effect sizes.

---

## 2. Methods

### 2.1 Signature identification and the three-stratum framework

Signatures were collected in four waves, and every step is reported at the level of individual records so that the funnel closes exactly.

Wave 1 (W1) comprised 13 signatures identified in a preregistered search.

Wave 2 (W2) expanded the metabolism stratum with six Europe PMC queries restricted to 2020 onward, anchored on amino acid, cholesterol, glycolysis, glutamine and energy metabolism: 767 records, 538 after de-duplication, 431 after applying the eligibility filter (PMCID present, absent from the existing list), 99 retained as tier A on the basis of a readable gene list in the abstract, and 62 sent to full-text extraction. Of these 62, 28 were included, 26 were excluded with an individually recorded reason, and 8 were never scored: six because the full text yielded no verifiable gene list and two because the record could not be retrieved (HTTP 500).

Wave 3 (W3) targeted the inflammation stratum with eight inflammation-specific queries: 1,139 hits, of which 851 were retrieved because one query returned 788 hits and the API returned only the first 500; 605 new candidates; 74 retained as tier A; 20 excluded because the title anchored an immune construct and 6 because the construct lay outside inflammation; 48 sent to full-text extraction. These 48 partition exactly: 2 included; 6 excluded with individually recorded reasons, three of which had a verifiable gene list and three of which failed the stated-count check as well as the construct check; 4 that had a verifiable list but no documented inflammation or chemokine provenance; and 36 that yielded no scorable record (25 without a verifiable gene list, 4 without a stated gene count, 5 with a count mismatch, 2 unretrievable). Only 16 of the 48 carried documentary evidence that the candidate pool derived from an inflammation or chemokine gene set, and only 9 produced a verifiable list with a matching stated count, which is why the stratum is small.

Wave 4 (W4) re-examined the 20 immune-anchored candidates that the construct rule had excluded from the main analysis, 16 of which had a verifiable gene list; 6 were included as a separate immune stratum and 10 excluded with recorded reasons.

Two rules not present in the original eligibility criteria were introduced in W2 and are disclosed here. First, signatures whose title anchored an explicitly immune-related construct (immune-related genes, immune-related signature, immune checkpoint) were **not** admitted to the main analysis; in W4 these were adjudicated separately as the immune stratum rather than discarded, because adjusting an immune signature for immune content is expected to attenuate it by construction and pooling such signatures with metabolism signatures would bias the pooled estimate. Second, signatures whose construct fell outside inflammation and nutrition/metabolism (ferroptosis, cuproptosis, disulfidptosis, pyroptosis, autophagy, m6A, hypoxia, ECM/CAF, lysosome, drug tolerance, DNA repair) were excluded regardless of statistical performance.

The final set comprised 49 signatures: 11 inflammation-related, 32 metabolism-related and 6 immune-related (W1 13, W2 28, W3 2, W4 6; Fig 1). The registry, the four amendment documents and the SHA-256 hashes of each wave are provided as supplementary material. Every candidate that entered the analysis has a registry row, every exclusion has a recorded reason, and every record that could not be scored is listed with its extraction status. No signature was removed on the basis of its result.

### 2.2 Cohorts

The primary cohort was TCGA lung adenocarcinoma (LUAD). Expression (log2 FPKM, Entrez-mapped to gene symbols) was combined with TCGA-CDR clinical data [2] and ESTIMATE scores [1]. The analysis set comprised 482 patients with 177 deaths after requiring non-missing age, sex, stage and (for M3) ImmuneScore and StromalScore.

Two independent cohorts were used for directional replication: GSE68465 (Director's Challenge, LUAD, HG-U133A; 12,545 genes; 439 patients with complete covariate data, 235 deaths) and GSE31210 (Japanese LUAD, stages IA–II, HG-U133 Plus 2; 16,504 genes; 226 patients with complete data, 35 deaths; age unavailable), both obtained through ArrayExpress mirrors of the GEO series [7–9]. A histology check used TCGA lung squamous cell carcinoma (LUSC; 482 patients, 209 deaths).

### 2.3 Signature scoring

For each signature, genes present on the expression platform were selected and z-scored across samples; the score was the mean of the z-scored values across genes, standardised across samples. Signatures with fewer than two available genes were not scored. Gene coverage was recorded for every signature–cohort pair.

### 2.4 Adjustment models

Four nested Cox proportional hazards models were fitted for overall survival: M0, signature score alone; M1, score + age + sex + stage (early vs late); M2, M1 + tumour purity; M3, M1 + ImmuneScore + StromalScore. Because purity and the ESTIMATE scores were collinear (Spearman ρ = −1.000) and the M2 and M3 signature coefficients were nearly identical (r = 0.998), M2 carried no independent information; the primary contrast was therefore M1 versus M3, and M2 results are reported in supplementary material only. Proportional hazards were checked with scaled Schoenfeld residuals.

### 2.5 Attenuation estimator and its validity boundary

Attenuation was defined on the log-hazard scale as A = 1 − β<sub>M3</sub>/β<sub>M1</sub>, so that A > 0 indicates shrinkage after adjustment, A < 0 indicates enlargement, and A > 1 or A < 0 with |A| > 100% indicates that β<sub>M1</sub> was close to zero and the ratio is numerically unstable. Nine signatures met the |A| > 100% criterion; all had a non-significant M1 coefficient. For this reason every summary is reported twice: over all signatures, and over the subset with a nominally significant M1 coefficient (P < 0.05), the only subset in which the ratio is interpretable. The latter subset is labelled exploratory throughout.

### 2.6 Uncertainty and multiplicity

Confidence intervals for medians were obtained by cluster bootstrap over signatures (5,000 resamples, seed 20261020) rather than over patients, because the unit of inference is the signature. Differences between strata were bootstrapped with the same scheme. Benjamini–Hochberg correction was applied across all 49 signatures for each model, so that q values reflect the full set.

### 2.7 Control axes and the random-panel null

Three control axes were scored with the same pipeline: a lymphocyte/cytotoxic axis (16 genes, expected to be a strong proxy for ImmuneScore), a neutrophil axis (16 genes) and a nutrition/liver-derived axis (16 genes, expected to be unrelated to immune content). For each of the 49 signatures, 200 random gene panels of identical size were drawn from the expressed transcriptome, the full M1→M3 pipeline was re-run on each panel, and the observed attenuation was placed in the resulting distribution as a percentile. Deviation from uniformity was assessed by Kolmogorov–Smirnov distance against U(0,1); because the 49 percentiles derive from the same patient sample and are therefore not independent, the associated P value is reported as descriptive only.

### 2.8 External directional replication

The full pipeline (same genes, same scoring, same M1 and M3, same estimator) was re-run in GSE68465 and GSE31210. Because the ESTIMATE package was not obtainable at the time of analysis, stromal and immune scores in these cohorts were recomputed from the 141-gene stromal and 141-gene immune signatures published with the original method [1], using single-sample GSEA as implemented in the GSVA package [3,4]. A linear rescaling of a covariate does not change the coefficient of the signature in a Cox model, so the attenuation ratio is invariant to this substitution; only the absolute scale of the scores differs from the official implementation. Age was available in GSE68465 but not in GSE31210 (all 226 records missing), so the GSE31210 M1 model omitted age; one signature (SIG26) had no gene present on either array platform and therefore no external value. This analysis was post hoc and is labelled as such.

### 2.9 Endpoint sensitivity and histology check

Overall survival was the primary endpoint. Progression-free interval (PFI) and disease-specific survival (DSS) were analysed with identical models in 482 and 447 patients respectively. The full pipeline was also re-run in TCGA-LUSC to ask whether the directionality was histology-specific; this analysis was post hoc and its interpretation is restricted as described in Results.

### 2.10 Software

Analyses used R 4.6.1 (survival, data.table, GSVA) and Python 3.13 for reporting. Every number in this manuscript is generated by scripts that read the result tables directly; no value was transcribed by hand. Scripts, registries and hash records are provided as supplementary material.

---

## 3. Results

### 3.1 The signature set, and why the inflammation stratum is small

The 49 included signatures contributed 11 inflammation, 32 metabolism and 6 immune signatures. The inflammation stratum proved resistant to expansion: eight inflammation-specific queries produced 1,139 hits, 851 retrieved records, 605 new candidates and 74 with a readable gene list, but only 16 had documentary evidence that the candidate pool derived from an inflammation or chemokine gene set, and 2 survived manual adjudication. The structural reason is that, in the NSCLC signature literature, inflammation is almost always published as *immune-related*: of the 74 tier-A candidates, 20 were excluded because their titles anchored an immune construct, and 16 of those 20 nevertheless contained a usable gene list. Inflammation as a construct has been absorbed by the immune construct in this literature, which is precisely why immune confounding is easily mistaken for an inflammation effect.

### 3.2 Main result: no systematic direction

Across the 49 signatures, median attenuation was −1.4% (95% CI −7.4% to +1.8%). Direction was close to balanced: 21 signatures shrank, 25 grew and 3 reversed sign. Restricting to the 32 signatures with a nominally significant unadjusted coefficient, the split was 17 shrinking versus 15 growing (exact two-sided binomial P = 0.86). Adjustment therefore had no systematic direction. The same picture holds at the level of false discovery: 31 signatures were significant at q < 0.05 before adjustment and 31 after, with two losses and two gains. The widely expected pattern of a class-wide loss of significance did not occur, which argues against the "unmasking" narrative in which published signatures are predominantly TME artefacts (Fig 2).

### 3.3 The control axes confirm that the pipeline has power

The lymphocyte axis correlated with ImmuneScore at ρ = 0.859 and attenuated by 68.9% after adjustment (P changed from 0.053 to 0.775): the same correction that leaves the 49 published signatures essentially unchanged abolishes a known immune proxy. The nutrition/liver-derived axis, which is unrelated to immune content (ρ = −0.098), attenuated by only 15.3%. The neutrophil axis had a near-zero unadjusted coefficient (β<sub>M1</sub> = 0.002, P = 0.976) and its attenuation (−5,831%) is a ratio artefact; it is listed but not interpreted (Fig 3b). The contrast between axes indicates that the pipeline detects confounding where confounding exists. A post hoc cross-check against ten alternative immune-deconvolution tools in the same cohort was concordant: the lymphocyte axis correlated significantly and positively with all ten (median Spearman ρ = 0.69) and the neutrophil axis with all ten (median ρ = 0.28), whereas the nutrition/liver-derived axis showed the weakest and least consistent association (median ρ = 0.16, 7 of 10 tools positive and significant); the control pattern is therefore not an artefact of the ESTIMATE implementation [5,6].

### 3.4 Published signatures attenuate less than random gene panels

For each signature, 200 random panels of matching size were passed through the same pipeline. The observed attenuation percentiles were not uniformly distributed but compressed toward the middle (median 0.440; 9 below 0.20, 37 between 0.20 and 0.80, 3 at or above 0.80; Kolmogorov–Smirnov D = 0.267, P = 0.001). In other words, the published signatures moved *less* under adjustment than randomly drawn gene sets did (Fig 3a). The three signatures in the upper tail (SIG01, SIG03, SIG26) all had non-significant unadjusted coefficients (P > 0.6), so their large attenuation values are ratio artefacts rather than evidence of strong confounding. Two caveats apply and are stated in the figure: the percentiles derive from a single patient sample and are not independent, so the P value is descriptive; and the dispersion of the random-panel distribution is itself large: a small random panel can appear to shrink or grow dramatically.

### 3.5 Displacement magnitude tracks correlation; its direction does not

The preregistered hypothesis was that attenuation would increase with the correlation between a signature and ImmuneScore. The data support this for magnitude but not for direction: |Δβ| correlated with |ρ| at 0.445, remained associated after controlling for |β<sub>M1</sub>| (partial ρ = 0.373, P = 0.009), while the signed attenuation was unrelated to ρ (ρ = −0.125, P = 0.391). The correlation between |Δβ| and |β<sub>M1</sub>| was negative (−0.355), excluding the trivial explanation that large effects move more. Adjustment therefore pushes effect sizes by an amount proportional to how strongly a signature tracks immune content, but pushes them in inconsistent directions. This is why the correct statement is not "adjustment weakens signatures" but "adjustment moves signatures according to their immune correlation, in either direction".

### 3.6 Three strata, and two preregistered predictions

Two directional predictions were registered before analysis. P1 anticipated greater attenuation in the inflammation stratum than in the metabolism stratum, on the grounds that inflammation signatures correlate more strongly with ImmuneScore (median |ρ| 0.482 vs 0.146). P1 was not supported and in fact ran in the opposite direction: the difference (inflammation minus metabolism) was −10.3% (95% CI −109.6% to +22.5%). P2 anticipated greater attenuation in the immune stratum than in either other stratum. The point estimates were in the predicted direction (immune minus inflammation +18.3%, 95% CI −75.0% to +118.8%; immune minus metabolism +8.1%, 95% CI −83.8% to +21.8%), but both intervals included zero, so P2 is not confirmed. The immune stratum itself (6 signatures) had median attenuation +8.4% (95% CI −84.8% to +20.5%), with 4 of 6 signatures attenuating, consistent in direction with the positive control but with an interval too wide to support inference at n = 6.

### 3.7 Directionality replicates externally, but only at the level of strata

In GSE68465 the immune stratum attenuated by 14.2% (4 of 6 signatures) and in GSE31210 by 5.3% (4 of 6); the inflammation stratum was negative in both (−5.5% and −2.4%) and the metabolism stratum was near zero (+1.3% and −1.5%). The ordering immune > metabolism ≈ inflammation therefore reproduced in all three cohorts (Fig 4). Individual values, however, did not: the Spearman correlation between a signature's attenuation in TCGA and in GSE68465 was −0.014, and in GSE31210 0.006. The reproducible unit is the stratum-level pattern, not the per-signature number, which is consistent with most per-signature ratios being noise-dominated.

### 3.8 Endpoint sensitivity

Re-running the pipeline with PFI and DSS gave median attenuation of +0.3% (95% CI −4.0% to +2.4%) and +2.0% (95% CI −2.1% to +5.6%), against −1.3% for overall survival in the identical endpoint pipeline; counting reversals separately, the direction split was 24 shrinking, 22 growing and 3 reversing for PFI and 25, 21 and 3 for DSS, and the three-stratum ordering was identical across endpoints (immune +8.5%, +6.6%, +19.7%; metabolism near zero; inflammation −9.9%, −11.1%, −0.8%). Attenuation values correlated across endpoints (ρ = 0.418 for PFI vs OS, 0.639 for DSS vs OS) but not across cohorts, indicating that the variation in A is driven mainly by the cohort rather than by the endpoint.

### 3.9 Histology check: the denominator fails in LUSC

In TCGA-LUSC (482 patients, 209 deaths) none of the 49 signatures reached significance at q < 0.05 and only five had a nominally significant unadjusted coefficient, with a Spearman correlation of −0.17 between the unadjusted coefficients in LUAD and LUSC. Because the denominator of A approaches zero in this cohort, the attenuation ratio is not defined and the cohort cannot be used for a directional test. Its direction counts (24 shrinking, 18 growing, 7 reversed) are reported as a side observation only. This is a limitation of the signature set's transportability, not evidence against the primary finding: the primary analysis concerns how much an existing effect changes, and requires that an effect exist.

---

## 4. Discussion

### 4.1 Principal finding

In 49 published NSCLC signatures spanning inflammation, metabolism and immune constructs, adjustment for tumour purity and immune/stromal content changed effect sizes by a median of −1.4%, with a direction split close to 1:1 and no class-wide loss of significance. The same adjustment abolished a known immune proxy, so the null result is not a power failure. Immune and stromal content, in other words, is not the variable that explains why these signatures carry prognostic information.

### 4.2 Published signatures attenuate less than random panels

The most informative result of this study is not the median but the distribution. Random gene panels of matched size, pushed through the same pipeline, produced a wide spread of apparent attenuation, with frequent large apparent shrinkage or growth. The published signatures occupied the middle of that distribution more densely than chance would predict. Two readings are compatible with this. The favourable reading is that published signatures are not arbitrary gene lists: they behave more stably under adjustment than random sets, as expected if their association with outcome is partly independent of TME content. The methodological reading is that a single signature's attenuation carries little information on its own, because the sampling variability of the ratio is large; only the aggregate, and specifically the direction split, is interpretable. Both readings argue against the practice of judging a signature by whether it "survives" TME adjustment in a single cohort.

### 4.3 The three-stratum picture, and why the inflammation stratum is thin

The three-stratum ordering (immune signatures attenuating most, metabolism signatures not at all, inflammation signatures in between and slightly negative) reproduced in both external cohorts and across all three endpoints, although no single contrast reached significance and the immune stratum has only six members. The ordering is stable across cohorts even though per-signature values correlate at approximately zero; stratum-level behaviour is therefore a property of the construct, whereas per-signature attenuation is a property of the cohort.

The inflammation stratum could not be expanded beyond 11 signatures, and the reason is informative rather than incidental. In this literature, inflammation is nearly always published as an immune construct: four-fifths of the inflammation-specific candidates that the construct rule excluded (16 of 20) carried a usable gene list but an immune title. Any study that anchors on "inflammation signatures" is therefore, in practice, sampling immune signatures, which is exactly the mechanism by which immune confounding is misattributed to inflammation.

### 4.4 Relation to previous work

Our results are consistent with, but distinct from, the literature on signature reproducibility. That literature documents that many published signatures do not validate in independent cohorts; we did not test validation and our data do not speak to it directly. What we tested is whether one specific, measurable quantity, TME content, accounts for the effect sizes that published signatures do report, and the answer is that it accounts for very little. The two findings are compatible: a signature can be non-reproducible across cohorts and still not be explained by TME content. Our histology check illustrates the distinction concretely: in LUSC the signatures largely failed to show an effect at all, which is a transportability problem, not an attenuation problem.

### 4.5 Implications and limitations

The practical implication is a reporting one. Because attenuation magnitude tracks the signature–ImmuneScore correlation while its direction does not, reporting a single adjusted hazard ratio without the unadjusted one, or reporting only the adjusted one, discards exactly the information that distinguishes a TME-driven signature from an independent one. We suggest that signature studies report the unadjusted and TME-adjusted coefficients together, and state the direction of change.

Several limitations bound these conclusions. The primary cohort is a single TCGA cohort; the external cohorts replicate directionality but not per-signature values, and both are array-based with partial gene coverage. The post hoc components (signature expansion, external replication, endpoint sensitivity and the histology check) are labelled as such, and the expansion rules, although fixed before analysis, were written after the first wave of results was known. The signature pool is also bounded by what the searches could retrieve: one inflammation-specific query returned 788 hits of which the API exposed only the first 500, and 48 records that reached full-text extraction could not be converted into a scoreable signature (8 in W2; in W3, 36 that yielded no scorable record plus 4 with no documented pool provenance) for the reasons itemised in Methods. The inflammation stratum may therefore be incomplete: a limitation on the completeness of the construct, not on the estimates for the signatures that were included. The immune stratum contains six signatures, so its interval is wide by construction. The attenuation ratio is undefined when the unadjusted coefficient is near zero, which affects nine signatures and one control axis; summaries are therefore given both over all signatures and over the nominally significant subset. Purity and the ESTIMATE scores are collinear, so the purity model adds nothing and is relegated to supplementary material. Finally, the percentile analysis derives from a single patient sample and its P values are descriptive. These limitations qualify the precision of the estimates, not the direction of the conclusion.

---

## 5. Conclusion

Across 49 published NSCLC prognostic signatures, adjustment for immune and stromal content produced no systematic change in effect size, while the same adjustment abolished a known immune proxy. Immune and stromal content is therefore not an explanatory variable for these signatures, and published signatures attenuate even less than random gene panels. Signature studies should report TME-adjusted effect sizes alongside unadjusted ones and state the direction of change, rather than treating TME content as a confounder presumed to inflate their results.

## Data and code availability

All analysis code (R and Python), the four hash-locked signature registries, all result tables and the manuscript source are archived at https://github.com/ximeng007cmx/nsclc_pathB_calibration (v1.0.0) and mirrored on Zenodo: **10.5281/zenodo.23537008**. The Zenodo deposit additionally contains the calibrated TCGA-LUAD/LUSC expression matrix and survival tables (pa_data.rds) required to re-run the pipeline. SHA-256 checksums of all archived files are recorded in the repository (MANIFEST.sha256; docs/pathA_dependencies.md). Expression and clinical data originate from TCGA/TCGA-CDR, the external cohorts from GEO/ArrayExpress (E-GEOD-68465, E-GEOD-31210) and the full-text extraction corpus from PubMed Central; these third-party public data are subject to their original terms and are not redistributed beyond the derived files listed above.

---

## References

> All references (9 in-text methodological/data-source references + 49 signature references) passed dual-source verification (Europe PMC + CrossRef) on 2026-10-21, 0 issues; itemised records in `references_verified.md`; the signature references are also listed in Table 1.

1. Yoshihara K, Shahmoradgoli M, Martínez E, Vegesna R, Kim H, Torres-Garcia W, et al. Inferring tumour purity and stromal and immune cell admixture from expression data. Nat Commun;4:2612. PMID 24113773. doi:10.1038/ncomms3612. **[verified]**
2. Liu J, Lichtenberg T, Hoadley KA, Poisson LM, Lazar AJ, Cherniack AD, et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell;173(2):400-416.e11. PMID 29625055. doi:10.1016/j.cell.2018.02.052. **[verified; DOI corrected, the earlier record 10.1016/j.cell.2018.03.046 was a mis-quote]**
3. Barbie DA, Tamayo P, Boehm JS, Kim SY, Moody SE, Dunn IF, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature;462(7269):108-112. PMID 19847166. doi:10.1038/nature08460. **[verified]**
4. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics;14:7. PMID 23323831. doi:10.1186/1471-2105-14-7. **[verified]**
5. Becht E, Giraldo NA, Lacroix L, Buttard B, Elarouci N, Petitprez F, et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biol;17(1):218. PMID 27765066. doi:10.1186/s13059-016-1070-5. **[verified]**
6. Aran D, Hu Z, Butte AJ. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biol;18(1):220. PMID 29141660. doi:10.1186/s13059-017-1349-1. **[verified]**
7. Director's Challenge Consortium for the Molecular Classification of Lung Adenocarcinoma, Shedden K, Taylor JM, Enkemann SA, Tsao MS, Yeatman TJ, et al. Gene expression-based survival prediction in lung adenocarcinoma: a multi-site, blinded validation study. Nat Med;14(8):822-827. PMID 18641660. doi:10.1038/nm.1790. **[verified]**
8. Okayama H, Kohno T, Ishii Y, Shimada Y, Shiraishi K, Iwakawa R, et al. Identification of genes upregulated in ALK-positive and EGFR/KRAS/ALK-negative lung adenocarcinomas. Cancer Res;72(1):100-111. PMID 22080568. doi:10.1158/0008-5472.can-11-1403. **[verified]**
9. Athar A, Füllgrabe A, George N, Iqbal H, Huerta L, Ali A, et al. ArrayExpress update - from bulk to single-cell expression data. Nucleic Acids Res;47(D1):D711-D715. PMID 30357387. doi:10.1093/nar/gky964. **[verified]**

The 49 source references for SIG01–SIG49: PMID/DOI for each are given in Table 1 and `references_verified.md` (generated automatically from the four wave registries; all verified).

---

## Figures and tables

- **Fig 1** (`figures/Fig1_flow.svg`) Four-wave signature inclusion flow: the funnel closes level by level, including the W3 retrieval truncation (788→500) and the disposition of the unscored records
- **Fig 2** (`figures/Fig2_forest.svg`) M1→M3 forest plot of the 49 signatures (three stratum colours; extracted from the project report, legend anglicised, xmlns added)
- **Fig 3** (`figures/Fig3_null_controls.svg`) a: random-panel percentile distribution (median 0.440, KS D = 0.267, P = 0.001); b: three control axes (lymphocyte 68.9%, nutrition 15.3%, neutrophil −5,831% ratio artefact)
- **Fig 4** (`figures/Fig4_cohorts_strata.svg`) Median A by three cohorts × three strata (immune > metabolism ≈ inflammation reproduced in all three)
- Captions: `manuscript/figure_captions.md` (English)

- **Table 1** (`tables.xlsx` / `tables.md`) List of the 49 signatures (SIG/PMID/DOI/construct/wave/gene count/pool provenance): generated from the registries; all PMID/DOI verified
- **Table 2** (`tables.xlsx` / `tables.md`) Multi-calibre summary (ALL/W1–W4/three strata, with bootstrap CIs): generated
- **Table 3** (`tables.xlsx` / `tables.md`) External re-estimation and endpoint sensitivity summary: generated
- **Supplementary** per-signature full-field table, exclusion list, LUSC check table, low-coverage list
