# Figure legends（Manuscript: MS_draft_v1）

> 图文件位于 `manuscript/figures/`，全部由 `tools/26_make_figures.py` 生成（Fig 2 自主报告提取后英化）。
> 所有数字与 `results/*.tsv` 直读一致，勾稽见 `06_amendment_numbers.md`。

---

**Fig 1. Identification of the 49 signatures.** Signatures were collected in four waves under rules fixed before
analysis (amendments 01–04). Boxes give the number of records at each step; every exclusion carries an individually
recorded reason and every record that could not be scored is listed with its extraction status. In Wave 3 one query
returned 788 hits of which the API exposed only the first 500, so 851 of 1,139 hits were screened; this is disclosed
as a limitation rather than repaired post hoc. Stratum colours are used consistently in Figs 2 and 4.

**Fig 2. Effect-size displacement from M1 to M3 for the 49 signatures.** Each row is one signature; lines connect
the unadjusted (open circle) and adjusted (filled circle) log-hazard ratio per SD of the signature score in
TCGA-LUAD (n = 482, 177 deaths). Orange-red: inflammation stratum (11); blue: metabolism stratum (32); purple:
immune stratum (6). Nine signatures with |A| > 100% (all with a non-significant M1 coefficient) and three with a
direction reversal are retained and shown as plotted; the ratio is interpretable only in the 32 signatures with a
nominally significant M1 coefficient.

**Fig 3. The random-panel null and the control axes.** **(a)** For each signature, 200 random gene panels of
matching size were drawn from the expressed transcriptome and pushed through the full M1→M3 pipeline; the observed
attenuation was placed in the resulting distribution as a percentile. The 49 percentiles compress toward the middle
(median 0.440; 9 below 0.20, 37 between 0.20 and 0.80, 3 at or above 0.80; Kolmogorov–Smirnov D = 0.267,
P = 0.001). The P value is descriptive: the percentiles derive from one patient sample and are not independent.
**(b)** The three control axes, scored and adjusted exactly like a signature. The lymphocyte/cytotoxic axis
(ρ = 0.859 with ImmuneScore) attenuates by 68.9% (P 0.053 → 0.775), showing that the pipeline detects confounding
where confounding exists; the nutrition/liver-derived axis (ρ = −0.098) attenuates by 15.3%; the neutrophil axis
has a near-zero unadjusted coefficient (β = 0.002, P = 0.976), so its attenuation ratio (−5,831%) is a numerical
artefact and is not interpreted.

**Fig 4. Stratum-level medians reproduce across cohorts; per-signature values do not.** Median attenuation by
stratum in TCGA-LUAD (OS; dark) and the two external cohorts (GSE68465, GSE31210; lighter shades of the same hue).
The ordering immune > metabolism ≈ inflammation reproduces in all three cohorts, whereas the per-signature
attenuation values correlate at approximately zero across cohorts (Spearman ρ = −0.014 and +0.006; n = 48 pairs,
SIG26 having no gene on either array platform). One signature (SIG26) is not scoreable externally, so the external
metabolism stratum contains 31 signatures.
