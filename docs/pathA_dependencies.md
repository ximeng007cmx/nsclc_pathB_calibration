# Companion (pathA) dependency checksums

The pathB pipeline reads three files from the companion `pathA` working directory
and uses the raw ArrayExpress downloads under `pathA/data/external/`.
These files accompany the publication (available from the corresponding author)
or are obtainable from the original public sources. Verify any obtained copies
against the SHA-256 checksums below before running the `R/` pipeline.

## Key derived files

```
# pa_data.rds (TCGA-LUAD/LUSC expression + survival; companion pathA pipeline)
a85a3e15af4af51925c3b3a150b337f7414f41e46c2c17dc3d711a00ba28ad40 *../pathA/work/pa_data.rds
# estimate_scores_matched.tsv (pathA ESTIMATE output; pathB uses work/estimate_matched.tsv dumped from it)
220933360251b8c67260828ef860a048dda898edfad09803cc7ef269eeabb380 *../pathA/results/estimate_scores_matched.tsv
# immune_tool_crosscheck_summary.tsv (ten-tool cross-check input)
3873fc28cd8c698acaf0922b7997d2645e3c236ae5f631bda64708e63bf91cff *../pathA/results/immune_tool_crosscheck_summary.tsv
```

## Raw external arrays (E-GEOD-68465, E-GEOD-31210)

716 files, ~594 MB total. Full checksum list: `pathA_external_data.sha256`
(same directory as this document; paths are relative to `pathA/`).

Sources:
- E-GEOD-68465 (Director's Challenge LUAD, HG-U133A): EBI BioStudies / ArrayExpress
- E-GEOD-31210 (Japanese LUAD IA-II, HG-U133 Plus 2): EBI BioStudies / ArrayExpress
- ESTIMATE scores in `pathA/results/estimate_scores_matched.tsv` derive from
  TCGA-LUAD/LUSC via the ESTIMATE method (Yoshihara et al., Nat Commun 2013).
