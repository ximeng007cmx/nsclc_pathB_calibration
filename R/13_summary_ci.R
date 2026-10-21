# =====================================================================
# 13_summary_ci.R — 汇总估计量的不确定性（bootstrap 百分位 CI）与比例 CI
#   只读 10/11/12 的结果文件，不做新模型。
# =====================================================================
suppressPackageStartupMessages({ library(data.table) })
PB <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
set.seed(20261020)

att <- fread(file.path(PB, "results", "attenuation_luad.tsv"))
lusc <- fread(file.path(PB, "results", "attenuation_lusc.tsv"))
ext <- fread(file.path(PB, "results", "external_vs_tcga.tsv"))
cox <- fread(file.path(PB, "results", "cox_models_luad.tsv"))

# 主估计量只对 |β_M1| 足够大的签名有意义（比值在 β≈0 处不稳定）——预先在此声明为描述性
att[, stable := abs(beta_M1) >= 0.10]
boot <- function(x, B = 5000) {
  x <- x[is.finite(x)]
  m <- replicate(B, median(sample(x, replace = TRUE)))
  c(median = as.numeric(median(x)), lo = as.numeric(quantile(m, .025)),
    hi = as.numeric(quantile(m, .975)), n = length(x))
}
b_all <- boot(att$att_immune)
b_st  <- boot(att[stable == TRUE, att_immune])
b_lu  <- boot(lusc$att_immune)

n1 <- att[, sum(q_M1 < .05)]; n3 <- att[, sum(q_M3 < .05)]
ci_prop <- function(k, n) {
  p <- k / n
  se <- sqrt(p * (1 - p) / n)
  c(k = k, n = n, p = p, lo = max(0, p - 1.96 * se), hi = min(1, p + 1.96 * se))
}

out <- list(
  median_attenuation_all = b_all,
  median_attenuation_stable_beta = b_st,
  median_attenuation_lusc = b_lu,
  prop_fdr_sig_M1 = ci_prop(n1, nrow(att)),
  prop_fdr_sig_M3 = ci_prop(n3, nrow(att)),
  n_stable = sum(att$stable),
  stable_ids = paste(att[stable == TRUE, sig_id], collapse = ","),
  ext_same_dir_GSE68465 = ext[endpoint %like% "68465", sum(same_direction)],
  ext_same_dir_GSE31210 = ext[endpoint %like% "31210", sum(same_direction)],
  ext_nominal_p05_GSE68465 = ext[endpoint %like% "68465", sum(p < .05)],
  ext_nominal_p05_GSE31210 = ext[endpoint %like% "31210", sum(p < .05)],
  median_coverage_GSE68465 = median(ext[endpoint %like% "68465", coverage]),
  median_coverage_GSE31210 = median(ext[endpoint %like% "31210", coverage])
)
for (nm in names(out)) {
  cat("\n##", nm, "\n"); print(out[[nm]])
}

# 落盘为一张人读表
rows <- data.table(
  item = c("median attenuation, all 13 (LUAD, M1->M3)",
           "median attenuation, |beta_M1|>=0.10 subset",
           "median attenuation, TCGA-LUSC",
           "FDR<0.05 in M1 (LUAD)", "FDR<0.05 in M3 (LUAD)",
           "external GSE68465: same direction", "external GSE68465: P<0.05",
           "external GSE31210: same direction", "external GSE31210: P<0.05",
           "median gene coverage GSE68465", "median gene coverage GSE31210"),
  value = c(
    sprintf("%.1f%% (%.1f to %.1f), n=%d", 100*b_all["median"], 100*b_all["lo"], 100*b_all["hi"], b_all["n"]),
    sprintf("%.1f%% (%.1f to %.1f), n=%d [%s]", 100*b_st["median"], 100*b_st["lo"], 100*b_st["hi"], b_st["n"], out$stable_ids),
    sprintf("%.1f%% (%.1f to %.1f), n=%d", 100*b_lu["median"], 100*b_lu["lo"], 100*b_lu["hi"], b_lu["n"]),
    sprintf("%d/13 (%.0f%%)", n1, 100*n1/13),
    sprintf("%d/13 (%.0f%%)", n3, 100*n3/13),
    sprintf("%d/13", out$ext_same_dir_GSE68465), sprintf("%d/13", out$ext_nominal_p05_GSE68465),
    sprintf("%d/13", out$ext_same_dir_GSE31210), sprintf("%d/13", out$ext_nominal_p05_GSE31210),
    sprintf("%.0f%%", 100*out$median_coverage_GSE68465), sprintf("%.0f%%", 100*out$median_coverage_GSE31210))
)
fwrite(rows, file.path(PB, "results", "headline_numbers.tsv"), sep = "\t")
print(rows)
