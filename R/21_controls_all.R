# =====================================================================
# 21_controls_all.R — 扩容后（41 个签名）的阳性对照与随机面板零分布
#   与 14_controls.R 同一逻辑，仅把签名集换成 13+28 合并，并按 wave/stratum 汇总。
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })
PB <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
PA <- "D:/WorkBuddy_Y02/out_ali_directions/pathA"
set.seed(20261020)

pa <- readRDS(file.path(PA, "work", "pa_data.rds")); expr <- pa$expr_log
est <- fread(file.path(PB, "work", "estimate_matched.tsv"))
regfiles <- list.files(PB, pattern = "^signature_registry.*\\.csv$", full.names = TRUE)
regfiles <- regfiles[!grepl("excluded", regfiles)]
WAVE_OF <- c("signature_registry.csv" = "W1",
             "signature_registry_expansion.csv" = "W2",
             "signature_registry_expansion_w3.csv" = "W3",
             "signature_registry_expansion_w4.csv" = "W4")
reg <- rbindlist(lapply(regfiles, function(f) {
  d <- fread(f); w <- WAVE_OF[[basename(f)]]; if (is.null(w)) w <- "W?"
  d[, .(sig_id, n_genes, stratum, wave = w)]
}), use.names = TRUE)
setorder(reg, wave, sig_id)

lu <- as.data.table(pa$surv_luad)
d <- merge(lu, est[, .(sample, Stromal, Immune, purity)], by = "sample", all.x = TRUE)
d <- d[sample %in% colnames(expr)][!is.na(time) & !is.na(status)]
d[, stage_late := factor(stage_late, levels = c("early", "late"))]
d[, gender := factor(gender)][, age := as.numeric(age)]
d <- d[!is.na(age) & !is.na(gender) & !is.na(stage_late) & !is.na(Immune) & !is.na(Stromal)]
message(sprintf("完整病例 n=%d，事件 %d；签名 %d 个", nrow(d), sum(d$status), nrow(reg)))

score_of <- function(genes) {
  g <- intersect(genes, rownames(expr))
  as.numeric(scale(colMeans(t(scale(t(expr[g, d$sample, drop = FALSE]))), na.rm = TRUE)))
}
att_of <- function(s) {
  dd <- copy(d); dd[, sig := as.numeric(scale(s))]
  f1 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late, data = dd)
  f3 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + Immune + Stromal, data = dd)
  1 - coef(f3)["sig"] / coef(f1)["sig"]
}

AXES <- list(
  lymphocyte = c("CD3D","CD3E","CD3G","CD8A","CD8B","LCK","GZMB","GZMK","PRF1","IL7R",
                 "CCL5","CXCL10","CXCR6","CD27","CXCL13","HLA-DRA"),
  neutrophil = c("ELANE","MPO","LTF","CEACAM8","S100A8","S100A9","S100A12","MMP9","PADI4",
                 "CXCR2","FPR1","DEFA4","CAMP","LCN2","ARG1","OLFM4"),
  nutrition_liver = c("ALB","TF","RBP4","APOA1","APOA2","APOB","GC","TTR","HPX","SERPINA1",
                      "AHSG","ORM1","FGA","FGG","ITIH2","PLG")
)
ctrl <- rbindlist(lapply(names(AXES), function(nm) {
  s <- score_of(AXES[[nm]]); rho <- cor(s, d$Immune, method = "spearman")
  dd <- copy(d); dd[, sig := s]
  f1 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late, data = dd)
  f3 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + Immune + Stromal, data = dd)
  data.table(axis = nm, n_genes = length(AXES[[nm]]), rho_immune = rho,
             beta_M1 = coef(f1)["sig"], p_M1 = summary(f1)$coefficients["sig", "Pr(>|z|)"],
             beta_M3 = coef(f3)["sig"], p_M3 = summary(f3)$coefficients["sig", "Pr(>|z|)"],
             attenuation = 1 - coef(f3)["sig"] / coef(f1)["sig"])
}))
print(ctrl); fwrite(ctrl, file.path(PB, "results", "control_axes_all.tsv"), sep = "\t")

pool <- rownames(expr)[rowMeans(expr[, d$sample] > 0) > 0.5]
B <- 200
t0 <- Sys.time()
null <- rbindlist(lapply(seq_len(nrow(reg)), function(i) {
  n <- reg$n_genes[i]
  a <- replicate(B, att_of(score_of(sample(pool, n))))
  data.table(sig_id = reg$sig_id[i], n_genes = n, attenuation = a)
}))
message(sprintf("零分布完成，用时 %.1f 分钟", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
fwrite(null, file.path(PB, "results", "null_random_panels_all.tsv"), sep = "\t")

nn <- null[, .(null_median = median(attenuation, na.rm = TRUE),
               null_q025 = quantile(attenuation, .025, na.rm = TRUE),
               null_q975 = quantile(attenuation, .975, na.rm = TRUE),
               n_draws = sum(is.finite(attenuation))), by = .(sig_id, n_genes)]

att <- fread(file.path(PB, "results", "attenuation_all.tsv"))
tmp <- merge(att[, .(sig_id, att_real = att_immune, wave, stratum)], null, by = "sig_id")
emp <- tmp[, .(pctile = mean(attenuation <= att_real, na.rm = TRUE),
               n = sum(is.finite(attenuation))), by = .(sig_id, wave, stratum)]
emp <- merge(emp, att[, .(sig_id, att_real = att_immune, p_M1, q_M1, q_M3, ratio_valid)],
             by = "sig_id")
emp <- merge(emp, nn, by = "sig_id")
setorder(emp, wave, sig_id)
fwrite(emp, file.path(PB, "results", "attenuation_pctile_all.tsv"), sep = "\t")

cat("\n=== 百分位是否服从 U(0,1) ===\n")
for (s in c("ALL", "W1", "W2", "W3")) {
  x <- if (s == "ALL") emp$pctile else emp[wave == s, pctile]
  ks <- suppressWarnings(ks.test(x, "punif", 0, 1))
  cat(sprintf("%-4s n=%2d  mean=%.3f median=%.3f  KS D=%.3f P=%.3f  |  <0.20:%d  0.20-0.80:%d  >=0.80:%d\n",
              s, length(x), mean(x), median(x), as.numeric(ks$statistic), ks$p.value,
              sum(x < 0.20), sum(x >= 0.20 & x < 0.80), sum(x >= 0.80)))
}
cat("\n=== 按分层 ===\n")
for (s in c("inflammation", "metabolism", "immune")) {
  x <- emp[stratum == s, pctile]
  ks <- suppressWarnings(ks.test(x, "punif", 0, 1))
  cat(sprintf("%-13s n=%2d  mean=%.3f  KS D=%.3f P=%.3f  >=0.80:%d\n",
              s, length(x), mean(x), as.numeric(ks$statistic), ks$p.value, sum(x >= 0.80)))
}
cat("\n=== 高百分位（>=0.80）明细 ===\n")
print(emp[pctile >= 0.80, .(sig_id, wave, stratum, att_real, pctile, p_M1, q_M3, null_median)])
message("[21] 完成")
