# =====================================================================
# 14_controls.R — 阳性对照与零分布
#   阳性对照：pathA 的淋巴细胞轴（与 ImmuneScore ρ=0.85，已知为免疫含量代理）
#            若校正方法有效，该轴应出现大幅衰减。
#   零分布：按各签名基因数抽取表达匹配随机基因面板，重复 200 次，
#            得到"任意基因面板"的衰减分布，用于判断已发表签名的衰减是否可区分。
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })
PB <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
PA <- "D:/WorkBuddy_Y02/out_ali_directions/pathA"
set.seed(20261020)

pa <- readRDS(file.path(PA, "work", "pa_data.rds")); expr <- pa$expr_log
est <- fread(file.path(PB, "work", "estimate_matched.tsv"))
reg <- fread(file.path(PB, "signature_registry.csv"))
lu  <- as.data.table(pa$surv_luad)
d <- merge(lu, est[, .(sample, Stromal, Immune, purity)], by = "sample", all.x = TRUE)
d <- d[sample %in% colnames(expr)][!is.na(time) & !is.na(status)]
d[, stage_late := factor(stage_late, levels = c("early", "late"))]
d[, gender := factor(gender)]; d[, age := as.numeric(age)]
d <- d[!is.na(age) & !is.na(gender) & !is.na(stage_late) & !is.na(Immune) & !is.na(Stromal)]
message(sprintf("完整病例 n=%d，事件 %d", nrow(d), sum(d$status)))

AXES <- list(
  lymphocyte = c("CD3D","CD3E","CD3G","CD8A","CD8B","LCK","GZMB","GZMK","PRF1","IL7R",
                 "CCL5","CXCL10","CXCR6","CD27","CXCL13","HLA-DRA"),
  neutrophil = c("ELANE","MPO","LTF","CEACAM8","S100A8","S100A9","S100A12","MMP9","PADI4",
                 "CXCR2","FPR1","DEFA4","CAMP","LCN2","ARG1","OLFM4"),
  nutrition_liver = c("ALB","TF","RBP4","APOA1","APOA2","APOB","GC","TTR","HPX","SERPINA1",
                      "AHSG","ORM1","FGA","FGG","ITIH2","PLG")
)

att_of <- function(score) {
  dd <- copy(d); dd[, sig := as.numeric(scale(score))]
  f1 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late, data = dd)
  f3 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + Immune + Stromal, data = dd)
  1 - coef(f3)["sig"] / coef(f1)["sig"]
}
score_of <- function(genes) {
  g <- intersect(genes, rownames(expr))
  as.numeric(scale(colMeans(t(scale(t(expr[g, d$sample, drop = FALSE]))), na.rm = TRUE)))
}

# ---- 阳性对照 ----
ctrl <- rbindlist(lapply(names(AXES), function(nm) {
  s <- score_of(AXES[[nm]])
  rho <- cor(s, d$Immune, method = "spearman")
  dd <- copy(d); dd[, sig := s]
  f1 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late, data = dd)
  f3 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + Immune + Stromal, data = dd)
  data.table(axis = nm, n_genes = length(AXES[[nm]]), rho_immune = rho,
             beta_M1 = coef(f1)["sig"], p_M1 = summary(f1)$coefficients["sig", "Pr(>|z|)"],
             beta_M3 = coef(f3)["sig"], p_M3 = summary(f3)$coefficients["sig", "Pr(>|z|)"],
             attenuation = 1 - coef(f3)["sig"] / coef(f1)["sig"])
}))
print(ctrl)
fwrite(ctrl, file.path(PB, "results", "control_axes.tsv"), sep = "\t")

# ---- 零分布：随机面板 ----
pool <- rownames(expr)[rowMeans(expr[, d$sample] > 0) > 0.5]  # 表达可检出基因
B <- 200
null <- rbindlist(lapply(seq_len(nrow(reg)), function(i) {
  n <- reg$n_genes[i]
  a <- replicate(B, att_of(score_of(sample(pool, n))))
  data.table(sig_id = reg$sig_id[i], n_genes = n, attenuation = a)
}))
fwrite(null, file.path(PB, "results", "null_random_panels.tsv"), sep = "\t")
nn <- null[, .(median = median(attenuation, na.rm = TRUE),
               q025 = quantile(attenuation, .025, na.rm = TRUE),
               q975 = quantile(attenuation, .975, na.rm = TRUE),
               n_draws = sum(is.finite(attenuation))), by = .(sig_id, n_genes)]
print(nn)
fwrite(nn, file.path(PB, "results", "null_summary.tsv"), sep = "\t")

# ---- 已发表签名的衰减落在随机面板分布的哪个分位 ----
att <- fread(file.path(PB, "results", "attenuation_luad.tsv"))
tmp <- merge(att[, .(sig_id, att_real = att_immune)], null, by = "sig_id")
emp <- tmp[, .(pctile = mean(attenuation <= att_real, na.rm = TRUE), n = sum(is.finite(attenuation))), by = sig_id]
emp <- merge(emp, att[, .(sig_id, att_real = att_immune, q_M1, q_M3)], by = "sig_id")
print(emp)
fwrite(emp, file.path(PB, "results", "attenuation_empirical_pctile.tsv"), sep = "\t")
message("[14] 完成")
