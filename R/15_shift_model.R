# =====================================================================
# 15_shift_model.R — 校正引起的位移是否可被"协变量结构"解释
#   线性近似：β_adj ≈ β_unadj − Σ_k γ_k · δ_k
#     γ_k = 协变量 k 对结局的（校正后）效应；δ_k = 签名分对协变量 k 的回归系数
#   若两者吻合，说明校正行为符合设计，位移方向由 ρ(签名,协变量)×协变量自身预后方向决定。
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })
PB <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
PA <- "D:/WorkBuddy_Y02/out_ali_directions/pathA"
pa <- readRDS(file.path(PA, "work", "pa_data.rds")); expr <- pa$expr_log
est <- fread(file.path(PB, "work", "estimate_matched.tsv"))
reg <- fread(file.path(PB, "signature_registry.csv"))
lu  <- as.data.table(pa$surv_luad)
d <- merge(lu, est[, .(sample, Stromal, Immune, purity)], by = "sample", all.x = TRUE)
d <- d[sample %in% colnames(expr)][!is.na(time) & !is.na(status)]
d[, stage_late := factor(stage_late, levels = c("early", "late"))]
d[, gender := factor(gender)]; d[, age := as.numeric(age)]
d <- d[!is.na(age) & !is.na(gender) & !is.na(stage_late) & !is.na(Immune) & !is.na(Stromal)]

# 协变量的"自身预后效应"（不含签名，含临床协变量；Immune/Stromal 按 SD 标准化，避免量纲不一致）
d[, `:=`(zI0 = as.numeric(scale(Immune)), zS0 = as.numeric(scale(Stromal)))]
gI <- coef(coxph(Surv(time, status) ~ zI0 + zS0 + age + gender + stage_late, data = d))
cat("gamma per SD (Immune, Stromal) =", round(gI[c("zI0", "zS0")], 4), "\n")

att <- fread(file.path(PB, "results", "attenuation_luad.tsv"))
rows <- list()
for (i in seq_len(nrow(reg))) {
  r <- reg[i]; g <- intersect(trimws(strsplit(r$genes, ",")[[1]]), rownames(expr))
  s <- as.numeric(scale(colMeans(t(scale(t(expr[g, d$sample, drop = FALSE]))), na.rm = TRUE)))
  dd <- copy(d); dd[, sig := s]
  dd[, `:=`(zI = as.numeric(scale(Immune)), zS = as.numeric(scale(Stromal)))]
  m_exp <- lm(sig ~ zI + zS, data = dd)          # δ
  delta <- coef(m_exp)[c("zI", "zS")]
  pred_shift <- -(gI[["zI0"]] * delta[["zI"]] + gI[["zS0"]] * delta[["zS"]])
  f1 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late, data = dd)
  f3 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + Immune + Stromal, data = dd)
  rows[[i]] <- data.table(sig_id = r$sig_id,
    delta_Immune = delta[["zI"]], delta_Stromal = delta[["zS"]],
    obs_shift = coef(f3)["sig"] - coef(f1)["sig"], pred_shift = pred_shift)
}
out <- rbindlist(rows)
out[, ratio := obs_shift / pred_shift]
print(out[order(-abs(obs_shift))])
ct <- cor.test(out$obs_shift, out$pred_shift, method = "spearman")
cat(sprintf("\nSpearman(observed, predicted) = %.3f, P = %.3g, n = %d\n", ct$estimate, ct$p.value, nrow(out)))
fwrite(out, file.path(PB, "results", "shift_model.tsv"), sep = "\t")
