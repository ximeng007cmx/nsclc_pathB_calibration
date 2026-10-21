# =====================================================================
# 11_diagnostics.R — 方案 B 的诊断与敏感性分析
#   ① 协变量本身的预后效应（若 Immune/Stromal 不预测 OS，则校正不可能造成衰减）
#   ② Immune/Stromal 共线性
#   ③ 绝对变化 |Δβ| 与权重敏感性（仅对符号明确的签名）
#   ④ 次要终点 TCGA-LUSC（对 LUAD 训练的签名属样本外）
#   ⑤ 比例风险假定汇总
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })

PB <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
PA <- "D:/WorkBuddy_Y02/out_ali_directions/pathA"
pa <- readRDS(file.path(PA, "work", "pa_data.rds")); expr <- pa$expr_log
est <- fread(file.path(PB, "work", "estimate_matched.tsv"))
reg <- fread(file.path(PB, "signature_registry.csv"))
lu  <- as.data.table(pa$surv_luad); ls_ <- as.data.table(pa$surv_lusc)

prep <- function(sv) {
  d <- merge(sv, est[, .(sample, Stromal, Immune, ESTIMATE, purity)], by = "sample", all.x = TRUE)
  d <- d[sample %in% colnames(expr)][!is.na(time) & !is.na(status)]
  d[, stage_late := factor(stage_late, levels = c("early", "late"))]
  d[, gender := factor(gender)]; d[, age := as.numeric(age)]
  d
}
d1 <- prep(lu); d2 <- prep(ls_)

# ---------- ① 协变量的预后效应 ----------
covs <- c("Immune", "Stromal", "purity", "ESTIMATE")
rows <- list()
for (nm in covs) {
  for (tag in c("LUAD", "LUSC")) {
    d <- if (tag == "LUAD") d1 else d2
    dd <- d[!is.na(get(nm))]
    z <- scale(dd[[nm]])[, 1]
    f <- coxph(Surv(time, status) ~ z, data = dd)
    rows[[length(rows) + 1]] <- data.table(endpoint = tag, covariate = nm, n = f$n, events = f$nevent,
      beta = coef(f)[1], HR = exp(coef(f)[1]),
      lo = exp(confint(f)[1]), hi = exp(confint(f)[2]), p = summary(f)$coefficients[1, "Pr(>|z|)"])
  }
}
cov_eff <- rbindlist(rows)
fwrite(cov_eff, file.path(PB, "results", "diag_covariate_effects.tsv"), sep = "\t")
print(cov_eff)

# ---------- ② 共线性 ----------
cm <- cor(d1[, .(Immune, Stromal, ESTIMATE, purity)], use = "pairwise.complete.obs", method = "spearman")
print(round(cm, 3))
fwrite(as.data.table(reshape2::melt(cm), keep.rownames = TRUE),
       file.path(PB, "results", "diag_covariate_correlation.tsv"), sep = "\t")

# ---------- ③ 绝对变化 + 加权敏感性 ----------
att <- fread(file.path(PB, "results", "attenuation_luad.tsv"))
att[, delta_beta := beta_M3 - beta_M1]
att[, delta_pct_of_M1 := delta_beta / abs(beta_M1)]
print(att[, .(sig_id, n_genes, beta_M1, beta_M3, att_immune, delta_beta, q_M1, q_M3)])

# 权重（仅当原文符号明确）：SIG05、SIG11
W <- list(
  SIG05 = c(CCL20 = 0.1241, CCR2 = -0.6192, GNAI3 = 0.3660, ITGA5 = 0.1240, NMI = 0.2050,
            PCDH7 = 0.1330, PSEN1 = 0.430, SLC11A2 = -0.3382),
  SIG11 = c(LDHA = 0.269752653768969, NSDHL = 0.114984371234971, TP53INP2 = 0.0310569577143764,
            FLT1 = 0.17593670859331, IRS1 = 0.103773854013971, ELOVL7 = 0.133274629300893,
            AGPS = 0.0763835280754936, FHL2 = 0.00688026723004499, MED6 = 0.549703365193899,
            PLIN3 = 0.109916205934083, VDAC1 = 0.0437103013689035, SULT2B1 = 0.160609136698149)
)
wrows <- list()
for (sid in names(W)) {
  w <- W[[sid]]; g <- intersect(names(w), rownames(expr))
  sub <- expr[g, d1$sample, drop = FALSE]
  s <- as.numeric(crossprod(w[g], sub))
  s <- as.numeric(scale(s))
  dd <- copy(d1); dd[, sig := s]
  f1 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late, data = dd)
  f3 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + Immune + Stromal, data = dd)
  wrows[[sid]] <- data.table(sig_id = sid, scoring = "published_weights",
    beta_M1 = coef(f1)["sig"], beta_M3 = coef(f3)["sig"],
    p_M1 = summary(f1)$coefficients["sig", "Pr(>|z|)"],
    p_M3 = summary(f3)$coefficients["sig", "Pr(>|z|)"],
    att = 1 - coef(f3)["sig"] / coef(f1)["sig"])
}
wt <- rbindlist(wrows)
fwrite(wt, file.path(PB, "results", "diag_weighted_scores.tsv"), sep = "\t")
print(wt)

# 等权（对照）
ew <- att[sig_id %in% names(W), .(sig_id, scoring = "unweighted", beta_M1, beta_M3, att = att_immune)]
print(rbind(ew, wt, fill = TRUE))

# ---------- ④ 次要终点 LUSC ----------
srows <- list()
for (i in seq_len(nrow(reg))) {
  r <- reg[i]; genes <- trimws(strsplit(r$genes, ",")[[1]])
  g <- intersect(genes, rownames(expr))
  s <- as.numeric(scale(colMeans(t(scale(t(expr[g, d2$sample, drop = FALSE]))), na.rm = TRUE)))
  dd <- copy(d2); dd[, sig := s]
  cc <- dd[!is.na(sig) & !is.na(age) & !is.na(gender) & !is.na(stage_late)]
  f1 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late, data = cc)
  cc3 <- cc[!is.na(Immune) & !is.na(Stromal)]
  f3 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + Immune + Stromal, data = cc3)
  srows[[i]] <- data.table(endpoint = "TCGA-LUSC, OS", sig_id = r$sig_id, n = f1$n, events = f1$nevent,
    beta_M1 = coef(f1)["sig"], beta_M3 = coef(f3)["sig"],
    p_M1 = summary(f1)$coefficients["sig", "Pr(>|z|)"],
    p_M3 = summary(f3)$coefficients["sig", "Pr(>|z|)"],
    att_immune = 1 - coef(f3)["sig"] / coef(f1)["sig"])
}
lusc <- rbindlist(srows)
lusc[, `:=`(q_M1 = p.adjust(p_M1, "BH"), q_M3 = p.adjust(p_M3, "BH"))]
fwrite(lusc, file.path(PB, "results", "attenuation_lusc.tsv"), sep = "\t")
print(lusc[, .(sig_id, beta_M1, beta_M3, att_immune, q_M1, q_M3)])

sm2 <- data.table(endpoint = "TCGA-LUSC, OS", n_signatures = nrow(lusc),
  median_attenuation = median(lusc$att_immune, na.rm = TRUE),
  n_FDR_sig_M1 = sum(lusc$q_M1 < 0.05, na.rm = TRUE),
  n_FDR_sig_M3 = sum(lusc$q_M3 < 0.05, na.rm = TRUE))
print(sm2)
fwrite(sm2, file.path(PB, "results", "summary_lusc.tsv"), sep = "\t")

# ---------- ⑤ PH 假定 ----------
print(att[, .(sig_id, ph_p_M3)])
message("[11] 完成")
