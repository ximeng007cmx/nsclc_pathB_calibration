# =====================================================================
# 10_calibration.R — 方案 B 主分析：已发表签名的纯度/免疫含量校准式重评估
#
# 预注册：../00_criteria_preregistered.md（sha256 见 work/LOCK.txt）
# 输入：pathA/work/pa_data.rds（表达+生存）、work/estimate_matched.tsv、signature_registry.csv
# 输出：results/signature_scores.tsv, signature_correlations.tsv,
#       cox_models.tsv, attenuation.tsv, summary.tsv
# 纪律：不修约、不改数值；所有模型与协变量按预注册文件执行。
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(survival)
})

PB <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
PA <- "D:/WorkBuddy_Y02/out_ali_directions/pathA"
dir.create(file.path(PB, "results"), showWarnings = FALSE, recursive = TRUE)

message("[10] 读取数据 …")
pa   <- readRDS(file.path(PA, "work", "pa_data.rds"))
expr <- pa$expr_log
est  <- fread(file.path(PB, "work", "estimate_matched.tsv"))
reg  <- fread(file.path(PB, "signature_registry.csv"))
lu   <- as.data.table(pa$surv_luad)
ls_  <- as.data.table(pa$surv_lusc)

# ---- 统一协变量 ----
prep <- function(sv, est_) {
  d <- merge(sv, est_[, .(sample, Stromal, Immune, ESTIMATE, purity)],
             by = "sample", all.x = TRUE)
  d <- d[sample %in% colnames(expr)]              # 仅保留表达矩阵内样本（与 pathA 口径一致）
  d <- d[!is.na(time) & !is.na(status)]
  d[, stage_late := factor(stage_late, levels = c("early", "late"))]
  d[, gender := factor(gender)]
  d[, age := as.numeric(age)]
  d
}
d_luad <- prep(lu, est)
d_lusc <- prep(ls_, est)
message(sprintf("     LUAD n=%d（事件 %d）；LUSC n=%d（事件 %d）",
                nrow(d_luad), sum(d_luad$status), nrow(d_lusc), sum(d_lusc$status)))

# ---- 打分：z 分数等权平均 → 队列内标准化 ----
score_set <- function(d, genes) {
  g <- intersect(genes, rownames(expr))
  sub <- expr[g, d$sample, drop = FALSE]
  z <- t(scale(t(sub)))
  s <- colMeans(z, na.rm = TRUE)
  s <- as.numeric(scale(s))
  list(score = s, used = g, coverage = length(g) / length(genes))
}

fit_models <- function(d, s, tag) {
  dd <- copy(d); dd[, sig := s]
  cc1 <- dd[!is.na(sig) & !is.na(age) & !is.na(gender) & !is.na(stage_late)]
  out <- list()
  f0 <- coxph(Surv(time, status) ~ sig, data = dd)
  out$M0 <- data.table(model = "M0", n = f0$n, events = f0$nevent,
                       beta = coef(f0)["sig"], se = sqrt(vcov(f0)["sig", "sig"]))
  f1 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late, data = cc1)
  out$M1 <- data.table(model = "M1", n = f1$n, events = f1$nevent,
                       beta = coef(f1)["sig"], se = sqrt(vcov(f1)["sig", "sig"]))
  cc2 <- cc1[!is.na(purity)]
  f2 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + purity, data = cc2)
  out$M2 <- data.table(model = "M2", n = f2$n, events = f2$nevent,
                       beta = coef(f2)["sig"], se = sqrt(vcov(f2)["sig", "sig"]))
  cc3 <- cc1[!is.na(Immune) & !is.na(Stromal)]
  f3 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + Immune + Stromal, data = cc3)
  out$M3 <- data.table(model = "M3", n = f3$n, events = f3$nevent,
                       beta = coef(f3)["sig"], se = sqrt(vcov(f3)["sig", "sig"]))
  # PH 假定
  zph <- tryCatch(cox.zph(f3)$table["sig", "p"], error = function(e) NA_real_)
  res <- rbindlist(out)
  res[, `:=`(sig_id = tag,
             HR = exp(beta),
             lo = exp(beta - 1.959964 * se),
             hi = exp(beta + 1.959964 * se),
             p = 2 * pnorm(-abs(beta / se)),
             ph_p_M3 = zph)]
  res[]
}

score_rows <- list(); cor_rows <- list(); cox_rows <- list(); att_rows <- list()

for (i in seq_len(nrow(reg))) {
  r <- reg[i]
  genes <- trimws(strsplit(r$genes, ",")[[1]])
  sc <- score_set(d_luad, genes)
  score_rows[[i]] <- data.table(sig_id = r$sig_id, sample = d_luad$sample,
                               score = sc$score, coverage = sc$coverage,
                               genes_used = paste(sc$used, collapse = ","),
                               genes_missing = paste(setdiff(genes, rownames(expr)), collapse = ","))
  dd <- d_luad
  # 相关性（Spearman，全样本，不完整者剔除）
  for (v in c("Immune", "Stromal", "purity")) {
    ok <- !is.na(dd[[v]]) & !is.na(sc$score)
    ct <- suppressWarnings(cor.test(sc$score[ok], dd[[v]][ok], method = "spearman"))
    cor_rows[[length(cor_rows) + 1]] <- data.table(sig_id = r$sig_id, variable = v,
                                                   rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
  }
  cm <- fit_models(dd, sc$score, r$sig_id)
  cox_rows[[i]] <- cm
  b1 <- cm[model == "M1", beta]; b3 <- cm[model == "M3", beta]; b2 <- cm[model == "M2", beta]
  att_rows[[i]] <- data.table(
    sig_id = r$sig_id, stratum = r$stratum, n_genes = r$n_genes,
    coverage = sc$coverage,
    beta_M0 = cm[model == "M0", beta], beta_M1 = b1, beta_M2 = b2, beta_M3 = b3,
    att_purity = 1 - b2 / b1,          # M1 -> M2
    att_immune = 1 - b3 / b1,          # M1 -> M3（主估计量）
    p_M1 = cm[model == "M1", p], p_M2 = cm[model == "M2", p], p_M3 = cm[model == "M3", p],
    ph_p_M3 = cm[model == "M3", ph_p_M3]
  )
  message(sprintf("     %s %2d genes  cov=%.2f  M1 β=%.3f  M3 β=%.3f  衰减=%.1f%%",
                  r$sig_id, r$n_genes, sc$coverage, b1, b3, 100 * (1 - b3 / b1)))
}

scores <- rbindlist(score_rows); cors <- rbindlist(cor_rows)
coxs <- rbindlist(cox_rows);    att  <- rbindlist(att_rows)

# BH 校正（按终点分别；此处为 LUAD）
coxs[, q := p.adjust(p, method = "BH"), by = model]
att[, `:=`(q_M1 = p.adjust(p_M1, "BH"), q_M3 = p.adjust(p_M3, "BH"))]
# 剂量-反应：衰减 vs 签名分×Immune 的相关
imm <- cors[variable == "Immune", .(sig_id, rho_immune = rho)]
att <- merge(att, imm, by = "sig_id")
att[, endpoint := "TCGA-LUAD, OS"]

fwrite(scores, file.path(PB, "results", "signature_scores_luad.tsv"), sep = "\t")
fwrite(cors,   file.path(PB, "results", "signature_correlations_luad.tsv"), sep = "\t")
fwrite(coxs,   file.path(PB, "results", "cox_models_luad.tsv"), sep = "\t")
fwrite(att,    file.path(PB, "results", "attenuation_luad.tsv"), sep = "\t")

# ---- 汇总 ----
med_att   <- median(att$att_immune, na.rm = TRUE)
n_sig_M1  <- sum(att$q_M1 < 0.05, na.rm = TRUE)
n_sig_M3  <- sum(att$q_M3 < 0.05, na.rm = TRUE)
ct_dose   <- suppressWarnings(cor.test(att$rho_immune, att$att_immune, method = "spearman"))
sm <- data.table(
  endpoint = "TCGA-LUAD, OS", n_signatures = nrow(att),
  median_attenuation = med_att,
  median_attenuation_purity_only = median(att$att_purity, na.rm = TRUE),
  n_FDR_sig_M1 = n_sig_M1, n_FDR_sig_M3 = n_sig_M3,
  n_lost_significance = n_sig_M1 - n_sig_M3,
  dose_response_rho = unname(ct_dose$estimate), dose_response_p = ct_dose$p.value,
  median_rho_immune = median(att$rho_immune, na.rm = TRUE),
  median_rho_purity = median(cors[variable == "purity", rho], na.rm = TRUE)
)
fwrite(sm, file.path(PB, "results", "summary_luad.tsv"), sep = "\t")
print(sm)
message("[10] 完成")
