# =====================================================================
# 23_lusc_calibration.R — LUSC 跨组织学复算（修订 04，post-hoc 稳健性检验）
#
# 纪律：与 20_calibration_all.R 逐项一致（同一表达矩阵、同一 score_set、同一
#       M1/M3 定义、同一估计量、同一 bootstrap 参数），仅把队列换成 LUSC。
#       不改签名名单、不改注册表、不改已锁定哈希。
# 修订文件：04_amendment_lusc.md（先于结果写定）
# 输出：results/attenuation_lusc49.tsv、summary_lusc49.tsv、control_axes_lusc.tsv
#       logs/21_lusc_calibration.log
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })

PB <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
PA <- "D:/WorkBuddy_Y02/out_ali_directions/pathA"

log_con <- file(file.path(PB, "logs", "21_lusc_calibration.log"), open = "wt")
sink(log_con, split = TRUE); on.exit({ sink(); close(log_con) }, add = TRUE)
cat("== 23_lusc_calibration.R（修订 04）==\n")

pa   <- readRDS(file.path(PA, "work", "pa_data.rds"))
expr <- pa$expr_log
est  <- fread(file.path(PB, "work", "estimate_matched.tsv"))

regfiles <- list.files(PB, pattern = "^signature_registry.*\\.csv$", full.names = TRUE)
regfiles <- regfiles[!grepl("excluded", regfiles)]
WAVE_OF <- c("signature_registry.csv" = "W1",
             "signature_registry_expansion.csv" = "W2",
             "signature_registry_expansion_w3.csv" = "W3",
             "signature_registry_expansion_w4.csv" = "W4")
keep <- c("sig_id", "pmcid", "pmid", "year", "journal", "stratum", "n_genes", "genes", "verification")
regs <- lapply(regfiles, function(f) {
  d <- fread(f); w <- WAVE_OF[[basename(f)]]; if (is.null(w)) w <- "W?"
  d2 <- d[, ..keep]; d2[, wave := w]; d2
})
reg <- rbindlist(regs, use.names = TRUE)
stopifnot(!anyDuplicated(reg$sig_id))
setorder(reg, wave, sig_id)
cat("签名注册：", nrow(reg), "个\n")

# ---- LUSC 队列（与 20 号 prep() 相同的清洗）----
sv <- as.data.table(pa$surv_lusc)
d <- merge(sv, est[, .(sample, Stromal, Immune, ESTIMATE, purity)], by = "sample", all.x = TRUE)
d <- d[sample %in% colnames(expr)][!is.na(time) & !is.na(status)]
d[, stage_late := factor(stage_late, levels = c("early", "late"))]
d[, gender := factor(gender)][, age := as.numeric(age)]
d <- d[!is.na(age) & !is.na(gender) & !is.na(stage_late) & !is.na(Immune) & !is.na(Stromal)]
cat(sprintf("LUSC 完整病例 n=%d，事件 %d（排除：stage 缺失/Discrepancy、age 缺失、无 ESTIMATE）\n",
            nrow(d), sum(d$status)))
stopifnot(nrow(d) >= 300, sum(d$status) >= 100)   # 修订 04 的 R4 停止条件

score_set <- function(dd, genes) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 2) return(list(score = rep(NA_real_, nrow(dd)), used = g,
                                 coverage = length(g) / length(genes)))
  sub <- expr[g, dd$sample, drop = FALSE]
  z <- t(scale(t(sub)))
  list(score = as.numeric(scale(colMeans(z, na.rm = TRUE))), used = g,
       coverage = length(g) / length(genes))
}

fit4 <- function(dd, s) {
  dd <- copy(dd); dd[, sig := s]
  cc <- dd[!is.na(sig) & !is.na(age) & !is.na(gender) & !is.na(stage_late)]
  f0 <- coxph(Surv(time, status) ~ sig, data = dd)
  f1 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late, data = cc)
  f2 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + purity, data = cc[!is.na(purity)])
  f3 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + Immune + Stromal,
              data = cc[!is.na(Immune) & !is.na(Stromal)])
  g <- function(f) {
    b <- coef(f)[["sig"]]; se <- sqrt(vcov(f)["sig", "sig"])
    list(beta = b, se = se, n = f$n, ev = f$nevent, p = 2 * pnorm(-abs(b / se)))
  }
  l <- list(M0 = g(f0), M1 = g(f1), M2 = g(f2), M3 = g(f3))
  ph <- tryCatch(cox.zph(f3)$table["sig", "p"], error = function(e) NA_real_)
  list(l = l, ph = ph)
}

out <- list(); cors <- list()
for (i in seq_len(nrow(reg))) {
  r <- reg[i]
  genes <- trimws(strsplit(r$genes, ",")[[1]])
  sc <- score_set(d, genes)
  ft <- fit4(d, sc$score)
  b1 <- ft$l$M1$beta; b3 <- ft$l$M3$beta; b2 <- ft$l$M2$beta
  out[[i]] <- data.table(
    sig_id = r$sig_id, wave = r$wave, stratum = r$stratum, n_genes = r$n_genes,
    coverage = sc$coverage, n_used = length(sc$used),
    beta_M0 = ft$l$M0$beta, beta_M1 = b1, beta_M2 = b2, beta_M3 = b3,
    se_M1 = ft$l$M1$se, se_M3 = ft$l$M3$se,
    att_purity = 1 - b2 / b1, att_immune = 1 - b3 / b1,
    p_M1 = ft$l$M1$p, p_M2 = ft$l$M2$p, p_M3 = ft$l$M3$p, ph_p_M3 = ft$ph,
    n = ft$l$M1$n, events = ft$l$M1$ev)
  ok <- !is.na(d$Immune) & !is.na(sc$score)
  ct <- suppressWarnings(cor.test(sc$score[ok], d$Immune[ok], method = "spearman"))
  cors[[i]] <- data.table(sig_id = r$sig_id, rho_immune = unname(ct$estimate))
  message(sprintf("  %-6s %-12s bM1=%+.3f (P=%.3g) bM3=%+.3f A=%+.1f%%",
                  r$sig_id, r$stratum, b1, out[[i]]$p_M1, b3, 100 * (1 - b3 / b1)))
}

att <- rbindlist(out)
att <- merge(att, rbindlist(cors), by = "sig_id")
att[, `:=`(q_M1 = p.adjust(p_M1, "BH"), q_M3 = p.adjust(p_M3, "BH"))]
att[, dir_flip := sign(beta_M1) != sign(beta_M3)]
att[, ratio_valid := p_M1 < 0.05]
fwrite(att, file.path(PB, "results", "attenuation_lusc49.tsv"), sep = "\t")
cat("逐签名输出：", nrow(att), "行\n")

boot_med <- function(x, B = 5000, seed = 20261020) {
  set.seed(seed); x <- x[is.finite(x)]
  if (length(x) < 3) return(c(median = median(x), lo = NA, hi = NA, n = length(x)))
  m <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  c(median = as.numeric(median(x)), lo = as.numeric(quantile(m, .025)),
    hi = as.numeric(quantile(m, .975)), n = length(x))
}

summ <- list()
add <- function(tag, sub) {
  s <- sub$att_immune * 100; b <- boot_med(s)
  bq <- boot_med(sub[ratio_valid == TRUE, att_immune * 100])
  summ[[length(summ) + 1]] <<- data.table(
    set = tag, n = nrow(sub), n_ratio_valid = nrow(sub[ratio_valid == TRUE]),
    median_A = b[["median"]], lo = b[["lo"]], hi = b[["hi"]],
    median_A_valid = bq[["median"]], lo_valid = bq[["lo"]], hi_valid = bq[["hi"]],
    n_abs_gt100 = sum(abs(s) > 100), n_dir_flip = sum(sub$dir_flip),
    n_FDR_M1 = sum(sub$q_M1 < 0.05), n_FDR_M3 = sum(sub$q_M3 < 0.05),
    median_abs_dbeta = median(abs(sub$beta_M3 - sub$beta_M1)))
}
add("ALL", att)
add("W1", att[wave == "W1"]); add("W2", att[wave == "W2"])
add("W3", att[wave == "W3"]); add("W4", att[wave == "W4"])
add("inflammation", att[stratum == "inflammation"])
add("metabolism", att[stratum == "metabolism"])
add("immune", att[stratum == "immune"])
sm <- rbindlist(summ)
fwrite(sm, file.path(PB, "results", "summary_lusc49.tsv"), sep = "\t")
print(sm)

# ---- 三条对照轴（管线功效的 LUSC 核对）----
AXES <- list(
  lymphocyte = c("CD3D","CD3E","CD3G","CD8A","CD8B","LCK","GZMB","GZMK","PRF1","IL7R",
                 "CCL5","CXCL10","CXCR6","CD27","CXCL13","HLA-DRA"),
  neutrophil = c("ELANE","MPO","LTF","CEACAM8","S100A8","S100A9","S100A12","MMP9","PADI4",
                 "CXCR2","FPR1","DEFA4","CAMP","LCN2","ARG1","OLFM4"),
  nutrition_liver = c("ALB","TF","RBP4","APOA1","APOA2","APOB","GC","TTR","HPX","SERPINA1",
                      "AHSG","ORM1","FGA","FGG","ITIH2","PLG"))
ctrl <- rbindlist(lapply(names(AXES), function(nm) {
  sc <- score_set(d, AXES[[nm]])
  ft <- fit4(d, sc$score)
  data.table(axis = nm, n_genes = length(AXES[[nm]]), coverage = sc$coverage,
             rho_immune = cor(sc$score, d$Immune, method = "spearman"),
             beta_M1 = ft$l$M1$beta, p_M1 = ft$l$M1$p,
             beta_M3 = ft$l$M3$beta, p_M3 = ft$l$M3$p,
             attenuation = 1 - ft$l$M3$beta / ft$l$M1$beta)
}))
fwrite(ctrl, file.path(PB, "results", "control_axes_lusc.tsv"), sep = "\t")
cat("\n三条对照轴（LUSC）：\n"); print(ctrl)
message("\n[23] 完成 -> results/attenuation_lusc49.tsv、summary_lusc49.tsv、control_axes_lusc.tsv")
