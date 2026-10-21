# =====================================================================
# 24_endpoint_sensitivity.R — 终点敏感性（PFI / DSS vs OS 主分析）
#
# 修订文件：05_amendment_endpoints.md（先于结果写定，E1–E3 判读规则）
# 纪律：与 20_calibration_all.R 逐项一致（同一队列、同一得分、同一 M1/M3、
#       同一估计量、同一 bootstrap 参数），仅把 Surv(time,status) 换成各终点。
# 输出：results/attenuation_endpoints.tsv、summary_endpoints.tsv
#       logs/22_endpoint_sensitivity.log
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })

PB <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
PA <- "D:/WorkBuddy_Y02/out_ali_directions/pathA"
log_con <- file(file.path(PB, "logs", "22_endpoint_sensitivity.log"), open = "wt")
sink(log_con, split = TRUE); on.exit({ sink(); close(log_con) }, add = TRUE)
cat("== 24_endpoint_sensitivity.R（修订 05）==\n")

pa <- readRDS(file.path(PA, "work", "pa_data.rds"))
expr <- pa$expr_log
est <- fread(file.path(PB, "work", "estimate_matched.tsv"))

regfiles <- list.files(PB, pattern = "^signature_registry.*\\.csv$", full.names = TRUE)
regfiles <- regfiles[!grepl("excluded", regfiles)]
WAVE_OF <- c("signature_registry.csv" = "W1",
             "signature_registry_expansion.csv" = "W2",
             "signature_registry_expansion_w3.csv" = "W3",
             "signature_registry_expansion_w4.csv" = "W4")
keep <- c("sig_id", "stratum", "n_genes", "genes")
regs <- lapply(regfiles, function(f) {
  d <- fread(f); w <- WAVE_OF[[basename(f)]]; if (is.null(w)) w <- "W?"
  d2 <- d[, ..keep]; d2[, wave := w]; d2
})
reg <- rbindlist(regs, use.names = TRUE)
setorder(reg, wave, sig_id)
cat("签名：", nrow(reg), "个\n")

base <- as.data.table(pa$surv_luad)
base <- merge(base, est[, .(sample, Stromal, Immune)], by = "sample", all.x = TRUE)
base <- base[sample %in% colnames(expr)]
base[, stage_late := factor(stage_late, levels = c("early", "late"))]
base[, gender := factor(gender)][, age := as.numeric(age)]

score_set <- function(dd, genes) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 2) return(rep(NA_real_, nrow(dd)))
  sub <- expr[g, dd$sample, drop = FALSE]
  as.numeric(scale(colMeans(t(scale(t(sub))), na.rm = TRUE)))
}

boot_med <- function(x, B = 5000, seed = 20261020) {
  set.seed(seed); x <- x[is.finite(x)]
  if (length(x) < 3) return(c(median = median(x), lo = NA, hi = NA, n = length(x)))
  m <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  c(median = as.numeric(median(x)), lo = as.numeric(quantile(m, .025)),
    hi = as.numeric(quantile(m, .975)), n = length(x))
}

ENDPOINTS <- c("OS", "PFI", "DSS")
all_rows <- list()
for (ep in ENDPOINTS) {
  d <- copy(base)
  d[, `:=`(time = as.numeric(get(paste0(ep, "_time"))),
           status = as.numeric(get(ep)))]
  d <- d[!is.na(time) & !is.na(status) & time > 0 &
           !is.na(age) & !is.na(gender) & !is.na(stage_late)]
  cat(sprintf("\n%s: n=%d，事件 %d\n", ep, nrow(d), sum(d$status)))
  if (nrow(d) < 300 || sum(d$status) < 100) { cat("  样本/事件不足，跳过\n"); next }

  for (i in seq_len(nrow(reg))) {
    r <- reg[i]
    s <- score_set(d, trimws(strsplit(r$genes, ",")[[1]]))
    if (all(is.na(s))) next
    dd <- copy(d); dd[, sig := s]
    # 与 20_calibration_all.R 完全一致：M1 不要求 Immune/Stromal 非缺失，M3 才要求
    cc1 <- dd[!is.na(sig)]
    cc3 <- dd[!is.na(sig) & !is.na(Immune) & !is.na(Stromal)]
    f1 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late, data = cc1)
    f3 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + Immune + Stromal, data = cc3)
    b1 <- coef(f1)[["sig"]]; b3 <- coef(f3)[["sig"]]
    p1 <- summary(f1)$coefficients["sig", "Pr(>|z|)"]
    p3 <- summary(f3)$coefficients["sig", "Pr(>|z|)"]
    all_rows[[length(all_rows) + 1]] <- data.table(
      endpoint = ep, sig_id = r$sig_id, wave = r$wave, stratum = r$stratum,
      n = f1$n, events = f1$nevent,
      beta_M1 = b1, beta_M3 = b3, p_M1 = p1, p_M3 = p3,
      A = 1 - b3 / b1, dir_flip = sign(b1) != sign(b3))
  }
}

out <- rbindlist(all_rows)
out[, q_M1 := p.adjust(p_M1, "BH"), by = endpoint]
out[, q_M3 := p.adjust(p_M3, "BH"), by = endpoint]
fwrite(out, file.path(PB, "results", "attenuation_endpoints.tsv"), sep = "\t")

summ <- list()
for (ep in unique(out$endpoint)) {
  sub <- out[endpoint == ep]
  for (st in c("ALL", "inflammation", "metabolism", "immune")) {
    ss <- if (st == "ALL") sub else sub[stratum == st]
    b <- boot_med(ss$A * 100)
    summ[[length(summ) + 1]] <- data.table(
      endpoint = ep, set = st, n = nrow(ss),
      n_ratio_valid = sum(ss$p_M1 < 0.05),
      median_A = b[["median"]], lo = b[["lo"]], hi = b[["hi"]],
      n_pos = sum(ss$A > 0), n_neg = sum(ss$A < 0), n_flip = sum(ss$dir_flip),
      median_abs_dbeta = median(abs(ss$beta_M3 - ss$beta_M1)),
      n_FDR_M1 = sum(ss$q_M1 < 0.05), n_FDR_M3 = sum(ss$q_M3 < 0.05))
  }
}
sm <- rbindlist(summ)
fwrite(sm, file.path(PB, "results", "summary_endpoints.tsv"), sep = "\t")
print(sm)

# 与 OS 的 A 值一致性
w <- dcast(out[, .(sig_id, endpoint, A)], sig_id ~ endpoint, value.var = "A")
for (ep in c("PFI", "DSS")) {
  if (all(c("OS", ep) %in% names(w))) {
    ok <- is.finite(w$OS) & is.finite(w[[ep]])
    ct <- suppressWarnings(cor.test(w$OS[ok], w[[ep]][ok], method = "spearman"))
    cat(sprintf("\nA(OS) vs A(%s): Spearman %.3f (P=%.3g, n=%d)\n",
                ep, unname(ct$estimate), ct$p.value, sum(ok)))
  }
}
message("\n[24] 完成 -> results/attenuation_endpoints.tsv、summary_endpoints.tsv")
