# =====================================================================
# 20_calibration_all.R — 方案 B 扩容后主分析（13 原锁定 + 28 扩容）
#
# 预注册：00_criteria_preregistered.md（主估计量不变）
# 修订  ：01_amendment_expansion_to20.md（扩容为 post-hoc，三口径并列报告）
# 纪律  ：不做探索性重跑；W1/W2/ALL 三口径都报；分层必报。
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })

PB <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
PA <- "D:/WorkBuddy_Y02/out_ali_directions/pathA"
dir.create(file.path(PB, "results"), showWarnings = FALSE, recursive = TRUE)

pa   <- readRDS(file.path(PA, "work", "pa_data.rds"))
expr <- pa$expr_log
est  <- fread(file.path(PB, "work", "estimate_matched.tsv"))

# 自动读取全部纳入登记表（W1 主表 + 各波次扩容表），排除 *_excluded.csv
regfiles <- list.files(PB, pattern = "^signature_registry.*\\.csv$", full.names = TRUE)
regfiles <- regfiles[!grepl("excluded", regfiles)]
WAVE_OF <- c("signature_registry.csv" = "W1",
             "signature_registry_expansion.csv" = "W2",
             "signature_registry_expansion_w3.csv" = "W3",
             "signature_registry_expansion_w4.csv" = "W4")
keep <- c("sig_id", "pmcid", "pmid", "year", "journal", "stratum", "n_genes", "genes",
          "verification")
regs <- lapply(regfiles, function(f) {
  d <- fread(f)
  w <- WAVE_OF[[basename(f)]]
  if (is.null(w)) w <- "W?"
  d[, wave := w]
  d[, ..keep] -> d2
  d2[, wave := w]
  d2
})
reg <- rbindlist(regs, use.names = TRUE)
stopifnot(!anyDuplicated(reg$sig_id), !anyDuplicated(reg$pmcid))
setorder(reg, wave, sig_id)
message(sprintf("签名总数 %d（W1 %d / W2 %d；炎症 %d / 代谢 %d）",
                nrow(reg), sum(reg$wave == "W1"), sum(reg$wave == "W2"),
                sum(reg$stratum == "inflammation"), sum(reg$stratum == "metabolism")))

prep <- function(sv, est_) {
  d <- merge(sv, est_[, .(sample, Stromal, Immune, ESTIMATE, purity)], by = "sample", all.x = TRUE)
  d <- d[sample %in% colnames(expr)]
  d <- d[!is.na(time) & !is.na(status)]
  d[, stage_late := factor(stage_late, levels = c("early", "late"))]
  d[, gender := factor(gender)][, age := as.numeric(age)]
  d
}
d_luad <- prep(as.data.table(pa$surv_luad), est)
d_lusc <- prep(as.data.table(pa$surv_lusc), est)

score_set <- function(d, genes) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 2) return(list(score = rep(NA_real_, nrow(d)), used = g,
                                 coverage = length(g) / length(genes)))
  sub <- expr[g, d$sample, drop = FALSE]
  z <- t(scale(t(sub)))
  s <- as.numeric(scale(colMeans(z, na.rm = TRUE)))
  list(score = s, used = g, coverage = length(g) / length(genes))
}

fit4 <- function(d, s) {
  dd <- copy(d); dd[, sig := s]
  cc <- dd[!is.na(sig) & !is.na(age) & !is.na(gender) & !is.na(stage_late)]
  f0 <- coxph(Surv(time, status) ~ sig, data = dd)
  f1 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late, data = cc)
  f2 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + purity,
              data = cc[!is.na(purity)])
  f3 <- coxph(Surv(time, status) ~ sig + age + gender + stage_late + Immune + Stromal,
              data = cc[!is.na(Immune) & !is.na(Stromal)])
  g <- function(f, nm) {
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
  sc <- score_set(d_luad, genes)
  ft <- fit4(d_luad, sc$score)
  b1 <- ft$l$M1$beta; b3 <- ft$l$M3$beta; b2 <- ft$l$M2$beta
  rr <- data.table(
    sig_id = r$sig_id, wave = r$wave, stratum = r$stratum, n_genes = r$n_genes,
    coverage = sc$coverage, n_used = length(sc$used),
    beta_M0 = ft$l$M0$beta, beta_M1 = b1, beta_M2 = b2, beta_M3 = b3,
    se_M1 = ft$l$M1$se, se_M3 = ft$l$M3$se,
    att_purity = 1 - b2 / b1, att_immune = 1 - b3 / b1,
    p_M1 = ft$l$M1$p, p_M2 = ft$l$M2$p, p_M3 = ft$l$M3$p, ph_p_M3 = ft$ph,
    n = ft$l$M1$n, events = ft$l$M1$ev
  )
  for (v in c("Immune", "Stromal", "purity")) {
    ok <- !is.na(d_luad[[v]]) & !is.na(sc$score)
    ct <- suppressWarnings(cor.test(sc$score[ok], d_luad[[v]][ok], method = "spearman"))
    cors[[length(cors) + 1]] <- data.table(sig_id = r$sig_id, variable = v,
                                           rho = unname(ct$estimate), p = ct$p.value)
  }
  out[[i]] <- rr
  message(sprintf("  %-6s %s %-12s %2d g cov=%.2f  bM1=%+.3f (P=%.3g)  bM3=%+.3f  A=%+.1f%%",
                  r$sig_id, r$wave, r$stratum, r$n_genes, sc$coverage, b1, rr$p_M1, b3,
                  100 * rr$att_immune))
}

att <- rbindlist(out)
cors <- rbindlist(cors)
att <- merge(att, cors[variable == "Immune", .(sig_id, rho_immune = rho)], by = "sig_id")
att[, `:=`(q_M1 = p.adjust(p_M1, "BH"), q_M3 = p.adjust(p_M3, "BH"))]
att[, dir_flip := sign(beta_M1) != sign(beta_M3)]
att[, ratio_valid := p_M1 < 0.05]
fwrite(att, file.path(PB, "results", "attenuation_all.tsv"), sep = "\t")
fwrite(cors, file.path(PB, "results", "correlations_all.tsv"), sep = "\t")

# ---------- bootstrap：按签名重抽样（cluster bootstrap） ----------
boot_med <- function(x, B = 5000, seed = 20261020) {
  set.seed(seed)
  x <- x[is.finite(x)]
  if (length(x) < 3) return(c(median = median(x), lo = NA, hi = NA, n = length(x)))
  m <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  c(median = as.numeric(median(x)), lo = as.numeric(quantile(m, .025)),
    hi = as.numeric(quantile(m, .975)), n = length(x))
}

summ <- list()
add <- function(tag, sub) {
  s <- sub$att_immune * 100
  b <- boot_med(s)
  q <- sub[ratio_valid == TRUE, att_immune * 100]
  bq <- boot_med(q)
  summ[[length(summ) + 1]] <<- data.table(
    set = tag, n = nrow(sub), n_ratio_valid = nrow(sub[ratio_valid == TRUE]),
    median_A = b[["median"]], lo = b[["lo"]], hi = b[["hi"]],
    median_A_valid = bq[["median"]], lo_valid = bq[["lo"]], hi_valid = bq[["hi"]],
    n_abs_gt100 = sum(abs(s) > 100), n_dir_flip = sum(sub$dir_flip),
    n_FDR_M1 = sum(sub$q_M1 < 0.05), n_FDR_M3 = sum(sub$q_M3 < 0.05),
    median_abs_dbeta = median(abs(sub$beta_M3 - sub$beta_M1))
  )
}
add("ALL", att)
add("W1", att[wave == "W1"])
add("W2", att[wave == "W2"])
add("W3", att[wave == "W3"])
add("inflammation", att[stratum == "inflammation"])
add("metabolism", att[stratum == "metabolism"])
add("W1-inflammation", att[wave == "W1" & stratum == "inflammation"])
add("W1-metabolism", att[wave == "W1" & stratum == "metabolism"])
add("W2-inflammation", att[wave == "W2" & stratum == "inflammation"])
add("W2-metabolism", att[wave == "W2" & stratum == "metabolism"])
add("W3-inflammation", att[wave == "W3" & stratum == "inflammation"])
add("W4", att[wave == "W4"])
add("immune", att[stratum == "immune"])
# 分层差值 bootstrap（修订 02 的 P1 + 修订 03 的 P2）
boot_diff <- function(x, y, B = 5000, seed = 20261020) {
  set.seed(seed)
  d <- replicate(B, median(sample(x, length(x), TRUE)) - median(sample(y, length(y), TRUE)))
  c(est = as.numeric(median(x) - median(y)),
    lo = as.numeric(quantile(d, .025)), hi = as.numeric(quantile(d, .975)))
}
pair_rows <- list()
mk <- function(a, b, tag) {
  xa <- att[stratum == a, att_immune]; xb <- att[stratum == b, att_immune]
  if (length(xa) >= 3 && length(xb) >= 3) {
    r <- boot_diff(xa, xb)
    pair_rows[[length(pair_rows) + 1]] <<- data.table(
      comparison = tag, est = r[["est"]], lo = r[["lo"]], hi = r[["hi"]],
      n_a = length(xa), n_b = length(xb))
  }
}
mk("inflammation", "metabolism", "P1: inflammation minus metabolism")
mk("immune", "inflammation", "P2a: immune minus inflammation")
mk("immune", "metabolism", "P2b: immune minus metabolism")
if (length(pair_rows)) fwrite(rbindlist(pair_rows),
                              file.path(PB, "results", "stratum_diff_boot.tsv"), sep = "\t")

sm <- rbindlist(summ)
fwrite(sm, file.path(PB, "results", "summary_all.tsv"), sep = "\t")
print(sm)
message("[20] 完成")
