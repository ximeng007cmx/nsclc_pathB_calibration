# =====================================================================
# 22_external_m1m3.R — 外部队列（GSE68465 / GSE31210）的 M1→M3 复算
#
# 性质：预览性分析（用户要求先看结果再决定是否纳入论文）。
#       未预注册；若纳入论文，必须作为 post-hoc 外部复算披露，且主注册表不改。
#
# 数据来源（复用 pathA 已构建并归档的外部队列）：
#   pathA/data/external/gse68465/  E-GEOD-68465（Director's Challenge LUAD, HG-U133A）442 例/236 事件
#   pathA/data/external/gse31210/  E-GEOD-31210（日本 LUAD IA–II, HG-U133 Plus 2）226 例/35 事件
#   表达矩阵 = log2(x+1) symbol×样本（tools/build_external_matrix.py，max-mean 探针合并）
#
# ESTIMATE 复算说明（官方包不可得，处理方式）：
#   estimate 包 Bioconductor 3.23 无、mdanderson 不再提供下载、GitHub 不可达。
#   官方签名集取自原始论文补充材料（Yoshihara et al. 2013 Nat Commun,
#   Supplementary Data 1, ESM_488.xlsx 'signatures' 表）：
#     Stromal141_UP (141 基因), Immune141_UP (141 基因)
#   用 GSVA::ssgsea 对两集合算单样本 ES 作为 StromalScore/ImmuneScore。
#   数学事实：Cox 模型中协变量的线性缩放不改变其余系数，故 sig 系数
#   （即 A = 1 - beta_M3/beta_M1）与官方分数严格等价；只在分数绝对量级上与官方差常数。
#
# 分期映射（先声明）：
#   GSE68465：stage 为 pN+pT 拼接 → late = (pT3|pT4|pN1|pN2)；pNX/pp → NA
#   GSE31210：I 期(IA/IB) = early；II = late（队列本身只有 IA–II，无 III/IV）
# M1 协变量按队列可得性：
#   GSE68465：sig + age + sex + stage_late
#   GSE31210：sig + sex + stage_late（SDRF 中 age 全部 "not specified"，无法纳入，须披露）
#
# 输出：
#   results/external_m1m3_long.tsv      队列 × 49 签名逐行
#   results/external_m1m3_summary.tsv   队列 × wave/stratum 中位 A 汇总
#   results/external_vs_tcga_rho.tsv    外部 A 与 TCGA A 的 Spearman 相关
#   logs/20_external_m1m3.log
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })
PB <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
PA <- "D:/WorkBuddy_Y02/out_ali_directions/pathA"

log_con <- file(file.path(PB, "logs", "20_external_m1m3.log"), open = "wt")
sink(log_con, split = TRUE); on.exit({ sink(); close(log_con) }, add = TRUE)
cat("== 22_external_m1m3.R ==\n")

# ---- 签名注册表（与 20_calibration_all.R 同一逻辑）----
regfiles <- list.files(PB, pattern = "^signature_registry.*\\.csv$", full.names = TRUE)
regfiles <- regfiles[!grepl("excluded", regfiles)]
WAVE_OF <- c("signature_registry.csv" = "W1",
             "signature_registry_expansion.csv" = "W2",
             "signature_registry_expansion_w3.csv" = "W3",
             "signature_registry_expansion_w4.csv" = "W4")
reg <- rbindlist(lapply(regfiles, function(x) {
  d <- fread(x); w <- WAVE_OF[[basename(x)]]; if (is.null(w)) w <- "W?"
  d[, wave := w]; d
}), use.names = TRUE, fill = TRUE)
reg <- reg[!duplicated(sig_id)]
stopifnot(!anyDuplicated(reg$sig_id))
cat("签名注册：", nrow(reg), "个（", paste(reg[, .N, by = wave][, paste0(wave, ":", N)], collapse = " "), "）\n")

# ---- ESTIMATE 官方签名集（Yoshihara 2013 Supplementary Data 1）----
sig_xlsx <- file.path(PB, "work", "estimate_ext", "ESM_488.xlsx")
stopifnot(file.exists(sig_xlsx))
# 读 xlsx 的 'signatures' 表：set 列 + 横排基因（用 python 已核验 141+141）
# R 侧用 openxlsx 不一定装了；改为读 python 预转的 tsv（下方 python 步骤先行）
est_sets <- fread(file.path(PB, "work", "estimate_ext", "estimate_signatures.tsv"))
stopifnot(all(c("Stromal141_UP", "Immune141_UP") %in% est_sets$set))
SIGS <- list(stromal = est_sets[set == "Stromal141_UP"]$gene,
             immune  = est_sets[set == "Immune141_UP"]$gene)
cat("ESTIMATE 签名集：stromal", length(SIGS$stromal), "/ immune", length(SIGS$immune), "基因\n")

score_of <- function(genes, expr, samples) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 1) return(NULL)
  m <- expr[g, samples, drop = FALSE]
  gm <- t(scale(t(m)))                    # 基因内标准化（与 TCGA 管线一致）
  as.numeric(scale(colMeans(gm, na.rm = TRUE)))
}

COHORTS <- list(
  list(id = "GSE68465", dir = "gse68465",
       m1_covars = c("age", "sex", "stage_late"), has_age = TRUE),
  list(id = "GSE31210", dir = "gse31210",
       m1_covars = c("sex", "stage_late"), has_age = FALSE)
)

long <- list()
for (co in COHORTS) {
  base <- file.path(PA, "data", "external", co$dir)
  ex <- as.matrix(fread(file.path(base, "expr_symbol.tsv.gz")), rownames = 1)
  cl <- fread(file.path(base, "clinical.tsv"), skip = 1)
  cl <- cl[!is.na(os_time_days) & os_time_days > 0]
  keep <- intersect(cl$sample, colnames(ex))
  cl <- cl[match(keep, sample)]
  ex <- ex[, keep, drop = FALSE]
  cl[, sex := factor(sex)]
  # 分期映射
  if (co$id == "GSE68465") {
    cl[, stage_late := fcase(grepl("pT3|pT4|pN1|pN2", stage), "late",
                             grepl("pT1|pT2", stage) & grepl("pN0", stage), "early",
                             default = NA_character_)]
  } else {
    cl[, stage_late := fcase(stage %in% c("II", "IIA", "IIB"), "late",
                             stage %in% c("IA", "IB"), "early",
                             default = NA_character_)]
  }
  cl[, stage_late := factor(stage_late, levels = c("early", "late"))]
  cl <- cl[!is.na(stage_late) & !is.na(sex)]
  if (co$has_age) {
    cl[, age := suppressWarnings(as.numeric(age))]
    cl <- cl[!is.na(age)]
  }
  keep <- cl$sample; ex <- ex[, keep, drop = FALSE]
  cat(sprintf("\n%s: %d 基因 x %d 样本（%d 事件）；stage late=%d early=%d\n",
              co$id, nrow(ex), ncol(ex), sum(cl$os_event == 1),
              sum(cl$stage_late == "late"), sum(cl$stage_late == "early")))

  # ESTIMATE 分数
  es <- GSVA::gsva(GSVA::ssgseaParam(
    exprData = ex, geneSets = list(StromalScore = SIGS$stromal, ImmuneScore = SIGS$immune)),
    verbose = FALSE)
  est <- data.table(sample = colnames(es), Stromal = es["StromalScore", ], Immune = es["ImmuneScore", ])
  d <- merge(cl, est, by = "sample")
  d <- d[order(d$sample)]; ex <- ex[, d$sample, drop = FALSE]
  cat(sprintf("ESTIMATE: stromal 中位 %.0f，immune 中位 %.0f（ES x 量级，与官方差常数倍，不影响 A）\n",
              median(d$Stromal), median(d$Immune)))

  for (i in seq_len(nrow(reg))) {
    r <- reg[i]
    s <- score_of(unlist(strsplit(r$genes, ",\\s*")), ex, d$sample)
    if (is.null(s)) { long[[length(long)+1]] <- data.table(
        cohort = co$id, sig_id = r$sig_id, wave = r$wave, stratum = r$stratum,
        n_genes = r$n_genes, n_used = 0L, coverage = 0,
        beta_M1 = NA_real_, p_M1 = NA_real_, beta_M3 = NA_real_, p_M3 = NA_real_,
        A = NA_real_, ratio_valid = FALSE, note = "no gene hits"); next }
    dd <- copy(d); dd[, sig := s]
    f1 <- do.call(coxph, list(formula = as.formula(paste(
      "Surv(os_time_days, os_event) ~ sig +", paste(co$m1_covars, collapse = " + "))), data = dd))
    f3 <- do.call(coxph, list(formula = as.formula(paste(
      "Surv(os_time_days, os_event) ~ sig +", paste(co$m1_covars, collapse = " + "),
      "+ Immune + Stromal")), data = dd))
    b1 <- coef(f1)["sig"]; b3 <- coef(f3)["sig"]
    p1 <- summary(f1)$coefficients["sig", "Pr(>|z|)"]
    p3 <- summary(f3)$coefficients["sig", "Pr(>|z|)"]
    long[[length(long)+1]] <- data.table(
      cohort = co$id, sig_id = r$sig_id, wave = r$wave, stratum = r$stratum,
      n_genes = r$n_genes, n_used = length(intersect(unlist(strsplit(r$genes, ",\\s*")), rownames(ex))),
      coverage = length(intersect(unlist(strsplit(r$genes, ",\\s*")), rownames(ex))) / r$n_genes,
      beta_M1 = b1, p_M1 = p1, beta_M3 = b3, p_M3 = p3,
      A = 1 - b3 / b1,
      ratio_valid = (p1 < 0.05) && (b1 != 0), note = "")
  }
}

out <- rbindlist(long)
fwrite(out, file.path(PB, "results", "external_m1m3_long.tsv"), sep = "\t")
cat("\n逐行输出：", nrow(out), "行\n")

# ---- 汇总（按队列 x wave / stratum）----
summ <- out[, .(
  n = sum(!is.na(A)),
  median_A = if (sum(!is.na(A)) >= 3) median(A, na.rm = TRUE) else NA_real_,
  n_pos = sum(A > 0, na.rm = TRUE),
  n_neg = sum(A < 0, na.rm = TRUE),
  n_valid = sum(ratio_valid, na.rm = TRUE),
  median_A_valid = if (sum(ratio_valid) >= 3) median(A[ratio_valid], na.rm = TRUE) else NA_real_,
  n_sig_M1 = sum(p_M1 < 0.05, na.rm = TRUE),
  n_sig_M3 = sum(p_M3 < 0.05, na.rm = TRUE)
), by = .(cohort, stratum)]
fwrite(summ, file.path(PB, "results", "external_m1m3_summary.tsv"), sep = "\t")
print(summ)

# ---- 与 TCGA 的 A 对比 ----
tcga <- fread(file.path(PB, "results", "attenuation_all.tsv"))[, .(sig_id, A_tcga = att_immune)]
cmp <- merge(out[!is.na(A), .(cohort, sig_id, A_ext = A)], tcga, by = "sig_id")
rho_tab <- cmp[, {
  ct <- suppressWarnings(cor.test(A_ext, A_tcga, method = "spearman"))
  .(n = .N, rho = unname(ct$estimate), p = ct$p.value)
}, by = cohort]
fwrite(rho_tab, file.path(PB, "results", "external_vs_tcga_rho.tsv"), sep = "\t")
cat("\n外部 A 与 TCGA A 的 Spearman：\n"); print(rho_tab)

message("\n[22] 完成 -> results/external_m1m3_*.tsv、external_vs_tcga_rho.tsv")
