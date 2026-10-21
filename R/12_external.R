# =====================================================================
# 12_external.R — 样本外验证：13 个签名在 GSE68465 / GSE31210（均为 LUAD）中重算
#   目的：区分"纯度伪影"与"过拟合/队列特异性"——若是前者，样本外应仍显著且被校正削弱；
#         若是后者，样本外应普遍不显著。
#   口径：为跨平台可比，主对比指标统一用 M0（签名分每 SD 的单变量 Cox），不引入分期解析噪声。
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })

PB <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
PA <- "D:/WorkBuddy_Y02/out_ali_directions/pathA"
reg <- fread(file.path(PB, "signature_registry.csv"))

load_cohort <- function(dir, tag) {
  e <- fread(cmd = sprintf("gzip -dc %s/expr_symbol.tsv.gz", dir), header = TRUE)
  setnames(e, 1, "symbol")
  m <- as.matrix(e[, -1]); rownames(m) <- e$symbol
  cl <- fread(file.path(dir, "clinical.tsv"), skip = 1)
  cl <- cl[!is.na(os_time_days) & os_time_days > 0]
  cl <- cl[sample %in% colnames(m)]
  list(expr = m, clin = cl, tag = tag)
}

score_in <- function(co, genes) {
  g <- intersect(genes, rownames(co$expr))
  sub <- co$expr[g, co$clin$sample, drop = FALSE]
  z <- t(scale(t(sub)))
  s <- as.numeric(scale(colMeans(z, na.rm = TRUE)))
  list(score = s, coverage = length(g) / length(genes), used = g)
}

run <- function(co) {
  rows <- list()
  for (i in seq_len(nrow(reg))) {
    r <- reg[i]; genes <- trimws(strsplit(r$genes, ",")[[1]])
    sc <- score_in(co, genes)
    d <- copy(co$clin); d[, sig := sc$score]
    f <- coxph(Surv(os_time_days, os_event) ~ sig, data = d)
    rows[[i]] <- data.table(endpoint = co$tag, sig_id = r$sig_id, n = f$n, events = f$nevent,
      coverage = sc$coverage, missing_genes = paste(setdiff(genes, rownames(co$expr)), collapse = ","),
      beta = coef(f)["sig"], HR = exp(coef(f)["sig"]),
      lo = exp(confint(f)[1]), hi = exp(confint(f)[2]),
      p = summary(f)$coefficients["sig", "Pr(>|z|)"])
  }
  out <- rbindlist(rows)
  out[, q := p.adjust(p, "BH")]
  out
}

c1 <- load_cohort(file.path(PA, "data", "external", "gse68465"), "GSE68465 (LUAD, n=442)")
c2 <- load_cohort(file.path(PA, "data", "external", "gse31210"), "GSE31210 (LUAD, n=226)")
message(sprintf("%s: %d 基因 × %d 样本；%s: %d × %d",
                c1$tag, nrow(c1$expr), ncol(c1$expr), c2$tag, nrow(c2$expr), ncol(c2$expr)))

res <- rbind(run(c1), run(c2))
fwrite(res, file.path(PB, "results", "external_replication.tsv"), sep = "\t")

# 与 TCGA-LUAD 的 M0 合并，判断方向一致性
tcga <- fread(file.path(PB, "results", "cox_models_luad.tsv"))[model == "M0", .(sig_id, beta_tcga = beta, p_tcga = p)]
cmp <- merge(res, tcga, by = "sig_id")
cmp[, same_direction := sign(beta) == sign(beta_tcga)]
fwrite(cmp, file.path(PB, "results", "external_vs_tcga.tsv"), sep = "\t")

print(cmp[, .(endpoint, sig_id, coverage, beta_tcga, beta, HR, p, q, same_direction)])
summ <- cmp[, .(n = .N,
                n_same_direction = sum(same_direction),
                n_p05 = sum(p < 0.05),
                n_q05 = sum(q < 0.05)), by = endpoint]
print(summ)
fwrite(summ, file.path(PB, "results", "external_summary.tsv"), sep = "\t")
message("[12] 完成")
