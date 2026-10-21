# =====================================================================
# 00_dump.R — 从 pathA 的已校验中间件导出 pathB 所需的最小数据集
#   输入（只读）：pathA/work/pa_data.rds               表达矩阵 + 生存表
#                pathA/results/estimate_scores_matched.tsv  ESTIMATE/纯度（已匹配）
#   输出：pathB/work/genes.tsv, surv_luad.tsv, surv_lusc.tsv, estimate_matched.tsv
# 说明：不重算任何数值，仅做格式转换与落盘。
# =====================================================================
suppressPackageStartupMessages(library(data.table))

PA  <- "D:/WorkBuddy_Y02/out_ali_directions/pathA"
PB  <- "D:/WorkBuddy_Y02/out_ali_directions/pathB"
dir.create(file.path(PB, "work"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PB, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PB, "logs"), showWarnings = FALSE, recursive = TRUE)

message("[00] 读取 pa_data.rds …")
pa <- readRDS(file.path(PA, "work", "pa_data.rds"))
expr <- pa$expr_log
message(sprintf("     expr_log: %d 基因 × %d 样本", nrow(expr), ncol(expr)))
message(sprintf("     基因名示例: %s", paste(head(rownames(expr), 5), collapse = ", ")))

fwrite(data.table(symbol = rownames(expr)), file.path(PB, "work", "genes.tsv"), sep = "\t")

for (nm in c("surv_luad", "surv_lusc")) {
  s <- as.data.table(pa[[nm]])
  message(sprintf("     %s: %d 行，列 = %s", nm, nrow(s), paste(names(s), collapse = ",")))
  fwrite(s, file.path(PB, "work", paste0(nm, ".tsv")), sep = "\t")
}

est <- fread(file.path(PA, "results", "estimate_scores_matched.tsv"))
message(sprintf("     ESTIMATE 匹配表: %d 行；matched=TRUE %d 行",
                nrow(est), sum(est$matched)))
fwrite(est, file.path(PB, "work", "estimate_matched.tsv"), sep = "\t")

# 样本清单（cohort 归属）：从 surv 表与 expr 列名交集确认
lu <- as.data.table(pa$surv_luad); ls <- as.data.table(pa$surv_lusc)
idcol <- if ("sample" %in% names(lu)) "sample" else names(lu)[1]
cohort <- rbind(
  data.table(sample = lu[[idcol]], cohort = "TCGA-LUAD"),
  data.table(sample = ls[[idcol]], cohort = "TCGA-LUSC")
)
cohort[, in_expr := sample %in% colnames(expr)]
fwrite(cohort, file.path(PB, "work", "cohort_map.tsv"), sep = "\t")
message(sprintf("     队列映射: LUAD %d / LUSC %d；在表达矩阵内 %d",
                sum(cohort$cohort == "TCGA-LUAD"), sum(cohort$cohort == "TCGA-LUSC"),
                sum(cohort$in_expr)))
message("[00] 完成")
