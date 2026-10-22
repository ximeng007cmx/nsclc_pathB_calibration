# Data and code availability statements (submission variants)

提交系统（Wiley / Cancer Medicine）里不同表单项要求长度不同，以下按用途分版，逐条可直接粘贴。
**核心事实**（与正文一致，勿改数字）：
- Zenodo DOI：`10.5281/zenodo.23537008`（记录页 https://zenodo.org/records/23537008 ，v1.0.0，cc-by-4.0）
- 代码仓库：https://github.com/ximeng007cmx/nsclc_pathB_calibration
- Zenodo 含两个文件：仓库快照 zip（346 条目）+ `pa_data.rds`（TCGA-LUAD/LUSC 表达矩阵与生存表，145,327,565 B）

---

## 1. 短版（表单 "Data availability statement"，约 55 词）

> The data that support the findings of this study are openly available in Zenodo at https://doi.org/10.5281/zenodo.23537008. The deposit contains the calibrated TCGA-LUAD/LUSC expression matrix and survival tables (pa_data.rds) together with a snapshot of the analysis code and all result tables. Source data are public (TCGA, GEO/ArrayExpress, PubMed Central).

## 2. 中版（Wiley "Data Availability Statement" 字段常用长度，约 110 词）

> All analysis code, the four hash-locked signature registries, all result tables and the calibrated TCGA-LUAD/LUSC expression matrix and survival tables are openly available in Zenodo at https://doi.org/10.5281/zenodo.23537008 (v1.0.0, CC-BY-4.0), mirrored as a code repository at https://github.com/ximeng007cmx/nsclc_pathB_calibration. SHA-256 checksums of every archived file are recorded in the repository (MANIFEST.sha256). Expression and clinical data derive from TCGA/TCGA-CDR, the external cohorts from GEO/ArrayExpress (E-GEOD-68465, E-GEOD-31210) and the full-text corpus from PubMed Central; these third-party data are subject to their original terms and are not redistributed beyond the derived files named above.

## 3. 长版（与正文 "Data and code availability" 一节完全一致，可直接引用）

> All analysis code (R and Python), the four hash-locked signature registries, all result tables and the manuscript source are archived at https://github.com/ximeng007cmx/nsclc_pathB_calibration (v1.0.0) and mirrored on Zenodo: 10.5281/zenodo.23537008. The Zenodo deposit additionally contains the calibrated TCGA-LUAD/LUSC expression matrix and survival tables (pa_data.rds) required to re-run the pipeline. SHA-256 checksums of all archived files are recorded in the repository (MANIFEST.sha256; docs/pathA_dependencies.md). Expression and clinical data originate from TCGA/TCGA-CDR, the external cohorts from GEO/ArrayExpress (E-GEOD-68465, E-GEOD-31210) and the full-text extraction corpus from PubMed Central; these third-party public data are subject to their original terms and are not redistributed beyond the derived files listed above.

## 4. Code availability（单独字段时用，约 60 词）

> All analysis code is openly available at https://github.com/ximeng007cmx/nsclc_pathB_calibration and archived at https://doi.org/10.5281/zenodo.23537008. The pipeline consists of 13 R scripts (signature scoring, Cox models M1 to M3, controls, external replication, endpoints, LUSC check) and 20 Python tools (registry building, screening and extraction, figures, tables, the 146-check reconciliation and the exact binomial check). Reproduction steps are given in the repository README.

## 5. 数据限制说明（如系统追问 "restrictions on data"）

> No restrictions apply to the derived files deposited in Zenodo. Raw third-party data are not redistributed: TCGA data are available through the Genomic Data Commons under its data use policy, and the external cohorts through ArrayExpress/BioStudies (E-GEOD-68465, E-GEOD-31210). Checksums of the raw array files used in this study are provided in the repository so that independently downloaded copies can be verified before re-analysis.

## 6. Ethics statement（Wiley 必填字段）

> This study is a secondary analysis of publicly available, de-identified data (TCGA, GEO/ArrayExpress, PubMed Central full texts). No human participants were recruited and no identifiable information was accessed, so ethics approval and informed consent were not required.

## 7. Funding statement（按实际情况二选一）

> This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

或（有资助时）：

> This work was supported by [funder] under grant [number].

## 8. Conflict of interest statement

> The author declares no competing interests.

## 9. Author contributions（单作者版；如有多作者请改为 CRediT 分项）

> Mingxiang Chen: conceptualisation, methodology, software, formal analysis, data curation, writing of the original draft, review and editing, and project administration.
