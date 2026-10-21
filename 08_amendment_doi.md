# 修订 08 —— Zenodo DOI 落稿（发布后回写）

- 日期：2026-10-22
- 触发：Zenodo deposit 23537008 发布成功（state: done，公开访问 HTTP 200）
- 性质：新增 "Data and code availability" 一节；无任何统计数字、结果或结论变动。

## 内容

**新增一节**（置于 5. Conclusion 与 References 之间）：

> ## Data and code availability
> All analysis code (R and Python), the four hash-locked signature registries, all result tables and the manuscript source are archived at https://github.com/ximeng007cmx/nsclc_pathB_calibration (v1.0.0) and mirrored on Zenodo: **10.5281/zenodo.23537008**. The Zenodo deposit additionally contains the calibrated TCGA-LUAD/LUSC expression matrix and survival tables (pa_data.rds) required to re-run the pipeline. SHA-256 checksums of all archived files are recorded in the repository (MANIFEST.sha256; docs/pathA_dependencies.md). Expression and clinical data originate from TCGA/TCGA-CDR, the external cohorts from GEO/ArrayExpress (E-GEOD-68465, E-GEOD-31210) and the full-text extraction corpus from PubMed Central; these third-party public data are subject to their original terms and are not redistributed beyond the derived files listed above.

**联动更新**：
- `README.md`（GitHub）：加 Zenodo DOI 徽章；Companion data 一节改为指向 DOI
- `tools/29_render_manuscript_html.py`：哈希断言更新为新版哈希
- `manuscript/MS_draft_v1.html` / `.pdf`：重渲染
- `MANIFEST.sha256`：重算（335+ 文件）
- `work/LOCK.txt`：追加发布与 DOI 记录

## Zenodo 侧事实（发布回执）

- DOI：**10.5281/zenodo.23537008**（conceptrecid 23537007）
- 记录页：https://zenodo.org/records/23537008 （state: done；未带 token 公开访问 200）
- 文件：nsclc_pathB_calibration-v1.0.0.zip（6,168,064 B，md5 023dc2c3…）+ pa_data.rds（145,327,565 B，md5 46c831c8…），两者与本地一致
- 作者：Chen, Mingxiang（单作者，用户指定）；许可 cc-by-4.0；v1.0.0

## 修订 08 后哈希

- MS_draft_v1.md SHA-256 前 16 位：（落盘后由 LOCK.txt 记录）
