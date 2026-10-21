# 修订 04：LUSC 跨组织学复算（post-hoc 稳健性检验）

> 写定时间：2026-10-21 10:40，**早于 LUSC 任何新结果产出**。
> 目的：回应"主分析仅 LUAD"这一审稿必问点；检验三层图景是否跨组织学成立。
> 性质：post-hoc 稳健性检验，**不改动**主注册表、主估计量、签名名单与已锁定哈希（W1–W4 四波 sha256 不变）。

---

## 1. 要回答的问题（且只回答这些）

1. 在 LUSC 中，49 个已发表签名的 M1→M3 衰减是否同样**无系统性方向**（收缩/放大约各半、中位 A 接近 0）？
2. 三层（炎症/代谢/免疫）的**中位 A 排序**是否与 LUAD 一致？
3. 三条对照轴（淋巴细胞/中性粒/营养肝源）在 LUSC 上是否复现"已知免疫代理会崩、无关轴不动"的方法功效证据？

**不回答**：不检验任何签名在 LUSC 的预后价值（FDR 显著数只作描述）；不重新挑选签名；不因 LUSC 结果好坏调整主分析。

## 2. 数据与模型（与 LUAD 主分析逐项一致，不做任何改动）

| 项 | 定义 |
|---|---|
| 队列 | TCGA-LUSC，`pa$surv_lusc`（541 例），完整病例 = 有 OS、age、gender、stage_late、ImmuneScore、StromalScore 且样本在表达矩阵内 |
| 表达 | 同一 `pa$expr_log`（log2），同一 `score_set()`（基因内 z 分数 → 列均值 → 标准化） |
| 分期 | 沿用 `pathA` 的 `stage_late`（Stage I/II = early，III/IV = late；`[Discrepancy]` 视为缺失） |
| M1 | `Surv(time, status) ~ sig + age + gender + stage_late` |
| M3 | M1 + `Immune + Stromal`（ESTIMATE，来自 `work/estimate_matched.tsv`，与 LUAD 同一份匹配表） |
| 估计量 | A = 1 − β_M3/β_M1（log-HR 尺度），**与主分析完全相同** |
| 不确定性 | 按签名 cluster bootstrap 5,000 次（seed 20261020），与主分析相同 |
| FDR | BH，按 LUSC 内 49 个签名校正（与 LUAD 分开校正，不合并） |

## 3. 预声明的判读规则（先写死，结果出来照此判）

- **R1（稳健性成立）**：LUSC 中位 A 落在 ±10% 内、收缩/放大比介于 1:2 与 2:1、且三层排序与 LUAD 相同（免疫 > 代谢 ≈ 炎症）→ 主结论升级为"跨组织学稳健"。
- **R2（部分复现）**：方向计数或三层排序只有一项吻合 → 报"部分复现"，主结论限定为 LUAD，LUSC 作组织学异质性证据。
- **R3（不一致）**：LUSC 出现系统性方向（如收缩:放大 > 2:1 或 < 1:2）或三层排序反转 → <strong>照实报不一致，不解释为"失败"</strong>；讨论组织学差异与签名来源（多数签名文献为 LUAD）作为解释。
- **R4（无效队列）**：完整病例 < 300 或事件 < 100 → 停止，只报样本量不足，不做推测。
- 任何情况下：**不因 LUSC 结果修改或剔除签名**；三层各报 n；|A|>100% 者照旧只列不解读。

## 4. 附加（成本低、同批跑）

- 三条对照轴在 LUSC 重跑（同基因集、同 M1/M3 定义）→ 输出 `results/control_axes_lusc.tsv`，用于回答"管线在 LUSC 是否有功效"。
- 记录 LUSC 的 ESTIMATE 与 purity 相关（核对 M2 冗余是否同样成立）。

## 5. 输出

- `results/attenuation_lusc49.tsv`（49 行 × 全字段）
- `results/summary_lusc49.tsv`（ALL / 三层 / 对照）
- `results/control_axes_lusc.tsv`
- `logs/21_lusc_calibration.log`

## 6. 披露义务（无论结果如何都要写进正文）

1. LUSC 队列与 LUAD 队列**样本不重叠**，但同属 TCGA，人群与平台一致，不构成完全独立验证。
2. 49 个签名的文献来源以 LUAD 为主，跨组织学迁移性本身有限——因此本分析定性为**稳健性检验**，不是外部验证。
3. LUSC 的 `[Discrepancy]`（4 例）与 age 缺失（5 例）已排除，须在 Methods 注明。
