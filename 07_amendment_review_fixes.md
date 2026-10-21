# 修订 07 —— 投稿前评审（pre-submission review）驱动修正

- 日期：2026-10-21
- 触发：`manuscript/review_presubmission_2026-10-21.md`（1 CRITICAL + 4 MAJOR + 6 MINOR）
- 性质：纯文字与记号修正 + 1 项数值更正（C1）。**不改任何统计结果、不改数据、不重跑管线。**
- 正文对象：`manuscript/MS_draft_v1.md`（修订 06 版，哈希 d0711d03a3641be3 → 修订 07 新哈希见本文件末尾）

---

## C1（CRITICAL）3.2 节二项检验 P 值：更正并补计算记录

**问题**：正文原句 "the split was 17 shrinking versus 15 growing — a binomial test gives P ≈ 1."
① results/work/logs 检索 "binomial" 零命中，无计算记录，违反本稿 2.10 节自述纪律；
② 数值错误：17:15（n=32, p=0.5）的双侧精确 P = **0.8601**，不是 ≈1。

**处置**：
- 新增 `tools/28_binomial_test.py`（精确二项检验，R binom.test 口径：pmf ≤ pmf(obs) 求和），
  独立锚点 1 − C(32,16)/2³² 交叉断言，输出 `results/binomial_direction.tsv`
  （SHA-256 前 16 位 70c27bcf922c78f1），并与 `manuscript_reconciliation.tsv` 的 3.2 节
  两行（17 / 15）交叉核对通过。
- 正文改为 "…the split was 17 shrinking versus 15 growing (exact two-sided binomial P = 0.86)."

**对结论的影响**：无。原意（方向接近 1:1、无系统方向）不变，且 0.86 比 "≈1" 更准确地支持同一结论。

## M1 英文正文 em-dash 19 处 → 逐处改写（逗号/冒号/括号/分号）

| 行号（修订 06 版） | 原文片段 | 改法 |
|---|---|---|
| 29×2 | TME 插入语 — … — | 括号 "(TME; the fraction of …)" |
| 31 | "…whole set — a form of…" | 逗号 |
| 45 | "8 were never scored — six because…" | 冒号 |
| 47 | "…stated count — which is why…" | 逗号 |
| 103 | "…15 growing — a binomial test…" | 随 C1 一并消除 |
| 107 | "…significant) — the control pattern…" | 分号 |
| 111 | "which is the point — a small random panel…" | 随 M3 改写消除 |
| 119×2 | P2 点估计 — … — | 括号对 |
| 141 | 4.2 标题 "— and what that means" | 删除破折号尾 |
| 147×2 | 4.3 三分层插入语 | 括号对 |
| 149 | "…immune signatures — which is exactly…" | 逗号 |
| 153×2 | "quantity — TME content — accounts" | 逗号对 |
| 159×3 | post hoc 列举 —…—；"incomplete — a limitation" | 括号对 / 冒号 |

## M2 投稿文件中文残留 → 全部英化（信息保留，不删除）

- 头部引用说明（原 5–9 行）→ English provenance note（内容逐条对应）
- References 前中文注（原 171 行）→ English
- SIG01–SIG49 中文说明（原 183 行）→ English
- "图表清单" 整节（原 187–198 行）→ "Figures and tables"（English）
- 投稿转换时上述内部说明节仍须整体剥离（已在评审报告注明）
- `figure_captions.md` / `tables.md` 内的中文 provenance 注属内部文件，本轮不动

## M3 口语化 3 处 → 中性表述

1. 3.4 "…is itself large, which is the point — a small random panel…" → "…is itself large: a small random panel can appear to shrink or grow dramatically."
2. 4.2 "…than random sets, which is what one expects if…" → "…than random sets, as expected if…"
3. 4.3 "The stability of the ordering across cohorts while per-signature values correlate at approximately zero is itself the message: stratum-level behaviour…" → "The ordering is stable across cohorts even though per-signature values correlate at approximately zero; stratum-level behaviour is therefore a property of the construct, whereas per-signature attenuation is a property of the cohort."

## M4 "P 0.053 to 0.775" 记号歧义（摘要 + 3.3）→ 明确表述

- 两处 "(P 0.053 to 0.775)" → "(P changed from 0.053 to 0.775)"

## 顺带（MINOR 1.3，同句低风险）

- 3.2 "which rules out the 'unmasking' narrative" → "which argues against the 'unmasking' narrative"

---

## 对勾稽与锁定的影响

- 数字层面：新增 1 个有留档的 P 值（0.8601）；其余所有数值未动。
- `25_reconcile_manuscript.py` 不受影响（其 146 项不含二项检验；如需可将新 P 值纳入，非必须）。
- LOCK.txt 追加修订 07 条目。

## 修订 07 后正文哈希

- MS_draft_v1.md SHA-256 前 16 位：**67dc50ec12860396**（修订 06 版为 d0711d03a3641be3）
- 修订后验收：em-dash 全文 0；CJK 仅剩头部说明中的真实文件名（8 字符）；勾稽重跑 146 项 = 145 对 + 1 人工核对（与修订 06 后一致），数字层无回归
- 本修订文件自身 SHA-256 前 16 位：3c45909633ce70b6（登记于 LOCK.txt 时再核）
