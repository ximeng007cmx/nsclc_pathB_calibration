# -*- coding: utf-8 -*-
"""对 Tier A 候选逐条套用 I1/I3/E3/E4 的规则判定，并打印理由。

设计原则：判定必须**可由第三人复核**，所以用显式规则表，而不是逐篇"看着像就留"。
凡规则表覆盖不到的情形，一律标 REVIEW（人工判），不默认纳入。

构念边界（I3）的来源：00_criteria_preregistered.md 第 2 节 I3
  "候选基因池明确来自'炎症反应相关基因'或'营养/代谢（脂质、脂肪酸、糖、氨基酸）相关基因'"
因此：
  - 炎症类词根：inflammat*, cytokine, IL-, TNF, NLRP3 等 -> I3 PASS(炎症)
  - 营养/代谢类词根：lipid, fatty acid, amino acid, glucose, glycolys, cholesterol,
    glutamine, energy metabolism, iron metabolism, vitamin, nutrient, diet
    -> I3 PASS(代谢)
  - 明确**不属于**上述两类（ferroptosis / cuproptosis / disulfidptosis / pyroptosis /
    autophagy / m6A / hypoxia / ECM / CAF / lysosome / drug-tolerance / DNA repair ...）
    -> I3 FAIL（构念外）
  - 免疫相关（immune-related / immunogenic cell death / immune signature / MHC / T-cell）
    -> 单独规则 CIRCULAR：用免疫签名去校正免疫含量属循环论证，不纳入
"""
import csv, re, sys

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
rows_all = list(csv.DictReader(open(f"{PB}/work/candidates_wave2_screened.csv", encoding="utf-8-sig")))
rows = [r for r in rows_all
        if r["sig"] == "True" and r["anchor"] == "True" and r["lung"] == "True"
        and 2 <= int(r["n_in_abstract"]) <= 15]
print(f"输入 Tier A 候选：{len(rows)}（Tier B/其余 {len(rows_all) - len(rows)} 条不进入本轮）\n")

INF = re.compile(r"inflammat|cytokin|interleukin|\bIL-?\d|tumou?r necrosis|NLRP3|"
                 r"systemic immune|neutrophil|lymphocyte|macrophage|CRP|prognostic nutritional", re.I)
MET = re.compile(r"lipid|fatty acid|amino acid|glucose|glycoly|cholesterol|glutamine|"
                 r"energy metabolism|metabolic|metabolism|iron metabolism|vitamin|nutrient|"
                 r"dietary|propionate|butyrate|choline|nicotinamide|folic acid|nitrogen metabolism|"
                 r"ammonia|glutathione|oxidative phosphorylation|lactate|ketone|bile acid|"
                 r"polyamine|spermidine|taurine|serine|methionine|histidine|tryptophan|"
                 r"one-carbon|carbohydrate|triglyceride|sphingolipid|arachidonic", re.I)
OUT_OF_CONSTRUCT = re.compile(r"ferroptosis|cuproptosis|disulfidptosis|pyroptosis|autophagy|"
                              r"\bm6A\b|methylation|hypoxia|extracellular matrix|\bECM\b|"
                              r"fibroblast|CAF\b|lysosom|drug.toleran|persister|DNA repair|"
                              r"PANoptosis|apoptosis|necroptosis|senescence|ubiquitin|"
                              r"focal adhesion|proliferation|mitochondrial biogenesis|"
                              r"immune exclusion|M2 macrophage|TEAD4|RN7SK|pseudogene", re.I)
IMMUNE_CIRC = re.compile(r"immune.related|immunogenic cell death|immune signature|"
                         r"immune gene|immunotherapy response signature|MHC|HLA|T.cell|"
                         r"CD8|immune landscape|immuno", re.I)
LNCRNA = re.compile(r"lncRNA|long non.coding|miRNA|microRNA|circRNA|circular RNA|"
                    r"ceRNA|pseudogene|snRNA|snoRNA", re.I)


def judge(r):
    blob = (r["title"] + " " + r["abstract_gene_list"] + " " + r["n_gene_cue"])
    title = r["title"]
    reasons = []
    # I1：2020 起（检索式已保证，冗余确认）
    if not (str(r["year"]).isdigit() and int(r["year"]) >= 2020):
        return "EXCL", "E_I1", "发表年 <2020"
    # E3：非编码主体
    if LNCRNA.search(title):
        return "EXCL", "E3", f"标题含非编码词：{LNCRNA.search(title).group(0)}"
    # I3：构念外
    ooc = OUT_OF_CONSTRUCT.search(title)
    inf = INF.search(blob)
    met = MET.search(blob)
    if ooc and not inf:
        return "EXCL", "I3", f"构念外锚点：{ooc.group(0)}（非炎症、非营养/代谢）"
    # 循环论证
    if IMMUNE_CIRC.search(title) and not inf:
        return "EXCL", "CIRC", f"免疫相关签名：{IMMUNE_CIRC.search(title).group(0)}（循环论证）"
    if inf:
        return "INCL", "I3-inflammation", inf.group(0)
    if met:
        return "INCL", "I3-metabolism", met.group(0)
    return "REVIEW", "I3-unclear", "无炎症/代谢词根，需人工读摘要"


out = []
for r in rows:
    st, rule, why = judge(r)
    out.append({**{k: r[k] for k in ("pmcid", "pmid", "doi", "year", "journal", "title",
                                     "abstract_gene_list", "n_in_abstract", "n_gene_cue",
                                     "queries")},
                "decision": st, "rule": rule, "reason": why})
    print(f"{st:6s} {rule:18s} {r['pmcid']:14s} {r['year']} n={int(r['n_in_abstract']):2d}  "
          f"{r['title'][:58]:58s} | {why[:36]}")

with open(f"{PB}/work/wave2_rule_judgement.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)

from collections import Counter
c = Counter((o["decision"], o["rule"]) for o in out)
print("\n=== 汇总 ===")
for k, v in sorted(c.items(), key=lambda x: -x[1]):
    print(f"  {k[0]:7s} {k[1]:18s} {v}")
print(f"\n-> work/wave2_rule_judgement.csv")
