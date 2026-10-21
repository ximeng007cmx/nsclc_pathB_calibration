# -*- coding: utf-8 -*-
"""修订 02 的第 3 层判定：对炎症定向检索的 Tier A 候选套用 INF-STRICT / CIRC / I3-OUT。

判定原则（与 08_classify_anchor.py 一致，但更严）：
  * 题名出现 CIRC 词（巨噬/中性粒/T/NK/免疫细胞浸润/ICD/免疫相关）-> 直接排除；
  * 题名出现构念外词 -> 直接排除；
  * 摘要里出现 CIRC/OUT 词但题名干净 -> 标 REVIEW，交人工读摘要；
  * 摘要无"炎症基因池"措辞 -> I3-INF 不过（INF-STRICT 的摘要代理）；
  * 其余 -> INCL，送全文抽取。
"""
import csv, re

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
rows = list(csv.DictReader(open(f"{PB}/work/wave3_screened.csv", encoding="utf-8-sig")))
rows = [r for r in rows if 2 <= int(r["n_in_abstract"]) <= 15]
print(f"Tier A 候选：{len(rows)}\n")

CIRC = re.compile(r"macrophage|neutrophil|\bNK\b|\bNK[- ]cell|\bT[- ]cell|immune cell|immunogen|"
                  r"immune[- ]related|immunotherapy|CD8|regulatory T|Treg", re.I)
OUT = re.compile(r"ferroptosis|cuproptosis|disulfidptosis|pyroptosis|autophagy|\bm6A\b|hypoxia|"
                 r"fibroblast|\bCAF\b|lysosom|DNA repair|necroptosis|PANoptosis|senescence|"
                 r"lncRNA|miRNA|microRNA|circRNA|single[- ]cell|exosome|mucin|glycosylation", re.I)
INFPOOL = re.compile(r"inflammation[- ]related genes?|inflammatory response[- ]?related genes?|"
                     r"\bIRGs?\b|\bDEIRGs?\b|inflammation[- ]associated genes?|inflammatory genes?|"
                     r"cytokine[- ]?related genes?|inflammat\w* gene set", re.I)

out = []
for r in rows:
    title = r["title"]; ab = r.get("abstract_gene_list", "")
    q = r.get("queries", "")
    c_t, o_t = CIRC.search(title), OUT.search(title)
    c_a, o_a = CIRC.search(ab), OUT.search(ab)
    if o_t:
        st, rule, why = "EXCL", "I3-OUT", f"题名构念外：{o_t.group(0)}"
    elif c_t:
        st, rule, why = "EXCL", "CIRC", f"题名免疫细胞类型/免疫签名：{c_t.group(0)}"
    elif False:
        pass
    elif o_a or c_a:
        st, rule, why = "REVIEW", "TXT-FLAG", f"摘要出现 {(o_a or c_a).group(0)}，需人工确认构念"
    else:
        has = bool(INFPOOL.search(title) or INFPOOL.search(ab) or r["inf_pool"] == "True")
        st = "INCL"
        rule = "I3-inflammation" if has else "I3-inflammation(全文待核)"
        why = ("摘要含炎症基因池措辞" if has else "摘要无炎症措辞，需全文核验基因池来源")
    out.append({**{k: r[k] for k in ("pmcid", "pmid", "doi", "year", "journal", "title",
                                     "abstract_gene_list", "n_in_abstract", "inf_pool_term", "queries")},
                "decision": st, "rule": rule, "reason": why})
    print(f"{st:6s} {rule:16s} {r['pmcid']:14s} {r['year']} n={int(r['n_in_abstract']):2d}  "
          f"{title[:60]:60s} | {why[:40]}")

with open(f"{PB}/work/wave3_rule_judgement.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)

from collections import Counter
print("\n=== 汇总 ===")
for k, v in Counter((o["decision"], o["rule"]) for o in out).most_common():
    print(f"  {k[0]:7s} {k[1]:16s} {v}")
print(f"\n-> work/wave3_rule_judgement.csv（{len(out)} 条）")
