# -*- coding: utf-8 -*-
"""构建炎症波次（W3）登记表。

人工定稿原则：INF-STRICT 要求正文明确写出**候选基因池来自炎症/趋化因子基因集**。
正则匹配到"KEGG 富集结果里出现 cytokine 通路"这类句子**不算**——那是结果注释，不是基因池来源。
本脚本把 8 条自动通过者逐条复核后的结论写死，并记录理由。
"""
import csv, hashlib, json, os

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
meta = {r["pmcid"]: r for r in
        csv.DictReader(open(f"{PB}/work/wave3_screened.csv", encoding="utf-8-sig"))}
universe = set(l.strip() for l in
               open(os.path.join(PB, "work", "genes.tsv"), encoding="utf-8").read().splitlines()[1:])
ext = {x["pmcid"]: x for x in json.load(open(f"{PB}/work/extraction_wave3.json", encoding="utf-8"))}

DECISIONS = [
    ("PMC12475510", "inflammation", 7,
     "CCR2,CCR4,CCR6,XCR1,CCL20,CXCL17,XCL2",
     "First, 45 chemokine genes and 24 chemokine receptor genes (Supplementary Table 2) were retrieved from the literature.",
     "INCL", "趋化因子系统基因池（45 趋化因子 + 24 受体，文献汇总）→ 明确炎症构念；7 基因签名，TCGA-LUAD"),
    ("PMC9575749", "inflammation", 7,
     "SHC1,PRKCD,VEGFC,RPS6KA1,CAT,CDC25C,GPI",
     "A novel gene signature derived from the CXC subfamily of chemokine receptors predicts the prognosis and immune infiltration of patients with lung adenocarcinoma.",
     "INCL", "以 CXC 趋化因子受体亚家族为候选池构建的 7 基因签名；炎症构念（趋化因子受体）"),
]

EXCLUDED = [
    ("PMC7802904", "I3", "锚点为肿瘤微环境（TME），非炎症基因集；来源句只是 KEGG 通路注释，且自述 6 基因与抽到 4 个不符"),
    ("PMC10319060", "I3", "目的在于找诊断标志物 + 免疫浸润景观；来源句是 GO/KEGG 富集注释，非基因池来源"),
    ("PMC10800788", "I3-OUT", "黏着斑（focal adhesion）构念，I3 失败；且自述 10 基因与抽到 13 个不符"),
    ("PMC12187644", "CIRC", "PD-L1(CD274) 相关签名，属免疫检查点/免疫构念"),
    ("PMC12541842", "I3-OUT", "奥希替尼耐药相关签名，构念外"),
    ("PMC11985211", "I3-OUT", "DNA 复制相关签名，构念外"),
]

rows = []
for i, (pmcid, stratum, n, genes, ev, verdict, reason) in enumerate(DECISIONS, start=42):
    gl = [g.strip() for g in genes.split(",")]
    x = ext.get(pmcid, {})
    assert len(gl) == n, f"{pmcid} 声称 {n} 实际 {len(gl)}"
    assert 2 <= n <= 15
    assert x.get("n_extracted") == n, f"{pmcid} 抽取数 {x.get('n_extracted')} != {n}"
    miss = [g for g in gl if g not in universe]
    m = meta.get(pmcid, {})
    rows.append(dict(
        sig_id=f"SIG{i:02d}", pmcid=pmcid, pmid=m.get("pmid", ""), doi=m.get("doi", ""),
        journal=m.get("journal", ""), year=m.get("year", ""), stratum=stratum,
        n_genes=n, genes=",".join(gl), n_listed=n, verification="body_provenance",
        evidence=ev, title=m.get("title", ""), wave="W3",
        in_universe=n - len(miss), missing=",".join(miss), verdict=verdict, reason=reason,
    ))

out = f"{PB}/signature_registry_expansion_w3.csv"
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
exout = f"{PB}/signature_registry_expansion_w3_excluded.csv"
with open(exout, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f); w.writerow(["pmcid", "rule", "reason"]); w.writerows(EXCLUDED)

print(f"W3 纳入 {len(rows)}；本轮排除 {len(EXCLUDED)}")
for r in rows:
    print(f"  {r['sig_id']} {r['pmcid']} {r['n_genes']} genes 覆盖 {r['in_universe']}/{r['n_genes']}")

h = hashlib.sha256(open(out, "rb").read()).hexdigest()
ha = hashlib.sha256(open(f"{PB}/02_amendment_inflammation_wave3.md", "rb").read()).hexdigest()
with open(f"{PB}/work/LOCK.txt", "a", encoding="utf-8") as f:
    f.write("\n--- inflammation wave 3 (2026-10-20 22:05) ---\n")
    f.write(f"amendment  02_amendment_inflammation_wave3.md sha256 {ha}\n")
    f.write(f"registry   signature_registry_expansion_w3.csv sha256 {h}\n")
    f.write(f"n_included {len(rows)}  n_excluded {len(EXCLUDED)}\n")
print(f"\nsha256 = {h}")
