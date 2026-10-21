# -*- coding: utf-8 -*-
"""修订 03：免疫分层（W4）登记表。逐条人工核验后写死，含排除理由。"""
import csv, hashlib, json, os

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
meta = {r["pmcid"]: r for r in
        csv.DictReader(open(f"{PB}/work/wave3_screened.csv", encoding="utf-8-sig"))}
universe = set(l.strip() for l in
               open(os.path.join(PB, "work", "genes.tsv"), encoding="utf-8").read().splitlines()[1:])
ext = {x["pmcid"]: x for x in json.load(open(f"{PB}/work/extraction_wave4.json", encoding="utf-8"))}

DECISIONS = [
    ("PMC12813068", "immune", 8,
     "PSMD1,ANGPTL4,LTBP3,MIF,NFATC3,NR1D2,PLXNB3,SP1",
     "IRGs were obtained from the ImmPort database (https://www.immport.org/).",
     "INCL", "免疫：IRG 池来自 ImmPort，8/8 与自述一致，LUSC"),
    ("PMC8182055", "immune", 8,
     "IKBKB,LTBR,MIF,PPARD,PPIA,PSME3,S100A6,SEMA4B",
     "A total of 2498 IRGs data were obtained from the ImmPort database (Bhattacharya et al., 2018).",
     "INCL", "免疫：IRG 池来自 ImmPort（2,498 基因），8/8 一致，LUAD 淋巴结转移"),
    ("PMC8643179", "immune", 6,
     "SREBF2,GP2,BMX,NR1H4,DDX41,GOPC",
     "IRGs from the InnateDB database.",
     "INCL", "免疫（先天免疫）：IRG 池来自 InnateDB，6/6 一致，早期 LUSC"),
    ("PMC7485220", "immune", 7,
     "RELA,NFKB1,SP1,VDR,FOS,PPARG,STAT3",
     "To investigate the potential ability of TFs in regulating these clinically relevant IRGs, a total of 318 TFs were downloaded from Cistrome.",
     "INCL", "免疫：免疫相关转录因子签名（TF 池来自 Cistrome，调控 IRG）；池为 TF 而非 IRG 本身，已标注 anchor=edge"),
    ("PMC9136055", "immune", 4,
     "DDX56,CTSL,ZC3H12D,PSMC5",
     "1,793 unique IRGs were extracted from the Immunology Database and Analysis Portal (ImmPort) database.",
     "INCL", "免疫：免疫相关 RNA 结合蛋白（RBP）签名，IRG 池来自 ImmPort，4/4 一致"),
    ("PMC7781839", "immune", 5,
     "CD40LG,TNFRSF6B,TNFSF13,TNFRSF13C,TNFRSF19",
     "Identification of a costimulatory molecule-based signature ... (共刺激分子基因集；5 个基因全部为 TNFRSF/TNFSF 家族共刺激分子)",
     "INCL", "免疫：共刺激分子签名，5/5 一致；列表自身即共刺激分子家族，构念自明"),
]

EXCLUDED = [
    ("PMC13313863", "I4", "自述 1 基因与抽到 4 个不符；且池来源句只是 GO 富集注释"),
    ("PMC7368195", "I4-POOL", "列表可读（8 基因）但基因池来源不可核验（仅 GO 富集句）；与 PMC11094586 有 4/8 基因重叠"),
    ("PMC7493716", "I4", "自述 4 基因与抽到 9 个不符"),
    ("PMC7372779", "I4", "自述 6 基因与抽到 11 个不符（BTLA/CD27/TNFRSF 家族，列表本身像免疫共刺激分子，但最终构成不可定）"),
    ("PMC11094586", "I4-POOL", "无自述基因数，池来源仅 GO/KEGG 富集句"),
    ("PMC8554550", "I4", "抽取含非基因符号 HR（假阳性）；自述 3 基因与抽到 4 个不符"),
    ("PMC8222816", "I2", "研究对象为肝细胞癌（HCC），非 NSCLC"),
    ("PMC7795632", "I3-OUT", "细胞死亡（自噬/凋亡/坏死）基因池，构念外；且自述 21 基因超出 I6"),
    ("PMC11745784", "I4", "免疫原性细胞死亡（ICD）签名，自述 17 基因与抽到 14 个不符"),
    ("PMC8317604", "I4", "题名为脂质代谢、抽取出的却是 CXCL9/CXCL10/IFNG 等干扰素响应基因，签名构成与基因池均不可核验"),
]

rows = []
for i, (pmcid, stratum, n, genes, ev, verdict, reason) in enumerate(DECISIONS, start=44):
    gl = [g.strip() for g in genes.split(",")]
    x = ext.get(pmcid, {})
    m = meta.get(pmcid, {})
    assert len(gl) == n, f"{pmcid}: {len(gl)} != {n}"
    assert x.get("n_extracted") == n, f"{pmcid}: 抽取 {x.get('n_extracted')} != {n}"
    assert all(g in universe for g in gl), f"{pmcid}: 基因不在全集 {missing}"
    rows.append(dict(
        sig_id=f"SIG{i:02d}", pmcid=pmcid, pmid=m.get("pmid", ""), doi=m.get("doi", ""),
        journal=m.get("journal", ""), year=m.get("year", ""), stratum=stratum,
        n_genes=n, genes=",".join(gl), n_listed=x.get("n_stated") or n,
        verification="body_provenance", evidence=ev, title=m.get("title", ""), wave="W4",
        in_universe=n, missing="", verdict=verdict, reason=reason,
    ))

out = f"{PB}/signature_registry_expansion_w4.csv"
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
exout = f"{PB}/signature_registry_expansion_w4_excluded.csv"
with open(exout, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f); w.writerow(["pmcid", "rule", "reason"]); w.writerows(EXCLUDED)

print(f"W4 纳入 {len(rows)}；排除 {len(EXCLUDED)}")
for r in rows:
    print(f"  {r['sig_id']} {r['pmcid']} {r['stratum']:6s} {r['n_genes']} g  {r['title'][:58]}")

h = hashlib.sha256(open(out, "rb").read()).hexdigest()
ha = hashlib.sha256(open(f"{PB}/03_amendment_immune_stratum.md", "rb").read()).hexdigest()
with open(f"{PB}/work/LOCK.txt", "a", encoding="utf-8") as f:
    f.write("\n--- immune stratum wave 4 (2026-10-20 22:35) ---\n")
    f.write(f"amendment  03_amendment_immune_stratum.md sha256 {ha}\n")
    f.write(f"registry   signature_registry_expansion_w4.csv sha256 {h}\n")
    f.write(f"n_included {len(rows)}  n_excluded {len(EXCLUDED)}\n")
print(f"\nsha256 = {h}")
