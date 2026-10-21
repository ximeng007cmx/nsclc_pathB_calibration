# -*- coding: utf-8 -*-
"""构建扩容签名登记表（wave2）。

纪律：
  * 每条的决定与理由写在 DECISIONS 里，可逐条复核；
  * assert 基因数 == 列表长度，任何不一致立即中止（不静默修补）；
  * 不在本脚本里"为了结果好看"增删任何一条，只执行上文已判定的决定。
"""
import csv, hashlib, json, os

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
meta = {r["pmcid"]: r for r in
        csv.DictReader(open(f"{PB}/work/candidates_wave2_screened.csv", encoding="utf-8-sig"))}
universe = set(l.strip() for l in
               open(os.path.join(PB, "work", "genes.tsv"), encoding="utf-8").read().splitlines()[1:])

# (pmcid, stratum, verification, genes, evidence, verdict, reason)
DECISIONS = [
    # ---------------- A_OK：自述数与列表一致 ----------------
    ("PMC9047773", "metabolism", "abstract",
     "AGER,AHSG,ALDH2,CIDEC,CYP17A1,FBP1,GNB3,GZMB,IGFBP1,SORD,SOX2,TRH,TYMS",
     "prognostic signature comprising 13 EM-related genes (AGER, AHSG, ALDH2, CIDEC, CYP17A1, FBP1, GNB3, GZMB, IGFBP1, SORD, SOX2, TRH and TYMS)",
     "INCL", "能量代谢 13 基因，TCGA-LUAD，摘要给出完整列表"),
    ("PMC8806683", "metabolism", "abstract",
     "HAVCR1,SPN,GAPDH,ANGPTL4,PRSS3,KRT8,LDHA,HMMR,SLC2A1,CYP24A1,LOXL2,TIMP1",
     "We constructed a prediction signature with 12 genes (HAVCR1, SPN, GAPDH, ANGPTL4, PRSS3, KRT8, LDHA, HMMR, SLC2A1, CYP24A1, LOXL2, TIMP1)",
     "INCL", "铁代谢 12 基因，LASSO-Cox，列表完整"),
    ("PMC13190771", "metabolism", "abstract",
     "HILPDA,ACHE,PPARD,ETNK1,ST3GAL2,NCOA2,HACD1,MBTPS2,UGT8,LPIN2,SPHK2,AKR1C4",
     "A prognostic signature comprising 12 lipid metabolism-related genes (HILPDA, ACHE, PPARD, ETNK1, ST3GAL2, NCOA2, HACD1, MBTPS2, UGT8, LPIN2, SPHK2, and AKR1C4)",
     "INCL", "脂质代谢 12 基因"),
    ("PMC12688223", "metabolism", "abstract",
     "LTBP1,MFGE8,ACTN1,CD59,CDC25C,SAAL1,SFXN4,PTTG1",
     "An 8-gene OXPHOS-related signature (LTBP1, MFGE8, ACTN1, CD59, CDC25C, SAAL1, SFXN4, PTTG1) was identified",
     "INCL", "氧化磷酸化 8 基因，LUSC"),
    ("PMC8034925", "metabolism", "title+body",
     "PLAUR,F2,UGT2B17,GNG7,IDO2,ST3GAL6,PIK3CG,GLS2",
     "A signature model of eight genes (PLAUR, F2, UGT2B17, GNG7, IDO2, ST3GAL6, PIK3CG, and GLS2)",
     "INCL", "代谢 8 基因签名"),
    ("PMC13559650", "metabolism", "abstract",
     "GALT,GANC,PGM1,GALM,B4GALT1,PGM2",
     "A six-gene prognostic model (GALT, GANC, PGM1, GALM, B4GALT1, PGM2) was developed",
     "INCL", "半乳糖代谢 6 基因"),
    ("PMC7557101", "metabolism", "abstract",
     "CYP2C9,DGAT1,UGT1A6,INS,HPGDS,LPL",
     "the six survival-related hub genes (CYP2C9, DGAT1, UGT1A6, INS, HPGDS, and LPL)",
     "INCL", "脂质代谢 6 hub 基因 + nomogram"),
    ("PMC7718790", "metabolism", "abstract",
     "PFKP,PKM,TPI1,LDHA,PTGES,TYMS",
     "The six-gene molecular signature (PFKP, PKM, TPI1, LDHA, PTGES, and TYMS)",
     "INCL", "代谢 6 基因"),
    ("PMC13073725", "metabolism", "abstract",
     "YWHAZ,HSPD1,NOTCH3,PGK1,PPARG",
     "we established a five-gene prognostic signature (YWHAZ, HSPD1, NOTCH3, PGK1, and PPARG)",
     "INCL", "线粒体代谢 5 基因"),
    ("PMC13261843", "metabolism", "abstract",
     "ACSL3,MCAT,NDUFAB1,OLAH,ACSL4",
     "A five-gene FAM-related risk signature (ACSL3, MCAT, NDUFAB1, OLAH, ACSL4) was identified",
     "INCL", "脂肪酸代谢 5 基因"),
    ("PMC11963273", "metabolism", "abstract",
     "S100P,GPX2,PRC1,ARNTL2,RGS20",
     "an optimal gene signature comprising five metabolism-stemness-related DEGs (msrDEGs) (S100P, GPX2, PRC1, ARNTL2, and RGS20)",
     "INCL", "代谢×干性交集 5 基因；构念为部分代谢，已标注"),
    ("PMC9174756", "metabolism", "abstract",
     "FKBP4,HMMR,B4GALT1,SLC2A1,STC1",
     "We constructed a 5-gene signature (FKBP4, HMMR, B4GALT1, SLC2A1, STC1)",
     "INCL", "糖酵解 5 基因"),
    ("PMC12955334", "metabolism", "abstract",
     "ACSL4,ACSS3,CHPT1,ACSM5",
     "the four model genes (ACSL4, ACSS3, CHPT1, and ACSM5)",
     "INCL", "线粒体脂质代谢 4 基因"),
    ("PMC12859070", "metabolism", "abstract",
     "ERCC1,KYNU,AKR1A1,PPIA",
     "multiomics prognostic model",
     "INCL", "脂肪酸/谷氨酰胺 4 基因模型（PPIA 为主，另有 3 基因）"),

    # ---------------- B_NOSTATED：无自述基因数，但证据句明确给出列表 ----------------
    ("PMC10619909", "metabolism", "body_enumerated",
     "SLC2A1,PKP2,FAM83A,TCN1,MS4A1,CLIC6,UBE2S,RRM2,CDC45,IGF2BP1,ANGPTL4,CD109",
     "A total of 12 key genes (SLC2A1, PKP2, FAM83A, TCN1, MS4A1, CLIC6, UBE2S, RRM2, CDC45, IGF2BP1, ANGPTL4, and CD109) were screened and used to develop a Risk",
     "INCL", "脂肪酸代谢 12 基因，机器学习筛选"),
    ("PMC12101037", "metabolism", "body_enumerated",
     "ANLN,PLK1,DLGAP5,PRC1,CYP4B1,MKI67,KIF23,BIRC5,TK1",
     "nine prognostic genes (ANLN, PLK1, DLGAP5, PRC1, CYP4B1, MKI67, KIF23, BIRC5, TK1)",
     "INCL", "叶酸代谢 9 基因"),
    ("PMC12098254", "metabolism", "evidence_sentence",
     "SCGB3A2,PGC,ADH7,RPS4Y1,KRT6A",
     "5 key FAMGs (SCGB3 A2, PGC, ADH7, RPS4Y1, and KRT6 A) were identified",
     "INCL", "脂肪酸代谢 5 基因；原文排版把 SCGB3A2/KRT6A 拆成两段，已合并"),
    ("PMC11404934", "metabolism", "body_enumerated",
     "ALG3,COL7A1,KL,MST1,SLC52A1",
     "5 prognostic related genes (ALG3, COL7A1, KL, MST1, and SLC52A1) were used to establish a LUAD polygenic",
     "INCL", "代谢相关 5 基因"),
    ("PMC8419226", "inflammation", "body_enumerated",
     "IL7R,IL5RA,IL20RB,IL11,IL22RA1",
     "A signature based on five identified IL(R)s (IL7R, IL5RA, IL20RB, IL11, IL22RA1) was constructed using the TCGA dataset",
     "INCL", "IL 受体 5 基因，炎症构念（本波唯一炎症类）"),
    ("PMC11084143", "metabolism", "body_enumerated",
     "MTHFD1,PDGFB,PIK3R3,CHKB",
     "4 choline metabolism-related genes (MTHFD1, PDGFB, PIK3R3, CHKB) were screened and developed a risk sign",
     "INCL", "胆碱代谢 4 基因"),
    ("PMC8361602", "metabolism", "body_enumerated",
     "TYMS,NME4,LDHA,SMOX",
     "The MGRS model comprising 4 genes (TYMS, NME4, LDHA, and SMOX) was developed",
     "INCL", "代谢-免疫整合模型 4 基因（以代谢为主，已标注）"),

    # ---------------- C_MISMATCH 中被判定可救者 ----------------
    ("PMC11726977", "metabolism", "evidence_sentence",
     "CPS1,ABCC2,CDKN3,SLC2A1,IGFBP1",
     "5 BMR-DEGs (CPS1, ABCC2, CDKN3, SLC2A1, and IGFBP1) were identified as individual prognosti",
     "INCL", "丁酸代谢 5 基因；自述数被误读为 1"),
    ("PMC12474791", "metabolism", "evidence_sentence",
     "LDHA,KYNU,SLC2A1,CFTR,MAOB",
     "A risk model was constructed using five characteristic genes (LDHA, KYNU, SLC2A1, CFTR, MAOB)",
     "INCL", "丙酸代谢 5 特征基因"),
    ("PMC9262104", "metabolism", "evidence_sentence",
     "DTYMK,GCDH,HEMK1,NEU1,NNT,NT5C3A,POLR3G,PPOX,SMS,TKFC,TRDMT1,UAP1,WARS2",
     "Metabolism Score was calculated based on 14 hub genes (DTYMK, GCDH, HEMK1, JMJD7.PLA2G4B, NEU1, NNT, NT5C3A, POLR3G, PPOX, SMS, TKFC, TRDMT1, UAP1 and WARS2)",
     "INCL", "代谢 13 基因；第 14 个 JMJD7.PLA2G4B 为读通融合符号，非单一基因，已剔除并记录"),
    ("PMC7778529", "metabolism", "evidence_sentence",
     "PKP2,LDHA,HMMR,COL5A1,B3GNT3",
     "Five genes (PKP2, LDHA, HMMR, COL5A1 and B3GNT3) were eventually identified to calculate",
     "INCL", "糖酵解 5 基因，NSCLC"),
    ("PMC7335383", "metabolism", "evidence_sentence",
     "GNPNAT1,LPGAT1,TYMS,LDHA,PTGES",
     "a novel 5-gene metabolic prognostic signature (including GNPNAT1, LPGAT1, TYMS, LDHA, and PTGES)",
     "INCL", "代谢 5 基因"),
    ("PMC13397370", "metabolism", "evidence_sentence",
     "ANLN,CIP2A,MEST,WDR76",
     "A gene signature of 4 genes (ANLN, CIP2A, MEST, WDR76) divided the patients into high-risk and",
     "INCL", "谷氨酰胺代谢 4 基因"),
    ("PMC7646164", "metabolism", "evidence_sentence",
     "TYMS,ALDH2,PKM,GNPNAT1,LDHA,ENTPD2,NT5E,MAOB",
     "the prognostic indexes based on 8 metabolism-related genes (TYMS, ALDH2, PKM, GNPNAT1, LDHA, ENTPD2, NT5E, MAOB) were established",
     "INCL", "代谢 8 基因预后指数（题名虽为 TME，候选池为代谢基因）"),
]

EXCLUDED = [
    # (pmcid, rule, reason)
    ("PMC7066120", "I3", "RNA 结合蛋白（RBP）构念，非炎症非营养/代谢"),
    ("PMC12816179", "I3", "锚点为 GPRC5A+ 细胞状态模块，非炎症/代谢"),
    ("PMC11283104", "CIRC", "NK 细胞基因签名，属免疫含量代理，循环论证"),
    ("PMC11876514", "CIRC", "巨噬细胞（TAM）基因签名，属免疫含量代理"),
    ("PMC10788202", "I2", "泛癌签名，训练集非 NSCLC 专属"),
    ("PMC9926575", "E4", "9 个基因为驱动代谢亚型的 hub 基因，未构成预后签名"),
    ("PMC9459509", "I3", "焦亡（pyroptosis）相关基因集，构念外（规则表 INF 命中 NLRP3 而漏判，人工纠正）"),
    ("PMC8028741", "CIRC", "TME 相关签名，等价于免疫/基质含量"),
    ("PMC13513166", "I3", "自噬分泌签名，构念外"),
    ("PMC11987467", "I3", "线粒体质量控制签名，构念外"),
    ("PMC8034940", "I3", "血小板相关基因，构念外"),
    ("PMC10864718", "E4", "meta 分析得到的顺铂耐药相关基因，非预后签名"),
    ("PMC8330534", "E4", "抽到的是 TP53/KRAS/EGFR 等驱动突变基因，抽取假阳性"),
    ("PMC9315097", "I4", "仅能从公式抽出 14 个符号，与自述数（2）不符，无法核验"),
    ("PMC7520091", "I6", "23 个基因，超出 2–15 上限"),
    ("PMC12290164", "I3", "肺癌与肺部感染共享血签名，构念不清"),
    ("PMC13414110", "I2", "6 个基因来自文献调研而非 TCGA 建模"),
    ("PMC10281056", "I3", "柚皮素网络药理学研究，构念外"),
    ("PMC7720456", "I4", "证据句引述的是他人 6 基因签名，本研究的 4 基因列表未能核验"),
    ("PMC8209576", "E2", "题名 6 基因、正文 5 基因，前后不一致"),
    ("PMC11341317", "I3", "免疫×氨基酸代谢混合签名，且 17 个中仅 4 个一致，无法定列表"),
    ("PMC12851720", "I4", "自述 7 基因、正文仅列 5 个验证基因，最终列表不可定"),
    ("PMC9784226", "I4", "仅给出 7 个'负相关'基因，未给出签名最终构成"),
    ("PMC8531531", "I3", "FDX1 单基因机制研究，铜死亡构念"),
    ("PMC13633477", "I4", "自述 3 基因、抽到 6 个，列表不可定"),
    ("PMC12595263", "E4", "6 个糖酵解基因为预后相关基因，非风险模型构成"),
]

rows = []
for i, (pmcid, stratum, verif, genes, ev, verdict, reason) in enumerate(DECISIONS, start=14):
    gl = [g.strip() for g in genes.split(",") if g.strip()]
    m = meta.get(pmcid, {})
    n_stated = None
    for x in json.load(open(f"{PB}/work/extraction_wave2.json", encoding="utf-8")):
        if x["pmcid"] == pmcid:
            n_stated = x.get("n_stated")
    assert len(gl) == len(set(gl)), f"{pmcid} 有重复基因"
    assert 2 <= len(gl) <= 15, f"{pmcid} 基因数 {len(gl)} 越界"
    miss = [g for g in gl if g not in universe]
    rows.append(dict(
        sig_id=f"SIG{i:02d}", pmcid=pmcid, pmid=m.get("pmid", ""), doi=m.get("doi", ""),
        journal=m.get("journal", ""), year=m.get("year", ""), stratum=stratum,
        n_genes=len(gl), genes=",".join(gl), n_listed=n_stated if n_stated else len(gl),
        verification=verif, evidence=ev, title=m.get("title", ""),
        wave="W2", in_universe=len(gl) - len(miss), missing=",".join(miss), verdict=verdict,
        reason=reason,
    ))

out = f"{PB}/signature_registry_expansion.csv"
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

exout = f"{PB}/signature_registry_expansion_excluded.csv"
with open(exout, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f); w.writerow(["pmcid", "rule", "reason"]); w.writerows(EXCLUDED)

print(f"纳入 {len(rows)} 条；排除 {len(EXCLUDED)} 条")
cov = [(r["sig_id"], r["pmcid"], r["n_genes"], r["in_universe"]) for r in rows if r["in_universe"] < r["n_genes"]]
print("基因全集覆盖不足者：", cov if cov else "无")
print("\n分层：", {s: sum(1 for r in rows if r["stratum"] == s) for s in {r["stratum"] for r in rows}})
print("验证方式：", {v: sum(1 for r in rows if r["verification"] == v) for v in {r["verification"] for r in rows}})

h = hashlib.sha256(open(out, "rb").read()).hexdigest()
he = hashlib.sha256(open(exout, "rb").read()).hexdigest()
ha = hashlib.sha256(open(f"{PB}/01_amendment_expansion_to20.md", "rb").read()).hexdigest()
with open(f"{PB}/work/LOCK.txt", "a", encoding="utf-8") as f:
    f.write("\n--- expansion wave (2026-10-20 21:45) ---\n")
    f.write(f"amendment  01_amendment_expansion_to20.md  sha256 {ha}\n")
    f.write(f"registry   signature_registry_expansion.csv sha256 {h}\n")
    f.write(f"excluded   signature_registry_expansion_excluded.csv sha256 {he}\n")
    f.write(f"n_included {len(rows)}  n_excluded {len(EXCLUDED)}\n")
print(f"\nsha256 registry = {h}\nsha256 excluded = {he}")
print(f"-> {out}\n-> {exout}")
