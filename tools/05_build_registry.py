# -*- coding: utf-8 -*-
"""生成签名登记表（signature_registry.csv）并锁定哈希。

基因列表来自正文/摘要的**逐篇人工核验**（evidence 字段保留原文片段），
不来自自动抽取的默认结果——自动抽取只用于定位，最终列表以人工核验为准。
"""
import csv, hashlib, json, os, re

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
CAND = os.path.join(PB, "work", "candidates_dedup.csv")

# 人工核验结果：sig_id, pmcid, pmid, 分层, 基因数, 基因列表, 核验依据, 依据原文
VERIFIED = [
 ("SIG01","PMC8464410","34580602","inflammation",9,
  "BTG2,CCL20,CD69,IL7R,MMP14,NMUR1,PCDH7,RNF144B,TPBG","abstract",
  "Thirteen inflammation-related genes ... LASSO regression further screened nine candidates (BTG2, CCL20, CD69, IL7R, MMP14, NMUR1, PCDH7, RNF144B, and TPBG)."),
 ("SIG02","PMC8766344","35069695","inflammation",5,
  "EREG,GPC3,IL7R,LAMP3,NMUR1","abstract+formula",
  "a five-IRG risk signature consisting of EREG, GPC3, IL7R, LAMP3, and NMUR1. Risk score = EREG x 0.076861 - GPC3 x 0.042023 - IL7R x 0.150135 - LAMP3 x 0.009471 - NMUR1 x 0.072508 (separator glyphs ambiguous; unweighted used)"),
 ("SIG03","PMC8964604","35372503","inflammation",7,
  "MMP14,BTG2,LAMP3,CCL20,TLR2,IL7R,PCDH7","abstract",
  "We constructed a novel inflammatory response-related 7-gene signature (MMP14, BTG2, LAMP3, CCL20, TLR2, IL7R, and PCDH7)."),
 ("SIG04","PMC9203558","35710585","inflammation",11,
  "ADM,GPC3,IL7R,NMI,NMUR1,PSEN1,PTPRE,PVR,SEMA4D,SERPINE1,SPHK1","abstract",
  "This established IRG signature including 11 IRGs (ADM, GPC3, IL7R, NMI, NMURI, PSEN1, PTPRE, PVR, SEMA4D, SERPINE1, SPHK1)"),
 ("SIG05","PMC9427773","36042374","inflammation",8,
  "CCL20,CCR2,GNAI3,ITGA5,NMI,PCDH7,PSEN1,SLC11A2","formula",
  "Risk Score = 0.1241*CCL20 - 0.6192*CCR2 + 0.3660*GNAI3 + 0.1240*ITGA5 + 0.2050*NMI + 0.1330*PCDH7 + 0.430*PSEN1 - 0.3382*SLC11A2"),
 ("SIG06","PMC9972108","36056909","inflammation",2,
  "KLF6,SGMS2","abstract",
  "A two-IRG signature consisting of KLF6 and SGMS2 was identified (LUSC)"),
 ("SIG07","PMC12259757","40053252","inflammation",2,
  "MS4A1,MS4A2","abstract",
  "a new model was developed based on two mRNAs (MS4A1 and MS4A2) in TCGA dataset"),
 ("SIG08","PMC13589103","42760759","inflammation",7,
  "BAK1,BMF,CYCS,FADD,IL1A,TFDP1,YWHAG","abstract",
  "The seven-gene prognostic model (comprising BAK1, BMF, CYCS, FADD, IL1A, TFDP1, and YWHAG)"),
 ("SIG09","PMC9125243","35615380","metabolism",10,
  "FBP2,ADH6,DHDH,PRKCB,INPP5J,ABAT,HK2,GNPNAT1,PLCB3,ACAT2","abstract",
  "ten genes ... including FBP2, ADH6, DHDH, PRKCB, INPP5J, ABAT, HK2, GNPNAT1, PLCB3, and ACAT2"),
 ("SIG10","PMC11521851","39478862","metabolism",10,
  "GPR115,SOAT2,CDH17,MOGAT2,COL11A1,TCN1,LGR5,SLC34A2,RHOV,DKK1","abstract",
  "a prognostic signature comprising 10 FAM-related genes (GPR115, SOAT2, CDH17, MOGAT2, COL11A1, TCN1, LGR5, SLC34A2, RHOV, and DKK1)"),
 ("SIG11","PMC8867215","35222371","metabolism",12,
  "LDHA,NSDHL,TP53INP2,FLT1,IRS1,ELOVL7,AGPS,FHL2,MED6,PLIN3,VDAC1,SULT2B1","formula",
  "Risk score = (0.269752653768969*LDHA)+(0.114984371234971*NSDHL)+(0.0310569577143764*TP53INP2)+(0.17593670859331*FLT1)+(0.103773854013971*IRS1)+(0.133274629300893*ELOVL7)+(0.0763835280754936*AGPS)+(0.00688026723004499*FHL2)+(0.549703365193899*MED6)+(0.109916205934083*PLIN3)+(0.0437103013689035*VDAC1)+(0.160609136698149*SULT2B1)"),
 ("SIG12","PMC13511999","42650107","metabolism",6,
  "CYP4B1,ACOXL,DPEP2,HPGDS,CA4,ALOX15","abstract",
  "The six-gene signature comprised CYP4B1, ACOXL, DPEP2, HPGDS, CA4, and ALOX15."),
 ("SIG13","PMC13128028","42003422","metabolism",7,
  "WIF1,GATA2,CD69,ID1,C4BPA,WFDC2,CCL20","abstract",
  "Seven histidine metabolism-related epithelial cell-specific genes ... (WIF1, GATA2, CD69, ID1, C4BPA, WFDC2, and CCL20)"),
]

EXCLUDED = [
 ("PMC8647082","34881236","8-gene signature; gene names only in Figure/Supplementary, not retrievable from text","E1"),
 ("PMC9736863","36497225","8-gene signature; construct shown in Figure S1 only","E1"),
 ("PMC9455632","36091041","states 10 genes but text contains inconsistent 4- and 7-gene lists; final list unverifiable","E2"),
 ("PMC10124980","37083272","no explicit final gene list in text (model described by score only)","E1"),
 ("PMC12102027","40407957","INF scoring system from subtypes; no fixed final gene list (MMP14 validated experimentally)","E1"),
 ("PMC10493788","37701116","final gene list not stated in text/abstract","E1"),
]

# 元数据（题名/期刊/DOI）从检索结果取，避免凭记忆填写
meta = {}
with open(CAND, encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r.get("pmcid"):
            meta[r["pmcid"]] = r

rows = []
for sid, pmc, pmid, stratum, n, genes, how, ev in VERIFIED:
    m = meta.get(pmc, {})
    rows.append(dict(
        sig_id=sid, pmcid=pmc, pmid=pmid, doi=m.get("doi", ""), journal=m.get("journal", ""),
        year=m.get("year", ""), stratum=stratum, n_genes=n, genes=genes,
        n_listed=len([g for g in genes.split(",") if g.strip()]),
        verification=how, evidence=ev, title=m.get("title", ""),
    ))

out = os.path.join(PB, "signature_registry.csv")
with open(out, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

ex_out = os.path.join(PB, "signature_registry_excluded.csv")
with open(ex_out, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["pmcid", "pmid", "exclusion_reason", "rule"])
    w.writerows(EXCLUDED)

# 一致性断言
bad = [r["sig_id"] for r in rows if r["n_genes"] != r["n_listed"]]
assert not bad, "基因数与列表长度不一致: %s" % bad
assert all(2 <= r["n_genes"] <= 15 for r in rows), "基因数越界"

h = hashlib.sha256(open(out, "rb").read()).hexdigest()
hc = hashlib.sha256(open(os.path.join(PB, "00_criteria_preregistered.md"), "rb").read()).hexdigest()
with open(os.path.join(PB, "work", "LOCK.txt"), "w", encoding="utf-8") as f:
    f.write("locked_at=2026-10-20T20:35+08:00\n")
    f.write("signature_registry.csv sha256=%s\n" % h)
    f.write("00_criteria_preregistered.md sha256=%s\n" % hc)
    f.write("n_included=%d n_excluded=%d\n" % (len(rows), len(EXCLUDED)))

print("纳入 %d 个签名；排除 %d 个" % (len(rows), len(EXCLUDED)))
print("registry sha256 =", h)
print("criteria sha256 =", hc)
for r in rows:
    print("  %s %-13s %s | %s | %d genes" % (r["sig_id"], r["pmcid"], r["journal"][:28], r["stratum"], r["n_genes"]))
