# -*- coding: utf-8 -*-
"""修订 02：炎症层定向检索 + 摘要级初筛。

与修订 01 的差别：
  * 8 组检索式全部围绕"炎症基因集/炎症反应/细胞因子/趋化因子/急性期/补体/炎症小体"，
    并**显式排除**细胞类型签名（macrophage/neutrophil/T/NK）与构念外词（ferroptosis/pyroptosis/
    autophagy/m6A/hypoxia/ECM/CAF），这两类分别对应 CIRC 与 I3-OUT；
  * 摘要级初筛除了"2–15 个基因符号"，额外要求出现**炎症基因池**的明确来源措辞（INF-STRICT 的摘要代理）。
"""
import csv, json, math, os, re, time, urllib.parse, urllib.request, html

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
OUT = f"{PB}/work/candidates_wave3.json"

LUNG = ('(TITLE:"non-small cell lung" OR TITLE:"lung adenocarcinoma" OR TITLE:"lung squamous" '
        'OR TITLE:NSCLC OR ABSTRACT:"lung adenocarcinoma" OR ABSTRACT:"non-small cell lung" OR ABSTRACT:NSCLC)')
SIG = ('(ABSTRACT:"gene signature" OR ABSTRACT:"prognostic signature" OR ABSTRACT:"risk score" '
       'OR ABSTRACT:"risk model" OR ABSTRACT:"prognostic model" OR TITLE:signature OR TITLE:model)')
# 细胞类型签名 -> CIRC；构念外 -> I3-OUT
EXCL = (' AND NOT (TITLE:macrophage OR TITLE:neutrophil OR TITLE:"T cell" OR TITLE:"NK cell" '
        'OR TITLE:ferroptosis OR TITLE:pyroptosis OR TITLE:cuproptosis OR TITLE:autophagy '
        'OR TITLE:hypoxia OR TITLE:fibroblast OR TITLE:lncRNA OR TITLE:miRNA OR TITLE:circRNA '
        'OR ABSTRACT:"single-cell" OR ABSTRACT:"single cell")')

QUERIES = {
    "W3Q1_irg": f'(ABSTRACT:"inflammation-related genes" OR ABSTRACT:"inflammation related genes" '
                 f'OR ABSTRACT:"inflammatory response-related genes" OR ABSTRACT:"inflammatory response related genes" '
                 f'OR ABSTRACT:IRGs) AND {LUNG} AND {SIG}{EXCL} AND SRC:MED AND PUB_YEAR:[2020 TO 2026]',
    "W3Q2_infresp": f'ABSTRACT:"inflammatory response" AND {LUNG} AND {SIG}{EXCL} AND SRC:MED AND PUB_YEAR:[2020 TO 2026]',
    "W3Q3_cytokine": f'(ABSTRACT:cytokine OR ABSTRACT:cytokines OR ABSTRACT:chemokine OR ABSTRACT:chemokines) '
                     f'AND {LUNG} AND {SIG}{EXCL} AND SRC:MED AND PUB_YEAR:[2020 TO 2026]',
    "W3Q4_acutephase": f'(ABSTRACT:"acute phase" OR ABSTRACT:complement OR ABSTRACT:inflammasome '
                       f'OR ABSTRACT:prostaglandin OR ABSTRACT:"NF-kB" OR ABSTRACT:"nuclear factor-kappa B" '
                       f'OR ABSTRACT:STAT3 OR ABSTRACT:"interleukin") AND {LUNG} AND {SIG}{EXCL} '
                       f'AND SRC:MED AND PUB_YEAR:[2020 TO 2026]',
    "W3Q5_sysinf": f'(ABSTRACT:"systemic inflammation" OR ABSTRACT:"systemic inflammatory" '
                   f'OR ABSTRACT:"inflammation-based" OR ABSTRACT:"inflammatory index") '
                   f'AND {LUNG} AND {SIG}{EXCL} AND SRC:MED AND PUB_YEAR:[2020 TO 2026]',
    "W3Q6_title": ('(TITLE:inflammation OR TITLE:inflammatory OR TITLE:"inflammation-related") '
                   'AND (TITLE:signature OR TITLE:model OR TITLE:score OR TITLE:index) '
                   'AND (TITLE:lung OR TITLE:pulmonary OR TITLE:NSCLC) '
                   f'AND SRC:MED AND PUB_YEAR:[2020 TO 2026]{EXCL}'),
    "W3Q7_lusc_inf": ('(ABSTRACT:"lung squamous" OR TITLE:"lung squamous" OR ABSTRACT:"squamous cell lung") '
                      'AND (ABSTRACT:inflammation OR ABSTRACT:inflammatory OR ABSTRACT:cytokine) '
                      f'AND {SIG}{EXCL} AND SRC:MED AND PUB_YEAR:[2020 TO 2026]'),
    "W3Q8_prognostic_inf": f'({LUNG}) AND (ABSTRACT:prognosis OR ABSTRACT:prognostic OR ABSTRACT:survival) '
                           f'AND (ABSTRACT:"inflammation" OR ABSTRACT:"inflammatory") AND ABSTRACT:TCGA '
                           f'{EXCL} AND SRC:MED AND PUB_YEAR:[2020 TO 2026]',
}


def fetch(q, page_size=100, cursor="*"):
    params = {"query": q, "format": "json", "pageSize": str(page_size),
              "resultType": "lite", "cursorMark": cursor}
    req = urllib.request.Request(BASE + "?" + urllib.parse.urlencode(params),
                                 headers={"User-Agent": "mailto:research@example.org"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode("utf-8"))


hits = {}
for k, q in QUERIES.items():
    cursor, got, rows, total = "*", 0, [], None
    while True:
        d = fetch(q, 100, cursor)
        total = d.get("hitCount")
        res = d.get("resultList", {}).get("result", [])
        if not res:
            break
        for it in res:
            rows.append({"pmid": it.get("pmid"), "pmcid": it.get("pmcid"), "doi": it.get("doi"),
                         "title": html.unescape(it.get("title") or ""),
                         "journal": it.get("journalTitle"), "year": it.get("pubYear")})
        got += len(res)
        nc = d.get("nextCursorMark")
        if not nc or nc == cursor or got >= (total or 0) or got >= 500:
            break
        cursor = nc
        time.sleep(0.35)
    hits[k] = {"query": q, "hitCount": total, "fetched": len(rows), "rows": rows}
    print(f"{k:22s} total={total:<6} fetched={len(rows)}", flush=True)

os.makedirs(f"{PB}/work", exist_ok=True)
json.dump(hits, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)

# ---------- 去重 ----------
known = set()
for f in ["signature_registry.csv", "signature_registry_expansion.csv"]:
    for r in csv.DictReader(open(f"{PB}/{f}", encoding="utf-8-sig")):
        known.add(r["pmcid"])
for f in ["signature_registry_excluded.csv", "signature_registry_expansion_excluded.csv"]:
    for r in csv.DictReader(open(f"{PB}/{f}", encoding="utf-8-sig")):
        known.add(r["pmcid"])
for f in ["work/extraction_raw.json", "work/extraction_wave2.json"]:
    for x in json.load(open(f"{PB}/{f}", encoding="utf-8")):
        known.add(x.get("pmcid"))

seen, new = {}, []
for k, blk in hits.items():
    for r in blk["rows"]:
        key = r["pmcid"] or r["pmid"]
        if not key:
            continue
        if key in seen:
            seen[key]["queries"].append(k); continue
        r["queries"] = [k]; seen[key] = r
for r in seen.values():
    if r["pmcid"] and str(r["year"]).isdigit() and int(r["year"]) >= 2020 and r["pmcid"] not in known:
        new.append(r)
print(f"\n去重 {len(seen)} 条；排除已知与 2020 前、无 PMCID 后：{len(new)} 条")

# ---------- 取摘要 ----------
def fetch_abs(pmids):
    out = {}
    for i in range(0, len(pmids), 25):
        chunk = pmids[i:i + 25]
        q = " OR ".join(f"EXT_ID:{p}" for p in chunk)
        params = {"query": q, "format": "json", "pageSize": "100",
                  "resultType": "core", "cursorMark": "*"}
        req = urllib.request.Request(BASE + "?" + urllib.parse.urlencode(params),
                                     headers={"User-Agent": "mailto:research@example.org"})
        try:
            with urllib.request.urlopen(req, timeout=120) as rr:
                d = json.loads(rr.read().decode("utf-8"))
        except Exception as e:
            print("  ! batch", i, e, flush=True); time.sleep(3); continue
        for it in d.get("resultList", {}).get("result", []):
            ab = it.get("abstractText") or ""
            out[it.get("pmid")] = html.unescape(re.sub(r"<[^>]+>", " ", ab))
        time.sleep(0.4)
    return out

abs_map = fetch_abs([r["pmid"] for r in new if r["pmid"]])
print(f"拿到摘要 {len(abs_map)} 篇")

GENELIKE = re.compile(r"^[A-Z][A-Z0-9]{1,8}(?:[-/][A-Z0-9]{1,6})?$")
STOP = set("TCGA NSCLC LUAD LUSC OS PFS ROC AUC HR CI GO KEGG LASSO RNA DNA MRNA COX DEG DEGS GSEA "
           "GDSC IC50 TME TIME EMT TMB MSI FDA NCCN AJCC TNM IHC PCR QPCR SD SE IQR BMI CRP ALB "
           "NLR PLR LMR PNI CONUT GPS SII SIRI ALI HALP NGS WGCNA CIBERSORT ESTIMATE GSVA SSGSEA "
           "FDR BH IRG IRGS DEIRG DEIRGS DEMG DEMGS".split())
INF_POOL = re.compile(r"inflammation[- ]related genes?|inflammatory response[- ]related genes?|"
                      r"\bIRGs?\b|\bDEIRGs?\b|cytokine[- ]related genes?|inflammat\w* gene set|"
                      r"inflammation[- ]associated genes?|inflammatory genes?", re.I)
CIRC_TXT = re.compile(r"macrophage|neutrophil|\bNK cell|\bT[- ]cell|immune cell|infiltrat|"
                      r"immunogenic cell death|immune[- ]related|immunotherapy response", re.I)
OUT_TXT = re.compile(r"ferroptosis|cuproptosis|disulfidptosis|pyroptosis|autophagy|\bm6A\b|"
                     r"hypoxia|fibroblast|\bCAF\b|lysosom|DNA repair|necropoptosis|necroptosis|"
                     r"PANoptosis|senescence|lncRNA|miRNA|circRNA|single[- ]cell", re.I)


def gene_list(ab):
    if not ab:
        return [], 0
    t = re.sub(r"[()\[\]]", " ", re.sub(r"\s+", " ", ab))
    toks = re.split(r"[,;]| and ", t)
    best, cur = [], []
    for tok in toks:
        c = tok.strip().strip(".").strip()
        if len(c.split()) == 1 and GENELIKE.match(c) and c.upper() not in STOP and not re.match(r"^\d", c):
            cur.append(c)
        else:
            if len(cur) > len(best): best = cur
            cur = []
    if len(cur) > len(best): best = cur
    return best, len(best)


rows = []
for r in new:
    ab = abs_map.get(r["pmid"], "")
    gl, n = gene_list(ab)
    title = r["title"]
    rows.append({
        "pmcid": r["pmcid"], "pmid": r["pmid"], "doi": r["doi"], "year": r["year"],
        "journal": r["journal"], "title": title,
        "abstract_gene_list": ",".join(gl), "n_in_abstract": n,
        "inf_pool": bool(INF_POOL.search(ab) or INF_POOL.search(title)),
        "inf_pool_term": (INF_POOL.search(ab) or INF_POOL.search(title)).group(0) if (INF_POOL.search(ab) or INF_POOL.search(title)) else "",
        "circ_hit": (CIRC_TXT.search(title) or CIRC_TXT.search(ab[:600])).group(0) if (CIRC_TXT.search(title) or CIRC_TXT.search(ab[:600])) else "",
        "out_hit": (OUT_TXT.search(title) or OUT_TXT.search(ab[:600])).group(0) if (OUT_TXT.search(title) or OUT_TXT.search(ab[:600])) else "",
        "queries": "|".join(r["queries"]),
    })

with open(f"{PB}/work/wave3_screened.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

tier = [r for r in rows if 2 <= r["n_in_abstract"] <= 15 and r["inf_pool"]]
print(f"\n【Tier A】摘要含 2–15 基因符号 且 有炎症基因池措辞：{len(tier)}")
for r in sorted(tier, key=lambda x: -int(x["year"])):
    flag = []
    if r["circ_hit"]: flag.append("CIRC:" + r["circ_hit"])
    if r["out_hit"]: flag.append("OUT:" + r["out_hit"])
    print(f"  {r['pmcid']:14s} {r['year']} n={int(r['n_in_abstract']):2d} {r['inf_pool_term'][:26]:26s} "
          f"{r['title'][:62]:62s} {' '.join(flag)}")
print(f"\n-> work/wave3_screened.csv（{len(rows)} 行）")
