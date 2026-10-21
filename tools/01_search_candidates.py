# -*- coding: utf-8 -*-
"""W1：候选签名检索（Europe PMC），只输出可核验的元数据，不做任何结论。
检索目标：NSCLC 中"以营养/炎症构念为锚点"的预后基因签名（2020-01-01 起，有 PMC 全文）。
"""
import json, re, time, urllib.parse, urllib.request, html, csv

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
OUT = r"D:/WorkBuddy_Y02/out_ali_directions/pathB/work/candidates_raw.json"

QUERIES = {
    "Q1_title_nutinf_sig": '(TITLE:"non-small cell lung" OR TITLE:"lung adenocarcinoma" OR TITLE:"lung squamous" OR TITLE:NSCLC) '
                           'AND (TITLE:nutrition* OR TITLE:"nutritional" OR TITLE:inflammat* OR TITLE:"immune-inflammation") '
                           'AND (TITLE:signature OR TITLE:model OR TITLE:score) AND (TITLE:prognos* OR TITLE:survival)',
    "Q2_abs_nutinf_sig": 'ABSTRACT:("nutrition" OR "nutritional" OR "nutritional index" OR "prognostic nutritional index") '
                         'AND ABSTRACT:("inflammation" OR "inflammatory" OR "systemic immune-inflammation") '
                         'AND ABSTRACT:("prognostic signature" OR "gene signature" OR "risk score" OR "prognostic model") '
                         'AND ABSTRACT:("non-small cell lung" OR "lung adenocarcinoma" OR "lung squamous" OR NSCLC)',
    "Q3_abs_index_sig": 'ABSTRACT:("advanced lung cancer inflammation index" OR "prognostic nutritional index" OR "controlling nutritional status" '
                        'OR "geriatric nutritional risk index" OR "systemic immune-inflammation index" OR "neutrophil-to-lymphocyte ratio") '
                        'AND ABSTRACT:("lung adenocarcinoma" OR "non-small cell lung" OR NSCLC) '
                        'AND ABSTRACT:(transcriptom* OR "gene expression" OR "gene signature" OR TCGA OR bioinformatic*)',
    "Q4_abs_nutrition_metab_sig": 'ABSTRACT:("amino acid" OR "lipid metabolism" OR "fatty acid" OR "glucose metabolism") '
                                  'AND ABSTRACT:("prognostic signature" OR "risk model") '
                                  'AND ABSTRACT:("lung adenocarcinoma" OR "non-small cell lung")',
}

def fetch(q, page_size=100, cursor="*"):
    params = {"query": q, "format": "json", "pageSize": str(page_size), "resultType": "lite", "cursorMark": cursor}
    url = BASE + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "mailto:research@example.org"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode("utf-8"))

allhits = {}
for k, q in QUERIES.items():
    cursor = "*"
    got = 0
    total = None
    rows = []
    while True:
        d = fetch(q, 100, cursor)
        total = d.get("hitCount")
        res = d.get("resultList", {}).get("result", [])
        if not res:
            break
        for it in res:
            rows.append({
                "pmid": it.get("pmid"), "pmcid": it.get("pmcid"), "doi": it.get("doi"),
                "title": html.unescape(it.get("title") or ""), "journal": it.get("journalTitle"),
                "year": it.get("pubYear"), "oa": it.get("isOpenAccess"),
                "has_ft": it.get("hasTextMinedTerms"), "in_epmc": it.get("inEPMC"),
            })
        got += len(res)
        nc = d.get("nextCursorMark")
        if not nc or nc == cursor or got >= (total or 0) or got >= 400:
            break
        cursor = nc
        time.sleep(0.4)
    allhits[k] = {"query": q, "hitCount": total, "fetched": len(rows), "rows": rows}
    print("%-30s total=%-5s fetched=%d" % (k, total, len(rows)), flush=True)

import os
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(allhits, f, ensure_ascii=False, indent=1)
print("->", OUT)
