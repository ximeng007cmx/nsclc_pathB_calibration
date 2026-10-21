# -*- coding: utf-8 -*-
"""扩容检索（第二波，目标把签名从 13 补到 20）。

纪律说明
--------
本脚本只负责"找候选"，不参与任何纳入判断。纳入仍按 00_criteria_preregistered.md 的
I1-I6 / E1-E5 执行，扩容波次在 01_amendment_expansion_to20.md 里声明为 post-hoc。

与第一波（01_search_candidates.py）的差别：
  1. 扩到 I3 构念的更多子锚点（氨基酸/胆固醇/糖酵解/谷氨酰胺/维生素/微量元素），
     但仍限定"营养或代谢"与"炎症"两大类；
  2. **不检索 immune-related / 免疫相关签名**——用免疫签名去校正免疫含量是循环论证，
     这类签名即便命中也不纳入，理由记入修订文件；
  3. 增加"摘要里直接出现基因列表"的加权检索（"-gene" + comprising/consisting），
     因为 I4 要求列表可在正文层面核验，摘要级列表最容易核。
"""
import json, os, re, time, urllib.parse, urllib.request, html

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
OUT = r"D:/WorkBuddy_Y02/out_ali_directions/pathB/work/candidates_wave2.json"

LUNG = ('(TITLE:"non-small cell lung" OR TITLE:"lung adenocarcinoma" OR TITLE:"lung squamous" '
        'OR TITLE:NSCLC OR ABSTRACT:"lung adenocarcinoma" OR ABSTRACT:"non-small cell lung" OR ABSTRACT:NSCLC)')

ANCHOR_NUT = ('(ABSTRACT:"lipid metabolism" OR ABSTRACT:"fatty acid metabolism" OR ABSTRACT:"amino acid metabolism" '
              'OR ABSTRACT:"glucose metabolism" OR ABSTRACT:"glycolysis" OR ABSTRACT:"cholesterol metabolism" '
              'OR ABSTRACT:"glutamine metabolism" OR ABSTRACT:"nutrition" OR ABSTRACT:"nutritional" '
              'OR ABSTRACT:"dietary" OR ABSTRACT:"vitamin" OR ABSTRACT:"trace element" OR ABSTRACT:"metabolic reprogramming")')
ANCHOR_INF = ('(ABSTRACT:"inflammation-related" OR ABSTRACT:"inflammatory response" OR ABSTRACT:"inflammation" '
              'OR ABSTRACT:"systemic inflammation" OR ABSTRACT:"cytokine")')
SIG = ('(ABSTRACT:"prognostic signature" OR ABSTRACT:"gene signature" OR ABSTRACT:"risk score" '
       'OR ABSTRACT:"prognostic model" OR ABSTRACT:"risk model" OR ABSTRACT:"prognostic gene")')

QUERIES = {
    # 代谢子锚点 × 签名（第一波只覆盖了 lipid/fatty acid/glucose/amino acid 四个词）
    "W2Q1_metab_sub": f'{LUNG} AND {ANCHOR_NUT} AND {SIG} AND SRC:MED AND PUB_YEAR:[2020 TO 2026]',
    # 炎症 × 签名（放宽到 cytokine/systemic inflammation）
    "W2Q2_inf_wide": f'{LUNG} AND {ANCHOR_INF} AND {SIG} AND SRC:MED AND PUB_YEAR:[2020 TO 2026]',
    # 摘要里直接出现 "N-gene signature" + comprising/consisting（I4 最易核的一类）
    "W2Q3_abs_genelist": (f'{LUNG} AND ABSTRACT:"-gene signature" AND '
                          '(ABSTRACT:"comprising" OR ABSTRACT:"consisting" OR ABSTRACT:"including") '
                          'AND SRC:MED AND PUB_YEAR:[2020 TO 2026]'),
    # LUAD 专属 + 代谢/炎症（第一波的 TITLE 检索可能漏掉只写 LUAD 的）
    "W2Q4_luad_metab": ('(TITLE:"lung adenocarcinoma" OR ABSTRACT:"lung adenocarcinoma") AND '
                        '(TITLE:"metabolism" OR TITLE:"metabolic" OR TITLE:"inflammation" OR TITLE:"inflammatory") AND '
                        '(TITLE:signature OR TITLE:model OR TITLE:score) AND SRC:MED AND PUB_YEAR:[2020 TO 2026]'),
    # LUSC 专属（第一波 LUSC 只命中 1 个签名，样本偏少）
    "W2Q5_lusc": ('(TITLE:"lung squamous" OR ABSTRACT:"lung squamous" OR ABSTRACT:"squamous cell lung") AND '
                  '(ABSTRACT:"inflammation" OR ABSTRACT:"inflammatory" OR ABSTRACT:"metabolism" OR ABSTRACT:"metabolic") AND '
                  f'{SIG} AND SRC:MED AND PUB_YEAR:[2020 TO 2026]'),
    # TCGA 训练集明示（I2 的直接证据）
    "W2Q6_tcga_explicit": (f'{LUNG} AND (ABSTRACT:TCGA OR ABSTRACT:"The Cancer Genome Atlas") AND '
                           '(ABSTRACT:"inflammation" OR ABSTRACT:"metabolism" OR ABSTRACT:"nutritional") AND '
                           f'{SIG} AND SRC:MED AND PUB_YEAR:[2020 TO 2026]'),
}


def fetch(q, page_size=100, cursor="*"):
    params = {"query": q, "format": "json", "pageSize": str(page_size),
              "resultType": "lite", "cursorMark": cursor}
    url = BASE + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "mailto:research@example.org"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode("utf-8"))


allhits = {}
for k, q in QUERIES.items():
    cursor, got, rows, total = "*", 0, [], None
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
            })
        got += len(res)
        nc = d.get("nextCursorMark")
        if not nc or nc == cursor or got >= (total or 0) or got >= 600:
            break
        cursor = nc
        time.sleep(0.4)
    allhits[k] = {"query": q, "hitCount": total, "fetched": len(rows), "rows": rows}
    print("%-22s total=%-6s fetched=%d" % (k, total, len(rows)), flush=True)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(allhits, f, ensure_ascii=False, indent=1)
print("->", OUT)
