# -*- coding: utf-8 -*-
"""把四组检索去重，按"是否 NSCLC / 是否签名类 / 是否有 PMC 全文"分档输出，供人工判定。"""
import json, re, csv, os

RAW = r"D:/WorkBuddy_Y02/out_ali_directions/pathB/work/candidates_raw.json"
OUT_CSV = r"D:/WorkBuddy_Y02/out_ali_directions/pathB/work/candidates_dedup.csv"

d = json.load(open(RAW, encoding="utf-8"))
seen = {}
for qk, blk in d.items():
    for r in blk["rows"]:
        key = r.get("pmid") or r.get("doi") or r.get("title")
        if key in seen:
            seen[key]["queries"].append(qk)
        else:
            r["queries"] = [qk]
            seen[key] = r

LUNG = re.compile(r"lung|NSCLC|LUAD|LUSC|adenocarcinoma", re.I)
SIG = re.compile(r"signature|model|score|index|nomogram|panel", re.I)
ANCHOR = re.compile(r"nutrition|inflammatory|inflammation|immune|metabol|lipid|amino", re.I)

rows = []
for r in seen.values():
    t = r["title"]
    rows.append({
        "pmid": r.get("pmid"), "pmcid": r.get("pmcid"), "doi": r.get("doi"), "year": r.get("year"),
        "journal": r.get("journal"), "title": t,
        "lung": bool(LUNG.search(t)), "sig": bool(SIG.search(t)), "anchor": bool(ANCHOR.search(t)),
        "has_pmcid": bool(r.get("pmcid")),
        "queries": ";".join(sorted(set(r["queries"]))),
    })

# 优先：有 PMCID 且 题名含肺 + 签名 + 营养/炎症锚点
rows.sort(key=lambda x: (not (x["has_pmcid"] and x["lung"] and x["sig"] and x["anchor"]), x["year"] or ""))

os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
with open(OUT_CSV, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print("去重后 %d 条；有 PMCID %d 条；三条件全中且有 PMCID %d 条"
      % (len(rows), sum(r["has_pmcid"] for r in rows),
         sum(1 for r in rows if r["has_pmcid"] and r["lung"] and r["sig"] and r["anchor"])))
print("\n=== 优先档（有 PMCID + 题名含肺 + 签名 + 营养/炎症锚点）===")
n = 0
for r in rows:
    if r["has_pmcid"] and r["lung"] and r["sig"] and r["anchor"]:
        n += 1
        print("%2d. [%s] %s | %s | PMID %s | %s" % (n, r["year"], r["title"][:120], r["journal"], r["pmid"], r["pmcid"]))
