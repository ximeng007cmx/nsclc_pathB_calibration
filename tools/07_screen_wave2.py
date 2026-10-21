# -*- coding: utf-8 -*-
"""扩容波次初筛：去重 -> 只保留 2020+ 且有 PMCID -> 从摘要里找"基因列表"证据。

不判断"是否该纳入"，只把 I4（可核验的基因列表）在摘要层面能证实的候选挑出来，
减少后续下载全文的数量。真正的纳入仍靠 03/04/05 的全文抽取与人工核验。
"""
import csv, json, re, time, urllib.parse, urllib.request, html, os

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# ---------- 1. 读入三处既有记录 ----------
reg = {r["pmcid"] for r in csv.DictReader(open(f"{PB}/signature_registry.csv", encoding="utf-8-sig"))}
exc = {r["pmcid"] for r in csv.DictReader(open(f"{PB}/signature_registry_excluded.csv", encoding="utf-8-sig"))}
ext = {x["pmcid"] for x in json.load(open(f"{PB}/work/extraction_raw.json", encoding="utf-8"))}
wave1 = {r["pmcid"] for r in csv.DictReader(open(f"{PB}/work/candidates_dedup.csv", encoding="utf-8-sig"))}
known = reg | exc | ext
print(f"已知：纳入 {len(reg)} / 排除 {len(exc)} / 已抽取 {len(ext)}；第一波候选池 {len(wave1)}")

# ---------- 2. 合并第二波并去重 ----------
w2 = json.load(open(f"{PB}/work/candidates_wave2.json", encoding="utf-8"))
seen, cand = {}, []
for qk, blk in w2.items():
    for r in blk["rows"]:
        key = r["pmcid"] or r["pmid"]
        if not key:
            continue
        if key in seen:
            seen[key]["queries"].append(qk)
            continue
        rec = dict(r); rec["queries"] = [qk]; seen[key] = rec
print(f"第二波去重后 {len(seen)} 条")

new = [r for r in seen.values()
       if r["pmcid"] and str(r["year"]).isdigit() and int(r["year"]) >= 2020
       and r["pmcid"] not in known]
print(f"排除已知 + 2020 起 + 有 PMCID 后：{len(new)} 条（其中第一波候选池内 {sum(1 for r in new if r['pmcid'] in wave1)} 条，"
      f"真正新增 {sum(1 for r in new if r['pmcid'] not in wave1)} 条）")
new.sort(key=lambda r: (-int(r["year"]), r["journal"] or ""))


# ---------- 3. 批量取摘要（一次 25 个 EXT_ID） ----------
def fetch_abstracts(pmids):
    out = {}
    for i in range(0, len(pmids), 25):
        chunk = pmids[i:i + 25]
        q = " OR ".join(f"EXT_ID:{p}" for p in chunk)
        params = {"query": q, "format": "json", "pageSize": "100",
                  "resultType": "core", "cursorMark": "*"}
        url = BASE + "?" + urllib.parse.urlencode(params)
        req = urllib.request.Request(url, headers={"User-Agent": "mailto:research@example.org"})
        try:
            with urllib.request.urlopen(req, timeout=120) as rr:
                d = json.loads(rr.read().decode("utf-8"))
        except Exception as e:
            print("  ! batch failed", i, e, flush=True)
            time.sleep(3); continue
        for it in d.get("resultList", {}).get("result", []):
            ab = it.get("abstractText") or ""
            out[it.get("pmid")] = html.unescape(re.sub(r"<[^>]+>", " ", ab))
        time.sleep(0.5)
    return out


pmids = [r["pmid"] for r in new if r["pmid"]]
print(f"取摘要：{len(pmids)} 篇 →", flush=True)
abs_map = fetch_abstracts(pmids)
print(f"拿到摘要 {len(abs_map)} 篇")

# ---------- 4. 摘要里的基因列表探测 ----------
GENELIKE = re.compile(r"^[A-Z][A-Z0-9]{1,8}(?:[-/][A-Z0-9]{1,6})?$")
STOP = {"TCGA", "NSCLC", "LUAD", "LUSC", "OS", "PFS", "ROC", "AUC", "HR", "CI", "GO", "KEGG",
        "LASSO", "RNA", "DNA", "MRNA", "MIRNA", "LNRNa".upper(), "COX", "DEGS", "DEG",
        "GSEA", "GDSC", "IC50", "TME", "TIME", "EMT", "TMB", "MSI", "FDA", "NCCN",
        "AJCC", "TNM", "IHC", "PCR", "QPCR", "SD", "SE", "IQR", "BMI", "CRP", "ALB",
        "NLR", "PLR", "LMR", "PNI", "CONUT", "GPS", "SII", "SIRI", "ALI", "HALP",
        "NGS", "WGCNA", "CIBERSORT", "ESTIMATE", "GSVA", "SSGSEA", "FDR", "BH"}
NGENE_CUE = re.compile(r"(\d+|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|"
                       r"thirteen|fourteen|fifteen)[\s-]*(?:gene|mRNA|IRG|DEMG|hub gene)", re.I)


def gene_list_from_abs(ab):
    """返回 (best_list, n_tokens, cue)。best_list 为最长的一段连续逗号/and 连接的符号串。"""
    if not ab:
        return [], 0, None
    t = re.sub(r"[()\[\]]", " ", ab)
    t = re.sub(r"\s+", " ", t)
    tokens = re.split(r"[,;]| and ", t)
    best, cur = [], []
    for tok in tokens:
        cand = tok.strip().strip(".").strip()
        words = cand.split()
        if len(words) == 1 and GENELIKE.match(cand) and cand.upper() not in STOP \
                and not re.match(r"^\d", cand):
            cur.append(cand)
        else:
            if len(cur) > len(best):
                best = cur
            cur = []
    if len(cur) > len(best):
        best = cur
    cue = NGENE_CUE.search(ab)
    return best, len(best), (cue.group(0) if cue else None)


rows = []
for r in new:
    ab = abs_map.get(r["pmid"], "")
    lst, n, cue = gene_list_from_abs(ab)
    has_sig = bool(re.search(r"signature|risk score|risk model|prognostic model|nomogram", ab, re.I))
    has_anchor = bool(re.search(r"inflammat|nutrition|metabol|lipid|fatty acid|amino acid|glucose|glycolys", ab, re.I))
    has_lung = bool(re.search(r"lung|NSCLC|LUAD|LUSC|pulmonary", ab, re.I))
    rows.append({
        "pmcid": r["pmcid"], "pmid": r["pmid"], "doi": r["doi"], "year": r["year"],
        "journal": r["journal"], "title": r["title"],
        "abstract_gene_list": ",".join(lst), "n_in_abstract": n,
        "n_gene_cue": cue or "",
        "sig": has_sig, "anchor": has_anchor, "lung": has_lung,
        "abs_len": len(ab), "queries": "|".join(r["queries"]),
    })

rows.sort(key=lambda x: (-x["n_in_abstract"], -int(x["year"])))
os.makedirs(f"{PB}/work", exist_ok=True)
with open(f"{PB}/work/candidates_wave2_screened.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

tier_a = [r for r in rows if r["sig"] and r["anchor"] and r["lung"] and 2 <= r["n_in_abstract"] <= 15]
tier_b = [r for r in rows if r["sig"] and r["anchor"] and r["lung"] and r["n_in_abstract"] < 2]
print(f"\n三条件全中且摘要含 2–15 个基因符号（Tier A，最可能满足 I4）：{len(tier_a)}")
for r in tier_a:
    print(f"  {r['pmcid']} {r['year']} n={r['n_in_abstract']:2d} cue={r['n_gene_cue'][:22]:22s} {r['title'][:64]}")
print(f"三条件全中但摘要无列表（Tier B，需查全文）：{len(tier_b)}")
print(f"\n-> work/candidates_wave2_screened.csv（{len(rows)} 行）")
