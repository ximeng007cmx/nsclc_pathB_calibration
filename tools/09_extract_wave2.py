# -*- coding: utf-8 -*-
"""扩容波次全文抽取：复用 03_extract_signatures.py 的抽取函数（不另写一套，避免两套逻辑漂移）。

输入：work/wave2_rule_judgement.csv 中 decision IN (INCL, REVIEW)
输出：work/extraction_wave2.json
"""
import csv, importlib.util, json, os, re, sys, time
import xml.etree.ElementTree as ET

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"

# ---- 复用第一波的抽取实现 ----
spec = importlib.util.spec_from_file_location("ex03", os.path.join(PB, "tools", "03_extract_signatures.py"))
ex03 = importlib.util.module_from_spec(spec)
sys.modules["ex03"] = ex03
spec.loader.exec_module(ex03)   # 模块级只定义函数，main() 在 __main__ 保护内

judg = list(csv.DictReader(open(f"{PB}/work/wave2_rule_judgement.csv", encoding="utf-8-sig")))
todo = [r for r in judg if r["decision"] in ("INCL", "REVIEW")]
print(f"待抽取：{len(todo)} 篇（INCL {sum(1 for r in todo if r['decision']=='INCL')}，"
      f"REVIEW {sum(1 for r in todo if r['decision']=='REVIEW')}）\n")

universe = set(l.strip() for l in
               open(os.path.join(PB, "work", "genes.tsv"), encoding="utf-8").read().splitlines()[1:])
print(f"基因全集 {len(universe)} 个\n")

recs = []
for i, r in enumerate(todo, 1):
    pmcid, pmid = r["pmcid"], r["pmid"]
    try:
        xml = ex03.fetch_xml(pmcid)
    except Exception as e:
        recs.append(dict(pmcid=pmcid, pmid=pmid, status="FETCH_FAIL", err=str(e)[:200],
                         decision=r["decision"], rule=r["rule"], title=r["title"]))
        print(f"{i:3d} {pmcid:14s} FETCH_FAIL {str(e)[:60]}", flush=True)
        continue

    try:
        t = ET.fromstring(xml.encode("utf-8"))
    except Exception as e:
        recs.append(dict(pmcid=pmcid, pmid=pmid, status="XML_FAIL", err=str(e)[:200],
                         decision=r["decision"], rule=r["rule"], title=r["title"]))
        print(f"{i:3d} {pmcid:14s} XML_FAIL", flush=True)
        continue

    def txt(xp):
        el = t.find(xp)
        return "".join(el.itertext()) if el is not None else ""

    title = txt(".//title-group/article-title")
    abstract = txt(".//abstract")
    body = txt(".//body")
    alltext = ex03.strip_tags(xml)

    n_stated, ev_n = ex03.stated_n(title, abstract, body)
    cands = ex03.gene_lists(abstract, universe) + ex03.gene_lists(body, universe)
    fg = ex03.formula_genes(alltext, universe)

    best = None
    if n_stated:
        for n, g, ev in cands:
            if n == n_stated:
                best = (n, g, ev); break
    if best is None and cands:
        best = cands[0]
    if best is None and fg:
        best = (len(fg), fg, "formula")

    if best:
        n_get, glist, ev = best
        status = "OK" if (n_stated and n_get == n_stated) else ("MISMATCH" if n_stated else "NO_STATED_N")
    else:
        n_get, glist, ev, status = 0, [], "", "NO_LIST"

    # I6：基因数 2–15
    i6_ok = 2 <= n_get <= 15
    recs.append(dict(pmcid=pmcid, pmid=pmid, year=r["year"], journal=r["journal"],
                     decision=r["decision"], rule=r["rule"], title=title,
                     n_stated=n_stated, n_stated_evidence=ev_n, n_extracted=n_get,
                     genes=";".join(glist), status=status, i6_ok=i6_ok,
                     evidence=re.sub(r"\s+", " ", ev)[:400], n_alt_lists=len(cands)))
    flag = "" if i6_ok else "  <I6 FAIL>"
    print(f"{i:3d} {pmcid:14s} stated={str(n_stated):4s} got={n_get:2d} {status:12s} "
          f"{','.join(glist)[:78]}{flag}", flush=True)
    time.sleep(0.15)

json.dump(recs, open(f"{PB}/work/extraction_wave2.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)

from collections import Counter
print("\n=== 状态汇总 ===")
for k, v in Counter((x.get("status"), x.get("i6_ok")) for x in recs).most_common():
    print(f"  {k[0]:12s} I6={k[1]}  {v}")
print(f"\n-> work/extraction_wave2.json（{len(recs)} 条）")
