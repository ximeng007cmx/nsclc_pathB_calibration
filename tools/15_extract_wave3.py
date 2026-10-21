# -*- coding: utf-8 -*-
"""修订 02 全文抽取：基因列表 + **炎症基因池来源证据**（INF-STRICT 的全文核验）。

INF-STRICT 判据（可复核）：正文里存在一句话同时满足
  (a) 含炎症词根（inflammation / inflammatory / cytokine / chemokine / acute phase / complement ...）
  (b) 含基因集来源词（were obtained from / retrieved from / downloaded from / extracted from /
      derived from / collected from / MSigDB / Molecular Signatures / GeneCards / GO term /
      Gene Ontology / KEGG / 文献汇总）
命中即记 inf_pass=True，并把证据句写进 json，供人工复核。
"""
import csv, importlib.util, json, os, re, sys, time
import xml.etree.ElementTree as ET

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
spec = importlib.util.spec_from_file_location("ex03", os.path.join(PB, "tools", "03_extract_signatures.py"))
ex03 = importlib.util.module_from_spec(spec); sys.modules["ex03"] = ex03
spec.loader.exec_module(ex03)

judg = list(csv.DictReader(open(f"{PB}/work/wave3_rule_judgement.csv", encoding="utf-8-sig")))
todo = [r for r in judg if r["decision"] in ("INCL", "REVIEW")]
print(f"待抽取 {len(todo)} 篇（INCL {sum(1 for r in todo if r['decision']=='INCL')}，"
      f"REVIEW {sum(1 for r in todo if r['decision']=='REVIEW')}）\n")

universe = set(l.strip() for l in
               open(os.path.join(PB, "work", "genes.tsv"), encoding="utf-8").read().splitlines()[1:])

INF_W = re.compile(r"inflammat\w*|cytokin\w*|chemokin\w*|acute[- ]phase|complement|"
                   r"interleukin|prostaglandin|inflammasome", re.I)
SRC_W = re.compile(r"(?:were|was)\s+(?:obtained|retrieved|downloaded|extracted|derived|collected|"
                   r"acquired|screened|identified)\s+(?:from|in)|MSigDB|Molecular Signatures|"
                   r"GeneCards|Gene Ontology|\bGO\b|KEGG|obtained from|retrieved from|"
                   r"downloaded from|extracted from|derived from", re.I)


def inf_provenance(body):
    """返回 (pass, evidence)。按句子切，找同时含炎症词与来源词的句子。"""
    if not body:
        return False, ""
    sents = re.split(r"(?<=[.;])\s+", re.sub(r"\s+", " ", body))
    hits = [s for s in sents if INF_W.search(s) and SRC_W.search(s) and len(s) < 600]
    # 优先含 "gene" 的句子（说明是基因池而不是通路讨论）
    hits.sort(key=lambda s: (0 if re.search(r"gene", s, re.I) else 1, len(s)))
    if hits:
        return True, hits[0][:400]
    return False, ""


recs = []
for i, r in enumerate(todo, 1):
    pmcid = r["pmcid"]
    try:
        xml = ex03.fetch_xml(pmcid)
        t = ET.fromstring(xml.encode("utf-8"))
    except Exception as e:
        recs.append(dict(pmcid=pmcid, pmid=r["pmid"], status="FETCH_FAIL", err=str(e)[:150],
                         decision=r["decision"], title=r["title"]))
        print(f"{i:3d} {pmcid:14s} FETCH_FAIL", flush=True); continue

    def txt(xp):
        el = t.find(xp); return "".join(el.itertext()) if el is not None else ""

    title, abstract, body = txt(".//title-group/article-title"), txt(".//abstract"), txt(".//body")
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

    ip, iev = inf_provenance(body)
    recs.append(dict(pmcid=pmcid, pmid=r["pmid"], year=r["year"], journal=r["journal"],
                     decision=r["decision"], rule=r["rule"], title=title,
                     n_stated=n_stated, n_stated_evidence=ev_n, n_extracted=n_get,
                     genes=";".join(glist), status=status, i6_ok=2 <= n_get <= 15,
                     evidence=re.sub(r"\s+", " ", ev)[:400],
                     inf_pass=ip, inf_evidence=iev))
    print(f"{i:3d} {pmcid:14s} stated={str(n_stated):4s} got={n_get:2d} {status:12s} "
          f"INF={'Y' if ip else 'n'}  {','.join(glist)[:58]}", flush=True)
    time.sleep(0.1)

json.dump(recs, open(f"{PB}/work/extraction_wave3.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)

from collections import Counter
print("\n=== 汇总 ===")
for k, v in Counter((x.get("status"), x.get("inf_pass")) for x in recs).most_common():
    print(f"  {str(k[0]):12s} INF={k[1]}  {v}")
ok = [x for x in recs if x.get("inf_pass") and x.get("status") in ("OK", "NO_STATED_N", "MISMATCH")
      and x.get("i6_ok")]
print(f"\nINF 通过且列表可抽（待人工定稿）：{len(ok)}")
for x in ok:
    print(f"  {x['pmcid']} {x['year']} {x['status']:12s} n={x['n_extracted']:2d} {x['title'][:66]}")
print(f"\n-> work/extraction_wave3.json（{len(recs)} 条）")
