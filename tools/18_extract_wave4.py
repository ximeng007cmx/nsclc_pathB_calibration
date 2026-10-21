# -*- coding: utf-8 -*-
"""修订 03：对 20 条 CIRC 候选做带证据的全文抽取（含基因池来源句，供分层归属判断）。"""
import csv, importlib.util, json, os, re, sys, time
import xml.etree.ElementTree as ET

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
spec = importlib.util.spec_from_file_location("ex03", os.path.join(PB, "tools", "03_extract_signatures.py"))
ex03 = importlib.util.module_from_spec(spec); sys.modules["ex03"] = ex03
spec.loader.exec_module(ex03)

judg = [r for r in csv.DictReader(open(f"{PB}/work/wave3_rule_judgement.csv", encoding="utf-8-sig"))
        if r["rule"] == "CIRC"]
print(f"CIRC 候选 {len(judg)} 篇\n")
universe = set(l.strip() for l in
               open(os.path.join(PB, "work", "genes.tsv"), encoding="utf-8").read().splitlines()[1:])

POOL_SRC = re.compile(r"(?:were|was)\s+(?:obtained|retrieved|downloaded|extracted|derived|collected|"
                      r"acquired|screened|selected|identified)\s+from|MSigDB|Molecular Signatures|"
                      r"GeneCards|Gene Ontology|\bGO terms?\b|KEGG|Harmonosome|from the literature", re.I)
GENESET_W = re.compile(r"gene[s]?\s+(?:set|list|pool)|genes?\b.{0,60}(?:set|list)|candidate genes?"
                       r"|hub genes?|\bIRGs?\b|\bIRDEGs?\b|\bDEIRGs?\b|immune[- ]related genes?"
                       r"|inflammat\w*|cytokin\w*|chemokin\w*|costimulat\w*|immunogenic", re.I)


def pool_evidence(body):
    if not body:
        return ""
    sents = re.split(r"(?<=[.;])\s+", re.sub(r"\s+", " ", body))
    hits = [x for x in sents if GENESET_W.search(x) and POOL_SRC.search(x) and len(x) < 600]
    hits.sort(key=len)
    return (hits[0] if hits else "")[:420]


recs = []
for i, r in enumerate(judg, 1):
    pmcid = r["pmcid"]
    try:
        xml = ex03.fetch_xml(pmcid)
        t = ET.fromstring(xml.encode("utf-8"))
    except Exception as e:
        recs.append(dict(pmcid=pmcid, pmid=r["pmid"], status="FETCH_FAIL", err=str(e)[:120],
                         title=r["title"]))
        print(f"{i:2d} {pmcid:14s} FETCH_FAIL", flush=True); continue

    def tx(xp):
        el = t.find(xp); return "".join(el.itertext()) if el is not None else ""

    ti, ab, bo = tx(".//title-group/article-title"), tx(".//abstract"), tx(".//body")
    n_st, ev_n = ex03.stated_n(ti, ab, bo)
    cands = ex03.gene_lists(ab, universe) + ex03.gene_lists(bo, universe)
    fg = ex03.formula_genes(ex03.strip_tags(xml), universe)
    best = None
    if n_st:
        for n, g, ev in cands:
            if n == n_st:
                best = (n, g, ev); break
    if best is None and cands:
        best = cands[0]
    if best is None and fg:
        best = (len(fg), fg, "formula")
    if best:
        n_get, glist, ev = best
        status = "OK" if (n_st and n_get == n_st) else ("MISMATCH" if n_st else "NO_STATED_N")
    else:
        n_get, glist, ev, status = 0, [], "", "NO_LIST"
    pe = pool_evidence(bo) or pool_evidence(ex03.strip_tags(xml))
    recs.append(dict(pmcid=pmcid, pmid=r["pmid"], year=r["year"], journal=r["journal"],
                     title=ti, n_stated=n_st, n_stated_evidence=ev_n, n_extracted=n_get,
                     genes=";".join(glist), status=status, i6_ok=2 <= n_get <= 15,
                     evidence=re.sub(r"\s+", " ", ev)[:400], pool_evidence=pe))
    print(f"{i:2d} {pmcid:14s} stated={str(n_st):4s} got={n_get:2d} {status:12s} "
          f"{','.join(glist)[:54]}", flush=True)
    time.sleep(0.1)

json.dump(recs, open(f"{PB}/work/extraction_wave4.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
from collections import Counter
print("\n=== 汇总 ===")
for k, v in Counter(x.get("status") for x in recs).most_common():
    print(f"  {k:12s} {v}")
usable = [x for x in recs if x.get("i6_ok") and x.get("status") in ("OK", "NO_STATED_N", "MISMATCH")]
print(f"\n可抽且 2–15 基因：{len(usable)}")
for x in usable:
    print(f"\n{x['pmcid']}  {x.get('year')}  {x.get('journal')}")
    print(f"  {x.get('title','')[:130]}")
    print(f"  stated={x.get('n_stated')} got={x.get('n_extracted')} {x.get('status')}")
    print(f"  GENES: {x.get('genes','')}")
    print(f"  池来源: {(x.get('pool_evidence') or '')[:240]}")
print(f"\n-> work/extraction_wave4.json")
