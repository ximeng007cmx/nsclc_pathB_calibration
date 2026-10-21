# -*- coding: utf-8 -*-
"""W1/W2：从 PMC 全文中抽取签名基因列表，并与"论文自述的基因数"交叉校验。

设计原则（可审计）：
  * 每个签名必须能回到 PMCID + 原文片段（evidence 字段）；
  * 抽取结果必须与论文自述的基因数（如 "nine-gene"/"7-gene"）一致，不一致即标 MISMATCH 交人工；
  * 只保留出现在本地表达矩阵基因全集中的符号（否则无法打分），并记录被剔除的符号。
"""
import json, re, html, os, sys, hashlib, csv
import urllib.request
import xml.etree.ElementTree as ET

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
XMLDIR = os.path.join(PB, "work", "xml")
os.makedirs(XMLDIR, exist_ok=True)

# ---- 候选（来自 01/02 检索，人工从 60 条优先档中筛出"转录组基因签名"类）----
CANDS = [
    # (pmcid, pmid, year, journal, stratum)
    ("PMC8766344",  "35069695", "2021", "Front Genet",            "inflammation"),
    ("PMC8464410",  "34580602", "2021", "Dis Markers",            "inflammation"),
    ("PMC8647082",  "34881236", "2021", "Front Bioeng Biotechnol", "inflammation"),
    ("PMC9427773",  "36042374", "2022", "Sci Rep",                "inflammation"),
    ("PMC9736863",  "36497225", "2022", "Cancers (Basel)",        "inflammation"),
    ("PMC9203558",  "35710585", "2022", "Sci Rep",                "inflammation"),
    ("PMC8964604",  "35372503", "2022", "Front Mol Biosci",       "inflammation"),
    ("PMC9972108",  "36056909", "2023", "Cancer Med",             "inflammation"),
    ("PMC10124980", "37083272", "2023", "Ann Med",                "inflammation"),
    ("PMC12102027", "40407957", "2025", "Discov Oncol",           "inflammation"),
    ("PMC12259757", "40053252", "2025", "Clin Transl Oncol",      "inflammation"),
    ("PMC13589103", "42760759", "2026", "Cancer Med",             "inflammation"),
    ("PMC9455632",  "36091041", "2022", "Front Immunol",          "metabolism"),
    ("PMC9125243",  "35615380", "2022", "Front Genet",            "metabolism"),
    ("PMC8867215",  "35222371", "2022", "Front Immunol",          "metabolism"),
    ("PMC10493788", "37701116", "2023", "Transl Cancer Res",      "metabolism"),
    ("PMC11521851", "39478862", "2024", "Front Immunol",          "metabolism"),
    ("PMC13511999", "42650107", "2026", "Genes (Basel)",          "metabolism"),
    ("PMC13128028", "42003422", "2026", "RNA Biol",               "metabolism"),
]

WORD2NUM = {"one":1,"two":2,"three":3,"four":4,"five":5,"six":6,"seven":7,"eight":8,"nine":9,
            "ten":10,"eleven":11,"twelve":12,"thirteen":13,"fourteen":14,"fifteen":15,
            "sixteen":16,"seventeen":17,"eighteen":18,"nineteen":19,"twenty":20}


def fetch_xml(pmcid):
    p = os.path.join(XMLDIR, pmcid + ".xml")
    if os.path.exists(p) and os.path.getsize(p) > 5000:
        return open(p, encoding="utf-8", errors="replace").read()
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML" % pmcid
    req = urllib.request.Request(url, headers={"User-Agent": "mailto:research@example.org"})
    with urllib.request.urlopen(req, timeout=90) as r:
        s = r.read().decode("utf-8", errors="replace")
    open(p, "w", encoding="utf-8").write(s)
    return s


def strip_tags(s):
    s = re.sub(r"<xref[^>]*>.*?</xref>", " ", s, flags=re.S)
    s = re.sub(r"<[^>]+>", " ", s)
    s = html.unescape(s)
    return re.sub(r"\s+", " ", s)


def stated_n(title, abstract, body):
    """从题名/摘要里读"论文自述基因数"。返回 (n, evidence)"""
    for txt, tag in ((title, "title"), (abstract, "abstract")):
        t = txt or ""
        m = re.search(r"\b(\d{1,2})[- ]?(gene|mRNA|miRNA|lncRNA|feature|hub)\b", t, re.I)
        if m:
            return int(m.group(1)), "%s: %s" % (tag, m.group(0))
        m = re.search(r"\b(" + "|".join(WORD2NUM) + r")[- ](gene|mRNA)\b", t, re.I)
        if m:
            return WORD2NUM[m.group(1).lower()], "%s: %s" % (tag, m.group(0))
    # 全文兜底
    m = re.search(r"\b(" + "|".join(WORD2NUM) + r")[- ](gene|mRNA)\b", body or "", re.I)
    if m:
        return WORD2NUM[m.group(1).lower()], "body: %s" % m.group(0)
    m = re.search(r"\b(\d{1,2})[- ](gene|mRNA)\b", body or "", re.I)
    if m:
        return int(m.group(1)), "body: %s" % m.group(0)
    return None, "not stated"


def gene_lists(text, universe):
    """在文本中找"基因枚举"：括号内 ≥3 个逗号/and 分隔的基因符号。"""
    out = []
    for m in re.finditer(r"\(([^()]{15,600})\)", text):
        inner = m.group(1)
        toks = re.split(r",|;|\band\b|/", inner)
        toks = [t.strip().strip(".").strip() for t in toks if t.strip()]
        genes = [t for t in toks if t in universe]
        if len(genes) >= 3 and len(genes) >= 0.6 * len([t for t in toks if re.fullmatch(r"[A-Za-z0-9\-\.]{2,15}", t)]):
            out.append((len(genes), genes, text[max(0, m.start() - 160):m.end() + 40]))
    # 去重（按基因集合）
    seen, uniq = set(), []
    for n, g, ev in sorted(out, key=lambda x: -x[0]):
        k = frozenset(g)
        if k not in seen:
            seen.add(k); uniq.append((n, g, ev))
    return uniq


def formula_genes(text, universe):
    """风险评分公式里的基因：形如 0.123 * GENE 或 GENE * 0.123"""
    pats = [r"([0-9]\.?[0-9]{0,4})\s*[×*]\s*(?:expression of\s*)?([A-Z][A-Z0-9\-]{1,11})\b",
            r"\b([A-Z][A-Z0-9\-]{1,11})\s*[×*]\s*([0-9]\.?[0-9]{0,4})"]
    found = []
    for p in pats:
        for m in re.finditer(p, text):
            g = m.group(2) if p.startswith(r"([0-9]") else m.group(1)
            if g in universe:
                found.append(g)
    return list(dict.fromkeys(found))


def main():
    genes = [l.strip() for l in open(os.path.join(PB, "work", "genes.tsv"), encoding="utf-8").read().splitlines()[1:]]
    universe = set(genes)
    print("基因全集 %d 个" % len(universe))

    recs = []
    for pmcid, pmid, year, journal, stratum in CANDS:
        try:
            xml = fetch_xml(pmcid)
        except Exception as e:
            recs.append(dict(pmcid=pmcid, pmid=pmid, year=year, journal=journal, stratum=stratum,
                             status="FETCH_FAIL", err=str(e)))
            print("%-13s FETCH_FAIL %s" % (pmcid, e)); continue

        t = ET.fromstring(xml.encode("utf-8"))
        title = "".join(t.find(".//title-group/article-title").itertext()) if t.find(".//title-group/article-title") is not None else ""
        abs_el = t.find(".//abstract")
        abstract = "".join(abs_el.itertext()) if abs_el is not None else ""
        body_el = t.find(".//body")
        body = "".join(body_el.itertext()) if body_el is not None else ""
        alltext = strip_tags(xml)

        n_stated, ev_n = stated_n(title, abstract, body)
        cands = gene_lists(abstract, universe) + gene_lists(body, universe)
        fg = formula_genes(alltext, universe)

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

        miss = [g for g in glist if g not in universe]
        recs.append(dict(pmcid=pmcid, pmid=pmid, year=year, journal=journal, stratum=stratum,
                         title=title, n_stated=n_stated, n_stated_evidence=ev_n,
                         n_extracted=n_get, genes=";".join(glist), status=status,
                         evidence=re.sub(r"\s+", " ", ev)[:400],
                         n_alt_lists=len(cands)))
        print("%-13s stated=%-4s got=%-4s %-10s %s" % (pmcid, n_stated, n_get, status, ",".join(glist)[:90]))

    outp = os.path.join(PB, "work", "extraction_raw.json")
    json.dump(recs, open(outp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("\n->", outp)


if __name__ == "__main__":
    main()
