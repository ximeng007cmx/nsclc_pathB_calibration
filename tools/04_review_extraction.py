# -*- coding: utf-8 -*-
"""人工复核用：对每篇打印候选基因列表 + 原文片段，供逐条判定（不改数据）。"""
import json, re, html, os, sys
import xml.etree.ElementTree as ET

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
XMLDIR = os.path.join(PB, "work", "xml")
TARGETS = sys.argv[1:] if len(sys.argv) > 1 else []

genes = [l.strip() for l in open(os.path.join(PB, "work", "genes.tsv"), encoding="utf-8").read().splitlines()[1:]]
universe = set(genes)

def strip_tags(s):
    s = re.sub(r"<[^>]+>", " ", s)
    return re.sub(r"\s+", " ", html.unescape(s))

for pmcid in TARGETS:
    p = os.path.join(XMLDIR, pmcid + ".xml")
    xml = open(p, encoding="utf-8", errors="replace").read()
    t = ET.fromstring(xml.encode("utf-8"))
    title = "".join(t.find(".//title-group/article-title").itertext()) if t.find(".//title-group/article-title") is not None else ""
    abs_el = t.find(".//abstract")
    abstract = "".join(abs_el.itertext()) if abs_el is not None else ""
    print("\n" + "=" * 100)
    print(pmcid, "|", re.sub(r"\s+", " ", title)[:150])
    print("-" * 100)
    print("ABSTRACT:", re.sub(r"\s+", " ", abstract)[:1400])

    body = strip_tags(xml)
    for kw in ["risk score", "riskScore", "prognostic model was", "signature consisting", "constructed a"]:
        for m in list(re.finditer(kw, body, re.I))[:2]:
            print("\n  [%s] ...%s..." % (kw, body[max(0, m.start() - 200):m.end() + 450]))
