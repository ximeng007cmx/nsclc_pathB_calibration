# -*- coding: utf-8 -*-
"""
29_render_manuscript_html.py — render the current manuscript (MS_draft_v1.md,
revision 07, sha256[:16] 67dc50ec12860396) into a single readable HTML file.

Read-only rendering for review; the .md master is NOT modified.

Appendices embedded for reading convenience:
  - Figures: the four SVGs (inlined) with their English legends from
    manuscript/figure_captions.md (CJK provenance note lines dropped)
  - Tables: manuscript/tables.md rendered (CJK provenance lines dropped)
"""
import io
import os
import re
import hashlib

import markdown

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MS = os.path.join(BASE, "manuscript", "MS_draft_v1.md")
CAPS = os.path.join(BASE, "manuscript", "figure_captions.md")
TBLS = os.path.join(BASE, "manuscript", "tables.md")
FIGDIR = os.path.join(BASE, "manuscript", "figures")
OUT = os.path.join(BASE, "manuscript", "MS_draft_v1.html")

md_text = io.open(MS, encoding="utf-8").read()

# integrity: render the revision-07 version, fail loudly otherwise
h = hashlib.sha256(io.open(MS, "rb").read()).hexdigest()[:16]
assert h == "67dc50ec12860396", "MS_draft_v1.md changed since revision 07: %s" % h

body = markdown.markdown(md_text, extensions=["tables"])

# --- figure legends (drop CJK provenance lines) ------------------------------
caps_raw = io.open(CAPS, encoding="utf-8").read()
caps_paras = [
    p for p in re.split(r"\n\s*\n", caps_raw)
    if p.startswith("**Fig") and not re.search(r"[\u4e00-\u9fff]", p)
]
assert len(caps_paras) == 4, len(caps_paras)
caps_html = [markdown.markdown(p, extensions=["tables"]) for p in caps_paras]

# --- figures inlined ----------------------------------------------------------
fig_files = ["Fig1_flow.svg", "Fig2_forest.svg", "Fig3_null_controls.svg",
             "Fig4_cohorts_strata.svg"]
fig_html = []
for f in fig_files:
    svg = io.open(os.path.join(FIGDIR, f), encoding="utf-8").read()
    assert "<svg" in svg and f.startswith("Fig")
    fig_html.append('<div class="figbox">%s</div>' % svg)

# --- tables (drop CJK provenance lines) ---------------------------------------
tbl_lines = [l for l in io.open(TBLS, encoding="utf-8").read().split("\n")
             if not re.search(r"[\u4e00-\u9fff]", l)]
tables_html = markdown.markdown("\n".join(tbl_lines), extensions=["tables"])

page = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>MS_draft_v1 (revision 07) — Immune and stromal content does not explain …</title>
<style>
  body { font-family: "Times New Roman", Times, serif; font-size: 12.5pt; line-height: 1.55;
         color: #1a1a1a; background: #ffffff; margin: 0; }
  .page { max-width: 880px; margin: 0 auto; padding: 48px 56px 80px; }
  h1 { font-size: 17pt; line-height: 1.35; margin: 0 0 6px; }
  h2 { font-size: 14pt; border-bottom: 1px solid #999; padding-bottom: 3px; margin-top: 34px; }
  h3 { font-size: 12.5pt; margin-top: 24px; }
  blockquote { margin: 10px 0; padding: 8px 14px; border-left: 3px solid #b8b8b8;
               background: #f6f6f6; font-size: 10.5pt; color: #444; }
  blockquote p { margin: 3px 0; }
  .meta { color: #555; font-size: 10.5pt; margin-bottom: 22px; }
  table { border-collapse: collapse; font-size: 9.5pt; margin: 12px 0; width: 100%; }
  th, td { border: 1px solid #bbb; padding: 3px 6px; text-align: left; vertical-align: top; }
  th { background: #efefef; }
  .tblscroll { overflow-x: auto; margin: 8px 0 20px; }
  .figbox { margin: 20px 0 6px; text-align: center; }
  .figbox svg { max-width: 100%; height: auto; }
  .figcap { font-size: 10.5pt; color: #333; margin: 0 0 26px; }
  .figcap p { margin: 4px 0; }
  .appendix h2 { margin-top: 46px; }
  hr { border: 0; border-top: 1px solid #ccc; margin: 30px 0; }
  code { font-family: Consolas, monospace; font-size: 0.9em; background: #f2f2f2; padding: 0 3px; }
  @page { size: A4; margin: 18mm 16mm; }
  @media print {
    body { font-size: 10.5pt; }
    .page { max-width: none; padding: 0; }
    h2 { margin-top: 22px; }
    .figbox { margin: 12px 0 4px; }
    .figcap { margin: 0 0 14px; font-size: 9pt; }
    table { font-size: 7.5pt; }
    th, td { padding: 2px 4px; }
    .tblscroll { overflow: visible; }
    blockquote { font-size: 9pt; }
    a { color: inherit; text-decoration: none; }
  }
</style>
</head>
<body><div class="page">
__BODY__
<div class="appendix">
<h2>Appendix A. Figures</h2>
__FIGS__
<h2>Appendix B. Tables</h2>
<div class="tblscroll">
__TBLS__
</div>
</div>
</div></body></html>"""

figs_block = "".join(
    '<div class="figbox">%s</div><div class="figcap">%s</div>' % (svg, cap)
    for svg, cap in zip(fig_html, caps_html)
)
page = page.replace("__BODY__", body).replace("__FIGS__", figs_block).replace("__TBLS__", tables_html)

io.open(OUT, "w", encoding="utf-8").write(page)
print("written:", OUT, "(%d bytes)" % os.path.getsize(OUT))
print("source hash verified:", h)
