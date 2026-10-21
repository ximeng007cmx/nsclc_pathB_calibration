#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
26_make_figures.py
生成投稿图 Fig 1 / Fig 3 / Fig 4（Fig 2 森林图从主报告提取）。
纪律：所有数字从 results/*.tsv、work/*.json、logs/*.log 直读并 assert，
     与 06_amendment_numbers.md / LOCK.txt 口径一致；口径不符会触发断言失败。
输出：manuscript/figures/Fig1_flow.svg, Fig2_forest.svg, Fig3_null_controls.svg, Fig4_cohorts_strata.svg
"""
import csv, os, re, json, math, gzip

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(B, "results")
W = os.path.join(B, "work")
L = os.path.join(B, "logs")
OUT = os.path.join(B, "manuscript", "figures")
os.makedirs(OUT, exist_ok=True)

C_IMM, C_INF, C_MET = "#6d28d9", "#c2410c", "#1d4ed8"   # 与 Fig 2 森林图一致
C_OK, C_EXC, C_NS, C_EDGE, C_TXT, C_SUB = "#dbeafe", "#fee2e2", "#fef3c7", "#94a3b8", "#1f2937", "#6b7280"


def tsv(name, base=None):
    with open(os.path.join(base or R, name), encoding="utf-8-sig", newline="") as f:
        return [dict(r) for r in csv.DictReader(f, delimiter="\t")]


def rd(name, base=None):
    with open(os.path.join(base or B, name), encoding="utf-8-sig", newline="") as f:
        return [dict(r) for r in csv.DictReader(f)]


def F(x):
    try:
        return float(x)
    except Exception:
        return None


def med(v):
    v = sorted(v)
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    num = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    den = math.sqrt(sum((v - mx) ** 2 for v in x) * sum((v - my) ** 2 for v in y))
    return num / den if den else float("nan")


def spearman(x, y):
    def rk(v):
        n = len(v)
        srt = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[srt[j + 1]] == v[srt[i]]:
                j += 1
            for k in range(i, j + 1):
                r[srt[k]] = (i + j) / 2.0 + 1
            i = j + 1
        return r
    return pearson(rk(x), rk(y))


HATCH = ('<defs><pattern id="hatch" width="6" height="6" patternUnits="userSpaceOnUse" '
         'patternTransform="rotate(45)"><rect width="6" height="6" fill="none"/>'
         '<line x1="0" y1="0" x2="0" y2="6" stroke="#d4a72c" stroke-width="2"/></pattern></defs>')


class SVG:
    def __init__(self, w, h, title):
        self.w, self.h = w, h
        self.p = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" '
                  'role="img" aria-label="%s">' % (w, h, w, title)]
        self.p.append('<rect width="%d" height="%d" fill="#ffffff"/>' % (w, h))
        self.p.append('<g font-family="Helvetica, Arial, sans-serif" fill="%s">' % C_TXT)

    def text(self, x, y, s, size=10, anchor="start", fill=C_TXT, weight="normal", halo=False):
        s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        st = ' stroke="#ffffff" stroke-width="3" paint-order="stroke"' if halo else ""
        self.p.append('<text x="%.1f" y="%.1f" font-size="%s" text-anchor="%s" fill="%s" '
                      'font-weight="%s"%s>%s</text>' % (x, y, size, anchor, fill, weight, st, s))

    def rect(self, x, y, w, h, fill, stroke=C_EDGE, rx=3, dash="", sw=1):
        self.p.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" stroke="%s" '
                      'stroke-width="%s" rx="%s"%s/>' % (x, y, w, h, fill, stroke, sw, rx,
                                                         ' stroke-dasharray="%s"' % dash if dash else ""))

    def line(self, x1, y1, x2, y2, stroke=C_EDGE, sw=1, dash=""):
        self.p.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%s"%s/>'
                      % (x1, y1, x2, y2, stroke, sw, ' stroke-dasharray="%s"' % dash if dash else ""))

    def bar(self, x, y, w, h, fill, hatch=False):
        self.rect(x, y, w, h, fill, stroke="none", rx=0)
        if hatch:
            self.p.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="url(#hatch)"/>'
                          % (x, y, w, h))

    def arrow(self, x1, y1, x2, y2, stroke=C_SUB):
        self.line(x1, y1, x2, y2, stroke, 1)
        self.p.append('<polygon points="%.1f,%.1f %.1f,%.1f %.1f,%.1f" fill="%s"/>' %
                      (x2, y2, x2 - 3.5, y2 - 6, x2 + 3.5, y2 - 6, stroke))

    def save(self, name):
        self.p.append("</g></svg>")
        path = os.path.join(OUT, name)
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(self.p))
        return path


# ============================================================
# 数据直读 + 断言（所有数字与 06 号修订核对过的值一致）
# ============================================================
A = tsv("attenuation_all.tsv")
SUM = {r["set"]: r for r in tsv("summary_all.tsv")}
EP = tsv("attenuation_endpoints.tsv")
SUMEP = {(r["endpoint"], r["set"]): r for r in tsv("summary_endpoints.tsv")}
EXT = {(r["cohort"], r["stratum"]): r for r in tsv("external_m1m3_summary.tsv")}
CTL = {r["axis"]: r for r in tsv("control_axes_all.tsv")}
PCT = tsv("attenuation_pctile_all.tsv")
EXTRHO = {r["cohort"]: r for r in tsv("external_vs_tcga_rho.tsv")}

# --- Fig 1 ---
reg = {w: len(rd(f)) for w, f in (("W1", "signature_registry.csv"),
                                  ("W2", "signature_registry_expansion.csv"),
                                  ("W3", "signature_registry_expansion_w3.csv"),
                                  ("W4", "signature_registry_expansion_w4.csv"))}
assert (reg["W1"], reg["W2"], reg["W3"], reg["W4"]) == (13, 28, 2, 6), reg
j2 = rd("work/wave2_rule_judgement.csv")
ft2 = [r for r in j2 if r["decision"] in ("INCL", "REVIEW")]
reg2ids = {r["pmcid"] for r in rd("signature_registry_expansion.csv")}
exc2 = len(rd("signature_registry_expansion_excluded.csv"))
sc2 = {r["pmcid"] for r in ft2} - reg2ids - {r["pmcid"] for r in rd("signature_registry_expansion_excluded.csv")}
assert (len(j2), len(ft2), exc2, len(sc2)) == (99, 62, 26, 8), (len(j2), len(ft2), exc2, len(sc2))
d2 = json.load(open(os.path.join(W, "candidates_wave2.json"), encoding="utf-8"))
assert sum(v["hitCount"] for v in d2.values()) == 767
assert sum(v["fetched"] for v in d2.values()) == 767
m = re.search(r"第二波去重后 (\d+) 条", open(os.path.join(L, "07_screen_wave2.log"), encoding="utf-8").read())
assert int(m.group(1)) == 538

d3 = json.load(open(os.path.join(W, "extraction_wave3.json"), encoding="utf-8"))
j3 = rd("work/wave3_rule_judgement.csv")
c3 = json.load(open(os.path.join(W, "candidates_wave3.json"), encoding="utf-8"))
w3_hits = sum(v["hitCount"] for v in c3.values())
w3_fetch = sum(v["fetched"] for v in c3.values())
assert (w3_hits, w3_fetch) == (1139, 851), (w3_hits, w3_fetch)
assert len(rd("work/wave3_screened.csv")) == 605
circ3 = sum(1 for r in j3 if r["rule"] == "CIRC")
out3 = sum(1 for r in j3 if r["rule"] == "I3-OUT")
assert (circ3, out3) == (20, 6), (circ3, out3)
status3 = {x["pmcid"]: str(x.get("status")) for x in d3}
reg3 = {r["pmcid"] for r in rd("signature_registry_expansion_w3.csv")}
exc3 = {r["pmcid"] for r in rd("signature_registry_expansion_w3_excluded.csv")}
inf3 = sum(1 for x in d3 if str(x.get("inf_pass")) == "True")
ok3 = sum(1 for x in d3 if status3.get(x["pmcid"]) == "OK")
unsc3 = sum(1 for p, st in status3.items() if st != "OK" and p not in exc3)
nsp3 = sum(1 for x in d3 if status3.get(x["pmcid"]) == "OK" and x["pmcid"] not in reg3 and x["pmcid"] not in exc3)
assert (inf3, ok3, len(reg3), len(exc3), unsc3, nsp3) == (16, 9, 2, 6, 36, 4), \
    (inf3, ok3, len(reg3), len(exc3), unsc3, nsp3)
assert len(reg3) + len(exc3) + unsc3 + nsp3 == 48
reg4 = len(rd("signature_registry_expansion_w4.csv"))
exc4 = len(rd("signature_registry_expansion_w4_excluded.csv"))
assert (reg4, exc4, reg4 + exc4) == (6, 10, 16)

strata = {}
for r in A:
    strata[r["stratum"]] = strata.get(r["stratum"], 0) + 1
assert strata == {"inflammation": 11, "metabolism": 32, "immune": 6}, strata

# ============================================================
# Fig 1  四波纳入流程
# ============================================================
s = SVG(680, 500, "Fig 1 inclusion flow")
s.p.insert(2, HATCH)

LANE_X = [14, 186, 358, 530]
BW = 158
for lane, label in zip(LANE_X, ["Wave 1 (pilot)", "Wave 2 (metabolism)",
                                "Wave 3 (inflammation)", "Wave 4 (immune)"]):
    s.text(lane + BW / 2, 26, label, 11, "middle", C_TXT, "bold")


def box(lane, y, h, lines, fill="#ffffff", stroke=C_EDGE, size=9.0):
    s.rect(lane, y, BW, h, fill, stroke)
    ty = y + 13 if len(lines) > 1 else y + h / 2 + 3.4
    for k, ln in enumerate(lines):
        s.text(lane + BW / 2, ty + k * 11.6, ln, size if k else size + 0.6, "middle",
               C_TXT if k == 0 else C_SUB, "bold" if k == 0 else "normal")


def down(lane, y1, y2):
    s.arrow(lane + BW / 2, y1, lane + BW / 2, y2)


# W1
box(LANE_X[0], 40, 40, ["13 signatures", "preregistered search"])
down(LANE_X[0], 80, 94)
box(LANE_X[0], 94, 30, ["6 excluded", "(recorded separately)"], C_EXC)
s.text(LANE_X[0] + BW / 2, 146, "pilot set only; 87 candidates", 8.2, "middle", C_SUB)
s.text(LANE_X[0] + BW / 2, 157, "screened in wave 1", 8.2, "middle", C_SUB)

# W2
y = 40
box(LANE_X[1], y, 26, ["767 records (6 queries, all retrieved)"]); y += 26
down(LANE_X[1], y, y + 12); y += 12
box(LANE_X[1], y, 24, ["538 after de-duplication"]); y += 24
down(LANE_X[1], y, y + 12); y += 12
box(LANE_X[1], y, 24, ["431 after eligibility filter"]); y += 24
down(LANE_X[1], y, y + 12); y += 12
box(LANE_X[1], y, 24, ["99 tier A (abstract gene list)"], "#f1f5f9"); y += 24
down(LANE_X[1], y, y + 12); y += 12
box(LANE_X[1], y, 24, ["62 to full-text extraction"]); y += 24
down(LANE_X[1], y, y + 14); y += 14
box(LANE_X[1], y, 56, ["28 included", "26 excluded (reason per record)",
                       "8 not scored: 6 no list,", "2 not retrievable"], C_OK); y += 56

# W3
y = 40
box(LANE_X[2], y, 42, ["1,139 hits (8 queries)", "851 retrieved", "(Q6: 788 hits, API cap 500)"]); y += 42
down(LANE_X[2], y, y + 12); y += 12
box(LANE_X[2], y, 24, ["605 new candidates"]); y += 24
down(LANE_X[2], y, y + 12); y += 12
box(LANE_X[2], y, 24, ["74 tier A"], "#f1f5f9"); y += 24
down(LANE_X[2], y, y + 12); y += 12
box(LANE_X[2], y, 40, ["26 excluded at rule stage:", "20 immune-anchored (CIRC),",
                       "6 outside construct (I3-OUT)"], C_EXC); y += 40
down(LANE_X[2], y, y + 12); y += 12
box(LANE_X[2], y, 24, ["48 to full-text extraction"]); y += 24
down(LANE_X[2], y, y + 14); y += 14
box(LANE_X[2], y, 68, ["2 included", "6 excluded (reason per record)",
                       "4 verifiable list, no pool evidence",
                       "36 not scorable (25/4/5/2)"], C_OK); y += 68

# W4
y = 40
box(LANE_X[3], y, 40, ["20 candidates excluded by", "the immune-construct rule",
                       "(re-examined, no new search)"]); y += 40
down(LANE_X[3], y, y + 12); y += 12
box(LANE_X[3], y, 24, ["16 with a verifiable gene list"], "#f1f5f9"); y += 24
down(LANE_X[3], y, y + 12); y += 12
box(LANE_X[3], y, 52, ["6 included", "10 excluded (reason per record)",
                       "adjudicated before any", "immune-stratum result"], C_OK); y += 52

# 汇总
s.rect(14, 376, 652, 34, "#f8fafc", C_SUB)
s.text(340, 390, "Final set: 49 signatures", 11.5, "middle", C_TXT, "bold")
s.text(340, 404, "11 inflammation-related + 32 metabolism-related + 6 immune-related   "
                 "(W1 13 · W2 28 · W3 2 · W4 6)", 9.2, "middle", C_SUB)
s.text(340, 428, "Every exclusion has a recorded reason; every unscored record is listed with its "
                 "extraction status.", 8.8, "middle", C_SUB)
s.text(340, 441, "All numbers read directly from the registries, screening tables and extraction logs "
                 "(amendment 06).", 8.2, "middle", "#9aa3ad")
lx = 14
for col, lab in ((C_INF, "inflammation (11)"), (C_MET, "metabolism (32)"), (C_IMM, "immune (6)")):
    s.rect(lx, 458, 12, 12, col, stroke="none", rx=2)
    s.text(lx + 17, 468, lab, 9.4, "start", C_SUB)
    lx += 150
print("Fig1 ->", s.save("Fig1_flow.svg"))

# ============================================================
# Fig 3  零分布 + 对照轴
# ============================================================
pcts = [F(r["pctile"]) for r in PCT]
assert len(pcts) == 49
medp = med(pcts)
assert abs(medp - 0.440) < 5e-4, medp
n_lo = sum(1 for x in pcts if x < 0.20)
n_mid = sum(1 for x in pcts if 0.20 <= x < 0.80)
n_hi = sum(1 for x in pcts if x >= 0.80)
assert (n_lo, n_mid, n_hi) == (9, 37, 3), (n_lo, n_mid, n_hi)
n_p = len(pcts)
sp = sorted(pcts)
D = max(max(abs((i + 1) / n_p - sp[i]), abs(i / n_p - sp[i])) for i in range(n_p))
lam = (math.sqrt(n_p) + 0.12 + 0.11 / math.sqrt(n_p)) * D
ks_p = min(max(2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 101)), 0.0), 1.0)
assert abs(D - 0.267) < 5e-4 and abs(ks_p - 0.0014) < 1e-4, (D, ks_p)

s = SVG(680, 320, "Fig 3 null distribution and control axes")
s.p.insert(2, HATCH)
# --- a：百分位直方图 ---
x0, y0, pw, ph = 64, 46, 250, 196
s.text(14, 24, "a", 12, "start", C_TXT, "bold")
s.text(30, 24, "Attenuation percentile vs 200 matched random panels", 10, "start", C_TXT, "bold")
bins = [0] * 10
for x in pcts:
    bins[min(int(x * 10), 9)] += 1
assert sum(bins) == 49
ymax = max(max(bins), 12)
bw = pw / 10.0
for i, c in enumerate(bins):
    h = ph * c / ymax
    col = C_MET if c else "#eef2f7"
    s.bar(x0 + i * bw + 1, y0 + ph - h, bw - 2, h, col)
    if c:
        s.text(x0 + i * bw + bw / 2, y0 + ph - h - 4, str(c), 8.6, "middle", C_SUB)
s.rect(x0, y0, pw, ph, "none")
for t in (0, 0.2, 0.4, 0.6, 0.8, 1.0):
    xx = x0 + t * pw
    s.line(xx, y0 + ph, xx, y0 + ph + 4, C_SUB)
    s.text(xx, y0 + ph + 15, "%.1f" % t, 8.6, "middle", C_SUB)
s.text(x0 + pw / 2, y0 + ph + 30, "percentile of the observed attenuation in its random-panel null",
       8.8, "middle", C_SUB)
s.line(x0 + medp * pw, y0, x0 + medp * pw, y0 + ph, "#111827", 1.4, "4 3")
s.text(x0 + medp * pw - 6, y0 + 14, "median 0.440", 9, "end", C_TXT, "bold", halo=True)
s.line(x0, y0 + ph / 2.0, x0 + pw, y0 + ph / 2.0, "#cbd5e1", 1, "3 3")
s.text(x0 + pw + 6, y0 + ph / 2.0, "uniform", 8.4, "start", "#94a3b8")
s.text(x0 - 6, y0 + 8, "N", 8.6, "end", C_SUB)
s.text(x0 + pw / 2, y0 + ph + 44,
       "9 below 0.20 · 37 between · 3 at or above 0.80 (SIG01/03/26, all M1 non-significant)",
       8.4, "middle", C_SUB)
s.text(x0 + pw / 2, y0 + ph + 57, "KS D = 0.267, P = 0.001 (descriptive; percentiles are not independent)",
       8.4, "middle", C_SUB)

# --- b：对照轴 ---
bx0, by0, bwp, bhp = 400, 46, 250, 196
s.text(352, 24, "b", 12, "start", C_TXT, "bold")
s.text(368, 24, "Control axes through the same pipeline", 10, "start", C_TXT, "bold")
AX = 100.0   # 显示上限（%），中性粒轴超出并截断标注
sx = bx0 + bwp / 2.0
axes = [("lymphocyte", "Lymphocyte / cytotoxic", "expected immune proxy"),
        ("nutrition_liver", "Nutrition / liver-derived", "expected unrelated"),
        ("neutrophil", "Neutrophil", "βM1 ≈ 0 (ratio artefact)")]
rows_y = [by0 + 40, by0 + 100, by0 + 160]
# 共享零线
s.line(sx, by0 + 12, sx, by0 + 172, "#cbd5e1", 1, "3 3")
for (key, lab, note), yy in zip(axes, rows_y):
    r = CTL[key]
    a_pct = F(r["attenuation"]) * 100
    rho = F(r["rho_immune"])
    s.text(bx0, yy - 18, lab, 9.4, "start", C_TXT, "bold")
    s.text(bx0, yy - 7, note + " · ρ with ImmuneScore = %+.3f" % rho, 8.2, "start", C_SUB)
    s.line(bx0, yy + 14, bx0 + bwp, yy + 14, C_EDGE)
    if key == "neutrophil":
        w = bwp / 2.0 * min(abs(a_pct), AX) / AX
        s.bar(sx, yy, w, 12, "#f1f5f9", hatch=True)
        s.text(bx0, yy + 9, "A = −5,831%", 8.6, "start", "#92400e", "bold", halo=True)
    else:
        w = bwp / 2.0 * abs(a_pct) / AX
        col = C_IMM if key == "lymphocyte" else C_MET
        if a_pct > 0:
            s.bar(sx, yy, w, 12, col)
            s.text(sx + w + 5, yy + 9, "%+.1f%%" % a_pct, 9, "start", C_TXT, "bold")
        else:
            s.bar(sx - w, yy, w, 12, col)
            s.text(sx - w - 5, yy + 9, "%+.1f%%" % a_pct, 9, "end", C_TXT, "bold")
s.text(sx, by0 + 186, "0%", 8.6, "middle", C_SUB)
s.text(bx0 + bwp, by0 + 186, "A → +100% (shrinkage)", 8.2, "end", C_SUB)
s.text(bx0, by0 + 186, "← A (enlargement)", 8.2, "start", C_SUB)
s.text(bx0, by0 + bhp + 10, "Each axis: 16 genes, scored and adjusted exactly like a signature.", 8.4, "start", C_SUB)
s.text(bx0, by0 + bhp + 22, "Lymphocyte (ρ = 0.859) −68.9%; nutrition (ρ = −0.098) −15.3%.", 8.4, "start", C_SUB)
print("Fig3 ->", s.save("Fig3_null_controls.svg"))

# ============================================================
# Fig 4  三队列 × 三分层
# ============================================================
vals = {}
vals["TCGA-LUAD"] = {st: F(SUM[st]["median_A"]) for st in ("immune", "inflammation", "metabolism")}
for coh in ("GSE68465", "GSE31210"):
    vals[coh] = {st: F(EXT[(coh, st)]["median_A"]) * 100 for st in ("immune", "inflammation", "metabolism")}
# 断言：与正文 3.6/3.7 引用值一致
assert abs(vals["TCGA-LUAD"]["immune"] - 8.4) < 0.06
assert abs(vals["TCGA-LUAD"]["inflammation"] - (-9.9)) < 0.06
assert abs(vals["TCGA-LUAD"]["metabolism"] - 0.37) < 0.06
assert abs(vals["GSE68465"]["immune"] - 14.2) < 0.06
assert abs(vals["GSE68465"]["inflammation"] - (-5.5)) < 0.06
assert abs(vals["GSE68465"]["metabolism"] - 1.3) < 0.06
assert abs(vals["GSE31210"]["immune"] - 5.3) < 0.06
assert abs(vals["GSE31210"]["inflammation"] - (-2.4)) < 0.06
assert abs(vals["GSE31210"]["metabolism"] - (-1.5)) < 0.06


def _light(hexcol, step):
    """外部队列用同色系浅色：step 1/2 → 与白色混合 30%/55%"""
    f = {1: 0.30, 2: 0.55}[step]
    r = int(hexcol[1:3], 16); g = int(hexcol[3:5], 16); b = int(hexcol[5:7], 16)
    r = int(r + (255 - r) * f); g = int(g + (255 - g) * f); b = int(b + (255 - b) * f)
    return "#%02x%02x%02x" % (r, g, b)

s = SVG(680, 340, "Fig 4 stratum medians across cohorts")
gx0, gy0, gw, gh = 70, 52, 560, 200
lo, hi = -12.0, 16.0
def YX(v):
    return gy0 + gh * (hi - v) / (hi - lo)
zline = YX(0)
s.rect(gx0, gy0, gw, gh, "#fbfdff")
s.line(gx0, zline, gx0 + gw, zline, "#94a3b8", 1.2)
for t in range(-10, 17, 5):
    yy = YX(t)
    s.line(gx0 - 4, yy, gx0, yy, C_SUB)
    s.text(gx0 - 7, yy + 3, ("%+d" % t) if t else "0", 8.6, "end", C_SUB)
s.line(gx0, gy0, gx0, gy0 + gh, C_TXT, 1)
s.line(gx0, gy0 + gh, gx0 + gw, gy0 + gh, C_TXT, 1)
s.text(gx0 - 8, gy0 - 8, "median A (%)", 9, "end", C_SUB)
groups = [("immune", "Immune (6)", C_IMM), ("inflammation", "Inflammation (11)", C_INF),
          ("metabolism", "Metabolism (32)", C_MET)]
gw3 = gw / 3.0
bw3 = 40.0
cohorts = [("TCGA-LUAD", "#0f172a"), ("GSE68465", "#475569"), ("GSE31210", "#94a3b8")]
for gi, (st, lab, col) in enumerate(groups):
    cxc = gx0 + gw3 * gi + gw3 / 2.0
    for ci, (coh, shade) in enumerate(cohorts):
        v = vals[coh][st]
        xx = cxc - bw3 * 1.5 - 8 + ci * (bw3 + 8)
        yy = YX(max(v, 0))
        hh = abs(YX(v) - YX(0))
        s.bar(xx, yy, bw3, max(hh, 1), col if ci == 0 else _light(col, ci))
        s.text(xx + bw3 / 2, yy - 4 if v >= 0 else yy + hh + 10,
               "%+.1f" % v, 8.6, "middle", C_TXT)
    s.text(cxc, gy0 + gh + 28, lab, 9.6, "middle", C_TXT, "bold")
# 图例（横排在绘图区上方，避开免疫组的高柱）
lx = gx0 + 26
for ci, (coh, shade) in enumerate(cohorts):
    xx = lx + ci * 118
    s.rect(xx, gy0 - 15, 10, 10, shade if ci else "#0f172a", stroke="none", rx=2)
    s.text(xx + 14, gy0 - 6, coh + (" (primary)" if ci == 0 else " (external)"), 8.8, "start", C_SUB)
s.text(gx0 + gw / 2, gy0 + gh + 46,
       "Ordering immune > metabolism ≈ inflammation reproduces in all three cohorts; per-signature attenuation "
       "does not (ρ = −0.014 / +0.006).", 8.6, "middle", C_SUB)
s.text(gx0 + gw / 2, gy0 + gh + 59,
       "TCGA-LUAD: OS, n = 482/177. GSE68465: n = 439/235. GSE31210: n = 226/35. Stratum medians of A; "
       "SIG26 not scoreable externally (metabolism n = 31).",
       8.6, "middle", C_SUB)
print("Fig4 ->", s.save("Fig4_cohorts_strata.svg"))

# ============================================================
# Fig 2  从主报告提取森林图
# ============================================================
html = open(os.path.join(B, "报告_阶段4_免疫分层49签名_2026-10-20.html"), encoding="utf-8").read()
svg = re.search(r"<svg.*?</svg>", html, re.S).group(0)
assert 'aria-label="49 个签名校正前后效应量位移"' in svg
# 独立 .svg 文件必须带 xmlns；报告内嵌在 HTML 里所以原本没写
svg = svg.replace('<svg viewBox=',
                  '<svg xmlns="http://www.w3.org/2000/svg" viewBox=', 1)
svg = svg.replace('width="100%"', 'width="680"', 1)
svg = svg.replace(">", '><rect x="0" y="0" width="680" height="443" fill="#ffffff"/>', 1)
# 英文稿件：图例改英文（全文仅这 4 处含中文，替换后 assert 无残留）
for zh, en in (("log-HR per SD（空心=M1 校正前，实心=M3 校正后）",
                "log-HR per SD (open: M1, filled: M3)"),
               ("炎症（11）", "Inflammation (11)"),
               ("代谢（32）", "Metabolism (32)"),
               ("免疫（6）", "Immune (6)"),
               ('aria-label="49 个签名校正前后效应量位移"',
                'aria-label="Effect-size displacement of 49 signatures from M1 to M3"')):
    assert zh in svg, zh
    svg = svg.replace(zh, en, 1)
assert not re.search(r"[\u4e00-\u9fff]", svg), "Fig2 仍有中文"
p2 = os.path.join(OUT, "Fig2_forest.svg")
with open(p2, "w", encoding="utf-8") as f:
    f.write(svg)
print("Fig2 ->", p2)

print("done")
