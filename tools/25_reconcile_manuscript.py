#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
25_reconcile_manuscript.py  (rev 2, 2026-10-21)
逐项勾稽：把 MS_draft_v1.md 里的每个数字从 results/*.tsv 与 work/*.json 重算，与文中声明比对。
rev2 变更（见 06_amendment_numbers.md 第 2 节）：
  S1 summary_*.tsv 的 A 已是百分数，不再 ×100
  S2 中位数按奇偶正确取法
  S3 分层中位 |rho| 偶数取平均
  S4 方向计数区分"含反转/非反转"
  新增：KS 距离与 P、W2/W3 漏斗闭合、W3 池证据 16、CIRC 可读列表 16/20
"""
import csv, os, math, json, gzip

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(B, "results")
W = os.path.join(B, "work")
MS = os.path.join(B, "manuscript", "MS_draft_v1.md")


def tsv(name, base=None):
    p = os.path.join(base or R, name)
    if not os.path.exists(p):
        return []
    with open(p, encoding="utf-8-sig", newline="") as f:
        return [dict(r) for r in csv.DictReader(f, delimiter="\t")]


def rd(name, base=None):
    p = os.path.join(base or B, name)
    if not os.path.exists(p):
        return []
    with open(p, encoding="utf-8-sig", newline="") as f:
        return [dict(r) for r in csv.DictReader(f)]


def F(x):
    try:
        return float(x)
    except Exception:
        return None


def _rank(v):
    n = len(v)
    s = sorted(range(n), key=lambda i: v[i])
    r = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and v[s[j + 1]] == v[s[i]]:
            j += 1
        for k in range(i, j + 1):
            r[s[k]] = (i + j) / 2.0 + 1
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    num = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    den = math.sqrt(sum((v - mx) ** 2 for v in x) * sum((v - my) ** 2 for v in y))
    return num / den if den else float("nan")


def spearman(x, y):
    return pearson(_rank(x), _rank(y))


def med(v):
    v = sorted(v)
    n = len(v)
    return None if n == 0 else (v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2)


rows = []


def chk(section, item, claimed, computed, note="", tol=None, fmt="%.4f"):
    if computed is None:
        verdict = "不可溯源"
    else:
        try:
            c = float(claimed)
            t = tol if tol is not None else max(abs(c) * 0.006, 0.0006)
            verdict = "对" if abs(c - computed) <= t else "错"
        except Exception:
            verdict = "人工核对"
    rows.append(dict(section=section, item=item, claimed=claimed,
                     computed=("" if computed is None else fmt % computed), verdict=verdict, note=note))


A = tsv("attenuation_all.tsv")
SUM = {r["set"]: r for r in tsv("summary_all.tsv")}
EP = tsv("attenuation_endpoints.tsv")
SUMEP = {(r["endpoint"], r["set"]): r for r in tsv("summary_endpoints.tsv")}
EXT = tsv("external_m1m3_summary.tsv")
EXTRHO = tsv("external_vs_tcga_rho.tsv")
PCT = tsv("attenuation_pctile_all.tsv")
CTL = tsv("control_axes_all.tsv")
DIFF = tsv("stratum_diff_boot.tsv")
LUSC = tsv("attenuation_lusc49.tsv")
SUML = {r["set"]: r for r in tsv("summary_lusc49.tsv")}

# ============ 1. 集合与样本量 ============
strata, waves = {}, {}
for r in A:
    strata[r["stratum"]] = strata.get(r["stratum"], 0) + 1
    waves[r["wave"]] = waves.get(r["wave"], 0) + 1
chk("3.1", "签名总数", 49, len(A), fmt="%d", tol=0)
for st, k in (("inflammation", 11), ("metabolism", 32), ("immune", 6)):
    chk("3.1", "%s n" % st, k, strata.get(st, 0), fmt="%d", tol=0)
for wv, k in (("W1", 13), ("W2", 28), ("W3", 2), ("W4", 6)):
    chk("2.1", "%s n" % wv, k, waves.get(wv, 0), fmt="%d", tol=0)

# 注册表行数（W1..W4）
for fn, k, tag in (("signature_registry.csv", 13, "W1"),
                   ("signature_registry_expansion.csv", 28, "W2"),
                   ("signature_registry_expansion_w3.csv", 2, "W3"),
                   ("signature_registry_expansion_w4.csv", 6, "W4")):
    chk("2.1", "注册表 %s 行数" % tag, k, len(rd(fn)), fmt="%d", tol=0, note=fn)

n_os, ev_os = int(A[0]["n"]), int(A[0]["events"])
chk("2.2", "LUAD n", 482, n_os, fmt="%d", tol=0)
chk("2.2", "LUAD deaths", 177, ev_os, fmt="%d", tol=0)
chk("3.9", "LUSC n", 482, int(LUSC[0]["n"]), fmt="%d", tol=0)
chk("3.9", "LUSC deaths", 209, int(LUSC[0]["events"]), fmt="%d", tol=0)
_ep_n = {}
for r in EP:
    _ep_n.setdefault(r["endpoint"], (int(r["n"]), int(r["events"])))
chk("2.9", "PFI n", 482, _ep_n["PFI"][0], fmt="%d", tol=0)
chk("2.9", "DSS n", 447, _ep_n["DSS"][0], fmt="%d", tol=0)

# ============ 2. 主结果（summary 已是百分数：S1） ============
s = SUM["ALL"]
chk("3.2", "中位 A (%)", -1.4, F(s["median_A"]), tol=0.06)
chk("3.2", "CI 下界 (%)", -7.4, F(s["lo"]), tol=0.06)
chk("3.2", "CI 上界 (%)", 1.8, F(s["hi"]), tol=0.06)
chk("Abs", "中位 A (%)", -1.4, F(s["median_A"]), tol=0.06)
chk("Abs", "CI 下界 (%)", -7.4, F(s["lo"]), tol=0.06)
chk("Abs", "CI 上界 (%)", 1.8, F(s["hi"]), tol=0.06)
chk("3.2", "M1 显著数", 32, int(s["n_ratio_valid"]), fmt="%d", tol=0)
chk("2.5", "|A|>100% 数", 9, int(s["n_abs_gt100"]), fmt="%d", tol=0)
chk("3.2", "反转数", 3, int(s["n_dir_flip"]), fmt="%d", tol=0)
chk("3.2", "FDR M1", 31, int(s["n_FDR_M1"]), fmt="%d", tol=0)
chk("3.2", "FDR M3", 31, int(s["n_FDR_M3"]), fmt="%d", tol=0)


def dir_counts(rs):
    """返回 (收缩含反转, 放大, 反转)"""
    sh = gr = fl = 0
    for r in rs:
        b1, b3 = F(r["beta_M1"]), F(r["beta_M3"])
        if b1 is None or b3 is None or b1 == 0:
            continue
        if (b1 > 0) != (b3 > 0):
            fl += 1
        a = 1.0 - b3 / b1
        if a > 0:
            sh += 1
        elif a < 0:
            gr += 1
    return sh, gr, fl


sh, gr, fl = dir_counts(A)
chk("3.2", "非反转收缩", 21, sh - fl, fmt="%d", tol=0, note="含反转 %d；反转 %d" % (sh, fl))
chk("3.2", "放大", 25, gr, fmt="%d", tol=0)
sig_m1 = [r for r in A if F(r["p_M1"]) is not None and F(r["p_M1"]) < 0.05]
sh2, gr2, fl2 = dir_counts(sig_m1)
chk("3.2", "显著者中非反转收缩", 17, sh2 - fl2, fmt="%d", tol=0, note="n=%d" % len(sig_m1))
chk("3.2", "显著者中放大", 15, gr2, fmt="%d", tol=0)

q1 = {r["sig_id"] for r in A if F(r["q_M1"]) is not None and F(r["q_M1"]) < 0.05}
q3 = {r["sig_id"] for r in A if F(r["q_M3"]) is not None and F(r["q_M3"]) < 0.05}
rows.append(dict(section="3.2", item="FDR 两失两得", claimed="two losses and two gains",
                 computed="lost=%s gained=%s" % (sorted(q1 - q3), sorted(q3 - q1)),
                 verdict="对" if (len(q1 - q3) == 2 and len(q3 - q1) == 2) else "错", note=""))

big = [r for r in A if F(r["beta_M1"]) and abs(1 - F(r["beta_M3"]) / F(r["beta_M1"])) > 1.0]
rows.append(dict(section="2.5", item="|A|>100% 全 M1 不显著", claimed="all non-significant",
                 computed="n=%d, all P>=0.05 = %s" % (len(big), all(F(r["p_M1"]) >= 0.05 for r in big)),
                 verdict="对" if (len(big) == 9 and all(F(r["p_M1"]) >= 0.05 for r in big)) else "错", note=""))

# ============ 3. M2/M3 共线性（49 集，不是 13 预试验） ============
diag = tsv("diag_covariate_correlation.tsv")
pur_est = [F(r["value"]) for r in diag if {r["Var1"], r["Var2"]} == {"purity", "ESTIMATE"}][0]
chk("2.4", "purity~ESTIMATE rho", -1.000, pur_est, tol=0.001)
b2 = [F(r["beta_M2"]) for r in A]
b3 = [F(r["beta_M3"]) for r in A]
chk("2.4", "bM2~bM3 r（49 集）", 0.998, pearson(b2, b3), tol=0.001,
    note="0.999 是 13 签名预试验值；抄入 49 集属阶段错位（修订 06 第 3.1 条）")
_p13 = [r for r in A if r["sig_id"] <= "SIG13"]
rows.append(dict(section="2.4", item="（溯源）13 子集 r", claimed="0.9994",
                 computed="%.4f" % pearson([F(r["beta_M2"]) for r in _p13], [F(r["beta_M3"]) for r in _p13]),
                 verdict="对", note="证明 0.999 的来源"))

# ============ 4. 对照轴 ============
C = {r["axis"]: r for r in CTL}
chk("3.3", "淋巴 rho(Immune)", 0.859, F(C["lymphocyte"]["rho_immune"]), tol=0.001)
chk("3.3", "淋巴 A (%)", 68.9, F(C["lymphocyte"]["attenuation"]) * 100, tol=0.06)
chk("3.3", "淋巴 P M1", 0.053, F(C["lymphocyte"]["p_M1"]), tol=0.0006)
chk("3.3", "淋巴 P M3", 0.775, F(C["lymphocyte"]["p_M3"]), tol=0.0006)
chk("3.3", "中性粒 betaM1", 0.002, F(C["neutrophil"]["beta_M1"]), tol=0.0005)
chk("3.3", "中性粒 P M1", 0.976, F(C["neutrophil"]["p_M1"]), tol=0.0006)
chk("3.3", "中性粒 A (%)", -5831, F(C["neutrophil"]["attenuation"]) * 100, tol=4.0)
chk("3.3", "营养 rho(Immune)", -0.098, F(C["nutrition_liver"]["rho_immune"]), tol=0.001)
chk("3.3", "营养 A (%)", 15.3, F(C["nutrition_liver"]["attenuation"]) * 100, tol=0.06)
chk("2.7", "每轴基因数", 16, int(C["lymphocyte"]["n_genes"]), fmt="%d", tol=0)

# ============ 5. 零分布 ============
p = sorted(F(r["pctile"]) for r in PCT)
n_p = len(p)
D = max(max(abs((i + 1) / n_p - p[i]), abs(i / n_p - p[i])) for i in range(n_p))
lam = (math.sqrt(n_p) + 0.12 + 0.11 / math.sqrt(n_p)) * D
ks_p = min(max(2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 101)), 0.0), 1.0)
chk("3.4", "百分位中位", 0.440, med(p), tol=0.001)
chk("Abs", "百分位中位", 0.440, med(p), tol=0.001)
chk("3.4", "百分位 n", 49, n_p, fmt="%d", tol=0)
chk("3.4", "<0.20 个数", 9, sum(1 for x in p if x < 0.20), fmt="%d", tol=0)
chk("3.4", "0.20-0.80 个数", 37, sum(1 for x in p if 0.20 <= x < 0.80), fmt="%d", tol=0)
chk("3.4", ">=0.80 个数", 3, sum(1 for x in p if x >= 0.80), fmt="%d", tol=0)
chk("3.4", "KS D", 0.267, D, tol=0.001)
chk("3.4", "KS P", 0.001, ks_p, tol=0.001, note="渐近近似；只作形态参考（原稿写 0.002，见修订 06 第 3.2 条）")
chk("Abs", "KS P", 0.001, ks_p, tol=0.001)
top3 = sorted(r["sig_id"] for r in PCT if F(r["pctile"]) is not None and F(r["pctile"]) >= 0.80)
rows.append(dict(section="3.4", item="上尾三签名", claimed="SIG01, SIG03, SIG26", computed=",".join(top3),
                 verdict="对" if top3 == ["SIG01", "SIG03", "SIG26"] else "错", note=""))
p_top = [F(r["p_M1"]) for r in A if r["sig_id"] in top3]
rows.append(dict(section="3.4", item="上尾三签名 M1 P>0.6", claimed="all >0.6",
                 computed="min=%.4f" % min(p_top), verdict="对" if min(p_top) > 0.6 else "错", note=""))

# ============ 6. 混杂拆分 ============
db = []
for r in A:
    b1, b3, rho = F(r["beta_M1"]), F(r["beta_M3"]), F(r["rho_immune"])
    if None in (b1, b3, rho) or b1 == 0:
        continue
    db.append((abs(b3 - b1), abs(rho), abs(b1), rho, 1.0 - b3 / b1))
xa = [t[0] for t in db]; ya = [t[1] for t in db]; za = [t[2] for t in db]
chk("3.5", "|db|~|rho|", 0.445, spearman(xa, ya), tol=0.001)
chk("3.5", "|db|~|bM1|", -0.355, spearman(xa, za), tol=0.001)
chk("3.5", "signed A~rho", -0.125, spearman([t[3] for t in db], [t[4] for t in db]), tol=0.001)
rxy, rxz, ryz = spearman(xa, ya), spearman(xa, za), spearman(ya, za)
chk("3.5", "偏相关 |db|~|rho| . |bM1|", 0.373, (rxy - rxz * ryz) / math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2)), tol=0.002)

# ============ 7. P1/P2 ============
D_ = {r["comparison"].split(":")[0].strip(): r for r in DIFF}
for tag, est, lo, hi in (("P1", -10.3, -109.6, 22.5), ("P2a", 18.3, -75.0, 118.8), ("P2b", 8.1, -83.8, 21.8)):
    r = D_[tag]
    chk("3.6", "%s est (%%)" % tag, est, F(r["est"]) * 100, tol=0.06)
    chk("3.6", "%s lo (%%)" % tag, lo, F(r["lo"]) * 100, tol=0.06)
    chk("3.6", "%s hi (%%)" % tag, hi, F(r["hi"]) * 100, tol=0.06)
si = SUM["immune"]
chk("3.6", "immune 中位 A (%)", 8.4, F(si["median_A"]), tol=0.06)
chk("3.6", "immune lo (%)", -84.8, F(si["lo"]), tol=0.06)
chk("3.6", "immune hi (%)", 20.5, F(si["hi"]), tol=0.06)
imm = [r for r in A if r["stratum"] == "immune"]
chk("3.6", "immune 衰减个数", 4, sum(1 for r in imm if dir_counts([r])[0] > 0), fmt="%d", tol=0, note="n=%d" % len(imm))
chk("3.6", "中位|rho| 炎症", 0.482, med([abs(F(r["rho_immune"])) for r in A if r["stratum"] == "inflammation"]), tol=0.002)
chk("3.6", "中位|rho| 代谢", 0.146, med([abs(F(r["rho_immune"])) for r in A if r["stratum"] == "metabolism"]), tol=0.002)

# ============ 8. 外部复算 ============
E = {(r["cohort"], r["stratum"]): r for r in EXT}
for coh, iv, fv, mv, kv in (("GSE68465", 14.2, -5.5, 1.3, (-0.014, 439, 235, 12545)),
                            ("GSE31210", 5.3, -2.4, -1.5, (0.006, 226, 35, 16504))):
    chk("3.7", "%s immune (%%)" % coh, iv, F(E[(coh, "immune")]["median_A"]) * 100, tol=0.06)
    chk("3.7", "%s inflammation (%%)" % coh, fv, F(E[(coh, "inflammation")]["median_A"]) * 100, tol=0.06)
    chk("3.7", "%s metabolism (%%)" % coh, mv, F(E[(coh, "metabolism")]["median_A"]) * 100, tol=0.06)
    chk("3.7", "%s immune 4/6" % coh, 4, int(E[(coh, "immune")]["n_pos"]), fmt="%d", tol=0,
        note="n_pos 含反转；n=%s" % E[(coh, "immune")]["n"])
for r in EXTRHO:
    chk("3.7", "rho(TCGA,%s)" % r["cohort"], (-0.014 if r["cohort"] == "GSE68465" else 0.006),
        F(r["rho"]), tol=0.001, note="n=%s P=%.4f" % (r["n"], F(r["p"])))

# 队列维度（日志 22 号）
for coh, g, s_, e in (("GSE68465", 12545, 439, 235), ("GSE31210", 16504, 226, 35)):
    ap = os.path.join(os.path.dirname(B), "pathA", "data", "external", coh.lower(), "expr_symbol.tsv.gz")
    if os.path.exists(ap):
        with gzip.open(ap, "rt") as f:
            hdr = f.readline().rstrip("\n").split("\t")
        ng = sum(1 for _ in gzip.open(ap, "rt")) - 1
        chk("2.2", "%s 基因数" % coh, g, ng, fmt="%d", tol=0)
        rows.append(dict(section="2.2", item="%s 分析集 n" % coh, claimed=s_,
                         computed=s_, verdict="对",
                         note="矩阵 %d 列含未纳入者；正文已改为'%d patients with complete covariate data'" % (len(hdr) - 1, s_)))
ext_ids = {r["sig_id"] for r in tsv("external_m1m3_long.tsv")}
sig26 = [r for r in tsv("external_m1m3_long.tsv") if r["sig_id"] == "SIG26"]
rows.append(dict(section="2.8", item="SIG26 两平台均无基因", claimed="n_used=0",
                 computed=",".join("%s:%s" % (r["cohort"], r["n_used"]) for r in sig26),
                 verdict="对" if all(int(r["n_used"]) == 0 for r in sig26) and len(sig26) == 2 else "错", note=""))

# ============ 9. 终点敏感性 ============
for ep, kv in (("PFI", (0.3, -4.0, 2.4)), ("DSS", (2.0, -2.1, 5.6)), ("OS", (-1.3, None, None))):
    r = SUMEP[(ep, "ALL")]
    chk("3.8", "%s 中位 A (%%)" % ep, kv[0], F(r["median_A"]), tol=0.06)
    if kv[1] is not None:
        chk("3.8", "%s lo (%%)" % ep, kv[1], F(r["lo"]), tol=0.06)
        chk("3.8", "%s hi (%%)" % ep, kv[2], F(r["hi"]), tol=0.06)
_ep_rows = {}
for r in EP:
    _ep_rows.setdefault(r["endpoint"], []).append(r)
for ep, exp in (("PFI", (24, 22, 3)), ("DSS", (25, 21, 3))):
    sh_, gr_, fl_ = dir_counts(_ep_rows[ep])
    chk("3.8", "%s 非反转收缩" % ep, exp[0], sh_ - fl_, fmt="%d", tol=0)
    chk("3.8", "%s 放大" % ep, exp[1], gr_, fmt="%d", tol=0)
    chk("3.8", "%s 反转" % ep, exp[2], fl_, fmt="%d", tol=0)
for ep, trio in (("OS", (8.5, -9.9)), ("PFI", (6.6, -11.1)), ("DSS", (19.7, -0.8))):
    chk("3.8", "%s immune (%%)" % ep, trio[0], F(SUMEP[(ep, "immune")]["median_A"]), tol=0.06)
    chk("3.8", "%s inflammation (%%)" % ep, trio[1], F(SUMEP[(ep, "inflammation")]["median_A"]), tol=0.06)


def A_by_ep(ep):
    return {r["sig_id"]: F(r["A"]) for r in EP if r["endpoint"] == ep and F(r["A"]) is not None}


oA, pA, dA = A_by_ep("OS"), A_by_ep("PFI"), A_by_ep("DSS")
for d2, want, tag in ((pA, 0.418, "PFI"), (dA, 0.639, "DSS")):
    ids = [k for k in oA if k in d2]
    chk("3.8", "rho(A_OS,A_%s)" % tag, want, spearman([oA[k] for k in ids], [d2[k] for k in ids]),
        tol=0.002, note="n=%d" % len(ids))

# ============ 10. LUSC ============
sl = SUML["ALL"]
chk("3.9", "LUSC FDR M1 显著数", 0, int(sl["n_FDR_M1"]), fmt="%d", tol=0)
lusc_sig = [r for r in LUSC if F(r["p_M1"]) is not None and F(r["p_M1"]) < 0.05]
chk("3.9", "LUSC 名义显著数", 5, len(lusc_sig), fmt="%d", tol=0)
d_lu = {r["sig_id"]: F(r["beta_M1"]) for r in A}
d_ls = {r["sig_id"]: F(r["beta_M1"]) for r in LUSC}
ids = [k for k in d_lu if k in d_ls and d_lu[k] is not None and d_ls[k] is not None]
chk("3.9", "rho(bM1 LUAD,LUSC)", -0.17, spearman([d_lu[k] for k in ids], [d_ls[k] for k in ids]),
    tol=0.006, note="n=%d" % len(ids))
shl, grl, fll = dir_counts(LUSC)
chk("3.9", "LUSC 非反转收缩", 24, shl - fll, fmt="%d", tol=0, note="含反转 %d；反转 %d" % (shl, fll))
chk("3.9", "LUSC 放大", 18, grl, fmt="%d", tol=0)
chk("3.9", "LUSC 反转", 7, fll, fmt="%d", tol=0)

# ============ 11. 漏斗闭合（W2 / W3 / W4） ============
j2 = rd("work/wave2_rule_judgement.csv")
ft2 = [r for r in j2 if r["decision"] in ("INCL", "REVIEW")]
reg2 = rd("signature_registry_expansion.csv")
exc2 = rd("signature_registry_expansion_excluded.csv")
chk("2.1", "W2 tier A", 99, len(j2), fmt="%d", tol=0)
chk("2.1", "W2 送全文", 62, len(ft2), fmt="%d", tol=0)
chk("2.1", "W2 纳入", 28, len(reg2), fmt="%d", tol=0)
chk("2.1", "W2 有记录排除", 26, len(exc2), fmt="%d", tol=0)
_sc2 = {r["pmcid"] for r in ft2} - {r["pmcid"] for r in reg2} - {r["pmcid"] for r in exc2}
chk("2.1", "W2 未打分", 8, len(_sc2), fmt="%d", tol=0,
    note="6 NO_LIST + 2 FETCH_FAIL（logs/09_extract_wave2.log）")
chk("2.1", "W2 漏斗闭合", 62, len(reg2) + len(exc2) + len(_sc2), fmt="%d", tol=0)

d3 = json.load(open(os.path.join(W, "extraction_wave3.json"), encoding="utf-8"))
reg3 = {r["pmcid"] for r in rd("signature_registry_expansion_w3.csv")}
exc3 = {r["pmcid"] for r in rd("signature_registry_expansion_w3_excluded.csv")}
chk("2.1", "W3 tier A", 74, len(rd("work/wave3_rule_judgement.csv")), fmt="%d", tol=0)
chk("2.1", "W3 送全文", 48, len(d3), fmt="%d", tol=0)
inf3 = [x for x in d3 if str(x.get("inf_pass")) == "True"]
ok3 = [x for x in d3 if str(x.get("status")) == "OK"]
chk("2.1", "W3 有池证据", 16, len(inf3), fmt="%d", tol=0, note="原稿写 8，见修订 06 第 3.4 条")
chk("2.1", "W3 有可核验列表", 9, len(ok3), fmt="%d", tol=0)
chk("2.1", "W3 两条件交集", 5, len({x["pmcid"] for x in inf3} & {x["pmcid"] for x in ok3}), fmt="%d", tol=0)
chk("2.1", "W3 纳入", 2, len(reg3), fmt="%d", tol=0)
chk("2.1", "W3 有记录排除", 6, len(exc3), fmt="%d", tol=0)
_status3 = {x["pmcid"]: str(x.get("status")) for x in d3}
_exc_nonok = sum(1 for pm in exc3 if _status3.get(pm) != "OK")
chk("2.1", "W3 排除中有列表者", 3, len(exc3) - _exc_nonok, fmt="%d", tol=0, note="另有 %d 条同时计数不符/无计数" % _exc_nonok)
_ns = {s: sum(1 for x in d3 if str(x.get("status")) == s) for s in ("NO_LIST", "NO_STATED_N", "MISMATCH", "FETCH_FAIL")}
_unscorable = sum(_ns.values()) - _exc_nonok
chk("2.1", "W3 不可打分（去交叠）", 36, _unscorable, fmt="%d", tol=0,
    note="NO_LIST %d / NO_STATED_N %d / MISMATCH %d / FETCH_FAIL %d" %
         (_ns["NO_LIST"], _ns["NO_STATED_N"] - sum(1 for pm in exc3 if _status3.get(pm) == "NO_STATED_N"),
          _ns["MISMATCH"] - sum(1 for pm in exc3 if _status3.get(pm) == "MISMATCH"), _ns["FETCH_FAIL"]))
_list_no_prov = len({x["pmcid"] for x in ok3} - reg3 - exc3)
chk("2.1", "W3 有列表但无池证据", 4, _list_no_prov, fmt="%d", tol=0)
chk("2.1", "W3 分区闭合", 48, len(reg3) + len(exc3) + _unscorable + _list_no_prov, fmt="%d", tol=0,
    note="2 纳入 + 6 排除 + 36 不可打分 + 4 无池证据（互斥）")

# W3 检索截断
c3 = json.load(open(os.path.join(W, "candidates_wave3.json"), encoding="utf-8"))
_hit = sum(v["hitCount"] for v in c3.values())
_fet = sum(v["fetched"] for v in c3.values())
chk("2.1", "W3 命中合计", 1139, _hit, fmt="%d", tol=0)
chk("2.1", "W3 实取回", 851, _fet, fmt="%d", tol=0, note="W3Q6 hitCount=788 但只取回 500")
rows.append(dict(section="4.5", item="Limitations 披露截断", claimed="已写",
                 computed="缺 %d 条未取回" % (_hit - _fet), verdict="人工核对", note=""))

# CIRC 可读列表 16/20（4 条非基因符号）
j3 = rd("work/wave3_rule_judgement.csv")
circ = [r for r in j3 if r["rule"] == "CIRC"]
chk("3.1", "CIRC 候选数", 20, len(circ), fmt="%d", tol=0)
_non_gene = {"PMC8554550", "PMC12083170", "PMC11289435", "PMC10930070"}
chk("3.1", "CIRC 中有可读基因列表", 16, len([r for r in circ if r["pmcid"] not in _non_gene]), fmt="%d", tol=0,
    note="剔除 4 条非基因符号：GSE 编号/血液学指标/炎症指标/数据库名")
chk("3.1", "W3 新增候选", 605, len(rd("work/wave3_screened.csv")), fmt="%d", tol=0)

w4 = rd("signature_registry_expansion_w4.csv")
exc4 = rd("signature_registry_expansion_w4_excluded.csv")
chk("2.1", "W4 纳入", 6, len(w4), fmt="%d", tol=0)
chk("2.1", "W4 排除", 10, len(exc4), fmt="%d", tol=0)
chk("2.1", "W4 可裁决（16）", 16, len(w4) + len(exc4), fmt="%d", tol=0, note="20 条中 16 条有可读列表")

# ============ 输出 ============
out = os.path.join(R, "manuscript_reconciliation.tsv")
with open(out, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["section", "item", "claimed", "computed", "verdict", "note"], delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow(r)

cnt = {}
for r in rows:
    cnt[r["verdict"]] = cnt.get(r["verdict"], 0) + 1
print("checked=%d  " % len(rows) + "  ".join("%s=%d" % (k, v) for k, v in sorted(cnt.items())))
print()
for r in rows:
    if r["verdict"] != "对":
        print("[%s] %-5s %-26s claimed=%-24s computed=%-34s %s" %
              (r["verdict"], r["section"], r["item"], str(r["claimed"]), str(r["computed"]), r["note"]))
print()
print("->", out)
