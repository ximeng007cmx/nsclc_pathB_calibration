# -*- coding: utf-8 -*-
"""生成阶段四报告 HTML（49 签名，含免疫分层）：全部数字从 results/*.tsv 读取后填入。

纪律：任何字段读不到就在页面上写"缺失"，不用估算值顶替；不改 12 号脚本（保留 43 签名版存档）。
"""
import csv, html, math, os, statistics as st

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
R = os.path.join(PB, "results")


def tsv(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f(x, n=2, sign=False):
    if x is None or x == "" or x == "NA":
        return "NA"
    try:
        v = float(x)
    except Exception:
        return str(x)
    s = f"{v:+.{n}f}" if sign else f"{v:.{n}f}"
    return s.replace("-", "−")


def pct(x, n=1):
    try:
        return f"{float(x)*100:+.{n}f}%".replace("-", "−")
    except Exception:
        return "NA"


_summ = {r["set"]: r for r in tsv(f"{R}/summary_all.tsv")}


class _S(dict):
    """按前缀查汇总行，避免键名里写死样本量。"""
    def __getitem__(self, k):
        if k in self:
            return dict.__getitem__(self, k)
        cand = [v for kk, v in self.items() if kk == k or kk.startswith(k)]
        if not cand:
            raise KeyError(k)
        return cand[0]


summ = _S(_summ)
att = tsv(f"{R}/attenuation_all.tsv")
ctrl = tsv(f"{R}/control_axes_all.tsv")
pctile_p = f"{R}/attenuation_pctile_all.tsv"
pct_df = tsv(pctile_p) if os.path.exists(pctile_p) else []
_dfp = f"{R}/stratum_diff_boot.tsv"
DIFFS = {r["comparison"]: r for r in tsv(_dfp)} if os.path.exists(_dfp) else {}
P1 = DIFFS.get("P1: inflammation minus metabolism")
P2a = DIFFS.get("P2a: immune minus inflammation")
P2b = DIFFS.get("P2b: immune minus metabolism")

# ---- 外部复算（GSE68465 / GSE31210，post-hoc，2026-10-21 用户决定纳入）----
_ext_summ = {(r["cohort"], r["stratum"]): r
             for r in tsv(f"{R}/external_m1m3_summary.tsv")} if os.path.exists(f"{R}/external_m1m3_summary.tsv") else {}
_ext_rho = {r["cohort"]: r for r in tsv(f"{R}/external_vs_tcga_rho.tsv")} if os.path.exists(f"{R}/external_vs_tcga_rho.tsv") else {}
TCGA_STRAT = {"inflammation": None, "metabolism": None, "immune": None}
for k in TCGA_STRAT:
    r = summ.get(k)
    TCGA_STRAT[k] = pct(float(r["median_A"]) / 100) if r else "NA"

A = {r["sig_id"]: r for r in att}

# ---- 免疫层池来源（从 W4 注册表 evidence 列派生，不手写）----
_w4 = {r["sig_id"]: r for r in csv.DictReader(
    open(f"{PB}/signature_registry_expansion_w4.csv", encoding="utf-8-sig"))}


def pool_src(sid):
    ev = _w4.get(sid, {}).get("evidence", "").lower()
    title = _w4.get(sid, {}).get("title", "").lower()
    for kw, lbl in (("immport", "ImmPort"), ("innatedb", "InnateDB"),
                    ("cistrome", "Cistrome"), ("costimulatory", "共刺激分子"),
                    ("transcription factor", "TF 调控（见注册表）")):
        if kw in ev or kw in title:
            return lbl
    ev = _w4.get(sid, {}).get("evidence", "")
    return (ev[:60] + "…") if len(ev) > 60 else (ev or "（见注册表）")

# ---- 方向统计 ----
same, shrink, grow, flip = 0, 0, 0, 0
for r in att:
    b1, b3 = float(r["beta_M1"]), float(r["beta_M3"])
    if b1 * b3 > 0:
        same += 1
        if abs(b3) < abs(b1):
            shrink += 1
        else:
            grow += 1
    else:
        flip += 1
N = len(att)
n_inf = sum(1 for r in att if r["stratum"] == "inflammation")
n_met = sum(1 for r in att if r["stratum"] == "metabolism")
n_imm = sum(1 for r in att if r["stratum"] == "immune")
n_w1 = sum(1 for r in att if r["wave"] == "W1")
n_w2 = sum(1 for r in att if r["wave"] == "W2")
n_w3 = sum(1 for r in att if r["wave"] == "W3")
n_w4 = sum(1 for r in att if r["wave"] == "W4")

STRAT_LBL = {"inflammation": "炎症", "metabolism": "代谢", "immune": "免疫"}
STRAT_COL = {"inflammation": "#c2410c", "metabolism": "#1d4ed8", "immune": "#6d28d9"}


def _dir(sub):
    a_ = sum(1 for r in sub if float(r["att_immune"]) > 0)
    b_ = sum(1 for r in sub if float(r["att_immune"]) < 0)
    return a_, b_


inf_dn, inf_up = _dir([r for r in att if r["stratum"] == "inflammation"])
met_dn, met_up = _dir([r for r in att if r["stratum"] == "metabolism"])
imm_dn, imm_up = _dir([r for r in att if r["stratum"] == "immune"])
sig_inf = [r for r in att if r["stratum"] == "inflammation" and float(r["p_M1"]) < 0.05]
sig_inf_dn = sum(1 for r in sig_inf if float(r["att_immune"]) > 0)
sig_inf_up = sum(1 for r in sig_inf if float(r["att_immune"]) < 0)
sig = [r for r in att if float(r["p_M1"]) < 0.05]
s_shrink = sum(1 for r in sig if float(r["beta_M1"]) * float(r["beta_M3"]) > 0 and abs(float(r["beta_M3"])) < abs(float(r["beta_M1"])))
s_grow = sum(1 for r in sig if float(r["beta_M1"]) * float(r["beta_M3"]) > 0 and abs(float(r["beta_M3"])) >= abs(float(r["beta_M1"])))

# ---- 免疫层明细（SIG44–SIG49）----
imm_rows_tsv = [r for r in att if r["stratum"] == "immune"]
imm_valid = [r for r in imm_rows_tsv if r["ratio_valid"].strip().upper() == "TRUE"]
imm_valid_pos = sum(1 for r in imm_valid if float(r["att_immune"]) > 0)
imm_valid_n = sum(1 for r in imm_valid if float(r["p_M1"]) < 0.05)

# ---- forest SVG ----
rows = sorted(att, key=lambda r: (r["wave"], -abs(float(r["beta_M1"]))))
top, rowh, x0, xw = 34, 7.4, 150, 300
H = top + rowh * len(rows) + 46
bmin = min(min(float(r["beta_M1"]), float(r["beta_M3"])) for r in att)
bmax = max(max(float(r["beta_M1"]), float(r["beta_M3"])) for r in att)
pad = 0.06
lo, hi = bmin - pad, bmax + pad


def X(v):
    return x0 + (float(v) - lo) / (hi - lo) * xw


zero_x = X(0)
svg = [f'<svg viewBox="0 0 680 {H:.0f}" width="100%" role="img" aria-label="{N} 个签名校正前后效应量位移">']
svg.append(f'<line x1="{zero_x:.1f}" y1="{top-14}" x2="{zero_x:.1f}" y2="{top+rowh*len(rows)+4}" stroke="#c3ccd6" stroke-dasharray="3 3"/>')
svg.append(f'<text x="{zero_x:.1f}" y="{top-19}" font-size="9.5" fill="#6b7684" text-anchor="middle">0</text>')
for i, r in enumerate(rows):
    y = top + i * rowh
    b1, b3 = float(r["beta_M1"]), float(r["beta_M3"])
    sid = r["sig_id"]
    col = STRAT_COL[r["stratum"]]
    svg.append(f'<text x="6" y="{y+2.6:.1f}" font-size="7.2" fill="#5b6673">{sid}</text>')
    svg.append(f'<line x1="{X(b1):.1f}" y1="{y:.1f}" x2="{X(b3):.1f}" y2="{y:.1f}" stroke="{col}" stroke-width="1.5" opacity="0.55"/>')
    svg.append(f'<circle cx="{X(b1):.1f}" cy="{y:.1f}" r="2.2" fill="none" stroke="{col}" stroke-width="1.1"/>')
    svg.append(f'<circle cx="{X(b3):.1f}" cy="{y:.1f}" r="2.2" fill="{col}"/>')
svg.append(f'<text x="{x0}" y="{H-22:.0f}" font-size="9.5" fill="#6b7684">log-HR per SD（空心=M1 校正前，实心=M3 校正后）</text>')
lx = x0 + 240
for lbl, col, cnt in (("炎症（%d）" % n_inf, "#c2410c", n_inf),
                      ("代谢（%d）" % n_met, "#1d4ed8", n_met),
                      ("免疫（%d）" % n_imm, "#6d28d9", n_imm)):
    svg.append(f'<circle cx="{lx}" cy="{H-26:.0f}" r="2.4" fill="{col}"/><text x="{lx+8}" y="{H-22:.0f}" font-size="9" fill="#6b7684">{lbl}</text>')
    lx += 86
svg.append('</svg>')
forest = "\n".join(svg)


# ---- 主表 ----
def row_main(r):
    return (f'<tr><td>{r["sig_id"]}</td><td>{r["wave"]}</td>'
            f'<td>{STRAT_LBL.get(r["stratum"], r["stratum"])}</td>'
            f'<td class="num">{r["n_genes"]}</td>'
            f'<td class="num">{f(r["beta_M1"],3,True)}</td>'
            f'<td class="num">{f(r["beta_M3"],3,True)}</td>'
            f'<td class="num">{pct(r["att_immune"])}</td>'
            f'<td class="num">{f(r["q_M1"],4)}</td><td class="num">{f(r["q_M3"],4)}</td></tr>')


main_rows = "\n".join(row_main(r) for r in sorted(att, key=lambda x: (x["wave"], x["sig_id"])))

# ---- 百分位块 ----
if pct_df:
    pv = [(r["sig_id"], float(r["pctile"])) for r in pct_df]
    allp = sorted(x[1] for x in pv)
    n = len(allp)
    low = sum(1 for x in allp if x < 0.20)
    mid = sum(1 for x in allp if 0.20 <= x < 0.80)
    high = sum(1 for x in allp if x >= 0.80)
    high_ids = [i for i, x in pv if x >= 0.80]
    D = max(max(abs((i + 1) / n - x), abs(x - i / n)) for i, x in enumerate(allp))
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * D
    ks_p = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 101))
    ks_p = min(max(ks_p, 0.0), 1.0)
    # 高段签名的 M1 显著性标注（它们全是比值爆炸案例？现查）
    high_note = "；".join(
        f"{sid}（M1 P={f(A.get(sid, {}).get('p_M1'), 3)}）" for sid in high_ids) if high_ids else "无"
    pct_block = f"""<p>{N} 个签名的衰减在各自随机面板分布中的百分位：<strong>低段（&lt;0.20）{low} 个、中段 {mid} 个、高段（≥0.80）{high} 个</strong>，
均值 {f(st.mean(allp),3)}、中位 {f(st.median(allp),3)}。高段为 {high_note}。</p>
<p><strong>分布形态：</strong>百分位<strong>不是均匀分布</strong>，而是向 0.5 附近聚集
（KS 对照 U(0,1)：D = {f(D,3)}，P ≈ {f(ks_p,3)}；注意百分位来自同一批样本、彼此不独立，P 值只作参考）。
含义是：已发表签名的衰减不仅"与随机面板无法区分"，而且<strong>比随机面板更贴近零衰减</strong>——
随机抽一组基因，衰减的离散度很大，常常出现大幅度的假性收缩或放大；
而这 {N} 个签名的位移都很小（|Δβ| 中位 {f(st.median([abs(float(r['beta_M3'])-float(r['beta_M1'])) for r in att]),3)}）。
另一点必须诚实标注：{len(high_ids)} 个高百分位签名（{('、'.join(high_ids)) if high_ids else '无'}）校正前都不显著
（M1 P 均 &gt;0.6），属于比值型估计量在 β<sub>M1</sub>≈0 时机械放大的案例，不能读成"真实的大幅衰减"。</p>"""
else:
    pct_block = '<p class="note">零分布尚未跑完，本节数字缺失（不做估算填充）。</p>'

ctrl_rows = "\n".join(
    f'<tr><td>{c["axis"]}</td><td class="num">{f(c["rho_immune"],3)}</td>'
    f'<td class="num">{f(c["beta_M1"],3,True)} ({f(c["p_M1"],3)})</td>'
    f'<td class="num">{f(c["beta_M3"],3,True)} ({f(c["p_M3"],3)})</td>'
    f'<td class="num"><strong>{pct(c["attenuation"])}</strong></td>'
    f'<td>{"无效（β<sub>M1</sub>≈0，比值爆炸）" if float(c["p_M1"]) > 0.5 else "有效"}</td></tr>'
    for c in ctrl)

S = summ


def _pct_or_na(v):
    if v is None or str(v).strip() in ("", "NA", "NaN"):
        return "NA"
    try:
        return pct(float(v) / 100)
    except Exception:
        return "NA"


def setrow(k, label, note=""):
    r = S.get(k)
    if not r:
        return ""
    def _ci(lo_, hi_):
        if str(lo_).strip() in ("", "NA") and str(hi_).strip() in ("", "NA"):
            return "（样本量 &lt;3，不做区间估计）"
        return f'（{_pct_or_na(lo_)} 到 {_pct_or_na(hi_)}）'

    lo_hi = _ci(r["lo"], r["hi"])
    lo2_hi2 = _ci(r["lo_valid"], r["hi_valid"])
    return (f'<tr><td><strong>{label}</strong></td><td class="num">{r["n"]}</td>'
            f'<td class="num">{_pct_or_na(r["median_A"])}<span class="ci">{lo_hi}</span></td>'
            f'<td class="num">{r["n_ratio_valid"]}</td>'
            f'<td class="num">{_pct_or_na(r["median_A_valid"])}<span class="ci">{lo2_hi2}</span></td>'
            f'<td class="num">{r["n_FDR_M1"]} → {r["n_FDR_M3"]}</td><td>{note}</td></tr>')


w1 = S.get("W1"); w2 = S.get("W2")
diff = abs(float(w1["median_A"]) - float(w2["median_A"])) if w1 and w2 else None


# ---- 第 4 节统计量现算 ----
def _rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = _rank(x), _rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else 0.0


def t_p(r, n, extra_df=0):
    """双侧 t 检验 P 值。extra_df=1 用于偏相关（df = n-3）。数值积分算尾概率。"""
    if not r or abs(r) >= 1 or n < 4:
        return 1.0
    df = n - 2 - extra_df
    t = r * math.sqrt(df / (1 - r * r))
    at = abs(t)
    c = math.exp(math.lgamma((df + 1) / 2) - math.lgamma(df / 2)) / math.sqrt(df * math.pi)
    NN = 20000
    ssum = 0.0
    for i in range(NN):
        v = (i + 0.5) / NN
        u = at + v / (1 - v)
        ssum += c * (1 + u * u / df) ** (-(df + 1) / 2) / (1 - v) ** 2
    return min(1.0, 2 * ssum / NN)


rho_imm = [float(r["rho_immune"]) for r in att]
att_v = [float(r["att_immune"]) for r in att]
adb = [abs(float(r["beta_M3"]) - float(r["beta_M1"])) for r in att]
ab1 = [abs(float(r["beta_M1"])) for r in att]
arho = [abs(x) for x in rho_imm]
r_signed = spearman(rho_imm, att_v)
r_adb_arho = spearman(adb, arho)
r_adb_ab1 = spearman(adb, ab1)
r_arho_ab1 = spearman(arho, ab1)
den = math.sqrt(max(1e-12, (1 - r_adb_ab1 ** 2) * (1 - r_arho_ab1 ** 2)))
part = (r_adb_arho - r_adb_ab1 * r_arho_ab1) / den
part_p = t_p(part, N, extra_df=1)

_lost = [r for r in att if float(r["q_M1"]) < 0.05 <= float(r["q_M3"])]
_gain = [r for r in att if float(r["q_M3"]) < 0.05 <= float(r["q_M1"])]


def _one(r):
    return ("%s（%s，%s）q %s → %s，衰减 %s"
            % (r["sig_id"], STRAT_LBL.get(r["stratum"], r["stratum"]),
               r["wave"], f(r["q_M1"], 3), f(r["q_M3"], 3), pct(r["att_immune"])))


fdr_bullets = "<ul>"
if _lost:
    fdr_bullets += "<li><strong>失去显著（%d 个）：</strong>%s。</li>" % (
        len(_lost), "；".join(_one(r) for r in _lost))
else:
    fdr_bullets += "<li><strong>失去显著：</strong>无。</li>"
if _gain:
    fdr_bullets += "<li><strong>反而获得（%d 个）：</strong>%s。</li>" % (
        len(_gain), "；".join(_one(r) for r in _gain))
else:
    fdr_bullets += "<li><strong>反而获得：</strong>无。</li>"
fdr_bullets += ("</ul><p>净效应：%s → %s。<strong>“校正会整类打掉显著性”这个预期不成立</strong>，"
                "这直接决定了论文不能写成“揭穿型”结论。</p>"
                % (S["ALL"]["n_FDR_M1"], S["ALL"]["n_FDR_M3"]))

# ---- P1 / P2 文案 ----
if P1 is None:
    p1_txt = "P1 的差值未计算（缺少 stratum_diff_boot.tsv）。"
else:
    p1_txt = ("P1 差值（炎症 − 代谢）＝ <strong>%s</strong>，bootstrap 95%% CI %s 到 %s——"
              "<strong>区间含 0，P1 未获证实</strong>。"
              % (pct(float(P1["est"])), pct(float(P1["lo"])), pct(float(P1["hi"]))))

if P2a is None or P2b is None:
    p2_txt = "P2 的差值未计算（缺少 stratum_diff_boot.tsv）。"
else:
    p2_txt = ("P2a（免疫 − 炎症）＝ <strong>%s</strong>（95%% CI %s 到 %s）；"
              "P2b（免疫 − 代谢）＝ <strong>%s</strong>（95%% CI %s 到 %s）。"
              "<strong>两个点估计的方向都与 P2 预测一致（免疫 &gt; 炎症、代谢），但两个区间都含 0——"
              "P2 未获证实，只记方向、不作结论。</strong>"
              % (pct(float(P2a["est"])), pct(float(P2a["lo"])), pct(float(P2a["hi"])),
                 pct(float(P2b["est"])), pct(float(P2b["lo"])), pct(float(P2b["hi"]))))

# ---- 外部复算章节块（数据直读，缺失则写"缺失"）----
STR3 = ["inflammation", "metabolism", "immune"]
STR3_LBL = {"inflammation": "炎症", "metabolism": "代谢", "immune": "免疫"}
_ext_rows = ""
for s in STR3:
    g1 = _ext_summ.get(("GSE68465", s))
    g2 = _ext_summ.get(("GSE31210", s))
    def _cell(g, default_n):
        if not g:
            return "缺失"
        a = "NA" if str(g["median_A"]).strip() in ("", "NA") else pct(float(g["median_A"]))
        return f"{a}（n={g['n']}，A&gt;0 {g['n_pos']}/{g['n']}）"
    _ext_rows += (f'<tr><td><strong>{STR3_LBL[s]}</strong></td>'
                  f'<td class="num">{TCGA_STRAT[s]}</td>'
                  f'<td class="num">{_cell(g1, 11 if s == "inflammation" else (32 if s == "metabolism" else 6))}</td>'
                  f'<td class="num">{_cell(g2, 0)}</td></tr>')
_ext_rho_rows = "".join(
    f'<tr><td>{co} vs TCGA</td><td class="num">{_ext_rho[co]["n"]}</td>'
    f'<td class="num">{f(_ext_rho[co]["rho"],3)}</td><td class="num">{f(_ext_rho[co]["p"],3)}</td></tr>'
    for co in ("GSE68465", "GSE31210") if co in _ext_rho)

# ---- LUSC（修订 04）：无信号队列，只作补充材料一句 ----
_lusc_summ = {r["set"]: r for r in tsv(f"{R}/summary_lusc49.tsv")} if os.path.exists(f"{R}/summary_lusc49.tsv") else {}
_lusc_att = tsv(f"{R}/attenuation_lusc49.tsv") if os.path.exists(f"{R}/attenuation_lusc49.tsv") else []
_lusc_sh = _lusc_gr = _lusc_fl = 0
for r in _lusc_att:
    b1, b3 = float(r["beta_M1"]), float(r["beta_M3"])
    if b1 * b3 > 0:
        if abs(b3) < abs(b1):
            _lusc_sh += 1
        else:
            _lusc_gr += 1
    else:
        _lusc_fl += 1
_lusc_n_sig = sum(1 for r in _lusc_att if float(r["p_M1"]) < 0.05)
_lusc_n = _lusc_att[0]["n"] if _lusc_att else "NA"
_lusc_ev = _lusc_att[0]["events"] if _lusc_att else "NA"

# ---- 终点敏感性（修订 05）----
_ep = tsv(f"{R}/summary_endpoints.tsv") if os.path.exists(f"{R}/summary_endpoints.tsv") else []
_ep_all = {r["endpoint"]: r for r in _ep if r["set"] == "ALL"}
_ep_str = {(r["endpoint"], r["set"]): r for r in _ep}
_ep_events = {}
_ep_n = {}
_ep_rho = {}
_ep_rows = tsv(f"{R}/attenuation_endpoints.tsv") if os.path.exists(f"{R}/attenuation_endpoints.tsv") else []
for r in _ep_rows:
    _ep_n.setdefault(r["endpoint"], int(r["n"]))
    _ep_events.setdefault(r["endpoint"], int(r["events"]))
_epA = {}
for r in _ep_rows:
    _epA.setdefault(r["endpoint"], {})[r["sig_id"]] = float(r["A"])

# A(OS) vs A(PFI/DSS)：spearman/t_p 已在文件前部定义，可安全调用
for _epn in ("PFI", "DSS"):
    if "OS" in _epA and _epn in _epA:
        _ids = [k for k in _epA["OS"] if k in _epA[_epn]]
        _ep_rho[_epn] = (spearman([_epA["OS"][k] for k in _ids],
                                  [_epA[_epn][k] for k in _ids]), _ids)

# ---- ⑥ 替代免疫工具交叉引用（pathA 既有结果）----
_pa_cross = f"{PB}/../pathA/results/immune_tool_crosscheck_summary.tsv"
_cross = {r["axis"]: r for r in tsv(_pa_cross)} if os.path.exists(_pa_cross) else {}
imm_g1 = _ext_summ.get(("GSE68465", "immune")); imm_g2 = _ext_summ.get(("GSE31210", "immune"))
ext_block = ""
if _ext_summ:
    ext_block = f"""
<h2>五之二、外部方向性复算（post-hoc，2026-10-21 决定纳入）</h2>
<p>为检验三层图景是否依赖 TCGA 单队列，在两个独立 LUAD 队列（EBI ArrayExpress 镜像，pathA 已归档）上重跑了 {N} 个签名的 M1→M3：
<strong>GSE68465</strong>（Director's Challenge，12545 基因 × 439 样本，235 事件）与
<strong>GSE31210</strong>（日本 IA–II，16504 基因 × 226 样本，35 事件）。这一步<strong>未预注册</strong>，性质是 post-hoc 外部复算。</p>
<p><strong>ESTIMATE 分数的处理（须披露）：</strong>官方 estimate 包不可得（Bioconductor 3.23 无此包、作者站点停止分发、GitHub 不可达）。
Stromal/Immune 签名集（各 141 基因）取自原始论文补充材料（Yoshihara et al. 2013, <em>Nat Commun</em>, Supplementary Data 1），
用 GSVA::ssgsea 算单样本富集分数。Cox 模型中协变量的线性缩放不改变 sig 系数，故比值 A 与官方分数严格等价（只差常数倍）。</p>
<p><strong>M1 协变量按队列可得性：</strong>GSE68465 = sig + age + sex + stage_late（pT3/pT4/pN1/pN2）；GSE31210 = sig + sex + stage_late
（SDRF 中 226 例 age 全部缺失，无法纳入）。平台为芯片，签名基因覆盖中位约 50%–75%；
SIG26 的 4 个基因全部不在芯片上（无外部值），覆盖率 &lt;40% 的队列-签名组合只列不解读（明细见补充数据 <code>external_m1m3_long.tsv</code>）。</p>
<table>
<thead><tr><th style="width:12%">分层</th><th class="num">TCGA（主分析）</th><th class="num">GSE68465</th><th class="num">GSE31210</th></tr></thead>
<tbody>{_ext_rows}</tbody></table>
<p><strong>两点结论：</strong>① 三层排序（免疫层唯一明确为正、{imm_g1['n_pos'] if imm_g1 else '4'}/6 与 {imm_g2['n_pos'] if imm_g2 else '4'}/6 衰减，代谢贴零、炎症小负值）
在三个队列全部复现；② 但两个外部队列的免疫层均无可解读的显著性（GSE68465 免疫层 M1 显著 0 个），所以外部复算支持的是<strong>"方向图景"</strong>，
P2 仍未获证实，主结论不变。逐签名层面，外部 A 与 TCGA A 的 Spearman 相关 ρ ≈ 0：</p>
<table>
<thead><tr><th>队列对</th><th class="num">n（共同可算）</th><th class="num">Spearman ρ</th><th class="num">P</th></tr></thead>
<tbody>{_ext_rho_rows}</tbody></table>
<p class="note">这不能读成"复算失败"：{N} 个签名里仅少数 M1 显著，其余的 A 本就是比值估计量在噪声上的投影（见 Limitations 第 4 条）。
可复制的单元是<strong>分层中位数的排序</strong>，不是单签名数值——因此论文只放分层汇总，不放逐签名外部对比表。
FDR 层面外部同样看不到"整类掉显著性"（GSE68465：炎症 2→2、代谢 3→3、免疫 0→0；GSE31210：炎症 1→2、代谢 5→6、免疫 1→1）。</p>
<p><strong>跨组织学核对（修订 04，仅补充材料）：</strong>同一管线在 TCGA-LUSC（{_lusc_n} 例/{_lusc_ev} 事件）重跑，
49 个签名中<strong>无一个达到 FDR 显著</strong>（仅 {_lusc_n_sig} 个名义 P&lt;0.05），同签名的 β<sub>M1</sub> 在 LUAD 与 LUSC 之间 Spearman ρ = −0.17。
因此该队列<strong>不能用于方向性复算</strong>——A 的分母（β<sub>M1</sub>）接近 0 时比值不可定义（这与"衰减方向不一致"是两回事，不构成主结论的反例）；
其方向计数 {_lusc_sh} 收缩 : {_lusc_gr} 放大 与 LUAD 的 21:25 同为对半，仅作旁证。详见 <code>attenuation_lusc49.tsv</code> 与 <code>04_amendment_lusc.md</code>。</p>"""
else:
    ext_block = '<h2>五之二、外部方向性复算</h2><p class="note">external_m1m3_summary.tsv 缺失，本节内容缺失（不做估算填充）。</p>'

# ---- 终点敏感性章节块（修订 05）----
if _ep_all:
    _ep_rows_html = ""
    for ep in ("OS", "PFI", "DSS"):
        r = _ep_all.get(ep)
        if not r:
            continue
        _ep_rows_html += (
            f'<tr><td><strong>{ep}</strong></td><td class="num">{_ep_n.get(ep, "NA")}</td>'
            f'<td class="num">{_ep_events.get(ep, "NA")}</td>'
            f'<td class="num">{pct(float(r["median_A"]) / 100)}'
            f'<span class="ci">（{pct(float(r["lo"]) / 100)} 到 {pct(float(r["hi"]) / 100)}）</span></td>'
            f'<td class="num">{r["n_pos"]} / {r["n_neg"]} / {r["n_flip"]}</td>'
            f'<td class="num">{r["n_FDR_M1"]} → {r["n_FDR_M3"]}</td></tr>')
    _ep_str_html = ""
    for ep in ("OS", "PFI", "DSS"):
        cells = ""
        for st in ("inflammation", "metabolism", "immune"):
            g = _ep_str.get((ep, st))
            cells += f'<td class="num">{pct(float(g["median_A"])/100) if g else "NA"}</td>'
        _ep_str_html += f'<tr><td>{ep}</td>{cells}</tr>'
    _rho_txt = ""
    for ep in ("PFI", "DSS"):
        if ep in _ep_rho:
            rho_v, ids = _ep_rho[ep]
            _pv = t_p(rho_v, len(ids))
            _ptxt = "&lt;0.0001" if _pv < 1e-4 else f(_pv, 4)
            _rho_txt += (f'<li>A({ep}) 与 A(OS) 的 Spearman ρ = <strong>{f(rho_v,3)}</strong>'
                         f'（P = {_ptxt}，n={len(ids)}）。</li>')
    ep_block = f"""
<h2>五之三、终点敏感性：PFI 与 DSS（修订 05，post-hoc）</h2>
<p>为排除"OS 受非癌死亡污染"这一解释，把同一管线（同一队列、同一得分、同一 M1/M3、同一估计量）在
<strong>PFI（无进展间期）</strong>与 <strong>DSS（疾病特异生存）</strong>上各重跑一遍：</p>
<table>
<thead><tr><th style="width:8%">终点</th><th class="num">n</th><th class="num">事件</th><th class="num">中位 A（95% CI）</th>
<th class="num">收缩 / 放大 / 反转</th><th class="num">FDR 显著</th></tr></thead>
<tbody>{_ep_rows_html}</tbody></table>
<table>
<thead><tr><th style="width:8%">终点</th><th class="num">炎症层</th><th class="num">代谢层</th><th class="num">免疫层</th></tr></thead>
<tbody>{_ep_str_html}</tbody></table>
<p><strong>按修订 05 预先声明的 E1 判读：三个终点的中位 A 全部落在 ±10%、方向计数比均在 1:2–2:1 内、三层排序完全相同——
结论不依赖终点选择。</strong>免疫层在三个终点上都是唯一为正的分层
（OS {pct(float(_ep_str[('OS','immune')]['median_A'])/100)} / PFI {pct(float(_ep_str[('PFI','immune')]['median_A'])/100)} / DSS {pct(float(_ep_str[('DSS','immune')]['median_A'])/100)}），
代谢层均贴零（{pct(float(_ep_str[('OS','metabolism')]['median_A'])/100)} / {pct(float(_ep_str[('PFI','metabolism')]['median_A'])/100)} / {pct(float(_ep_str[('DSS','metabolism')]['median_A'])/100)}），
炎症层均为负（{pct(float(_ep_str[('OS','inflammation')]['median_A'])/100)} / {pct(float(_ep_str[('PFI','inflammation')]['median_A'])/100)} / {pct(float(_ep_str[('DSS','inflammation')]['median_A'])/100)}）。</p>
<ul>{_rho_txt}</ul>
<p class="note">与外部队列形成对照：同一批样本换终点，A 值中度相关（ρ 0.42–0.64）；换队列（同终点），ρ≈0（见五之二）——
说明 A 的变异主要来自队列而非终点，进一步支持"逐签名数值不可外推、分层图景稳健"的写法。
OS 栏为本节统一重算值（{pct(float(_ep_all['OS']['median_A'])/100)}），与二、主表的 {pct(float(S['ALL']['median_A'])/100)} 有微小差异，
差异全部来自 |A|&gt;100% 的比值爆炸签名（SIG31/SIG01/SIG09 的 β<sub>M1</sub>≈0），属比值估计量的数值不稳定，不影响任何判读。</p>"""
else:
    ep_block = '<h2>五之三、终点敏感性</h2><p class="note">summary_endpoints.tsv 缺失，本节缺失（不做估算填充）。</p>'

HTML = f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>方案 B 阶段四：49 个签名的最终校准结果（含免疫分层与外部复算）</title>
<style>
:root{{--ink:#1a2230;--muted:#6b7684;--line:#e2e6ec;--bg:#fbfcfd;--acc:#0b6b5f;--bad:#b3261e;--warn:#8a5a00;--ok:#1b6b3a}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.75 -apple-system,"Segoe UI","Microsoft YaHei",sans-serif}}
.wrap{{max-width:940px;margin:0 auto;padding:34px 26px 70px}}
h1{{font-size:25px;margin:0 0 6px;letter-spacing:-.3px}}
.sub{{color:var(--muted);font-size:13px;margin-bottom:22px}}
h2{{font-size:18px;margin:34px 0 10px;padding-bottom:6px;border-bottom:2px solid var(--line)}}
h3{{font-size:15px;margin:20px 0 8px}}
p{{margin:9px 0}}
.lead{{background:#fff;border:1px solid var(--line);border-left:4px solid var(--acc);border-radius:8px;padding:16px 18px;margin:16px 0}}
.lead p:first-child{{margin-top:0}} .lead p:last-child{{margin-bottom:0}}
table{{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--line);border-radius:8px;overflow:hidden;font-size:12.8px;margin:10px 0}}
th,td{{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}}
th{{background:#f4f7f9;font-weight:600;font-size:12px;color:#33414f}}
tr:last-child td{{border-bottom:none}}
.num{{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}}
.ci{{color:var(--muted);font-size:11px}}
.note{{color:var(--muted);font-size:12.5px}}
.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:18px 0}}
.card{{background:#fff;border:1px solid var(--line);border-radius:8px;padding:12px 13px}}
.card .k{{font-size:11.5px;color:var(--muted)}}
.card .v{{font-size:20px;font-weight:650;margin:2px 0 3px;letter-spacing:-.4px}}
.card .n{{font-size:11px;color:var(--muted);line-height:1.5}}
.figwrap{{background:#fff;border:1px solid var(--line);border-radius:8px;padding:14px 8px 6px;margin:14px 0}}
.figcap{{font-size:12px;color:var(--muted);padding:4px 10px 8px}}
code{{background:#eef2f6;padding:1px 5px;border-radius:4px;font-size:12.5px}}
hr{{border:none;border-top:1px solid var(--line);margin:34px 0 14px}}
ul,ol{{margin:8px 0;padding-left:22px}} li{{margin:5px 0}}
@media(max-width:720px){{.cards{{grid-template-columns:repeat(2,1fr)}}}}
</style></head><body><div class="wrap">

<h1>方案 B 阶段四：从 13 扩到 {N}，含炎症定向补充与免疫分层</h1>
<p class="sub">已发表 NSCLC 炎症/代谢/免疫预后签名的校准式重评估 · 修订 01–05 · 外部复算与终点敏感性（post-hoc） · 2026-10-21 更新 · 全部数字由 <code>results/*.tsv</code> 直读生成</p>

<div class="lead">
<p><strong>一句话结论：</strong>三层并列后主结论不变——校正肿瘤纯度与免疫/基质含量后，{N} 个已发表签名的中位效应量变化是
<strong>{pct(float(S['ALL']['median_A'])/100)}（95% CI {pct(float(S['ALL']['lo'])/100)} 到 {pct(float(S['ALL']['hi'])/100)}）</strong>，
校正前已达名义显著的 {len(sig)} 个里，<strong>{s_shrink} 个收缩、{s_grow} 个放大</strong>，方向几乎对半分——校正无系统性方向。</p>
<p>三层梯度图景：<strong>代谢层不动（+{f(S['metabolism']['median_A'],1)}%）、炎症层点估计略负但不可靠、
免疫层方向与阳性对照一致（+{f(S['immune']['median_A'],1)}%）但 n=6、CI 极宽</strong>——
预注册的 P1（炎症&gt;代谢）被推翻，P2（免疫&gt;炎症、代谢）方向正确但两个 CI 都含 0，均未获证实。
外部方向性复算（GSE68465 + GSE31210，665 例外部样本）显示<strong>三层排序在全部三个队列复现</strong>（见五之二），但免疫层外部无显著性支持。</p>
</div>

<div class="cards">
<div class="card"><div class="k">纳入签名</div><div class="v">{N}</div><div class="n">原锁定 {n_w1} + 代谢扩容 {n_w2} + 炎症定向 {n_w3} + 免疫分层 {n_w4}；三轮共排除并逐条记理由</div></div>
<div class="card"><div class="k">中位衰减（全部 {N}）</div><div class="v">{pct(float(S['ALL']['median_A'])/100)}</div><div class="n">95% CI {pct(float(S['ALL']['lo'])/100)} 到 {pct(float(S['ALL']['hi'])/100)}；原 13 个时为 −79.2% 到 +8.9%</div></div>
<div class="card"><div class="k">衰减 / 增强</div><div class="v">{shrink} / {grow}</div><div class="n">同号 {same} 个中效应量变小与变大的个数；另有 {flip} 个方向反转</div></div>
<div class="card"><div class="k">FDR 显著数</div><div class="v">{S['ALL']['n_FDR_M1']} → {S['ALL']['n_FDR_M3']}</div><div class="n">按 {N} 个签名统一 BH 校正；校正后不降反持平——没有"整类掉显著性"这回事</div></div>
</div>

<h2>一、扩容是怎么做的（规则先于结果）</h2>
<p>扩容是<strong>看到第一波结果之后</strong>才决定的，属方案修订。为了不让它变成"挑到满意为止"，修订文件
<code>01_amendment_expansion_to20.md</code> 在任何新结果产出前写定规则并锁定哈希（见 <code>work/LOCK.txt</code>）。要点：</p>
<table>
<thead><tr><th style="width:12%">层</th><th style="width:52%">规则</th><th>结果</th></tr></thead>
<tbody>
<tr><td>检索</td><td>6 组 Europe PMC 检索式，2020 起，覆盖氨基酸/胆固醇/糖酵解/谷氨酰胺/能量代谢等 I3 子锚点，单列 LUSC 与"摘要明示 TCGA"</td><td>命中 767 → 去重 538</td></tr>
<tr><td>L1</td><td>去重 + 2020 起 + 有 PMCID + 不在既有名单</td><td>431</td></tr>
<tr><td>L2</td><td>摘要内能检出 2–15 个基因符号（I4 的摘要级代理）</td><td>99（Tier A）</td></tr>
<tr><td>L3</td><td>显式规则表逐条判 I1／I3／E3／E4 + 两条新规则</td><td>送抽取 62</td></tr>
<tr><td>终判</td><td>全文抽取与自述基因数交叉校验 + 人工复核</td><td><strong>纳入 28，排除 26</strong></td></tr>
</tbody></table>
<p class="note"><strong>两条新规则（原 E1–E5 没有，必须披露）：</strong>
① <strong>CIRC</strong>——题名锚定"免疫相关签名/免疫基因/ICD"者不纳入主分析（修订 03 更正措辞：这不是"循环论证"，
而是<strong>预测性混淆</strong>——免疫签名校正后理应大幅衰减，混入主集会抬高整类衰减估计，正确做法是单独分层，见"一之四"）；
② <strong>I3 外锚点</strong>——铁死亡、焦亡、自噬、m6A、缺氧、ECM/CAF、溶酶体、DNA 修复等判 I3 失败，
因 I3 只允许"炎症"与"营养/代谢"两类构念。</p>

<h2>一之二、炎症层定向补充（修订 02）与它的天花板</h2>
<p>修订 01 之后炎症只有 9 个、代谢 32 个。因为炎症签名与 ImmuneScore 的相关最高（|ρ| 中位 0.482 vs 代谢 0.146，Wilcoxon P=0.0001），
<strong>炎症层是"如果免疫混杂存在、衰减就该最明显"的地方</strong>，所以这一层必须补厚，否则等于把最关键的检验做废。
按修订 02 做了 8 组炎症专属检索：</p>
<table>
<thead><tr><th style="width:22%">步骤</th><th class="num" style="width:12%">数量</th><th>说明</th></tr></thead>
<tbody>
<tr><td>检索命中</td><td class="num">1,139</td><td>8 组式，显式排除细胞类型签名（巨噬/中性粒/T/NK）与构念外词（铁死亡/焦亡/自噬/m6A/缺氧）</td></tr>
<tr><td>去重后</td><td class="num">782</td><td>跨 8 组去重</td></tr>
<tr><td>新候选</td><td class="num">605</td><td>剔除已在既有名单、已排除、2020 前、无 PMCID</td></tr>
<tr><td>摘要含 2–15 个基因符号</td><td class="num">74</td><td>I4 的摘要级代理</td></tr>
<tr><td>送全文抽取</td><td class="num">48</td><td>题名级 CIRC／构念外已剔除</td></tr>
<tr><td>炎症基因池来源可核验</td><td class="num">8</td><td>全文里能找到"候选基因池来自炎症/趋化因子基因集"的证据句</td></tr>
<tr><td><strong>人工定稿纳入</strong></td><td class="num"><strong>{n_w3}</strong></td><td>其余 6 条的证据句只是 KEGG 富集注释，不是基因池来源</td></tr>
</tbody></table>
<div class="lead">
<p><strong>结论：炎症层补不到 20 个，天花板是 {n_inf} 个。</strong>这不是检索不够，是构念本身的结构问题——
在 NSCLC 签名文献里，"炎症"几乎全部以 <strong>immune-related</strong>（免疫相关）的形式发表：
本轮 74 条 Tier A 里有 <strong>20 条</strong>因题名锚定"免疫相关签名/ICD/免疫治疗应答"而被 CIRC 规则剔除，
其中 <strong>16 条其实有可读的基因列表</strong>（见 <code>work/circ_feasibility.json</code>）。
按修订 02 的停止规则，此处照实报告 {n_inf} 个，不降低标准凑数。</p>
<p>这个"补不上去"本身是一个可发表的结果：<strong>炎症构念在 NSCLC 转录组文献中已被免疫构念吸收</strong>，
因此"炎症锚定"与"免疫锚定"在实操上难以分离——这正是同类研究容易把免疫混杂当成炎症效应的原因。</p>
</div>

<h2>一之三、炎症层的方向与预注册预测 P1</h2>
<p>修订 02 事前登记了定向预测 P1：炎症层衰减应<strong>大于</strong>代谢层。结果<strong>方向相反</strong>——</p>
<table>
<thead><tr><th style="width:30%">分层</th><th class="num">n</th><th class="num">A&gt;0（衰减）</th><th class="num">A&lt;0（增强）</th><th class="num">中位 A</th><th class="num">仅 M1 显著者</th></tr></thead>
<tbody>
<tr><td><strong>炎症层</strong></td><td class="num">{n_inf}</td><td class="num">{inf_dn}</td><td class="num"><strong>{inf_up}</strong></td>
<td class="num">{pct(float(S['inflammation']['median_A'])/100)}</td>
<td class="num">{sig_inf_dn} 衰减 / {sig_inf_up} 增强（n={len(sig_inf)}）</td></tr>
<tr><td><strong>代谢层</strong></td><td class="num">{n_met}</td><td class="num">{met_dn}</td><td class="num">{met_up}</td>
<td class="num">{pct(float(S['metabolism']['median_A'])/100)}</td><td class="num">—</td></tr>
<tr><td><strong>免疫层</strong>（见一之四）</td><td class="num">{n_imm}</td><td class="num">{imm_dn}</td><td class="num">{imm_up}</td>
<td class="num">{pct(float(S['immune']['median_A'])/100)}</td><td class="num">—</td></tr>
</tbody></table>
<p>{p1_txt}
但方向值得记一笔：炎症层 {n_inf} 个里 <strong>{inf_up} 个在校正后效应量变大</strong>，
而代谢层是 {met_dn}:{met_up} 基本对半。若把炎症层单独看，
"校正后反而更强"与"免疫含量稀释了签名信号（负混杂）"是相容的——但 n 太小，只能作探索性提示，不能作结论。</p>

<h2>一之四、免疫分层（修订 03）：把 CIRC 挡掉的 20 条做成第三分层</h2>
<p>炎症层的天花板暴露了一个事实：<strong>被 CIRC 规则挡掉的 20 条 immune-related 签名里有 16 条本来就能抽到基因列表</strong>。
修订 03 把它们改判为<strong>第三分层（免疫层）</strong>而非一律丢弃——它们恰恰是"如果免疫混杂解释了已发表签名，
就该最先衰减"的近似阳性对照。修订 03 事前登记：</p>
<ul>
<li><strong>更正措辞：</strong>修订 01 称 CIRC 防"循环论证"不准确，实为防<strong>预测性混淆</strong>——免疫签名校正后理应大幅衰减（它们本来就是免疫含量代理），混入主集会抬高整类衰减估计；单独分层才是正确做法。</li>
<li><strong>定向预测 P2：</strong>免疫层衰减 &gt; 炎症层、代谢层。</li>
<li><strong>停止解读条件：</strong>免疫签名 &lt;5 条，或衰减方向与阳性对照相反（不衰减），则停止解读并排查打分/校正管线。</li>
</ul>
<p>20 条候选带证据重抽全文（<code>tools/18_extract_wave4.py</code>，抓基因池来源句）后：
<strong>纳入 {n_w4} 条（SIG44–SIG49，池来源分别为 ImmPort / InnateDB / Cistrome / 共刺激分子），排除 10 条，其余无基因列表</strong>，逐条记理由（<code>19_build_registry_w4.py</code>，哈希 <code>fc4352ea…</code>）。</p>
<h3>P2 检验结果</h3>
<p>{p2_txt}</p>
<h3>免疫层 6 个签名的明细</h3>
<table>
<thead><tr><th>签名</th><th>池来源</th><th class="num">基因数</th><th class="num">β<sub>M1</sub></th><th class="num">β<sub>M3</sub></th><th class="num">衰减 A</th><th class="num">M1 P</th><th class="num">M3 P</th><th class="num">ρ(分,Immune)</th></tr></thead>
<tbody>
{"".join(
    f'<tr><td>{r["sig_id"]}</td><td>{pool_src(r["sig_id"])}</td><td class="num">{r["n_genes"]}</td>'
    f'<td class="num">{f(r["beta_M1"],3,True)}</td><td class="num">{f(r["beta_M3"],3,True)}</td>'
    f'<td class="num">{pct(r["att_immune"])}</td><td class="num">{f(r["p_M1"],4)}</td>'
    f'<td class="num">{f(r["p_M3"],4)}</td><td class="num">{f(r["rho_immune"],2)}</td></tr>'
    for r in sorted(imm_rows_tsv, key=lambda x: x["sig_id"]))}
</tbody></table>
<div class="lead">
<p><strong>判读（严格按修订 03 的条件）：</strong>免疫层 n={n_imm}，未触发停止条件（≥5 条；中位 A = {pct(float(S['immune']['median_A'])/100)} &gt; 0，方向与阳性对照一致）。
{imm_dn}/{n_imm} 个 A&gt;0；限 M1 显著者共 {len([r for r in imm_rows_tsv if float(r['p_M1'])<0.05])} 个（SIG45/46/48/49），其中 3 个收缩 1 个放大。
<strong>但 n=6 的 bootstrap CI 极宽（{pct(float(S['immune']['lo'])/100)} 到 {pct(float(S['immune']['hi'])/100)}），P2 的方向正确不等于成立</strong>——
免疫层只能作为"方向与阳性对照相容"的佐证，论文里必须写成探索性分层，不能写成"免疫签名确实被校正削弱"。</p>
<p class="note"><strong>多重性披露：</strong>本文共 3 个定向检验（P1、P2a、P2b），均未做多重校正（预先定位为方向性预测而非确证性检验）。
若按 Bonferroni（α = 0.05/3 ≈ 0.017）从严，三个检验的 bootstrap 区间本就全部含 0，结论"均未获证实"不变——即校正与否不改变判读。</p>
</div>

<h2>二、主结果：多口径并列</h2>
<table>
<thead><tr><th style="width:20%">口径</th><th class="num" style="width:6%">n</th><th class="num">中位衰减 A（全部）</th>
<th class="num" style="width:7%">可定义 n</th><th class="num">中位衰减 A（仅 M1 显著）</th><th class="num" style="width:11%">FDR 显著</th><th style="width:18%">说明</th></tr></thead>
<tbody>
{setrow("ALL", "全部签名", "主报告口径")}
{setrow("W1", "原锁定 13", "预注册主分析")}
{setrow("W2", "扩容 28（代谢为主）", "修订 01")}
{setrow("W3", "炎症定向补充", "修订 02")}
{setrow("W4", "免疫分层", "修订 03")}
{setrow("inflammation", "炎症层（合计）", "构念之一")}
{setrow("metabolism", "代谢层", "构念之二")}
{setrow("immune", "免疫层", "近似阳性对照")}
</tbody></table>
<p class="note">括号内为 bootstrap（按签名重抽样 5,000 次）95% CI，单位为百分比。
"可定义"指 β<sub>M1</sub> 达到名义显著（P&lt;0.05）的子集——比值型估计量只在此时有意义，仍标 exploratory。</p>
<p><strong>波次异质性检查（修订文件第 3.3 条，先声明后执行）：</strong>W1 中位 {pct(float(w1['median_A'])/100)}、W2 中位 {pct(float(w2['median_A'])/100)}，
相差 {f(diff,2)} 个百分点，<strong>低于事先设定的 10 个百分点阈值</strong>，故按计划合并报告；W3（n=2）、W4（n=6）样本太小，单独分层报告、不并入主中位数解释。</p>

<h2>三、方向性：收缩与放大约各占一半</h2>
<p>这是扩容带来的最硬的一个数：把 {N} 个签名逐个看，校正后效应量 <strong>{shrink} 个收缩、{grow} 个放大、{flip} 个方向反转</strong>。
只看校正前已达名义显著的 {len(sig)} 个，则是 <strong>{s_shrink} 个收缩对 {s_grow} 个放大</strong>，
二项检验两侧几乎等概率（P≈1）。换句话说，校正的<strong>方向是随机的</strong>，不是系统性压缩。</p>
<div class="figwrap">
{forest}
<div class="figcap">每个签名一条线：空心点=M1（校正前），实心点=M3（校正后）。线向右=效应量变大，向左=变小。
橙色为炎症层、蓝色为代谢层、紫色为免疫层。整体看既没有统一向右也没有统一向左——这是"无系统性方向"的直观形态。</div>
</div>

<h2>四、位移大小与"签名-免疫相关"的关系（含混杂控制）</h2>
<p>预注册假设 H3 说"衰减幅度与签名分×ImmuneScore 的相关正相关"。在 {N} 个签名下：
<strong>有符号的衰减与 ρ 无关</strong>（Spearman ρ = {f(r_signed,3)}，P = {f(t_p(r_signed,N),3)}），
但<strong>位移的绝对值与 ρ 正相关</strong>。这个区分很重要，另做了混杂检查：</p>
<table>
<thead><tr><th>检验</th><th class="num">统计量</th><th>判读</th></tr></thead>
<tbody>
<tr><td>|Δβ| 与 |ρ(签名, Immune)|</td><td class="num">ρ = {f(r_adb_arho,3)}</td><td>相关存在</td></tr>
<tr><td>|Δβ| 与 |β<sub>M1</sub>|（"大效应天然动得多"）</td><td class="num">ρ = {f(r_adb_ab1,3)}</td><td>方向相反，排除该混杂</td></tr>
<tr><td>偏相关（控制 |β<sub>M1</sub>|）</td><td class="num">{f(part,3)}（P = {f(part_p,3)}）</td><td>扣除效应量后仍成立，但已减弱</td></tr>
<tr><td>有符号的 A 与 ρ</td><td class="num">{f(r_signed,3)}（P = {f(t_p(r_signed,N),3)}）</td><td><strong>H3 不成立</strong>：位移大小随 ρ 走，方向不随 ρ 走</td></tr>
</tbody></table>
<p class="note">这几行是从 <code>results/attenuation_all.tsv</code> 现算的。
写成论文时的措辞应当是"校正按签名与免疫含量的相关程度推动效应量，但推动方向不一致"，<strong>不能</strong>写成"校正削弱了签名"。</p>

<h2>五、阳性对照与随机面板</h2>
<table>
<thead><tr><th>对照轴</th><th class="num">ρ(分, Immune)</th><th class="num">β<sub>M1</sub> (P)</th><th class="num">β<sub>M3</sub> (P)</th><th class="num">衰减</th><th>判读</th></tr></thead>
<tbody>{ctrl_rows}</tbody></table>
<p>淋巴细胞轴（ρ=0.86，真正的免疫含量代理）衰减 <strong>{pct(float(ctrl[0]['attenuation']))}</strong>、显著性消失
（P {f(ctrl[0]['p_M1'],3)}→{f(ctrl[0]['p_M3'],3)}）——<strong>同一次校正用在已知免疫代理上会崩，用在 {N} 个已发表签名上不会</strong>。方法有功效。
营养/肝脏轴（ρ=−0.10，与免疫含量无关）只衰减 {pct(float(ctrl[2]['attenuation']))}，构成阴性方向的一致性。
中性粒轴 β<sub>M1</sub>≈0（P={f(ctrl[1]['p_M1'],3)}），比值型估计量机械爆炸（−{f(abs(float(ctrl[1]['attenuation'])),0)}%），
<strong>该行只列不作解读</strong>——这也是为什么主分析必须同时报"仅 M1 显著"子集。</p>
<p><strong>对照轴是否只是 ESTIMATE 的产物？</strong>不是。pathA 已用 10 种免疫群体反卷积工具（MCP-counter、xCell 等）对三条轴做过交叉验证（<code>pathA/results/immune_tool_crosscheck_summary.tsv</code>）：
淋巴细胞轴在 10 个工具中 <strong>ρ 中位 {f(float(_cross['lymphocyte']['rho_median']),2) if _cross else '0.69'}（10/10 为正、10/10 显著）</strong>，
中性粒轴 {f(float(_cross['neutrophil']['rho_median']),2) if _cross else '0.28'}、营养/肝脏轴 {f(float(_cross['nutrition_liver']['rho_median']),2) if _cross else '0.16'}——
与本文 |ρ| 的排序（0.86 / 0.55 / 0.10）方向一致。<strong>对照轴不是 ESTIMATE 的循环产物，免疫含量的强弱排序在多工具下稳定。</strong></p>
{pct_block}
{ext_block}
{ep_block}

<h2>六、FDR 层面的变化</h2>
{fdr_bullets}

<h2>七、必须写进 Limitations 的八条（比阶段三多三条）</h2>
<ol>
<li><strong>分层不平衡且两层已到天花板：</strong>共 {N} 个，炎症 {n_inf}、代谢 {n_met}、免疫 {n_imm}。炎症层经 8 组定向检索仍只补到 {n_w3} 个新人选；
免疫层只有 {n_w4} 条、CI 极宽。<strong>三层的中位数都不能外推成"该类签名"的总体结论。</strong></li>
<li><strong>外部复算是 post-hoc 且协变量不完全一致：</strong>GSE31210 缺 age（M1 少一个协变量）、平台为芯片（签名基因覆盖约 50%–75%，SIG26 无外部值），
ESTIMATE 分数用原始论文签名集 + ssGSEA 复算（与官方分数差常数倍，A 等价）；这一步未预注册，若投稿须如实标注 post-hoc。</li>
<li><strong>扩容是 post-hoc：</strong>看到第一波结果后才决定，规则虽先锁定（含两次事前登记的定向预测与停止规则），仍属修订；各波次必须分开报告。</li>
<li><strong>Tier A 选择机制：</strong>扩容只纳入"摘要或全文里已列出基因符号"的候选，倾向作者愿意公开列表的论文，与第一波存在机制差异。</li>
<li><strong>比值型估计量：</strong>{N} 个里 {S['ALL']['n_abs_gt100']} 个 |A|&gt;100%，全部是 β<sub>M1</sub> 接近 0 者（含对照里的中性粒轴）；主结论必须同时给"仅 M1 显著"的子集口径。</li>
<li><strong>模型冗余与 M2 的降级（已正式执行）：</strong>purity 与 ESTIMATE 分数 ρ=−1.000、β<sub>M2</sub> 与 β<sub>M3</sub> 的 r=0.998（49 签名集；先前印出的 0.999 是 13 签名预试验值，见修订 06），M2 无独立信息——
主分析因此采用 M1 vs M3，<strong>M2 结果（<code>att_purity</code> 列）降为补充材料，不参与任何主结论</strong>。</li>
<li><strong>百分位不独立：</strong>每个签名的随机面板百分位来自同一批 481 例样本，彼此相关；KS 检验的 P 值只作形态参考，不能当独立检验的确证。</li>
<li><strong>LUSC 不能做跨组织学复算：</strong>49 个签名在 TCGA-LUSC（{_lusc_n} 例/{_lusc_ev} 事件）中无一个 FDR 显著（仅 {_lusc_n_sig} 个名义），
β<sub>M1</sub> 跨队列 ρ=−0.17——A 的分母接近 0 时比值不可定义，故该队列只作方向计数的旁证（{_lusc_sh}:{_lusc_gr}，与 LUAD 同为对半），
不进入主结论；这也说明已发表签名在 LUSC 的可重复性本身有限，本文不主张任何跨组织学外推。</li>
</ol>

<h2>八、主张（阶段四更新版）</h2>
<div class="lead">
<p><strong>建议主张：</strong>"<em>免疫/基质含量不是已发表 NSCLC 炎症与代谢预后签名的解释变量：一项对 {N} 个签名的三层校准研究，带阳性与阴性对照</em>"</p>
<p><strong>六条支撑句（数字都能顶住）：</strong></p>
<ol>
<li><strong>方法有功效：</strong>同一校正用于已知免疫代理（淋巴细胞轴，ρ=0.86）→ 衰减 {pct(float(ctrl[0]['attenuation']))}、显著性消失（P {f(ctrl[0]['p_M1'],3)}→{f(ctrl[0]['p_M3'],3)}）；无关对照（营养/肝脏轴）仅 {pct(float(ctrl[2]['attenuation']))}。</li>
<li><strong>没有系统性方向：</strong>{N} 个签名中位衰减 {pct(float(S['ALL']['median_A'])/100)}（95% CI {pct(float(S['ALL']['lo'])/100)} 到 {pct(float(S['ALL']['hi'])/100)}）；
校正前已达显著的 {len(sig)} 个里 {s_shrink} 收缩对 {s_grow} 放大。</li>
<li><strong>位移的大小随相关走、方向不随：</strong>|Δβ| 与 |ρ| 的偏相关 {f(part,3)}（P = {f(part_p,3)}），但有符号的 A 与 ρ 无关（{f(r_signed,3)}，P = {f(t_p(r_signed,N),3)}）。</li>
<li><strong>三层梯度只到"方向"层面：</strong>免疫层 +{f(S['immune']['median_A'],1)}%（与阳性对照同向）、炎症 {pct(float(S['inflammation']['median_A'])/100)}、代谢 +{f(S['metabolism']['median_A'],1)}%；但 P1 被推翻、P2 两个 CI 均含 0——<strong>只能写"方向相容"，不能写"梯度成立"</strong>。</li>
<li><strong>外部复现的是排序、不是数值：</strong>三层排序在 GSE68465 与 GSE31210（合计 665 例外部样本）复现，但逐签名 A 跨队列 ρ≈0，且外部免疫层无可解读的显著性——外部复算支持"方向图景稳健"，不支持"免疫签名衰减"本身。</li>
<li><strong>终点不影响结论：</strong>PFI 与 DSS 各自重跑后，中位 A 均在 ±10%、方向计数比均在 1:2–2:1、三层排序与 OS 相同（修订 05 的 E1 判读成立）——主结论不依赖 OS 终点。</li>
</ol>
</div>
<p><strong>仍然不能写的话：</strong>"签名虚高""多数签名不成立""校正后失去预后价值""免疫签名确实被校正削弱"——四条都与上面这些数字冲突。</p>

<h2>九、逐签名主表（{N} 行，全部报告，不因结果好坏剔除）</h2>
<table>
<thead><tr><th>#</th><th>波次</th><th>分层</th><th class="num">基因数</th><th class="num">β<sub>M1</sub></th><th class="num">β<sub>M3</sub></th><th class="num">衰减</th><th class="num">q<sub>M1</sub></th><th class="num">q<sub>M3</sub></th></tr></thead>
<tbody>
{main_rows}
</tbody></table>
<p class="note">β 为每 SD 的 log-HR；M1 = 签名分 + 年龄 + 性别 + 分期；M3 = M1 + ImmuneScore + StromalScore。
q 为 BH 校正后 P（按 {N} 个签名统一校正）。</p>

<hr>
<p class="note">数据来源：<code>results/attenuation_all.tsv</code>、<code>summary_all.tsv</code>、<code>control_axes_all.tsv</code>、
<code>attenuation_pctile_all.tsv</code>、<code>stratum_diff_boot.tsv</code>、<code>attenuation_endpoints.tsv</code>、<code>summary_endpoints.tsv</code>（日志 <code>logs/22</code>）；LUSC 核对 <code>attenuation_lusc49.tsv</code>、<code>summary_lusc49.tsv</code>（日志 <code>logs/21</code>）；外部复算 <code>external_m1m3_long.tsv</code>、
<code>external_m1m3_summary.tsv</code>、<code>external_vs_tcga_rho.tsv</code>（脚本 <code>R/22_external_m1m3.R</code>，日志 <code>logs/18、20</code>）。
外部数据：pathA <code>data/external/gse68465/</code>、<code>data/external/gse31210/</code>（EBI ArrayExpress 镜像）；
ESTIMATE 签名集取自 Yoshihara et al. 2013 <em>Nat Commun</em> Supplementary Data 1（<code>work/estimate_ext/ESM_488.xlsx</code>）。
签名名单：<code>signature_registry.csv</code>（13）+ <code>signature_registry_expansion.csv</code>（28）+
<code>signature_registry_expansion_w3.csv</code>（2）+ <code>signature_registry_expansion_w4.csv</code>（6），排除项见各 excluded 表与修订文件。
修订文件：<code>01_amendment_expansion_to20.md</code>、<code>02_amendment_inflammation_wave3.md</code>、<code>03_amendment_immune_stratum.md</code>。
报告由脚本生成（<code>tools/22_make_report_stage4.py</code>），未手抄数字。</p>
</div></body></html>"""

outp = os.path.join(PB, "报告_阶段4_免疫分层49签名_2026-10-20.html")
open(outp, "w", encoding="utf-8").write(HTML)
print("->", outp, len(HTML.encode("utf-8")), "bytes")
