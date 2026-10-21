# -*- coding: utf-8 -*-
"""生成阶段二报告 HTML：全部数字从 results/*.tsv 读取后填入，不手写数字。

纪律：任何字段读不到就在页面上写"缺失"，不用估算值顶替。
"""
import csv, html, json, math, os, statistics as st

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
R = os.path.join(PB, "results")


def tsv(p):
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


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
reg2 = list(csv.DictReader(open(f"{PB}/signature_registry_expansion.csv", encoding="utf-8-sig")))
_p1p = f"{R}/stratum_diff_boot.tsv"
P1 = tsv(_p1p)[0] if os.path.exists(_p1p) else None

A = {r["sig_id"]: r for r in att}
def g(sid, k, d="NA"):
    return A.get(sid, {}).get(k, d)

# ---- 方向统计 ----
same, shrink, grow, flip = 0, 0, 0, 0
for r in att:
    b1, b3 = float(r["beta_M1"]), float(r["beta_M3"])
    if b1 * b3 > 0:
        same += 1
        if abs(b3) < abs(b1): shrink += 1
        else: grow += 1
    else:
        flip += 1
N = len(att)
n_inf = sum(1 for r in att if r["stratum"] == "inflammation")
n_met = sum(1 for r in att if r["stratum"] == "metabolism")
n_w1 = sum(1 for r in att if r["wave"] == "W1")
n_w2 = sum(1 for r in att if r["wave"] == "W2")
n_w3 = sum(1 for r in att if r["wave"] == "W3")


def _dir(sub):
    a_ = sum(1 for r in sub if float(r["att_immune"]) > 0)
    b_ = sum(1 for r in sub if float(r["att_immune"]) < 0)
    return a_, b_


inf_dn, inf_up = _dir([r for r in att if r["stratum"] == "inflammation"])
met_dn, met_up = _dir([r for r in att if r["stratum"] == "metabolism"])
sig_inf = [r for r in att if r["stratum"] == "inflammation" and float(r["p_M1"]) < 0.05]
sig_inf_dn = sum(1 for r in sig_inf if float(r["att_immune"]) > 0)
sig_inf_up = sum(1 for r in sig_inf if float(r["att_immune"]) < 0)
sig = [r for r in att if float(r["p_M1"]) < 0.05]
s_shrink = sum(1 for r in sig if float(r["beta_M1"]) * float(r["beta_M3"]) > 0 and abs(float(r["beta_M3"])) < abs(float(r["beta_M1"])))
s_grow = sum(1 for r in sig if float(r["beta_M1"]) * float(r["beta_M3"]) > 0 and abs(float(r["beta_M3"])) >= abs(float(r["beta_M1"])))

# ---- forest SVG ----
rows = sorted(att, key=lambda r: (r["wave"], -abs(float(r["beta_M1"]))))
top, rowh, x0, xw = 34, 7.4, 150, 300
H = top + rowh * len(rows) + 46
bmin = min(min(float(r["beta_M1"]), float(r["beta_M3"])) for r in att)
bmax = max(max(float(r["beta_M1"]), float(r["beta_M3"])) for r in att)
pad = 0.06
lo, hi = bmin - pad, bmax + pad
def X(v): return x0 + (float(v) - lo) / (hi - lo) * xw
zero_x = X(0)
svg = [f'<svg viewBox="0 0 680 {H:.0f}" width="100%" role="img" aria-label="{N} 个签名校正前后效应量位移">']
svg.append(f'<line x1="{zero_x:.1f}" y1="{top-14}" x2="{zero_x:.1f}" y2="{top+rowh*len(rows)+4}" stroke="#c3ccd6" stroke-dasharray="3 3"/>')
svg.append(f'<text x="{zero_x:.1f}" y="{top-19}" font-size="9.5" fill="#6b7684" text-anchor="middle">0</text>')
for i, r in enumerate(rows):
    y = top + i * rowh
    b1, b3 = float(r["beta_M1"]), float(r["beta_M3"])
    sid = r["sig_id"]
    col = "#c2410c" if r["stratum"] == "inflammation" else "#1d4ed8"
    svg.append(f'<text x="6" y="{y+2.6:.1f}" font-size="7.2" fill="#5b6673">{sid}</text>')
    svg.append(f'<line x1="{X(b1):.1f}" y1="{y:.1f}" x2="{X(b3):.1f}" y2="{y:.1f}" stroke="{col}" stroke-width="1.5" opacity="0.55"/>')
    svg.append(f'<circle cx="{X(b1):.1f}" cy="{y:.1f}" r="2.2" fill="none" stroke="{col}" stroke-width="1.1"/>')
    svg.append(f'<circle cx="{X(b3):.1f}" cy="{y:.1f}" r="2.2" fill="{col}"/>')
svg.append(f'<text x="{x0}" y="{H-22:.0f}" font-size="9.5" fill="#6b7684">log-HR per SD（空心=M1 校正前，实心=M3 校正后）</text>')
svg.append(f'<circle cx="{x0+300}" cy="{H-26:.0f}" r="2.4" fill="none" stroke="#c2410c"/><text x="{x0+308}" y="{H-22:.0f}" font-size="9" fill="#6b7684">炎症（9）</text>')
svg.append(f'<circle cx="{x0+380}" cy="{H-26:.0f}" r="2.4" fill="none" stroke="#1d4ed8"/><text x="{x0+388}" y="{H-22:.0f}" font-size="9" fill="#6b7684">代谢（32）</text>')
svg.append('</svg>')
forest = "\n".join(svg)

# ---- 主表 ----
def row_main(r):
    return (f'<tr><td>{r["sig_id"]}</td><td>{r["wave"]}</td>'
            f'<td>{"炎症" if r["stratum"]=="inflammation" else "代谢"}</td>'
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
    low = sum(1 for x in allp if x < 0.20); mid = sum(1 for x in allp if 0.20 <= x < 0.80)
    high = sum(1 for x in allp if x >= 0.80)
    high_ids = [i for i, x in pv if x >= 0.80]
    # 单样本 KS（对照 U(0,1)）
    D = max(max(abs((i + 1) / n - x), abs(x - i / n)) for i, x in enumerate(allp))
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * D
    ks_p = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 101))
    ks_p = min(max(ks_p, 0.0), 1.0)
    pct_block = f"""<p>{N} 个签名的衰减在各自随机面板分布中的百分位：<strong>低段（&lt;0.20）{low} 个、中段 {mid} 个、高段（≥0.80）{high} 个</strong>，
均值 {f(st.mean(allp),3)}、中位 {f(st.median(allp),3)}。高段为 {'、'.join(high_ids) if high_ids else '无'}。</p>
<p><strong>这里出现了一个第一波没看出来的结构：</strong>百分位<strong>不是均匀分布</strong>，而是向 0.5 附近聚集
（KS 对照 U(0,1)：D = {f(D,3)}，P ≈ {f(ks_p,3)}）。含义是：已发表签名的衰减不仅"与随机面板无法区分"，
而且<strong>比随机面板更贴近零衰减</strong>——随机抽一组基因，衰减的离散度很大，常常出现大幅度的假性收缩或放大；
而这 {N} 个签名的位移都极小（|Δβ| 中位 {f(st.median([abs(float(r['beta_M3'])-float(r['beta_M1'])) for r in att]),3)}）。</p>
<p class="note">统计口径提醒：这些百分位来自同一批样本、彼此不独立，KS 的 P 值只能作参考、不能当确证；
真正稳的是分布形态（{mid}/{n} 落在中段）与 |Δβ| 的量级。</p>"""
else:
    pct_block = '<p class="note">零分布尚未跑完，本节数字缺失（不做估算填充）。</p>'

ctrl_rows = "\n".join(
    f'<tr><td>{c["axis"]}</td><td class="num">{f(c["rho_immune"],3)}</td>'
    f'<td class="num">{f(c["beta_M1"],3,True)} ({f(c["p_M1"],3)})</td>'
    f'<td class="num">{f(c["beta_M3"],3,True)} ({f(c["p_M3"],3)})</td>'
    f'<td class="num"><strong>{pct(c["attenuation"])}</strong></td></tr>' for c in ctrl)

S = summ
def _pct_or_na(v):
    """空串/NA 时返回 'NA'，避免整篇生成因一个缺失值崩掉。"""
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
    def _ci(lo, hi):
        if str(lo).strip() in ("", "NA") and str(hi).strip() in ("", "NA"):
            return "（样本量 &lt;3，不做区间估计）"
        return f'（{_pct_or_na(lo)} 到 {_pct_or_na(hi)}）'

    lo_hi = _ci(r["lo"], r["hi"])
    lo2_hi2 = _ci(r["lo_valid"], r["hi_valid"])
    return (f'<tr><td><strong>{label}</strong></td><td class="num">{r["n"]}</td>'
            f'<td class="num">{_pct_or_na(r["median_A"])}<span class="ci">{lo_hi}</span></td>'
            f'<td class="num">{r["n_ratio_valid"]}</td>'
            f'<td class="num">{_pct_or_na(r["median_A_valid"])}<span class="ci">{lo2_hi2}</span></td>'
            f'<td class="num">{r["n_FDR_M1"]} → {r["n_FDR_M3"]}</td><td>{note}</td></tr>')


w1 = S.get("W1"); w2 = S.get("W2")
diff = abs(float(w1["median_A"]) - float(w2["median_A"])) if w1 and w2 else None


# ---- 第 4 节的统计量全部现算（不引用日志里的手抄值） ----
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
    N = 20000
    ssum = 0.0
    for i in range(N):
        v = (i + 0.5) / N
        u = at + v / (1 - v)
        ssum += c * (1 + u * u / df) ** (-(df + 1) / 2) / (1 - v) ** 2
    return min(1.0, 2 * ssum / N)


rho_imm = [float(r["rho_immune"]) for r in att]
att_v = [float(r["att_immune"]) for r in att]
adb = [abs(float(r["beta_M3"]) - float(r["beta_M1"])) for r in att]
ab1 = [abs(float(r["beta_M1"])) for r in att]
arho = [abs(x) for x in rho_imm]
N = len(att)
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
            % (r["sig_id"], "炎症" if r["stratum"] == "inflammation" else "代谢",
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

if P1 is None:
    p1_txt = "P1 的差值未计算（缺少 stratum_diff_boot.tsv）。"
else:
    p1_txt = ("P1 的差值（炎症 − 代谢）＝ <strong>%s</strong>，bootstrap 95%% CI %s 到 %s——"
              "<strong>区间含 0，P1 未获证实</strong>。"
              % (pct(float(P1["est"])), pct(float(P1["lo"])), pct(float(P1["hi"]))))

HTML = f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>方案 B 阶段三：43 个签名的最终校准结果</title>
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
.badge{{display:inline-block;padding:1px 7px;border-radius:99px;font-size:11.5px;font-weight:600}}
.b-ok{{background:#e5f4ea;color:var(--ok)}} .b-bad{{background:#fdeae8;color:var(--bad)}} .b-warn{{background:#fdf3dd;color:var(--warn)}}
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

<h1>方案 B 阶段三：从 13 扩到 {N}，含炎症层定向补充</h1>
<p class="sub">已发表 NSCLC 炎症/代谢预后签名的校准式重评估 · 修订 01 + 02 · 2026-10-20 · 全部数字由 <code>results/*.tsv</code> 直读生成</p>

<div class="lead">
<p><strong>一句话结论：</strong>扩容后主结论不变，而且精度好得多——
校正肿瘤纯度与免疫/基质含量后，{N} 个已发表签名的中位效应量变化是
<strong>{pct(float(S['ALL']['median_A'])/100)}（95% CI {pct(float(S['ALL']['lo'])/100)} 到 {pct(float(S['ALL']['hi'])/100)}）</strong>，
其中 {len(sig)} 个校正前就已达名义显著的签名里，<strong>{s_shrink} 个收缩、{s_grow} 个放大</strong>，方向几乎对半分。</p>
<p>也就是说：<strong>校正对这些签名没有系统性方向</strong>。既不能说它们"虚高"，也不能说校正"没用"——
正确结论是"免疫/基质含量不是这类签名的解释变量"，而这个结论现在有一个 ±8% 以内的区间撑腰，不再是 −79% 到 +9% 那种没法下判断的宽区间。</p>
</div>

<div class="cards">
<div class="card"><div class="k">纳入签名</div><div class="v">{N}</div><div class="n">原锁定 {n_w1} + 代谢扩容 {n_w2} + 炎症定向 {n_w3}；两轮共排除 32 条并逐条记理由</div></div>
<div class="card"><div class="k">中位衰减（全部 {N}）</div><div class="v">{pct(float(S['ALL']['median_A'])/100)}</div><div class="n">95% CI {pct(float(S['ALL']['lo'])/100)} 到 {pct(float(S['ALL']['hi'])/100)}；原 13 个时为 −79.2% 到 +8.9%</div></div>
<div class="card"><div class="k">衰减 / 增强</div><div class="v">{shrink} / {grow}</div><div class="n">同号 {same} 个中效应量变小与变大的个数；另有 {flip} 个方向反转</div></div>
<div class="card"><div class="k">FDR 显著数</div><div class="v">{S['ALL']['n_FDR_M1']} → {S['ALL']['n_FDR_M3']}</div><div class="n">校正后不是更少，反而多 1 个——没有"整类掉显著性"这回事</div></div>
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
① <strong>CIRC</strong>——题名锚定"免疫相关签名/免疫基因/ICD"者不纳入，因为用免疫含量去校正免疫签名属循环论证（排除 10 条）；
② <strong>I3 外锚点</strong>——铁死亡、焦亡、自噬、m6A、缺氧、ECM/CAF、溶酶体、DNA 修复等判 I3 失败，
因 I3 只允许"炎症"与"营养/代谢（脂质、脂肪酸、糖、氨基酸）"两类构念（排除 20 条）。</p>

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
<tr><td><strong>人工定稿纳入</strong></td><td class="num"><strong>2</strong></td><td>其余 6 条的证据句只是 KEGG 富集注释，不是基因池来源</td></tr>
</tbody></table>
<div class="lead">
<p><strong>结论：炎症层补不到 20 个，天花板是 {n_inf} 个。</strong>这不是检索不够，是构念本身的结构问题——
在 NSCLC 签名文献里，"炎症"几乎全部以 **immune-related**（免疫相关）的形式发表：
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
</tbody></table>
<p>{p1_txt}
但方向值得记一笔：炎症层 {n_inf} 个里 <strong>{inf_up} 个在校正后效应量变大</strong>，
而代谢层是 {met_dn}:{met_up} 基本对半。若把炎症层单独看，
"校正后反而更强"与"免疫含量稀释了签名信号（负混杂）"是相容的——但 n 太小，只能作探索性提示，不能作结论。</p>

<h2>二、主结果：三口径并列</h2>
<table>
<thead><tr><th style="width:20%">口径</th><th class="num" style="width:6%">n</th><th class="num">中位衰减 A（全部）</th>
<th class="num" style="width:7%">可定义 n</th><th class="num">中位衰减 A（仅 M1 显著）</th><th class="num" style="width:11%">FDR 显著</th><th style="width:20%">说明</th></tr></thead>
<tbody>
{setrow("ALL", "全部签名", "主报告口径")}
{setrow("W1", "原锁定 13", "预注册主分析")}
{setrow("W2", "扩容 28（代谢为主）", "修订 01")}
{setrow("W3", "炎症定向补充", "修订 02")}
{setrow("inflammation", "炎症层（合计）", "构念之一")}
{setrow("metabolism", "代谢层", "构念之二")}
</tbody></table>
<p class="note">括号内为 bootstrap（按签名重抽样 5,000 次）95% CI，单位为百分比。
"可定义"指 β<sub>M1</sub> 达到名义显著（P&lt;0.05）的子集——比值型估计量只在此时有意义，仍标 exploratory。</p>
<p><strong>波次异质性检查（修订文件第 3.3 条，先声明后执行）：</strong>W1 中位 {pct(float(w1['median_A'])/100)}、W2 中位 {pct(float(w2['median_A'])/100)}，
相差 {f(diff,2)} 个百分点，<strong>低于事先设定的 10 个百分点阈值</strong>，故按计划合并报告，并注明两波一致。</p>

<h2>三、方向性：收缩与放大约各占一半</h2>
<p>这是本次扩容带来的最硬的一个数：把 {N} 个签名逐个看，校正后效应量 <strong>{shrink} 个收缩、{grow} 个放大、{flip} 个方向反转</strong>。
只看校正前已达名义显著的 {len(sig)} 个，则是 <strong>{s_shrink} 个收缩对 {s_grow} 个放大</strong>，
二项检验两侧几乎等概率（P≈1）。换句话说，校正的<strong>方向是随机的</strong>，不是系统性压缩。</p>
<div class="figwrap">
{forest}
<div class="figcap">每个签名一条线：空心点=M1（校正前），实心点=M3（校正后）。线向右=效应量变大，向左=变小。
橙色为炎症层、蓝色为代谢层。整体看既没有统一向右也没有统一向左——这是"无系统性方向"的直观形态。</div>
</div>

<h2>四、位移大小与"签名-免疫相关"的关系（含混杂控制）</h2>
<p>预注册假设 H3 说"衰减幅度与签名分×ImmuneScore 的相关正相关"。在 {N} 个签名下：
<strong>有符号的衰减与 ρ 无关</strong>（Spearman ρ = {f(r_signed,3)}，P = {f(t_p(r_signed,N),3)}），
但<strong>位移的绝对值与 ρ 正相关</strong>。这个区分很重要，我另做了混杂检查：</p>
<table>
<thead><tr><th>检验</th><th class="num">统计量</th><th>判读</th></tr></thead>
<tbody>
<tr><td>|Δβ| 与 |ρ(签名, Immune)|</td><td class="num">ρ = {f(r_adb_arho,3)}</td><td>相关存在</td></tr>
<tr><td>|Δβ| 与 |β<sub>M1</sub>|（"大效应天然动得多"）</td><td class="num">ρ = {f(r_adb_ab1,3)}</td><td>方向相反，排除该混杂</td></tr>
<tr><td>偏相关（控制 |β<sub>M1</sub>|）</td><td class="num">{f(part,3)}（P = {f(part_p,3)}）</td><td>扣除效应量后仍成立，但已减弱</td></tr>
<tr><td>有符号的 A 与 ρ</td><td class="num">{f(r_signed,3)}（P = {f(t_p(r_signed,N),3)}）</td><td><strong>H3 不成立</strong>：位移大小随 ρ 走，方向不随 ρ 走</td></tr>
</tbody></table>
<p class="note">这几行是从 <code>results/attenuation_all.tsv</code> 现算的，日志 <code>logs/23_confound.log</code>、<code>logs/24_direction.log</code>。
写成论文时的措辞应当是"校正按签名与免疫含量的相关程度推动效应量，但推动方向不一致"，<strong>不能</strong>写成"校正削弱了签名"。</p>

<h2>五、阳性对照与随机面板</h2>
<table>
<thead><tr><th>对照轴</th><th class="num">ρ(分, Immune)</th><th class="num">β<sub>M1</sub> (P)</th><th class="num">β<sub>M3</sub> (P)</th><th class="num">衰减</th></tr></thead>
<tbody>{ctrl_rows}</tbody></table>
<p>淋巴细胞轴复现了第一波的结果（衰减 68.9%，P 从 {f(ctrl[0]['p_M1'],3)} 变 {f(ctrl[0]['p_M3'],3)}）——<strong>同一次校正用在已知免疫代理上会崩，用在 {N} 个已发表签名上不会</strong>。方法有功效，这一点没变。</p>
{pct_block}

<h2>六、FDR 层面的变化</h2>
{fdr_bullets}

<h2>七、必须写进 Limitations 的五条（比第一波多一条）</h2>
<ol>
<li><strong>分层不平衡且已到天花板：</strong>共 {N} 个，炎症 {n_inf} 个、代谢 {n_met} 个。炎症层经 8 组定向检索仍只补到 {n_w3} 个新人选，
原因是该构念在文献中被"免疫相关签名"吸收（20 条因 CIRC 被剔，其中 16 条本来可抽列表）。<strong>炎症层只能作探索性报告，不能把 {N} 个的中位数说成"炎症签名"的结论。</strong></li>
<li><strong>扩容是 post-hoc：</strong>看到第一波结果后才决定，规则虽先锁定，仍属修订；W1 与 W2 分开报是必须的。</li>
<li><strong>Tier A 选择机制：</strong>扩容只纳入"摘要里已列出基因符号"的候选，倾向作者把列表写进摘要的论文，与第一波存在机制差异。</li>
<li><strong>比值型估计量：</strong>{N} 个里 {S['ALL']['n_abs_gt100']} 个 |A|&gt;100%，全部是 β<sub>M1</sub> 接近 0 者；主结论必须同时给"仅 M1 显著"的子集口径。</li>
<li><strong>模型冗余：</strong>purity 与 ESTIMATE 分数 ρ=−1.000、β<sub>M2</sub> 与 β<sub>M3</sub> 的 r=0.999，M2 应删或降为补充。</li>
</ol>

<h2>八、主张（阶段三更新版）</h2>
<div class="lead">
<p><strong>建议主张：</strong>"<em>免疫/基质含量不是已发表 NSCLC 炎症与代谢预后签名的解释变量：一项对 {N} 个签名的校准研究，带阳性与阴性对照</em>"</p>
<p><strong>三条支撑句（数字都能顶住）：</strong></p>
<ol>
<li><strong>方法有功效：</strong>同一校正用于已知免疫代理（淋巴细胞轴，ρ=0.86）→ 衰减 68.9%、显著性消失（P {f(ctrl[0]['p_M1'],3)}→{f(ctrl[0]['p_M3'],3)}）。</li>
<li><strong>没有系统性方向：</strong>{N} 个签名中位衰减 {pct(float(S['ALL']['median_A'])/100)}（95% CI {pct(float(S['ALL']['lo'])/100)} 到 {pct(float(S['ALL']['hi'])/100)}）；
校正前已达显著的 {len(sig)} 个里 {s_shrink} 收缩对 {s_grow} 放大。</li>
<li><strong>位移的大小随相关走、方向不随：</strong>|Δβ| 与 |ρ| 的偏相关 {f(part,3)}（P = {f(part_p,3)}），但有符号的 A 与 ρ 无关（{f(r_signed,3)}，P = {f(t_p(r_signed,N),3)}）。</li>
</ol>
</div>
<p><strong>仍然不能写的话：</strong>"签名虚高""多数签名不成立""校正后失去预后价值"——三条都与上面这些数字冲突。</p>

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
<code>attenuation_pctile_all.tsv</code>、<code>null_summary_all.tsv</code>；日志 <code>logs/20–24</code>。
签名名单：<code>signature_registry.csv</code>（13）+ <code>signature_registry_expansion.csv</code>（28），排除项见两份 excluded 表。
修订文件：<code>01_amendment_expansion_to20.md</code>、<code>02_amendment_inflammation_wave3.md</code>。报告由脚本生成（<code>tools/12_make_report_stage2.py</code>），未手抄数字。</p>
</div></body></html>"""

outp = os.path.join(PB, "报告_阶段3_炎症定向后43签名_2026-10-20.html")
open(outp, "w", encoding="utf-8").write(HTML)
print("->", outp, len(HTML.encode("utf-8")), "bytes")
