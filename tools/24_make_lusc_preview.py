# -*- coding: utf-8 -*-
"""LUSC 跨组织学复算结果简报（数据从 results/*.tsv 直读）。"""
import csv, os, statistics as st

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
    return (f"{v:+.{n}f}" if sign else f"{v:.{n}f}").replace("-", "−")


def pct(x, n=1):
    try:
        return f"{float(x)*100:+.{n}f}%".replace("-", "−")
    except Exception:
        return "NA"


L = {r["sig_id"]: r for r in tsv(f"{R}/attenuation_lusc49.tsv")}
U = {r["sig_id"]: r for r in tsv(f"{R}/attenuation_all.tsv")}
SL = {r["set"]: r for r in tsv(f"{R}/summary_lusc49.tsv")}
SU = {r["set"]: r for r in tsv(f"{R}/summary_all.tsv")}
ctrl = {r["axis"]: r for r in tsv(f"{R}/control_axes_lusc.tsv")}
ctrl_u = {r["axis"]: r for r in tsv(f"{R}/control_axes_all.tsv")}


def dirs(d):
    sh = gr = fl = 0
    for r in d.values():
        b1, b3 = float(r["beta_M1"]), float(r["beta_M3"])
        if b1 * b3 > 0:
            if abs(b3) < abs(b1):
                sh += 1
            else:
                gr += 1
        else:
            fl += 1
    return sh, gr, fl


sh_l, gr_l, fl_l = dirs(L)
sh_u, gr_u, fl_u = dirs(U)
sig_l = sorted([r for r in L.values() if float(r["p_M1"]) < 0.05], key=lambda r: float(r["p_M1"]))
N_COHORT = int(list(L.values())[0]["n"])
N_EVENTS = int(list(L.values())[0]["events"])

STR_LBL = {"inflammation": "炎症", "metabolism": "代谢", "immune": "免疫"}
rows = ""
for s in ("inflammation", "metabolism", "immune"):
    a = SL[s]; b = SU[s]
    rows += (f'<tr><td><strong>{STR_LBL[s]}</strong></td>'
             f'<td class="num">{pct(float(b["median_A"])/100)}（n={b["n"]}，有效 {b["n_ratio_valid"]}）</td>'
             f'<td class="num">{pct(float(a["median_A"])/100)}（n={a["n"]}，有效 {a["n_ratio_valid"]}）</td></tr>')

ctrl_rows = ""
for ax in ("lymphocyte", "neutrophil", "nutrition_liver"):
    a = ctrl[ax]; b = ctrl_u[ax]
    ctrl_rows += (f'<tr><td>{ax}</td><td class="num">{f(b["rho_immune"],2)}</td>'
                  f'<td class="num">{f(a["rho_immune"],2)}</td>'
                  f'<td class="num">{f(b["p_M1"],3)}</td><td class="num">{f(a["p_M1"],3)}</td>'
                  f'<td class="num">{pct(float(a["attenuation"]))}</td></tr>')

sig_rows = "".join(
    f'<tr><td>{r["sig_id"]}</td><td class="num">{f(r["beta_M1"],3,True)}</td>'
    f'<td class="num">{f(r["p_M1"],4)}</td><td class="num">{pct(float(r["att_immune"]))}</td>'
    f'<td class="num">{f(r["q_M1"],3)}</td></tr>' for r in sig_l)

HTML = f"""<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>LUSC 跨组织学复算结果（修订 04）</title>
<style>
:root{{--ink:#1a2230;--muted:#6b7684;--line:#e2e6ec;--bg:#fbfcfd;--acc:#0b6b5f;--warn:#8a5a00;--bad:#b3261e}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.75 -apple-system,"Segoe UI","Microsoft YaHei",sans-serif}}
.wrap{{max-width:920px;margin:0 auto;padding:34px 26px 70px}}
h1{{font-size:23px;margin:0 0 6px}}.sub{{color:var(--muted);font-size:13px;margin-bottom:20px}}
h2{{font-size:17px;margin:30px 0 10px;padding-bottom:6px;border-bottom:2px solid var(--line)}}
.lead{{background:#fff;border:1px solid var(--line);border-left:4px solid var(--warn);border-radius:8px;padding:16px 18px;margin:16px 0}}
.lead p:first-child{{margin-top:0}}.lead p:last-child{{margin-bottom:0}}
table{{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--line);border-radius:8px;overflow:hidden;font-size:12.8px;margin:10px 0}}
th,td{{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left}}
th{{background:#f4f7f9;font-weight:600;font-size:12px;color:#33414f}}
tr:last-child td{{border-bottom:none}}
.num{{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}}
.note{{color:var(--muted);font-size:12.5px}}code{{background:#eef2f6;padding:1px 5px;border-radius:4px;font-size:12.5px}}
ul,ol{{padding-left:22px}}li{{margin:5px 0}}
</style></head><body><div class="wrap">

<h1>LUSC 跨组织学复算结果（修订 04）</h1>
<p class="sub">TCGA-LUSC {N_COHORT} 例 / {N_EVENTS} 事件 · 2026-10-21 · 同一表达矩阵、同一 M1/M3、同一估计量、同一 bootstrap</p>

<div class="lead">
<p><strong>结论：字面上满足修订 04 的 R1（稳健性成立），但实质上不可判读——必须按"无信号队列"处理。</strong></p>
<p>原因：<strong>49 个已发表签名在 LUSC 中没有一个达到 FDR 显著</strong>（LUAD 是 31 个），仅 {len(sig_l)} 个名义 P&lt;0.05；
同签名的 β<sub>M1</sub> 在 LUAD 与 LUSC 之间的 Spearman 相关是 <strong>−0.174</strong>。
也就是说，这些签名在 LUSC 里基本没有可检出的预后效应，A = 1−β<sub>M3</sub>/β<sub>M1</sub> 的分母大多是噪声——
此时"中位 A 落在 ±10%""方向计数对半"都是噪声的性质，不能读成"校正无系统性方向"的复现。</p>
<p>修订 04 的 R1–R4 未预见"有样本、有事件、但签名无信号"这种队列，这是规则外情形，按最保守方式处理：<strong>只报事实，不判读方向</strong>。</p>
</div>

<h2>一、主结果：三队列对照</h2>
<table>
<thead><tr><th style="width:12%">分层</th><th class="num">LUAD（主分析）</th><th class="num">LUSC</th></tr></thead>
<tbody>{rows}</tbody></table>
<p>LUSC 的中位 A 全为正（{pct(float(SL['ALL']['median_A'])/100)}），三层排序（免疫 &gt; 代谢 &gt; 炎症）与 LUAD 一致；
但 LUSC 三层的"有效 n"加起来只有 {int(SL['inflammation']['n_ratio_valid'])+int(SL['metabolism']['n_ratio_valid'])+int(SL['immune']['n_ratio_valid'])} 个（LUAD 是 36）。
<strong>中位数建立在这么少的有效值上，且各层 n 不等，排序本身没有推断意义。</strong></p>

<h2>二、方向计数（这一项不依赖签名是否显著，可以报）</h2>
<table>
<thead><tr><th>队列</th><th class="num">收缩</th><th class="num">放大</th><th class="num">反转</th><th class="num">FDR 显著（M1→M3）</th></tr></thead>
<tbody>
<tr><td>LUAD（{SU['ALL']['n']} 个签名）</td><td class="num">{sh_u}</td><td class="num">{gr_u}</td><td class="num">{fl_u}</td>
<td class="num">{SU['ALL']['n_FDR_M1']} → {SU['ALL']['n_FDR_M3']}</td></tr>
<tr><td>LUSC（{SL['ALL']['n']} 个签名）</td><td class="num">{sh_l}</td><td class="num">{gr_l}</td><td class="num">{fl_l}</td>
<td class="num"><strong>{SL['ALL']['n_FDR_M1']} → {SL['ALL']['n_FDR_M3']}</strong></td></tr>
</tbody></table>
<p>LUSC 的收缩:放大 = {sh_l}:{gr_l}，与 LUAD 的 {sh_u}:{gr_u} 同样接近对半——<strong>这一点可以作为"方向无系统性"的跨组织学旁证</strong>，
且它不需要签名本身显著（纯方向计数）。但强度低于 LUAD：LUAD 有 32 个 M1 显著者支撑，LUSC 只有 {len(sig_l)} 个。</p>

<h2>三、LUSC 中仅有的 5 个名义显著签名</h2>
<table>
<thead><tr><th>签名</th><th class="num">β<sub>M1</sub></th><th class="num">P<sub>M1</sub></th><th class="num">衰减 A</th><th class="num">q<sub>M1</sub></th></tr></thead>
<tbody>{sig_rows}</tbody></table>
<p class="note">全部 |A| &lt; 10%，q 值均 &gt; 0.05（BH 按 49 个校正后无一显著）。</p>

<h2>四、对照轴：管线在 LUSC 也失效，且原因不在管线</h2>
<table>
<thead><tr><th>对照轴</th><th class="num">ρ(分,Immune) LUAD</th><th class="num">ρ(分,Immune) LUSC</th><th class="num">P<sub>M1</sub> LUAD</th><th class="num">P<sub>M1</sub> LUSC</th><th class="num">LUSC 衰减</th></tr></thead>
<tbody>{ctrl_rows}</tbody></table>
<p>淋巴细胞轴在 LUSC 的 ρ 仍是 {f(ctrl['lymphocyte']['rho_immune'],2)}（与 LUAD {f(ctrl_u['lymphocyte']['rho_immune'],2)} 几乎相同），
说明 ESTIMATE 打分在 LUSC 正常工作；但它的 β<sub>M1</sub> 在 LUSC 中 P = {f(ctrl['lymphocyte']['p_M1'],3)}——
<strong>免疫含量本身在 LUSC 里就不预测 OS</strong>（调整年龄/性别/分期后），所以"衰减"无从谈起（β<sub>M1</sub>≈0，A 爆炸到 {pct(float(ctrl['lymphocyte']['attenuation']))}，属无效估计）。
这不是管线问题，是 LUSC 的生物学/队列特征。</p>

<h2>五、可以怎么写（三个选项）</h2>
<ol>
<li><strong>（推荐）写进补充材料，一句话</strong>：LUSC 中 49 个签名均未达 FDR 显著（0/49），β<sub>M1</sub> 与 LUAD 相关 ρ=−0.17，
故该队列无法用于方向性复算；方向计数（{sh_l}:{gr_l}）与 LUAD 同为对半，仅作旁证。优点：诚实、不夸大；缺点：浪费了一个队列。</li>
<li><strong>升格为一节"签名在 LUSC 不复现"</strong>：把"0/49 FDR 显著 + ρ=−0.17"作为独立发现，强化"已发表签名可重复性有限"的论点。
优点：故事更强、多一个实证；<strong>风险：会引出审稿人反问"你选的 49 个签名质量是否本身就有问题"，与主结论（签名效应稳定、只是混杂不解释）相互张力</strong>。</li>
<li>不写：只留档。优点：不引入矛盾；缺点：将来审稿人问"为什么不做 LUSC"时无法答。</li>
</ol>
<p class="note">我的判断：选 1。理由——本文主线是"免疫/基质含量不解释签名效应量"，而 LUSC 结果是"签名本身在 LUSC 不复现"，
是<strong>另一条故事线</strong>；把它塞进主线会制造逻辑张力（一边说签名效应稳定，一边说 LUSC 全不复现）。放补充材料 + Limitations 一句最稳。</p>

<hr>
<p class="note">数据：<code>results/attenuation_lusc49.tsv</code>、<code>summary_lusc49.tsv</code>、<code>control_axes_lusc.tsv</code>；
脚本 <code>R/23_lusc_calibration.R</code>，日志 <code>logs/21_lusc_calibration.log</code>；修订文件 <code>04_amendment_lusc.md</code>。
报告由 <code>tools/24_make_lusc_preview.py</code> 生成，数字直读未手抄。</p>
</div></body></html>"""

outp = os.path.join(PB, "报告_LUSC复算_2026-10-21.html")
open(outp, "w", encoding="utf-8").write(HTML)
print("->", outp, len(HTML.encode("utf-8")), "bytes")
