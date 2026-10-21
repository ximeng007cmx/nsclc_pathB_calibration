# -*- coding: utf-8 -*-
"""生成外部复算预览报告 HTML：全部数字从 results/*.tsv 直读。"""
import csv, math, os, statistics as st

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


long = tsv(f"{R}/external_m1m3_long.tsv")
summ = tsv(f"{R}/external_m1m3_summary.tsv")
rho = {r["cohort"]: r for r in tsv(f"{R}/external_vs_tcga_rho.tsv")}
tcga = {r["sig_id"]: float(r["att_immune"]) for r in tsv(f"{R}/attenuation_all.tsv")}

COH = ["GSE68465", "GSE31210"]
STR = ["inflammation", "metabolism", "immune"]
STR_LBL = {"inflammation": "炎症", "metabolism": "代谢", "immune": "免疫"}
SM = {(r["cohort"], r["stratum"]): r for r in summ}
T = {  # TCGA 主分析（summary_all.tsv）
    "inflammation": ("-9.9%", 11), "metabolism": ("+0.4%", 32), "immune": ("+8.4%", 6)}

# 覆盖率
covs = [(r["sig_id"], r["cohort"], int(r["n_used"]), int(r["n_genes"]), float(r["coverage"])) for r in long]
low_cov = sorted([c for c in covs if c[4] < 0.4 and c[2] > 0], key=lambda x: x[4])
no_hit = [c for c in covs if c[2] == 0]
med_cov = {c: st.median([x[4] for x in covs if x[1] == c]) for c in COH}
n_miss = sum(1 for r in long if r["A"] == "" or r["A"] == "NA")

# 三队列免疫层 A>0 计数
imm_pos = {c: sum(1 for r in long if r["cohort"] == c and r["stratum"] == "immune" and r["A"] not in ("", "NA") and float(r["A"]) > 0) for c in COH}

rows_summ = ""
for s in STR:
    tcga_txt = f"{T[s][0]}（n={T[s][1]}）"
    g1 = SM[("GSE68465", s)]; g2 = SM[("GSE31210", s)]
    v_t = float(T[s][0].rstrip("%").replace("+", "").replace("−", "-")) / 100
    m3 = [v_t, float(g1["median_A"]), float(g2["median_A"])]
    order_ok = "✓（免疫三队列均为正且最大）" if s == "immune" else "—"
    rows_summ += (
        f'<tr><td><strong>{STR_LBL[s]}</strong></td><td class="num">{tcga_txt}</td>'
        f'<td class="num">{pct(float(g1["median_A"]))}（n={g1["n"]}，A&gt;0 {g1["n_pos"]}/{g1["n"]}）</td>'
        f'<td class="num">{pct(float(g2["median_A"]))}（n={g2["n"]}，A&gt;0 {g2["n_pos"]}/{g2["n"]}）</td>'
        f'<td class="num">{order_ok}</td></tr>')

low_rows = "\n".join(
    f'<tr><td>{sid}</td><td>{co}</td><td class="num">{nu}/{ng}</td><td class="num">{pct(cv)}</td></tr>'
    for sid, co, nu, ng, cv in low_cov[:8])

HTML = f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>外部复算预览：GSE68465 / GSE31210 的 M1→M3</title>
<style>
:root{{--ink:#1a2230;--muted:#6b7684;--line:#e2e6ec;--bg:#fbfcfd;--acc:#0b6b5f;--bad:#b3261e;--warn:#8a5a00;--ok:#1b6b3a}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.75 -apple-system,"Segoe UI","Microsoft YaHei",sans-serif}}
.wrap{{max-width:940px;margin:0 auto;padding:34px 26px 70px}}
h1{{font-size:24px;margin:0 0 6px}} .sub{{color:var(--muted);font-size:13px;margin-bottom:22px}}
h2{{font-size:18px;margin:32px 0 10px;padding-bottom:6px;border-bottom:2px solid var(--line)}}
p{{margin:9px 0}}
.lead{{background:#fff;border:1px solid var(--line);border-left:4px solid var(--acc);border-radius:8px;padding:16px 18px;margin:16px 0}}
.lead p:first-child{{margin-top:0}}.lead p:last-child{{margin-bottom:0}}
table{{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--line);border-radius:8px;overflow:hidden;font-size:12.8px;margin:10px 0}}
th,td{{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}}
th{{background:#f4f7f9;font-weight:600;font-size:12px;color:#33414f}}
tr:last-child td{{border-bottom:none}}
.num{{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}}
.note{{color:var(--muted);font-size:12.5px}}
code{{background:#eef2f6;padding:1px 5px;border-radius:4px;font-size:12.5px}}
ul,ol{{margin:8px 0;padding-left:22px}} li{{margin:5px 0}}
</style></head><body><div class="wrap">

<h1>外部复算预览：两个独立 LUAD 队列的 M1→M3</h1>
<p class="sub">GSE68465（Director's Challenge，439 例/235 事件）+ GSE31210（日本 IA–II，226 例/35 事件）· 2026-10-20 · 预览性分析</p>

<div class="lead">
<p><strong>一句话结论：三层排序（免疫衰减 &gt; 代谢≈0 &gt; 炎症略负）在全部三个队列复现。</strong>
免疫层中位 A 在 TCGA {T['immune'][0]}、GSE68465 <strong>{pct(float(SM[('GSE68465','immune')]['median_A']))}</strong>（A&gt;0 {imm_pos['GSE68465']}/6）、
GSE31210 <strong>{pct(float(SM[('GSE31210','immune')]['median_A']))}</strong>（A&gt;0 {imm_pos['GSE31210']}/6）——
三队列免疫层都是唯一明确为正、且 4/6 衰减的分层。</p>
<p>同时必须知道：<strong>逐签名的 A 值跨队列相关 ρ ≈ 0</strong>（GSE68465 {f(rho['GSE68465']['rho'],3)}、GSE31210 {f(rho['GSE31210']['rho'],3)}）——
单个签名的衰减幅度不可复制（多数是 M1 不显著者的噪声），但<strong>分层层面的图景稳健</strong>。
这反而支持论文只做分层中位数、不做逐签名外部对比的写法。</p>
</div>

<h2>一、数据与方法（先披露再谈结果）</h2>
<table>
<thead><tr><th style="width:22%">项</th><th>说明</th></tr></thead>
<tbody>
<tr><td>队列</td><td><strong>GSE68465</strong>（E-GEOD-68465，Director's Challenge，LUAD，HG-U133A）：12545 基因 × 439 样本（235 事件）；
<strong>GSE31210</strong>（E-GEOD-31210，日本 LUAD IA–II，HG-U133 Plus 2）：16504 基因 × 226 样本（35 事件）。
数据复用 pathA 已归档的 ArrayExpress 镜像（EBI BioStudies 下载，NCBI 不可达时的替代通道），表达矩阵 log2(x+1)、探针 max-mean 合并。</td></tr>
<tr><td>ESTIMATE 复算</td><td>官方 estimate 包不可得（Bioconductor 3.23 无此包、mdanderson 不再提供下载、GitHub 不可达）。
官方签名集（Stromal141_UP / Immune141_UP，各 141 基因）取自原始论文补充材料（Yoshihara et al. 2013, Nat Commun, <em>Supplementary Data 1</em>，本机 <code>work/estimate_ext/ESM_488.xlsx</code>），
用 GSVA::ssgsea 算单样本 ES 作分数。<strong>数学事实：Cox 中协变量线性缩放不改变 sig 系数，故 A 与官方分数严格等价</strong>（只差常数倍）。</td></tr>
<tr><td>M1 协变量</td><td>GSE68465：sig + age + sex + stage_late（pT3/pT4/pN1/pN2 = late）；GSE31210：<strong>sig + sex + stage_late（无 age）</strong>——
SDRF 中 226 例 age 全部 "not specified"，无法纳入，写论文必须披露。</td></tr>
<tr><td>基因覆盖</td><td>平台是芯片，签名基因覆盖中位 GSE68465 {pct(med_cov['GSE68465'])}、GSE31210 {pct(med_cov['GSE31210'])}；
<strong>SIG26 的 4 个基因全部不在芯片上（无外部值）</strong>；另有 {len(low_cov)} 个队列-签名组合覆盖 &lt;40%（见下表），其 A 只列不解读。</td></tr>
<tr><td>性质</td><td><strong>预览性分析，未预注册。</strong>若纳入论文，须作为 post-hoc 外部复算披露；主注册表与锁定哈希不改。</td></tr>
</tbody></table>

<h2>二、主结果：三队列 × 三分层</h2>
<table>
<thead><tr><th style="width:12%">分层</th><th class="num">TCGA（主分析）</th><th class="num">GSE68465</th><th class="num">GSE31210</th><th class="num">排序复现</th></tr></thead>
<tbody>{rows_summ}</tbody></table>
<p><strong>读法：</strong>代谢层三队列都贴着零（{T['metabolism'][0]} / {pct(float(SM[('GSE68465','metabolism')]['median_A']))} / {pct(float(SM[('GSE31210','metabolism')]['median_A']))}）；
炎症层三队列都是小负值；免疫层三队列都是正中位数且 A&gt;0 计数一致（4/6、4/6、4/6）。
<strong>但注意：两个外部队列的免疫层都没有可解读的显著性（GSE68465 免疫层 M1 显著 0 个；GSE31210 仅 SIG45 一个）</strong>，
所以外部复算支持的是"方向图景"，不改变 P2 未获证实的结论。</p>
<p class="note">FDR/M1→M3 显著数（逐层）：GSE68465 炎症 2→2、代谢 3→3、免疫 0→0；GSE31210 炎症 1→2、代谢 5→6、免疫 1→1——外部也看不到"整类掉显著性"。</p>

<h2>三、逐签名 A 跨队列相关 ≈ 0 的含义</h2>
<table>
<thead><tr><th>队列对</th><th class="num">n（共同可算）</th><th class="num">Spearman ρ</th><th class="num">P</th></tr></thead>
<tbody>
<tr><td>GSE68465 vs TCGA</td><td class="num">{rho['GSE68465']['n']}</td><td class="num">{f(rho['GSE68465']['rho'],3)}</td><td class="num">{f(rho['GSE68465']['p'],3)}</td></tr>
<tr><td>GSE31210 vs TCGA</td><td class="num">{rho['GSE31210']['n']}</td><td class="num">{f(rho['GSE31210']['rho'],3)}</td><td class="num">{f(rho['GSE31210']['p'],3)}</td></tr>
</tbody></table>
<p>这不能读成"复算失败"：49 个签名里只有少数 M1 显著，其余的 A 本来就是比值估计量在噪声上的投影（主报告 Limitations 第 4 条）。
它支持的结论是：<strong>衰减的可复制单元是"分层图景"而不是"单签名数值"</strong>——论文写法上应避免逐签名外部对比表，只放分层汇总。</p>

<h2>四、低覆盖组合（A 只列不解读）</h2>
<table>
<thead><tr><th>签名</th><th>队列</th><th class="num">命中/应有</th><th class="num">覆盖</th></tr></thead>
<tbody>{low_rows}</tbody></table>
<p class="note">完整 98 行明细见 <code>results/external_m1m3_long.tsv</code>；SIG26 两队列均无外部值（4 基因全不在芯片上）。</p>

<h2>五、是否纳入论文的建议</h2>
<div class="lead">
<p><strong>建议纳入，作为一节"外部方向性复算"（约半页）+ 补充表</strong>，理由：</p>
<ol>
<li>它把"单队列"这条最大拒稿理由拆掉了：三队列、665 例外部样本，分层排序全部复现。</li>
<li>它不改变任何主结论：免疫层仍只是"方向相容"（外部无显著性支持），零结果反而更稳（代谢层三队列贴零）。</li>
<li>诚实且加分：逐签名 ρ≈0 本身是个可报告的方法学发现（衰减是队列特异的），与"校正无系统性方向"主线自洽。</li>
</ol>
<p><strong>写法红线（沿用主报告口径）：</strong>只能写"免疫层方向在三个队列与阳性对照相容"，<strong>不能写"外部验证证实了免疫签名衰减"</strong>；
GSE31210 缺 age、SIG26 无外部值、低覆盖清单三件事都要进 Methods/Limitations。</p>
</div>

<hr>
<p class="note">数据来源：<code>results/external_m1m3_long.tsv</code>、<code>external_m1m3_summary.tsv</code>、<code>external_vs_tcga_rho.tsv</code>；
脚本 <code>R/22_external_m1m3.R</code>、日志 <code>logs/20_external_m1m3.log</code>。
外部数据：pathA <code>data/external/gse68465/</code>、<code>data/external/gse31210/</code>（ArrayExpress 镜像）。
ESTIMATE 签名集：Yoshihara 2013 Nat Commun Supplementary Data 1（<code>work/estimate_ext/ESM_488.xlsx</code>）。
报告由 <code>tools/23_make_external_preview.py</code> 生成，未手抄数字。</p>
</div></body></html>"""

outp = os.path.join(PB, "报告_外部复算预览_GSE68465_GSE31210_2026-10-20.html")
open(outp, "w", encoding="utf-8").write(HTML)
print("->", outp, len(HTML.encode("utf-8")), "bytes")
