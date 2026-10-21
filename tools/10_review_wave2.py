# -*- coding: utf-8 -*-
"""扩容波次人工复核辅助：把抽取结果按"需人工判断的可能"分档打印，附证据窗口。

分档（决定是否需要人工读原文）：
  A_OK        stated == extracted，且基因数 2–15            -> 可直接进入名单（仍抽检）
  B_NOSTATED  无自述基因数但抽到列表                       -> 需看证据句是否真是"签名构成"
  C_MISMATCH  自述数与抽到数不符                           -> 必须人工
  D_FAIL      NO_LIST / FETCH_FAIL / I6 越界               -> 记录排除理由
"""
import json, re, sys

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
recs = json.load(open(f"{PB}/work/extraction_wave2.json", encoding="utf-8"))

def tier(x):
    st = x.get("status")
    if st == "OK" and x.get("i6_ok"):
        return "A_OK"
    if st == "NO_STATED_N" and x.get("i6_ok"):
        return "B_NOSTATED"
    if st == "MISMATCH":
        return "C_MISMATCH"
    return "D_FAIL"

buckets = {}
for x in recs:
    buckets.setdefault(tier(x), []).append(x)

for k in ("A_OK", "B_NOSTATED", "C_MISMATCH", "D_FAIL"):
    print(f"\n{'='*100}\n{k}  n={len(buckets.get(k, []))}\n{'='*100}")
    for x in buckets.get(k, []):
        ev = (x.get("evidence") or "").replace("\n", " ")
        print(f"\n{x['pmcid']}  PMID {x.get('pmid')}  {x.get('year')}  {x.get('journal')}")
        print(f"  {x.get('title','')[:150]}")
        print(f"  rule={x.get('rule')}  stated={x.get('n_stated')}  got={x.get('n_extracted')}  "
              f"status={x.get('status')}  i6={x.get('i6_ok')}")
        print(f"  GENES: {x.get('genes','')}")
        print(f"  stated-evidence: {str(x.get('n_stated_evidence'))[:120]}")
        print(f"  list-evidence  : {ev[:300]}")
