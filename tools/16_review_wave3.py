# -*- coding: utf-8 -*-
"""炎症波次复核辅助：只打印"INF 通过且列表可抽"的条目，附基因池来源证据句。"""
import json

PB = r"D:/WorkBuddy_Y02/out_ali_directions/pathB"
recs = json.load(open(f"{PB}/work/extraction_wave3.json", encoding="utf-8"))
ok = [x for x in recs if x.get("inf_pass") and x.get("status") in ("OK", "NO_STATED_N", "MISMATCH")
      and x.get("i6_ok")]
bad = [x for x in recs if not (x.get("inf_pass") and x.get("i6_ok"))]
print(f"INF 通过 + 可抽：{len(ok)}；其余：{len(bad)}\n")
for x in sorted(ok, key=lambda y: str(y.get("year"))):
    print(f"\n{x['pmcid']}  PMID {x.get('pmid')}  {x.get('year')}  {x.get('journal')}")
    print(f"  {x.get('title','')[:140]}")
    print(f"  stated={x.get('n_stated')} got={x.get('n_extracted')} status={x.get('status')}")
    print(f"  GENES: {x.get('genes','')}")
    print(f"  INF 证据: {(x.get('inf_evidence') or '')[:280]}")
print("\n\n===== INF 未通过（供核对是否漏判）=====")
for x in bad:
    print(f"  {x['pmcid']} {x.get('year')} {x.get('status'):12s} inf={x.get('inf_pass')} "
          f"n={x.get('n_extracted')} {x.get('title','')[:72]}")
