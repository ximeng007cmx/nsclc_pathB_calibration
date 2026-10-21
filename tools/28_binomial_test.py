# -*- coding: utf-8 -*-
"""
28_binomial_test.py — C1 fix (pre-submission review 2026-10-21)

Manuscript MS_draft_v1.md §3.2 previously stated "a binomial test gives P ≈ 1"
for the direction split 17 shrinking vs 15 growing among the 32 signatures with
a nominally significant unadjusted coefficient. That number had no
computational record and is numerically wrong.

This script computes the exact two-sided binomial test (R binom.test default
method: sum of probabilities <= P(observed)) and writes the result to
results/binomial_direction.tsv. Every downstream number is asserted.

Input counts are hard-coded FROM the pipeline output and cross-checked against
results/manuscript_reconciliation.tsv (item: direction split 17:15, n=32).
"""
import csv
import hashlib
import io
import math
import os

# --- inputs ( asserted against reconciliation output below ) -----------------
N_SIG = 32          # signatures with nominally significant M1 coefficient
K_SHRINK = 17       # attenuated (A > 0) among the 32
K_GROW = 15         # grew (A < 0) among the 32
P_NULL = 0.5

assert K_SHRINK + K_GROW == N_SIG
assert N_SIG % 2 == 0  # two-sided symmetric around N/2

# cross-check against the reconciliation table: §3.2 rows record 17 and 15
RECON = os.path.join(os.path.dirname(__file__), "..", "results",
                     "manuscript_reconciliation.tsv")
vals = {}
if os.path.exists(RECON):
    with io.open(RECON, encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            item = str(row.get("item", ""))
            if str(row.get("section", "")) == "3.2" and "显著者中" in item:
                vals[item] = row.get("claimed")
assert any(v == "17" for v in vals.values()), "shrink=17 row not found: %s" % vals
assert any(v == "15" for v in vals.values()), "grow=15 row not found: %s" % vals


def binom_pmf(k, n, p):
    return comb(n, k) * p ** k * (1 - p) ** (n - k) if False else \
        math.comb(n, k) * p ** k * (1 - p) ** (n - k)


def binom_test_two_sided(k, n, p=P_NULL):
    """R binom.test default: two-sided = sum of pmf values <= pmf(k)."""
    pk = binom_pmf(k, n, p)
    tail = sum(binom_pmf(i, n, p) for i in range(n + 1) if binom_pmf(i, n, p) <= pk + 1e-15)
    return min(1.0, tail)


p_exact = binom_test_two_sided(K_SHRINK, N_SIG)
# anchor value computed independently: 1 - P(X=16) for n=32, p=0.5
anchor = 1.0 - math.comb(32, 16) / 2 ** 32
assert abs(p_exact - anchor) < 1e-12, (p_exact, anchor)
assert abs(p_exact - 0.8601) < 5e-5, p_exact   # 0.86008...

out_path = os.path.join(os.path.dirname(__file__), "..", "results",
                        "binomial_direction.tsv")
rows = [
    ["comparison", "n", "k_shrink", "k_grow", "test", "method", "p_value"],
    ["direction split among 32 nominally significant signatures (M1 P<0.05), "
     "TCGA-LUAD, M1 vs M3 attenuation sign",
     N_SIG, K_SHRINK, K_GROW,
     "exact binomial two-sided vs p=0.5",
     "R binom.test convention (sum of pmf <= pmf(observed))",
     "%.4f" % p_exact],
]
with io.open(out_path, "w", encoding="utf-8", newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows(rows)

print("exact two-sided binomial P (17 vs 15, n=32) = %.6f" % p_exact)
print("anchor (1 - P(X=16))                        = %.6f" % anchor)
print("written:", os.path.normpath(out_path))

h = hashlib.sha256(io.open(out_path, "rb").read()).hexdigest()[:16]
print("tsv sha256[:16] =", h)
