#!/usr/bin/env python3
"""coloc combine.abf ported to Python — PP.H0-H4 from two per-variant log-ABF vectors.

Port of coloc's combine.abf / logsum / logdiff -- bit-exact vs R wherever R is finite, with
one deliberate, documented divergence on non-finite H3 terms (see `combine_abf`) (R/claudia.R,
github.com/chr1swallace/coloc @ master, fetched 2026-08-10), the computation behind
coloc.abf (Giambartolomei et al. 2014) once per-variant lABFs exist. coloc-susie
(Wallace 2021) runs the same combine.abf on per-signal SuSiE lbf vectors, so for
single-signal ABF inputs (our 1-causal benchmark) the two coincide.

Inputs here are the simulation ABF `lbf` columns (natural-log Bayes factors), restricted
to the variants non-missing on BOTH sides -- coloc.abf's internal merge keeps shared
SNPs only.

This module is the single implementation of the port. `scripts/analysis/coloc_similarity.py`
re-exports it so part 1's scripts and the validation harness keep working unchanged; there is
deliberately no second copy, because a numerically validated port that exists twice will drift.

Validated against coloc's R source (not a transcription of it) by
`PYTHONHASHSEED=0 python3 scripts/analysis/validate_coloc_port.py`, which downloads R/claudia.R, sources
logsum/logdiff/combine.abf verbatim, and compares 10 deterministic cases: worst
|py - R| = 3.6e-15 over all five posteriors. Per-case table:
results/simulations/similarity_algorithm_selection/coloc_port_validation.tsv.
Method write-up: .../similarity_algorithms.tex section 5.

Priors (coloc defaults): p1 = p2 = 1e-4, p12 = 1e-5.
"""
from __future__ import annotations

import numpy as np


def logsum(x: np.ndarray) -> float:
    """R coloc logsum: max(x) + log(sum(exp(x - max(x))))."""
    m = np.max(x)
    return float(m + np.log(np.sum(np.exp(x - m))))


def logdiff(x: float, y: float) -> float:
    """R coloc logdiff: max(x,y) + log(exp(x - max) - exp(y - max))."""
    m = max(x, y)
    d = np.exp(x - m) - np.exp(y - m)
    with np.errstate(divide="ignore"):
        return float(m + np.log(d))


def combine_abf(l1: np.ndarray, l2: np.ndarray,
                p1: float = 1e-4, p2: float = 1e-4, p12: float = 1e-5) -> np.ndarray:
    """PP.H0..PP.H4 for two aligned log-ABF vectors over the same variants.

    ONE deliberate divergence from R's combine.abf (2026-08-21 review): when `logdiff`'s
    argument is non-positive (S1*S2 <= S12 -- exactly so for a single shared variant, or by
    float rounding), R produces NaN/-Inf that PROPAGATES through the denominator, turning all
    five posteriors NaN. This port instead drops the non-finite H3 term from the denominator,
    i.e. treats "H3 impossible" as H3 mass 0 and returns well-defined posteriors. Everywhere
    both are finite the port is bit-exact vs R (validate_coloc_port.py, worst 3.6e-15).
    """
    lsum = l1 + l2
    ls1, ls2, ls12 = logsum(l1), logsum(l2), logsum(lsum)
    lh = np.array([
        0.0,
        np.log(p1) + ls1,
        np.log(p2) + ls2,
        np.log(p1) + np.log(p2) + logdiff(ls1 + ls2, ls12),
        np.log(p12) + ls12,
    ])
    # logdiff is -inf when S1*S2 == S12 exactly (single shared variant): H3 mass is 0
    denom = logsum(lh[np.isfinite(lh)])
    pp = np.exp(lh - denom)
    pp[~np.isfinite(pp)] = 0.0
    return pp
