#!/usr/bin/env python3
"""coloc combine.abf ported to Python — PP.H0-H4 from two per-variant log-ABF vectors.

Port of coloc's combine.abf / logsum / logdiff. The computation computes 
PP.H0-H4 given per-variant natural-log Bayes factors (lABFs), restricted 
to the variants non-missing on BOTH sides.

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

    Note on divergence from R's combine.abf: when `logdiff`'s argument is 
    non-positive (S1*S2 <= S12 -- exactly so for a single shared variant, or by
    float rounding), R produces NaN/-Inf that propagates through the denominator, 
    turning all five posteriors NaN. This implementation instead drops the 
    non-finite H3 term from the denominator (treats "H3 impossible" as H3 mass 0) 
    and returns well-defined posteriors.
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
