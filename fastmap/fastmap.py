"""FastMap: combine per-cohort SuSiE fine-mapping results across cohorts into a single
set of posterior inclusion probabilities (PIPs) per variant, without joint modeling.

Each cohort contributes one or more SuSiE "single effect" components (alpha columns) per
genomic region. FastMap greedily merges components across cohorts that colocalise on one
pairwise similarity score (`PP.H4` from coloc-susie, or weighted Jaccard), then reports one
final PIP per variant from the L strongest resulting components by evidence (`component_f1`).

The original coupled "production" algorithm (PIP-coverage merge criterion with Jaccard gate
and escalating-`div` relaxation, plus a plain-Jaccard dedup at selection) was REMOVED on
2026-08-21 (Ran): every forward analysis uses the single-score paths, and keeping a frozen
second algorithm invited drift. Reproducing pre-removal outputs requires the snapshot
`scripts/fastmap_algorithm/fastmap_v08212026_pre_production_removal.py` (the last state
with both algorithms; the repo is not under git).

This module contains only the combination algorithm. Reading per-cohort SuSiE output
(.rds files) or real-data summary-stat feathers into the `region_df`/`pips_df` inputs
below is data-source-specific plumbing that belongs in a separate script.
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pandas as pd


def _ustack_idxmax(df):
    """Equivalent of `df.infer_objects(copy=False).fillna(-2).unstack().idxmax()`:
    the `(column_name, index_name)` of the FIRST maximum in column-major order.

    Tie-breaking is the load-bearing detail. `df.unstack()` yields a Series ordered
    column-major, so `idxmax` returns the first occurrence *in that order*; flattening the
    TRANSPOSE in C order reproduces exactly that traversal, and `argmax` likewise returns the
    first occurrence. Getting this wrong would not raise -- it would silently let a different
    pair win the argmax, which then cascades through every later merge, the same failure mode
    as the `PYTHONHASHSEED` non-determinism that made the stored 2025 grid one random draw.
    """
    a = df.to_numpy(dtype=float)
    a = np.where(np.isnan(a), -2.0, a)
    n_rows = a.shape[0]
    k = int(a.T.reshape(-1).argmax())      # C-order over transpose == column-major over a
    ci, ri = divmod(k, n_rows)
    return df.columns[ci], df.index[ri]


def _cohort_of(column_name) -> str:
    """The cohort a `{cohort}_alpha{i}` column (or single component name) belongs to:
    everything before the final underscore-separated token.

    THE ONLY PLACE this convention is parsed. It is load-bearing for A1 (a merged component
    must draw on distinct cohorts) and for the backfill's cohort lookup, and it used to be
    inlined at four call sites -- a convention shift fixed in one place would have silently
    NaN'ed `lbf_fastmap` and degraded the F1 tie-break (2026-08-22 review finding #10)."""
    return "_".join(str(column_name).split("_")[:-1])


def _cohort_names(component_name: str) -> set:
    """The set of cohort names referenced by a component name, stripping each name's
    trailing alpha-index suffix. A component produced by combining earlier components has
    a comma-joined name, e.g. "EUR_sim9_..._alpha1,AFR_sim1_..._alpha2"."""
    return {_cohort_of(name) for name in component_name.split(",")}


def shares_cohort(name_a: str, name_b: str) -> bool:
    """True if two components share an underlying cohort. Used to avoid combining two
    alpha columns that came from the same cohort's SuSiE run."""
    return len(_cohort_names(name_a) & _cohort_names(name_b)) > 0


def weighted_jaccard_matrix(df, reverse_colname_map, top_n=500):
    """Pairwise weighted Jaccard `W[a,b] = sum_v min(p_av, p_bv) / sum_v max(p_av, p_bv)`,
    computed on top-N-masked columns with missing values treated as 0.

    This is part 1's winning similarity (`v08012025_weighted_jaccard` in
    `scripts/analysis/similarity_metrics.py`), reproduced here so FastMap can use it as a
    merge criterion. Two fidelity points carried over deliberately from that benchmark:

      - the top-N mask is `~col.isin(col.nlargest(top_n))`, which is VALUE-based, not a rank
        cutoff, so ties can admit more than `top_n` entries;
      - `mx == 0` (both columns all-zero over the masked set) yields 0.0, not NaN.

    Symmetric by construction. Entries for pairs sharing a cohort, and the diagonal, are left
    NaN so the caller's argmax can never select them.
    """
    names = list(df.columns)
    masked = df.mask(df.apply(lambda col: ~col.isin(col.nlargest(int(top_n))), axis=0))
    d = masked.fillna(0).to_numpy(dtype=float)

    out = pd.DataFrame(np.nan, index=names, columns=names, dtype=float)
    for i, row in enumerate(names):
        a = d[:, i][:, None]
        mn = np.minimum(a, d).sum(axis=0)
        mx = np.maximum(a, d).sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            vals = np.where(mx > 0, mn / mx, 0.0)
        for j, col in enumerate(names):
            if j <= i:
                continue
            if shares_cohort(reverse_colname_map[row], reverse_colname_map[col]):
                continue
            out.loc[row, col] = out.loc[col, row] = vals[j]
    return out


def coloc_susie_matrix(lbf_df, reverse_colname_map, p1=1e-4, p2=1e-4, p12=1e-5,
                        overlap_min=0.5, trim_by_posterior=True):
    """Pairwise coloc-susie `PP.H4` between components, from their per-variant log-Bayes-factor
    vectors (SuSiE's ``lbf_variable`` row for each single-effect component).

    This is coloc-susie's all-pairs-of-credible-sets computation (Wallace 2021) evaluated on
    FastMap's own component-pair granularity, so no reduction step is needed. Each pair is
    scored by `coloc.combine_abf` over the variants **measured on both sides**, which is what
    ``coloc.abf``'s internal merge does. One consequence recorded in part 1 and carried here:
    a pair sharing a single variant makes H3 impossible (``logdiff`` -> -inf, handled inside
    `combine_abf`) -- which inflates `PP.H4` exactly when the overlap is thinnest, the case the
    trim below exists to kill.

    ``trim_by_posterior`` / ``overlap_min`` (added 2026-08-21, Ran) port ``coloc.bf_bf``'s
    identically-named gate at coloc's own defaults (TRUE, 0.5): a pair is admissible only if the
    shared variants capture at least ``overlap_min`` of EACH side's posterior signal mass.
    coloc computes that mass via ``logbf_to_pp`` with the scalar per-SNP prior and takes
    ``rowSums(pp[, isnps]) / rowSums(pp[, all non-null snps])``; for a scalar prior the prior
    cancels in that ratio (and the null column is excluded from both sums), so it reduces
    EXACTLY to ``prop = exp(logsum(lbf[shared]) - logsum(lbf[measured]))`` -- an algebraic
    identity, not an approximation. Per-SNP prior VECTORS (where nothing would cancel) are not
    supported here; we only ever use scalars. An inadmissible pair is left NaN, like a
    same-cohort pair: the caller's argmax can never select it, but nothing is remembered --
    the matrix is recomputed after every accepted merge, so a trimmed pair is re-evaluated
    each round and can become admissible later (a merged component's measured set is its
    parents' union, which can only grow the overlap). A pair sharing NO variants is prop=0 and
    NaN under the trim; with ``trim_by_posterior=False`` it scores 0.0 as before (that flag
    reproduces the pre-2026-08-21 MATRIX bit-exactly; the stored part-2 sweep TABLES also
    predate the same-day removal of the C1 lbf fill, so reproducing them end-to-end requires
    the pre-removal code state -- see `merge_lbf`).

    Note the missing-variant fill convention does **not** enter the pair SCORE -- restricting
    to shared variants makes a first-round pair's score independent of it, whichever convention
    is in force. It matters only for the merged log-BF vectors later rounds consume; see
    `merge_lbf` and ``documentation/coloc_susie_merged_lbf.pdf`` §4. The trim's denominators
    likewise read only MEASURED entries, never fills.

    Symmetric; same-cohort pairs and the diagonal are left NaN so the caller's argmax can never
    select them.
    """
    from fastmap.coloc import combine_abf, logsum

    names = list(lbf_df.columns)
    out = pd.DataFrame(np.nan, index=names, columns=names, dtype=float)
    # trim denominators: each component's total posterior signal mass, over what it measured
    log_total = {c: logsum(lbf_df[c].dropna().to_numpy(dtype=float)) for c in names}
    for i, row in enumerate(names):
        for j, col in enumerate(names):
            if j <= i:
                continue
            if shares_cohort(reverse_colname_map[row], reverse_colname_map[col]):
                continue
            la, lb = lbf_df[row], lbf_df[col]
            both = la.notna() & lb.notna()
            if not both.any():
                if not trim_by_posterior:
                    out.loc[row, col] = out.loc[col, row] = 0.0
                continue
            a = la[both].to_numpy(dtype=float)
            b = lb[both].to_numpy(dtype=float)
            if trim_by_posterior:
                prop_row = np.exp(logsum(a) - log_total[row])
                prop_col = np.exp(logsum(b) - log_total[col])
                if prop_row < overlap_min or prop_col < overlap_min:
                    continue
            pp = combine_abf(a, b, p1=p1, p2=p2, p12=p12)
            out.loc[row, col] = out.loc[col, row] = pp[4]
    return out


def merge_lbf(lbf_df: pd.DataFrame) -> pd.Series:
    """Merged per-variant log-Bayes-factor for a set of components: ``sum of lbf_k``.

    The SUM is derived, not chosen: conditional on the shared causal variant, disjoint samples
    factorise the joint Bayes factor, so log-BFs add -- the same summand coloc's own ``S12``
    term integrates (``documentation/coloc_susie_merged_lbf.pdf`` §3). It requires every
    component in the set to come from a distinct cohort (assumption A1), which the caller
    asserts. **A variant measured on every side is unaffected by anything below.**

    A cohort's MISSING entry (§4's case (ii)) takes ``lbf := log(S_k / n_k)``, the log of the
    **arithmetic mean Bayes factor** over the variants column *k* did measure -- "as likely
    causal as a typical measured variant in this cohort" (the C2 convention; the sole rule
    since 2026-08-21, see below). A variant measured by NO component in the set stays NaN, so
    the support is always the union of what was actually measured.

    This fill is the same value `combine_pips` imputes at the same merge: its coverage scaling
    plus the ``1/M'`` fill is algebraically ``BF := S_k/n_k`` (§4; verified to 1.96e-13 on all
    155,235 benchmark pairs by `scripts/analysis/lbf_sum_vs_combine_pips.py`). Using it here
    makes a merged component's two representations ONE object: ``alpha == softmax(lbf)`` holds
    for a merged component just as it does for an original one, inductively across merge rounds
    -- worst deviation 1.9e-15 (`scripts/analysis/check_lbf_alpha_consistency.py`).
    Posterior-neutral, and §4 proves it is the *unique* rule with that property.

    Two costs of this rule, recorded rather than hidden. (1) Bayes-factor distributions are
    heavily right-skewed, so a cohort with a sharp signal at some *other* variant has
    ``S_k/n_k >> 1`` and its unmeasured variants inherit a large Bayes factor. (2) The sum is
    not associative across nestings -- ``(a+b)+c`` differs from a flat 3-way merge, because the
    merged column is re-filled from its own ``S/n``. That is not a regression relative to
    `combine_pips`, which has never been associative under incomplete coverage (measured 0.96
    apart on the same trials), and the merge loop only ever merges pairwise in a deterministic
    greedy order, so no run is ambiguous.

    HISTORY: a ``fill`` parameter used to select between this rule ("C2") and ``lbf := 0``
    ("C1", likelihood-neutral, the default through 2026-08-19, what every stored pre-2026-08-21
    part-2 sweep table was produced with). Ran removed the option entirely on 2026-08-21:
    nothing in the process needs C1, and a knob with one legal value is a trap. Reproducing the
    stored tables now requires the pre-removal code state, not a switch. The C1-vs-C2
    measurements (§8: C1 mildly better on accuracy axes tried, top-ranked variant differs on
    16.3% of pairs; the decision rests on the consistency identity above) remain in
    ``coloc_susie_merged_lbf.pdf`` and the code-structure doc.

    The fill has NO effect on the score of a pair of two ORIGINAL components, because
    `coloc_susie_matrix` restricts to variants shared by both sides. It bites from the first
    accepted merge onward, on every pair that has a merged side. Note that does NOT mean "only
    when a merged component is merged again": the merged component sits in every subsequent
    score matrix, so its values feed both the argmax and the `(sim >= threshold).any()` stop
    test, and a run where it is never selected again can still differ because a different pair
    won or the loop stopped earlier.
    """
    from fastmap.coloc import logsum

    all_missing = lbf_df.isna().all(axis=1)
    filled = lbf_df.copy()
    for col in filled.columns:
        v = filled[col]
        measured = v.notna()
        n_k = int(measured.sum())
        if n_k == 0:
            continue
        # log(S_k / n_k) computed via logsum, not log(sum(exp(.))): SuSiE lbf reaches |lbf| ~
        # 800 in this data (the regime validate_coloc_port.py deliberately exercises), where
        # a naive exp() overflows to inf and the fill silently becomes inf.
        filled[col] = v.fillna(logsum(v[measured].to_numpy(dtype=float)) - np.log(n_k))
    merged = filled.sum(axis=1, min_count=1)
    # Re-impose "measured by nobody stays NaN": the per-column fill above would otherwise
    # promote such a variant to a real value, inventing support that no cohort provided.
    merged[all_missing] = np.nan
    return merged


def component_f1(lbf: pd.Series) -> float:
    """A component's evidence score: ``logsum(lbf_measured) - log(n_measured)``, the log of the
    arithmetic MEAN Bayes factor over the variants the component measured. NaN if it measured
    nothing.

    Under the uniform within-component prior this IS the single-effect regression's model
    log Bayes factor -- the quantity SuSiE maximises when fitting each effect's prior variance
    ``V`` and stores per component as ``lbf`` (verified equal to 6e-4 on 1,520 real components,
    2026-08-21; ``V -> 0`` forces it to 0). Unlike ``V`` it is defined identically for ORIGINAL
    and MERGED components, since a merged component carries a merged lbf vector (`merge_lbf`).
    The same quantity already appears twice elsewhere: it is `merge_lbf`'s missing-entry fill
    value ``log(S_k/n_k)`` and coloc's per-trait ``ls1``/``ls2`` term.

    Since 2026-08-21 (Ran) this is FastMap's component PRIORITY SCORE on every lbf-carrying
    path, and it is reported per shipped component as the ``comp_{i}_f1`` output column --
    exposed to the user exactly the way SuSiE exposes ``V``. Measured on 76 staged tuning
    regions at the production setting: merged components rank first in every region that has
    one (41/45 rank #1, the rest #2 behind another merge), because log-BFs ADD at a
    colocalised peak -- cross-cohort corroboration is more evidence, and this score says so.

    Overflow-safe at the |lbf| ~ 800 regime via `coloc.logsum`.
    """
    from fastmap.coloc import logsum

    v = lbf.dropna().to_numpy(dtype=float)
    if len(v) == 0:
        return float("nan")
    return float(logsum(v) - np.log(len(v)))


def combine_pips(df: pd.DataFrame) -> pd.DataFrame:
    """Combine two or more columns of per-variant PIPs into one.

    Rows missing in EVERY column are dropped (case (iii): no cohort measured the variant, so it
    is not reportable from this set). Each remaining column is scaled by its **coverage
    fraction** `c_k = n_k / M`, where `n_k` is the number of variants that column measured and
    `M = len(df)` is the region's variant count; missing entries are filled with the uniform
    prior `1 / M`; the columns are multiplied and renormalised to sum to 1.

    ``M``, not ``M'``, in both places (corrected 2026-08-20, Ran)
    ---------------------------------------------------------
    The code previously used `M' = len(df.dropna(how="all"))` -- the union support of the columns
    being merged -- for both the scaling and the fill, and the derivation note called `n_k/M'` a
    "coverage fraction", which it is not: coverage means the share of the REGION's variants that
    column measured, so the denominator is `M`.

    **The two halves must move together, and when they do the estimator is unchanged.** Only the
    RATIO of a measured entry's scale to the fill value survives renormalisation:

        c_k / f  =  (n_k / D) / (1 / D)  =  n_k        for either D in {M, M'},

    so `(n_k/M, 1/M)` and `(n_k/M', 1/M')` are the same estimator -- verified to 8.3e-17 over 300
    random pairs -- and both imply the same missing-variant Bayes factor `S_k/n_k`, the arithmetic
    mean over the variants that column measured (`coloc_susie_merged_lbf.pdf` §4, whose remark
    "x does not depend on M'" is exactly this invariance). Changing the scaling to `n_k/M` while
    leaving the fill at `1/M'` would NOT be neutral: the ratio becomes `n_k M'/M`, every fill
    Bayes factor is inflated by `M/M' >= 1`, and posterior-neutrality is lost (measured: PIPs move
    by up to 0.076, fills by a factor 1.5 at M'/M = 0.67). So this rewrite is a naming and clarity
    fix, not a numerical one.

    Note `M' = M` whenever the columns' supports cover the region -- which is always true for a
    TWO-cohort dataset, since a component's support is its cohort's variant set. So the M/M'
    distinction cannot have affected Figure 3 or Figure 4 at all; it is only ever visible with
    three or more cohorts (Figure 2's 10-cohort ABF arm, where M'/M ran ~0.53-0.59 in a synthetic
    check).
    """
    df_no_na_rows = df.dropna(how="all").copy()
    m_region = len(df)
    coverage_scale = (len(df_no_na_rows) - df_no_na_rows.isnull().sum()) / m_region
    df_scaled = df_no_na_rows * coverage_scale
    df_scaled.fillna(1 / m_region, inplace=True)
    product = df_scaled.prod(axis=1)
    total = product.sum()
    if not total > 0:
        # Reachable when the merged columns have DISJOINT exact-zero supports: SuSiE alphas
        # underflow to exactly 0.0 once a component's lbf spread exceeds ~745 (this data
        # reaches ~800), so every row's product can be 0. Renormalising would return an
        # all-NaN component that propagates through the merge loop without raising
        # (2026-08-22 review finding). Refuse loudly instead.
        raise ValueError(
            "combine_pips: the merged columns' product is zero (or non-finite) on every "
            "variant -- disjoint exact-zero supports. This merge cannot be normalised; "
            "an all-NaN component would propagate silently.")
    result = pd.DataFrame({"prob_fastmap": product / total})
    return result.reindex(df.index)


def get_components(sorted_priority_score, L=10):
    """The final selection: the first L keys of `sorted_priority_score` (F1-descending on
    lbf-carrying paths -- see `component_f1`).

    There is deliberately nothing else here. The plain-Jaccard dedup this function used to run
    was dropped on 2026-08-21 (redundant by construction on the single-score paths -- the merge
    loop only stops once every surviving pair is BELOW `similarity_threshold` on the very
    metric a same-metric dedup would test), and the coupled production path that still needed
    it was removed the same day. So `L` alone decides the size of the reported set and the
    priority-score ORDER alone decides its membership; `L` binds in all 759 staged tuning
    regions.
    """
    if not sorted_priority_score:
        # Explicit rather than a bare next(iter(...)): StopIteration names nothing about the
        # cause, and inside a generator PEP 479 would convert it to an unrelated RuntimeError.
        raise ValueError("get_components: empty ranking -- no components survived the merge")
    return list(sorted_priority_score)[:L]


def get_unused_cohorts(used_keys, cohorts, reverse_colname_map):
    """Cohorts not contributing to any of the final selected components."""
    used_cohorts = set().union(*(_cohort_names(reverse_colname_map[k]) for k in used_keys))
    return set(cohorts) - used_cohorts


def prepare_region(region_df):
    """`combine_region`'s preamble: cohort names, and the int-keyed column relabelling.

    Components are relabelled to integers because merged components are named by joining their
    parents' keys with a comma, and doing that on the original long cohort strings would make
    the keys unbounded. Returns `(renamed_df, colname_map, reverse_colname_map, cohorts)`.

    Public so that a caller needing to reuse one merge across several
    (top_n, L) settings can compose
    `prepare_region -> _merge_by_similarity -> get_components -> finalize_region` without
    duplicating this relabelling. The part-2 sweep driver does exactly that, and keeps its
    timing instrumentation outside this module (plan §6.4.2).
    """
    # sorted, not list(set(...)): set iteration order is PYTHONHASHSEED-dependent, and this
    # list reaches the backfill machinery (review finding #2, fixed 2026-08-21).
    cohorts = sorted({_cohort_of(col) for col in region_df.columns})
    colname_map = {colname: i for i, colname in enumerate(region_df.columns)}
    reverse_colname_map = {i: colname for colname, i in colname_map.items()}
    return region_df.rename(columns=colname_map), colname_map, reverse_colname_map, cohorts


class MergeState(NamedTuple):
    """Everything the greedy merge produces, before any component is *selected*.

    Split out because `L` does not participate in the merge -- it is read only by the final
    `get_components` call. So one merge can serve every `L`, which is what the part-2 sweep
    exploits: the merge is the expensive part (each accepted merge triggers a quadratic
    recompute of the similarity matrix).

    `sorted_priority_score` and `component_similarity` are read-only for consumers; `all_cols_df`
    and `all_lbf_df` must not be mutated either, since they are shared across the derived settings.
    Since 2026-08-21 the priority is `component_f1` (log mean BF) for every component whenever the
    merge ran with an lbf frame; the legacy max-alpha / 1-over-signal-size pair survives only on
    lbf-less runs. See `_merge_by_similarity`.

    `component_similarity` is the FINAL round's matrix of the SAME similarity the merge was run
    on -- `W` for weighted Jaccard, `PP.H4` for coloc-susie -- evaluated on the surviving
    components. It is DIAGNOSTIC only: no consumer selects components with it, because on these
    paths the merge already guarantees every surviving pair scores below `similarity_threshold`
    (see `_merge_by_similarity`), so a dedup pass on the same metric could never fire. Through
    2026-08-20 this field was the plain unweighted Jaccard matrix and `get_components` deduped
    on it at `component_similarity_threshold`; that mixed two metrics in one algorithm (Ran,
    2026-08-21) and the dedup stage was dropped rather than rescaled.

    `all_lbf_df` mirrors `all_cols_df` on the evidence scale: one column per component ever
    created, original or merged, holding that component's per-variant log-BF vector. It is
    `None` when the merge ran without `lbf_df` (the weighted-Jaccard path). Kept because
    `finalize_region` needs it to report `lbf_fastmap`; before 2026-08-20 the merged log-BFs were
    discarded when the loop ended, so the evidence scale existed only inside the merge.
    """
    sorted_priority_score: dict
    component_similarity: pd.DataFrame
    all_cols_df: pd.DataFrame
    all_lbf_df: pd.DataFrame = None


def _merge_by_similarity(region_df, reverse_colname_map, *, similarity, top_n,
                          similarity_threshold, lbf_df=None, coloc_priors=(1e-4, 1e-4, 1e-5),
                          overlap_min=0.5, trim_by_posterior=True):
    """Merge components greedily on a single pairwise similarity score.

    Used for the part-2 hyperparameter sweep, where the merge criterion is one threshold on
    one symmetric score (weighted Jaccard `W`, or coloc-susie `PP.H4`) rather than
    the coupled coverage/Jaccard/`div` rules of the removed production algorithm. The highest-scoring admissible pair is
    merged, scores are recomputed, and the loop repeats until nothing clears
    `similarity_threshold`.

    ONE similarity governs the whole run: `score` is bound once, below, and every round's
    stopping test and argmax read it. Nothing else enters a merge decision.

    There is deliberately NO final dedup pass here. Until 2026-08-21 `get_components` deduped
    the survivors on the plain unweighted Jaccard matrix at `component_similarity_threshold`,
    which meant a coloc-susie run merged on `PP.H4` and then selected on Jaccard -- two metrics
    in one algorithm (Ran). Rescaling that threshold onto the swept metric would not have helped:
    the loop exits only when NO surviving pair scores >= `similarity_threshold`, so a same-metric
    dedup can only ever fire in the band below the merge threshold, and is vacuous whenever
    `component_similarity_threshold >= similarity_threshold`. The stage was dropped instead, and
    `get_components` now takes the top `L` by priority score on these paths.

    Merged PIPs come from `combine_pips` (see
    `documentation/coloc_susie_merged_lbf.pdf` corollary 1, which shows the
    product-and-renormalise rule is the Bayes-correct merged posterior).
    """
    if similarity == "weighted_jaccard":
        def score(df, lbf):
            return weighted_jaccard_matrix(df, reverse_colname_map, top_n=top_n)
    elif similarity == "coloc_susie":
        if lbf_df is None:
            raise ValueError(
                "similarity='coloc_susie' requires lbf_df: per-variant log-Bayes-factor "
                "columns (SuSiE lbf_variable), one per component, aligned with region_df. "
                "The absolute lbf scale is NOT recoverable from alpha -- softmax discards it "
                "-- so it must be supplied, not derived."
            )
        if list(lbf_df.columns) != list(region_df.columns):
            raise ValueError("lbf_df columns must match region_df columns, in the same order")
        p1, p2, p12 = coloc_priors

        def score(df, lbf):
            return coloc_susie_matrix(lbf, reverse_colname_map, p1=p1, p2=p2, p12=p12,
                                      overlap_min=overlap_min,
                                      trim_by_posterior=trim_by_posterior)
    else:
        raise ValueError(f"unknown similarity={similarity!r}")

    updated_df = region_df.copy()
    all_cols_df = updated_df.copy()
    updated_lbf = None if lbf_df is None else lbf_df.copy()
    all_lbf_df = None if lbf_df is None else lbf_df.copy()
    # Priority score (decides which L components ship, since the L cap binds in practice):
    # F1 = log mean BF, one statistic for original and merged components alike, whenever an
    # lbf frame is available (Ran, 2026-08-21; see `component_f1`). Without lbf (a
    # weighted-jaccard run that wasn't given one) the legacy two-branch score remains: max
    # alpha for originals, 1/signal_size for merged -- two incomparable statistics, kept only
    # because there is nothing better to compute from alpha alone.
    if all_lbf_df is not None:
        priority_score = {col: component_f1(all_lbf_df[col]) for col in all_cols_df.columns}
    else:
        priority_score = {col: all_cols_df[col].max() for col in all_cols_df.columns}

    sim = score(updated_df, updated_lbf)
    # `cand` is where pairs get ruled out; `sim` stays a pristine matrix of the metric on the
    # CURRENT columns, so it is still meaningful when the loop exits (MergeState reports it).
    cand = sim.copy()

    while updated_df.shape[1] > 1:
        if not (cand.to_numpy(dtype=float) >= similarity_threshold).any():
            break
        col_key, row_key = _ustack_idxmax(cand)
        new_name = f"{row_key},{col_key}"

        if new_name in all_cols_df.columns:
            # Already generated this exact merge -- rule the pair out and take the next best.
            cand.loc[row_key, col_key] = cand.loc[col_key, row_key] = -1
            continue

        # A1 from the derivation: every component in a merged set must come from a distinct
        # cohort, else the same samples are counted twice and the log-BFs no longer add.
        # `weighted_jaccard_matrix` already leaves same-cohort pairs NaN; assert rather than
        # trust, because A1 is load-bearing for correctness, not a modelling preference.
        assert not shares_cohort(reverse_colname_map[row_key], reverse_colname_map[col_key]), (
            f"A1 violated: {reverse_colname_map[row_key]} and {reverse_colname_map[col_key]} "
            "share a cohort"
        )

        combined_pips = combine_pips(updated_df.loc[:, [row_key, col_key]].copy())
        updated_df.drop([row_key, col_key], axis=1, inplace=True)
        reverse_colname_map[new_name] = (
            f"{reverse_colname_map[row_key]},{reverse_colname_map[col_key]}"
        )
        updated_df[new_name] = combined_pips["prob_fastmap"]
        all_cols_df[new_name] = combined_pips["prob_fastmap"]

        if updated_lbf is not None:
            # Carry the log-BF scale forward alongside alpha (corollary 3): softmax discards
            # log S_12, the evidence strength that separates H4 from H1/H2/H0, so a merged
            # component whose lbf was not tracked could not be scored in a later round.
            merged = merge_lbf(updated_lbf.loc[:, [row_key, col_key]])
            updated_lbf.drop([row_key, col_key], axis=1, inplace=True)
            updated_lbf[new_name] = merged
            # Mirror of `all_cols_df[new_name]`: keep every component's log-BF vector, not just
            # the surviving ones, so `finalize_region` can report the lbf of whichever component
            # a variant's PIP came from.
            all_lbf_df[new_name] = merged
            priority_score[new_name] = component_f1(merged)
        else:
            signal_size = (combined_pips["prob_fastmap"] > 1 / len(combined_pips)).sum()
            priority_score[new_name] = 1 / signal_size

        sim = score(updated_df, updated_lbf)
        cand = sim.copy()

    sorted_priority_score = {k: v for k, v in
                             sorted(priority_score.items(), key=lambda item: item[1], reverse=True)
                             if k in updated_df.columns}
    return MergeState(sorted_priority_score=sorted_priority_score, component_similarity=sim,
                      all_cols_df=all_cols_df, all_lbf_df=all_lbf_df)


def _select_components_by_similarity(region_df, reverse_colname_map, *, similarity, top_n,
                                      similarity_threshold, L,
                                      lbf_df=None, coloc_priors=(1e-4, 1e-4, 1e-5),
                                      overlap_min=0.5, trim_by_posterior=True):
    """Merge, then take the top `L` components by priority score.

    No `component_similarity_threshold`: the dedup pass it controlled was dropped on
    2026-08-21 because it tested a different metric than the merge, and rescaling it onto the
    merge metric would have made it vacuous. See `_merge_by_similarity` and `get_components`.
    """
    state = _merge_by_similarity(
        region_df, reverse_colname_map, similarity=similarity, top_n=top_n,
        similarity_threshold=similarity_threshold, lbf_df=lbf_df, coloc_priors=coloc_priors,
        overlap_min=overlap_min, trim_by_posterior=trim_by_posterior)
    used_keys = get_components(state.sorted_priority_score, L=L)
    return used_keys, state.all_cols_df, state.all_lbf_df


def combine_region(region: str, region_df: pd.DataFrame, pips_df: pd.DataFrame, *,
                    similarity, similarity_threshold, top_n=500, L=10, lbf_df=None,
                    coloc_priors=(1e-4, 1e-4, 1e-5),
                    overlap_min=0.5, trim_by_posterior=True,
                    lbf_backfill="max_component") -> pd.DataFrame:
    """Run FastMap's combination algorithm for a single genomic region.

    The original coupled ``similarity="production"`` algorithm and its five parameters
    (`max_div`, `max_threshold`, `initial_reverse_max_threshold`, `jaccard_threshold`,
    `pip_threshold`), plus `component_similarity_threshold`, were REMOVED on 2026-08-21 (Ran)
    -- passing ``similarity="production"`` now raises. See the module docstring for the
    snapshot that reproduces pre-removal outputs.

    Parameters
    ----------
    region : label for this region, stored in the output's ``region`` column.
    region_df : per-variant alpha (single-effect posterior) columns, one set of L columns
        per cohort, indexed by variant. Column names must encode their cohort as
        ``{cohort}_alpha{i}`` (everything before the final underscore-separated token is
        taken as the cohort name).
    pips_df : per-variant marginal PIP columns, one per cohort, named ``{cohort}_prob``,
        same index as `region_df`. Used to backfill variants whose region a cohort simply
        wasn't run on (see "unused cohorts" handling below).
    similarity : merge criterion, required: ``"coloc_susie"`` (production, setting 107) or
        ``"weighted_jaccard"``. One threshold on one symmetric score; the highest-scoring
        admissible pair merges until nothing clears `similarity_threshold`.
    similarity_threshold : required acceptance threshold (`PP.H4` for coloc-susie, `W` for
        weighted Jaccard).
    top_n, L : `top_n` feeds only `weighted_jaccard_matrix`'s value-based top-N mask (inert
        on the coloc-susie path, kept for grid-row fidelity); `L` caps the number of final
        reported components -- selection is the top L by `component_f1` (see
        `get_components`).
    lbf_df : required for ``similarity="coloc_susie"``; strongly recommended for
        ``"weighted_jaccard"`` too, where the SCORE ignores it but the F1 priority score, the
        ``comp_{i}_f1`` columns and ``lbf_fastmap`` all need it (without it the legacy
        two-branch priority applies -- see `_merge_by_similarity`). Per-variant natural-log
        Bayes factors (SuSiE's ``lbf_variable``), one column per component, same columns and
        order as `region_df`, same index. It cannot be derived from `region_df`:
        ``alpha = softmax(lbf)`` discards the absolute scale, which is precisely the quantity
        coloc needs (see `merge_lbf`).
    coloc_priors : ``(p1, p2, p12)`` for coloc-susie, at coloc's defaults. Part 1 measured the
        ranking to be exactly `p12`-invariant and `p1 = p2 = 1e-4` optimal. Note `p1`/`p2` are
        per-*study* priors and drift in meaning once one side is itself a merged two-cohort
        component; keeping them fixed for all pairs is a recorded choice, not a derivation.
        A merged component's log-BF vector fills a cohort's MISSING entries with that cohort's
        mean Bayes factor -- the same fill `combine_pips` applies to the alpha side of the very
        same merge, so the merged component's alpha and lbf describe one object. This was the
        ``lbf_fill`` parameter ("C2", vs "C1" = ``lbf := 0``, the default through 2026-08-19)
        until Ran removed the option entirely on 2026-08-21; see `merge_lbf`'s HISTORY note.

    overlap_min, trim_by_posterior : coloc-susie pair-admissibility gate, ignored unless
        ``similarity="coloc_susie"`` (added 2026-08-21, Ran). Exact port of ``coloc.bf_bf``'s
        parameters of the same names at coloc's own defaults (0.5, TRUE): a pair may merge only
        if the variants measured on BOTH sides capture at least `overlap_min` of each side's
        posterior signal mass. See `coloc_susie_matrix` for the reduction to
        ``exp(logsum(lbf[shared]) - logsum(lbf[measured]))`` and the NaN semantics: an
        inadmissible pair is skipped THIS round only -- the matrix is recomputed after every
        accepted merge, so the pair is reconsidered each round (merging grows a component's
        measured set, so admissibility can be gained, never permanently lost).
        ``trim_by_posterior=False`` restores the pre-trim pair scores; note the stored part-2
        sweep tables are no longer reproducible from this code anyway, since they were built
        with the C1 lbf fill that was removed on 2026-08-21 (see `merge_lbf`).

    lbf_backfill : which of a cohort's L components supplies the log-BF of a BACKFILLED variant
        -- ``"max_component"`` (default) or ``"pip_argmax_component"``. See `backfill_lbf`. Read
        only when `lbf_df` is supplied, i.e. on the coloc-susie path.

    Returns
    -------
    One row per variant with a ``prob_fastmap`` column (the combined PIP), a ``note``
    column recording provenance (which component(s) it came from, or which cohort it was
    backfilled from), ``comp_i_val``/``comp_i_name`` columns for each of the up-to-L final
    components, and a ``region`` column. When `lbf_df` is supplied there are also:
    ``lbf_fastmap`` -- the log-BF of the component named in ``note`` for a covered variant,
    and the same-cohort copy described in `backfill_lbf` for a backfilled one, so PIP and
    evidence always come from the same place (Ran, 2026-08-20); and one ``comp_{i}_f1`` column
    per shipped component -- its F1 evidence score (`component_f1`), the same statistic the
    priority score selected it by, exposed the way SuSiE exposes ``V`` (Ran, 2026-08-21).
    """
    if similarity == "production":
        raise ValueError(
            "similarity='production' was removed on 2026-08-21 (Ran): the coupled "
            "coverage/Jaccard/div algorithm is no longer in this module. To reproduce a "
            "pre-removal output, use scripts/fastmap_algorithm/"
            "fastmap_v08212026_pre_production_removal.py, the last state carrying both "
            "algorithms.")
    if similarity_threshold is None:
        raise ValueError(f"similarity={similarity!r} requires similarity_threshold")

    region_df, colname_map, reverse_colname_map, cohorts = prepare_region(region_df)
    used_keys, all_cols_df, all_lbf_df = _select_components_by_similarity(
        region_df, reverse_colname_map, similarity=similarity, top_n=top_n,
        similarity_threshold=similarity_threshold, L=L,
        lbf_df=None if lbf_df is None else lbf_df.rename(columns=colname_map),
        coloc_priors=coloc_priors,
        overlap_min=overlap_min, trim_by_posterior=trim_by_posterior,
    )
    return finalize_region(region, all_cols_df, used_keys, pips_df, reverse_colname_map, cohorts,
                            all_lbf_df=all_lbf_df, cohort_lbf_df=lbf_df,
                            lbf_backfill=lbf_backfill)


def _cohort_of_note(note: str) -> str:
    """`filling_from_EUR_sim9_chr3_Omni25_HRC` -> the cohort name. Inverse of the `note` string
    `finalize_region` writes for a backfilled variant."""
    return note[len("filling_from_"):]


def backfill_lbf(cohort_lbf_df: pd.DataFrame, cohorts: pd.Series,
                  rule: str = "max_component", cohort_alpha_df: pd.DataFrame = None):
    """Per-variant log-BF copied from ONE named cohort per variant, for backfilled variants.

    A backfilled variant's PIP is that cohort's **marginal** PIP -- SuSiE's `pip`, aggregated
    over all L single effects -- and there is no single log-BF that corresponds to an aggregated
    PIP: `lbf_variable` is per component. So which of that cohort's L components to read is a
    convention, and this is where it is written down.

    ``"max_component"`` (default)
        ``max_i lbf[{cohort}_alpha{i}, j]`` -- the strongest single-effect evidence that cohort
        has for the variant. A real Bayes factor from a real component, on the absolute scale, and
        it needs no reference to alpha.
    ``"pip_argmax_component"``
        the lbf of the component with the largest ``alpha`` at that variant -- "the component the
        PIP mostly came from". Requires `alpha` as well, so the caller passes it; kept because it
        is the more literal reading of "whichever component's lbf the PIPs are filled from".

    Not offered: a log-sum-exp across the cohort's components. The L single effects are fitted
    jointly on one dataset, so adding their Bayes factors would count the same samples L times --
    the very thing assumption A1 exists to prevent.

    `cohorts` is one cohort NAME per row (index-aligned); rows whose cohort is missing, or whose
    cohort has no columns here, come back NaN rather than raising -- the `filling_from_unknown`
    branch legitimately has no cohort.
    """
    if rule not in ("max_component", "pip_argmax_component"):
        raise ValueError(f"unknown lbf backfill rule {rule!r}; expected 'max_component' or "
                         f"'pip_argmax_component'")
    if rule == "pip_argmax_component" and cohort_alpha_df is None:
        raise ValueError("rule='pip_argmax_component' needs cohort_alpha_df: the per-component "
                         "alpha frame with the same {cohort}_alpha{i} column names")

    named = cohorts.dropna()
    out = pd.Series(np.nan, index=cohorts.index, dtype=float)
    for cohort, rows in named.groupby(named).groups.items():
        cols = [c for c in cohort_lbf_df.columns if _cohort_of(c) == cohort]
        if not cols:
            continue
        lbf_block = cohort_lbf_df.loc[rows, cols]
        if rule == "max_component":
            out.loc[rows] = lbf_block.max(axis=1)
            continue
        alpha_block = cohort_alpha_df.reindex(index=rows, columns=cols)
        # A variant this cohort did not measure is NaN in every one of its components; skip it
        # rather than let `idxmax` raise, and leave the row NaN.
        has_alpha = alpha_block.notna().any(axis=1)
        if not has_alpha.any():
            continue
        pos = np.array([cols.index(k) for k in alpha_block[has_alpha].idxmax(axis=1)])
        vals = np.take_along_axis(
            lbf_block.loc[has_alpha[has_alpha].index].to_numpy(dtype=float), pos[:, None], axis=1)
        out.loc[has_alpha[has_alpha].index] = vals.ravel()
    return out


def finalize_region(region, all_cols_df, used_keys, pips_df, reverse_colname_map, cohorts, *,
                     all_lbf_df=None, cohort_lbf_df=None, lbf_backfill="max_component"):
    """Turn a selected component set into the reported per-variant table.

    Extracted verbatim from `combine_region`'s tail so that the part-2 sweep can reuse ONE
    merge across many (top_n, L) settings -- see `MergeState`. Keeping
    it as a function rather than reimplementing it in the sweep driver matters because the
    backfill rule and the `note` provenance semantics live here, and a second copy would drift.

    Does not mutate `all_cols_df` (it is shared across derived settings): every write goes to a
    `.copy()` or to a fresh frame.

    `all_lbf_df` / `cohort_lbf_df` are optional and add the `lbf_fastmap` column plus one
    `comp_{i}_f1` column per shipped component (the F1 evidence score, `component_f1`); every
    other column is unchanged whether they are passed or not. Both are needed because the two
    strata of the output get their evidence from different places, exactly mirroring how they get
    their PIP (Ran, 2026-08-20):

      * a variant some SELECTED component covers takes the log-BF of the component named in
        `note` -- i.e. the same component its `prob_fastmap` is dominated by. `all_lbf_df` is
        `MergeState.all_lbf_df`, keyed like `all_cols_df`.
      * a BACKFILLED variant takes its log-BF from the same cohort its copied marginal PIP came
        from, per `backfill_lbf`. `cohort_lbf_df` is the ORIGINAL per-component lbf frame, with
        `{cohort}_alpha{i}` column names still intact (not the int-relabelled one).

    Backfill provenance tie-break (Ran, 2026-08-21): when several unused cohorts tie on the
    marginal PIP (exact 1.0 ties are realistic), the reported source cohort is the one with the
    highest best-component F1 (`component_f1` over `cohort_lbf_df`), then alphabetical; without
    `cohort_lbf_df`, alphabetical alone. Deterministic regardless of PYTHONHASHSEED -- before
    2026-08-21 this tie broke on Python-set iteration order, so `note`/`lbf_fastmap` (never
    `prob_fastmap`) varied with the hash seed on tied rows.

    Passing neither omits the columns entirely rather than filling them with NaN, so lbf-less
    runs keep their old schema.
    """
    real_names = [reverse_colname_map[k] for k in used_keys]

    covered = all_cols_df[used_keys].dropna(how="all")
    result = covered.copy()
    # Explicit DataFrame.prod (skipna=True): a component missing at a variant contributes
    # nothing (its 1-alpha factor is skipped). np.prod used to dispatch here implicitly.
    result["prob_fastmap"] = 1 - (1 - covered[used_keys]).prod(axis=1)
    winner = covered[used_keys].idxmax(axis=1)
    result["note"] = winner.map(reverse_colname_map)

    if all_lbf_df is not None:
        # The lbf of the winning component, row by row. Positional gather rather than a
        # per-row lookup: `idxmax` already named the column, so this is one take_along_axis.
        lbf_block = all_lbf_df.reindex(index=covered.index, columns=used_keys).to_numpy(dtype=float)
        pos = np.array([used_keys.index(k) for k in winner])
        result["lbf_fastmap"] = np.take_along_axis(lbf_block, pos[:, None], axis=1).ravel()

    if len(result) < len(all_cols_df):
        # Variants no cohort in `used_keys` was run on at all -- backfill with the highest
        # marginal PIP among the cohorts that didn't make it into the final components.
        uncovered = all_cols_df[used_keys][all_cols_df[used_keys].isna().all(axis=1)].copy()
        unused_cohorts = get_unused_cohorts(used_keys, cohorts, reverse_colname_map)
        # Backfill tie-break (Ran, 2026-08-21): when two unused cohorts offer the same marginal
        # PIP (they saturate at exactly 1.0, so exact ties are realistic), prefer the cohort
        # with the STRONGER evidence -- its best component F1 -- then alphabetical. Implemented
        # by ordering the columns, since `idxmax` takes the first maximum in column order.
        # Before this, the order was Python-set iteration order, i.e. the winning cohort's name
        # in `note` (and hence `lbf_fastmap`) depended on PYTHONHASHSEED (review finding #2).
        if cohort_lbf_df is not None:
            def _cohort_f1(c):
                f1s = [component_f1(cohort_lbf_df[col]) for col in cohort_lbf_df.columns
                       if _cohort_of(col) == c]
                finite = [v for v in f1s if pd.notna(v)]
                return max(finite) if finite else float("-inf")
            ordered_cohorts = sorted(unused_cohorts, key=lambda c: (-_cohort_f1(c), c))
        else:
            ordered_cohorts = sorted(unused_cohorts)
        unused_pips = pips_df[[f"{c}_prob" for c in ordered_cohorts]].reindex(uncovered.index)
        if unused_pips.shape[1] == 0:
            # DEFENSIVE ONLY -- this should be unreachable, and reaching it means the input
            # is malformed rather than that the algorithm went wrong.
            #
            # Getting here needs every cohort to be in `used_keys` *and* some variant to be
            # NaN across all of them. With well-formed input that cannot happen:
            # `fastmap_input_processing.get_all_results` builds region_df by an outer merge
            # of each cohort's own .z variants, so every indexed variant was measured by at
            # least one cohort; a variant measured only by cohort k stays non-NaN in whichever
            # used component contains k (`combine_pips` is NaN only where *both* parents are);
            # and if k's component is not selected, `get_unused_cohorts` reports k as unused,
            # which is the normal backfill path below.
            #
            # The guard exists because `idxmax(axis=1)` raises "argmax of an empty sequence"
            # on a zero-column frame, which would take down the whole region with an error
            # that says nothing about the real cause. Falling through to the
            # `filling_from_unknown` value the code already defines is both safer and more
            # legible. Added 2026-08-14 after a synthetic test manufactured the state; no
            # real-data occurrence is known.
            uncovered["prob_fastmap"] = np.nan
            uncovered["note"] = "filling_from_unknown"
            if all_lbf_df is not None:
                uncovered["lbf_fastmap"] = np.nan
        else:
            uncovered["prob_fastmap"] = unused_pips.max(axis=1).values
            # `idxmax` on an all-NaN row relies on deprecated behaviour (pandas 2.3
            # FutureWarning: it will RAISE) -- review finding #4. Split the rows explicitly:
            # measured rows take the first maximum in the F1-then-alphabetical column order
            # above; all-NaN rows are `filling_from_unknown` without touching idxmax.
            has_pip = unused_pips.notna().any(axis=1)
            note = pd.Series("filling_from_unknown", index=uncovered.index)
            # removesuffix, not replace: replace('_prob', '') fires ANYWHERE in the string,
            # so a cohort name containing '_prob' would corrupt the note -- and the note is
            # decoded back to a cohort by _cohort_of_note for the lbf backfill.
            note[has_pip] = unused_pips[has_pip].idxmax(axis=1).apply(
                lambda x: f"filling_from_{str(x).removesuffix('_prob')}")
            uncovered["note"] = note
            if all_lbf_df is not None:
                # Same cohort the copied PIP came from, so the two reported quantities describe
                # the same evidence. `cohort_lbf_df` is required for this; without it the column
                # would silently be NaN on exactly the rows Ran asked to fill.
                if cohort_lbf_df is None:
                    raise ValueError(
                        "finalize_region: all_lbf_df was given but cohort_lbf_df was not, so the "
                        "backfilled rows could not take an lbf from their source cohort. Pass "
                        "the original {cohort}_alpha{i} lbf frame, or pass neither and omit the "
                        "lbf_fastmap column.")
                src = uncovered["note"].where(uncovered["note"] != "filling_from_unknown").dropna()
                cohort_alpha_df = None
                if lbf_backfill == "pip_argmax_component":
                    # The ORIGINAL alpha columns, under their real names: a merged component's
                    # name contains a comma, an original's does not.
                    orig = [k for k in all_cols_df.columns
                            if "," not in str(reverse_colname_map[k])]
                    cohort_alpha_df = all_cols_df[orig].rename(columns=reverse_colname_map)
                uncovered["lbf_fastmap"] = backfill_lbf(
                    cohort_lbf_df, src.map(_cohort_of_note).reindex(uncovered.index),
                    rule=lbf_backfill, cohort_alpha_df=cohort_alpha_df)
        result = pd.concat([result, uncovered], axis=0)

    result["region"] = region
    result.rename(columns={key: f"comp_{i+1}_val" for i, key in enumerate(used_keys)}, inplace=True)
    for i, name in enumerate(real_names):
        result[f"comp_{i+1}_name"] = name
    if all_lbf_df is not None:
        # Per-component evidence, exposed the way SuSiE exposes V: one scalar per shipped
        # component, constant within the region, next to its comp_{i}_name. This is the same
        # statistic the priority score ranks on (Ran, 2026-08-21; see `component_f1`).
        for i, key in enumerate(used_keys):
            result[f"comp_{i+1}_f1"] = component_f1(all_lbf_df[key])
    return result


def fastmap(results_dict: dict, pips_dict: dict, **kwargs) -> pd.DataFrame:
    """Run FastMap across every region in `results_dict`, concatenating the per-region
    results (see `combine_region` for the input format and hyperparameters)."""
    return pd.concat(
        [combine_region(region, results_dict[region], pips_dict[region], **kwargs)
         for region in results_dict],
        axis=0,
    )
