"""FastMap: combine per-cohort SuSiE fine-mapping results across cohorts into a single
set of posterior inclusion probabilities (PIPs) per variant, without joint modeling.

Each cohort contributes one or more SuSiE "single effect" components (alpha columns) per
genomic region. FastMap greedily merges components across cohorts that colocalise on one
pairwise similarity score (`PP.H4` from coloc-susie), then reports per ICS PIPs and one
final PIP per variant from the L strongest resulting components by evidence.
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pandas as pd


def _ustack_idxmax(df):
    """Equivalent of `df.infer_objects(copy=False).fillna(-2).unstack().idxmax()`:
    the `(column_name, index_name)` of the FIRST maximum in column-major order.
    
    Tie-breaking ensures the first occurrence in column-major order is returned consistently.
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
    """
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

      - the top-N mask is `~col.isin(col.nlargest(top_n))`, which is VALUE-based, not a rank
        cutoff, so ties can admit more than `top_n` entries.
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

    Each pair is scored by `coloc.combine_abf` over the variants measured on both sides.

    ``trim_by_posterior`` / ``overlap_min`` ports ``coloc.bf_bf``'s gates at defaults (TRUE, 0.5): 
    a pair is admissible only if the shared variants capture at least ``overlap_min`` of EACH 
    side's posterior signal mass. An inadmissible pair is left NaN, but the matrix is recomputed 
    after every accepted merge, so a trimmed pair can become admissible later.

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

    Requires every component in the set to come from a distinct cohort. A variant measured 
    on every side is unaffected.

    A cohort's MISSING entry takes ``lbf := log(S_k / n_k)``, the log of the arithmetic 
    mean Bayes factor over the variants column *k* did measure. A variant measured by NO 
    component in the set stays NaN, so the support is always the union of what was actually 
    measured.
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
        filled[col] = v.fillna(logsum(v[measured].to_numpy(dtype=float)) - np.log(n_k))
    merged = filled.sum(axis=1, min_count=1)
    # Re-impose "measured by nobody stays NaN"
    merged[all_missing] = np.nan
    return merged


def component_f1(lbf: pd.Series) -> float:
    """A component's evidence score: ``logsum(lbf_measured) - log(n_measured)``, the log of the
    arithmetic MEAN Bayes factor over the variants the component measured. NaN if it measured
    nothing.
    """
    from fastmap.coloc import logsum

    v = lbf.dropna().to_numpy(dtype=float)
    if len(v) == 0:
        return float("nan")
    return float(logsum(v) - np.log(len(v)))


def combine_pips(df: pd.DataFrame) -> pd.DataFrame:
    """Combine two or more columns of per-variant PIPs into one.

    Rows missing in EVERY column are dropped. Each remaining column is scaled by its 
    coverage fraction `c_k = n_k / M`, where `n_k` is the number of variants that column 
    measured and `M = len(df)` is the region's variant count; missing entries are filled 
    with the uniform prior `1 / M`; the columns are multiplied and renormalised to sum to 1.
    """
    df_no_na_rows = df.dropna(how="all").copy()
    m_region = len(df)
    coverage_scale = (len(df_no_na_rows) - df_no_na_rows.isnull().sum()) / m_region
    df_scaled = df_no_na_rows * coverage_scale
    df_scaled.fillna(1 / m_region, inplace=True)
    product = df_scaled.prod(axis=1)
    total = product.sum()
    if not total > 0:
        raise ValueError(
            "combine_pips: the merged columns' product is zero (or non-finite) on every "
            "variant -- disjoint exact-zero supports. This merge cannot be normalised.")
    result = pd.DataFrame({"prob_fastmap": product / total})
    return result.reindex(df.index)


def get_components(sorted_priority_score, L=10):
    """The final selection: the first L keys of `sorted_priority_score`."""
    if not sorted_priority_score:
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
    """
    cohorts = sorted({_cohort_of(col) for col in region_df.columns})
    colname_map = {colname: i for i, colname in enumerate(region_df.columns)}
    reverse_colname_map = {i: colname for colname, i in colname_map.items()}
    return region_df.rename(columns=colname_map), colname_map, reverse_colname_map, cohorts


class MergeState(NamedTuple):
    """Everything the greedy merge produces, before any component is selected.
    
    `sorted_priority_score` and `component_similarity` are read-only for consumers; `all_cols_df`
    and `all_lbf_df` must not be mutated either, since they are shared across derived settings.
    """
    sorted_priority_score: dict
    component_similarity: pd.DataFrame
    all_cols_df: pd.DataFrame
    all_lbf_df: pd.DataFrame = None


def _merge_by_similarity(region_df, reverse_colname_map, *, similarity, top_n,
                         similarity_threshold, lbf_df=None, coloc_priors=(1e-4, 1e-4, 1e-5),
                         overlap_min=0.5, trim_by_posterior=True):
    """Merge components greedily on a single pairwise similarity score.

    The highest-scoring admissible pair is merged, scores are recomputed, and the loop 
    repeats until nothing clears `similarity_threshold`.
    """
    if similarity == "weighted_jaccard":
        def score(df, lbf):
            return weighted_jaccard_matrix(df, reverse_colname_map, top_n=top_n)
    elif similarity == "coloc_susie":
        if lbf_df is None:
            raise ValueError(
                "similarity='coloc_susie' requires lbf_df: per-variant log-Bayes-factor "
                "columns (SuSiE lbf_variable), one per component, aligned with region_df."
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
    
    if all_lbf_df is not None:
        priority_score = {col: component_f1(all_lbf_df[col]) for col in all_cols_df.columns}
    else:
        priority_score = {col: all_cols_df[col].max() for col in all_cols_df.columns}

    sim = score(updated_df, updated_lbf)
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
            merged = merge_lbf(updated_lbf.loc[:, [row_key, col_key]])
            updated_lbf.drop([row_key, col_key], axis=1, inplace=True)
            updated_lbf[new_name] = merged
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
    """Merge, then take the top `L` components by priority score."""
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

    Parameters
    ----------
    region : label for this region, stored in the output's ``region`` column.
    region_df : per-variant alpha (single-effect posterior) columns, one set of L columns
        per cohort, indexed by variant. Column names must encode their cohort as
        ``{cohort}_alpha{i}``.
    pips_df : per-variant marginal PIP columns, one per cohort, named ``{cohort}_prob``,
        same index as `region_df`. Used to backfill variants whose region a cohort simply
        wasn't run on.
    similarity : merge criterion, required: ``"coloc_susie"`` or ``"weighted_jaccard"``.
    similarity_threshold : required acceptance threshold.
    top_n, L : `top_n` feeds only `weighted_jaccard_matrix`'s value-based top-N mask; 
        `L` caps the number of final reported components.
    lbf_df : required for ``similarity="coloc_susie"``; recommended for
        ``"weighted_jaccard"`` too. Per-variant natural-log Bayes factors.
    coloc_priors : ``(p1, p2, p12)`` for coloc-susie, at coloc's defaults.
    overlap_min, trim_by_posterior : coloc-susie pair-admissibility gate, ignored unless
        ``similarity="coloc_susie"``.
    lbf_backfill : which of a cohort's L components supplies the log-BF of a BACKFILLED variant
        -- ``"max_component"`` (default) or ``"pip_argmax_component"``.

    Returns
    -------
    One row per variant with a ``prob_fastmap`` column (the combined PIP), a ``note``
    column recording provenance, ``comp_i_val``/``comp_i_name`` columns for each of the 
    up-to-L final components, and a ``region`` column. When `lbf_df` is supplied there are also 
    ``lbf_fastmap`` and one ``comp_{i}_f1`` column per shipped component.
    """
    if similarity == "production":
        raise ValueError("similarity='production' is no longer supported.")
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

    ``"max_component"`` (default)
        ``max_i lbf[{cohort}_alpha{i}, j]`` -- the strongest single-effect evidence that cohort
        has for the variant.
    ``"pip_argmax_component"``
        the lbf of the component with the largest ``alpha`` at that variant. Requires `alpha`.
    """
    if rule not in ("max_component", "pip_argmax_component"):
        raise ValueError(f"unknown lbf backfill rule {rule!r}; expected 'max_component' or "
                         f"'pip_argmax_component'")
    if rule == "pip_argmax_component" and cohort_alpha_df is None:
        raise ValueError("rule='pip_argmax_component' needs cohort_alpha_df")

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
    """Turn a selected component set into the reported per-variant table."""
    real_names = [reverse_colname_map[k] for k in used_keys]

    covered = all_cols_df[used_keys].dropna(how="all")
    result = covered.copy()
    
    result["prob_fastmap"] = 1 - (1 - covered[used_keys]).prod(axis=1)
    winner = covered[used_keys].idxmax(axis=1)
    result["note"] = winner.map(reverse_colname_map)

    if all_lbf_df is not None:
        lbf_block = all_lbf_df.reindex(index=covered.index, columns=used_keys).to_numpy(dtype=float)
        pos = np.array([used_keys.index(k) for k in winner])
        result["lbf_fastmap"] = np.take_along_axis(lbf_block, pos[:, None], axis=1).ravel()

    if len(result) < len(all_cols_df):
        uncovered = all_cols_df[used_keys][all_cols_df[used_keys].isna().all(axis=1)].copy()
        unused_cohorts = get_unused_cohorts(used_keys, cohorts, reverse_colname_map)
        
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
            uncovered["prob_fastmap"] = np.nan
            uncovered["note"] = "filling_from_unknown"
            if all_lbf_df is not None:
                uncovered["lbf_fastmap"] = np.nan
        else:
            uncovered["prob_fastmap"] = unused_pips.max(axis=1).values
            has_pip = unused_pips.notna().any(axis=1)
            note = pd.Series("filling_from_unknown", index=uncovered.index)
            note[has_pip] = unused_pips[has_pip].idxmax(axis=1).apply(
                lambda x: f"filling_from_{str(x).removesuffix('_prob')}")
            uncovered["note"] = note
            if all_lbf_df is not None:
                if cohort_lbf_df is None:
                    raise ValueError(
                        "finalize_region: all_lbf_df was given but cohort_lbf_df was not. Pass "
                        "the original {cohort}_alpha{i} lbf frame, or pass neither and omit the "
                        "lbf_fastmap column.")
                src = uncovered["note"].where(uncovered["note"] != "filling_from_unknown").dropna()
                cohort_alpha_df = None
                if lbf_backfill == "pip_argmax_component":
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
