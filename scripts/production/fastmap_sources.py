#!/usr/bin/env python3
"""Input adapters: per-cohort fine-mapping output -> the three frames FastMap consumes.

One adapter per manuscript figure. All three return the SAME contract, so `run_fastmap.py`
carries no source-specific logic:

    alpha  per-variant single-effect posteriors, columns "{cohort}_alpha{i}"   -> region_df
    lbf    per-variant natural-log Bayes factors, SAME columns, SAME order     -> lbf_df
    pip    per-variant marginal PIP, one "{cohort}_prob" column per cohort     -> pips_df

all indexed by variant id, outer-merged across cohorts (so a variant measured by any cohort
appears exactly once and unmeasured entries are NaN), and column-aligned.

  sim_1causal      Figure 2   per-cohort ABF `.abf.snp`, ONE component per cohort, 10 cohorts
  sim_multicausal  Figure 3   per-cohort SuSiE `.susie.rds` + `.z`, 10 components, 2 cohorts
  real_data        Figure 4   FinnGen aggregated feather + UKBB per-region `.susie.rds` + `.z`

Why `lbf` is loaded everywhere and cannot be derived
----------------------------------------------------
The production setting is coloc-susie, whose score is `coloc.combine_abf` over per-variant
log-Bayes-factor vectors. `alpha = softmax(lbf)` discards the absolute scale, and that scale is
exactly the quantity coloc's `S12` term integrates -- so lbf must be read from source, never
reconstructed from alpha. `_select_components_by_similarity` refuses to run without it.

Cohort naming is load-bearing
-----------------------------
`fastmap._cohort_names()` takes the cohort of a component to be everything before the final
underscore-separated token, and A1 (a merged component must draw on distinct cohorts, else the
same samples are counted twice) is enforced by comparing those strings. So each adapter strips
the phenotype token out of the cohort name -- `EAS_sim1_chr3_pheno26_Omni25_1000G` becomes
`EAS_sim1_chr3_Omni25_1000G` -- matching production `fastmap_input_processing.get_all_results`,
and the real-data adapter renames FinnGen/UKBB columns into the same `{cohort}_alpha{i}` shape
rather than passing through the raw `alpha1` / `alpha_1` names (see `load_real_data`).

The rpy2 view bug
-----------------
`np.asarray()` on an rpy2 object returns a VIEW into R's memory (`owndata is False`). The R
handle is reassigned on every loop iteration, so R's GC reclaims the old buffer and any
DataFrame still wrapping it reads whatever now occupies those bytes -- observed as alpha frames
with negative "posteriors" (`lbf_variable` bleeding through) and a hard segfault. Every
`np.asarray` here is followed by `.copy()`. This is not defensive; removing it silently
corrupts output.
"""
from __future__ import annotations

import glob
import json
import multiprocessing as mp
import os
import re
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from functools import reduce

import numpy as np
import pandas as pd

#: Root GCS bucket for simulation/real-data inputs. Set FASTMAP_DATA_ROOT to your own bucket;
#: the paths below (private lab infrastructure) are examples of the expected layout, not
#: a public dataset.
DATA_ROOT = os.environ.get("FASTMAP_DATA_ROOT", "gs://YOUR-BUCKET")

SIMS = f"{DATA_ROOT}/simulations/sims"
CONFIGS = f"{DATA_ROOT}/simulations/scripts/configs"

# Figure 2 / Figure 3 production cohort sets, as actually run (cross-checked against the
# Cromwell inputs JSONs that produced the shipped results, not inferred).
SIM_SOURCES = {
    "sim_1causal": {
        "config_tsv": f"{CONFIGS}/All_random_configs10_Omni25_noTOPMed.tsv",
        "config_index": 1,                       # 1-based row, matching get_individual_cohorts.py
        "assoc_dir": f"{SIMS}/assoc_1causal",
        "meta_dir": f"{SIMS}/meta_analysis_1causal",
        "kind": "abf",
        "susie_l": 1,
    },
    "sim_multicausal": {
        "config_tsv": f"{CONFIGS}/All_random_configs2_Omni25_noTOPMed.tsv",
        "config_index": 1,                       # -> AFR_sim1_Omni25_1000G + EUR_sim9_Omni25_HRC
        "assoc_dir": f"{SIMS}/assoc_3causal_27h2g_fixed",
        "meta_dir": f"{SIMS}/meta_analysis_3causal_27h2g_fixed",
        "kind": "susie",
        "susie_l": 10,
    },
}

# Figure 4: FinnGen R12 endpoint <-> UKBB phenotype. Height is included (it is excluded at
# PLOTTING time, not here -- its 834 regions dominate any pooled count, so Figure 4 calls it
# out separately rather than dropping it from the run).
# The 15 phenotypes of the June-2026 analysis (Daly lab deck slide 8), FG endpoint -> UKBB
# name. Every FG key below was VERIFIED 2026-08-21 against the actual feathers in
# gs://.../FINNGEN_meta_regions/UKBB_FINNGEN_R12/June_2026/ -- four earlier keys (I9_MI,
# J10_ASTHMA, K11_CROHN, M13_RHEUMA) were remembered names that do not exist there and would
# have 404ed at staging; the real endpoints are I9_MI_STRICT, J10_ASTHMA_EXMORE,
# K11_CD_STRICT2, RHEUMA_SEROPOS_WIDE.
REAL_DATA_PHENOS = {
    "I9_CHD": "CAD", "I9_MI_STRICT": "MI", "I9_AF": "AFib", "C_STROKE": "IS",
    "E4_HYTHY_AI_STRICT": "Hypothyroidism", "T2D": "T2D", "T1D": "T1D",
    "K11_CD_STRICT2": "Crohn", "RHEUMA_SEROPOS_WIDE": "RA",
    "H7_CATARACTSENILE": "Cataract",
    "C3_BRONCHUS_LUNG_EXALLC": "LuC", "C3_SKIN_EXALLC": "SkC",
    "G6_AD_WIDE": "Alzheimer_LTFH",
    "BMI_IRN": "BMI", "HEIGHT_IRN": "Height",
}
# Present in GCS but NOT in the 15-phenotype analysis: UKBB SuSiE dirs + FG .snp exist for
# these, but the FG June_2026 FEATHER (what stage_real_data reads) was never written --
# regenerate it from the .snp if any is ever added. C3_BREAST_EXALLC additionally has no
# UKBB counterpart at all.
REAL_DATA_CANDIDATES_NO_FG_FEATHER = {
    "J10_ASTHMA_EXMORE": "Asthma_Combined", "H7_GLAUCOMA": "Glaucoma_Combined",
    "C3_PROSTATE_EXALLC": "PrC_M",
}
FG_FEATHER = f"{DATA_ROOT}/results/FINNGEN_meta_regions/UKBB_FINNGEN_R12/June_2026/{{fg}}.susie.all_snps.feather"
UKB_SUSIE_DIR = f"{DATA_ROOT}/results/UKBB_meta_regions/UKBB_FINNGEN_R12/{{ukb}}"
UKB_Z_DIR = f"{DATA_ROOT}/data/FINNGEN_UKBB_meta/R12/finemap_inputs/{{ukb}}"
REGION_MAP = f"{DATA_ROOT}/data/FINNGEN_UKBB_meta/R12/meta_generated_bed_files/June_2026/{{fg}}_38_to_37_mapping.tsv"
#: Per-variant b37->b38 maps precomputed by an internal liftover pipeline (BCFtools/+liftover,
#: pyliftover-gated; not included in this repo). The Fig 4 cross-cohort join keys on lifted
#: b38 (chrom,pos,ref,alt), NOT on rsid: 4.5% of FG variants are unnamed (the rsid index
#: crashed the loader), rsid matched only 11.7% of indels, and mispaired 643 T1D multi-allelics
#: with the wrong allele (quantified 2026-08-24).
UKB_LIFTOVER_MAP = f"{DATA_ROOT}/data/FINNGEN_UKBB_meta/R12/variant_liftover_b37_to_b38/{{ukb}}.b37_to_b38.tsv.gz"


def norm_chrom_b37(c) -> str:
    """.z chromosome value ('01'..'22', '23'/'X') -> '1'..'22'/'X' (no chr prefix)."""
    s = str(c).strip()
    if s.startswith("chr"):
        s = s[3:]
    if s.upper() in ("X", "23"):
        return "X"
    n = int(s)
    if not 1 <= n <= 22:
        raise ValueError(f"unexpected chromosome {c!r}")
    return str(n)


def z_src_keys(z: pd.DataFrame) -> pd.Series:
    """The liftover map's `src_key` ('1:12345:A:G'), reproduced from .z columns. Must stay
    byte-identical to `build_ukb_liftover_maps.read_z_variants`."""
    return (z["chromosome"].map(norm_chrom_b37) + ":" + z["position"].astype(str) + ":"
            + z["allele1"].str.upper() + ":" + z["allele2"].str.upper())


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------

def _run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def _gsutil_cp(patterns, dest, quiet=True):
    """Multi-file copy into `dest`. Uses `gcloud storage cp` rather than `gsutil -m cp`: same
    semantics (recursive multi-pattern copy, still `check=True`-guarded), but the rewritten
    transfer tool has much better throughput on the many-small-file glob patterns
    (`*.susie.rds`/`*.z`) this staging path pulls -- confirmed empirically 2026-09-10 on
    C_STROKE/IS's 19-region staging."""
    os.makedirs(dest, exist_ok=True)
    cmd = ["gcloud", "storage", "cp"] + (["-q"] if quiet else []) + list(patterns) + [dest]
    _run(cmd)


STAGING_MARKER = ".staging_complete"


def _ensure_staged(dest: str, patterns: list) -> None:
    """Download `patterns` into `dest` unless a prior COMPLETE download is recorded there.

    The re-entrancy check is a marker file written only AFTER gsutil returns success -- NOT
    "any matching file exists". An interrupted `gsutil -m cp` leaves a partial file set
    behind; a glob-based check would then skip re-staging forever and the loaders would
    silently drop the missing regions from production outputs while the run is recorded as
    done (2026-08-22 review finding #1). `gsutil -m cp` exits nonzero on any per-file
    failure, and `_run` uses check=True, so reaching the marker write means every file
    matched by `patterns` was copied. Re-running the cp over an existing partial download
    just overwrites -- idempotent, and cheap next to the correctness.
    """
    marker = os.path.join(dest, STAGING_MARKER)
    if os.path.exists(marker):
        return
    if os.path.isdir(dest) and os.listdir(dest):
        print(f"  {dest}: files present but no {STAGING_MARKER} -- prior download was "
              "interrupted; re-staging", flush=True)
    _gsutil_cp(patterns, dest)
    n_files = len([f for f in os.listdir(dest) if f != STAGING_MARKER])
    with open(marker, "w") as fh:
        json.dump({"patterns": list(patterns), "n_files": n_files,
                   "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}, fh)


def _gsutil_cp_file(src: str, local: str) -> None:
    """Download one file atomically: fetch to a temp name, os.replace into place.

    `os.path.exists(local)` is the re-entrancy check for single files, so `local` must never
    name a partially written download.
    """
    if os.path.exists(local):
        return
    tmp = local + ".staging_tmp"
    _run(["gsutil", "-q", "cp", src, tmp])
    os.replace(tmp, local)


def _outer(a, b):
    return pd.merge(a, b, how="outer", left_index=True, right_index=True)


def cohorts_for_config(config_tsv: str, config_index: int, local_root: str) -> list:
    """Cohort names for one config row. `config_index` is 1-BASED, as in the WDL inputs.

    Reproduces `get_individual_cohorts.py`: read the whitespace-delimited config with no
    header, take row `config_index - 1`. Off-by-one here would silently fine-map a different
    cohort set, so the 1-based convention is kept rather than "cleaned up".
    """
    os.makedirs(local_root, exist_ok=True)
    local = os.path.join(local_root, os.path.basename(config_tsv))
    if not os.path.exists(local):
        _run(["gsutil", "-q", "cp", config_tsv, local])
    df = pd.read_csv(local, delimiter=r"\s+", header=None)
    if not 1 <= config_index <= len(df):
        raise ValueError(f"config_index {config_index} out of range for {config_tsv} "
                         f"({len(df)} rows)")
    return [str(x) for x in df.iloc[config_index - 1, :].to_list()]


def _with_pheno(cohort: str, pheno: int) -> str:
    """`EUR_sim9_chr3_Omni25_HRC` -> `EUR_sim9_chr3_pheno12_Omni25_HRC`, the on-disk form."""
    return cohort.replace("_chr3_", f"_chr3_pheno{pheno}_", 1)


def _stage_regions_bed(source: str, pheno: int, config_index: int, local_root: str) -> list:
    """Download and parse one phenotype's meta-analysis region BED -> region-name list.

    The shipped pipeline was DRIVEN by this bed (`get_cohort_region_result_names.py` turned
    its rows into `chr{c}.{start}-{end}` filename tokens and the WDL demanded a unique
    per-cohort match for each). The per-cohort assoc directories contain MORE `.abf.snp` /
    `.susie.rds` files than the bed names -- overlapping boundary variants of the same loci
    and non-meta-significant regions -- so globbing "whatever files exist" both over-counts
    regions and trips the all-cohorts completeness check on files the shipped run never read
    (pheno1 1-causal: 27 globbed names vs the bed's 11, which reproduce the shipped output's
    region set exactly).
    """
    spec = SIM_SOURCES[source]
    config_name = os.path.basename(spec["config_tsv"]).replace(".tsv", "")
    stem = f"meta.pheno{pheno}.config{config_index}"
    remote = f"{spec['meta_dir']}/{config_name}/{stem}/{stem}.bed"
    local = os.path.join(local_root, f"pheno{pheno}", f"{stem}.bed")
    os.makedirs(os.path.dirname(local), exist_ok=True)
    _gsutil_cp_file(remote, local)
    bed = pd.read_csv(local, header=None, delimiter="\t")
    regions = [f"chr{c}.{s}-{e}" for c, s, e in bed.itertuples(index=False)]
    if not regions:
        raise RuntimeError(f"{remote}: empty region bed")
    if len(set(regions)) != len(regions):
        raise RuntimeError(f"{remote}: duplicated region rows")
    return regions


def _restrict_to_bed(per_cohort: dict, bed_regions: list, accounting) -> dict:
    """Keep only the bed's regions in each cohort's file map; extra on-disk files are
    IGNORED by design (they are what the shipped run never read), but their count is
    recorded so the manifest shows what was set aside. A bed region with no file in some
    cohort stays absent from that cohort's map, so `_resolve_regions` names it loudly."""
    wanted = set(bed_regions)
    # A bed region with no file in ANY cohort would otherwise vanish from the union that
    # `_resolve_regions` iterates -- the one silent-drop path its strict mode cannot see.
    nowhere = wanted - set().union(*(set(f) for f in per_cohort.values()))
    if nowhere:
        raise RuntimeError(
            f"{len(nowhere)} bed region(s) have no staged file in any cohort: "
            f"{sorted(nowhere)} -- staging is incomplete or the bed does not match this "
            "assoc directory")
    ignored = {c: sorted(set(files) - wanted) for c, files in per_cohort.items()}
    if accounting is not None:
        accounting["regions_bed"] = len(bed_regions)
        accounting["files_on_disk_not_in_bed"] = {
            c: len(extra) for c, extra in ignored.items() if extra}
    n_extra = sum(len(v) for v in ignored.values())
    if n_extra:
        print(f"  {n_extra} per-cohort region file(s) not named by the meta bed are ignored "
              f"(as in the shipped pipeline); counts per cohort are in the manifest")
    return {c: {r: f for r, f in files.items() if r in wanted}
            for c, files in per_cohort.items()}


def _resolve_regions(per_cohort: dict, min_cohorts, accounting) -> list:
    """The `(region, present_cohorts)` work list, with completeness enforced LOUDLY.

    `per_cohort` maps cohort -> {region: ...files...}. Default (`min_cohorts=None`): every
    region must be present in EVERY cohort, and anything less RAISES with the offenders named
    -- this matches the legacy WDL path, whose `get_fname` demanded a unique match per
    cohort x region. Silent-skip modes are opt-in only: pass an integer `min_cohorts` to
    tolerate regions present in at least that many cohorts; the skips are then recorded in
    `accounting["regions_skipped"]` and printed, never dropped without trace
    (2026-08-22 review: silent region loss is the worst failure mode this pipeline has).
    An entirely absent cohort is always fatal in either mode.
    """
    cohorts = list(per_cohort)
    empty = [c for c in cohorts if not per_cohort[c]]
    if empty:
        raise RuntimeError(
            f"cohort(s) with ZERO regions on disk: {empty} -- staging is incomplete or the "
            "source is missing; refusing to proceed (every region would look "
            "cohort-incomplete and the output would be silently wrong)")
    all_regions = sorted(set().union(*(set(v) for v in per_cohort.values())))
    work, skipped = [], []
    for region in all_regions:
        present = [c for c in cohorts if region in per_cohort[c]]
        if len(present) < (len(cohorts) if min_cohorts is None else min_cohorts):
            skipped.append((region, len(present),
                            [c for c in cohorts if region not in per_cohort[c]]))
        else:
            work.append((region, present))
    if accounting is not None:
        accounting["regions_per_cohort"] = {c: len(per_cohort[c]) for c in cohorts}
        accounting["regions_union"] = len(all_regions)
        accounting["regions_processed"] = len(work)
        accounting["regions_skipped"] = [
            {"region": r, "n_present": n, "missing_from": miss} for r, n, miss in skipped]
    if skipped and min_cohorts is None:
        detail = "; ".join(f"{r} missing from {miss}" for r, _, miss in skipped[:5])
        raise RuntimeError(
            f"{len(skipped)} region(s) are not present in every cohort ({detail}"
            f"{' ...' if len(skipped) > 5 else ''}). This is either incomplete staging or "
            "missing upstream output. Refusing to skip silently; pass min_cohorts=<k> "
            "explicitly if partial-cohort regions are genuinely intended.")
    if skipped:
        print(f"  skipped {len(skipped)} region(s) present in < {min_cohorts} cohorts: "
              f"{[(r, n) for r, n, _ in skipped[:5]]}{' ...' if len(skipped) > 5 else ''}")
    return work


def _check_frames(region, alpha, lbf, pip):
    """The invariants every adapter must satisfy before FastMap sees the frames."""
    if list(alpha.columns) != list(lbf.columns):
        raise AssertionError(f"{region}: alpha/lbf columns differ")
    if not (alpha.index.equals(lbf.index) and alpha.index.equals(pip.index)):
        raise AssertionError(f"{region}: frame indices disagree")
    if alpha.index.duplicated().any():
        raise AssertionError(f"{region}: {int(alpha.index.duplicated().sum())} duplicated "
                             "variant ids -- an outer merge would have fanned out")
    # A variant present in the index must have been measured by someone, else `combine_pips`
    # and the backfill both have nothing to work from.
    if alpha.isna().all(axis=1).any():
        raise AssertionError(f"{region}: variants with no alpha in any cohort")


# ---------------------------------------------------------------------------------------------
# Figure 2 -- single-causal simulations, per-cohort ABF
# ---------------------------------------------------------------------------------------------

def stage_sim_1causal(pheno: int, local_root: str, config_index: int = 1) -> dict:
    """Download one phenotype's per-cohort ABF output plus the meta-region bed that names
    the regions to process. Re-entrant via `_ensure_staged`'s completion marker (never via
    "some files exist")."""
    spec = SIM_SOURCES["sim_1causal"]
    cohorts = cohorts_for_config(spec["config_tsv"], config_index, local_root)
    dest_root = os.path.join(local_root, f"pheno{pheno}")
    regions = _stage_regions_bed("sim_1causal", pheno, config_index, local_root)
    for cohort in cohorts:
        cp = _with_pheno(cohort, pheno)
        _ensure_staged(os.path.join(dest_root, cp),
                       [f"{spec['assoc_dir']}/pheno{pheno}/{cp}/*.abf.snp"])
    return {"cohorts": cohorts, "dir": dest_root, "bed_regions": regions}


def load_sim_1causal(pheno: int, staged: dict, min_cohorts: int = None, accounting: dict = None):
    """Yield `(region, alpha, lbf, pip)` for the 1-causal sims.

    ABF gives ONE component per cohort, so each cohort contributes a single `{cohort}_alpha1`
    column whose values are the ABF posterior (`prob`) and whose lbf is the ABF `lbf`. On
    single-component inputs coloc-susie's score reduces exactly to coloc's `combine.abf`, which
    is what part 1 benchmarked -- consistent, not a special case.

    The work list is the meta-region bed staged alongside the inputs (`staged["bed_regions"]`),
    exactly as in the shipped pipeline; per-cohort files the bed does not name are ignored.
    Region completeness is then enforced by `_resolve_regions`: by default a bed region
    missing from ANY cohort raises (matching the legacy WDL path); pass an integer
    `min_cohorts` only if partial-cohort regions are genuinely intended, and they are then
    recorded in `accounting`.
    """
    cohorts, root = staged["cohorts"], staged["dir"]
    per_cohort = {}
    for cohort in cohorts:
        cp = _with_pheno(cohort, pheno)
        per_cohort[cohort] = {
            os.path.basename(f)[len(cp) + 1:-len(".abf.snp")]: f
            for f in glob.glob(f"{root}/{cp}/*.abf.snp")}
    per_cohort = _restrict_to_bed(per_cohort, staged["bed_regions"], accounting)
    for region, present in _resolve_regions(per_cohort, min_cohorts, accounting):
        a_parts, l_parts, p_parts = [], [], []
        for cohort in present:
            d = pd.read_csv(per_cohort[cohort][region], delimiter=r"\s+",
                            usecols=["rsid", "prob", "lbf"])
            idx = pd.Index(d["rsid"], name="variant")
            if idx.duplicated().any():
                raise ValueError(f"{per_cohort[cohort][region]}: duplicated rsid")
            a_parts.append(pd.DataFrame({f"{cohort}_alpha1": d["prob"].values}, index=idx))
            l_parts.append(pd.DataFrame({f"{cohort}_alpha1": d["lbf"].values}, index=idx))
            p_parts.append(pd.DataFrame({f"{cohort}_prob": d["prob"].values}, index=idx))
        alpha = reduce(_outer, a_parts)
        lbf = reduce(_outer, l_parts)[alpha.columns]
        pip = reduce(_outer, p_parts)
        _check_frames(region, alpha, lbf, pip)
        yield region, alpha, lbf, pip


# ---------------------------------------------------------------------------------------------
# Figure 3 -- multi-causal simulations, per-cohort SuSiE
# ---------------------------------------------------------------------------------------------

def stage_sim_multicausal(pheno: int, local_root: str, config_index: int = 1) -> dict:
    spec = SIM_SOURCES["sim_multicausal"]
    cohorts = cohorts_for_config(spec["config_tsv"], config_index, local_root)
    dest_root = os.path.join(local_root, f"pheno{pheno}")
    for cohort in cohorts:
        cp = _with_pheno(cohort, pheno)
        src = f"{spec['assoc_dir']}/pheno{pheno}/{cp}"
        # One marker covers BOTH patterns: the old code checked only *.susie.rds, so an
        # interruption between the rds and z copies looked complete.
        _ensure_staged(os.path.join(dest_root, cp), [f"{src}/*.susie.rds", f"{src}/*.z"])
    regions = _stage_regions_bed("sim_multicausal", pheno, config_index, local_root)
    return {"cohorts": cohorts, "dir": dest_root, "bed_regions": regions}


def _read_susie_rds(rds_path: str, z_path: str, cohort: str, susie_l: int,
                    variant_col: str = "rsid", index_func=None):
    """One cohort x region SuSiE fit -> (alpha, lbf, pip) frames, plus variant positions.

    `index_func(z_df) -> pd.Index` overrides the default rsid index; the real-data loader
    uses it to key UKBB variants by their lifted b38 composite ids."""
    import rpy2.robjects as ro

    r = ro.r["readRDS"](rds_path)
    z = pd.read_csv(z_path, delimiter=r"\s+")
    if variant_col not in z.columns:
        raise ValueError(f"{z_path}: no '{variant_col}' column (has {list(z.columns)[:8]})")

    # `.copy()` is REQUIRED -- see the rpy2 view bug in the module docstring.
    alpha = np.asarray(r.rx2("alpha")).copy()            # L x n_variants
    lbf = np.asarray(r.rx2("lbf_variable")).copy()       # L x n_variants
    pip = np.asarray(r.rx2("pip"), dtype=float).copy()   # n_variants

    if alpha.shape != lbf.shape:
        raise ValueError(f"{rds_path}: alpha {alpha.shape} != lbf_variable {lbf.shape}")
    if alpha.shape[0] != susie_l:
        raise ValueError(f"{rds_path}: alpha has {alpha.shape[0]} rows, expected L={susie_l}")
    if alpha.shape[1] != len(z) or len(pip) != len(z):
        raise ValueError(f"{rds_path}: {alpha.shape[1]} alpha / {len(pip)} pip variants but "
                         f"{z_path} has {len(z)}")

    cols = [f"{cohort}_alpha{i}" for i in range(1, susie_l + 1)]
    idx = (index_func(z) if index_func is not None
           else pd.Index(z[variant_col], name="variant"))
    a = pd.DataFrame(alpha.T, index=idx, columns=cols)
    l = pd.DataFrame(lbf.T, index=idx, columns=cols)
    # SuSiE's OWN `pip`, not `alpha1`: alpha1 is one single effect while pip aggregates all L,
    # so in a 3-causal simulation every variant belonging to the 2nd or 3rd signal would get a
    # near-zero stand-in. Production uses `pip` and so do we.
    p = pd.DataFrame({f"{cohort}_prob": pip}, index=idx)
    pos = (pd.DataFrame({"position": z["position"].values}, index=idx)
           if "position" in z.columns else None)
    return a, l, p, pos


_worker_src_to_b38 = None  # set once per worker process by `_init_ukb_worker`


def _init_ukb_worker(src_to_b38: dict) -> None:
    """`ProcessPoolExecutor(initializer=...)` hook: stashes the liftover map in this worker's
    global once, rather than pickling and re-sending it (potentially millions of entries) on
    every one of thousands of per-region task submissions."""
    global _worker_src_to_b38
    _worker_src_to_b38 = src_to_b38


def _ukb_region_worker(region38: str, rds_path: str, z_path: str, cohort: str, susie_l: int):
    """One UKBB region's `_read_susie_rds`, run in a worker process (see `load_real_data`).

    Top-level and picklable-args-only by design: `ProcessPoolExecutor` must ship this
    function and its arguments to the worker, and `rpy2`'s embedded R is NOT fork-safe, so
    the pool uses the `spawn` start method -- each worker imports rpy2 and starts its own R
    interpreter fresh rather than inheriting the parent's (which may already have R loaded).
    Reproduces `ukb_index`'s b38-lookup and duprep/unlifted accounting exactly, just returning
    the counts instead of mutating a closed-over `nonlocal` (which can't cross a process
    boundary) so the caller can sum them after collecting all workers' results.
    """
    counts = {"unlifted": 0, "twin": 0}

    def index_func(z):
        keys = z_src_keys(z)
        mapped = keys.map(_worker_src_to_b38)
        unlifted = mapped.isna()
        counts["unlifted"] = int(unlifted.sum())
        counts["twin"] = int(mapped.str.startswith("b37_duprep::", na=False).sum())
        mapped[unlifted] = "b37_unlifted::" + keys[unlifted]
        return pd.Index(mapped, name="variant")

    a, l, p, _ = _read_susie_rds(rds_path, z_path, cohort, susie_l, index_func=index_func)
    return region38, a, l, p, counts["unlifted"], counts["twin"]


def load_sim_multicausal(pheno: int, staged: dict, susie_l: int = 10, min_cohorts: int = None,
                         accounting: dict = None):
    """Yield `(region, alpha, lbf, pip)` for the multi-causal sims. Work list and
    completeness semantics as in `load_sim_1causal`: the meta-region bed names the regions,
    and by default any cohort-incomplete bed region raises."""
    cohorts, root = staged["cohorts"], staged["dir"]
    per_cohort = {}
    for cohort in cohorts:
        cp = _with_pheno(cohort, pheno)
        per_cohort[cohort] = {
            os.path.basename(f)[len(cp) + 1:-len(".susie.rds")]: (
                f, f"{root}/{cp}/{cp}.{os.path.basename(f)[len(cp) + 1:-len('.susie.rds')]}.z")
            for f in glob.glob(f"{root}/{cp}/*.susie.rds")}
    per_cohort = _restrict_to_bed(per_cohort, staged["bed_regions"], accounting)
    for region, present in _resolve_regions(per_cohort, min_cohorts, accounting):
        a_parts, l_parts, p_parts = [], [], []
        for cohort in present:
            rds, zf = per_cohort[cohort][region]
            a, l, p, _ = _read_susie_rds(rds, zf, cohort, susie_l)
            a_parts.append(a)
            l_parts.append(l)
            p_parts.append(p)
        alpha = reduce(_outer, a_parts)
        lbf = reduce(_outer, l_parts)[alpha.columns]
        pip = reduce(_outer, p_parts)
        _check_frames(region, alpha, lbf, pip)
        yield region, alpha, lbf, pip


# ---------------------------------------------------------------------------------------------
# Figure 4 -- FinnGen R12 + UKBB real data
# ---------------------------------------------------------------------------------------------
# FinnGen: one aggregated feather per endpoint, already carrying alpha1..10 AND
# lbf_variable1..10, indexed by rsid with a `region_grch38` column.
#
# UKBB: read from the per-region `.susie.rds` + `.z`, NOT from the aggregated
# `{pheno}.merged_with_z.feather`. That feather carries alpha_1..alpha_10 but NO lbf_variable
# columns, so coloc-susie cannot be scored from it. The per-region SuSiE fits do have
# lbf_variable (`run_susieR.R` was invoked with `--write-lbf-variable` and `--save-susie-obj`),
# so this is a re-read of existing output -- SuSiE is NOT re-run.
#
# Region identity: UKBB regions are GRCh37, FinnGen's are GRCh38. The per-endpoint
# `{fg}_38_to_37_mapping.tsv` supplies the correspondence; the `start - 1` shift reproduces the
# half-open/closed convention difference the production notebook applied. Getting this wrong
# does not raise, it just produces an empty or wrong join, so `load_real_data` asserts a
# non-trivial match rate.

def stage_real_data(fg_pheno: str, ukb_pheno: str, local_root: str) -> dict:
    dest = os.path.join(local_root, f"{fg_pheno}__{ukb_pheno}")
    os.makedirs(dest, exist_ok=True)
    fg_local = os.path.join(dest, f"{fg_pheno}.susie.all_snps.feather")
    _gsutil_cp_file(FG_FEATHER.format(fg=fg_pheno), fg_local)
    map_local = os.path.join(dest, f"{fg_pheno}_38_to_37_mapping.tsv")
    _gsutil_cp_file(REGION_MAP.format(fg=fg_pheno), map_local)
    # ONE completion marker covers both the .susie.rds and the .z pattern; the old separate
    # any-file-exists checks were the staging re-entrancy hole (2026-08-22 review finding #1).
    ukb_dest = os.path.join(dest, "ukb")
    _ensure_staged(ukb_dest, [f"{UKB_SUSIE_DIR.format(ukb=ukb_pheno)}/*.susie.rds",
                              f"{UKB_Z_DIR.format(ukb=ukb_pheno)}/*.z"])
    lift_local = os.path.join(dest, f"{ukb_pheno}.b37_to_b38.tsv.gz")
    _gsutil_cp_file(UKB_LIFTOVER_MAP.format(ukb=ukb_pheno), lift_local)
    return {"dir": dest, "fg_feather": fg_local, "region_map": map_local, "ukb_dir": ukb_dest,
            "liftover_map": lift_local, "fg_pheno": fg_pheno, "ukb_pheno": ukb_pheno}


def _region_37_to_38(map_path: str) -> dict:
    """GRCh37 region label -> GRCh38 region label, per the production notebook's construction."""
    m = pd.read_csv(map_path, delimiter=r"\s+")
    collapsed = (m.groupby(["GRCh38_region", "GRCh37_chr"])
                  .agg({"GRCh37_start": "min", "GRCh37_end": "max"}).reset_index())
    g37 = ("chr" + collapsed["GRCh37_chr"].astype(str) + "."
           + collapsed["GRCh37_start"].astype(str) + "-"
           + collapsed["GRCh37_end"].astype(str))
    g38 = collapsed["GRCh38_region"].str.replace(":", ".", regex=False)

    def shift(r):
        chrom, coords = r.split(".")
        start, end = coords.split("-")
        return f"{chrom}.{int(start) - 1}-{end}"

    if g37.duplicated().any():
        # dict(zip(...)) would silently keep the LAST mapping for a duplicated GRCh37 label.
        raise ValueError(f"{map_path}: duplicated GRCh37 region label(s) after collapsing: "
                         f"{sorted(g37[g37.duplicated()].unique())}")
    return dict(zip(g37, g38.map(shift)))


def load_real_data(staged: dict, susie_l: int = 10, cohort_names=("FINNGEN", "UKBB"),
                   accounting: dict = None):
    """Yield `(region, alpha, lbf, pip)` per GRCh38 region for the FinnGen+UKBB pair.

    Column renaming, and why it is not cosmetic
    -------------------------------------------
    FinnGen's feather names its components `alpha1..alpha10` and UKBB's `alpha_1..alpha_10`.
    Passed through untouched, `fastmap._cohort_names()` would read those as cohorts `""` and
    `"alpha"` -- which happens to keep A1 satisfied but makes the `note` provenance strings
    meaningless and breaks `finalize_region`'s backfill, which looks up `f"{cohort}_prob"`.
    They are renamed to `FINNGEN_alpha{i}` / `UKBB_alpha{i}` and `FINNGEN_prob` / `UKBB_prob`.

    Outer merge, not inner
    ----------------------
    The production notebook inner-joined FinnGen to UKBB and then bolted on a separate
    "rescue cohort-specific SNPs" pass afterwards. This adapter outer-joins instead, which is
    the same convention the simulations use and routes cohort-specific variants through the
    algorithm's own designed paths: `combine_pips`' coverage-scaled fill while merging, and
    `finalize_region`'s `filling_from_{cohort}` backfill for variants no selected component
    covers. The rescue pass becomes unnecessary rather than being reimplemented. NOTE this
    changes real-data handling relative to the shipped `fastmap.v04212026` run.

    Join key: lifted b38 positional ids, not rsid (Ran, 2026-08-24)
    ---------------------------------------------------------------
    Both cohorts are indexed by FinnGen-style composite ids `chr{c}_{pos}_{ref}_{alt}` on
    GRCh38: FinnGen by its own `variant` column, UKBB by translating each .z variant through
    the precomputed BCFtools/+liftover map (see `build_ukb_liftover_maps.py`). rsid keying
    is retired: it crashed on FinnGen's 4.5% unnamed variants, matched only 11.7% of indels,
    and mispaired multi-allelics whose rsid the two datasets assigned to different alt
    alleles (dbSNP b142-era vs b155). A UKBB variant absent from the map (liftover reject;
    ~0.001%) keeps a `b37_unlifted::`-prefixed id so it flows through as cohort-specific and
    is counted in `accounting`, never silently dropped.
    """
    fg_name, ukb_name = cohort_names
    fg_alpha = [f"alpha{i}" for i in range(1, susie_l + 1)]
    fg_lbf = [f"lbf_variable{i}" for i in range(1, susie_l + 1)]
    fg_cols = fg_alpha + fg_lbf + ["variant", "prob", "region_grch38"]
    # The feather carries 61 columns (incl. 7 string columns: chromosome, allele1/2, SNP,
    # rsid, ...) but only the ones above are ever read below. On HEIGHT_IRN (21.5M rows) a
    # full pd.read_feather materializes those strings as pandas object arrays and blows past
    # 30GB RSS, OOM-killing the whole cgroup. Restricting columns at read time keeps pyarrow
    # from ever building the unused string columns.
    available = pd.read_feather(staged["fg_feather"], columns=[]).columns  # cheap: no data read
    missing = [c for c in fg_cols if c not in available]
    if missing:
        raise ValueError(f"{staged['fg_feather']}: missing columns {missing}")
    fg = pd.read_feather(staged["fg_feather"], columns=fg_cols)
    if not fg["variant"].str.match(r"chr[0-9X]+_[0-9]+_[A-Z]+_[A-Z]+").all():
        bad = fg.loc[~fg["variant"].str.match(r"chr[0-9X]+_[0-9]+_[A-Z]+_[A-Z]+"), "variant"]
        raise ValueError(f"{staged['fg_feather']}: unexpected variant id format, e.g. "
                         f"{bad.head(3).tolist()}")

    lift = pd.read_csv(staged["liftover_map"], sep="\t")
    b38_id = ("chr" + lift["b38_chrom"].str.replace("chr", "", regex=False) + "_"
              + lift["b38_pos"].astype(str) + "_" + lift["b38_ref"] + "_" + lift["b38_alt"])
    # Representation twins (two b37 encodings of one indel; `b38_dup_rank` > 0) must not
    # both carry the b38 id or the outer merge fans out; the non-canonical twin becomes a
    # cohort-specific row instead.
    twin = lift["b38_dup_rank"] > 0
    b38_id[twin] = "b37_duprep::" + lift.loc[twin, "src_key"]
    src_to_b38 = dict(zip(lift["src_key"], b38_id))
    n_unlifted_total = n_twin_total = 0
    # (the per-region b38 index build now happens inside `_ukb_region_worker`, in parallel
    # worker processes -- see the `pool.map` call below.)

    r37_to_38 = _region_37_to_38(staged["region_map"])
    ukb_prefix = f"UKB.{staged['ukb_pheno']}."
    ukb_files, unmapped = {}, []
    for rds in sorted(glob.glob(f"{staged['ukb_dir']}/*.susie.rds")):
        base = os.path.basename(rds)[:-len(".susie.rds")]
        if not base.startswith(ukb_prefix):
            raise ValueError(f"unexpected UKBB file name {base!r} (want {ukb_prefix}<region>)")
        region37 = base[len(ukb_prefix):]
        zf = f"{staged['ukb_dir']}/{base}.z"
        if not os.path.exists(zf):
            # The staging marker guarantees the local copy matches GCS, so a missing .z means
            # the SOURCE lacks it -- a data-integrity problem, not a staging hiccup. The old
            # warn-and-skip would silently drop the region from Figure 4.
            raise FileNotFoundError(
                f"{staged['ukb_pheno']}: UKBB region {region37} has a .susie.rds but no .z "
                f"({zf}). The source data is inconsistent; refusing to skip silently.")
        region38 = r37_to_38.get(region37)
        if region38 is None:
            # Expected for UKBB regions outside this FG endpoint's region universe (the map
            # is FG-derived); they could never pair with FinnGen anyway. Counted, not dropped
            # without trace.
            unmapped.append(region37)
            continue
        if region38 in ukb_files:
            # Two GRCh37 labels can map to one GRCh38 region when its liftover spans two
            # GRCh37 chromosomes. Silently keeping the last glob-order file would pair FinnGen
            # with an arbitrary partial UKBB region (2026-08-22 review finding #6).
            raise ValueError(
                f"{staged['ukb_pheno']}: GRCh38 region {region38} maps from TWO UKBB fits "
                f"({ukb_files[region38][2]} and {region37}). Ambiguous pairing -- decide "
                "explicitly how to handle split-liftover regions before running this "
                "endpoint.")
        ukb_files[region38] = (rds, zf, region37)

    mapped = len(ukb_files)
    if mapped == 0:
        raise AssertionError(
            f"{staged['ukb_pheno']}: zero UKBB regions mapped to GRCh38. The region-label "
            f"convention has changed -- check {staged['region_map']} and the start-1 shift.")
    if unmapped:
        print(f"  {len(unmapped)} UKBB region(s) not in this endpoint's 38<->37 map "
              f"(outside FinnGen's region universe): {unmapped[:3]}"
              f"{' ...' if len(unmapped) > 3 else ''}")

    fg_by_region = dict(tuple(fg.groupby("region_grch38")))
    shared = sorted(set(fg_by_region) & set(ukb_files))
    print(f"  FinnGen regions {len(fg_by_region)}, UKBB regions {mapped}, shared {len(shared)}")
    if not shared:
        raise AssertionError("no region shared between FinnGen and UKBB after mapping")
    if accounting is not None:
        accounting["fg_regions"] = len(fg_by_region)
        accounting["ukb_regions_mapped"] = mapped
        accounting["ukb_regions_unmapped"] = unmapped
        accounting["regions_processed"] = len(shared)
        accounting["fg_only_regions"] = sorted(set(fg_by_region) - set(ukb_files))
        accounting["ukb_only_regions"] = sorted(set(ukb_files) - set(fg_by_region))

    # The per-region UKB read (readRDS via rpy2 + z read + b38 index build) is ~0.5s/region
    # and was the dominant serial cost on a genome-wide trait like BMI (thousands of shared
    # regions, single CPU core pinned the whole run). Fan it out across a `spawn` process
    # pool -- `spawn`, not the default `fork`, because rpy2's embedded R is not fork-safe and
    # a forked worker can inherit a half-initialized R state from the parent. `pool.map`
    # preserves `shared`'s order in its output, so the per-region body below (FG prep, dup
    # checks, merge, `_check_frames`, yield) is otherwise UNCHANGED and still fully ordered.
    n_workers = min(len(shared), os.cpu_count() or 4, 8)
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx,
                             initializer=_init_ukb_worker, initargs=(src_to_b38,)) as pool:
        ukb_reads = pool.map(
            _ukb_region_worker, shared,
            [ukb_files[r][0] for r in shared], [ukb_files[r][1] for r in shared],
            [ukb_name] * len(shared), [susie_l] * len(shared))
        for region38, a_ukb, l_ukb, p_ukb, n_unlifted, n_twin in ukb_reads:
            n_unlifted_total += n_unlifted
            n_twin_total += n_twin

            g = fg_by_region[region38]
            idx = pd.Index(g["variant"], name="variant")
            if idx.duplicated().any():
                # Same variant appearing twice inside one region would fan out on merge.
                # FinnGen's feather is per-region unique in the runs seen; fail rather than
                # silently inflate.
                raise ValueError(f"{region38}: duplicated FinnGen variant id")
            a_fg = pd.DataFrame(g[fg_alpha].to_numpy(dtype=float), index=idx,
                                columns=[f"{fg_name}_alpha{i}" for i in range(1, susie_l + 1)])
            l_fg = pd.DataFrame(g[fg_lbf].to_numpy(dtype=float), index=idx,
                                columns=[f"{fg_name}_alpha{i}" for i in range(1, susie_l + 1)])
            p_fg = pd.DataFrame({f"{fg_name}_prob": g["prob"].to_numpy(dtype=float)}, index=idx)

            if a_ukb.index.duplicated().any():
                # FG's side is already checked above; UKB's was never checked because the old
                # three-way `merge(how="outer")` would fan out visibly and `_check_frames`
                # would catch it downstream. Union+reindex below does NOT fan out on a
                # duplicated key (it just picks one), so the check has to be explicit here.
                raise ValueError(f"{region38}: duplicated UKBB variant id")

            # alpha/lbf/pip's FG and UKB frames all share the SAME two index sets, so three
            # separate `merge(how="outer")` calls redo the identical join alignment three
            # times. Both indices are now known-unique, so union-once + reindex is exactly
            # equivalent to `how="outer"` merge (no possible fan-out) at ~1/3 the cost.
            union_idx = a_fg.index.union(a_ukb.index)
            alpha = pd.concat([a_fg.reindex(union_idx), a_ukb.reindex(union_idx)], axis=1)
            lbf = pd.concat([l_fg.reindex(union_idx), l_ukb.reindex(union_idx)],
                            axis=1)[alpha.columns]
            pip = pd.concat([p_fg.reindex(union_idx), p_ukb.reindex(union_idx)], axis=1)
            _check_frames(region38, alpha, lbf, pip)
            if accounting is not None:
                accounting["ukb_variants_unlifted"] = n_unlifted_total
                accounting["ukb_variants_duplicate_representation"] = n_twin_total
            yield region38, alpha, lbf, pip


# ---------------------------------------------------------------------------------------------
# dispatch
# ---------------------------------------------------------------------------------------------

STAGERS = {
    "sim_1causal": stage_sim_1causal,
    "sim_multicausal": stage_sim_multicausal,
}
LOADERS = {
    "sim_1causal": load_sim_1causal,
    "sim_multicausal": load_sim_multicausal,
}
