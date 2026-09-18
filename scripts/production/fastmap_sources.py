#!/usr/bin/env python3
"""Input adapters: per-cohort fine-mapping output -> the three frames FastMap consumes.

One adapter per manuscript figure. All three return the SAME contract:

    alpha  per-variant single-effect posteriors, columns "{cohort}_alpha{i}"   -> region_df
    lbf    per-variant natural-log Bayes factors, SAME columns, SAME order     -> lbf_df
    pip    per-variant marginal PIP, one "{cohort}_prob" column per cohort     -> pips_df

All are indexed by variant id, outer-merged across cohorts, and column-aligned.

  sim_1causal      Figure 2   per-cohort ABF `.abf.snp`, 1 component per cohort, 10 cohorts
  sim_multicausal  Figure 3   per-cohort SuSiE `.susie.rds` + `.z`, 10 components, 2 cohorts
  real_data        Figure 4   FinnGen aggregated feather + UKBB per-region `.susie.rds` + `.z`
"""
from __future__ import annotations

import glob
import json
import multiprocessing as mp
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from functools import reduce

import numpy as np
import pandas as pd

#: Root GCS bucket for simulation/real-data inputs.
DATA_ROOT = os.environ.get("FASTMAP_DATA_ROOT", "gs://YOUR-BUCKET")

SIMS = f"{DATA_ROOT}/simulations/sims"
CONFIGS = f"{DATA_ROOT}/simulations/scripts/configs"

# Figure 2 / Figure 3 cohort sets.
SIM_SOURCES = {
    "sim_1causal": {
        "config_tsv": f"{CONFIGS}/All_random_configs10_Omni25_noTOPMed.tsv",
        "config_index": 1,
        "assoc_dir": f"{SIMS}/assoc_1causal",
        "meta_dir": f"{SIMS}/meta_analysis_1causal",
        "kind": "abf",
        "susie_l": 1,
    },
    "sim_multicausal": {
        "config_tsv": f"{CONFIGS}/All_random_configs2_Omni25_noTOPMed.tsv",
        "config_index": 1,
        "assoc_dir": f"{SIMS}/assoc_3causal_27h2g_fixed",
        "meta_dir": f"{SIMS}/meta_analysis_3causal_27h2g_fixed",
        "kind": "susie",
        "susie_l": 10,
    },
}

# Figure 4: FinnGen R12 endpoint <-> UKBB phenotype mapping.
REAL_DATA_PHENOS = {
    "I9_CHD": "CAD", "I9_MI_STRICT": "MI", "I9_AF": "AFib", "C_STROKE": "IS",
    "E4_HYTHY_AI_STRICT": "Hypothyroidism", "T2D": "T2D", "T1D": "T1D",
    "K11_CD_STRICT2": "Crohn", "RHEUMA_SEROPOS_WIDE": "RA",
    "H7_CATARACTSENILE": "Cataract",
    "C3_BRONCHUS_LUNG_EXALLC": "LuC", "C3_SKIN_EXALLC": "SkC",
    "G6_AD_WIDE": "Alzheimer_LTFH",
    "BMI_IRN": "BMI", "HEIGHT_IRN": "Height",
}

# Present in GCS but missing the necessary FINNGEN feather output.
REAL_DATA_CANDIDATES_NO_FG_FEATHER = {
    "J10_ASTHMA_EXMORE": "Asthma_Combined", "H7_GLAUCOMA": "Glaucoma_Combined",
    "C3_PROSTATE_EXALLC": "PrC_M",
}

FG_FEATHER = f"{DATA_ROOT}/results/FINNGEN_meta_regions/UKBB_FINNGEN_R12/June_2026/{{fg}}.susie.all_snps.feather"
UKB_SUSIE_DIR = f"{DATA_ROOT}/results/UKBB_meta_regions/UKBB_FINNGEN_R12/{{ukb}}"
UKB_Z_DIR = f"{DATA_ROOT}/data/FINNGEN_UKBB_meta/R12/finemap_inputs/{{ukb}}"
REGION_MAP = f"{DATA_ROOT}/data/FINNGEN_UKBB_meta/R12/meta_generated_bed_files/June_2026/{{fg}}_38_to_37_mapping.tsv"
#: Per-variant b37->b38 liftover maps.
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
    """The liftover map's `src_key` reproduced from .z columns."""
    return (z["chromosome"].map(norm_chrom_b37) + ":" + z["position"].astype(str) + ":"
            + z["allele1"].str.upper() + ":" + z["allele2"].str.upper())


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------

def _run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def _gsutil_cp(patterns, dest, quiet=True):
    """Multi-file copy into `dest` using `gcloud storage cp`."""
    os.makedirs(dest, exist_ok=True)
    cmd = ["gcloud", "storage", "cp"] + (["-q"] if quiet else []) + list(patterns) + [dest]
    _run(cmd)


STAGING_MARKER = ".staging_complete"


def _ensure_staged(dest: str, patterns: list) -> None:
    """Download `patterns` into `dest` unless a prior COMPLETE download is recorded there."""
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
    """Download one file atomically: fetch to a temp name, os.replace into place."""
    if os.path.exists(local):
        return
    tmp = local + ".staging_tmp"
    _run(["gsutil", "-q", "cp", src, tmp])
    os.replace(tmp, local)


def _outer(a, b):
    return pd.merge(a, b, how="outer", left_index=True, right_index=True)


def cohorts_for_config(config_tsv: str, config_index: int, local_root: str) -> list:
    """Cohort names for one config row. `config_index` is 1-BASED."""
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
    """Download and parse one phenotype's meta-analysis region BED -> region-name list."""
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
    """Keep only the bed's regions in each cohort's file map; extra files are ignored."""
    wanted = set(bed_regions)
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
        print(f"  {n_extra} per-cohort region file(s) not named by the meta bed are ignored; "
              f"counts per cohort are in the manifest")
    return {c: {r: f for r, f in files.items() if r in wanted}
            for c, files in per_cohort.items()}


def _resolve_regions(per_cohort: dict, min_cohorts, accounting) -> list:
    """The `(region, present_cohorts)` work list, checking region presence across cohorts."""
    cohorts = list(per_cohort)
    empty = [c for c in cohorts if not per_cohort[c]]
    if empty:
        raise RuntimeError(
            f"cohort(s) with ZERO regions on disk: {empty} -- staging is incomplete or the "
            "source is missing.")
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
            f"{' ...' if len(skipped) > 5 else ''}). Pass min_cohorts=<k> "
            "explicitly if partial-cohort regions are intended.")
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
    if alpha.isna().all(axis=1).any():
        raise AssertionError(f"{region}: variants with no alpha in any cohort")


# ---------------------------------------------------------------------------------------------
# Figure 2 -- single-causal simulations, per-cohort ABF
# ---------------------------------------------------------------------------------------------

def stage_sim_1causal(pheno: int, local_root: str, config_index: int = 1) -> dict:
    """Download one phenotype's per-cohort ABF output plus the meta-region bed."""
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
    """Yield `(region, alpha, lbf, pip)` for the 1-causal sims."""
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
    """Download one phenotype's per-cohort SuSiE output plus the meta-region bed."""
    spec = SIM_SOURCES["sim_multicausal"]
    cohorts = cohorts_for_config(spec["config_tsv"], config_index, local_root)
    dest_root = os.path.join(local_root, f"pheno{pheno}")
    for cohort in cohorts:
        cp = _with_pheno(cohort, pheno)
        src = f"{spec['assoc_dir']}/pheno{pheno}/{cp}"
        _ensure_staged(os.path.join(dest_root, cp), [f"{src}/*.susie.rds", f"{src}/*.z"])
    regions = _stage_regions_bed("sim_multicausal", pheno, config_index, local_root)
    return {"cohorts": cohorts, "dir": dest_root, "bed_regions": regions}


def _read_susie_rds(rds_path: str, z_path: str, cohort: str, susie_l: int,
                    variant_col: str = "rsid", index_func=None):
    """One cohort x region SuSiE fit -> (alpha, lbf, pip) frames, plus variant positions."""
    import rpy2.robjects as ro

    r = ro.r["readRDS"](rds_path)
    z = pd.read_csv(z_path, delimiter=r"\s+")
    if variant_col not in z.columns:
        raise ValueError(f"{z_path}: no '{variant_col}' column (has {list(z.columns)[:8]})")

    # `.copy()` is required because rpy2 returns a view into R memory.
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
    p = pd.DataFrame({f"{cohort}_prob": pip}, index=idx)
    pos = (pd.DataFrame({"position": z["position"].values}, index=idx)
           if "position" in z.columns else None)
    return a, l, p, pos


_worker_src_to_b38 = None


def _init_ukb_worker(src_to_b38: dict) -> None:
    """Initialize a worker process with the global liftover map."""
    global _worker_src_to_b38
    _worker_src_to_b38 = src_to_b38


def _ukb_region_worker(region38: str, rds_path: str, z_path: str, cohort: str, susie_l: int):
    """Worker logic for reading and mapping one UKBB region."""
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
    """Yield `(region, alpha, lbf, pip)` for the multi-causal sims."""
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

def stage_real_data(fg_pheno: str, ukb_pheno: str, local_root: str) -> dict:
    dest = os.path.join(local_root, f"{fg_pheno}__{ukb_pheno}")
    os.makedirs(dest, exist_ok=True)
    fg_local = os.path.join(dest, f"{fg_pheno}.susie.all_snps.feather")
    _gsutil_cp_file(FG_FEATHER.format(fg=fg_pheno), fg_local)
    map_local = os.path.join(dest, f"{fg_pheno}_38_to_37_mapping.tsv")
    _gsutil_cp_file(REGION_MAP.format(fg=fg_pheno), map_local)
    ukb_dest = os.path.join(dest, "ukb")
    _ensure_staged(ukb_dest, [f"{UKB_SUSIE_DIR.format(ukb=ukb_pheno)}/*.susie.rds",
                              f"{UKB_Z_DIR.format(ukb=ukb_pheno)}/*.z"])
    lift_local = os.path.join(dest, f"{ukb_pheno}.b37_to_b38.tsv.gz")
    _gsutil_cp_file(UKB_LIFTOVER_MAP.format(ukb=ukb_pheno), lift_local)
    return {"dir": dest, "fg_feather": fg_local, "region_map": map_local, "ukb_dir": ukb_dest,
            "liftover_map": lift_local, "fg_pheno": fg_pheno, "ukb_pheno": ukb_pheno}


def _region_37_to_38(map_path: str) -> dict:
    """GRCh37 region label -> GRCh38 region label mapping."""
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
        raise ValueError(f"{map_path}: duplicated GRCh37 region label(s) after collapsing: "
                         f"{sorted(g37[g37.duplicated()].unique())}")
    return dict(zip(g37, g38.map(shift)))


def load_real_data(staged: dict, susie_l: int = 10, cohort_names=("FINNGEN", "UKBB"),
                   accounting: dict = None):
    """Yield `(region, alpha, lbf, pip)` per GRCh38 region for the FinnGen+UKBB pair."""
    fg_name, ukb_name = cohort_names
    fg_alpha = [f"alpha{i}" for i in range(1, susie_l + 1)]
    fg_lbf = [f"lbf_variable{i}" for i in range(1, susie_l + 1)]
    fg_cols = fg_alpha + fg_lbf + ["variant", "prob", "region_grch38"]

    available = pd.read_feather(staged["fg_feather"], columns=[]).columns
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

    twin = lift["b38_dup_rank"] > 0
    b38_id[twin] = "b37_duprep::" + lift.loc[twin, "src_key"]
    src_to_b38 = dict(zip(lift["src_key"], b38_id))
    n_unlifted_total = n_twin_total = 0

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
            raise FileNotFoundError(
                f"{staged['ukb_pheno']}: UKBB region {region37} has a .susie.rds but no .z "
                f"({zf}).")
        region38 = r37_to_38.get(region37)
        if region38 is None:
            unmapped.append(region37)
            continue
        if region38 in ukb_files:
            raise ValueError(
                f"{staged['ukb_pheno']}: GRCh38 region {region38} maps from TWO UKBB fits "
                f"({ukb_files[region38][2]} and {region37}). Ambiguous pairing.")
        ukb_files[region38] = (rds, zf, region37)

    mapped = len(ukb_files)
    if mapped == 0:
        raise AssertionError(
            f"{staged['ukb_pheno']}: zero UKBB regions mapped to GRCh38.")
    if unmapped:
        print(f"  {len(unmapped)} UKBB region(s) not in this endpoint's 38<->37 map: "
              f"{unmapped[:3]}{' ...' if len(unmapped) > 3 else ''}")

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
                raise ValueError(f"{region38}: duplicated FinnGen variant id")
            a_fg = pd.DataFrame(g[fg_alpha].to_numpy(dtype=float), index=idx,
                                columns=[f"{fg_name}_alpha{i}" for i in range(1, susie_l + 1)])
            l_fg = pd.DataFrame(g[fg_lbf].to_numpy(dtype=float), index=idx,
                                columns=[f"{fg_name}_alpha{i}" for i in range(1, susie_l + 1)])
            p_fg = pd.DataFrame({f"{fg_name}_prob": g["prob"].to_numpy(dtype=float)}, index=idx)

            if a_ukb.index.duplicated().any():
                raise ValueError(f"{region38}: duplicated UKBB variant id")

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
