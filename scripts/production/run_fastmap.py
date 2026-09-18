#!/usr/bin/env python3
"""Production FastMap runner -- ONE entry point, ONE operating point, three data sources.

This is the script that regenerates the FastMap arm of all three results figures at the
part-2-selected operating point (`production_settings.PRODUCTION` by default; see
`fastmap/production_settings.py` for its full provenance). It replaces four divergent invocation paths that had drifted:

  * `scripts/simulation_pipeline/fastmap_input_processing.py` (sims, coupled algorithm, CLI
    flags for the now-INAPPLICABLE production-only parameters)
  * `jupyter_scripts/Real_Data_Application/FastMap_for_FinnGen_UKBB_meta_v04212026.ipynb`
    (real data -- an inline REIMPLEMENTATION of the algorithm, with the phenotype hard-coded in
    a cell and its own separate "rescue cohort-specific SNPs" pass)
  * the six legacy copies under `scripts/fastmap_algorithm/`
  * the part-2 sweep driver

Every one of them called a different copy of the algorithm or a different parameter set. Here
the algorithm is imported from `fastmap/fastmap.py`, the parameters come from
`fastmap/production_settings.py`, and the per-source plumbing lives in `fastmap_sources.py`.
Nothing in this file hard-codes a hyperparameter.

Usage
-----
  # Figure 2 -- single-causal sims, one phenotype
  python3 scripts/production/run_fastmap.py --source sim_1causal --pheno 1

  # Figure 3 -- multi-causal sims, a range
  python3 scripts/production/run_fastmap.py --source sim_multicausal --phenos 1-100

  # Figure 4 -- real data, all 18 FinnGen/UKBB endpoint pairs
  python3 scripts/production/run_fastmap.py --source real_data --all-phenos

  # sensitivity: same code, a different vetted operating point (v2 grid ids; see ALTERNATIVES)
  python3 scripts/production/run_fastmap.py --source sim_multicausal --pheno 1 --setting 35

  # custom parameters: the three decision knobs default to the recommended values
  # (PP.H4 threshold 0.1, overlap_min 0.5, L 10) and can each be overridden; outputs are
  # tagged with the deviation (here '.setting107_L1') and the manifest records it
  python3 scripts/production/run_fastmap.py --source sim_1causal --phenos 1-400 --L 1

Outputs, per phenotype
----------------------
  {prefix}.fastmap.snp.feather   one row per (variant, region): `prob_fastmap`, `note`,
                                 `comp_i_val` / `comp_i_name`, `region`
  {prefix}.runtime.tsv           per-region cpu / wall / peak-RSS / n_snps / n_components.
                                 FastMap had NO per-region runtime record before this; Figure 3
                                 panel e is a CPU-time comparison, and the numbers it currently
                                 rests on came from resource logs that were partly clobbered by
                                 the Stage 5/6 `.resource.log` name collision.
  {prefix}.manifest.json         setting, conventions, input paths, code hashes, versions.
                                 Written so that "which parameters produced this file" is a
                                 property of the file, not of somebody's memory -- the exact
                                 failure that forced the part-2 re-analysis.

`PYTHONHASHSEED` is pinned to 0 and enforced. The legacy coupled metric iterated a Python set
of variant-id strings, making its output depend on per-process string hashing; the current
similarity paths are deterministic, but the run stays pinned so that any future comparison
against a legacy artifact is meaningful.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import fastmap_sources as src                                    # noqa: E402
from fastmap.fastmap import combine_region                       # noqa: E402
from fastmap.production_settings import (                         # noqa: E402
    PRODUCTION, assert_matches_grid, get_setting)


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _parse_phenos(spec: str) -> list:
    """"1", "1-100", "1,5,7-9" -> a sorted list of ints."""
    out = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-")
            out.update(range(int(lo), int(hi) + 1))
        else:
            out.add(int(part))
    return sorted(out)


def run_one(regions_iter, setting, out_prefix: str, out_dir: str) -> dict:
    """Run FastMap over one phenotype's regions; write the feather and the runtime table.

    Regions are streamed rather than collected into a dict first (which is what
    `fastmap.fastmap()` does) so that a 900-region real-data endpoint does not need every
    region's frames resident at once, and so a per-region failure names its region.
    """
    kwargs = setting.as_kwargs()
    per_region, runtime_rows = [], []
    for region, alpha, lbf, pip in regions_iter:
        t_cpu = time.process_time()
        t_wall = time.perf_counter()
        rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        try:
            res = combine_region(region, alpha, pip, lbf_df=lbf, **kwargs)
        except Exception as exc:                       # noqa: BLE001 -- want the region name
            raise RuntimeError(f"region {region} failed: {exc}") from exc
        cpu = time.process_time() - t_cpu
        wall = time.perf_counter() - t_wall
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        n_comp = sum(1 for c in res.columns if c.endswith("_name"))
        per_region.append(res)
        runtime_rows.append({
            "region": region, "n_snps": len(alpha), "n_components_in": alpha.shape[1],
            "n_components_out": n_comp, "cpu_seconds": cpu, "wall_seconds": wall,
            # ru_maxrss is a process high-water mark, so it is monotone and NOT a per-region
            # peak. Both columns are recorded and named accordingly; do not read
            # `peak_rss_kb_highwater` as this region's cost.
            "peak_rss_kb_highwater": rss, "peak_rss_kb_grew": max(0, rss - rss_before),
        })
        print(f"    {region}: {len(alpha)} variants, {alpha.shape[1]} -> {n_comp} components, "
              f"{cpu:.2f} cpu-s", flush=True)

    if not per_region:
        raise RuntimeError(f"{out_prefix}: no regions produced output")

    result = pd.concat(per_region, axis=0)
    os.makedirs(out_dir, exist_ok=True)
    feather = os.path.join(out_dir, f"{out_prefix}.fastmap.snp.feather")
    result.reset_index().to_feather(feather)
    runtime = pd.DataFrame(runtime_rows)
    runtime.to_csv(os.path.join(out_dir, f"{out_prefix}.runtime.tsv"), sep="\t", index=False)

    notes = result["note"].astype(str)
    return {
        "output_feather": feather,
        "n_regions": len(per_region),
        "n_rows": int(len(result)),
        "n_variants_distinct": int(result.index.nunique()),
        "total_cpu_seconds": float(runtime["cpu_seconds"].sum()),
        "total_wall_seconds": float(runtime["wall_seconds"].sum()),
        "peak_rss_kb": int(runtime["peak_rss_kb_highwater"].max()),
        # `filling_from_*` = variants no SELECTED component covered, backfilled from the best
        # marginal PIP among unused cohorts. Reported because Figure 4's headline must not
        # credit FastMap for variants whose PIP was copied rather than combined.
        "n_backfilled": int(notes.str.startswith("filling_from_").sum()),
        "frac_backfilled": float(notes.str.startswith("filling_from_").mean()),
        "prob_fastmap_na": int(result["prob_fastmap"].isna().sum()),
        # `lbf_fastmap` exists only on the coloc-susie path (the coupled path never loads lbf).
        "lbf_fastmap_reported": (int(result["lbf_fastmap"].notna().sum())
                                 if "lbf_fastmap" in result.columns else None),
        "lbf_fastmap_backfilled": (
            int(result.loc[notes.str.startswith("filling_from_"), "lbf_fastmap"].notna().sum())
            if "lbf_fastmap" in result.columns else None),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True,
                    choices=["sim_1causal", "sim_multicausal", "real_data"])
    ap.add_argument("--pheno", type=int, help="single simulation phenotype")
    ap.add_argument("--phenos", type=str, help='simulation phenotypes, e.g. "1-100" or "1,5,7-9"')
    ap.add_argument("--fg-pheno", type=str, help="FinnGen endpoint (real_data)")
    ap.add_argument("--ukb-pheno", type=str, help="UKBB phenotype (real_data)")
    ap.add_argument("--all-phenos", action="store_true",
                    help="real_data: every pair in fastmap_sources.REAL_DATA_PHENOS")
    ap.add_argument("--setting", type=int, default=PRODUCTION.setting_id,
                    help=f"vetted operating point id (default {PRODUCTION.setting_id} "
                         "= production)")
    ap.add_argument("--threshold", type=float, default=None,
                    help=f"PP.H4 similarity threshold override "
                         f"(default {PRODUCTION.threshold}, the recommended value)")
    ap.add_argument("--overlap-min", type=float, default=None,
                    help=f"pair-admissibility overlap trim override "
                         f"(default {PRODUCTION.overlap_min}, the recommended value)")
    ap.add_argument("--L", type=int, default=None,
                    help=f"max components in the combined output override "
                         f"(default {PRODUCTION.L}, the recommended value)")
    ap.add_argument("--config-index", type=int, default=None,
                    help="1-based config row; defaults to the source's production value")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--stage-dir", default=None, help="where inputs are downloaded")
    ap.add_argument("--upload-to", default=None, help="gs:// prefix to copy outputs to")
    ap.add_argument("--keep-inputs", action="store_true",
                    help="do not delete staged inputs after each phenotype")
    ap.add_argument("--redo", action="store_true",
                    help="recompute even if a manifest done-marker exists")
    a = ap.parse_args()

    if os.environ.get("PYTHONHASHSEED") != "0":
        # Re-exec rather than warn: a run that silently used a random hash seed is not
        # reproducible, and this script's whole purpose is to make the run auditable.
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable] + sys.argv)

    setting = get_setting(a.setting).customize(
        threshold=a.threshold, overlap_min=a.overlap_min, L=a.L)
    assert_matches_grid(setting)
    print(f"FastMap production run -- {setting.describe()}")
    if setting.overrides:
        print(f"  NOTE: this is a CUSTOMIZED run -- "
              + ", ".join(f"{n}={v}" for n, v in setting.overrides)
              + f" overriding setting {setting.setting_id}; outputs are tagged "
              + f"'setting{setting.file_tag}'.")
    elif setting.setting_id != PRODUCTION.setting_id:
        print(f"  NOTE: this is a SENSITIVITY run, not the production operating point "
              f"({PRODUCTION.setting_id}).")

    stage_dir = a.stage_dir or os.path.join("/home/jupyter", f"fastmap_stage_{a.source}")
    out_dir = a.out_dir or os.path.join(ROOT, "results", "production",
                                        f"{a.source}_setting{setting.file_tag}")
    os.makedirs(out_dir, exist_ok=True)

    # ---- build the work list -----------------------------------------------------------
    if a.source == "real_data":
        if a.all_phenos:
            jobs = sorted(src.REAL_DATA_PHENOS.items())
        elif a.fg_pheno and a.ukb_pheno:
            jobs = [(a.fg_pheno, a.ukb_pheno)]
        elif a.fg_pheno:
            jobs = [(a.fg_pheno, src.REAL_DATA_PHENOS[a.fg_pheno])]
        else:
            ap.error("real_data needs --fg-pheno (and optionally --ukb-pheno) or --all-phenos")
    else:
        if a.pheno is not None:
            jobs = [a.pheno]
        elif a.phenos:
            jobs = _parse_phenos(a.phenos)
        else:
            ap.error(f"{a.source} needs --pheno or --phenos")

    spec = src.SIM_SOURCES.get(a.source)
    config_index = a.config_index if a.config_index is not None else (
        spec["config_index"] if spec else None)

    code_hashes = {
        rel: _sha256(os.path.join(ROOT, rel)) for rel in [
            "fastmap/fastmap.py", "fastmap/coloc.py", "fastmap/production_settings.py",
            "scripts/production/fastmap_sources.py", "scripts/production/run_fastmap.py"]}

    failures = []
    for job in jobs:
        # ---- stage + load -------------------------------------------------------------
        if a.source == "real_data":
            fg_pheno, ukb_pheno = job
            label = f"{fg_pheno}/{ukb_pheno}"
            out_prefix = (f"FG_R12.UKBB.fastmap.setting{setting.file_tag}."
                          f"{fg_pheno}.{ukb_pheno}")
            inputs_desc = {"fg_pheno": fg_pheno, "ukb_pheno": ukb_pheno,
                           "fg_feather": src.FG_FEATHER.format(fg=fg_pheno),
                           "ukb_susie_dir": src.UKB_SUSIE_DIR.format(ukb=ukb_pheno),
                           "ukb_z_dir": src.UKB_Z_DIR.format(ukb=ukb_pheno),
                           "region_map": src.REGION_MAP.format(fg=fg_pheno)}
        else:
            pheno = job
            label = f"pheno{pheno}"
            config_name = os.path.basename(spec["config_tsv"]).replace(".tsv", "")
            out_prefix = (f"{config_name}.config{config_index}.pheno{pheno}."
                          f"setting{setting.file_tag}")
            inputs_desc = {"pheno": pheno, "config_tsv": spec["config_tsv"],
                           "config_index": config_index,
                           "assoc_dir": f"{spec['assoc_dir']}/pheno{pheno}"}

        done_marker = os.path.join(out_dir, f"{out_prefix}.manifest.json")
        if os.path.exists(done_marker) and not a.redo:
            print(f"\n[{label}] already done ({done_marker}); skipping (--redo to recompute)")
            continue

        print(f"\n[{label}] staging inputs into {stage_dir}", flush=True)
        t0 = time.time()
        # Filled by the loader: per-cohort region counts and every region it did NOT process,
        # so the manifest records exactly what went in -- a run whose inputs were quietly
        # short can no longer look identical to a complete one.
        accounting = {}
        try:
            if a.source == "real_data":
                staged = src.stage_real_data(fg_pheno, ukb_pheno, stage_dir)
                regions = src.load_real_data(staged, accounting=accounting)
            elif a.source == "sim_1causal":
                staged = src.stage_sim_1causal(pheno, stage_dir, config_index)
                regions = src.load_sim_1causal(pheno, staged, accounting=accounting)
            else:
                staged = src.stage_sim_multicausal(pheno, stage_dir, config_index)
                regions = src.load_sim_multicausal(
                    pheno, staged, susie_l=spec["susie_l"], accounting=accounting)
            stats = run_one(regions, setting, out_prefix, out_dir)
            if accounting.get("regions_processed") not in (None, stats["n_regions"]):
                raise RuntimeError(
                    f"loader planned {accounting['regions_processed']} regions but "
                    f"{stats['n_regions']} produced output -- refusing to write a manifest "
                    "for a partially processed phenotype")
        except Exception as exc:                       # noqa: BLE001
            print(f"[{label}] FAILED: {exc}", file=sys.stderr)
            failures.append((label, str(exc)))
            continue

        manifest = {
            "source": a.source,
            "label": label,
            "setting": setting._asdict(),
            "inapplicable_parameters": {
                "names": ["max_div", "max_threshold", "initial_reverse_max_threshold",
                          "jaccard_threshold", "pip_threshold"],
                "why": ("these belonged to the coupled similarity='production' algorithm, "
                        "REMOVED from fastmap.py on 2026-08-21; before removal they were "
                        "verified INAPPLICABLE to the single-score paths (bit-identical "
                        "output over 27 combinations). Methods text must say INAPPLICABLE, "
                        "not 'unswept'."),
            },
            "conventions": {
                "merge_rule": (
                    "log-BFs add: lbf_merged = sum_k lbf_k, equivalently alpha_merged is the "
                    "elementwise product renormalised. Derived (section 3 of "
                    "documentation/coloc_susie_merged_lbf.pdf), not a convention, and NOT what "
                    "C1/C2 name. A variant measured on every side is unaffected by any fill."),
                "merged_lbf_missing_entry_fill": (
                    "C2 (the only behaviour since 2026-08-21, when Ran removed the C1 option "
                    "from the code) -- lbf := log(S_k/n_k), the log arithmetic-mean Bayes "
                    "factor over the variants that cohort measured (posterior-neutral, and "
                    "section 4 proves it unique). This is the SAME value combine_pips imputes "
                    "on the alpha side of the same merge, so a merged component's alpha and "
                    "lbf are one object: alpha == softmax(lbf), verified to 1.9e-15 "
                    "inductively across rounds. It cannot affect the score of a pair of two "
                    "original components (coloc_susie_matrix restricts to variants shared by "
                    "both sides), and it acts from the first accepted merge onward on every "
                    "pair with a merged side -- including merges that are never accepted, "
                    "since those scores still drive the argmax and the stop test."),
                "merged_alpha_missing_entry_fill": (
                    "C2 (implicit, and unchanged). combine_pips scales each column by its "
                    "coverage fraction and fills missing entries with 1/M'; read together "
                    "those two operations give an unmeasured variant the ARITHMETIC MEAN "
                    "Bayes factor of the variants that cohort did measure, log(S_k/n_k). "
                    "Section 4 of the note proves this is the UNIQUE posterior-neutral rule. "
                    "It is in force for every PIP reported here."),
                "uncovered_variant_backfill": (
                    "finalize_region: a variant no SELECTED component covers takes the MAX "
                    "marginal PIP among the cohorts absent from the final component set, and "
                    "note = 'filling_from_{cohort}'. This is a PIP copy, not a combination; "
                    "n_backfilled below counts them so downstream figures need not credit "
                    "FastMap for them."),
                "lbf_fastmap_column": (
                    "Added 2026-08-20 (Ran). Every row carries the log-BF that goes with its "
                    "reported PIP: a covered variant gets the lbf of the component named in "
                    "`note` (the same component its prob_fastmap is dominated by), and a "
                    "backfilled variant gets it from the same cohort its copied marginal PIP "
                    "came from. Which of that cohort's L components supplies it is the "
                    f"`lbf_backfill` convention: {setting.lbf_backfill}"
                    + (" -- max_i lbf over that cohort's components, the strongest "
                       "single-effect evidence it has for the variant."
                       if setting.lbf_backfill == "max_component" else
                       " -- the lbf of the component with the largest alpha at that variant.")
                    + " A log-sum-exp across components is deliberately NOT used: the L single "
                      "effects are fitted jointly on one dataset, so adding their Bayes factors "
                      "would count the same samples L times."),
                "pip_source": ("SuSiE's own `pip` vector (not alpha1) for SuSiE sources; the "
                               "ABF posterior `prob` for the single-causal ABF source."),
                "cross_cohort_merge": "outer join across cohorts (union of measured variants)",
                "variant_join_key": (
                    "real_data only: lifted GRCh38 composite ids chr{c}_{pos}_{ref}_{alt} -- "
                    "FinnGen's own `variant` column; UKBB translated through the precomputed "
                    "BCFtools/+liftover map (build_ukb_liftover_maps.py, pyliftover-gated). "
                    "rsid keying RETIRED 2026-08-24: it crashed on FinnGen's unnamed variants "
                    "(4.5%), matched 11.7% of indels, and mispaired multi-allelics across "
                    "dbSNP vintages (UKBB ~b142 vs FinnGen b155). Unlifted UKBB variants "
                    "(~0.001%) keep a b37_unlifted:: id, flow through as cohort-specific, and "
                    "are counted in input_accounting.ukb_variants_unlifted. Sims join on the "
                    "shared simulation variant ids, unchanged."),
            },
            "inputs": inputs_desc,
            "input_accounting": accounting,
            "stats": stats,
            "runtime_tsv": os.path.join(out_dir, f"{out_prefix}.runtime.tsv"),
            "elapsed_seconds_including_staging": round(time.time() - t0, 1),
            "code_sha256": code_hashes,
            "environment": {
                "python": sys.version.split()[0], "numpy": np.__version__,
                "pandas": pd.__version__, "platform": platform.platform(),
                "PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED"),
            },
        }
        with open(done_marker, "w") as fh:
            json.dump(manifest, fh, indent=2, default=str)
        print(f"[{label}] {stats['n_regions']} regions, {stats['n_rows']:,} rows, "
              f"{stats['total_cpu_seconds']:.0f} cpu-s, "
              f"{stats['frac_backfilled']:.2%} backfilled -> {stats['output_feather']}")

        if a.upload_to:
            _dest = a.upload_to.rstrip("/") + "/"
            subprocess.run(["gsutil", "-m", "-q", "cp",
                            stats["output_feather"],
                            os.path.join(out_dir, f"{out_prefix}.runtime.tsv"),
                            done_marker, _dest], check=True)
            print(f"[{label}] uploaded to {_dest}")

        if not a.keep_inputs and a.source != "real_data":
            subprocess.run(["rm", "-rf", staged["dir"]], check=False)

    if failures:
        print(f"\n{len(failures)} phenotype(s) FAILED:", file=sys.stderr)
        for label, err in failures:
            print(f"  {label}: {err}", file=sys.stderr)
        sys.exit(1)
    print(f"\nAll {len(jobs)} job(s) complete. Outputs in {out_dir}")


if __name__ == "__main__":
    main()
