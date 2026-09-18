#!/usr/bin/env python3
"""Production FastMap runner -- ONE entry point, ONE operating point, three data sources.

This script regenerates the FastMap arm of results figures at the selected operating point.
The algorithm is imported from `fastmap/fastmap.py`, the parameters come from
`fastmap/production_settings.py`, and the per-source plumbing lives in `fastmap_sources.py`.

Usage
-----
  # Figure 2 -- single-causal sims, one phenotype
  python3 scripts/production/run_fastmap.py --source sim_1causal --pheno 1

  # Figure 3 -- multi-causal sims, a range
  python3 scripts/production/run_fastmap.py --source sim_multicausal --phenos 1-100

  # Figure 4 -- real data, all FinnGen/UKBB endpoint pairs
  python3 scripts/production/run_fastmap.py --source real_data --all-phenos

  # sensitivity: same code, a different vetted operating point
  python3 scripts/production/run_fastmap.py --source sim_multicausal --pheno 1 --setting 35

  # custom parameters: override the decision knobs (threshold, overlap_min, L)
  python3 scripts/production/run_fastmap.py --source sim_1causal --phenos 1-400 --L 1

Outputs, per phenotype
----------------------
  {prefix}.fastmap.snp.feather   one row per (variant, region): prob_fastmap, note,
                                 comp_i_val / comp_i_name, region
  {prefix}.runtime.tsv           per-region cpu / wall / peak-RSS / n_snps / n_components.
  {prefix}.manifest.json         setting, conventions, input paths, code hashes, versions.
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

import fastmap_sources as src                                        # noqa: E402
from fastmap.fastmap import combine_region                           # noqa: E402
from fastmap.production_settings import (                            # noqa: E402
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
    """Run FastMap over one phenotype's regions; write the feather and the runtime table."""
    kwargs = setting.as_kwargs()
    per_region, runtime_rows = [], []
    for region, alpha, lbf, pip in regions_iter:
        t_cpu = time.process_time()
        t_wall = time.perf_counter()
        rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        try:
            res = combine_region(region, alpha, pip, lbf_df=lbf, **kwargs)
        except Exception as exc:                                     # noqa: BLE001
            raise RuntimeError(f"region {region} failed: {exc}") from exc
        cpu = time.process_time() - t_cpu
        wall = time.perf_counter() - t_wall
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        n_comp = sum(1 for c in res.columns if c.endswith("_name"))
        per_region.append(res)
        runtime_rows.append({
            "region": region, "n_snps": len(alpha), "n_components_in": alpha.shape[1],
            "n_components_out": n_comp, "cpu_seconds": cpu, "wall_seconds": wall,
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
        "n_backfilled": int(notes.str.startswith("filling_from_").sum()),
        "frac_backfilled": float(notes.str.startswith("filling_from_").mean()),
        "prob_fastmap_na": int(result["prob_fastmap"].isna().sum()),
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
                    help=f"vetted operating point id (default {PRODUCTION.setting_id})")
    ap.add_argument("--threshold", type=float, default=None,
                    help=f"PP.H4 similarity threshold override (default {PRODUCTION.threshold})")
    ap.add_argument("--overlap-min", type=float, default=None,
                    help=f"pair-admissibility overlap trim override (default {PRODUCTION.overlap_min})")
    ap.add_argument("--L", type=int, default=None,
                    help=f"max components in combined output override (default {PRODUCTION.L})")
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
        except Exception as exc:                                       # noqa: BLE001
            print(f"[{label}] FAILED: {exc}", file=sys.stderr)
            failures.append((label, str(exc)))
            continue

        manifest = {
            "source": a.source,
            "label": label,
            "setting": setting._asdict(),
            "conventions": {
                "merge_rule": "log-BFs add: lbf_merged = sum_k lbf_k, equivalently alpha_merged is the elementwise product renormalised.",
                "merged_lbf_missing_entry_fill": "lbf := log(S_k/n_k), the log arithmetic-mean Bayes factor over the variants that cohort measured.",
                "merged_alpha_missing_entry_fill": "combine_pips scales each column by its coverage fraction and fills missing entries with 1/M'.",
                "uncovered_variant_backfill": "A variant no selected component covers takes the max marginal PIP among unused cohorts.",
                "lbf_fastmap_column": (
                    f"Every row carries the log-BF that goes with its reported PIP. "
                    f"lbf_backfill convention: {setting.lbf_backfill}"
                ),
                "pip_source": "SuSiE's own `pip` vector for SuSiE sources; the ABF posterior `prob` for ABF sources.",
                "cross_cohort_merge": "Outer join across cohorts (union of measured variants).",
                "variant_join_key": "real_data: lifted GRCh38 composite ids. Sims: shared simulation variant ids.",
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
