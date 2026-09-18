# FastMap

FastMap combines per-cohort SuSiE fine-mapping results across cohorts into a single set of
posterior inclusion probabilities (PIPs) per variant, without joint modeling. It greedily
merges single-effect components across cohorts that colocalize (by `PP.H4` from
coloc-susie, or weighted Jaccard), then reports one final PIP per variant from the `L`
strongest resulting components by evidence.

## Repository layout

- `fastmap/` — the current, published algorithm: the core combination logic
  (`fastmap.py`), the pairwise coloc scoring it calls (`coloc.py`), and the vetted
  production operating point plus grid/sensitivity settings (`production_settings.py`).
- `scripts/production/` — drivers that read per-cohort SuSiE output and real-data summary
  statistics into the `region_df` / `pips_df` inputs `fastmap.fastmap.combine_region`
  expects, and run FastMap at a given setting (`run_fastmap.py`, `fastmap_sources.py`).
- `scripts/`, `wdl/` (root-level) — earlier beta scripts and WDL, kept for history; superseded
  by `fastmap/` and `scripts/production/` above.

## Installation

```
pip install -r requirements.txt
```

Python 3.10+.

## Usage

`fastmap.fastmap.combine_region` is the core entry point; see its docstring for the expected
`region_df` / `pips_df` schema. `scripts/production/run_fastmap.py` is an example driver:

```
python scripts/production/run_fastmap.py --help
```

`scripts/production/fastmap_sources.py` reads inputs from GCS paths under an environment
variable, `FASTMAP_DATA_ROOT`, pointing at your own bucket with the layout documented in
that file's comments; it does not include any private data.

## Status

Method and code accompany the FastMap manuscript. See the manuscript's Code Availability
section for the version used to generate published results.

## License

MIT (see `LICENSE`).

For questions, contact rancui@broadinstitute.org.
