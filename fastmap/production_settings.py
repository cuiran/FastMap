"""The single source of truth for FastMap's recommended operating point.

Every analysis in the manuscript -- Figure 2 (single-causal simulations), Figure 3
(multi-causal SuSiEx comparison) and Figure 4 (FinnGen R12 + UKBB real data) -- must be run
at ONE operating point, and this module is where that point is written down. Nothing else in
the repository is allowed to hard-code these values; import `PRODUCTION` instead.

For USERS of the software the three decision knobs -- the PP.H4 similarity `threshold`
(recommended 0.1), the pair-admissibility `overlap_min` (recommended 0.5) and the component
budget `L` (recommended 10) -- are defaults, not constants: override any of them per run with
`FastMapSetting.customize()` or run_fastmap.py's --threshold / --overlap-min / --L flags.
A customized run is tagged as such in its output filenames (`file_tag`, e.g. `107_L1`) and
its manifest, so it can never be mistaken for a run at the vetted point.

Provenance of setting 107
-------------------------
`setting_id = 107` is a row of the thr-0.1 supplement to the part-2 grid-v2 sweep,
`results/simulations/part2_sweep_rerun/grid_thr01_supplement.tsv` (ids 96-111 extend the v2
grid's 0-95). `assert_matches_grid()` re-reads the v2 grid files and refuses to run if the
numbers below have drifted from them.

The parameter values are CONVENTION-FIRST, chosen by Ran (2026-08-21) and then verified
rather than tuned:

    PP.H4 >= 0.1     a round threshold (the v2 grid's AP peak is at 0.091; 0.1 is its
                     conventional rounding, avoiding the appearance of overfitting)
    overlap_min 0.5  coloc.bf_bf's own default for the pair-admissibility trim
    L = 10           the prevailing choice in FINNGEN and UKBB fine-mapping analyses

Because 0.1 is not a v2 grid point, a dedicated 16-setting supplement sweep (thresholds all
0.1, same 759 tuning regions, workflows 917f0f90/edd808c4, 2026-08-22) VERIFIED the rounding
instead of presuming it: setting 107 matches the v2 default cell 27 (PP.H4 >= 0.091, same
overlap_min and L) to four decimals on every decision statistic, with exactly 4 of 1,564
causal-variant PIPs differing (max |delta| 0.064, none crossing 0.5 or 0.9).

Measured performance, pooled over both tuning arms (1,564 in-region causal variants),
at PIP >= 0.9:

    n at or above 0.9      175
    causal among them      163      -> recall 163 (10.42%), FDR 6.9%
    mean PIP               0.98647
    proportion causal      0.93143
    calibration gap        +0.0550
    Brier skill score      0.1784
    standalone CPU         ~280 s over 759 regions (supplement run; the v2 twin 27
                           measured 375 s -- same order, run-to-run VM jitter applies)

Computed by `scripts/analysis/part2_joint_metrics.py`, whose formulas are regression-checked
against `results/simulations/part2_sweep_rerun/decision_joint_metrics.tsv`.

Why THIS setting and not another frontier point
-----------------------------------------------
The v2 frontier is tight -- gap +0.050..+0.065 against recall 147..165 -- so no candidate
dominates. R1*'s mechanical pick was AGAIN inert-objective + CPU-tie-break (id 93, strictly
dominated: recall 154, gap +0.0709) and was reported, not adopted, same as in the retired
sweep. Within the frontier the choice is convention, recorded above. The v2 alternatives on
or near the frontier are kept in `ALTERNATIVES` with their measured numbers so sensitivity
runs do not retype them.

HISTORY: retirement of setting 221 (2026-08-22)
-----------------------------------------------
Production was previously `221` = `coloc_susie, PP.H4 >= 0.091, top_n 100, cst 0.1, L 10`, a
row of the RETIRED 360-row grid (`results/simulations/part2_sweep/grid.tsv`) measured under
the C1 lbf fill (gap +0.0509, recall 161/1,564, 624 CPU-s, all C1 numbers). Four
sweep-invalidating code changes (C2-only lbf fill, component-dedup removal, the
overlap-trim port, the F1 priority score) forced the grid-v2 rerun that replaced it; the v2
grid also dropped the `top_n` and `component_similarity_threshold` dimensions (top_n is
unread on the coloc-susie path since the dedup removal; cst was removed from the code).
This record therefore no longer carries those fields. The retired numbering (221/253/230/
182) indexes only the retired grid and must not be used with the v2 files.
"""
from __future__ import annotations

import csv
import os
from typing import NamedTuple

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#: The v2 grid and its thr-0.1 supplement share one numbering (0-95 / 96-111); a setting id
#: is looked up in whichever file holds it.
GRID_TSVS = (
    os.path.join(_ROOT, "results", "simulations", "part2_sweep_rerun", "grid.tsv"),
    os.path.join(_ROOT, "results", "simulations", "part2_sweep_rerun",
                 "grid_thr01_supplement.tsv"),
)


#: The user-facing knobs. Each has a recommended default (the values in `PRODUCTION`) and can
#: be overridden per run via `FastMapSetting.customize()` / run_fastmap.py's --threshold /
#: --overlap-min / --L flags. Everything else in a setting is an algorithm convention, not a
#: tuning knob, and stays pinned.
CUSTOMIZABLE = ("threshold", "overlap_min", "L")


class FastMapSetting(NamedTuple):
    """One fully specified FastMap operating point.

    The field names match `fastmap.combine_region`'s keyword arguments exactly, except
    `threshold` -> `similarity_threshold`, so `as_kwargs()` can be splatted straight in.
    `top_n` is deliberately NOT a field: it is unread on the coloc-susie path (the Jaccard
    dedup that consumed it was removed 2026-08-21) and no weighted-jaccard operating point
    is vetted, so `combine_region`'s own default applies and cannot affect output.
    """

    setting_id: int
    algorithm: str          # -> combine_region(similarity=...)
    threshold: float        # -> combine_region(similarity_threshold=...)
    overlap_min: float      # coloc.bf_bf's pair-admissibility trim, see coloc_susie_matrix
    L: int
    coloc_priors: tuple     # (p1, p2, p12); ignored unless algorithm == "coloc_susie"
    trim_by_posterior: bool = True
    lbf_backfill: str = "max_component"   # which component supplies a BACKFILLED variant's lbf
    #: user overrides applied on top of the base setting, as ((name, value), ...); empty for
    #: a vetted grid point. Filled by `customize()`, recorded in manifests, and reflected in
    #: `file_tag` so a customized run can never masquerade as the vetted setting.
    overrides: tuple = ()

    def customize(self, threshold: float = None, overlap_min: float = None,
                  L: int = None) -> "FastMapSetting":
        """Return this setting with any of the user-facing knobs overridden.

        Passing a knob's current value is a no-op (no override is recorded), so scripted
        callers can always pass all three. Values are validated here so a typo fails at
        parse time, not deep inside `combine_region`.
        """
        if self.overrides:
            raise ValueError(f"{self.describe()} is already customized; customize the vetted "
                             f"base setting instead of stacking overrides")
        changes = []
        if threshold is not None and threshold != self.threshold:
            if not 0.0 < threshold <= 1.0:
                raise ValueError(f"threshold (PP.H4) must be in (0, 1], got {threshold}")
            changes.append(("threshold", float(threshold)))
        if overlap_min is not None and overlap_min != self.overlap_min:
            if not 0.0 <= overlap_min <= 1.0:
                raise ValueError(f"overlap_min must be in [0, 1], got {overlap_min}")
            changes.append(("overlap_min", float(overlap_min)))
        if L is not None and L != self.L:
            if not (isinstance(L, int) and L >= 1):
                raise ValueError(f"L must be an integer >= 1, got {L!r}")
            changes.append(("L", int(L)))
        if not changes:
            return self
        return self._replace(overrides=tuple(changes),
                             **{name: value for name, value in changes})

    @property
    def file_tag(self) -> str:
        """Token for output filenames/dirs: '107' for a vetted point, '107_L1' etc. for a
        customized one, so custom outputs never collide with (or impersonate) vetted ones."""
        tag = str(self.setting_id)
        for name, value in self.overrides:
            short = {"threshold": "thr", "overlap_min": "ov", "L": "L"}[name]
            tag += f"_{short}{value:g}"
        return tag

    def as_kwargs(self) -> dict:
        """Keyword arguments for `fastmap.fastmap()` / `fastmap.combine_region()`."""
        return {
            "similarity": self.algorithm,
            "similarity_threshold": self.threshold,
            "L": self.L,
            "coloc_priors": self.coloc_priors,
            "overlap_min": self.overlap_min,
            "trim_by_posterior": self.trim_by_posterior,
            "lbf_backfill": self.lbf_backfill,
        }

    def describe(self) -> str:
        custom = ""
        if self.overrides:
            custom = (" [CUSTOM: "
                      + ", ".join(f"{n}={v}" for n, v in self.overrides)
                      + f" overriding setting {self.setting_id}]")
        return (f"setting {self.file_tag}: {self.algorithm} "
                f"threshold={self.threshold} overlap_min={self.overlap_min} L={self.L} "
                f"coloc_priors={self.coloc_priors} "
                f"trim_by_posterior={self.trim_by_posterior} "
                f"lbf_fill=C2[only option since 2026-08-21] "
                f"lbf_backfill={self.lbf_backfill}" + custom)


PRODUCTION = FastMapSetting(
    setting_id=107,
    algorithm="coloc_susie",
    threshold=0.1,
    overlap_min=0.5,
    L=10,
    # coloc defaults for (p1, p2); p12 = 1e-5 is coloc.abf's default (coloc.bf_bf's own p12
    # default is 5e-6, a within-coloc discrepancy). Part 1 proved the pairwise RANKING is
    # exactly p12-invariant (PP.H4 is monotone in S12/D0) and the swept threshold absorbs the
    # scale, so the "combine.abf on per-signal SuSiE lbf" framing keeps 1e-5.
    coloc_priors=(1e-4, 1e-4, 1e-5),
    # lbf fill is NOT a field: C2 (a cohort's MISSING entry takes its mean Bayes factor,
    # log(S_k/n_k), the same fill `combine_pips` applies to the alpha side of the very same
    # merge) became the default on 2026-08-20 and the ONLY behaviour on 2026-08-21, when Ran
    # removed the C1 option (`lbf := 0`) from the code entirely. See `fastmap.merge_lbf`.
    # A backfilled variant's `lbf_fastmap` comes from the cohort its copied marginal PIP came
    # from (Ran, 2026-08-20). That cohort has L components and the marginal PIP aggregates all
    # of them, so which component's lbf to read is a convention: "max_component" takes the
    # strongest single-effect evidence that cohort has for the variant; the alternative,
    # "pip_argmax_component", disagrees on 55% of rows. See `fastmap.backfill_lbf`.
    lbf_backfill="max_component",
)

#: v2-grid alternatives on or near the (gap@0.9, recall@0.9) frontier, for `--setting`
#: overrides and sensitivity runs. These are NOT production. Numbers are pooled A+B from the
#: v2 rerun (gap@0.9 / recall@0.9 of 1,564 / FDR@0.9 / standalone CPU).
ALTERNATIVES = {
    #  27: the v2 grid's default cell, PP.H4 >= 0.091 -- production's un-rounded twin.
    #      +0.0550 / 163 / 6.9% / 375 s. Numerically ~= production (4 causal PIPs differ).
    27: PRODUCTION._replace(setting_id=27, threshold=0.091),
    #  19: same threshold, trim-vacuous anchor (overlap_min 0). +0.0509 / 161 / 6.4% / 491 s.
    19: PRODUCTION._replace(setting_id=19, threshold=0.091, overlap_min=0.0),
    #  35: recall plateau. +0.0596 / 165 / 7.3% / 444 s.
    35: PRODUCTION._replace(setting_id=35, threshold=0.3, overlap_min=0.0),
    #  51: recall plateau, best BSS of the 165-club. +0.0645 / 165 / 7.8% / 414 s.
    51: PRODUCTION._replace(setting_id=51, threshold=0.5, overlap_min=0.0),
    #  11: loosest threshold; under trim+F1 the old C1-era calibration cliff is GONE, so this
    #      is no longer a low-recall outlier. +0.0500 / 147 / 6.4% / 403 s.
    11: PRODUCTION._replace(setting_id=11, threshold=0.03),
}


def get_setting(setting_id: int = PRODUCTION.setting_id) -> FastMapSetting:
    if setting_id == PRODUCTION.setting_id:
        return PRODUCTION
    if setting_id in ALTERNATIVES:
        return ALTERNATIVES[setting_id]
    raise ValueError(
        f"setting {setting_id} is not one of the vetted operating points "
        f"{[PRODUCTION.setting_id] + sorted(ALTERNATIVES)}. To run at custom parameter "
        f"values, override the user-facing knobs on a vetted base instead -- "
        f"`get_setting().customize(threshold=..., overlap_min=..., L=...)`, or "
        f"run_fastmap.py's --threshold/--overlap-min/--L flags -- so the deviation is "
        f"recorded in the manifest and the file tag. New MANUSCRIPT operating points still "
        f"belong here with their provenance. (Retired-grid ids such as 221 index the "
        f"pre-rerun 360-row grid and have no meaning in the v2 numbering.)")


def assert_matches_grid(setting: FastMapSetting, grid_tsvs=GRID_TSVS) -> None:
    """Fail loudly if `setting` has drifted from the sweep grid it claims to come from.

    This is cheap insurance against the failure mode that produced the whole part-2 exercise:
    a production operating point whose parameter values nobody could trace back to a
    measurement. The id is looked up across the v2 grid files (main + thr-0.1 supplement,
    disjoint id ranges). If no grid file is present (e.g. a deployed container that carries
    only the algorithm), the check is skipped with a warning rather than blocking the run.

    A CUSTOMIZED setting (non-empty `overrides`) is deliberately off-grid in the overridden
    fields; those fields are skipped and the remaining fields are still checked against the
    base setting's grid row, so the un-customized part of the pedigree stays verified.
    """
    present = [p for p in grid_tsvs if os.path.exists(p)]
    if not present:
        print(f"WARNING: none of {grid_tsvs} found; cannot verify {setting.describe()} "
              f"against the sweep grid")
        return
    rows = {}
    for path in present:
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows[int(r["setting_id"])] = r
    if setting.setting_id not in rows:
        raise AssertionError(f"setting_id {setting.setting_id} absent from {present}")
    row = rows[setting.setting_id]
    overridden = {name for name, _ in setting.overrides}
    mismatches = []
    for name, got, want in [
        ("algorithm", setting.algorithm, row["algorithm"]),
        ("threshold", setting.threshold, float(row["threshold"])),
        ("overlap_min", setting.overlap_min, float(row["overlap_min"])),
        ("L", setting.L, int(row["L"])),
    ]:
        if name in overridden:
            continue
        if got != want:
            mismatches.append(f"{name}: module says {got!r}, grid says {want!r}")
    if mismatches:
        raise AssertionError(
            f"setting {setting.setting_id} does not match the sweep grid:\n  "
            + "\n  ".join(mismatches))
