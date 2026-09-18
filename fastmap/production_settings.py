"""The single source of truth for FastMap's recommended operating point.

Every analysis using FastMap should be run at ONE operating point, and this module
is where that point is defined.

For users of the software, the decision knobs -- the PP.H4 similarity `threshold`
(recommended 0.1), the pair-admissibility `overlap_min` (recommended 0.5) and the component
budget `L` (recommended 10) -- are defaults, not constants. Override any of them per run with
`FastMapSetting.customize()` or run_fastmap.py's --threshold / --overlap-min / --L flags.
A customized run is tagged as such in its output filenames (e.g. `107_L1`) and
its manifest.
"""
from __future__ import annotations

import csv
import os
from typing import NamedTuple

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: Paths to the grid sweep files used to verify settings.
GRID_TSVS = (
    os.path.join(_ROOT, "results", "simulations", "part2_sweep_rerun", "grid.tsv"),
    os.path.join(_ROOT, "results", "simulations", "part2_sweep_rerun",
                 "grid_thr01_supplement.tsv"),
)

#: The user-facing knobs that can be overridden per run.
CUSTOMIZABLE = ("threshold", "overlap_min", "L")


class FastMapSetting(NamedTuple):
    """One fully specified FastMap operating point.

    The field names match `fastmap.combine_region`'s keyword arguments exactly, except
    `threshold` -> `similarity_threshold`.
    """

    setting_id: int
    algorithm: str          # -> combine_region(similarity=...)
    threshold: float        # -> combine_region(similarity_threshold=...)
    overlap_min: float      # coloc.bf_bf's pair-admissibility trim, see coloc_susie_matrix
    L: int
    coloc_priors: tuple     # (p1, p2, p12); ignored unless algorithm == "coloc_susie"
    trim_by_posterior: bool = True
    lbf_backfill: str = "max_component"   # which component supplies a BACKFILLED variant's lbf
    
    #: User overrides applied on top of the base setting. Filled by `customize()`, 
    #: recorded in manifests, and reflected in `file_tag`.
    overrides: tuple = ()

    def customize(self, threshold: float = None, overlap_min: float = None,
                  L: int = None) -> "FastMapSetting":
        """Return this setting with any of the user-facing knobs overridden.

        Passing a knob's current value is a no-op (no override is recorded).
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
        customized one."""
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
                f"lbf_backfill={self.lbf_backfill}" + custom)


PRODUCTION = FastMapSetting(
    setting_id=107,
    algorithm="coloc_susie",
    threshold=0.1,
    overlap_min=0.5,
    L=10,
    coloc_priors=(1e-4, 1e-4, 1e-5),
    lbf_backfill="max_component",
)

#: Alternative settings on or near the frontier, for overrides and sensitivity runs.
ALTERNATIVES = {
    27: PRODUCTION._replace(setting_id=27, threshold=0.091),
    19: PRODUCTION._replace(setting_id=19, threshold=0.091, overlap_min=0.0),
    35: PRODUCTION._replace(setting_id=35, threshold=0.3, overlap_min=0.0),
    51: PRODUCTION._replace(setting_id=51, threshold=0.5, overlap_min=0.0),
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
        f"run_fastmap.py's --threshold/--overlap-min/--L flags."
    )


def assert_matches_grid(setting: FastMapSetting, grid_tsvs=GRID_TSVS) -> None:
    """Fail loudly if `setting` has drifted from the sweep grid it claims to come from.

    If no grid file is present, the check is skipped with a warning rather than blocking the run.

    A CUSTOMIZED setting skips the overridden fields and the remaining fields are checked against 
    the base setting's grid row.
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
