# FastMap (Beta)

**FastMap (Beta)** is a beta implementation of a fast, federated multi-biobank fine-mapping method.  
It takes SuSiE outputs from separate cohorts (ideally fine-mapped with in-sample LD), and combining them to approximate joint Posterior Inclusion Probabilities (PIPs).
FastMap (Beta) works in single-causal-variant setting, as well as multi-causal-variant setting.

---

## Overview

The main entry point is the function:
```python
fastmap_multi_causal(
    results_dict,
    pips_dict,
    max_thresh=0.1,
    initial_rev_max_thresh=0.001,
    max_div=3,
    pip_thresh=0.5,
    jaccard_thresh=0,
    L=10,
    top_n=100
)
```

## Arguments

| Argument                 | Type    | Default    | Description                                                     |
| ------------------------ | ------- | ---------- | --------------------------------------------------------------- |
| `results_dict`           | `dict`  | *required* | Dictionary of SuSiE alpha values from all cohorts.      |
| `pips_dict`              | `dict`  | *required* | Dictionary of posterior inclusion probabilities (PIPs).         |
| `max_thresh`             | `float` | `0.1`      | Maximum threshold for calling signal overlap.                    |
| `initial_rev_max_thresh` | `float` | `0.001`    | Initial reverse threshold for refining overlapping signals. |
| `max_div`                | `int`   | `3`        | Maximum denominator to `max_thresh`.       |
| `pip_thresh`             | `float` | `0.5`      | Minimum PIP threshold for remaining signals.                    |
| `jaccard_thresh`         | `float` | `0.0`      | Jaccard similarity threshold for merging overlapping signals.   |
| `L`                      | `int`   | `10`       | Number of causal components to consider.                |
| `top_n`                  | `int`   | `100`      | Number of top variants to include for checking signal overlap.             |

## Current optimal parameters

| Parameter Name       | Value |
|----------------------|-------|
| `max_div`            | `1`   |
| `max_threshold`      | `0.0` |
| `rev_max_threshold`  | `0.001` |
| `top_n`              | `500` |
| `jaccard_threshold`  | `0.0` |
| `pip_threshold`      | `0.1` |
| `L`                  | `10`  |
| `susie_l`            | `10`  |

## Status
In development.
For questions, suggestions please contact rancui@broadinstitute.org
