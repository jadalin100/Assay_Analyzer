# =============================================================================
# stats.py — Statistical analysis: ANOVA / Kruskal-Wallis + post-hoc tests.
# =============================================================================

import json
import numpy as np
from scipy import stats as scipy_stats
import config

try:
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False


def _group_values(results: list[dict], metric_key: str) -> dict[str, list[float]]:
    """Group metric values by concentration label."""
    groups: dict[str, list[float]] = {}
    for r in results:
        conc = r["concentration_cfu_ml"]
        label = f"{conc:.0e}" if conc > 0 else "0"
        val = r.get(metric_key)
        if val is not None and not (isinstance(val, float) and np.isnan(val)):
            groups.setdefault(label, []).append(float(val))
    return groups


def run_statistics(results: list[dict], output_path: str) -> dict:
    """
    Run one-way ANOVA (or Kruskal-Wallis fallback) on final RG ratio across groups.
    Adds Tukey HSD post-hoc if statsmodels is available.
    Saves results to a JSON file.
    """
    groups = _group_values(results, "final_RG_ratio")
    group_arrays = [np.array(v) for v in groups.values()]
    group_labels = list(groups.keys())

    stat_output = {"groups": {k: v for k, v in groups.items()}, "metric": "final_RG_ratio"}

    if len(group_arrays) < 2:
        stat_output["note"] = "Fewer than 2 groups — no statistical test performed."
        _save(stat_output, output_path)
        return stat_output

    # Normality test
    normal = True
    for arr in group_arrays:
        if len(arr) >= 3:
            _, p = scipy_stats.shapiro(arr)
            if p < config.ALPHA:
                normal = False
                break
        elif len(arr) < 3:
            normal = False

    # Levene test for equal variances
    if len(group_arrays) >= 2 and all(len(a) >= 2 for a in group_arrays):
        _, lev_p = scipy_stats.levene(*group_arrays)
    else:
        lev_p = 1.0

    if normal and lev_p >= config.ALPHA:
        stat, p_val = scipy_stats.f_oneway(*group_arrays)
        test_name = "One-way ANOVA"
    else:
        stat, p_val = scipy_stats.kruskal(*group_arrays)
        test_name = "Kruskal-Wallis"

    stat_output["test"] = test_name
    stat_output["statistic"] = float(stat)
    stat_output["p_value"] = float(p_val)
    stat_output["significant"] = bool(p_val < config.ALPHA)
    print(f"  {test_name}: statistic={stat:.4f}  p={p_val:.4f}  "
          f"{'significant' if p_val < config.ALPHA else 'not significant'}")

    # Post-hoc Tukey HSD
    if HAS_STATSMODELS and p_val < config.ALPHA:
        all_vals, all_labels = [], []
        for label, vals in groups.items():
            all_vals.extend(vals)
            all_labels.extend([label] * len(vals))
        tukey = pairwise_tukeyhsd(all_vals, all_labels, alpha=config.ALPHA)
        stat_output["tukey_hsd"] = str(tukey.summary())
        print(f"  Tukey HSD:\n{tukey.summary()}")

    _save(stat_output, output_path)
    return stat_output


def _save(data: dict, path: str):
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"  Statistics JSON saved: {path}")
