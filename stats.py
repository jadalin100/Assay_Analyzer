# =============================================================================
# stats.py — Statistical analysis: ANOVA / Kruskal-Wallis + post-hoc tests.
# =============================================================================
#
# Three continuous metrics are tested independently:
#
#   ratio_auc              — area under the normalized R/G ratio curve (ratio·min).
#                            Primary metric: captures both timing and magnitude of
#                            colour conversion.  Higher concentration → larger AUC.
#
#   time_to_ratio_threshold — first timepoint where normalized R/G ≥ RATIO_THRESHOLD.
#                            Higher concentration → earlier crossing.
#                            Groups that never cross are excluded from this test.
#
#   final_RG_ratio         — normalized R/G at the final timepoint.
#                            Simpler endpoint; useful for cross-trial comparison.
#
# Test selection (per metric, per group):
#   • Shapiro-Wilk normality check on each group's values.
#   • Levene's test for homogeneity of variance.
#   • One-way ANOVA if both pass; Kruskal-Wallis otherwise.
#   • Tukey HSD post-hoc (if statsmodels available) when ANOVA is significant.
#   • Dunn's test with Bonferroni correction (if scikit-posthocs available)
#     when Kruskal-Wallis is significant.
#
# Replication note:
#   ANOVA requires ≥ 3 independent observations per group.  Run 3 wells per
#   concentration in the same trial to satisfy this.  The pipeline automatically
#   pools multiple samples that share the same concentration label.
#
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

try:
    import pandas as pd
    import scikit_posthocs as sp
    HAS_POSTHOCS = True
except ImportError:
    HAS_POSTHOCS = False

# Metrics to test, with human-readable labels and expected direction of effect.
_METRICS = [
    {
        "key":       "total_detection_score",
        "label":     "Total detection score (SD + slope + Likert event counts)",
        "direction": "higher = stronger / earlier detection across all algorithms",
    },
    {
        "key":       "rg_auc",
        "label":     "AUC of normalized R/G ratio curve (ratio·min)",
        "direction": "higher = more/faster colour conversion",
    },
    {
        "key":       "rb_auc",
        "label":     "AUC of normalized R/B ratio curve (ratio·min)",
        "direction": "higher = more/faster colour conversion",
    },
    {
        "key":       "gb_auc",
        "label":     "AUC of normalized G/B ratio curve (ratio·min)",
        "direction": "lower = more conversion (G falls as resazurin converts)",
    },
    {
        "key":       "time_to_ratio_threshold",
        "label":     f"Time to R/G ≥ {getattr(config, 'RATIO_THRESHOLD', 1.10):.2f} (min)",
        "direction": "lower = faster conversion; groups that never cross are excluded",
    },
    {
        "key":       "final_RG_ratio",
        "label":     "Final normalized R/G ratio",
        "direction": "higher = more conversion by end of assay",
    },
]


def _group_values(results: list[dict], metric_key: str) -> dict[str, list[float]]:
    """
    Collect metric values grouped by concentration label.

    Multiple samples with the same concentration are pooled — this is the
    mechanism that handles replicate wells.  Values of None are excluded.
    """
    groups: dict[str, list[float]] = {}
    for r in results:
        conc = r["concentration_cfu_ml"]
        label = f"10^{int(round(np.log10(conc)))}" if conc > 0 else "0 (negative)"
        val = r.get(metric_key)
        if val is not None and not (isinstance(val, float) and np.isnan(val)):
            groups.setdefault(label, []).append(float(val))
    return groups


def _run_one_metric(results: list[dict], metric: dict) -> dict:
    """Run the full statistical pipeline for a single metric. Returns a result dict."""
    key   = metric["key"]
    label = metric["label"]
    groups = _group_values(results, key)

    out = {
        "metric":    key,
        "label":     label,
        "direction": metric["direction"],
        "groups":    {k: v for k, v in groups.items()},
    }

    # Minimum groups check
    valid_groups = {k: v for k, v in groups.items() if len(v) >= 1}
    if len(valid_groups) < 2:
        out["note"] = "Fewer than 2 groups with data — no test performed."
        return out

    # Replication warning
    low_n = [k for k, v in valid_groups.items() if len(v) < 3]
    if low_n:
        out["replication_warning"] = (
            f"Groups {low_n} have fewer than 3 replicates. "
            "Run 3+ wells per concentration for valid ANOVA. "
            "Test result is shown but has very low statistical power."
        )

    group_arrays  = [np.array(v) for v in valid_groups.values()]
    group_labels  = list(valid_groups.keys())

    # Normality (Shapiro-Wilk requires n ≥ 3; skip check for smaller groups)
    normal = True
    for arr in group_arrays:
        if len(arr) >= 3:
            _, p = scipy_stats.shapiro(arr)
            if p < config.ALPHA:
                normal = False
                break
        else:
            normal = False   # can't confirm normality with < 3 points

    # Levene's test for equal variances (requires ≥ 2 values per group)
    groups_for_levene = [a for a in group_arrays if len(a) >= 2]
    if len(groups_for_levene) >= 2:
        _, lev_p = scipy_stats.levene(*groups_for_levene)
    else:
        lev_p = 1.0

    # Main test
    # For time_to_ratio_threshold, groups that never cross are excluded (None),
    # so some groups may have 0 observations here — drop them before testing.
    testable = [a for a in group_arrays if len(a) >= 1]
    if len(testable) < 2:
        out["note"] = "Fewer than 2 groups with valid observations — no test performed."
        return out

    if normal and lev_p >= config.ALPHA and all(len(a) >= 3 for a in testable):
        stat, p_val = scipy_stats.f_oneway(*testable)
        test_name = "One-way ANOVA"
    else:
        stat, p_val = scipy_stats.kruskal(*testable)
        test_name = "Kruskal-Wallis"

    out["test"]        = test_name
    out["statistic"]   = float(stat)
    out["p_value"]     = float(p_val)
    out["significant"] = bool(p_val < config.ALPHA)

    sig_str = "SIGNIFICANT" if p_val < config.ALPHA else "not significant"
    print(f"    [{key}]  {test_name}: H/F={stat:.4f}  p={p_val:.4f}  {sig_str}")

    # Post-hoc tests when the overall test is significant
    if p_val < config.ALPHA:
        all_vals   = [v for vals in valid_groups.values() for v in vals]
        all_labels = [lbl for lbl, vals in valid_groups.items() for _ in vals]

        # Tukey HSD (ANOVA route, requires statsmodels)
        if HAS_STATSMODELS and test_name == "One-way ANOVA":
            tukey = pairwise_tukeyhsd(all_vals, all_labels, alpha=config.ALPHA)
            out["tukey_hsd"] = str(tukey.summary())
            print(f"      Tukey HSD:\n{tukey.summary()}")

        # Dunn's test with Bonferroni correction (Kruskal-Wallis route)
        elif HAS_POSTHOCS and test_name == "Kruskal-Wallis":
            dunn = sp.posthoc_dunn(
                pd.DataFrame({"val": all_vals, "group": all_labels}),
                val_col="val", group_col="group", p_adjust="bonferroni",
            )
            out["dunn_bonferroni"] = dunn.to_dict()
            print(f"      Dunn's (Bonferroni):\n{dunn.round(4)}")
        else:
            # Manual pairwise Mann-Whitney U with Bonferroni correction
            pairs = []
            n_pairs = len(group_labels) * (len(group_labels) - 1) // 2
            for i in range(len(group_labels)):
                for j in range(i + 1, len(group_labels)):
                    a, b = group_arrays[i], group_arrays[j]
                    if len(a) >= 1 and len(b) >= 1:
                        u, raw_p = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
                        adj_p = min(raw_p * n_pairs, 1.0)  # Bonferroni
                        pairs.append({
                            "group1": group_labels[i],
                            "group2": group_labels[j],
                            "U": float(u),
                            "p_raw": float(raw_p),
                            "p_bonferroni": float(adj_p),
                            "significant": bool(adj_p < config.ALPHA),
                        })
            out["pairwise_mannwhitney_bonferroni"] = pairs
            sig_pairs = [(p["group1"], p["group2"]) for p in pairs if p["significant"]]
            print(f"      Pairwise Mann-Whitney (Bonferroni): "
                  f"{len(sig_pairs)}/{len(pairs)} pairs significant")

    return out


def run_statistics(results: list[dict], output_path: str) -> dict:
    """
    Run statistical tests on all continuous kinetic metrics.
    Saves a unified JSON report and prints a summary to stdout.
    """
    print(f"  Testing {len(_METRICS)} metrics "
          f"(RATIO_THRESHOLD={getattr(config, 'RATIO_THRESHOLD', 1.10):.2f}):")

    stat_output = {"metrics": []}
    for metric in _METRICS:
        result = _run_one_metric(results, metric)
        stat_output["metrics"].append(result)

    _save(stat_output, output_path)
    return stat_output


def _save(data: dict, path: str):
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"  Statistics JSON saved: {path}")
