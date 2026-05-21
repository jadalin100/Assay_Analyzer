# =============================================================================
# metrics.py — Sensitivity, Specificity, PPV, NPV for each algorithm.
# =============================================================================


def _confusion(results: list[dict], detected_key: str) -> tuple[int, int, int, int]:
    """Return (TP, TN, FP, FN) for a given detection key."""
    TP = TN = FP = FN = 0
    for r in results:
        gt = r["ground_truth_positive"]
        det = r[detected_key]
        if gt and det:
            TP += 1
        elif not gt and not det:
            TN += 1
        elif not gt and det:
            FP += 1
        elif gt and not det:
            FN += 1
    return TP, TN, FP, FN


def _safe_div(num: float, denom: float) -> float:
    return num / denom if denom > 0 else float("nan")


def compute_metrics(results: list[dict]) -> dict:
    """
    Compute sensitivity, specificity, PPV, and NPV for SD and slope algorithms.

    Returns a dict with metrics for each algorithm.
    """
    metrics = {}
    for algo, key in [
        ("SD",     "sd_detected"),
        ("Slope",  "slope_detected"),
        ("Likert", "likert_detected"),
    ]:
        TP, TN, FP, FN = _confusion(results, key)
        sensitivity = _safe_div(TP, TP + FN) * 100
        specificity = _safe_div(TN, TN + FP) * 100
        ppv = _safe_div(TP, TP + FP) * 100
        npv = _safe_div(TN, TN + FN) * 100
        metrics[algo] = {
            "TP": TP, "TN": TN, "FP": FP, "FN": FN,
            "sensitivity_pct": sensitivity,
            "specificity_pct": specificity,
            "PPV_pct": ppv,
            "NPV_pct": npv,
        }
        print(f"  {algo:6s}  Sensitivity: {sensitivity:.1f}%  Specificity: {specificity:.1f}%  "
              f"PPV: {ppv:.1f}%  NPV: {npv:.1f}%")
    return metrics
