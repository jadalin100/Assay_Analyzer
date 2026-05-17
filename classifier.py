# =============================================================================
# classifier.py — Per-sample classification using SD and slope algorithms.
# =============================================================================

import numpy as np
import config
from sd_algorithm import run_sd_algorithm
from slope_algorithm import run_slope_algorithm


def _compute_ratio_auc(times: list, vals: list) -> float | None:
    """
    Area under a ratio channel's curve using the trapezoidal rule (ratio·min).

    After reagent-reference normalization a non-converting group stays near 1.0,
    giving AUC ≈ total_duration_min.  A converting group rises above 1.0,
    giving a proportionally larger AUC.

    Returns None if fewer than 2 valid timepoints are available.
    """
    clean = [(t, v) for t, v in zip(times, vals) if v is not None]
    if len(clean) < 2:
        return None
    t_arr = np.array([p[0] for p in clean], dtype=float)
    v_arr = np.array([p[1] for p in clean], dtype=float)
    return float(np.trapz(v_arr, t_arr))


def _rate_of_change(times: list, rg_vals: list) -> float | None:
    """
    Linear slope of the normalized R/G ratio over the last SLOPE_CURRENT_POINTS
    timepoints (units: ratio / minute).

    Positive → ratio is still rising at trial end (active conversion).
    Negative → ratio peaked and is falling (conversion complete or reversing).
    Near zero → plateau reached.

    Uses SLOPE_CURRENT_POINTS for consistency with the slope algorithm's window.
    Returns None if fewer than 2 valid points are available in the window.
    """
    window = getattr(config, "SLOPE_CURRENT_POINTS", 4)
    clean = [(t, v) for t, v in zip(times, rg_vals) if v is not None]
    if len(clean) < 2:
        return None
    recent = clean[-window:]
    if len(recent) < 2:
        return None
    t_arr = np.array([p[0] for p in recent], dtype=float)
    v_arr = np.array([p[1] for p in recent], dtype=float)
    if len(np.unique(t_arr)) < 2:
        return None
    coeffs = np.polyfit(t_arr, v_arr, 1)
    return float(coeffs[0])


def _concordance(sd_first: float | None, slope_first: float | None) -> dict:
    """
    Assess agreement between SD and slope first_detection_min.

    Status:
      'concordant_negative'  — neither algorithm detected
      'concordant_positive'  — both detected, timing within 2 timepoints (10 min)
      'timing_discordant'    — both detected, but timing differs by > 10 min;
                               suggests one algorithm caught gradual variance
                               before the other caught the directional slope
      'algorithm_discordant' — one detected, the other did not;
                               the detecting algorithm may need threshold adjustment

    discordance_min is |sd_first - slope_first| when both detected, else None.
    """
    tolerance = 2 * getattr(config, "CAPTURE_INTERVAL_MIN", 5)   # 10 min at 5-min intervals

    if sd_first is None and slope_first is None:
        return {"status": "concordant_negative", "discordance_min": None}

    if sd_first is not None and slope_first is not None:
        diff = abs(sd_first - slope_first)
        status = "concordant_positive" if diff <= tolerance else "timing_discordant"
        return {"status": status, "discordance_min": float(diff)}

    return {"status": "algorithm_discordant", "discordance_min": None}


def _time_to_ratio_threshold(times: list, rg_vals: list) -> float | None:
    """
    Return the first timepoint (minutes) at which the normalized R/G ratio
    exceeds config.RATIO_THRESHOLD.

    Unlike first_detection_min (which depends on rolling statistical windows),
    this metric uses a fixed threshold on the normalized ratio and is therefore
    directly comparable across trials and concentrations.

    Returns None if the threshold is never crossed.
    """
    threshold = getattr(config, "RATIO_THRESHOLD", 1.10)
    for t, v in zip(times, rg_vals):
        if v is not None and v >= threshold:
            return float(t)
    return None


def _assign_load_tier(sd_first: float | None, slope_first: float | None) -> str:
    """
    Assign a bacterial load tier based on the earliest detection time across
    both algorithms. Returns a tier label from config.LOAD_TIERS, or 'negative'.
    """
    candidates = [t for t in [sd_first, slope_first] if t is not None]
    if not candidates:
        return "negative"
    earliest = min(candidates)
    for label, cutoff in sorted(config.LOAD_TIERS.items(), key=lambda x: x[1]):
        if earliest <= cutoff:
            return label
    return "negative"


def classify_sample(sample: dict) -> dict:
    """
    Run both detection algorithms on one sample and return combined results.

    Parameters
    ----------
    sample : dict with keys 'name', 'concentration', 'rgb_data'

    Returns
    -------
    dict with classification results for the sample
    """
    rgb_data = sample["rgb_data"]
    concentration = sample["concentration"]

    sd_result = run_sd_algorithm(rgb_data)
    slope_result = run_slope_algorithm(rgb_data)

    ground_truth = concentration >= config.UTI_POSITIVE_THRESHOLD_CFU

    final_rg = rgb_data["RG_ratio"][-1] if rgb_data.get("RG_ratio") else None
    final_rb = rgb_data["RB_ratio"][-1] if rgb_data.get("RB_ratio") else None
    final_gb = rgb_data["GB_ratio"][-1] if rgb_data.get("GB_ratio") else None

    times   = rgb_data["times"]
    rg_vals = rgb_data.get("RG_ratio", [])
    rb_vals = rgb_data.get("RB_ratio", [])
    gb_vals = rgb_data.get("GB_ratio", [])

    # ── Continuous kinetic metrics ────────────────────────────────────────────
    # Computed on whatever values are in rgb_data — reagent-normalized if
    # USE_REAGENT_REFERENCE=True, raw otherwise.

    # AUC for each ratio channel (ratio·min, trapezoidal rule).
    # Primary ANOVA variable once replicates exist.
    rg_auc = _compute_ratio_auc(times, rg_vals)
    rb_auc = _compute_ratio_auc(times, rb_vals)
    gb_auc = _compute_ratio_auc(times, gb_vals)

    # Time to fixed-threshold crossing on R/G (most sensitive primary channel).
    time_to_thr = _time_to_ratio_threshold(times, rg_vals)

    # Rate of change: slope of normalized R/G in the last SLOPE_CURRENT_POINTS
    # timepoints (ratio/min).  Distinguishes active-converting groups from
    # groups that already plateaued.
    roc = _rate_of_change(times, rg_vals)

    # Concordance between the two algorithms' first_detection_min.
    concordance = _concordance(
        sd_result["first_detection_min"],
        slope_result["first_detection_min"],
    )

    # ── Combined detection score ──────────────────────────────────────────────
    # SD event count (boolean sum over timepoints) + slope total weighted score
    # (continuous 0–6 per timepoint).  Single number for ANOVA across trials.
    sd_event_count        = int(sum(sd_result.get("detection_flags", [])))
    slope_total           = float(sum(slope_result.get("weighted_scores", [])))
    total_detection_score = slope_total + sd_event_count

    load_tier = _assign_load_tier(
        sd_result["first_detection_min"],
        slope_result["first_detection_min"],
    )

    return {
        "name":                     sample["name"],
        "concentration_cfu_ml":     concentration,
        "ground_truth_positive":    ground_truth,
        # ── SD algorithm ──────────────────────────────────────────────────────
        "sd_detected":              sd_result["detected"],
        "sd_first_detection_min":   sd_result["first_detection_min"],
        "sd_event_count":           sd_event_count,
        "sd_result":                sd_result,
        # ── Slope algorithm ───────────────────────────────────────────────────
        "slope_detected":           slope_result["detected"],
        "slope_first_detection_min": slope_result["first_detection_min"],
        "slope_total_score":        slope_total,
        "slope_result":             slope_result,
        # ── Combined score ────────────────────────────────────────────────────
        "total_detection_score":    total_detection_score,
        "load_tier":                load_tier,
        "concordance_status":       concordance["status"],
        "concordance_discordance_min": concordance["discordance_min"],
        # ── Continuous kinetic metrics ────────────────────────────────────────
        "rg_auc":                   rg_auc,
        "rb_auc":                   rb_auc,
        "gb_auc":                   gb_auc,
        "time_to_ratio_threshold":  time_to_thr,
        "rate_of_change_rg":        roc,
        # ── Final ratio values ────────────────────────────────────────────────
        "final_RG_ratio":           final_rg,
        "final_RB_ratio":           final_rb,
        "final_GB_ratio":           final_gb,
        "rgb_data":                 rgb_data,
    }


def classify_all(samples: list[dict]) -> list[dict]:
    """Classify all samples and return a list of result dicts."""
    results = []
    for sample in samples:
        result = classify_sample(sample)
        print(f"  {sample['name']:30s}  SD: {result['sd_detected']}  "
              f"Slope: {result['slope_detected']}  [{result['load_tier']}]")
        results.append(result)
    return results
