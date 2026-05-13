# =============================================================================
# slope_algorithm.py — Slope detection algorithm (Ciaccheri et al. 2023)
# =============================================================================

import math
import numpy as np
import config


def _linear_fit(x: list, y: list) -> tuple[float, float]:
    """Return (slope, R_squared) for a linear fit of y on x."""
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    if len(x) < 2:
        return 0.0, 0.0
    coeffs = np.polyfit(x, y, 1)
    slope = coeffs[0]
    y_pred = np.polyval(coeffs, x)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return float(slope), float(r2)


def _threshold_angle(channel_mean: float, multi: bool = False) -> float:
    """
    Compute detection threshold angle in degrees: θ = a·exp(−mean/b) + c.
    Called with the mean value of the channel's current window.
    """
    if multi:
        a, b, c = config.SLOPE_A_MULTI, config.SLOPE_B_MULTI, config.SLOPE_C_MULTI
    else:
        a, b, c = config.SLOPE_A_SINGLE, config.SLOPE_B_SINGLE, config.SLOPE_C_SINGLE
    return float(a * np.exp(-channel_mean / b) + c)


def _cosine_between_slopes(slope1: float, slope2: float) -> float:
    """
    Cosine of the angle between two slope vectors (1, slope1) and (1, slope2).

    Each slope m defines a direction vector (1, m) in (time, value) space.
    Formula: cos(θ) = (1 + m1·m2) / (√(1+m1²) · √(1+m2²))

    Returns 1.0 (no angular change) if denominator is zero.
    """
    denom = math.sqrt(1.0 + slope1 ** 2) * math.sqrt(1.0 + slope2 ** 2)
    if denom == 0.0:
        return 1.0
    return (1.0 + slope1 * slope2) / denom


def _channel_score(
    baseline_slope: float,
    current_slope: float,
    current_r2: float,
    channel_mean: float,
) -> int:
    """
    Score a single channel: 0, 1, or 2.

    Detection requires BOTH conditions:
      1. The angle between baseline slope and current slope exceeds the
         threshold angle (cos(angle) < cos(threshold_angle_in_degrees)).
      2. R² of the current window fit > 0.67 → score 1; > 0.80 → score 2.
    """
    threshold_deg = _threshold_angle(channel_mean)
    cos_threshold = math.cos(math.radians(threshold_deg))
    cos_angle = _cosine_between_slopes(baseline_slope, current_slope)

    if cos_angle < cos_threshold:
        if current_r2 > config.SLOPE_R2_SCORE_2:
            return 2
        elif current_r2 > config.SLOPE_R2_SCORE_1:
            return 1
    return 0


def run_slope_algorithm(rgb_data: dict) -> dict:
    """
    Apply the slope algorithm to RGB and ratio time series.

    Baseline window: fixed, always the first SLOPE_BASELINE_POINTS timepoints.
    Current window:  sliding, last SLOPE_CURRENT_POINTS timepoints up to i.

    Scoring begins at SLOPE_RAPID_CHANGE_FALLBACK (i=3, t=30 min by default).
    All 6 channels scored: R, G, B, R/G, R/B, G/B.
    Raw channels are included because the slope algorithm compares rates of
    change between windows — lighting drift affects both windows equally and
    does not produce false positives. SD uses ratios-only because variance IS
    sensitive to absolute brightness noise.

    Each channel is scored 0/1/2 via cosine check + R².
    Total score = unweighted sum across all 6 channels (max 12 per timepoint).
    Detection when total score > SLOPE_WEIGHTED_SCORE_THRESHOLD.
    """
    times = rgb_data["times"]
    n = len(times)

    # All 6 channels scored: R, G, B + R/G, R/B, G/B.
    # Raw channels are included here (unlike SD) because the slope algorithm
    # compares rate-of-change between windows, not absolute values — global
    # brightness drift produces similar slopes in both baseline and current
    # windows and does not cause false positives.
    R = rgb_data["R"]
    G = rgb_data["G"]
    B = rgb_data["B"]
    channels: dict[str, list] = {"R": R, "G": G, "B": B}
    for key in ("RG_ratio", "RB_ratio", "GB_ratio"):
        vals = rgb_data.get(key, [])
        if vals and all(v is not None for v in vals):
            channels[key] = vals

    bp              = min(config.SLOPE_BASELINE_POINTS, n)
    cp              = min(config.SLOPE_CURRENT_POINTS, n)
    score_threshold = config.SLOPE_WEIGHTED_SCORE_THRESHOLD
    effective_min   = min(config.SLOPE_RAPID_CHANGE_FALLBACK, n - 1)

    weighted_scores: list[float] = []
    detected = False
    first_detection_min = None

    for i in range(n):
        if i < effective_min:
            weighted_scores.append(0.0)
            continue

        # Fixed baseline: always the first bp timepoints
        baseline_end = min(bp, i)
        b_times = times[:baseline_end]

        # Sliding current window: last cp points up to and including i
        cur_start = max(0, i + 1 - cp)
        c_times   = times[cur_start: i + 1]

        total_score = 0
        for ch_values in channels.values():
            b_vals = ch_values[:baseline_end]
            c_vals = ch_values[cur_start: i + 1]

            baseline_slope, _  = _linear_fit(b_times, b_vals)
            current_slope, r2  = _linear_fit(c_times, c_vals)
            ch_mean            = float(np.mean(c_vals)) if c_vals else 0.0

            total_score += _channel_score(baseline_slope, current_slope, r2, ch_mean)

        weighted_scores.append(float(total_score))

        if not detected and total_score > score_threshold:
            detected = True
            first_detection_min = times[i]

    return {
        "times": times,
        "weighted_scores": weighted_scores,
        "detected": detected,
        "first_detection_min": first_detection_min,
    }
