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

    Scoring begins at SLOPE_RAPID_CHANGE_FALLBACK (i=10, t=50 min by default).
    Ratio channels only: R/G, R/B, G/B.
    Raw R/G/B excluded — auto-exposure drift creates false slope changes in
    absolute brightness. Ratios cancel this drift by construction.
    SLOPE_A/B/C_SINGLE parameters are calibrated for ratio-scale (~1.0).

    Each channel is scored 0/1/2 via cosine check + R².
    Total score = unweighted sum across all 6 channels (max 12 per timepoint).
    Detection when total score > SLOPE_WEIGHTED_SCORE_THRESHOLD.
    """
    times = rgb_data["times"]
    n = len(times)

    # Ratio channels only: R/G, R/B, G/B.
    # Raw R, G, B excluded — auto-exposure drift produces a real slope change
    # in absolute brightness (fast initial drop then stabilisation) that the
    # cosine check detects as a false positive even without any color change.
    # Ratios cancel global brightness drift by construction.
    # NOTE: SLOPE_A/B/C_SINGLE are calibrated for ratio-scale (~1.0), not 0-255.
    channels: dict[str, list] = {}
    for key in ("RG_ratio", "RB_ratio", "GB_ratio"):
        vals = rgb_data.get(key, [])
        if vals and all(v is not None for v in vals):
            channels[key] = vals

    bp              = min(config.SLOPE_BASELINE_POINTS, n)
    cp              = min(config.SLOPE_CURRENT_POINTS, n)
    score_threshold = config.SLOPE_WEIGHTED_SCORE_THRESHOLD
    effective_min   = min(config.SLOPE_RAPID_CHANGE_FALLBACK, n - 1)

    min_consec: int = max(1, getattr(config, "SLOPE_MIN_CONSECUTIVE", 2))

    # Camera artifact handling.
    # When artifact timepoints exist (set by main.py from the reagent reference
    # single-frame jump detector), both the baseline and the scoring loop must
    # avoid those frames so the contaminated spike does not distort the angular
    # comparison that drives the slope algorithm.
    #
    # Baseline strategy:
    #   Normal (no artifact): use first SLOPE_BASELINE_POINTS timepoints as
    #   designed — a confirmed stable early-assay window.
    #   Artifact present: filter artifact indices OUT of the baseline so the
    #   baseline slope is computed only on clean, uncontaminated frames.
    #   If fewer than 3 clean baseline frames remain, extend the baseline window
    #   forward (past the artifact) until 3+ clean frames are collected.
    #
    # Detection loop: suppress scoring at any artifact timepoint and reset the
    # consecutive counter (same pattern as the SD algorithm window reset).
    artifact_times: set = set(rgb_data.get("artifact_timepoints", []))

    # Build the clean baseline index list, extending if necessary.
    def _clean_baseline_indices(times, artifact_times, bp):
        """
        Return a list of up to bp indices that are not in artifact_times,
        drawn from the earliest timepoints.  If the first bp timepoints don't
        yield enough clean frames, extend forward until bp clean indices exist.
        """
        clean = []
        for idx in range(len(times)):
            if times[idx] not in artifact_times:
                clean.append(idx)
            if len(clean) >= bp:
                break
        return clean

    baseline_indices = _clean_baseline_indices(times, artifact_times, bp)
    b_times_clean = [times[i] for i in baseline_indices]

    # Pre-compute baseline slopes once (fixed for all detection timepoints).
    baseline_slopes: dict[str, float] = {}
    for key, ch_values in channels.items():
        b_vals_clean = [ch_values[i] for i in baseline_indices]
        baseline_slopes[key], _ = _linear_fit(b_times_clean, b_vals_clean)

    weighted_scores: list[float] = []
    detected = False
    first_detection_min = None
    consecutive = 0          # frames in a row currently above threshold
    run_start_min = None     # time of the first frame in the current run

    for i in range(n):
        # Suppress scoring at artifact timepoints and reset consecutive counter.
        if times[i] in artifact_times:
            weighted_scores.append(0.0)
            consecutive = 0
            run_start_min = None
            continue

        if i < effective_min:
            weighted_scores.append(0.0)
            consecutive = 0
            run_start_min = None
            continue

        # Sliding current window: last cp points up to and including i,
        # excluding any artifact timepoints.
        cur_indices = [
            j for j in range(max(0, i + 1 - cp), i + 1)
            if times[j] not in artifact_times
        ]
        if len(cur_indices) < 2:
            weighted_scores.append(0.0)
            continue
        c_times = [times[j] for j in cur_indices]

        total_score = 0
        for key, ch_values in channels.items():
            c_vals = [ch_values[j] for j in cur_indices]

            baseline_slope = baseline_slopes.get(key, 0.0)
            current_slope, r2 = _linear_fit(c_times, c_vals)
            ch_mean = float(np.mean(c_vals)) if c_vals else 0.0

            total_score += _channel_score(baseline_slope, current_slope, r2, ch_mean)

        weighted_scores.append(float(total_score))

        # Consecutive-frame persistence check.
        # A single isolated spike above the threshold is treated as noise.
        # Detection only fires once the score has exceeded the threshold for
        # SLOPE_MIN_CONSECUTIVE frames in a row.  first_detection_min is
        # recorded as the start of the run, not the frame that completed it.
        if total_score > score_threshold:
            if consecutive == 0:
                run_start_min = times[i]   # mark the beginning of this run
            consecutive += 1
            if not detected and consecutive >= min_consec:
                detected = True
                first_detection_min = run_start_min
        else:
            consecutive = 0
            run_start_min = None

    return {
        "times": times,
        "weighted_scores": weighted_scores,
        "detected": detected,
        "first_detection_min": first_detection_min,
    }
