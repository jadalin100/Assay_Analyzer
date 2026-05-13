# =============================================================================
# sd_algorithm.py — Standard Deviation detection algorithm (Malanoski et al. 2016)
# =============================================================================

import numpy as np
import config


def _sd_ratio(values: list[float]) -> float:
    """Compute unweighted SD / mean for a channel window."""
    arr = np.array(values, dtype=float)
    if len(arr) < 2:
        return 0.0
    mean = np.mean(arr)
    if mean == 0:
        return 0.0
    return float(np.std(arr, ddof=1) / mean)


def _detection_triggered(ratios: list[float], thr_all: float, thr_two: float) -> bool:
    """
    Return True if:
      - ALL values exceed thr_all, OR
      - ANY two values exceed thr_two
    """
    all_exceed = all(v > thr_all for v in ratios)
    two_exceed = sum(v > thr_two for v in ratios) >= 2
    return all_exceed or two_exceed


def run_sd_algorithm(rgb_data: dict) -> dict:
    """
    Apply the SD algorithm to RGB and ratio time series.

    For each timepoint, uses all data up to that point (rolling window capped
    at SD_WINDOW) to compute SD/mean for R, G, B, R/G, R/B, and G/B.

    Primary channels (R, G, B) and ratio channels (R/G, R/B, G/B) are evaluated
    independently. Detection triggers if either set meets the criteria:
      - ALL channels in the set > SD_THRESHOLD_ALL, OR
      - ANY two channels in the set > SD_THRESHOLD_TWO
    """
    times = rgb_data["times"]
    R = rgb_data["R"]
    G = rgb_data["G"]
    B = rgb_data["B"]
    RG = rgb_data.get("RG_ratio", [])
    RB = rgb_data.get("RB_ratio", [])
    GB = rgb_data.get("GB_ratio", [])
    n = len(times)

    window = config.SD_WINDOW
    thr_all = config.SD_THRESHOLD_ALL
    thr_two = config.SD_THRESHOLD_TWO

    ratio_R, ratio_G, ratio_B = [], [], []
    ratio_RG, ratio_RB, ratio_GB = [], [], []
    detection_flags = []   # True at timepoint i if criteria met (regardless of prior detection)
    detected = False
    first_detection_min = None

    for i in range(n):
        start = max(0, i + 1 - window)

        # Primary channels — unweighted SD/mean
        rr = _sd_ratio(R[start: i + 1])
        rg = _sd_ratio(G[start: i + 1])
        rb = _sd_ratio(B[start: i + 1])
        ratio_R.append(rr)
        ratio_G.append(rg)
        ratio_B.append(rb)

        # Ratio channels (filter None values present when ratio is disabled)
        rrg = _sd_ratio([v for v in RG[start: i + 1] if v is not None]) if RG else 0.0
        rrb = _sd_ratio([v for v in RB[start: i + 1] if v is not None]) if RB else 0.0
        rgb_ = _sd_ratio([v for v in GB[start: i + 1] if v is not None]) if GB else 0.0
        ratio_RG.append(rrg)
        ratio_RB.append(rrb)
        ratio_GB.append(rgb_)

        # Detection driven by ratio channels only (R/G, R/B, G/B).
        # Raw R, G, B CVs are still computed for visualization but do not
        # contribute to detection — they are too sensitive to lighting and
        # JPEG compression noise.
        points_in_window = (i + 1) - start
        if points_in_window >= 3:
            event = _detection_triggered([rrg, rrb, rgb_], thr_all, thr_two)
        else:
            event = False

        detection_flags.append(event)
        if not detected and event:
            detected = True
            first_detection_min = times[i]

    return {
        "times": times,
        "ratio_R": ratio_R,
        "ratio_G": ratio_G,
        "ratio_B": ratio_B,
        "ratio_RG": ratio_RG,
        "ratio_RB": ratio_RB,
        "ratio_GB": ratio_GB,
        "detection_flags": detection_flags,
        "detected": detected,
        "first_detection_min": first_detection_min,
    }
