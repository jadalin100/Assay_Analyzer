# =============================================================================
# likert_algorithm.py — Colorimetric Likert detection algorithm.
#
# Matches the well's raw RGB colour at each timepoint to one of five reference
# levels derived from real trial photos (stored in config.LIKERT_REFERENCE_COLORS).
# Detection fires when the matched level exceeds LIKERT_POSITIVE_THRESHOLD for
# LIKERT_MIN_CONSECUTIVE consecutive frames.
#
# Why raw R, G, B?
#   The reference colours were sampled directly from JPEG well photos, so they
#   live on the same scale as rgb_data["R/G/B"] (white-balance corrected, 0-255).
#   The blue-purple → magenta-pink hue shift is large (~100 units in B, ~180 in G)
#   so small lighting offsets don't affect level assignment.
#
# !! RECALIBRATE config.LIKERT_REFERENCE_COLORS after the first locked-exposure
#    trial.  To recalibrate, run the pipeline and inspect the mean R, G, B values
#    in rgb_timeseries.csv for the negative control (Level 1) and 10^8 well at
#    t=120 (Level 5).  Interpolate Levels 2–4 between those anchors. !!
# =============================================================================

import math
import config


def _euclidean_rgb(r: float, g: float, b: float, ref: tuple) -> float:
    """Euclidean distance in RGB space between (r, g, b) and a reference colour."""
    dr = r - ref[0]
    dg = g - ref[1]
    db = b - ref[2]
    return math.sqrt(dr * dr + dg * dg + db * db)


def _closest_level(r: float, g: float, b: float,
                   ref_colors: list[tuple]) -> tuple[int, float]:
    """
    Return (level, distance) for the nearest reference colour.

    level is 1-indexed (1 = most negative, 5 = most positive by default).
    distance is the Euclidean RGB distance to the matched reference.
    """
    best_level = 1
    best_dist  = float("inf")
    for i, ref in enumerate(ref_colors):
        dist = _euclidean_rgb(r, g, b, ref)
        if dist < best_dist:
            best_dist  = dist
            best_level = i + 1
    return best_level, best_dist


def run_likert_algorithm(rgb_data: dict) -> dict:
    """
    Apply the Likert colour-matching algorithm to raw RGB time series.

    At each timepoint the mean well colour (rgb_data["R"], ["G"], ["B"]) is
    compared to LIKERT_REFERENCE_COLORS via Euclidean RGB distance.  The
    nearest reference is assigned as the Likert level (1–N, where N is the
    number of reference colours; typically 5).

    Detection logic mirrors the SD and slope algorithms:
      - A score >= LIKERT_POSITIVE_THRESHOLD at each frame is an "event".
      - Detection fires only after LIKERT_MIN_CONSECUTIVE consecutive events.
      - first_detection_min records the START of the first qualifying run.
      - Artifact timepoints (camera AE jumps) are treated as non-events and
        reset the consecutive counter, consistent with the other algorithms.

    Parameters
    ----------
    rgb_data : dict
        Must contain keys "times", "R", "G", "B".
        Optionally contains "artifact_timepoints" (set of times to suppress).

    Returns
    -------
    dict with keys:
        times            — list of timepoints (min)
        likert_scores    — list of int, one level (1-5) per timepoint
        match_distances  — list of float, RGB distance to matched reference
        event_flags      — list of bool, True where score >= threshold
        detected         — bool
        first_detection_min — float | None
    """
    times = rgb_data["times"]
    R     = rgb_data["R"]
    G     = rgb_data["G"]
    B     = rgb_data["B"]
    n     = len(times)

    ref_colors  = list(config.LIKERT_REFERENCE_COLORS)
    threshold   = int(getattr(config, "LIKERT_POSITIVE_THRESHOLD", 4))
    min_consec  = max(1, int(getattr(config, "LIKERT_MIN_CONSECUTIVE", 2)))
    artifact_times: set = set(rgb_data.get("artifact_timepoints", []))

    likert_scores:   list[int]   = []
    match_distances: list[float] = []
    event_flags:     list[bool]  = []

    detected           = False
    first_detection_min = None
    consecutive        = 0
    run_start_min      = None

    for i in range(n):
        t = times[i]

        # Compute nearest reference level regardless of artifact status
        # (we still want a score for the CSV/plot — just don't detect on it).
        level, dist = _closest_level(R[i], G[i], B[i], ref_colors)
        likert_scores.append(level)
        match_distances.append(round(dist, 2))

        # Artifact frames are suppressed: treat as a non-event and reset counter.
        if t in artifact_times:
            event_flags.append(False)
            consecutive   = 0
            run_start_min = None
            continue

        above = level >= threshold
        event_flags.append(above)

        if above:
            if consecutive == 0:
                run_start_min = t
            consecutive += 1
            if not detected and consecutive >= min_consec:
                detected            = True
                first_detection_min = run_start_min
        else:
            consecutive   = 0
            run_start_min = None

    return {
        "times":              times,
        "likert_scores":      likert_scores,
        "match_distances":    match_distances,
        "event_flags":        event_flags,
        "detected":           detected,
        "first_detection_min": first_detection_min,
    }
