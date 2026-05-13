# =============================================================================
# reporting.py — Print and save per-group summary tables for both algorithms.
# =============================================================================

import os
import csv
import config


def _col(text: str, width: int, align: str = "<") -> str:
    return f"{str(text):{align}{width}}"


def _det(flag: bool) -> str:
    return "YES" if flag else "no"


# ---------------------------------------------------------------------------
# Slope table
# ---------------------------------------------------------------------------

def print_slope_table(results: list[dict]):
    """Print per-group slope-algorithm total score at each timepoint."""
    W_GROUP = 36
    W_TIME  = 7
    W_DET   = 10
    W_FIRST = 16

    # Gather all timepoints from first result
    first_times = results[0]["slope_result"]["times"]
    time_headers = "".join(_col(f"t={int(t)}", W_TIME) for t in first_times)

    header = (
        _col("Group", W_GROUP)
        + time_headers
        + _col("Detected", W_DET)
        + _col("First detect (min)", W_FIRST)
    )
    sep = "-" * len(header)

    print()
    print("=" * len(sep))
    print("  SLOPE ALGORITHM — Total channel score at each timepoint")
    print(f"  Detection threshold: total score > {config.SLOPE_WEIGHTED_SCORE_THRESHOLD}")
    print("=" * len(sep))
    print(header)
    print(sep)

    for r in results:
        sl = r["slope_result"]
        score_cells = "".join(
            _col(f"{s:.0f}", W_TIME) for s in sl["weighted_scores"]
        )
        first = r.get("slope_first_detection_min")
        first_str = f"{int(first)} min" if first is not None else "—"
        row = (
            _col(r["name"], W_GROUP)
            + score_cells
            + _col(_det(r["slope_detected"]), W_DET)
            + _col(first_str, W_FIRST)
        )
        print(row)

    print(sep)
    print()


def save_slope_csv(results: list[dict], out_dir: str):
    """Save slope total scores per timepoint as CSV."""
    first_times = results[0]["slope_result"]["times"]
    time_cols   = [f"score_t{int(t)}" for t in first_times]

    fieldnames = (
        ["group", "concentration_cfu_ml"]
        + time_cols
        + ["detected", "first_detection_min"]
    )

    path = os.path.join(out_dir, "slope_scores.csv")
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            sl  = r["slope_result"]
            row = {
                "group":                 r["name"],
                "concentration_cfu_ml":  r["concentration_cfu_ml"],
                **{col: round(s, 4) for col, s in zip(time_cols, sl["weighted_scores"])},
                "detected":              r["slope_detected"],
                "first_detection_min":   r.get("slope_first_detection_min", ""),
            }
            writer.writerow(row)
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# SD table
# ---------------------------------------------------------------------------

def print_sd_table(results: list[dict]):
    """Print per-group SD/mean (CV) at final timepoint. Ratio channels drive detection."""
    thr_all = config.SD_THRESHOLD_ALL
    thr_two = config.SD_THRESHOLD_TWO

    diag_channels = ["R(diag)", "G(diag)", "B(diag)"]
    det_channels  = ["R/G", "R/B", "G/B"]
    diag_keys     = ["ratio_R", "ratio_G", "ratio_B"]
    det_keys      = ["ratio_RG", "ratio_RB", "ratio_GB"]

    W_GROUP = 36
    W_VAL   = 10
    W_FLAG  = 12
    W_DET   = 8

    channel_header = (
        "".join(_col(ch, W_VAL) for ch in diag_channels)
        + "".join(_col(ch, W_VAL) for ch in det_channels)
    )
    header = (
        _col("Group", W_GROUP)
        + channel_header
        + _col("All>thr_all", W_FLAG)
        + _col("Two>thr_two", W_FLAG)
        + _col("Detected", W_DET)
    )
    sep = "-" * len(header)

    print("=" * len(sep))
    print("  SD ALGORITHM — CV (SD/mean) at final timepoint")
    print(f"  Detection channels: R/G, R/B, G/B only  |  R, G, B shown for diagnostics")
    print(f"  thr_all={thr_all}  (ALL ratio channels exceed)  thr_two={thr_two}  (ANY 2 ratio channels exceed)")
    print("=" * len(sep))
    print(header)
    print(sep)

    for r in results:
        sd  = r["sd_result"]
        diag_vals = [sd[k][-1] if sd[k] else 0.0 for k in diag_keys]
        det_vals  = [sd[k][-1] if sd[k] else 0.0 for k in det_keys]

        all_exceed = all(v > thr_all for v in det_vals)
        two_exceed = sum(v > thr_two for v in det_vals) >= 2

        val_cells = (
            "".join(_col(f"{v:.5f}", W_VAL) for v in diag_vals)
            + "".join(_col(f"{v:.5f}", W_VAL) for v in det_vals)
        )
        row = (
            _col(r["name"], W_GROUP)
            + val_cells
            + _col("YES" if all_exceed else "no", W_FLAG)
            + _col("YES" if two_exceed else "no", W_FLAG)
            + _col(_det(r["sd_detected"]), W_DET)
        )
        print(row)

    print(sep)
    print()


def save_sd_csv(results: list[dict], out_dir: str):
    """Save SD/mean values at final timepoint as CSV."""
    keys = ["ratio_R", "ratio_G", "ratio_B", "ratio_RG", "ratio_RB", "ratio_GB"]
    col_names = ["cv_R", "cv_G", "cv_B", "cv_RG", "cv_RB", "cv_GB"]

    fieldnames = (
        ["group", "concentration_cfu_ml"]
        + col_names
        + ["all_exceed_thr_all", "two_exceed_thr_two", "detected", "first_detection_min"]
    )

    thr_all = config.SD_THRESHOLD_ALL
    thr_two = config.SD_THRESHOLD_TWO

    path = os.path.join(out_dir, "sd_channel_values.csv")
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            sd   = r["sd_result"]
            vals = [sd[k][-1] if sd[k] else 0.0 for k in keys]
            row  = {
                "group":                 r["name"],
                "concentration_cfu_ml":  r["concentration_cfu_ml"],
                **{col: round(v, 6) for col, v in zip(col_names, vals)},
                "all_exceed_thr_all":    all(v > thr_all for v in vals),
                "two_exceed_thr_two":    sum(v > thr_two for v in vals) >= 2,
                "detected":              r["sd_detected"],
                "first_detection_min":   r.get("sd_first_detection_min", ""),
            }
            writer.writerow(row)
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# Detection events table
# ---------------------------------------------------------------------------

def print_detection_events(results: list[dict]):
    """
    Per-group table of every timepoint where a detection event fired.

    SD event:    criteria met at that timepoint (regardless of whether a prior
                 timepoint already triggered overall detection).
    Slope event: total channel score > SLOPE_WEIGHTED_SCORE_THRESHOLD at that
                 timepoint.
    """
    threshold = config.SLOPE_WEIGHTED_SCORE_THRESHOLD

    W_GROUP = 36
    W_ALGO  = 10
    W_TIMES = 40
    W_COUNT = 7

    header = (
        _col("Group", W_GROUP)
        + _col("Algorithm", W_ALGO)
        + _col("Event timepoints (min)", W_TIMES)
        + _col("Count", W_COUNT)
    )
    sep = "-" * len(header)

    print()
    print("=" * len(sep))
    print("  DETECTION EVENTS — timepoints where criteria were met")
    print("=" * len(sep))
    print(header)
    print(sep)

    for r in results:
        sd  = r["sd_result"]
        sl  = r["slope_result"]

        sd_events    = [t for t, flag in zip(sd["times"], sd["detection_flags"]) if flag]
        slope_events = [t for t, s in zip(sl["times"], sl["weighted_scores"]) if s > threshold]

        for algo, events in [("SD", sd_events), ("Slope", slope_events)]:
            times_str = ", ".join(f"{int(t)}" for t in events) if events else "—"
            print(
                _col(r["name"] if algo == "SD" else "", W_GROUP)
                + _col(algo, W_ALGO)
                + _col(times_str, W_TIMES)
                + _col(str(len(events)), W_COUNT)
            )
        print(_col("", W_GROUP) + "-" * (W_ALGO + W_TIMES + W_COUNT))

    print(sep)
    print()


def save_detection_events_csv(results: list[dict], out_dir: str):
    """Save detection events (timepoints where criteria met) as CSV."""
    threshold = config.SLOPE_WEIGHTED_SCORE_THRESHOLD
    fieldnames = [
        "group", "concentration_cfu_ml",
        "sd_event_count", "sd_event_timepoints_min",
        "slope_event_count", "slope_event_timepoints_min",
        "slope_total_score",
    ]
    path = os.path.join(out_dir, "detection_events.csv")
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            sd = r["sd_result"]
            sl = r["slope_result"]
            sd_events    = [t for t, flag in zip(sd["times"], sd["detection_flags"]) if flag]
            slope_events = [t for t, s in zip(sl["times"], sl["weighted_scores"]) if s > threshold]
            writer.writerow({
                "group":                        r["name"],
                "concentration_cfu_ml":         r["concentration_cfu_ml"],
                "sd_event_count":               len(sd_events),
                "sd_event_timepoints_min":      ";".join(str(int(t)) for t in sd_events),
                "slope_event_count":            len(slope_events),
                "slope_event_timepoints_min":   ";".join(str(int(t)) for t in slope_events),
                "slope_total_score":            round(sum(sl["weighted_scores"]), 2),
            })
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# Combined entry point
# ---------------------------------------------------------------------------

def print_and_save_tables(results: list[dict], out_dir: str):
    """Print both tables to the console and save them as CSVs."""
    print_slope_table(results)
    save_slope_csv(results, out_dir)
    print_sd_table(results)
    save_sd_csv(results, out_dir)
    print_detection_events(results)
    save_detection_events_csv(results, out_dir)
