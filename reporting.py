# =============================================================================
# reporting.py — Print and save per-group summary tables for both algorithms.
#
# CSV structure (all time-series files use long format: one row per group × timepoint)
# =============================================================================

import os
import csv
import config


def _col(text: str, width: int, align: str = "<") -> str:
    return f"{str(text):{align}{width}}"


def _det(flag: bool) -> str:
    return "YES" if flag else "no"


def _fmt_val(v, decimals: int = 4) -> str:
    if v is None or v == "":
        return "—"
    try:
        return f"{float(v):.{decimals}f}"
    except (TypeError, ValueError):
        return str(v)


# =============================================================================
# Slope algorithm
# =============================================================================

def print_slope_table(results: list[dict]):
    """Print per-group slope algorithm total score at each timepoint."""
    W_GROUP = 36
    W_TIME  = 7
    W_DET   = 10
    W_FIRST = 18

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
    print(f"  Threshold: score > {config.SLOPE_WEIGHTED_SCORE_THRESHOLD}  "
          f"for ≥ {config.SLOPE_MIN_CONSECUTIVE} consecutive frames")
    print("=" * len(sep))
    print(header)
    print(sep)

    for r in results:
        sl = r["slope_result"]
        score_cells = "".join(_col(f"{s:.1f}", W_TIME) for s in sl["weighted_scores"])
        first = r.get("slope_first_detection_min")
        first_str = f"{int(first)} min" if first is not None else "—"
        print(
            _col(r["name"], W_GROUP)
            + score_cells
            + _col(_det(r["slope_detected"]), W_DET)
            + _col(first_str, W_FIRST)
        )

    print(sep)
    print()


def save_slope_csv(results: list[dict], out_dir: str):
    """
    Save slope scores — long format: one row per group × timepoint.
    Summary section appended after a blank row.
    """
    path = os.path.join(out_dir, "slope_scores.csv")
    threshold = config.SLOPE_WEIGHTED_SCORE_THRESHOLD

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)

        # ── Per-timepoint section ──────────────────────────────────────────────
        writer.writerow(["# SLOPE SCORES — PER TIMEPOINT"])
        writer.writerow([
            "group", "concentration_cfu_ml", "time_min",
            "weighted_score", "above_threshold",
        ])
        for r in results:
            sl = r["slope_result"]
            for t, score in zip(sl["times"], sl["weighted_scores"]):
                writer.writerow([
                    r["name"],
                    r["concentration_cfu_ml"],
                    int(t),
                    round(score, 4),
                    score > threshold,
                ])

        # ── Per-group summary ──────────────────────────────────────────────────
        writer.writerow([])
        writer.writerow(["# SLOPE SCORES — PER-GROUP SUMMARY"])
        writer.writerow([
            "group", "concentration_cfu_ml",
            "slope_total_score", "slope_event_count",
            "detected", "first_detection_min",
        ])
        for r in results:
            sl = r["slope_result"]
            event_count = sum(1 for s in sl["weighted_scores"] if s > threshold)
            writer.writerow([
                r["name"],
                r["concentration_cfu_ml"],
                round(r.get("slope_total_score", 0.0), 4),
                event_count,
                r["slope_detected"],
                int(r["slope_first_detection_min"]) if r.get("slope_first_detection_min") is not None else "",
            ])

    print(f"  Saved: {path}")


# =============================================================================
# SD algorithm
# =============================================================================

def print_sd_table(results: list[dict]):
    """
    Print SD algorithm CV values at every timepoint for ratio channels.
    Threshold lines shown for reference.
    """
    thr_all = config.SD_THRESHOLD_ALL
    thr_two = config.SD_THRESHOLD_TWO

    W_GROUP = 30
    W_TIME  = 6
    W_CH    = 8
    W_FLAG  = 7

    # Build timepoints from first result
    first_times = results[0]["sd_result"]["times"]

    header_top = (
        _col("Group", W_GROUP)
        + "".join(_col(f"t={int(t)}", W_TIME + W_CH * 3, "^") for t in first_times)
        + _col("Det", W_FLAG)
        + _col("1st(min)", 9)
    )
    header_ch = _col("", W_GROUP) + "".join(
        _col("R/G", W_TIME) + _col("R/B", W_CH) + _col("G/B", W_CH)
        for _ in first_times
    )
    sep = "-" * len(header_top)

    print()
    print("=" * len(sep))
    print("  SD ALGORITHM — CV (SD/mean) at every timepoint — ratio channels only")
    print(f"  thr_all={thr_all} (ALL channels exceed)  "
          f"thr_two={thr_two} (ANY 2 channels exceed)  "
          f"min_consec={config.SD_MIN_CONSECUTIVE}")
    print("=" * len(sep))
    print(header_top)
    print(header_ch)
    print(sep)

    for r in results:
        sd = r["sd_result"]
        cells = ""
        for i in range(len(sd["times"])):
            rg = sd["ratio_RG"][i] if i < len(sd["ratio_RG"]) else 0.0
            rb = sd["ratio_RB"][i] if i < len(sd["ratio_RB"]) else 0.0
            gb = sd["ratio_GB"][i] if i < len(sd["ratio_GB"]) else 0.0
            cells += (
                _col(f"{rg:.3f}", W_TIME)
                + _col(f"{rb:.3f}", W_CH)
                + _col(f"{gb:.3f}", W_CH)
            )
        first = r.get("sd_first_detection_min")
        print(
            _col(r["name"], W_GROUP)
            + cells
            + _col(_det(r["sd_detected"]), W_FLAG)
            + _col(f"{int(first)}" if first is not None else "—", 9)
        )

    print(sep)
    print()


def save_sd_csv(results: list[dict], out_dir: str):
    """
    Save SD algorithm CV values — long format: one row per group × timepoint.
    Summary section appended after a blank row.
    """
    thr_all = config.SD_THRESHOLD_ALL
    thr_two = config.SD_THRESHOLD_TWO
    path = os.path.join(out_dir, "sd_channel_values.csv")

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)

        # ── Per-timepoint section ──────────────────────────────────────────────
        writer.writerow([
            f"# SD ALGORITHM — CV (SD/mean) PER TIMEPOINT",
            f"thr_all={thr_all}", f"thr_two={thr_two}",
        ])
        writer.writerow([
            "group", "concentration_cfu_ml", "time_min",
            "cv_R", "cv_G", "cv_B",          # diagnostic only — not used for detection
            "cv_RG", "cv_RB", "cv_GB",        # detection channels
            "sd_flag",                         # 1 if criteria met at this timepoint
            "meets_thr_all", "meets_thr_two",  # per-timepoint threshold checks
        ])
        for r in results:
            sd = r["sd_result"]
            times = sd["times"]
            for i, t in enumerate(times):
                def _g(key, i=i):
                    lst = sd.get(key, [])
                    return round(lst[i], 6) if i < len(lst) else ""

                cv_rg = sd["ratio_RG"][i] if i < len(sd["ratio_RG"]) else 0.0
                cv_rb = sd["ratio_RB"][i] if i < len(sd["ratio_RB"]) else 0.0
                cv_gb = sd["ratio_GB"][i] if i < len(sd["ratio_GB"]) else 0.0
                writer.writerow([
                    r["name"],
                    r["concentration_cfu_ml"],
                    int(t),
                    _g("ratio_R"), _g("ratio_G"), _g("ratio_B"),
                    round(cv_rg, 6), round(cv_rb, 6), round(cv_gb, 6),
                    int(bool(sd["detection_flags"][i])) if i < len(sd["detection_flags"]) else "",
                    all(v > thr_all for v in [cv_rg, cv_rb, cv_gb]),
                    sum(v > thr_two for v in [cv_rg, cv_rb, cv_gb]) >= 2,
                ])

        # ── Per-group summary ──────────────────────────────────────────────────
        writer.writerow([])
        writer.writerow(["# SD ALGORITHM — PER-GROUP SUMMARY (final timepoint CVs)"])
        writer.writerow([
            "group", "concentration_cfu_ml",
            "final_cv_RG", "final_cv_RB", "final_cv_GB",
            "sd_event_count", "detected", "first_detection_min",
        ])
        for r in results:
            sd = r["sd_result"]
            def _last(key):
                lst = sd.get(key, [])
                return round(lst[-1], 6) if lst else ""
            writer.writerow([
                r["name"],
                r["concentration_cfu_ml"],
                _last("ratio_RG"), _last("ratio_RB"), _last("ratio_GB"),
                r.get("sd_event_count", ""),
                r["sd_detected"],
                int(r["sd_first_detection_min"]) if r.get("sd_first_detection_min") is not None else "",
            ])

    print(f"  Saved: {path}")


# =============================================================================
# Likert algorithm
# =============================================================================

def print_likert_table(results: list[dict]):
    """Print per-group Likert score (1–5) at each timepoint."""
    W_GROUP = 36
    W_TIME  = 7
    W_DET   = 10
    W_FIRST = 18

    first_times = results[0]["likert_result"]["times"]
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
    print("  LIKERT ALGORITHM — Closest colour reference level at each timepoint")
    print(f"  Threshold: level >= {config.LIKERT_POSITIVE_THRESHOLD}  "
          f"for >= {config.LIKERT_MIN_CONSECUTIVE} consecutive frames")
    print("=" * len(sep))
    print(header)
    print(sep)

    for r in results:
        lk = r["likert_result"]
        score_cells = "".join(_col(str(s), W_TIME) for s in lk["likert_scores"])
        first = r.get("likert_first_detection_min")
        first_str = f"{int(first)} min" if first is not None else "—"
        print(
            _col(r["name"], W_GROUP)
            + score_cells
            + _col(_det(r["likert_detected"]), W_DET)
            + _col(first_str, W_FIRST)
        )

    print(sep)
    print()


def save_likert_csv(results: list[dict], out_dir: str):
    """
    Save Likert scores — long format: one row per group × timepoint.
    Summary section appended after a blank row.
    """
    path = os.path.join(out_dir, "likert_scores.csv")
    threshold = config.LIKERT_POSITIVE_THRESHOLD

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)

        # ── Reference colour table ─────────────────────────────────────────────
        writer.writerow(["# LIKERT REFERENCE COLORS"])
        writer.writerow(["level", "label", "R", "G", "B"])
        labels = [
            "No change (Negative)",
            "Slight change (Likely negative)",
            "Moderate change (Ambiguous)",
            "Clear change (Likely positive)",
            "Definite change (Positive)",
        ]
        for i, (color, label) in enumerate(zip(config.LIKERT_REFERENCE_COLORS, labels)):
            writer.writerow([i + 1, label, color[0], color[1], color[2]])
        writer.writerow([])

        # ── Per-timepoint section ──────────────────────────────────────────────
        writer.writerow(["# LIKERT SCORES — PER TIMEPOINT"])
        writer.writerow([
            "group", "concentration_cfu_ml", "time_min",
            "likert_score", "match_distance_rgb",
            "above_threshold",
        ])
        for r in results:
            lk = r["likert_result"]
            for t, score, dist, flag in zip(
                lk["times"], lk["likert_scores"],
                lk["match_distances"], lk["event_flags"],
            ):
                writer.writerow([
                    r["name"],
                    r["concentration_cfu_ml"],
                    int(t),
                    score,
                    round(dist, 2),
                    flag,
                ])

        # ── Per-group summary ──────────────────────────────────────────────────
        writer.writerow([])
        writer.writerow(["# LIKERT SCORES — PER-GROUP SUMMARY"])
        writer.writerow([
            "group", "concentration_cfu_ml",
            "likert_event_count", "final_likert_score",
            "detected", "first_detection_min",
        ])
        for r in results:
            lk = r["likert_result"]
            writer.writerow([
                r["name"],
                r["concentration_cfu_ml"],
                r.get("likert_event_count", 0),
                lk["likert_scores"][-1] if lk["likert_scores"] else "",
                r["likert_detected"],
                int(r["likert_first_detection_min"])
                    if r.get("likert_first_detection_min") is not None else "",
            ])

    print(f"  Saved: {path}")


# =============================================================================
# Detection events
# =============================================================================

def print_detection_events(results: list[dict]):
    """Per-group table of every timepoint where a detection event fired."""
    threshold = config.SLOPE_WEIGHTED_SCORE_THRESHOLD
    W_GROUP = 36
    W_ALGO  = 10
    W_TIMES = 44
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
        sd = r["sd_result"]
        sl = r["slope_result"]
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
    """
    Save detection events — one row per group with event timepoints listed,
    plus total_detection_score and concordance status for quick comparison.
    """
    threshold = config.SLOPE_WEIGHTED_SCORE_THRESHOLD
    path = os.path.join(out_dir, "detection_events.csv")

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["# DETECTION EVENTS — timepoints where each algorithm's criteria were met"])
        writer.writerow([
            "group", "concentration_cfu_ml",
            "sd_event_count", "sd_event_timepoints_min",
            "sd_detected", "sd_first_detection_min",
            "slope_event_count", "slope_event_timepoints_min",
            "slope_detected", "slope_first_detection_min",
            "slope_total_score",
            "likert_event_count", "likert_event_timepoints_min",
            "likert_detected", "likert_first_detection_min",
            "total_detection_score",
            "concordance_status", "concordance_discordance_min",
        ])
        for r in results:
            sd = r["sd_result"]
            sl = r["slope_result"]
            lk = r["likert_result"]
            sd_events     = [t for t, flag in zip(sd["times"], sd["detection_flags"]) if flag]
            slope_events  = [t for t, s in zip(sl["times"], sl["weighted_scores"]) if s > threshold]
            likert_threshold = getattr(config, "LIKERT_POSITIVE_THRESHOLD", 4)
            likert_events = [t for t, flag in zip(lk["times"], lk["event_flags"]) if flag]
            writer.writerow([
                r["name"],
                r["concentration_cfu_ml"],
                len(sd_events),
                ";".join(str(int(t)) for t in sd_events),
                r["sd_detected"],
                int(r["sd_first_detection_min"]) if r.get("sd_first_detection_min") is not None else "",
                len(slope_events),
                ";".join(str(int(t)) for t in slope_events),
                r["slope_detected"],
                int(r["slope_first_detection_min"]) if r.get("slope_first_detection_min") is not None else "",
                round(r.get("slope_total_score", 0.0), 4),
                len(likert_events),
                ";".join(str(int(t)) for t in likert_events),
                r["likert_detected"],
                int(r["likert_first_detection_min"]) if r.get("likert_first_detection_min") is not None else "",
                round(r.get("total_detection_score", 0.0), 4),
                r.get("concordance_status", ""),
                r.get("concordance_discordance_min", ""),
            ])

    print(f"  Saved: {path}")


# =============================================================================
# Total detection score sheet
# =============================================================================

def print_total_detection_table(results: list[dict]):
    """Print per-group summary of total_detection_score and concordance."""
    W_GROUP = 36
    W_VAL   = 18
    W_TOTAL = 14
    W_CONC  = 24

    header = (
        _col("Group", W_GROUP)
        + _col("SD event count", W_VAL)
        + _col("Slope total score", W_VAL)
        + _col("Likert events", W_VAL)
        + _col("TOTAL score", W_TOTAL)
        + _col("Concordance", W_CONC)
        + _col("Load tier", 24)
    )
    sep = "-" * len(header)

    print()
    print("=" * len(sep))
    print("  TOTAL DETECTION SCORE  =  SD event count  +  slope total weighted score  +  Likert event count")
    print("=" * len(sep))
    print(header)
    print(sep)

    for r in results:
        conc_status = r.get("concordance_status", "—")
        disc = r.get("concordance_discordance_min")
        conc_str = conc_status if disc is None else f"{conc_status} ({int(disc)} min)"
        print(
            _col(r["name"], W_GROUP)
            + _col(str(r.get("sd_event_count", 0)), W_VAL)
            + _col(_fmt_val(r.get("slope_total_score"), 2), W_VAL)
            + _col(str(r.get("likert_event_count", 0)), W_VAL)
            + _col(_fmt_val(r.get("total_detection_score"), 2), W_TOTAL)
            + _col(conc_str, W_CONC)
            + _col(r.get("load_tier", "—"), 24)
        )

    print(sep)
    print()


def save_total_detection_csv(results: list[dict], out_dir: str):
    """
    Save total detection scores — two sections in one file:

    Section 1 — per-timepoint breakdown: one row per group × timepoint showing
                 SD flag, slope score, timepoint contribution, and cumulative total.
    Section 2 — per-group summary: totals, kinetic metrics, and concordance.
    """
    path = os.path.join(out_dir, "total_detection_score.csv")
    threshold = config.SLOPE_WEIGHTED_SCORE_THRESHOLD

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)

        # ── Section 1: Per-timepoint breakdown ────────────────────────────────
        writer.writerow(["# SECTION 1 — PER-TIMEPOINT BREAKDOWN"])
        writer.writerow([
            "# sd_flag=1 when SD criteria met at that timepoint",
            "# slope_score=weighted slope score at that timepoint (0–6)",
            "# contribution=sd_flag + slope_score",
            "# cumulative=running total from t=0",
        ])
        writer.writerow([
            "group", "concentration_cfu_ml", "time_min",
            "sd_flag", "slope_score", "contribution", "cumulative_total",
        ])

        for r in results:
            sd = r["sd_result"]
            sl = r["slope_result"]
            sd_flags     = list(sd.get("detection_flags", []))
            slope_scores = list(sl.get("weighted_scores", []))
            times        = sd["times"]
            n = len(times)
            sd_flags     += [0]   * (n - len(sd_flags))
            slope_scores += [0.0] * (n - len(slope_scores))

            cumulative = 0.0
            for i, t in enumerate(times):
                sd_val   = int(bool(sd_flags[i]))
                sl_val   = float(slope_scores[i])
                contrib  = sd_val + sl_val
                cumulative += contrib
                writer.writerow([
                    r["name"],
                    r["concentration_cfu_ml"],
                    int(t),
                    sd_val,
                    round(sl_val, 4),
                    round(contrib, 4),
                    round(cumulative, 4),
                ])

        # ── Section 2: Per-group summary ──────────────────────────────────────
        writer.writerow([])
        writer.writerow(["# SECTION 2 — PER-GROUP SUMMARY"])
        writer.writerow([
            "group", "concentration_cfu_ml",
            "sd_event_count", "slope_total_score", "likert_event_count",
            "total_detection_score",
            "sd_detected", "slope_detected", "likert_detected",
            "concordance_status", "concordance_discordance_min",
            "load_tier",
            "rg_auc", "rb_auc", "gb_auc",
            "time_to_ratio_threshold_min",
            "rate_of_change_rg",
        ])
        for r in results:
            writer.writerow([
                r["name"],
                r["concentration_cfu_ml"],
                r.get("sd_event_count", ""),
                round(r.get("slope_total_score", 0.0), 4) if r.get("slope_total_score") is not None else "",
                r.get("likert_event_count", ""),
                round(r.get("total_detection_score", 0.0), 4) if r.get("total_detection_score") is not None else "",
                r.get("sd_detected", ""),
                r.get("slope_detected", ""),
                r.get("likert_detected", ""),
                r.get("concordance_status", ""),
                r.get("concordance_discordance_min", ""),
                r.get("load_tier", ""),
                round(r["rg_auc"], 2) if r.get("rg_auc") is not None else "",
                round(r["rb_auc"], 2) if r.get("rb_auc") is not None else "",
                round(r["gb_auc"], 2) if r.get("gb_auc") is not None else "",
                int(r["time_to_ratio_threshold"]) if r.get("time_to_ratio_threshold") is not None else "never",
                round(r["rate_of_change_rg"], 6) if r.get("rate_of_change_rg") is not None else "",
            ])

    print(f"  Saved: {path}")


# =============================================================================
# Combined entry point
# =============================================================================

def print_and_save_tables(results: list[dict], out_dir: str):
    """Print all tables to the console and save them as CSVs."""
    print_slope_table(results)
    save_slope_csv(results, out_dir)
    print_sd_table(results)
    save_sd_csv(results, out_dir)
    print_likert_table(results)
    save_likert_csv(results, out_dir)
    print_detection_events(results)
    save_detection_events_csv(results, out_dir)
    print_total_detection_table(results)
    save_total_detection_csv(results, out_dir)
