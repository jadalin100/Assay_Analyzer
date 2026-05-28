"""
build_averaged_results.py
=========================
Aggregates all 4 trial result CSVs into a single averaged_results/ folder.

Outputs
-------
averaged_rgb_timeseries.csv       — mean ± SD of R, G, B, RG/RB/GB ratios per group × timepoint
averaged_sd_channel_values.csv    — mean ± SD of all CV channels + sd_flag rate per group × timepoint
averaged_slope_scores.csv         — mean ± SD of weighted_score + detection rate per group × timepoint
averaged_likert_scores.csv        — mean ± SD of likert_score + match_distance per group × timepoint
averaged_total_timeseries.csv     — mean ± SD of sd_flag, slope_score, cumulative per group × timepoint
averaged_results_summary.csv      — mean ± SD of all per-group metrics across 4 trials
averaged_detection_summary.csv    — detection counts, rates, mean first_detection_min per group

Plots
-----
avg_rg_ratio.png, avg_rb_ratio.png, avg_gb_ratio.png   — ratio time-series with ±1 SD shading
avg_sd_cv.png                                           — CV time-series per channel with ±1 SD
avg_slope_scores.png                                    — slope weighted score with ±1 SD
avg_likert_scores.png                                   — Likert level with ±1 SD
avg_total_detection.png                                 — cumulative detection score with ±1 SD
avg_performance_summary.png                             — sensitivity/specificity bar chart
"""

import os
import csv
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

# ── Config ────────────────────────────────────────────────────────────────────

BASE_DIR     = os.path.dirname(os.path.abspath(__file__))
TRIAL_DIRS   = [os.path.join(BASE_DIR, f"trial-{i}", "results") for i in range(1, 5)]
OUT_DIR      = os.path.join(BASE_DIR, "averaged_results")
N_TRIALS     = 4

# Display order for groups in plots / CSVs
GROUP_ORDER = [
    "0 (Sterile Negative Control)",
    "10^1", "10^2", "10^3", "10^4",
    "10^5", "10^6", "10^7",
    "10^8 (Positive Control)",
]

# Clinical positive threshold
UTI_POSITIVE_CFU = 1e5

# Colour palette — one colour per group for all plots
PALETTE = {
    "0 (Sterile Negative Control)": "#555555",
    "10^1":                          "#4477AA",
    "10^2":                          "#66CCEE",
    "10^3":                          "#228833",
    "10^4":                          "#CCBB44",
    "10^5":                          "#EE6677",
    "10^6":                          "#AA3377",
    "10^7":                          "#BB5566",
    "10^8 (Positive Control)":       "#000000",
}

os.makedirs(OUT_DIR, exist_ok=True)
print(f"Output → {OUT_DIR}")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _read_csv_skip_comments(path: str) -> list[dict]:
    """Read a CSV, skipping lines that start with '#'. Strips group name whitespace."""
    with open(path, newline="") as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    rows = []
    for row in reader:
        row["group"] = row["group"].strip()
        rows.append(row)
    return rows


def _read_csv_section(path: str, section_marker: str) -> list[dict]:
    """
    Read a CSV that has multiple sections separated by comment blocks.
    Returns rows from the section whose header line follows `section_marker`.
    """
    with open(path, newline="") as f:
        lines = f.readlines()

    # Find the section
    target = -1
    for i, line in enumerate(lines):
        if section_marker in line:
            target = i
            break
    if target == -1:
        return []

    # Find the header row after the marker
    header_idx = -1
    for i in range(target + 1, len(lines)):
        if not lines[i].startswith("#") and lines[i].strip():
            header_idx = i
            break
    if header_idx == -1:
        return []

    # Collect rows until the next comment block or EOF
    data_lines = []
    for i in range(header_idx, len(lines)):
        if lines[i].startswith("#") and i > header_idx:
            break
        data_lines.append(lines[i])

    reader = csv.DictReader(data_lines)
    rows = []
    for row in reader:
        if "group" in row:
            row["group"] = row["group"].strip()
        rows.append(row)
    return rows


def _pool(trial_rows_list: list[list[dict]],
          key_cols: list[str],
          numeric_cols: list[str]) -> dict:
    """
    Pool rows from N trials.

    Returns dict keyed by tuple(row[k] for k in key_cols) →
        { col: [val_trial1, val_trial2, ...] }

    Rows are skipped if any key column is missing or if numeric key columns
    (anything other than 'group') cannot be parsed as float — this filters
    out stray summary-section header lines embedded in CSVs.
    """
    pool = defaultdict(lambda: defaultdict(list))
    for rows in trial_rows_list:
        for row in rows:
            # Build key — skip row if a numeric key column is non-parseable
            try:
                key_parts = []
                for k in key_cols:
                    v = row.get(k, "")
                    if v is None:
                        raise ValueError(f"None for key {k}")
                    v = str(v).strip()
                    if k != "group":          # group is kept as string
                        float(v)             # raises ValueError if not numeric
                    key_parts.append(v if k == "group" else v)
                key = tuple(key_parts)
            except (ValueError, KeyError, TypeError):
                continue                     # skip summary / header ghost rows

            for col in numeric_cols:
                try:
                    val = row.get(col)
                    if val is not None and str(val).strip() != "":
                        pool[key][col].append(float(val))
                except (ValueError, KeyError, TypeError):
                    pass
    return pool


def _write_averaged_csv(path: str,
                        pool: dict,
                        key_cols: list[str],
                        numeric_cols: list[str],
                        key_order: list[tuple] | None = None):
    """Write mean ± SD rows to a CSV."""
    header = list(key_cols)
    for col in numeric_cols:
        header += [f"{col}_mean", f"{col}_sd", f"{col}_N"]

    keys = key_order if key_order else sorted(pool.keys())

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for key in keys:
            if key not in pool:
                continue
            row = list(key)
            for col in numeric_cols:
                vals = pool[key][col]
                if vals:
                    row += [round(np.mean(vals), 6),
                            round(np.std(vals, ddof=1) if len(vals) > 1 else 0.0, 6),
                            len(vals)]
                else:
                    row += ["", "", 0]
            writer.writerow(row)
    print(f"  Wrote {path}")


def _sort_key_for_group(name: str) -> int:
    try:
        return GROUP_ORDER.index(name)
    except ValueError:
        return len(GROUP_ORDER)


# ── 1. RGB time-series ────────────────────────────────────────────────────────
print("\n[1/7] Averaging RGB time-series...")

rgb_rows_per_trial = [
    _read_csv_skip_comments(os.path.join(d, "rgb_timeseries.csv"))
    for d in TRIAL_DIRS
]

rgb_numeric = ["R", "G", "B", "RG_ratio", "RB_ratio", "GB_ratio"]
rgb_pool = _pool(rgb_rows_per_trial, ["group", "concentration_cfu_ml", "time_min"], rgb_numeric)

# Build sorted key order: group order first, then timepoint
all_rgb_keys = list(rgb_pool.keys())
sorted_rgb_keys = sorted(all_rgb_keys,
    key=lambda k: (_sort_key_for_group(k[0]), float(k[2])))

_write_averaged_csv(
    os.path.join(OUT_DIR, "averaged_rgb_timeseries.csv"),
    rgb_pool,
    ["group", "concentration_cfu_ml", "time_min"],
    rgb_numeric,
    key_order=sorted_rgb_keys,
)


# ── 2. SD channel values ──────────────────────────────────────────────────────
print("[2/7] Averaging SD channel values...")

sd_rows_per_trial = [
    _read_csv_skip_comments(os.path.join(d, "sd_channel_values.csv"))
    for d in TRIAL_DIRS
]

sd_numeric = ["cv_R", "cv_G", "cv_B", "cv_RG", "cv_RB", "cv_GB", "sd_flag"]
sd_pool = _pool(sd_rows_per_trial, ["group", "concentration_cfu_ml", "time_min"], sd_numeric)

sorted_sd_keys = sorted(sd_pool.keys(),
    key=lambda k: (_sort_key_for_group(k[0]), float(k[2])))

_write_averaged_csv(
    os.path.join(OUT_DIR, "averaged_sd_channel_values.csv"),
    sd_pool,
    ["group", "concentration_cfu_ml", "time_min"],
    sd_numeric,
    key_order=sorted_sd_keys,
)


# ── 3. Slope scores ───────────────────────────────────────────────────────────
print("[3/7] Averaging slope scores...")

slope_rows_per_trial = [
    _read_csv_skip_comments(os.path.join(d, "slope_scores.csv"))
    for d in TRIAL_DIRS
]

slope_numeric = ["weighted_score"]
slope_pool = _pool(slope_rows_per_trial,
                   ["group", "concentration_cfu_ml", "time_min"], slope_numeric)

sorted_slope_keys = sorted(slope_pool.keys(),
    key=lambda k: (_sort_key_for_group(k[0]), float(k[2])))

_write_averaged_csv(
    os.path.join(OUT_DIR, "averaged_slope_scores.csv"),
    slope_pool,
    ["group", "concentration_cfu_ml", "time_min"],
    slope_numeric,
    key_order=sorted_slope_keys,
)


# ── 4. Likert scores ──────────────────────────────────────────────────────────
print("[4/7] Averaging Likert scores...")

likert_rows_per_trial = []
for d in TRIAL_DIRS:
    rows = _read_csv_section(
        os.path.join(d, "likert_scores.csv"),
        "LIKERT SCORES"
    )
    likert_rows_per_trial.append(rows)

likert_numeric = ["likert_score", "match_distance_rgb"]
likert_pool = _pool(likert_rows_per_trial,
                    ["group", "concentration_cfu_ml", "time_min"], likert_numeric)

sorted_likert_keys = sorted(likert_pool.keys(),
    key=lambda k: (_sort_key_for_group(k[0]), float(k[2])))

_write_averaged_csv(
    os.path.join(OUT_DIR, "averaged_likert_scores.csv"),
    likert_pool,
    ["group", "concentration_cfu_ml", "time_min"],
    likert_numeric,
    key_order=sorted_likert_keys,
)


# ── 5. Total detection score time-series ──────────────────────────────────────
print("[5/7] Averaging total detection timeseries...")

tds_rows_per_trial = [
    _read_csv_section(os.path.join(d, "total_detection_score.csv"), "SECTION 1")
    for d in TRIAL_DIRS
]

tds_numeric = ["sd_flag", "slope_score", "contribution", "cumulative_total"]
tds_pool = _pool(tds_rows_per_trial,
                 ["group", "concentration_cfu_ml", "time_min"], tds_numeric)

sorted_tds_keys = sorted(tds_pool.keys(),
    key=lambda k: (_sort_key_for_group(k[0]), float(k[2])))

_write_averaged_csv(
    os.path.join(OUT_DIR, "averaged_total_timeseries.csv"),
    tds_pool,
    ["group", "concentration_cfu_ml", "time_min"],
    tds_numeric,
    key_order=sorted_tds_keys,
)


# ── 6. Per-group results summary ──────────────────────────────────────────────
print("[6/7] Averaging per-group results summary...")

summary_rows_per_trial = [
    _read_csv_section(os.path.join(d, "total_detection_score.csv"), "SECTION 2")
    for d in TRIAL_DIRS
]

summary_numeric = [
    "sd_event_count", "slope_total_score", "likert_event_count", "total_detection_score",
    "rg_auc", "rb_auc", "gb_auc", "rate_of_change_rg",
]
summary_bool = ["sd_detected", "slope_detected", "likert_detected"]

summary_pool = _pool(summary_rows_per_trial,
                     ["group", "concentration_cfu_ml"], summary_numeric)

# Also pool boolean detection as 0/1 to get detection rate
bool_pool = defaultdict(lambda: defaultdict(list))
for rows in summary_rows_per_trial:
    for row in rows:
        key = (row["group"].strip(), row["concentration_cfu_ml"])
        for col in summary_bool:
            val = row.get(col, "").strip().lower()
            if val in ("true", "false"):
                bool_pool[key][col].append(1 if val == "true" else 0)

sorted_summary_keys = sorted(
    set(list(summary_pool.keys()) + list(bool_pool.keys())),
    key=lambda k: _sort_key_for_group(k[0])
)

out_path = os.path.join(OUT_DIR, "averaged_results_summary.csv")
header = ["group", "concentration_cfu_ml", "N_trials", "ground_truth_positive"]
for col in summary_numeric:
    header += [f"{col}_mean", f"{col}_sd"]
for col in summary_bool:
    header += [f"{col}_detection_rate"]
header += ["time_to_ratio_threshold_mean", "time_to_ratio_threshold_sd"]

# time_to_ratio_threshold is "never" for some rows — handle separately
ttr_pool = defaultdict(list)
for rows in summary_rows_per_trial:
    for row in rows:
        key = (row["group"].strip(), row["concentration_cfu_ml"])
        val = row.get("time_to_ratio_threshold_min", "").strip()
        if val and val.lower() != "never":
            try:
                ttr_pool[key].append(float(val))
            except ValueError:
                pass

with open(out_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    for key in sorted_summary_keys:
        if key not in summary_pool and key not in bool_pool:
            continue
        group, conc = key
        try:
            gt = float(conc) >= UTI_POSITIVE_CFU
        except ValueError:
            gt = False

        all_vals = summary_pool[key]
        n = max((len(v) for v in all_vals.values()), default=0)

        row_out = [group, conc, n, gt]
        for col in summary_numeric:
            vals = all_vals[col]
            if vals:
                row_out += [round(np.mean(vals), 4),
                            round(np.std(vals, ddof=1) if len(vals) > 1 else 0.0, 4)]
            else:
                row_out += ["", ""]
        for col in summary_bool:
            vals = bool_pool[key][col]
            row_out.append(round(np.mean(vals), 4) if vals else "")

        ttr_vals = ttr_pool[key]
        if ttr_vals:
            row_out += [round(np.mean(ttr_vals), 2),
                        round(np.std(ttr_vals, ddof=1) if len(ttr_vals) > 1 else 0.0, 2)]
        else:
            row_out += ["never", ""]

        writer.writerow(row_out)

print(f"  Wrote {out_path}")


# ── 7. Detection summary (first_detection_min, per-algorithm) ─────────────────
det_rows_per_trial = []
for d in TRIAL_DIRS:
    path = os.path.join(d, "detection_events.csv")
    with open(path, newline="") as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    rows = [{k: v.strip() for k, v in row.items()} for row in reader]
    det_rows_per_trial.append(rows)

det_numeric = ["sd_first_detection_min", "slope_first_detection_min",
               "likert_first_detection_min", "total_detection_score",
               "slope_total_score", "sd_event_count", "likert_event_count"]
det_bool = ["sd_detected", "slope_detected", "likert_detected"]

det_pool = defaultdict(lambda: defaultdict(list))
det_bool_pool = defaultdict(lambda: defaultdict(list))

for rows in det_rows_per_trial:
    for row in rows:
        key = (row["group"], row["concentration_cfu_ml"])
        for col in det_numeric:
            val = row.get(col, "").strip()
            if val:
                try:
                    det_pool[key][col].append(float(val))
                except ValueError:
                    pass
        for col in det_bool:
            val = row.get(col, "").strip().lower()
            if val in ("true", "false"):
                det_bool_pool[key][col].append(1 if val == "true" else 0)

sorted_det_keys = sorted(
    set(list(det_pool.keys()) + list(det_bool_pool.keys())),
    key=lambda k: _sort_key_for_group(k[0])
)

det_out = os.path.join(OUT_DIR, "averaged_detection_summary.csv")
det_header = ["group", "concentration_cfu_ml", "ground_truth_positive",
              "N_trials",
              "sd_detection_rate", "sd_first_detection_min_mean", "sd_first_detection_min_sd",
              "slope_detection_rate", "slope_first_detection_min_mean", "slope_first_detection_min_sd",
              "likert_detection_rate", "likert_first_detection_min_mean", "likert_first_detection_min_sd",
              "combined_detection_rate",
              "total_detection_score_mean", "total_detection_score_sd"]

with open(det_out, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(det_header)
    for key in sorted_det_keys:
        group, conc = key
        try:
            gt = float(conc) >= UTI_POSITIVE_CFU
        except ValueError:
            gt = False

        n = max((len(v) for v in det_bool_pool[key].values()), default=N_TRIALS)

        sd_rate   = round(np.mean(det_bool_pool[key]["sd_detected"]),   4) if det_bool_pool[key]["sd_detected"]    else 0
        slp_rate  = round(np.mean(det_bool_pool[key]["slope_detected"]), 4) if det_bool_pool[key]["slope_detected"] else 0
        lik_rate  = round(np.mean(det_bool_pool[key]["likert_detected"]),4) if det_bool_pool[key]["likert_detected"] else 0

        # combined = detected by any algorithm in that trial
        combined_detections = []
        for rows in det_rows_per_trial:
            for row in rows:
                if row["group"] == group:
                    sd  = row.get("sd_detected","").strip().lower() == "true"
                    slp = row.get("slope_detected","").strip().lower() == "true"
                    lik = row.get("likert_detected","").strip().lower() == "true"
                    combined_detections.append(1 if (sd or slp or lik) else 0)
        comb_rate = round(np.mean(combined_detections), 4) if combined_detections else 0

        def _mean_sd(vals):
            if not vals:
                return "", ""
            return round(np.mean(vals), 2), round(np.std(vals, ddof=1) if len(vals) > 1 else 0.0, 2)

        sd_fdm, sd_fds    = _mean_sd(det_pool[key]["sd_first_detection_min"])
        slp_fdm, slp_fds  = _mean_sd(det_pool[key]["slope_first_detection_min"])
        lik_fdm, lik_fds  = _mean_sd(det_pool[key]["likert_first_detection_min"])
        tds_m, tds_s      = _mean_sd(det_pool[key]["total_detection_score"])

        writer.writerow([
            group, conc, gt, n,
            sd_rate,  sd_fdm,  sd_fds,
            slp_rate, slp_fdm, slp_fds,
            lik_rate, lik_fdm, lik_fds,
            comb_rate,
            tds_m, tds_s,
        ])

print(f"  Wrote {det_out}")


# ── Plots ─────────────────────────────────────────────────────────────────────
print("[7/7] Generating averaged plots...")
DPI = 150
ALPHA_BAND = 0.18

groups_ordered = GROUP_ORDER

def _extract_series(pool, group, conc, key_col, mean_col, sd_col=None):
    """Pull time-series mean (and optional SD) for a specific group."""
    times, means, sds = [], [], []
    for key, vals in pool.items():
        if key[0] != group:
            continue
        try:
            t = float(key[2])
        except (IndexError, ValueError):
            continue
        m_vals = vals.get(mean_col, [])
        if not m_vals:
            continue
        times.append(t)
        means.append(np.mean(m_vals))
        if sd_col:
            sds.append(np.std(m_vals, ddof=1) if len(m_vals) > 1 else 0.0)
    order = np.argsort(times)
    t_arr = np.array(times)[order]
    m_arr = np.array(means)[order]
    s_arr = np.array(sds)[order] if sds else None
    return t_arr, m_arr, s_arr


# ── Plot helper ────────────────────────────────────────────────────────────────
def _ratio_plot(pool, ratio_col, ylabel, title, fname,
                positive_groups=None, threshold_label=None):
    fig, ax = plt.subplots(figsize=(10, 5))
    for group in groups_ordered:
        color = PALETTE.get(group, "#888888")
        ls = "--" if "Sterile" in group or "Positive" in group else "-"
        lw = 1.4 if "Sterile" in group or "Positive" in group else 1.8
        t_arr, m_arr, s_arr = _extract_series(rgb_pool, group, None, "time_min",
                                               ratio_col, ratio_col)
        if len(t_arr) == 0:
            continue
        ax.plot(t_arr, m_arr, color=color, linestyle=ls, linewidth=lw, label=group)
        if s_arr is not None:
            ax.fill_between(t_arr, m_arr - s_arr, m_arr + s_arr,
                            color=color, alpha=ALPHA_BAND)

    ax.axhline(1.0, color="#AAAAAA", linewidth=0.8, linestyle=":")
    if positive_groups:
        ax.axhline(1.5, color="#CC0000", linewidth=0.9, linestyle="--",
                   label="R/G = 1.50 threshold")
    ax.set_xlabel("Time (min)", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.legend(fontsize=7, ncol=2, loc="upper left")
    ax.set_xlim(0, 120)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, fname), dpi=DPI)
    plt.close(fig)
    print(f"  Saved {fname}")


_ratio_plot(rgb_pool, "RG_ratio", "Mean normalized R/G ratio ± 1 SD",
            "Averaged R/G Ratio — All 4 Trials", "avg_rg_ratio.png",
            positive_groups=True)

_ratio_plot(rgb_pool, "RB_ratio", "Mean normalized R/B ratio ± 1 SD",
            "Averaged R/B Ratio — All 4 Trials", "avg_rb_ratio.png")

_ratio_plot(rgb_pool, "GB_ratio", "Mean normalized G/B ratio ± 1 SD",
            "Averaged G/B Ratio — All 4 Trials", "avg_gb_ratio.png")


# ── SD CV plot ─────────────────────────────────────────────────────────────────
def _cv_plot(pool, cv_col, ylabel, title, fname):
    fig, ax = plt.subplots(figsize=(10, 5))
    for group in groups_ordered:
        color = PALETTE.get(group, "#888888")
        ls = "--" if "Sterile" in group or "Positive" in group else "-"
        t_arr, m_arr, s_arr = _extract_series(pool, group, None, "time_min",
                                               cv_col, cv_col)
        if len(t_arr) == 0:
            continue
        ax.plot(t_arr, m_arr, color=color, linestyle=ls, linewidth=1.6, label=group)
        if s_arr is not None:
            ax.fill_between(t_arr, m_arr - s_arr, m_arr + s_arr,
                            color=color, alpha=ALPHA_BAND)
    ax.axhline(0.040, color="#CC4400", linewidth=0.9, linestyle="--", label="thr_all = 0.040")
    ax.axhline(0.080, color="#880000", linewidth=0.9, linestyle=":",  label="thr_two = 0.080")
    ax.set_xlabel("Time (min)", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.legend(fontsize=7, ncol=2, loc="upper left")
    ax.set_xlim(0, 120)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, fname), dpi=DPI)
    plt.close(fig)
    print(f"  Saved {fname}")

_cv_plot(sd_pool, "cv_RG", "Mean CV (R/G) ± 1 SD",
         "Averaged SD CV — R/G channel — All 4 Trials", "avg_sd_cv_rg.png")
_cv_plot(sd_pool, "cv_RB", "Mean CV (R/B) ± 1 SD",
         "Averaged SD CV — R/B channel — All 4 Trials", "avg_sd_cv_rb.png")
_cv_plot(sd_pool, "cv_GB", "Mean CV (G/B) ± 1 SD",
         "Averaged SD CV — G/B channel — All 4 Trials", "avg_sd_cv_gb.png")


# ── Slope score plot ──────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))
for group in groups_ordered:
    color = PALETTE.get(group, "#888888")
    ls = "--" if "Sterile" in group or "Positive" in group else "-"
    t_arr, m_arr, s_arr = _extract_series(slope_pool, group, None, "time_min",
                                           "weighted_score", "weighted_score")
    if len(t_arr) == 0:
        continue
    ax.plot(t_arr, m_arr, color=color, linestyle=ls, linewidth=1.6, label=group)
    if s_arr is not None:
        ax.fill_between(t_arr, m_arr - s_arr, m_arr + s_arr,
                        color=color, alpha=ALPHA_BAND)
ax.axhline(4.0, color="#CC0000", linewidth=0.9, linestyle="--", label="threshold > 4.0")
ax.set_xlabel("Time (min)", fontsize=11)
ax.set_ylabel("Mean weighted slope score ± 1 SD", fontsize=11)
ax.set_title("Averaged Slope Scores — All 4 Trials", fontsize=13, fontweight="bold")
ax.legend(fontsize=7, ncol=2, loc="upper left")
ax.set_xlim(0, 120)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "avg_slope_scores.png"), dpi=DPI)
plt.close(fig)
print("  Saved avg_slope_scores.png")


# ── Likert score plot ─────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))
for group in groups_ordered:
    color = PALETTE.get(group, "#888888")
    ls = "--" if "Sterile" in group or "Positive" in group else "-"
    t_arr, m_arr, s_arr = _extract_series(likert_pool, group, None, "time_min",
                                           "likert_score", "likert_score")
    if len(t_arr) == 0:
        continue
    ax.plot(t_arr, m_arr, color=color, linestyle=ls, linewidth=1.6, label=group)
    if s_arr is not None:
        ax.fill_between(t_arr, m_arr - s_arr, m_arr + s_arr,
                        color=color, alpha=ALPHA_BAND)
ax.axhline(4.0, color="#CC0000", linewidth=0.9, linestyle="--", label="positive threshold ≥ 4")
ax.set_yticks([1, 2, 3, 4, 5])
ax.set_yticklabels(["1\nNegative", "2\nLikely−", "3\nAmbig.", "4\nLikely+", "5\nPositive"])
ax.set_xlabel("Time (min)", fontsize=11)
ax.set_ylabel("Mean Likert level ± 1 SD", fontsize=11)
ax.set_title("Averaged Likert Scores — All 4 Trials", fontsize=13, fontweight="bold")
ax.legend(fontsize=7, ncol=2, loc="upper left")
ax.set_xlim(0, 120)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "avg_likert_scores.png"), dpi=DPI)
plt.close(fig)
print("  Saved avg_likert_scores.png")


# ── Cumulative detection score plot ───────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))
for group in groups_ordered:
    color = PALETTE.get(group, "#888888")
    ls = "--" if "Sterile" in group or "Positive" in group else "-"
    t_arr, m_arr, s_arr = _extract_series(tds_pool, group, None, "time_min",
                                           "cumulative_total", "cumulative_total")
    if len(t_arr) == 0:
        continue
    ax.plot(t_arr, m_arr, color=color, linestyle=ls, linewidth=1.6, label=group)
    if s_arr is not None:
        ax.fill_between(t_arr, m_arr - s_arr, m_arr + s_arr,
                        color=color, alpha=ALPHA_BAND)
ax.set_xlabel("Time (min)", fontsize=11)
ax.set_ylabel("Mean cumulative detection score ± 1 SD", fontsize=11)
ax.set_title("Averaged Cumulative Detection Score — All 4 Trials", fontsize=13, fontweight="bold")
ax.legend(fontsize=7, ncol=2, loc="upper left")
ax.set_xlim(0, 120)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "avg_total_detection.png"), dpi=DPI)
plt.close(fig)
print("  Saved avg_total_detection.png")


# ── Performance summary bar chart ─────────────────────────────────────────────
# Load averaged detection summary for rates
det_data = {}
with open(det_out, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        det_data[row["group"]] = row

algos = ["SD", "Slope", "Likert", "Combined"]
algo_keys = {
    "SD":       "sd_detection_rate",
    "Slope":    "slope_detection_rate",
    "Likert":   "likert_detection_rate",
    "Combined": "combined_detection_rate",
}

tp_groups  = [g for g in GROUP_ORDER if float(det_data.get(g, {}).get("concentration_cfu_ml", 0) or 0) >= UTI_POSITIVE_CFU]
neg_groups = [g for g in GROUP_ORDER if float(det_data.get(g, {}).get("concentration_cfu_ml", 0) or 0) < UTI_POSITIVE_CFU]

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle("Averaged Detection Rates Across 4 Trials", fontsize=13, fontweight="bold")

bar_colors = ["#4477AA", "#228833", "#AA3377", "#CC3300"]

for ax_i, (ax, subset, subtitle) in enumerate(zip(
        axes,
        [tp_groups, neg_groups],
        ["True Positives (≥10⁵ CFU/mL)", "True Negatives (<10⁵ CFU/mL)"])):
    x = np.arange(len(subset))
    w = 0.18
    for ai, (algo, key) in enumerate(algo_keys.items()):
        vals = []
        for g in subset:
            row = det_data.get(g, {})
            try:
                vals.append(float(row.get(key, 0) or 0) * 100)
            except ValueError:
                vals.append(0.0)
        ax.bar(x + ai * w, vals, w, label=algo, color=bar_colors[ai], alpha=0.85)

    ax.set_xticks(x + 1.5 * w)
    ax.set_xticklabels([g.replace(" (", "\n(") for g in subset], fontsize=8)
    ax.set_ylim(0, 110)
    ax.set_ylabel("Detection rate (%)", fontsize=10)
    ax.set_title(subtitle, fontsize=10)
    ax.axhline(100, color="#777777", linewidth=0.5, linestyle=":")
    if ax_i == 0:
        ax.legend(fontsize=8, loc="lower right")

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "avg_performance_summary.png"), dpi=DPI)
plt.close(fig)
print("  Saved avg_performance_summary.png")


# ── AUC bar chart (mean ± SD per group) ───────────────────────────────────────
auc_data = {}
with open(os.path.join(OUT_DIR, "averaged_results_summary.csv"), newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        auc_data[row["group"]] = row

fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
fig.suptitle("Mean AUC per Concentration Group — Averaged Across 4 Trials",
             fontsize=12, fontweight="bold")

for ax, metric, label in zip(
        axes,
        ["rg_auc", "rb_auc", "gb_auc"],
        ["R/G AUC (ratio·min)", "R/B AUC (ratio·min)", "G/B AUC (ratio·min)"]):
    groups_plot = [g for g in GROUP_ORDER if g in auc_data]
    means = []
    errs  = []
    colors = []
    for g in groups_plot:
        row = auc_data[g]
        try:
            means.append(float(row.get(f"{metric}_mean") or 0))
            errs.append(float(row.get(f"{metric}_sd") or 0))
        except ValueError:
            means.append(0); errs.append(0)
        try:
            gt = float(row.get("concentration_cfu_ml", 0)) >= UTI_POSITIVE_CFU
        except ValueError:
            gt = False
        colors.append("#CC3300" if gt else "#4477AA")

    x = np.arange(len(groups_plot))
    ax.bar(x, means, yerr=errs, color=colors, alpha=0.82, capsize=4, width=0.65)
    ax.set_xticks(x)
    ax.set_xticklabels([g.replace(" (", "\n(") for g in groups_plot], fontsize=7, rotation=30, ha="right")
    ax.set_ylabel(label, fontsize=9)
    ax.set_title(label, fontsize=10)

pos_patch = mpatches.Patch(color="#CC3300", label="Clinical positive (≥10⁵)")
neg_patch = mpatches.Patch(color="#4477AA", label="Clinical negative (<10⁵)")
fig.legend(handles=[pos_patch, neg_patch], loc="lower center", ncol=2,
           fontsize=9, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig(os.path.join(OUT_DIR, "avg_auc_bars.png"), dpi=DPI)
plt.close(fig)
print("  Saved avg_auc_bars.png")


print(f"\n✓ Done. {len(os.listdir(OUT_DIR))} files in {OUT_DIR}")
