#!/usr/bin/env python3
# =============================================================================
# cross_trial_stats.py — Pool all trial results and run cross-trial statistics.
#
# Usage
# -----
#   python3 cross_trial_stats.py                  # trials 1–4
#   python3 cross_trial_stats.py --trials 1 2 3   # specific trials
#
# What it does
# ------------
# Reads results_summary.csv from each trial-N/results/ directory, pools
# rows by concentration group across all trials (each trial = one replicate),
# then runs the same ANOVA / Kruskal-Wallis pipeline as stats.py.
#
# Output
# ------
#   results/cross_trial_statistics.json   — full test results
#   results/cross_trial_summary.csv       — one row per group with N, mean ± SD
# =============================================================================

import argparse
import csv
import json
import os
import sys
import numpy as np

# Make sure the package root is on the path
sys.path.insert(0, os.path.dirname(__file__))
import stats as stats_module
import config


def _load_trial(trial_n: int) -> list[dict]:
    """Load results_summary.csv for a single trial. Returns list of row dicts."""
    path = os.path.join(f"trial-{trial_n}", "results", "results_summary.csv")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing: {path}")
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    # Cast numeric fields
    numeric = [
        "concentration_cfu_ml", "total_detection_score",
        "rg_auc", "rb_auc", "gb_auc",
        "time_to_ratio_threshold", "final_RG_ratio", "final_RB_ratio", "final_GB_ratio",
        "rate_of_change_rg",
    ]
    for r in rows:
        r["_trial"] = trial_n
        for col in numeric:
            raw = r.get(col, "")
            try:
                r[col] = float(raw) if raw not in ("", None) else None
            except ValueError:
                r[col] = None
        # Boolean
        for col in ("ground_truth_positive", "overall_positive",
                    "sd_detected", "slope_detected", "likert_detected"):
            r[col] = r.get(col, "").lower() == "true"
    return rows


def _save_summary_csv(groups: dict[str, list[dict]], metrics: list[str], out_path: str):
    """Save a wide-format summary table: group × metric (mean, SD, N)."""
    fieldnames = ["group", "N"] + [f"{m}_mean" for m in metrics] + [f"{m}_sd" for m in metrics]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for group, rows in sorted(groups.items()):
            row = {"group": group, "N": len(rows)}
            for m in metrics:
                vals = [r[m] for r in rows if r.get(m) is not None]
                row[f"{m}_mean"] = float(np.mean(vals)) if vals else None
                row[f"{m}_sd"]   = float(np.std(vals, ddof=1)) if len(vals) > 1 else None
            w.writerow(row)
    print(f"  Summary CSV saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Cross-trial statistical analysis — pool all trials and run ANOVA/KW."
    )
    parser.add_argument(
        "--trials", nargs="+", type=int, default=[1, 2, 3, 4], metavar="N",
        help="Trial numbers to include (default: 1 2 3 4)",
    )
    parser.add_argument(
        "--out-dir", default="cross_trial_results",
        help="Directory to write output files (default: cross_trial_results/)",
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Load and pool all trials
    all_rows: list[dict] = []
    print(f"\nLoading {len(args.trials)} trials...")
    for n in args.trials:
        try:
            rows = _load_trial(n)
            all_rows.extend(rows)
            print(f"  Trial {n}: {len(rows)} groups loaded")
        except FileNotFoundError as e:
            print(f"  WARNING: {e} — skipping trial {n}")

    if not all_rows:
        print("ERROR: No data loaded. Exiting.")
        sys.exit(1)

    # Group by concentration for summary CSV
    conc_groups: dict[str, list[dict]] = {}
    for r in all_rows:
        conc = r.get("concentration_cfu_ml")
        if conc is None:
            continue
        label = f"10^{int(round(np.log10(conc)))}" if conc > 0 else "0 (negative)"
        conc_groups.setdefault(label, []).append(r)

    print(f"\nConcentration groups (N = replicates across trials):")
    for label, rows in sorted(conc_groups.items(),
                               key=lambda x: float(x[1][0]["concentration_cfu_ml"])):
        print(f"  {label:12s}  N={len(rows)}  "
              f"total_score: {[r['total_detection_score'] for r in rows]}  "
              f"final_RG: {[round(r['final_RG_ratio'],3) if r['final_RG_ratio'] else None for r in rows]}")

    # Summary CSV
    summary_metrics = ["total_detection_score", "rg_auc", "rb_auc", "gb_auc",
                       "final_RG_ratio", "final_RB_ratio", "final_GB_ratio", "time_to_ratio_threshold"]
    _save_summary_csv(
        conc_groups, summary_metrics,
        os.path.join(args.out_dir, "cross_trial_summary.csv")
    )

    # Run statistical tests (reuses stats.py pipeline)
    print(f"\nRunning cross-trial statistical tests ({len(all_rows)} total observations):")
    stat_results = stats_module.run_statistics(
        all_rows,
        os.path.join(args.out_dir, "cross_trial_statistics.json")
    )

    # Print readable summary
    print(f"\n{'='*70}")
    print("  CROSS-TRIAL STATISTICAL SUMMARY")
    print(f"{'='*70}")
    for m in stat_results.get("metrics", []):
        test  = m.get("test", "—")
        stat  = m.get("statistic")
        pval  = m.get("p_value")
        sig   = m.get("significant")
        note  = m.get("note", "")
        warn  = "  ⚠ " + m.get("replication_warning", "").split(".")[0] if m.get("replication_warning") else ""
        label = m["label"]
        if pval is not None:
            sig_str = "✓ SIGNIFICANT" if sig else "✗ not significant"
            print(f"\n  {label}")
            print(f"    {test}: stat={stat:.4f}  p={pval:.4f}  {sig_str}{warn}")
        else:
            print(f"\n  {label}")
            print(f"    {note or 'no test performed'}")

    print(f"\n  Output written to: {args.out_dir}/")


if __name__ == "__main__":
    main()
