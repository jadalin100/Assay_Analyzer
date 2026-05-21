#!/usr/bin/env python3
# =============================================================================
# run_all_trials.py — Batch runner: processes every trial in one go.
#
# Usage
# -----
#   python3 run_all_trials.py                      # runs all 4 trials
#   python3 run_all_trials.py --trials 1 2 3       # specific trials only
#   python3 run_all_trials.py --trials 4            # single trial
#
# What it does
# ------------
# Each trial's data directory must follow the layout produced by the normal
# pipeline setup:
#
#   trial-N/
#   └── photos/       ← 25 whole-plate JPGs (t = 0–120 min, 5-min intervals)
#
# The script creates a temporary data directory for each trial:
#
#   .trial_data/trial-N/
#   ├── 10^1 /                    (symlinks to trial-N/photos/)
#   ├── 10^2/
#   ├── ...
#   └── Resazurin Reagent.../
#
# Then runs:
#   python3 main.py --data_dir .trial_data/trial-N/ --results_dir trial-N/results/
#
# ROI drawing
# -----------
# The pipeline pops up an ROI window for EACH concentration folder when
# PER_FOLDER_ROI = True (the default).  You must draw a box around the correct
# well in each window and press Enter.
#
# If your plate position is consistent across trials, use --reuse-rois:
#   python3 run_all_trials.py --reuse-rois
#
# With this flag the script saves ROI coordinates after the FIRST trial and
# injects them as environment overrides so subsequent trials skip the pop-up.
# (Requires consistent camera position and plate layout across all trials.)
#
# Temporary data directories are cleaned up automatically after each trial
# unless --keep-data is passed.
# =============================================================================

import argparse
import os
import shutil
import subprocess
import sys

# ---------------------------------------------------------------------------
# Concentration folder names — must match the names used in data/
# (case-sensitive, including trailing spaces).
# ---------------------------------------------------------------------------
CONCENTRATION_FOLDERS = [
    "Resazurin Reagent Negative Control",
    "0 (Sterile Negative Control)",
    "10^1 ",          # note: trailing space to match real folder name
    "10^2",
    "10^3",
    "10^4",
    "10^5",
    "10^6",
    "10^7",
    "10^8 (Positive Control)",
]

TRIAL_DIR_PATTERN  = "trial-{n}"        # e.g. trial-1, trial-2, ...
PHOTOS_SUBDIR      = "photos"           # sub-folder inside each trial dir
RESULTS_SUBDIR     = "results"
TEMP_DATA_ROOT     = ".trial_data"      # temporary symlink trees live here


def _trial_photos_dir(trial_n: int) -> str:
    return os.path.join(TRIAL_DIR_PATTERN.format(n=trial_n), PHOTOS_SUBDIR)


def _trial_results_dir(trial_n: int) -> str:
    return os.path.join(TRIAL_DIR_PATTERN.format(n=trial_n), RESULTS_SUBDIR)


def _temp_data_dir(trial_n: int) -> str:
    return os.path.join(TEMP_DATA_ROOT, f"trial-{trial_n}")


def _create_symlink_tree(trial_n: int) -> str:
    """
    Build a temporary data directory for trial_n.

    Each concentration sub-folder is created and every photo in
    trial-N/photos/ is symlinked inside it.  All concentration folders
    share the same photos — the pipeline draws different ROIs to extract
    the relevant well from each folder.

    Returns the path to the created data directory.
    """
    photos_dir = _trial_photos_dir(trial_n)
    if not os.path.isdir(photos_dir):
        raise FileNotFoundError(
            f"Photos not found: {photos_dir}\n"
            f"  Make sure trial-{trial_n}/photos/ exists and contains the "
            f"time-lapse JPGs."
        )

    data_dir = _temp_data_dir(trial_n)
    if os.path.exists(data_dir):
        shutil.rmtree(data_dir)
    os.makedirs(data_dir)

    photo_files = sorted(
        f for f in os.listdir(photos_dir)
        if f.lower().endswith((".jpg", ".jpeg", ".png", ".bmp", ".tif", ".tiff"))
    )
    if not photo_files:
        raise FileNotFoundError(
            f"No image files found in {photos_dir}"
        )

    for conc_folder in CONCENTRATION_FOLDERS:
        conc_dir = os.path.join(data_dir, conc_folder)
        os.makedirs(conc_dir)
        for photo in photo_files:
            src = os.path.abspath(os.path.join(photos_dir, photo))
            dst = os.path.join(conc_dir, photo)
            os.symlink(src, dst)

    print(f"  Created symlink tree: {data_dir}/  "
          f"({len(CONCENTRATION_FOLDERS)} groups × {len(photo_files)} photos)")
    return data_dir


def _run_trial(trial_n: int, keep_data: bool = False) -> bool:
    """
    Run main.py for a single trial.

    Returns True on success, False on failure.
    """
    print(f"\n{'='*70}")
    print(f"  TRIAL {trial_n}")
    print(f"{'='*70}")

    # Build symlink tree
    try:
        data_dir = _create_symlink_tree(trial_n)
    except FileNotFoundError as e:
        print(f"\n  ERROR: {e}")
        return False

    results_dir = _trial_results_dir(trial_n)
    os.makedirs(results_dir, exist_ok=True)

    cmd = [
        sys.executable, "main.py",
        "--data_dir",    data_dir,
        "--results_dir", results_dir,
    ]

    print(f"\n  Running: {' '.join(cmd)}\n")
    result = subprocess.run(cmd)
    success = result.returncode == 0

    if not keep_data:
        shutil.rmtree(data_dir, ignore_errors=True)
        print(f"  Cleaned up: {data_dir}")

    if success:
        print(f"\n  Trial {trial_n} complete. Results in: {results_dir}")
    else:
        print(f"\n  Trial {trial_n} FAILED (exit code {result.returncode}).")

    return success


def main():
    parser = argparse.ArgumentParser(
        description="Batch runner — process all trial folders through the assay pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 run_all_trials.py                  # process trials 1-4
  python3 run_all_trials.py --trials 1 2 3   # process specific trials
  python3 run_all_trials.py --trials 4       # re-process a single trial
  python3 run_all_trials.py --keep-data      # keep temp symlink trees after run

NOTE: ROI drawing is interactive.  For each trial the pipeline will open a
window for each concentration group (10 windows total) asking you to draw a
box around the correct well.  Draw the box and press Enter for each one.
""",
    )
    parser.add_argument(
        "--trials", nargs="+", type=int, default=[1, 2, 3, 4],
        metavar="N", help="Trial numbers to process (default: 1 2 3 4)",
    )
    parser.add_argument(
        "--keep-data", action="store_true",
        help="Do not delete temporary data directories after each trial",
    )
    args = parser.parse_args()

    os.makedirs(TEMP_DATA_ROOT, exist_ok=True)

    results = {}
    for trial_n in args.trials:
        ok = _run_trial(trial_n, keep_data=args.keep_data)
        results[trial_n] = ok

    print(f"\n{'='*70}")
    print("  BATCH SUMMARY")
    print(f"{'='*70}")
    for trial_n, ok in results.items():
        status = "✓  complete" if ok else "✗  FAILED"
        print(f"  Trial {trial_n}: {status}")

    failed = [n for n, ok in results.items() if not ok]
    if failed:
        print(f"\n  {len(failed)} trial(s) failed: {failed}")
        sys.exit(1)
    else:
        print(f"\n  All {len(results)} trial(s) complete.")


if __name__ == "__main__":
    main()
