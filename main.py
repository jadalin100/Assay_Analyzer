# =============================================================================
# main.py — Entry point for the Assay Analyzer pipeline.
# =============================================================================

import os
import re
import sys
import csv
import argparse

import config
from ingestion import load_sample_images, extract_rgb_series, extract_rgb_from_video
from classifier import classify_all
from metrics import compute_metrics
from stats import run_statistics
from visualization import generate_all_plots
from reporting import print_and_save_tables

VIDEO_EXTS = {".mp4", ".avi", ".mov", ".mkv", ".m4v"}
IMAGE_EXTS = {".jpg", ".jpeg", ".png", ".bmp", ".tif", ".tiff"}


def _should_skip(folder_name: str) -> bool:
    """Return True if this folder should be silently ignored (e.g. Reference Card)."""
    name = folder_name.lower()
    return any(skip in name for skip in config.SKIP_FOLDERS)


def _folder_sort_key(name: str) -> tuple:
    """
    Sort folders by their position in config.FOLDER_ORDER.
    Unrecognised folders sort to the end, then alphabetically.
    """
    lower = name.lower().strip()
    for i, prefix in enumerate(config.FOLDER_ORDER):
        if lower.startswith(prefix.lower()):
            return (i, lower)
    return (len(config.FOLDER_ORDER), lower)


def _infer_concentration(folder_name: str) -> float:
    """
    Infer concentration in CFU/mL from the folder name.

    Priority order:
      1. Numeric pattern: any '10^N' found anywhere in the name → 10^N CFU/mL.
         This correctly handles names like 'Low 10^5', 'High 10^7',
         'Positive Control 10^8' regardless of the leading keyword.
      2. Keyword prefix: longest matching key in config.CONCENTRATIONS.
         Used for controls with no numeric concentration (e.g. 'Sterile Negative Control').
    """
    name = folder_name.lower().strip()

    # Step 1 — extract 10^N from anywhere in the name
    match = re.search(r'10\^(\d+)', name)
    if match:
        return 10.0 ** int(match.group(1))

    # Step 2 — keyword prefix fallback (for controls)
    for key in sorted(config.CONCENTRATIONS.keys(), key=len, reverse=True):
        if name.startswith(key.lower()):
            return config.CONCENTRATIONS[key]

    raise ValueError(
        f"Cannot infer concentration from folder name: '{folder_name}'\n"
        f"  Include '10^N' anywhere in the name (e.g. 'Low 10^5', 'Positive Control 10^8'), or\n"
        f"  rename to start with one of: {list(config.CONCENTRATIONS.keys())}\n"
        f"  or add a new entry to CONCENTRATIONS in config.py"
    )


def _load_sample(sample_dir: str, folder_name: str, shared_roi=None, shared_ref_roi=None) -> dict:
    """Load images or video from a sample folder and extract RGB series."""
    files = os.listdir(sample_dir)
    video_files = [f for f in files if os.path.splitext(f)[1].lower() in VIDEO_EXTS]

    if video_files:
        video_path = os.path.join(sample_dir, video_files[0])
        print(f"  Loading video: {video_files[0]}")
        rgb_data = extract_rgb_from_video(
            video_path, roi=shared_roi, folder_name=folder_name, ref_roi=shared_ref_roi
        )
    else:
        images = load_sample_images(sample_dir)
        rgb_data = extract_rgb_series(
            images, roi=shared_roi, folder_name=folder_name, ref_roi=shared_ref_roi
        )

    concentration = _infer_concentration(folder_name)
    return {
        "name": folder_name,
        "concentration": concentration,
        "rgb_data": rgb_data,
    }


def _save_rgb_timeseries_csv(results: list[dict], path: str):
    """
    Save raw RGB values and ratios at every timepoint.
    One row per (group × timepoint): group, time_min, R, G, B, RG_ratio, RB_ratio, GB_ratio.
    """
    fieldnames = ["group", "concentration_cfu_ml", "time_min",
                  "R", "G", "B", "RG_ratio", "RB_ratio", "GB_ratio"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            rgb = r["rgb_data"]
            times = rgb["times"]
            for i, t in enumerate(times):
                writer.writerow({
                    "group":                r["name"],
                    "concentration_cfu_ml": r["concentration_cfu_ml"],
                    "time_min":             int(t),
                    "R":                    round(rgb["R"][i], 4) if i < len(rgb["R"]) else "",
                    "G":                    round(rgb["G"][i], 4) if i < len(rgb["G"]) else "",
                    "B":                    round(rgb["B"][i], 4) if i < len(rgb["B"]) else "",
                    "RG_ratio":             round(rgb["RG_ratio"][i], 6) if rgb.get("RG_ratio") and i < len(rgb["RG_ratio"]) else "",
                    "RB_ratio":             round(rgb["RB_ratio"][i], 6) if rgb.get("RB_ratio") and i < len(rgb["RB_ratio"]) else "",
                    "GB_ratio":             round(rgb["GB_ratio"][i], 6) if rgb.get("GB_ratio") and i < len(rgb["GB_ratio"]) else "",
                })
    print(f"  CSV saved: {path}")


def _save_csv(results: list[dict], path: str):
    """Save a summary CSV with one row per sample."""
    fieldnames = [
        "name", "concentration_cfu_ml", "ground_truth_positive",
        "overall_positive",
        "sd_detected", "sd_first_detection_min",
        "slope_detected", "slope_first_detection_min", "slope_total_score",
        "final_RG_ratio", "final_RB_ratio", "final_GB_ratio",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            slope_scores = r.get("slope_result", {}).get("weighted_scores", [])
            row = {k: r.get(k, "") for k in fieldnames}
            row["overall_positive"] = r.get("sd_detected") or r.get("slope_detected")
            row["slope_total_score"] = sum(slope_scores)
            writer.writerow(row)
    print(f"  CSV saved: {path}")


def main():
    parser = argparse.ArgumentParser(description="Smartphone Colorimetric Assay Analyzer")
    parser.add_argument("--data_dir", default="data/", help="Path to data folder")
    args = parser.parse_args()

    data_dir = args.data_dir
    results_dir = config.RESULTS_DIR
    os.makedirs(results_dir, exist_ok=True)

    # ── Stage 1: Load ─────────────────────────────────────────────────────────
    print("[1/5] Loading samples...")
    if not os.path.isdir(data_dir):
        print(f"ERROR: data directory not found: {data_dir}")
        sys.exit(1)

    all_subdirs = [
        d for d in os.listdir(data_dir)
        if os.path.isdir(os.path.join(data_dir, d))
    ]
    subdirs = sorted(all_subdirs, key=_folder_sort_key)
    if not subdirs:
        print(f"ERROR: No sample subdirectories found in: {data_dir}")
        sys.exit(1)

    # Print the order so the user can confirm before any pop-ups appear
    sample_folders = [f for f in subdirs if not _should_skip(f)]
    print("\n  Processing order:")
    for i, f in enumerate(sample_folders, 1):
        print(f"    {i}. {f}")
    print()

    per_folder = getattr(config, "PER_FOLDER_ROI", False)
    samples = []
    shared_roi = None
    shared_ref_roi = None   # reference card ROI — asked once, reused for all groups

    for folder in subdirs:
        if _should_skip(folder):
            print(f"\n  Skipping:  {folder}  (not a sample folder)")
            continue
        folder_path = os.path.join(data_dir, folder)
        print(f"\n  Loading: {folder}")
        if per_folder:
            # Each folder gets its own well ROI drawn on that folder's first image.
            roi_for_this_folder = None
        else:
            # Reuse the single shared well ROI drawn on the very first folder.
            roi_for_this_folder = shared_roi
        try:
            sample = _load_sample(
                folder_path, folder,
                shared_roi=roi_for_this_folder,
                shared_ref_roi=shared_ref_roi,
            )
            # In global-shared mode, capture the well ROI from the first folder only
            if not config.PER_IMAGE_ROI and not per_folder and shared_roi is None:
                shared_roi = sample["rgb_data"]["roi"]
            # Capture reference card ROI from the first group and reuse for all others
            if shared_ref_roi is None and sample["rgb_data"].get("ref_roi") is not None:
                shared_ref_roi = sample["rgb_data"]["ref_roi"]
            samples.append(sample)
        except ValueError as e:
            if "No ROI selected" in str(e):
                print(
                    f"\n  !! ERROR for group '{folder}': You pressed ENTER without drawing a box.\n"
                    f"     Please re-run the program and draw a rectangle around the correct well\n"
                    f"     when the '{folder}' pop-up appears, then press ENTER.\n"
                )
                sys.exit(1)
            else:
                print(f"  WARNING: Skipping '{folder}': {e}")
        except Exception as e:
            print(f"  WARNING: Skipping '{folder}': {e}")

    if not samples:
        print("ERROR: No samples loaded successfully.")
        sys.exit(1)

    print(f"\n  Loaded {len(samples)} samples.")

    # ── Stage 2: Classify ─────────────────────────────────────────────────────
    print("\n[2/5] Running detection algorithms...")
    results = classify_all(samples)

    # ── Stage 3: Metrics ──────────────────────────────────────────────────────
    print("\n[3/5] Computing performance metrics...")
    metrics = compute_metrics(results)

    # ── Stage 4: Statistics ───────────────────────────────────────────────────
    print("\n[4/5] Running statistical analysis...")
    stats_path = os.path.join(results_dir, "statistics.json")
    run_statistics(results, stats_path)

    # ── Stage 5: Plots & Save ─────────────────────────────────────────────────
    print("\n[5/5] Generating plots and saving results...")
    generate_all_plots(results, metrics, results_dir)
    csv_path = os.path.join(results_dir, "results_summary.csv")
    _save_csv(results, csv_path)
    rgb_csv_path = os.path.join(results_dir, "rgb_timeseries.csv")
    _save_rgb_timeseries_csv(results, rgb_csv_path)

    # ── Stage 6: Per-group tables ─────────────────────────────────────────────
    print("\n[+] Per-group algorithm tables:")
    print_and_save_tables(results, results_dir)

    print(f"\nDone. All results in: {results_dir}")


if __name__ == "__main__":
    main()
