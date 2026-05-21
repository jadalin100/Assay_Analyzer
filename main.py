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


def _is_reagent_reference(folder_name: str) -> bool:
    """
    Return True if this folder is the resazurin reagent control reference.
    When USE_REAGENT_REFERENCE is True this folder is loaded for its RGB series
    but excluded from classification — its values normalise all other groups.
    """
    if not getattr(config, "USE_REAGENT_REFERENCE", False):
        return False
    prefix = getattr(config, "REAGENT_REFERENCE_FOLDER", "").lower().strip()
    return bool(prefix) and folder_name.lower().strip().startswith(prefix)


def _normalize_to_reagent_reference(samples: list[dict], ref_rgb: dict) -> None:
    """
    Divide each sample's R/G, R/B, and G/B ratio values by the reagent control's
    value at the same timepoint.  Edits rgb_data in-place.

    After normalization:
      - A stable (non-converting) group stays near 1.0 throughout the assay.
      - A colour-changing group rises above 1.0 as bacteria reduce resazurin.
      - Lighting drift and auto-exposure changes cancel by construction because
        both the group and the reference are affected equally.

    Timepoints present in the group but missing from the reference are left
    unnormalized (rare edge case — both folders should share the same image set).
    """
    ref_by_time: dict[float, int] = {t: idx for idx, t in enumerate(ref_rgb["times"])}
    ratio_keys = [k for k in ("RG_ratio", "RB_ratio", "GB_ratio") if ref_rgb.get(k)]

    for sample in samples:
        rgb = sample["rgb_data"]
        for i, t in enumerate(rgb["times"]):
            ref_i = ref_by_time.get(t)
            if ref_i is None:
                continue
            for key in ratio_keys:
                if not rgb.get(key):
                    continue
                ref_val = ref_rgb[key][ref_i]
                if ref_val and ref_val != 0.0:
                    rgb[key][i] = rgb[key][i] / ref_val


def _save_reagent_reference_csv(ref_rgb: dict, ref_name: str, path: str):
    """Save the raw (pre-normalization) reagent reference time series."""
    fieldnames = ["group", "time_min", "R", "G", "B", "RG_ratio", "RB_ratio", "GB_ratio"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        times = ref_rgb["times"]
        for i, t in enumerate(times):
            writer.writerow({
                "group":    ref_name,
                "time_min": int(t),
                "R":        round(ref_rgb["R"][i], 4)          if i < len(ref_rgb.get("R", [])) else "",
                "G":        round(ref_rgb["G"][i], 4)          if i < len(ref_rgb.get("G", [])) else "",
                "B":        round(ref_rgb["B"][i], 4)          if i < len(ref_rgb.get("B", [])) else "",
                "RG_ratio": round(ref_rgb["RG_ratio"][i], 6)   if ref_rgb.get("RG_ratio") and i < len(ref_rgb["RG_ratio"]) else "",
                "RB_ratio": round(ref_rgb["RB_ratio"][i], 6)   if ref_rgb.get("RB_ratio") and i < len(ref_rgb["RB_ratio"]) else "",
                "GB_ratio": round(ref_rgb["GB_ratio"][i], 6)   if ref_rgb.get("GB_ratio") and i < len(ref_rgb["GB_ratio"]) else "",
            })
    print(f"  CSV saved: {path}")


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

    # Step 2 — folder starts with a bare "0" → 0 CFU/mL (e.g. "0 (Sterile Negative Control)")
    if re.match(r'^0\b', name):
        return 0.0

    # Step 3 — keyword prefix fallback (for controls)
    for key in sorted(config.CONCENTRATIONS.keys(), key=len, reverse=True):
        if name.startswith(key.lower()):
            return config.CONCENTRATIONS[key]

    raise ValueError(
        f"Cannot infer concentration from folder name: '{folder_name}'\n"
        f"  Include '10^N' anywhere in the name (e.g. 'Low 10^5', 'Positive Control 10^8'), or\n"
        f"  start with '0' for a sterile/zero-concentration control, or\n"
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


def _save_rgb_timeseries_csv(results: list[dict], path: str,
                             reagent_normalized: bool = False):
    """
    Save RGB values and ratios at every timepoint.
    One row per (group × timepoint).

    The 'reagent_normalized' column is True when USE_REAGENT_REFERENCE was applied —
    ratio columns then represent values relative to the reagent control (1.0 = same
    as reagent), not raw pixel ratios.
    """
    fieldnames = [
        "group", "concentration_cfu_ml", "time_min",
        "R", "G", "B",
        "RG_ratio", "RB_ratio", "GB_ratio",
        "reagent_normalized",
    ]
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
                    "R":       round(rgb["R"][i], 4)          if i < len(rgb.get("R", [])) else "",
                    "G":       round(rgb["G"][i], 4)          if i < len(rgb.get("G", [])) else "",
                    "B":       round(rgb["B"][i], 4)          if i < len(rgb.get("B", [])) else "",
                    "RG_ratio": round(rgb["RG_ratio"][i], 6)  if rgb.get("RG_ratio") and i < len(rgb["RG_ratio"]) else "",
                    "RB_ratio": round(rgb["RB_ratio"][i], 6)  if rgb.get("RB_ratio") and i < len(rgb["RB_ratio"]) else "",
                    "GB_ratio": round(rgb["GB_ratio"][i], 6)  if rgb.get("GB_ratio") and i < len(rgb["GB_ratio"]) else "",
                    "reagent_normalized": reagent_normalized,
                })
    print(f"  CSV saved: {path}")


def _save_csv(results: list[dict], path: str):
    """Save a summary CSV with one row per sample."""
    fieldnames = [
        # ── Identity ──────────────────────────────────────────────────────────
        "name", "concentration_cfu_ml", "ground_truth_positive", "overall_positive",
        # ── Combined score (primary ANOVA variable) ───────────────────────────
        "total_detection_score",        # SD event count + slope total weighted score
        "load_tier",
        "concordance_status",           # concordant_negative/positive, timing_discordant, algorithm_discordant
        "concordance_discordance_min",  # |sd_first - slope_first| when both detected
        # ── SD algorithm ──────────────────────────────────────────────────────
        "sd_detected", "sd_first_detection_min", "sd_event_count",
        # ── Slope algorithm ───────────────────────────────────────────────────
        "slope_detected", "slope_first_detection_min", "slope_total_score",
        # ── Likert algorithm ──────────────────────────────────────────────────
        "likert_detected", "likert_first_detection_min", "likert_event_count",
        # ── Continuous kinetic metrics ────────────────────────────────────────
        "rg_auc",                       # area under normalized R/G curve (ratio·min)
        "rb_auc",                       # area under normalized R/B curve
        "gb_auc",                       # area under normalized G/B curve
        "time_to_ratio_threshold",      # first t (min) where R/G ≥ RATIO_THRESHOLD
        "rate_of_change_rg",            # slope of R/G in last SLOPE_CURRENT_POINTS frames (ratio/min)
        # ── Final ratio values ────────────────────────────────────────────────
        "final_RG_ratio", "final_RB_ratio", "final_GB_ratio",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            row = {k: r.get(k, "") for k in fieldnames}
            row["overall_positive"] = (
                r.get("sd_detected") or r.get("slope_detected") or r.get("likert_detected")
            )
            # Round floats for readability
            for key in ("total_detection_score", "slope_total_score",
                        "rg_auc", "rb_auc", "gb_auc", "rate_of_change_rg"):
                if r.get(key) is not None:
                    row[key] = round(r[key], 4)
            if r.get("time_to_ratio_threshold") is not None:
                row["time_to_ratio_threshold"] = int(r["time_to_ratio_threshold"])
            writer.writerow(row)
    print(f"  CSV saved: {path}")


def main():
    parser = argparse.ArgumentParser(description="Smartphone Colorimetric Assay Analyzer")
    parser.add_argument("--data_dir",     default="data/",            help="Path to data folder")
    parser.add_argument("--results_dir",  default=config.RESULTS_DIR, help="Where to write results")
    args = parser.parse_args()

    data_dir    = args.data_dir
    results_dir = args.results_dir
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

    # Print the order so the user can confirm before any pop-ups appear.
    # Reagent reference is loaded but not classified — marked with (reference).
    sample_folders = [f for f in subdirs if not _should_skip(f)]
    print("\n  Processing order:")
    for i, f in enumerate(sample_folders, 1):
        tag = "  [reference — loaded for normalization, not classified]" \
              if _is_reagent_reference(f) else ""
        print(f"    {i}. {f}{tag}")
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

    # ── Stage 1a: Camera artifact detection ───────────────────────────────────
    # A single-frame R/G jump larger than ARTIFACT_JUMP_THRESHOLD almost certainly
    # reflects a camera auto-exposure or white-balance change, not bacterial conversion.
    # Real bacterial conversion is gradual (~0.01–0.05 ratio units per 5-min frame).
    #
    # Artifact timepoints are detected primarily from the reagent reference (most
    # reliable baseline — should stay flat).  Any detected artifact timepoints are
    # stamped onto all sample rgb_data dicts so the SD and slope algorithms can
    # reset their rolling windows and suppress false detections around those frames.
    artifact_threshold = getattr(config, "ARTIFACT_JUMP_THRESHOLD", 0.30)

    # Collect artifact timepoints from the reagent reference if available,
    # otherwise fall back to scanning each sample individually.
    ref_artifact_times: set = set()
    ref_sample_for_artifact = next(
        (s for s in samples if _is_reagent_reference(s["name"])), None
    )
    # If reagent reference was already removed from samples (post-normalization),
    # re-detect from the normalized samples using a tighter per-sample scan.
    scan_source = ref_sample_for_artifact or None

    if scan_source:
        rg_ref = [v for v in scan_source["rgb_data"].get("RG_ratio", []) if v is not None]
        t_ref  = scan_source["rgb_data"].get("times", [])

        # Find all frame-to-frame jumps exceeding the artifact threshold.
        jump_idxs = [i for i in range(1, len(rg_ref))
                     if abs(rg_ref[i] - rg_ref[i - 1]) >= artifact_threshold]

        # Mark the FULL contamination window, not just the endpoints.
        # A spike at index i_start and recovery at i_end leave ALL frames between
        # them elevated — the camera AE hasn't settled yet even though the
        # per-frame deltas look small.  Marking only the spike and recovery frames
        # (abs(delta) >= threshold) lets the intermediate frames contaminate the
        # slope baseline.
        # Strategy: for each consecutive spike→recovery pair, mark every frame from
        # the spike through the recovery (inclusive).  Without a following jump the
        # spike frame itself is still marked.
        for k, j_idx in enumerate(jump_idxs):
            ref_artifact_times.add(t_ref[j_idx] if j_idx < len(t_ref) else None)
            if k + 1 < len(jump_idxs):
                next_j = jump_idxs[k + 1]
                for between_idx in range(j_idx + 1, next_j):
                    t_between = t_ref[between_idx] if between_idx < len(t_ref) else None
                    ref_artifact_times.add(t_between)

    # Also scan each sample for jumps (catches cases where reference is missing)
    per_sample_artifact_times: dict[str, set] = {}
    artifact_warnings = []
    for s in samples:
        rg    = [v for v in s["rgb_data"].get("RG_ratio", []) if v is not None]
        times = s["rgb_data"].get("times", [])

        # Separate positive spikes (+) from negative recoveries (−).
        # A camera AE event is a positive spike followed by a negative recovery.
        # Real bacterial conversion is a large positive jump WITHOUT a following
        # negative recovery (the ratio stays elevated or keeps rising).
        #
        # Window expansion ONLY applies to spike(+) → recovery(−) pairs.
        # A positive jump that has no subsequent negative recovery is left unmarked
        # (it is genuine biology, not a camera event).
        # An unmatched recovery (no preceding spike) marks just that single frame.
        spike_idxs    = [i for i in range(1, len(rg)) if rg[i] - rg[i - 1] >=  artifact_threshold]
        recovery_idxs = [i for i in range(1, len(rg)) if rg[i] - rg[i - 1] <= -artifact_threshold]

        sample_artifacts: set = set()
        matched_recoveries: set = set()

        for spike_idx in spike_idxs:
            # Find the earliest unmatched recovery that follows this spike.
            next_rec = next(
                (r for r in recovery_idxs if r > spike_idx and r not in matched_recoveries),
                None,
            )
            if next_rec is not None:
                # Classic AE spike → recovery: mark the full contamination window.
                t_spike = times[spike_idx] if spike_idx < len(times) else None
                sample_artifacts.add(t_spike)
                artifact_warnings.append((s["name"], t_spike,
                                          abs(rg[spike_idx] - rg[spike_idx - 1])))
                matched_recoveries.add(next_rec)
                t_rec = times[next_rec] if next_rec < len(times) else None
                sample_artifacts.add(t_rec)
                for between_idx in range(spike_idx + 1, next_rec):
                    t_between = times[between_idx] if between_idx < len(times) else None
                    sample_artifacts.add(t_between)
            # If no recovery follows, this is likely real bacterial conversion.
            # Do NOT mark it — marking it would reset the SD rolling window and
            # suppress detection of genuine colour change.

        # Mark any recovery that had no paired spike (e.g. reference caught the spike,
        # but the per-sample delta was below threshold at the spike frame).
        for rec_idx in recovery_idxs:
            if rec_idx not in matched_recoveries:
                t_rec = times[rec_idx] if rec_idx < len(times) else None
                sample_artifacts.add(t_rec)
                artifact_warnings.append((s["name"], t_rec,
                                          abs(rg[rec_idx] - rg[rec_idx - 1])))

        per_sample_artifact_times[s["name"]] = sample_artifacts

    # Combine: reference-detected artifacts apply to everyone (camera-wide event).
    # Sample-specific artifacts apply only to that sample.
    all_detected = ref_artifact_times.union(*per_sample_artifact_times.values()) \
        if per_sample_artifact_times else ref_artifact_times

    for s in samples:
        combined = ref_artifact_times | per_sample_artifact_times.get(s["name"], set())
        s["rgb_data"]["artifact_timepoints"] = combined

    if all_detected:
        unique_warnings = sorted(set((t, d) for _, t, d in artifact_warnings), key=lambda x: x[0] or 0)
        print(
            "\n  !! CAMERA ARTIFACT WARNING — large single-frame R/G jumps detected.\n"
            "     These most likely reflect an auto-exposure or white-balance change,\n"
            "     NOT bacterial conversion.  Verify that AE/AF was locked before capture.\n"
            f"     Artifact timepoints: {sorted(t for t in all_detected if t is not None)}\n"
            "     SD rolling windows will be RESET at these timepoints to prevent\n"
            "     artifact-driven false detections.\n"
        )
        for gname, t, d in artifact_warnings[:10]:   # cap output to 10 lines
            print(f"     {gname:<35}  t={t} min  Δ R/G = {d:.3f}")
        if len(artifact_warnings) > 10:
            print(f"     ... and {len(artifact_warnings) - 10} more")
        print(
            f"\n     Threshold: ARTIFACT_JUMP_THRESHOLD = {artifact_threshold} "
            f"(edit in config.py if needed)\n"
        )
    else:
        print("  No single-frame artifact jumps detected  ✓")

    # ── Stage 1b: Timepoint count validation ──────────────────────────────────
    time_counts = {s["name"]: len(s["rgb_data"]["times"]) for s in samples}
    unique_counts = set(time_counts.values())
    if len(unique_counts) > 1:
        print("\n  !! WARNING: Groups have different numbers of timepoints.")
        print("     Normalization will be misaligned at timepoints missing from any folder.")
        print("     Check that no images are missing or extra:\n")
        for name, count in time_counts.items():
            marker = " <-- DIFFERENT" if count != max(unique_counts) else ""
            print(f"    {count:3d} timepoints  {name}{marker}")
        print()
    else:
        tp_count = next(iter(unique_counts))
        print(f"  All groups: {tp_count} timepoints  ✓")

    # ── Stage 1c: Reagent reference normalization ──────────────────────────────
    if getattr(config, "USE_REAGENT_REFERENCE", False):
        ref_sample = next(
            (s for s in samples if _is_reagent_reference(s["name"])), None
        )
        if ref_sample:
            ref_rgb_raw = ref_sample["rgb_data"]
            ref_csv = os.path.join(results_dir, "reagent_reference_timeseries.csv")
            _save_reagent_reference_csv(ref_rgb_raw, ref_sample["name"], ref_csv)
            samples = [s for s in samples if not _is_reagent_reference(s["name"])]
            _normalize_to_reagent_reference(samples, ref_rgb_raw)
            print(
                f"\n  Reagent reference: '{ref_sample['name']}'\n"
                f"  Ratios normalized for {len(samples)} samples. "
                f"(Raw reference saved to {ref_csv})"
            )
        else:
            print(
                f"\n  WARNING: USE_REAGENT_REFERENCE=True but no folder matching "
                f"'{config.REAGENT_REFERENCE_FOLDER}' was found. "
                f"Skipping normalization."
            )

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
    _save_rgb_timeseries_csv(
        results, rgb_csv_path,
        reagent_normalized=getattr(config, "USE_REAGENT_REFERENCE", False),
    )

    # ── Stage 6: Per-group tables ─────────────────────────────────────────────
    print("\n[+] Per-group algorithm tables:")
    print_and_save_tables(results, results_dir)

    print(f"\nDone. All results in: {results_dir}")


if __name__ == "__main__":
    main()
