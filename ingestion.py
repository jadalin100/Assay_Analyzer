# =============================================================================
# ingestion.py — Module 1: Load images, extract ROI, produce RGB time series.
# =============================================================================

import os
import math
import cv2
import numpy as np
import config


def _select_roi_interactively(
    image: np.ndarray,
    window_title: str = "Select ROI — press ENTER to confirm",
    prompt: str = "Draw a rectangle around the PAD strip area, then press ENTER or SPACE.",
) -> tuple[int, int, int, int]:
    """
    Open a window and let the user draw a ROI rectangle. Returns (x, y, w, h).
    After selection, prints the pixel coordinates and mean RGB so the user can
    verify they selected the correct region before the program continues.
    """
    print(f"\n{'─' * 60}")
    print(f"  POP-UP ACTION: {prompt}")
    print(f"{'─' * 60}")
    roi = cv2.selectROI(window_title, image, fromCenter=False)
    cv2.destroyAllWindows()
    if roi == (0, 0, 0, 0):
        raise ValueError("No ROI selected.")

    x, y, w, h = roi
    crop = image[y: y + h, x: x + w]
    mean_bgr = cv2.mean(crop)[:3]
    r_mean = mean_bgr[2]
    g_mean = mean_bgr[1]
    b_mean = mean_bgr[0]

    print(f"  Selected region: x={x}  y={y}  w={w}  h={h}  ({w * h} px)")
    print(f"  Mean color:      R={r_mean:.0f}  G={g_mean:.0f}  B={b_mean:.0f}")

    # Warn if the well ROI looks unexpectedly bright (likely the reference card)
    if r_mean > 210 and g_mean > 210 and b_mean > 210:
        print(f"  !! WARNING: This region is very bright (near-white).")
        print(f"             If you meant to draw around a colored WELL,")
        print(f"             re-run the program and select the well area instead.")

    return roi


def _extract_mean_rgb(
    image: np.ndarray,
    roi: tuple[int, int, int, int],
    label: str = "",
    verbose: bool = True,
) -> tuple[float, float, float]:
    """
    Crop to ROI, optionally blur, return mean (R, G, B).

    Every pixel inside the selected rectangle is included in the average —
    heterogeneous colour across a strip is handled correctly. Gaussian blur
    smooths sharp local noise before the mean is computed.

    If verbose is True, prints per-channel standard deviation across the ROI
    so you can confirm colour variation is being captured and averaged.
    """
    x, y, w, h = roi
    crop = image[y: y + h, x: x + w]

    if config.BLUR_BEFORE_EXTRACTION:
        k = config.BLUR_KERNEL_SIZE
        crop = cv2.GaussianBlur(crop, (k, k), 0)

    mean_bgr = cv2.mean(crop)[:3]
    r, g, b = mean_bgr[2], mean_bgr[1], mean_bgr[0]

    if verbose and label:
        crop_rgb = crop[:, :, ::-1].astype(np.float32)
        std_r = float(np.std(crop_rgb[:, :, 0]))
        std_g = float(np.std(crop_rgb[:, :, 1]))
        std_b = float(np.std(crop_rgb[:, :, 2]))
        n_px = crop.shape[0] * crop.shape[1]
        print(
            f"    [{label}] mean RGB = ({r:.1f}, {g:.1f}, {b:.1f})  "
            f"std = ({std_r:.1f}, {std_g:.1f}, {std_b:.1f})  "
            f"over {n_px} pixels"
        )

    return float(r), float(g), float(b)


def _normalize_by_reference_card(
    r: float, g: float, b: float,
    image: np.ndarray,
    ref_roi: tuple[int, int, int, int],
    verbose: bool = True,
) -> tuple[float, float, float]:
    """
    Normalize PAD RGB values using the white swatch of the reference card.

    Correction is multiplicative (white-balance style):
        scale_channel = 255 / white_swatch_channel
        corrected_pad = raw_pad x scale_channel

    Example: white swatch reads R=240, G=235, B=255 -> scales (1.063, 1.085, 1.000)
    applied to the pad values, removing the lighting bias.
    """
    ref_r, ref_g, ref_b = _extract_mean_rgb(image, ref_roi, label="white swatch", verbose=verbose)

    eps = 1e-6
    scale_r = 255.0 / (ref_r + eps)
    scale_g = 255.0 / (ref_g + eps)
    scale_b = 255.0 / (ref_b + eps)

    norm_r = min(r * scale_r, 255.0)
    norm_g = min(g * scale_g, 255.0)
    norm_b = min(b * scale_b, 255.0)

    if verbose:
        print(
            f"    [ref card] scale factors: R x{scale_r:.4f}  G x{scale_g:.4f}  B x{scale_b:.4f}"
        )
        print(
            f"    [ref card] pad corrected: ({r:.1f}->{norm_r:.1f}, "
            f"{g:.1f}->{norm_g:.1f}, {b:.1f}->{norm_b:.1f})"
        )

    return norm_r, norm_g, norm_b


def load_sample_images(image_dir: str) -> list[np.ndarray]:
    """Load all images from a directory in alphabetical order."""
    extensions = {".jpg", ".jpeg", ".png", ".bmp", ".tif", ".tiff"}
    files = sorted(
        f for f in os.listdir(image_dir)
        if os.path.splitext(f)[1].lower() in extensions
    )
    if not files:
        raise FileNotFoundError(f"No images found in: {image_dir}")
    images = [cv2.imread(os.path.join(image_dir, f)) for f in files]
    print(f"  Loaded {len(images)} images from '{image_dir}'")
    return images


def extract_rgb_series(
    images: list[np.ndarray],
    roi: tuple[int, int, int, int] | None = None,
    folder_name: str = "",
    ref_roi: tuple[int, int, int, int] | None = None,
) -> dict:
    """Extract mean RGB per frame from a list of images."""
    n = len(images)
    interval = config.CAPTURE_INTERVAL_MIN
    per_image = getattr(config, "PER_IMAGE_ROI", False)

    group = folder_name if folder_name else "this group"

    fixed_roi = None
    fixed_ref_roi = None

    if not per_image:
        if roi is not None:
            fixed_roi = roi
        elif config.ROI is not None:
            fixed_roi = config.ROI
        else:
            n_steps = 2 if (config.USE_REFERENCE_CARD and ref_roi is None
                            and config.REFERENCE_CARD_ROI is None) else 1
            print(f"\n  >>> STEP 1 of {n_steps}: Draw around the COLORED WELL for '{group}'")
            print(f"      (Do NOT draw around the reference card — that comes next)")
            fixed_roi = _select_roi_interactively(
                images[0],
                window_title=f"STEP 1/{n_steps} — WELL: {group}",
                prompt=f"Draw rectangle around the COLORED WELL for '{group}' — NOT the reference card.",
            )

        if config.USE_REFERENCE_CARD:
            if ref_roi is not None:
                fixed_ref_roi = ref_roi
                print(f"    [ref card] Reusing reference card ROI from first group: {fixed_ref_roi}")
            elif config.REFERENCE_CARD_ROI is not None:
                fixed_ref_roi = config.REFERENCE_CARD_ROI
            else:
                print(f"\n  >>> STEP 2 of 2: Draw around the WHITE REFERENCE CARD")
                print(f"      (This is the white swatch on your color reference card — NOT a well)")
                fixed_ref_roi = _select_roi_interactively(
                    images[0],
                    window_title="STEP 2/2 — WHITE REFERENCE CARD (not a well)",
                    prompt="Draw rectangle around the WHITE REFERENCE CARD swatch — NOT a well.",
                )
                print(f"    Reference card ROI: {fixed_ref_roi}  "
                      f"(copy into config.REFERENCE_CARD_ROI to skip next time)")

    times, R_vals, G_vals, B_vals = [], [], [], []
    RB_ratio_vals, RG_ratio_vals, GB_ratio_vals = [], [], []

    for i, img in enumerate(images):
        t = i * interval

        if per_image:
            print(f"\n  --- Image {i + 1} of {n}  (t = {t:.0f} min) ---")
            current_roi = _select_roi_interactively(
                img,
                window_title=f"[WELL ROI] Image {i + 1}/{n} — draw PAD ROI, press ENTER",
                prompt=f"WELL ROI for '{group}': draw rectangle around the assay pad (image {i + 1}/{n}, t={t:.0f} min).",
            )
            print(f"    Pad ROI: {current_roi}")
            current_ref_roi = None
            if config.USE_REFERENCE_CARD:
                current_ref_roi = _select_roi_interactively(
                    img,
                    window_title=f"[WHITE REFERENCE CARD] Image {i + 1}/{n} — press ENTER",
                    prompt=f"WHITE REFERENCE CARD: draw rectangle around the white swatch (image {i + 1}/{n}).",
                )
                print(f"    Reference ROI: {current_ref_roi}")
        else:
            current_roi = fixed_roi
            current_ref_roi = fixed_ref_roi

        print(f"\n  t = {t:.0f} min")
        r, g, b = _extract_mean_rgb(img, current_roi, label="assay pad", verbose=True)

        if config.USE_REFERENCE_CARD and current_ref_roi is not None:
            r, g, b = _normalize_by_reference_card(r, g, b, img, current_ref_roi, verbose=True)

        eps = 1e-6
        rb_ratio = r / (b + eps) if config.COMPUTE_RB_RATIO else None
        rg_ratio = r / (g + eps) if config.COMPUTE_RG_RATIO else None
        gb_ratio = g / (b + eps)

        times.append(t)
        R_vals.append(r)
        G_vals.append(g)
        B_vals.append(b)
        RB_ratio_vals.append(rb_ratio)
        RG_ratio_vals.append(rg_ratio)
        GB_ratio_vals.append(gb_ratio)

    return {
        "times": times,
        "R": R_vals,
        "G": G_vals,
        "B": B_vals,
        "RB_ratio": RB_ratio_vals,
        "RG_ratio": RG_ratio_vals,
        "GB_ratio": GB_ratio_vals,
        "roi": "per-image" if per_image else fixed_roi,
        "ref_roi": fixed_ref_roi,
        "normalized": config.USE_REFERENCE_CARD,
    }


def extract_rgb_from_video(video_path: str, roi: tuple[int, int, int, int] | None = None, folder_name: str = "", ref_roi: tuple[int, int, int, int] | None = None) -> dict:
    """Extract RGB time series from a video file."""
    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        raise IOError(f"Cannot open video: {video_path}")

    fps = cap.get(cv2.CAP_PROP_FPS)
    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    duration_min = total_frames / fps / 60
    print(f"  Video: {fps:.1f} fps, {total_frames} frames, ~{duration_min:.1f} min")

    interval_min = config.CAPTURE_INTERVAL_MIN
    n_points = math.floor(duration_min / interval_min) + 1
    sample_times = [i * interval_min for i in range(n_points)]

    ret, first_frame = cap.read()
    if not ret:
        raise IOError("Cannot read first frame.")

    group = folder_name if folder_name else "this group"

    per_image = getattr(config, "PER_IMAGE_ROI", False)
    fixed_roi = None
    fixed_ref_roi = None

    if not per_image:
        if roi is not None:
            fixed_roi = roi
        elif config.ROI is not None:
            fixed_roi = config.ROI
        else:
            fixed_roi = _select_roi_interactively(
                first_frame,
                window_title=f"[WELL ROI] {group} — draw rectangle, press ENTER",
                prompt=f"WELL ROI for '{group}': draw a rectangle around this group's well, then press ENTER.",
            )
            print(f"    ROI selected: {fixed_roi}  (copy this into config.ROI to skip next time)")

        if config.USE_REFERENCE_CARD:
            if ref_roi is not None:
                fixed_ref_roi = ref_roi
                print(f"    [ref card] Reusing reference card ROI from first group: {fixed_ref_roi}")
            elif config.REFERENCE_CARD_ROI is not None:
                fixed_ref_roi = config.REFERENCE_CARD_ROI
            else:
                print("  Now select the WHITE swatch on your reference card (asked only once).")
                fixed_ref_roi = _select_roi_interactively(
                    first_frame,
                    window_title="[WHITE REFERENCE CARD] — draw rectangle, press ENTER",
                    prompt="WHITE REFERENCE CARD: draw a rectangle around the white swatch, then press ENTER.",
                )
                print(f"    Reference card ROI: {fixed_ref_roi}  "
                      f"(copy into config.REFERENCE_CARD_ROI to skip next time)")

    times, R_vals, G_vals, B_vals = [], [], [], []
    RB_ratio_vals, RG_ratio_vals, GB_ratio_vals = [], [], []
    n_points = len(sample_times)

    for idx, t_min in enumerate(sample_times):
        frame_idx = int(t_min * 60 * fps)
        cap.set(cv2.CAP_PROP_POS_FRAMES, frame_idx)
        ret, frame = cap.read()
        if not ret:
            print(f"  Warning: could not read frame at t={t_min} min, skipping.")
            continue

        if per_image:
            print(f"\n  --- Frame {idx + 1} of {n_points}  (t = {t_min:.0f} min) ---")
            current_roi = _select_roi_interactively(
                frame,
                window_title=f"[WELL ROI] '{group}' Frame {idx + 1}/{n_points} — press ENTER",
                prompt=f"WELL ROI for '{group}': draw rectangle around the assay pad (t={t_min:.0f} min).",
            )
            current_ref_roi = None
            if config.USE_REFERENCE_CARD:
                current_ref_roi = _select_roi_interactively(
                    frame,
                    window_title=f"[WHITE REFERENCE CARD] Frame {idx + 1}/{n_points} — press ENTER",
                    prompt=f"WHITE REFERENCE CARD: draw rectangle around the white swatch (t={t_min:.0f} min).",
                )
        else:
            current_roi = fixed_roi
            current_ref_roi = fixed_ref_roi

        print(f"\n  t = {t_min:.0f} min")
        r, g, b = _extract_mean_rgb(frame, current_roi, label="assay pad", verbose=True)

        if config.USE_REFERENCE_CARD and current_ref_roi is not None:
            r, g, b = _normalize_by_reference_card(r, g, b, frame, current_ref_roi, verbose=True)

        eps = 1e-6
        rb_ratio = r / (b + eps) if config.COMPUTE_RB_RATIO else None
        rg_ratio = r / (g + eps) if config.COMPUTE_RG_RATIO else None
        gb_ratio = g / (b + eps)

        times.append(t_min)
        R_vals.append(r)
        G_vals.append(g)
        B_vals.append(b)
        RB_ratio_vals.append(rb_ratio)
        RG_ratio_vals.append(rg_ratio)
        GB_ratio_vals.append(gb_ratio)

    cap.release()
    print(f"  Extracted {len(times)} time points from video.")
    return {
        "times": times,
        "R": R_vals,
        "G": G_vals,
        "B": B_vals,
        "RB_ratio": RB_ratio_vals,
        "RG_ratio": RG_ratio_vals,
        "GB_ratio": GB_ratio_vals,
        "roi": "per-image" if per_image else fixed_roi,
        "ref_roi": fixed_ref_roi,
        "normalized": config.USE_REFERENCE_CARD,
    }
