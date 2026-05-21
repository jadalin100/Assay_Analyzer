# =============================================================================
# config.py — All tunable parameters in one place.
# Edit this file to adapt the pipeline to your experimental setup.
# =============================================================================

# -----------------------------------------------------------------------------
# Timing
# -----------------------------------------------------------------------------
CAPTURE_INTERVAL_MIN = 5        # Minutes between image captures
TOTAL_DURATION_MIN = 120        # Total assay duration in minutes

# -----------------------------------------------------------------------------
# Image / ROI
# -----------------------------------------------------------------------------
ROI = None                  # pop-up appears on first image — draw the pad box
                            # Set to (x, y, width, height) to skip pop-up.

PER_IMAGE_ROI = False       # True = draw ROI for every image (hand-held photos)

PER_FOLDER_ROI = True       # True = draw ROI once per sample folder (recommended
                            # when all folders contain the same photos and each
                            # group must be boxed individually in the first image)

BLUR_BEFORE_EXTRACTION = True
BLUR_KERNEL_SIZE = 5            # Must be odd

# -----------------------------------------------------------------------------
# Color reference card normalization
# -----------------------------------------------------------------------------
USE_REFERENCE_CARD = True
REFERENCE_CARD_ROI = None       # pop-up appears on first image — draw white swatch
                                # Set to (x, y, width, height) to skip pop-up.

# -----------------------------------------------------------------------------
# Reagent control reference normalization
# -----------------------------------------------------------------------------
# When True, the resazurin reagent control folder is loaded but NOT classified.
# Instead, every other group's R/G, R/B, and G/B ratio values are divided by
# the reagent control's values at the same timepoint.
#
# Effect: a stable (non-converting) group stays at ~1.0 throughout.
#         a colour-changing group rises above 1.0 as bacteria convert resazurin.
#         Auto-exposure drift and lighting changes are cancelled by construction.
#
# NOTE: SD and slope thresholds were calibrated on raw ratios ~1.0.
# After normalization the mean is still ~1.0 for negatives, so threshold
# formulas remain valid in sign — but the absolute noise floor will be lower
# (drift is removed).  Recalibrate SD_THRESHOLD_ALL/TWO after the first
# locked-exposure normalized trial.
USE_REAGENT_REFERENCE = True
# Folder name prefix (case-insensitive) that identifies the reagent control.
# Any folder whose name starts with this string is treated as the reference.
REAGENT_REFERENCE_FOLDER = "resazurin reagent"

# -----------------------------------------------------------------------------
# Colorimetric ratio metrics
# -----------------------------------------------------------------------------
COMPUTE_RG_RATIO = True         # R / G  (recommended primary metric)
COMPUTE_RB_RATIO = True         # R / B  (strong signal: red rises, blue falls simultaneously)
COMPUTE_GB_RATIO = True         # G / B  (complementary; G falls as resazurin converts)

# -----------------------------------------------------------------------------
# Camera artifact detection
# -----------------------------------------------------------------------------
# Real bacterial conversion changes the R/G ratio gradually (~0.01–0.05 per 5-min frame).
# A single-frame delta larger than this threshold is almost certainly a camera
# auto-exposure or white-balance event, not biology.
# The pipeline will print a warning listing the group and timepoint if this fires.
# Set to None to disable the check entirely.
ARTIFACT_JUMP_THRESHOLD = 0.30   # ratio units per 5-min frame

# -----------------------------------------------------------------------------
# Standard Deviation Algorithm (Malanoski et al., 2016)
# -----------------------------------------------------------------------------
# SD_WINDOW = 9 → rolling window of 45 min (9 × 5-min intervals at CAPTURE_INTERVAL_MIN=5).
# Full 120-min trial has 25 timepoints; window slides so recent variance is tracked,
# not all-time CV.
SD_WINDOW = 9

# Thresholds calibrated 2026-05-18 on UNLOCKED auto-exposure data (noisy).
#   Sweep of thr_all × thr_two pairs on current trial data:
#     thr_all=0.060, thr_two=0.080 → FP=[], FN=[10^5, 10^6]   (misses 10^6)
#     thr_all=0.040, thr_two=0.080 → FP=[], FN=[10^5]          ← OPTIMAL
#     thr_all=0.040, thr_two=0.075 → FP=[10^4], FN=[10^5]
#   thr_all=0.040 lets the all-three condition catch 10^6 (cv_RG≈0.040, cv_RB≈0.094 at t=115).
#   thr_two=0.080 eliminates the 10^4 false positive (cv_RB=0.075 no longer triggers).
#   10^5 is still missed by SD — caught only weakly by slope (borderline clinical threshold).
# !! RECALIBRATE after first locked-exposure + reagent-normalized trial !!
#   With locked exposure, negative CV should drop to ~0.018.
#   Lower SD_THRESHOLD_ALL toward 0.025–0.035 and SD_THRESHOLD_TWO toward 0.030–0.045
#   until negatives stop triggering, then verify true positives still detect.
SD_THRESHOLD_ALL = 0.040        # ALL ratio channels must exceed this, OR
SD_THRESHOLD_TWO = 0.080        # ANY two ratio channels must exceed this

# Minimum consecutive timepoints that must meet the SD threshold before detection
# is declared. Mirrors SLOPE_MIN_CONSECUTIVE — eliminates single-frame spikes.
# first_detection_min is set to the START of the run, not the confirming frame.
SD_MIN_CONSECUTIVE = 2

# -----------------------------------------------------------------------------
# Slope Algorithm (Ciaccheri et al., 2023)
# -----------------------------------------------------------------------------
# Calibrated from 10^5 baseline trial (2026-05-16, 5-min intervals, 0–55 min + t=120 final).
# Stable phase confirmed t=0–30 min (indices 0–6). Color change for 10^5 occurs ~60–120 min.
# Non-overlapping windows guaranteed: first score at i=10 (t=50 min) since bp+cp-1=10.
# After next full trial (all groups, 120 min): verify overlap artifact is absent and
# adjust SLOPE_RAPID_CHANGE_FALLBACK down if color change first appears earlier than t=50.
SLOPE_BASELINE_POINTS = 7       # Indices 0–6 = t=0–30 min (confirmed stable phase)
SLOPE_CURRENT_POINTS  = 4       # Rolling window of 4 points = 20-min window
SLOPE_RAPID_CHANGE_FALLBACK = 10 # Index 10 = t=50 min at CAPTURE_INTERVAL_MIN=5.
                                  # Formula: index × CAPTURE_INTERVAL_MIN = time.
                                  # If you change CAPTURE_INTERVAL_MIN, recalculate
                                  # the index so the first score still lands at ~t=50 min.

# Threshold angle formula: θ = A·exp(−mean/B) + C  (degrees)
# Calibrated for RATIO channels (mean ≈ 0.9–1.2 after reagent normalization), NOT raw 0–255.
# At mean=1.0: θ ≈ 0.20·exp(−0.5) + 0.05 ≈ 0.17°
# Signal angle for 10^5 colour change estimated at ~0.30° from baseline trial.
# !! RECALIBRATE after first full 120-min reagent-normalized trial across all concentrations !!
#   - If false positives appear: raise SLOPE_C_SINGLE (higher floor makes less sensitive)
#   - If 10^5 is missed: lower SLOPE_C_SINGLE
#   - Recalibrate SLOPE_A/B after you have angle measurements from multiple concentrations.
SLOPE_A_SINGLE = 0.20
SLOPE_B_SINGLE = 2.0
SLOPE_C_SINGLE = 0.05

SLOPE_A_MULTI = 0.20
SLOPE_B_MULTI = 2.0
SLOPE_C_MULTI = 0.05

SLOPE_R2_SCORE_1 = 0.67
SLOPE_R2_SCORE_2 = 0.80

# With 3 ratio channels (R/G, R/B, G/B) each scoring 0/1/2, max total = 6.
# Threshold > 4 requires strong multi-channel evidence — practically this means
# at least two channels must score 2 each (strong angular divergence + R²>0.80),
# which corresponds to clear, sustained resazurin-to-resorufin conversion.
# Calibrated from first trial (2026-05-18, unlocked AE, photo-converted reference):
#   negatives (10^1–10^4) peaked at 4 → threshold > 4 eliminates all FP
#   10^6/7/8 regularly scored 5–6 → correctly detected
#   10^5 (borderline) scored ≤4 → missed by slope, caught by SD instead
# !! RECALIBRATE after first locked-exposure reagent-normalized trial.
#    With clean data, negatives should score 0 → threshold can drop toward 2–3. !!
SLOPE_WEIGHTED_SCORE_THRESHOLD = 4.0

# Minimum consecutive timepoints exceeding SLOPE_WEIGHTED_SCORE_THRESHOLD
# before detection is declared.  2 frames = 10 min at 5-min capture intervals.
# first_detection_min records the START of the run, not the confirming frame.
SLOPE_MIN_CONSECUTIVE = 2

# -----------------------------------------------------------------------------
# Experimental groups
# -----------------------------------------------------------------------------
CONCENTRATIONS = {
    # Standard keyword prefixes
    "control":          0,
    "reagent_control":  0,
    "positive_control": 1e7,
    "low":              1e4,
    "medium":           1e5,
    "high":             1e6,
    # Numeric folder names (e.g. "10^1" through "10^8")
    "10^1":             1e1,
    "10^2":             1e2,
    "10^3":             1e3,
    "10^4":             1e4,
    "10^5":             1e5,
    "10^6":             1e6,
    "10^7":             1e7,
    "10^8":             1e8,
    # Descriptive negative control names
    "sterile":          0,
    "resazurin":        0,
    "negative":         0,
    "reagent":          0,
}

# Folders whose names contain any of these words are silently skipped
SKIP_FOLDERS = ["reference card", "reference", ".ds_store"]

# Desired processing order for sample folders (case-insensitive prefix match).
# Folders not matching any entry are processed last, in alphabetical order.
FOLDER_ORDER = [
    # Reagent reference — loaded first, used for normalization, not classified
    "resazurin reagent control",
    "resazurin reagent",
    "reagent control",
    # Sterile negative controls
    "sterile negative control",
    "sterile negative",
    "sterile",
    "negative",
    # Sub-threshold concentrations (clinical negatives; ground_truth_positive = False)
    # Note: 10^4 can grow ~6 doublings over 120 min and may show a late weak signal.
    "10^1",
    "10^2",
    "10^3",
    "10^4",
    "low",
    # UTI-positive threshold and above
    "10^5",
    "medium",
    "10^6",
    "high",
    "10^7",
    "positive control",
    "10^8",
]

UTI_POSITIVE_THRESHOLD_CFU = 1e5   # ≥ 10^5 CFU/mL = true positive (clinical UTI threshold)

# -----------------------------------------------------------------------------
# Bacterial load tiers — based on time-to-first-detection (minutes)
# Detection earlier = faster conversion = higher bacterial load.
# Cutoffs are PLACEHOLDERS — calibrate after first full 120-min trial.
# Use the earliest first_detection_min across SD and slope for assignment.
# -----------------------------------------------------------------------------
LOAD_TIERS = {
    # tier_label:       max first_detection_min to qualify for this tier
    # !! ALL CUTOFFS ARE PLACEHOLDERS except 120 min !!
    # Calibrate 60-min and 90-min boundaries after first full 120-min
    # reagent-normalized trial by reading first_detection_min for 10^6 and 10^7 groups.
    "heavy  (≥10^7)":   60,   # Placeholder — calibrate from 10^7 first_detection_min
    "moderate (≥10^6)": 90,   # Placeholder — calibrate from 10^6 first_detection_min
    "low  (≥10^5)":    120,   # Confirmed — 10^5 converts between t=60–120 min
    # Anything detected after 120 min or not at all → negative
}

# -----------------------------------------------------------------------------
# Likert Colour-Matching Algorithm
# -----------------------------------------------------------------------------
# Five reference colours derived from real trial photos (trials 1–4).
# Level 1 = unreacted resazurin (blue-purple, negative control at t=120 min).
# Level 5 = fully converted resorufin (bright magenta, 10^8 CFU/mL at t=120 min,
#           averaged across all four trials).
# Levels 2–4 interpolate the transition.
#
# These values were sampled directly from JPEG well photos using PIL, so they
# match the raw camera pixel scale (0-255).  The pipeline's white-balance
# correction applies multiplicative scale factors (~1.0–1.1) per channel —
# a small offset that does not affect level assignment given the large colour
# distance (~180 units) between Level 1 and Level 5.
#
# !! RECALIBRATE after the first locked-exposure trial !!
#    From rgb_timeseries.csv, read the mean R, G, B of the sterile negative
#    at t=120 min → Level 1, and of the 10^8 group at t=120 min → Level 5.
#    Linearly interpolate Levels 2–4 between those anchors (or sample
#    directly from intermediate concentration wells).
LIKERT_REFERENCE_COLORS = [
    ( 94,  86, 198),   # Level 1 — No change        (Negative)
    ( 80,  68, 175),   # Level 2 — Slight change    (Likely negative)
    ( 92,  30, 154),   # Level 3 — Moderate change  (Ambiguous)
    (210,  50, 158),   # Level 4 — Clear change     (Likely positive)
    (240,  60, 135),   # Level 5 — Definite change  (Positive) — bright hot pink (ref: resorufin fully converted)
]

# Minimum Likert level to count as a "positive event" at a timepoint.
# Score ≥ 4 = test positive (consistent with the printed Likert reference card).
LIKERT_POSITIVE_THRESHOLD = 4

# Minimum consecutive timepoints at >= LIKERT_POSITIVE_THRESHOLD before
# detection is declared.  Mirrors SD_MIN_CONSECUTIVE and SLOPE_MIN_CONSECUTIVE.
LIKERT_MIN_CONSECUTIVE = 2

# -----------------------------------------------------------------------------
# Statistical analysis
# -----------------------------------------------------------------------------
ALPHA = 0.05

# Fixed normalized-ratio threshold for time-to-threshold metric.
# After reagent-reference normalization, a stable negative stays at ~1.0.
# A value of 1.10 means the group's R/G ratio is 10% above the reagent control —
# a conservative threshold that should not fire on noise.
# Calibrate after first normalized trial: raise if negatives briefly touch 1.10,
# lower if true positives are slow to reach it.
RATIO_THRESHOLD = 1.50   # raised from 1.10 — artifact inflated baselines to ~1.4 for
                         # negatives in the unlocked-AE trial; re-evaluate after first
                         # locked-exposure trial where negatives should stay near 1.0

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------
RESULTS_DIR = "results"
PLOT_DPI = 150
