# =============================================================================
# config.py — All tunable parameters in one place.
# Edit this file to adapt the pipeline to your experimental setup.
# =============================================================================

# -----------------------------------------------------------------------------
# Timing
# -----------------------------------------------------------------------------
CAPTURE_INTERVAL_MIN = 10       # Minutes between image captures
TOTAL_DURATION_MIN = 90         # Total assay duration in minutes

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
# Colorimetric ratio metrics
# -----------------------------------------------------------------------------
COMPUTE_RB_RATIO = True         # R / B
COMPUTE_RG_RATIO = True         # R / G  (recommended primary metric)

# -----------------------------------------------------------------------------
# Standard Deviation Algorithm (Malanoski et al., 2016)
# -----------------------------------------------------------------------------
# SD_WINDOW = 9 matches the number of timepoints. The rolling window grows
# from 1 point up to 9, so at the final timepoint it uses the full series.
SD_WINDOW = 9
# PLACEHOLDER values calibrated from imaging noise floor (2026-05-12, first trial).
# Original paper value (0.00015) is too low for smartphone JPEG — fires on pure noise.
# After next successful trial: raise until true negatives no longer detect,
# then verify true positives still detect. Update both values together.
SD_THRESHOLD_ALL = 0.025        # ALL channels in a set must exceed this, OR
SD_THRESHOLD_TWO = 0.08         # ANY two channels in a set must exceed this

# -----------------------------------------------------------------------------
# Slope Algorithm (Ciaccheri et al., 2023)
# -----------------------------------------------------------------------------
# PLACEHOLDER values calibrated for 9 timepoints (0–80 min, 10-min intervals).
# TODO: After your next successful trial, note when color change first appears
#       visually and set SLOPE_BASELINE_POINTS = that timepoint index,
#       SLOPE_CURRENT_POINTS = 3-4, SLOPE_RAPID_CHANGE_FALLBACK = 2 or 3.
SLOPE_BASELINE_POINTS = 4       # First 4 timepoints (0–30 min) used as baseline
SLOPE_CURRENT_POINTS  = 4       # Rolling window of 4 recent timepoints
SLOPE_RAPID_CHANGE_FALLBACK = 3 # Skip first 3 indices; first scoring at i=3 (t=30 min)

SLOPE_A_SINGLE = 20
SLOPE_B_SINGLE = 130
SLOPE_C_SINGLE = 0.45

SLOPE_A_MULTI = 70
SLOPE_B_MULTI = 30
SLOPE_C_MULTI = 0.45

SLOPE_R2_SCORE_1 = 0.67
SLOPE_R2_SCORE_2 = 0.80

# Per spec: total score across all channels must be > 1 (i.e., ≥ 2 in integer terms).
# With 6 channels (R, G, B, R/G, R/B, G/B) each scoring 0/1/2, max total = 12.
SLOPE_WEIGHTED_SCORE_THRESHOLD = 1.0

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
    # Numeric folder names (e.g. "10^5", "10^6", "10^7", "10^8")
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
    "resazurin reagent control",
    "resazurin reagent",
    "reagent control",
    "sterile negative control",
    "sterile negative",
    "sterile",
    "negative",
    "low",
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
# Statistical analysis
# -----------------------------------------------------------------------------
ALPHA = 0.05

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------
RESULTS_DIR = "results"
PLOT_DPI = 150
