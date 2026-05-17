# Smartphone-Based Resazurin UTI Detection Pipeline

A Python pipeline for quantitative detection of urinary tract infections (UTIs) using smartphone colorimetric analysis of resazurin-based assays. The system captures the resazurin-to-resorufin color transition over a 120-minute incubation window and applies two independent detection algorithms to classify bacterial load — no specialized laboratory equipment required.

## Overview

Resazurin (Alamar Blue) shifts from blue to pink as metabolically active bacteria reduce it to resorufin. This pipeline extracts RGB pixel values from time-lapse smartphone images of assay wells, normalizes them against a co-incubated reagent control to cancel lighting drift, and runs two complementary algorithms to detect and quantify the conversion signal.

**Clinical threshold:** ≥ 10⁵ CFU/mL (standard UTI diagnostic criterion)

## Features

- **Dual-algorithm detection** — SD (variance-based) and Slope (angular divergence) algorithms operate independently and are combined into a unified Total Detection Score
- **Reagent reference normalization** — divides each sample's ratios by the reagent-only control at each timepoint, canceling auto-exposure drift and lighting changes
- **White-balance correction** — per-image multiplicative correction using a physical reference card
- **Three colorimetric channels** — R/G, R/B, and G/B ratios computed and analyzed in parallel
- **Consecutive-frame requirement** — both algorithms require ≥ 2 consecutive frames above threshold before firing, eliminating single-frame noise spikes
- **Quantitative kinetic metrics** — AUC (per channel), Time to Threshold, and Rate of Change for regression and ANOVA analysis
- **Bacterial load tiering** — Heavy / Moderate / Low burden classification based on time-to-first-detection
- **Concordance assessment** — agreement status between the two algorithms, used as a confidence qualifier and clinical tiebreaker
- **Statistical analysis** — Shapiro-Wilk normality test, Kruskal-Wallis / one-way ANOVA, Tukey HSD / Dunn's Bonferroni post-hoc
- **Power analysis** — Wilson score CI-based sample size estimation for sensitivity and specificity targets
- **Structured CSV exports** — long-format time-series tables, per-group summary sheets, total detection score breakdown
- **Visualization** — ratio time-series plots, SD CV plots, slope score heatmaps, stacked detection score bar charts

## Project Structure

```
assay_analyzer/
├── main.py               # Pipeline entry point — orchestrates all stages
├── config.py             # All tunable parameters (edit this before each run)
├── ingestion.py          # Image loading, ROI extraction, RGB ratio computation
├── classifier.py         # Per-sample classification, kinetic metrics, concordance
├── sd_algorithm.py       # Standard Deviation algorithm (Malanoski et al., 2016)
├── slope_algorithm.py    # Slope/angle algorithm (Ciaccheri et al., 2023)
├── reporting.py          # CSV export functions
├── stats.py              # Statistical testing (ANOVA, Kruskal-Wallis, post-hoc)
├── visualization.py      # Matplotlib plots
├── power_analysis.py     # Sample size / CI estimation (run standalone)
├── metrics.py            # Shared metric utilities
├── requirements.txt      # Python dependencies
└── Assay_Methodology.docx  # Full methodology document
```

## Setup

### Requirements

- Python 3.10+
- pip packages:

```bash
pip install -r requirements.txt
```

Dependencies: `numpy`, `scipy`, `matplotlib`, `opencv-python`, `statsmodels`

Optional (for post-hoc tests):
```bash
pip install pandas scikit-posthocs
```

### Image Folder Structure

Organize your data directory as follows, with one subfolder per sample group:

```
data/
├── resazurin reagent/     # Reagent-only control — used for normalization, not classified
├── sterile negative/      # True negative (0 CFU/mL)
├── 10^1/                  # Sub-threshold negatives
├── 10^2/
├── 10^3/
├── 10^4/
├── 10^5/                  # Clinical UTI threshold — first true positive
├── 10^6/
├── 10^7/
└── reference card/        # White-balance reference card (skipped from classification)
```

Folder names are matched by case-insensitive prefix against the entries in `config.py`. Any folder not matching a known prefix is processed last.

## Usage

### Running the full pipeline

```bash
python3 main.py
```

On the first run you will be prompted to draw ROI boxes around each well using your mouse. Coordinates are stored and reused for all subsequent images in the series.

### Running the power analysis standalone

```bash
python3 power_analysis.py
# Override defaults from the command line:
python3 power_analysis.py --sensitivity 0.85 --specificity 1.0 --target 0.15
```

## Configuration

All parameters are in `config.py`. Key settings to review before each run:

| Parameter | Default | Description |
|---|---|---|
| `CAPTURE_INTERVAL_MIN` | `5` | Minutes between image captures |
| `TOTAL_DURATION_MIN` | `120` | Total assay duration |
| `USE_REAGENT_REFERENCE` | `True` | Normalize against reagent control well |
| `USE_REFERENCE_CARD` | `True` | Apply white-balance correction via reference card |
| `PER_FOLDER_ROI` | `True` | Draw one ROI per sample folder (recommended) |
| `SD_THRESHOLD_ALL` | `0.060` | CV threshold — all three channels must exceed |
| `SD_THRESHOLD_TWO` | `0.075` | CV threshold — any two channels must exceed |
| `SD_MIN_CONSECUTIVE` | `2` | Consecutive frames required before SD detection fires |
| `SLOPE_MIN_CONSECUTIVE` | `2` | Consecutive frames required before slope detection fires |
| `RATIO_THRESHOLD` | `1.10` | Fixed R/G threshold for Time-to-Threshold metric |
| `UTI_POSITIVE_THRESHOLD_CFU` | `1e5` | Clinical positive cutoff (CFU/mL) |

> **Note:** SD thresholds and slope parameters (`SLOPE_A/B/C_SINGLE`) are currently calibrated for unlocked auto-exposure data. After the first locked-exposure reagent-normalized trial, lower the SD thresholds and recalibrate the slope parameters using the measured angle values per concentration group.

## Output

Results are saved to the `results/` directory:

| File | Contents |
|---|---|
| `rgb_timeseries.csv` | Normalized R/G, R/B, G/B ratios — one row per group × timepoint |
| `sd_results.csv` | CV values at every timepoint for all channels + detection summary |
| `slope_results.csv` | Per-timepoint slope scores + detection summary |
| `detection_events.csv` | First detection time, load tier, concordance — one row per group |
| `total_detection_score.csv` | Per-timepoint score breakdown + full kinetic metrics summary |
| `classification_results.csv` | Complete per-group result with all metrics |
| `anova_results.csv` | Statistical test results across concentration groups |
| `*.png` | Ratio plots, SD CV plots, slope heatmaps, detection score charts |

## Detection Algorithms

### SD Algorithm (Malanoski et al., 2016)
Monitors the rolling coefficient of variation (CV = SD / mean) over a 9-timepoint window. Detection fires when all three ratio channels exceed `SD_THRESHOLD_ALL`, or any two exceed `SD_THRESHOLD_TWO`, for ≥ 2 consecutive frames. A 30-minute guard period suppresses early-assay noise.

### Slope Algorithm (Ciaccheri et al., 2023)
Computes the cosine of the angle between a fixed 7-point baseline regression (t = 0–30 min) and a sliding 4-point current window. An adaptive threshold `θ = A·exp(−mean/B) + C` accounts for the mean ratio value. Each channel is scored 0–2 based on angular divergence and regression quality (R²). Detection fires when the sum exceeds a weighted threshold for ≥ 2 consecutive frames. Scoring begins at t = 50 min to ensure non-overlapping baseline and current windows.

### Total Detection Score
`Total Detection Score = SD Event Count + Slope Total Weighted Score`

A continuous, ANOVA-ready variable that captures sustained evidence of color conversion from two complementary perspectives.

## References

- Malanoski, A.P. et al. (2016). Simultaneous Identification of Multiple Analytes Using a Self-Reporting Spectroscopic Hybridization Assay. *Analytical Chemistry.*
- Ciaccheri, L. et al. (2023). A Machine Learning Approach for the Rapid Identification of Bacterial Contamination. *Sensors.*
- Ramirez, R. et al. Resazurin-based methods for evaluating bacterial viability. *Journal of Microbiological Methods.*
