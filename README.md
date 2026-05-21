# Smartphone-Based Resazurin UTI Detection Pipeline

A Python pipeline for quantitative detection of urinary tract infections (UTIs) using smartphone colorimetric analysis of resazurin-based assays. The system captures the resazurin-to-resorufin color transition over a 120-minute incubation window and applies three independent detection algorithms to classify bacterial load — no specialized laboratory equipment required.

## Overview

Resazurin (Alamar Blue) shifts from blue to pink as metabolically active bacteria reduce it to resorufin. This pipeline extracts RGB pixel values from time-lapse smartphone images of assay wells, normalizes them against a co-incubated reagent control to cancel lighting drift, and runs three complementary algorithms to detect and quantify the conversion signal.

**Clinical threshold:** ≥ 10⁵ CFU/mL (standard UTI diagnostic criterion)

## Experimental Data

Four trials have been completed. Each trial folder contains time-lapse photos (t = 0–120 min, 5-min intervals) and full processed results.

| Trial | Groups tested | Notes |
|---|---|---|
| `trial-1/` | 10¹–10⁸, sterile negative | Auto-exposure unlocked |
| `trial-2/` | 10¹–10⁸, sterile negative | Auto-exposure unlocked |
| `trial-3/` | 10¹–10⁸, sterile negative | Auto-exposure unlocked |
| `trial-4/` | 10¹–10⁸, sterile negative | Auto-exposure unlocked |

> All trials to date use unlocked auto-exposure. Algorithm thresholds (SD and slope) are calibrated for this condition. A locked-exposure trial is planned; parameters should be recalibrated afterward.

## Features

- **Triple-algorithm detection** — SD (variance-based), Slope (angular divergence), and Likert (RGB colour-matching) algorithms operate independently; a Combined positive fires when any one detects
- **Reagent reference normalization** — divides each sample's ratios by the reagent-only control at each timepoint, canceling auto-exposure drift and lighting changes
- **White-balance correction** — per-image multiplicative correction using a physical reference card
- **Three colorimetric channels** — R/G, R/B, and G/B ratios computed and analyzed in parallel
- **Consecutive-frame requirement** — all algorithms require ≥ 2 consecutive frames above threshold before firing, eliminating single-frame noise spikes
- **Quantitative kinetic metrics** — AUC (per channel), Time to Threshold, and Rate of Change for regression and ANOVA analysis
- **Bacterial load tiering** — Heavy / Moderate / Low burden classification based on time-to-first-detection
- **Concordance assessment** — agreement status between algorithms, used as a confidence qualifier and clinical tiebreaker
- **Per-trial statistics** — Shapiro-Wilk normality test, Kruskal-Wallis / one-way ANOVA, Tukey HSD / Dunn's Bonferroni post-hoc, saved per trial as `statistics.json`
- **Cross-trial statistics** — pools all 4 trials (N=4 per concentration) for group-level ANOVA/Kruskal-Wallis; results in `cross_trial_results/`
- **Power analysis** — Wilson score CI-based sample size estimation for sensitivity and specificity targets
- **Structured CSV exports** — long-format time-series tables, per-group summary sheets, total detection score breakdown
- **Visualization** — ratio time-series plots, SD CV plots, slope score heatmaps, Likert score step plots, stacked detection score bar charts, performance summary with confusion matrix
- **Likert color reference scale** — `likert_scale.png`, a 5-level colorimetric reference card derived from real well photos for human observer scoring

## Project Structure

```
assay_analyzer/
├── main.py                  # Pipeline entry point — orchestrates all stages
├── config.py                # All tunable parameters (edit this before each run)
├── ingestion.py             # Image loading, ROI extraction, RGB ratio computation
├── classifier.py            # Per-sample classification, kinetic metrics, concordance
├── sd_algorithm.py          # Standard Deviation algorithm (Malanoski et al., 2016)
├── slope_algorithm.py       # Slope/angle algorithm (Ciaccheri et al., 2023)
├── likert_algorithm.py      # Likert colour-matching algorithm (Euclidean RGB distance)
├── reporting.py             # CSV export functions
├── stats.py                 # Per-trial statistical testing (ANOVA, Kruskal-Wallis, post-hoc)
├── cross_trial_stats.py     # Cross-trial pooled statistical analysis
├── visualization.py         # Matplotlib plots
├── metrics.py               # Sensitivity, specificity, PPV, NPV per algorithm
├── power_analysis.py        # Sample size / CI estimation (run standalone)
├── run_all_trials.py        # Batch runner — processes all trial folders sequentially
├── requirements.txt         # Python dependencies
├── likert_scale.png         # Human observer color reference card (Levels 1–5)
├── Assay_Methodology.docx   # Full methodology document
├── cross_trial_results/
│   ├── cross_trial_statistics.json  # Pooled ANOVA/KW test results across all trials
│   └── cross_trial_summary.csv      # Per-group means ± SD for all key metrics
├── trial-1/
│   ├── photos/              # 25 time-lapse JPGs (t = 0–120 min)
│   └── results/             # CSVs, PNGs, statistics.json, groups/
├── trial-2/
│   ├── photos/
│   └── results/
├── trial-3/
│   ├── photos/
│   └── results/
└── trial-4/
    ├── photos/
    └── results/
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

### Running the cross-trial statistical analysis

After all trials have been processed and their results saved to `trial-N/results/`:

```bash
python3 cross_trial_stats.py                  # pool trials 1–4
python3 cross_trial_stats.py --trials 1 2 3   # specific trials only
```

Output is written to `cross_trial_results/`.

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
| `SD_THRESHOLD_ALL` | `0.040` | CV threshold — all three channels must exceed |
| `SD_THRESHOLD_TWO` | `0.080` | CV threshold — any two channels must exceed |
| `SD_MIN_CONSECUTIVE` | `2` | Consecutive frames required before SD detection fires |
| `SLOPE_MIN_CONSECUTIVE` | `2` | Consecutive frames required before slope detection fires |
| `LIKERT_POSITIVE_THRESHOLD` | `4` | Minimum Likert level to count as a positive event |
| `LIKERT_MIN_CONSECUTIVE` | `2` | Consecutive frames required before Likert detection fires |
| `RATIO_THRESHOLD` | `1.50` | Fixed R/G threshold for Time-to-Threshold metric |
| `UTI_POSITIVE_THRESHOLD_CFU` | `1e5` | Clinical positive cutoff (CFU/mL) |

> **Note:** SD thresholds and slope parameters (`SLOPE_A/B/C_SINGLE`) are currently calibrated for unlocked auto-exposure data. After the first locked-exposure reagent-normalized trial, lower the SD thresholds and recalibrate the slope parameters using the measured angle values per concentration group. Recalibrate `LIKERT_REFERENCE_COLORS` by reading the mean R, G, B values from `rgb_timeseries.csv` for the negative control (Level 1) and 10⁸ well at t=120 min (Level 5).

## Output

### Per-trial results (`trial-N/results/`)

| File | Contents |
|---|---|
| `rgb_timeseries.csv` | Normalized R/G, R/B, G/B ratios — one row per group × timepoint |
| `sd_channel_values.csv` | CV values at every timepoint for all channels + detection summary |
| `slope_scores.csv` | Per-timepoint slope scores + detection summary |
| `likert_scores.csv` | Per-timepoint Likert colour-match level (1–5) + detection summary |
| `detection_events.csv` | First detection time, load tier, concordance — one row per group |
| `total_detection_score.csv` | Per-timepoint score breakdown + full kinetic metrics summary |
| `results_summary.csv` | Complete per-group result with all metrics |
| `statistics.json` | Per-trial statistical test results (Kruskal-Wallis / ANOVA + post-hoc) |
| `performance_summary.png` | Sensitivity/specificity bar chart + confusion matrix for all algorithms |
| `*.png` | Ratio plots, SD CV plots, slope heatmaps, Likert score plots, detection score charts |
| `groups/*.png` | Per-group 4-panel plots: RGB trace, SD, slope, and Likert over time |

### Cross-trial results (`cross_trial_results/`)

| File | Contents |
|---|---|
| `cross_trial_statistics.json` | Pooled ANOVA/Kruskal-Wallis results across all 4 trials (N=4 per group) with Tukey HSD / Dunn's Bonferroni post-hoc where significant |
| `cross_trial_summary.csv` | Per-concentration-group mean ± SD for all key metrics |

## Detection Algorithms

### SD Algorithm (Malanoski et al., 2016)
Monitors the rolling coefficient of variation (CV = SD / mean) over a 9-timepoint window. Detection fires when all three ratio channels exceed `SD_THRESHOLD_ALL`, or any two exceed `SD_THRESHOLD_TWO`, for ≥ 2 consecutive frames. A 30-minute guard period suppresses early-assay noise.

### Slope Algorithm (Ciaccheri et al., 2023)
Computes the cosine of the angle between a fixed 7-point baseline regression (t = 0–30 min) and a sliding 4-point current window. An adaptive threshold `θ = A·exp(−mean/B) + C` accounts for the mean ratio value. Each channel is scored 0–2 based on angular divergence and regression quality (R²). Detection fires when the sum exceeds a weighted threshold for ≥ 2 consecutive frames. Scoring begins at t = 50 min to ensure non-overlapping baseline and current windows.

### Likert Colour-Matching Algorithm
At each timepoint the mean RGB value of the well is compared by Euclidean distance to five reference colours stored in `config.LIKERT_REFERENCE_COLORS`. These were sampled directly from actual trial photos: Level 1 from the sterile negative control at t = 120 min (unreacted blue-purple resazurin), Level 5 from the 10⁸ well at t = 120 min averaged across all four trials (fully converted bright pink resorufin). Detection fires when the matched level ≥ `LIKERT_POSITIVE_THRESHOLD` (default 4) for ≥ `LIKERT_MIN_CONSECUTIVE` consecutive frames. This algorithm is unique in that it operates on absolute RGB colour space rather than ratio variance or slope — it asks "does this well look pink?" directly.

| Level | Label | Interpretation | RGB |
|---|---|---|---|
| 1 | No change | Negative | (94, 86, 198) — blue-purple |
| 2 | Slight change | Likely negative | (80, 68, 175) |
| 3 | Moderate change | Ambiguous | (92, 30, 154) — deep purple |
| 4 | Clear change | Likely positive | (210, 50, 158) — purple-pink |
| 5 | Definite change | Positive | (240, 60, 135) — bright hot pink |

### Total Detection Score
`Total Detection Score = SD Event Count + Slope Total Weighted Score + Likert Event Count`

A continuous variable that captures sustained evidence of color conversion from three complementary perspectives: variance (SD), directional slope, and direct colour matching.

## Cross-Trial Statistical Analysis

With four trials completed, each concentration group has N=4 independent replicates — enough for group-level statistical testing. `cross_trial_stats.py` pools all trials and applies the same normality-adaptive test selection as `stats.py`: Shapiro-Wilk per group, Levene's test for equal variances, then one-way ANOVA if both pass or Kruskal-Wallis otherwise.

**Results across trials 1–4 (N=4 per group, 9 concentration groups):**

| Metric | Test | Statistic | p-value | Significant |
|---|---|---|---|---|
| Total detection score | Kruskal-Wallis | H = 17.60 | 0.025 | ✓ |
| R/G AUC | Kruskal-Wallis | H = 17.27 | 0.027 | ✓ |
| R/B AUC | One-way ANOVA | F = 8.55 | < 0.001 | ✓ |
| G/B AUC | Kruskal-Wallis | H = 8.41 | 0.394 | ✗ |
| Time to R/G ≥ 1.50 | Kruskal-Wallis | H = 3.67 | 0.598 | ✗ |
| Final R/G ratio | Kruskal-Wallis | H = 18.91 | 0.015 | ✓ |

R/B AUC (the only metric meeting normality and equal-variance assumptions) was analyzed by one-way ANOVA. Tukey HSD post-hoc showed 10⁸ significantly higher than all other groups (p ≤ 0.001). Kruskal-Wallis was used for the remaining metrics due to Shapiro-Wilk or Levene violations in one or more groups with N=4. G/B AUC was not significant, consistent with G and B falling together during conversion and producing a relatively stable ratio.

## Human Observer Scoring

`likert_scale.png` is a colorimetric reference card for human observer validation studies. It shows five levels of resazurin conversion derived from actual well photos across all four trials:

| Level | Label | Interpretation |
|---|---|---|
| 1 | No change | Negative |
| 2 | Slight change | Likely negative |
| 3 | Moderate change | Ambiguous |
| 4 | Clear change | Likely positive |
| 5 | Definite change | Positive |

Score ≥ 4 = test positive. Designed for inclusion in research poster presentations.

## References

- Malanoski, A.P. et al. (2016). Simultaneous Identification of Multiple Analytes Using a Self-Reporting Spectroscopic Hybridization Assay. *Analytical Chemistry.*
- Ciaccheri, L. et al. (2023). A Machine Learning Approach for the Rapid Identification of Bacterial Contamination. *Sensors.*
- Ramirez, R. et al. Resazurin-based methods for evaluating bacterial viability. *Journal of Microbiological Methods.*
