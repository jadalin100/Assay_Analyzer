# Assay Analyzer — Project Context for Claude

Smartphone-based colorimetric UTI detection pipeline. Bacteria reduce resazurin (blue) to resorufin (pink); R/G and R/B ratios rise with conversion. Images are taken every 5 min over 120 min from phone.

**Clinical threshold:** ≥ 10⁵ CFU/mL = UTI positive. Groups 10^1–10^4 are clinical negatives.

---

## Architecture

| File | Role |
|------|------|
| `main.py` | Pipeline entry point — orchestrates all stages |
| `config.py` | All tunable parameters — edit here, not in algorithm files |
| `ingestion.py` | Image loading, ROI selection, RGB extraction |
| `classifier.py` | Runs both algorithms, computes kinetic metrics, assembles result dict |
| `sd_algorithm.py` | SD (variance-based) detection algorithm |
| `slope_algorithm.py` | Slope (angular divergence) detection algorithm |
| `reporting.py` | CSV exports and results summary |
| `visualization.py` | Ratio time-series, CV plots, slope heatmaps |
| `stats.py` | Shapiro-Wilk, Kruskal-Wallis/ANOVA, Tukey/Dunn post-hoc |
| `metrics.py` | AUC, time-to-threshold, rate of change |

---

## Detection Pipeline (main.py stages)

**Stage 1a — Camera artifact scan:**
- Scans the reagent reference folder's R/G ratio for large frame-to-frame jumps (≥ `ARTIFACT_JUMP_THRESHOLD` = 0.30)
- Positive spike → negative recovery pair = auto-exposure (AE) event → marks the FULL window (spike through recovery) as artifact
- Unmatched positive spike = real bacterial conversion → NOT marked
- Unmatched negative recovery = single-frame glitch → marks that frame only
- Same directional logic applied per-sample to avoid marking bacterial conversion as artifact
- Artifact timepoints are excluded from SD rolling windows and slope baseline

**Stage 1b — Reference card normalization:** per-image white-balance correction using physical reference card swatch.

**Stage 2 — Reagent reference normalization:** each sample's ratios divided by the reagent-only control at matching timepoints. Negatives stay near 1.0; converting groups rise above 1.0.

**Stage 3 — Classification:** SD and slope algorithms run independently on normalized ratios.

**Stage 4 — Reporting and visualization.**

---

## SD Algorithm

Rolling coefficient of variation (CV = SD/mean) over a 9-timepoint window.

Detection fires when, for ≥ 2 consecutive frames, EITHER:
- ALL three channels (R/G, R/B, G/B) exceed `SD_THRESHOLD_ALL`, OR
- ANY two channels exceed `SD_THRESHOLD_TWO`

**Current thresholds (calibrated 2026-05-18, unlocked AE data):**
```python
SD_THRESHOLD_ALL = 0.040   # lowered from 0.060 — enables 10^6 detection
SD_THRESHOLD_TWO = 0.080   # raised from 0.075 — eliminates 10^4 false positive
SD_MIN_CONSECUTIVE = 2
```

**Calibration rationale:**
- `thr_all=0.040` catches 10^6 via all-three condition (cv_RG≈0.040, cv_RB≈0.094 at t=115–120)
- `thr_two=0.080` eliminates 10^4 FP (its cv_RB=0.075 was just breaching the old 0.075 threshold)
- 10^5 is still missed by SD — acceptable as borderline clinical threshold; caught weakly by slope

---

## Slope Algorithm

Compares a fixed baseline regression (indices 0–6, t=0–30 min) against a sliding 4-point current window. Scores each ratio channel 0/1/2 based on angular divergence and R²; max score per timepoint = 6.

Detection fires when total weighted score > `SLOPE_WEIGHTED_SCORE_THRESHOLD` for ≥ 2 consecutive frames.

**Current thresholds (calibrated 2026-05-18):**
```python
SLOPE_WEIGHTED_SCORE_THRESHOLD = 4.0   # raised from 1.0
SLOPE_MIN_CONSECUTIVE = 2
SLOPE_BASELINE_POINTS = 7   # indices 0–6 = t=0–30 min
SLOPE_CURRENT_POINTS  = 4   # 20-min rolling window
```

**Calibration rationale:**
- Clinical negatives (10^1–10^4) max out at score 4 at any timepoint
- Clinical positives (10^6–10^8) regularly score 5–6
- Threshold > 4 (strict greater-than) gives clean separation with no false positives
- 10^5 scores ≤ 4 — missed by slope, caught by SD instead

---

## Total Detection Score

```python
total_detection_score = slope_total_weighted_score + sd_event_count
```

**Note:** slope dominates this metric (~3:1 scale ratio). This score is kept for convenience in the classification report. For ANOVA, use `rg_auc` (already computed in `classifier.py`) — it is a direct physical measurement, threshold-independent, and stable across trials.

---

## Kinetic Metrics (classifier.py)

| Metric | Description |
|--------|-------------|
| `rg_auc` / `rb_auc` / `gb_auc` | Area under ratio curve (ratio·min, trapezoidal rule) — **primary ANOVA variable** |
| `time_to_ratio_threshold` | First timepoint where normalized R/G ≥ `RATIO_THRESHOLD` (1.50) |
| `rate_of_change_rg` | Linear slope of R/G over last 4 timepoints (ratio/min) |
| `concordance_status` | Agreement between SD and slope first_detection_min |

---

## Recalibration Checklist (do after first locked-exposure + covered-reference trial)

- [ ] **SD thresholds** — with locked AE, negative CV should drop to ~0.018. Lower `SD_THRESHOLD_ALL` toward 0.025–0.035 and `SD_THRESHOLD_TWO` toward 0.030–0.045.
- [ ] **Slope thresholds** — with clean data, negatives should score 0. `SLOPE_WEIGHTED_SCORE_THRESHOLD` can drop toward 2–3.
- [ ] **LOAD_TIERS** — update 60-min and 90-min cutoffs using actual `first_detection_min` for 10^6 and 10^7.
- [ ] **RATIO_THRESHOLD** — currently 1.50 (inflated due to unlocked AE artifacts). Lower toward 1.10 once negatives stay near 1.0.
- [ ] **SLOPE_RAPID_CHANGE_FALLBACK** — if color change first appears earlier than t=50 min, recalculate index.
- [ ] **ARTIFACT_JUMP_THRESHOLD** — currently 0.30 (calibrated for unlocked AE). May need adjustment with locked exposure.

---

## Known Issues / Design Decisions

- **10^5 is missed by both algorithms** — borderline clinical threshold; considered acceptable for now. Revisit after locked-exposure trial.
- **Reference card must be covered/stable** — photo-converted reference card inflated baselines in the 2026-05-18 trial. Cover the reference card during incubation.
- **`total_detection_score` is slope-dominated** — not suitable as primary ANOVA variable. Use `rg_auc` instead.
- **Android AE lock** — use the padlock icon in the camera app to lock exposure and white balance before starting the trial.

---

## Recent Commit History (as of 2026-05-18)

- `a61b2c1` Recalibrate SD and slope thresholds from 2026-05-18 trial data
- `311f1ab` Fix per-sample artifact scan to not flag real bacterial conversion as artifact
- `f62ec08` Fix contamination window: mark all frames between AE spike and recovery as artifacts
- `dc725d2` Artifact-aware SD rolling window reset
