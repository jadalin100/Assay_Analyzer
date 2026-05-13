# =============================================================================
# power_analysis.py — Sample size estimation for diagnostic assay validation.
#
# Run standalone:  python3 power_analysis.py
#
# Answers the question: "How many replicate trials do we need so that our
# estimated sensitivity and specificity have acceptably narrow confidence
# intervals?"
#
# Definitions for this experiment
# --------------------------------
# One "trial"  = one complete assay run: all 6 groups photographed at all
#                9 time points, then analysed by the pipeline.
# n trials     = n independent experiment days (or repeat runs).
#
# Sensitivity and specificity are proportions (0–1).  The Wilson score
# interval is used throughout — it is more accurate than the normal
# approximation for small n and proportions near 0 or 1.
# =============================================================================

import math
import sys


# ---------------------------------------------------------------------------
# Configuration — edit these to match your expectations
# ---------------------------------------------------------------------------

N_POSITIVE_GROUPS  = 4    # Groups that are true positives (10^5, 10^6, 10^7, 10^8)
N_NEGATIVE_GROUPS  = 2    # Groups that are true negatives (Sterile ctrl, Reagent ctrl)

# Your best guess at true sensitivity / specificity based on pilot data or
# literature.  If unknown, use 0.90 as a conservative assumption.
EXPECTED_SENSITIVITY = 0.90   # 90% — conservative; change to 1.0 if you expect perfect
EXPECTED_SPECIFICITY = 0.95   # 95%

# Desired half-width of the 95% confidence interval (e.g. 0.15 = ±15 %)
TARGET_CI_HALF_WIDTH = 0.20   # ±20 % is adequate for a pilot study

# Range of replicates to evaluate
N_MIN = 1
N_MAX = 30

ALPHA = 0.05   # 1 - confidence level

# ---------------------------------------------------------------------------
# Core statistics
# ---------------------------------------------------------------------------

def wilson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """
    Wilson score confidence interval for a proportion k/n.
    Returns (lower, upper).  Handles n=0 gracefully.
    """
    if n == 0:
        return 0.0, 1.0
    z = _z(alpha)
    p_hat = k / n
    denom = 1 + z ** 2 / n
    centre = (p_hat + z ** 2 / (2 * n)) / denom
    margin = (z * math.sqrt(p_hat * (1 - p_hat) / n + z ** 2 / (4 * n ** 2))) / denom
    return max(0.0, centre - margin), min(1.0, centre + margin)


def _z(alpha: float) -> float:
    """Two-sided z critical value via rational approximation (Abramowitz & Stegun)."""
    p = 1 - alpha / 2
    t = math.sqrt(-2 * math.log(1 - p))
    c0, c1, c2 = 2.515517, 0.802853, 0.010328
    d1, d2, d3 = 1.432788, 0.189269, 0.001308
    return t - (c0 + c1 * t + c2 * t ** 2) / (1 + d1 * t + d2 * t ** 2 + d3 * t ** 3)


def half_width(k: int, n: int, alpha: float = 0.05) -> float:
    lo, hi = wilson_ci(k, n, alpha)
    return (hi - lo) / 2


def min_n_for_target(p_expected: float, target_hw: float, alpha: float = 0.05,
                     n_max: int = 200) -> int:
    """
    Return the smallest n such that the expected CI half-width ≤ target_hw,
    assuming the observed proportion equals p_expected.
    """
    for n in range(1, n_max + 1):
        k = round(p_expected * n)
        if half_width(k, n, alpha) <= target_hw:
            return n
    return n_max


# ---------------------------------------------------------------------------
# Per-trial totals
# Each trial contributes N_POSITIVE_GROUPS binary outcomes to sensitivity
# and N_NEGATIVE_GROUPS binary outcomes to specificity.
# After T trials: n_pos = T × N_POSITIVE_GROUPS, n_neg = T × N_NEGATIVE_GROUPS
# ---------------------------------------------------------------------------

def trials_for_target(p_expected: float, n_groups: int,
                      target_hw: float, alpha: float = 0.05,
                      max_trials: int = 50) -> int:
    """Minimum number of *trials* to reach target CI half-width."""
    for t in range(1, max_trials + 1):
        n = t * n_groups
        k = round(p_expected * n)
        if half_width(k, n, alpha) <= target_hw:
            return t
    return max_trials


# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------

def _bar(value: float, width: int = 20, char: str = "█") -> str:
    filled = round(value * width)
    return char * filled + "·" * (width - filled)


def print_ci_table():
    z_val = _z(ALPHA)
    print()
    print("=" * 78)
    print("  POWER ANALYSIS — Smartphone Colorimetric Assay Validation")
    print("=" * 78)
    print(f"  Setup  : {N_POSITIVE_GROUPS} positive groups × {N_NEGATIVE_GROUPS} negative groups per trial")
    print(f"  Target : 95% CI half-width ≤ ±{TARGET_CI_HALF_WIDTH:.0%}")
    print(f"  Assumed sensitivity: {EXPECTED_SENSITIVITY:.0%}   Assumed specificity: {EXPECTED_SPECIFICITY:.0%}")
    print()

    # ── Sensitivity table ────────────────────────────────────────────────────
    print("  ── SENSITIVITY (positive groups) ──────────────────────────────")
    print(f"  {'Trials':>7}  {'n pos obs':>10}  {'k expected':>11}  "
          f"{'95% CI':>20}  {'Half-width':>11}  {'Goal met':>8}")
    print("  " + "-" * 74)

    for t in range(N_MIN, N_MAX + 1):
        n = t * N_POSITIVE_GROUPS
        k = round(EXPECTED_SENSITIVITY * n)
        lo, hi = wilson_ci(k, n, ALPHA)
        hw = (hi - lo) / 2
        met = "✓" if hw <= TARGET_CI_HALF_WIDTH else ""
        ci_str = f"[{lo:.3f}, {hi:.3f}]"
        print(f"  {t:>7}  {n:>10}  {k:>11}  {ci_str:>20}  {hw:>10.3f}   {met:>8}")
        if t >= 5 and met:
            break   # stop once goal is met and we've shown at least 5 rows

    t_sens = trials_for_target(EXPECTED_SENSITIVITY, N_POSITIVE_GROUPS,
                               TARGET_CI_HALF_WIDTH, ALPHA)
    print(f"\n  → Minimum trials for ±{TARGET_CI_HALF_WIDTH:.0%} CI on sensitivity: {t_sens}")

    # ── Specificity table ─────────────────────────────────────────────────────
    print()
    print("  ── SPECIFICITY (negative groups) ───────────────────────────────")
    print(f"  {'Trials':>7}  {'n neg obs':>10}  {'k expected':>11}  "
          f"{'95% CI':>20}  {'Half-width':>11}  {'Goal met':>8}")
    print("  " + "-" * 74)

    for t in range(N_MIN, N_MAX + 1):
        n = t * N_NEGATIVE_GROUPS
        k = round(EXPECTED_SPECIFICITY * n)
        lo, hi = wilson_ci(k, n, ALPHA)
        hw = (hi - lo) / 2
        met = "✓" if hw <= TARGET_CI_HALF_WIDTH else ""
        ci_str = f"[{lo:.3f}, {hi:.3f}]"
        print(f"  {t:>7}  {n:>10}  {k:>11}  {ci_str:>20}  {hw:>10.3f}   {met:>8}")
        if t >= 5 and met:
            break

    t_spec = trials_for_target(EXPECTED_SPECIFICITY, N_NEGATIVE_GROUPS,
                               TARGET_CI_HALF_WIDTH, ALPHA)
    print(f"\n  → Minimum trials for ±{TARGET_CI_HALF_WIDTH:.0%} CI on specificity: {t_spec}")

    # ── Summary ───────────────────────────────────────────────────────────────
    t_needed = max(t_sens, t_spec)
    print()
    print("=" * 78)
    print(f"  RECOMMENDATION")
    print(f"  {'─' * 74}")
    print(f"  To achieve ±{TARGET_CI_HALF_WIDTH:.0%} confidence intervals on BOTH metrics:")
    print(f"  Run at least  {t_needed}  complete replicate trial(s).")
    print()
    print(f"  Context:")
    print(f"    3 trials  — proof-of-concept / feasibility (wide CIs, ~±30–40%)")
    print(f"    5 trials  — adequate pilot study (CIs narrow enough to interpret)")
    print(f"   10 trials  — preliminary data suitable for a research report")
    print(f"   20+ trials — rigorous clinical-grade validation")
    print()
    print(f"  Note: specificity is harder to pin down with only {N_NEGATIVE_GROUPS} negative")
    print(f"  groups per trial. Adding a third negative control (e.g. buffer only)")
    print(f"  would reduce the number of trials needed.")
    print("=" * 78)
    print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Allow quick override from command line:
    #   python3 power_analysis.py --sensitivity 0.85 --specificity 1.0 --target 0.15
    import argparse
    parser = argparse.ArgumentParser(description="Power analysis for assay validation")
    parser.add_argument("--sensitivity", type=float, default=EXPECTED_SENSITIVITY)
    parser.add_argument("--specificity", type=float, default=EXPECTED_SPECIFICITY)
    parser.add_argument("--target",      type=float, default=TARGET_CI_HALF_WIDTH,
                        help="Target CI half-width, e.g. 0.15 for ±15%%")
    parser.add_argument("--n_pos",       type=int,   default=N_POSITIVE_GROUPS)
    parser.add_argument("--n_neg",       type=int,   default=N_NEGATIVE_GROUPS)
    args = parser.parse_args()

    EXPECTED_SENSITIVITY = args.sensitivity
    EXPECTED_SPECIFICITY = args.specificity
    TARGET_CI_HALF_WIDTH = args.target
    N_POSITIVE_GROUPS    = args.n_pos
    N_NEGATIVE_GROUPS    = args.n_neg

    print_ci_table()
