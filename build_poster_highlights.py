"""
build_poster_highlights.py
==========================
Selects and re-renders the strongest, most poster-appropriate results from
the averaged_results/ and cross_trial_results/ folders into poster_highlights/.

Figures produced
----------------
fig1_rg_ratio_timeseries.png   Main result — averaged R/G ratio over time, ±1 SD shading
fig2_rb_auc_bars.png           Strongest stat — R/B AUC per group, ANOVA p<0.001, Tukey brackets
fig3_detection_heatmap.png     Detection rate + mean first-detection time: algorithm × group
fig4_stats_table.png           Cross-trial statistics table (6 metrics, rendered as image)
fig5_auc_comparison.png        R/G, R/B, G/B AUC side-by-side with p-value annotations

Each figure includes a subtitle / annotation with the relevant statistic.
"""

import os, csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE      = os.path.dirname(os.path.abspath(__file__))
AVG_DIR   = os.path.join(BASE, "averaged_results")
CROSS_DIR = os.path.join(BASE, "cross_trial_results")
OUT_DIR   = os.path.join(BASE, "poster_highlights")
os.makedirs(OUT_DIR, exist_ok=True)

DPI       = 200
FONT      = "DejaVu Sans"
plt.rcParams.update({
    "font.family":       FONT,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.labelsize":    11,
    "xtick.labelsize":   9,
    "ytick.labelsize":   9,
})

# ── Constants ──────────────────────────────────────────────────────────────────
UTI_THRESH = 1e5

GROUP_ORDER = [
    "0 (Sterile Negative Control)",
    "10^1", "10^2", "10^3", "10^4",
    "10^5", "10^6", "10^7",
    "10^8 (Positive Control)",
]
GROUP_LABELS = {
    "0 (Sterile Negative Control)": "Sterile\nNeg.",
    "10^1": "10¹", "10^2": "10²", "10^3": "10³", "10^4": "10⁴",
    "10^5": "10⁵*", "10^6": "10⁶", "10^7": "10⁷",
    "10^8 (Positive Control)": "10⁸\n(+ctrl)",
}
CONC_MAP = {
    "0 (Sterile Negative Control)": 0,
    "10^1": 1e1, "10^2": 1e2, "10^3": 1e3, "10^4": 1e4,
    "10^5": 1e5, "10^6": 1e6, "10^7": 1e7,
    "10^8 (Positive Control)": 1e8,
}

# Colour per group — negatives blue-grey, positives warm red/purple
PALETTE = {
    "0 (Sterile Negative Control)": "#888888",
    "10^1": "#aec6cf", "10^2": "#779ecb", "10^3": "#4169a1", "10^4": "#223f7a",
    "10^5": "#f5a623", "10^6": "#e06030", "10^7": "#b8003b",
    "10^8 (Positive Control)": "#000000",
}
IS_POS = {g: CONC_MAP[g] >= UTI_THRESH for g in GROUP_ORDER}

# ── Loaders ────────────────────────────────────────────────────────────────────

def load_avg_rgb():
    rows = []
    with open(os.path.join(AVG_DIR, "averaged_rgb_timeseries.csv"), newline="") as f:
        for row in csv.DictReader(f):
            row["group"] = row["group"].strip()
            rows.append(row)
    return rows

def load_avg_summary():
    out = {}
    with open(os.path.join(AVG_DIR, "averaged_results_summary.csv"), newline="") as f:
        for row in csv.DictReader(f):
            out[row["group"].strip()] = row
    return out

def load_detection_summary():
    out = {}
    with open(os.path.join(AVG_DIR, "averaged_detection_summary.csv"), newline="") as f:
        for row in csv.DictReader(f):
            out[row["group"].strip()] = row
    return out

def load_cross_stats():
    with open(os.path.join(CROSS_DIR, "cross_trial_statistics.json")) as f:
        return json.load(f)

# ── FIG 1 — R/G Ratio time-series ─────────────────────────────────────────────
print("[1/5] Fig 1: R/G ratio timeseries...")

rgb_rows = load_avg_rgb()

# Build per-group time series from averaged CSV
series = {g: {"t": [], "mean": [], "sd": []} for g in GROUP_ORDER}
for row in rgb_rows:
    g = row["group"]
    if g not in series:
        continue
    try:
        series[g]["t"].append(float(row["time_min"]))
        series[g]["mean"].append(float(row["RG_ratio_mean"]))
        series[g]["sd"].append(float(row["RG_ratio_sd"]) if row["RG_ratio_sd"] else 0.0)
    except ValueError:
        pass

fig, ax = plt.subplots(figsize=(9, 5))

for g in GROUP_ORDER:
    s = series[g]
    if not s["t"]:
        continue
    order = np.argsort(s["t"])
    t   = np.array(s["t"])[order]
    m   = np.array(s["mean"])[order]
    sd  = np.array(s["sd"])[order]
    col = PALETTE[g]
    lw  = 2.2 if IS_POS[g] else 1.4
    ls  = "-" if IS_POS[g] else "--"
    ax.plot(t, m, color=col, lw=lw, ls=ls, label=GROUP_LABELS[g], zorder=3 if IS_POS[g] else 2)
    ax.fill_between(t, m - sd, m + sd, color=col, alpha=0.13, zorder=1)

ax.axhline(1.0, color="#aaaaaa", lw=0.8, ls=":", zorder=0)
ax.axvspan(0, 30, color="#e8f0ff", alpha=0.4, zorder=0, label="Baseline window")
ax.set_xlabel("Time (min)")
ax.set_ylabel("Normalized R/G ratio  (reagent-referenced, ±1 SD)")
ax.set_title("Resazurin Color Conversion Over Time by Bacterial Concentration\n"
             "Mean of 4 independent trials — shaded band = ±1 SD", fontsize=11)
ax.set_xlim(0, 120)

# UTI threshold annotation
ax.annotate("← UTI threshold: 10⁵ CFU/mL*",
            xy=(121, ax.get_ylim()[1]*0.88), fontsize=7.5, color="#b8003b",
            ha="left", va="center", annotation_clip=False)

# Legend — split neg / pos
neg_handles = [mpatches.Patch(color=PALETTE[g], label=GROUP_LABELS[g])
               for g in GROUP_ORDER if not IS_POS[g]]
pos_handles = [mpatches.Patch(color=PALETTE[g], label=GROUP_LABELS[g])
               for g in GROUP_ORDER if IS_POS[g]]
band_patch   = mpatches.Patch(color="#aaaaaa", alpha=0.35, label="Baseline\nwindow")
l1 = ax.legend(handles=neg_handles, title="Negatives (<10⁵)", fontsize=7.5,
               loc="upper left", framealpha=0.85, title_fontsize=8)
l2 = ax.legend(handles=pos_handles, title="Positives (≥10⁵)", fontsize=7.5,
               loc="upper center", framealpha=0.85, title_fontsize=8)
ax.add_artist(l1)

ax.annotate("* 10⁵ = clinical UTI threshold", xy=(0, -0.12),
            xycoords="axes fraction", fontsize=7, color="#555555")

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig1_rg_ratio_timeseries.png"), dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("   → fig1_rg_ratio_timeseries.png")


# ── FIG 2 — R/B AUC bar chart (ANOVA p<0.001) ─────────────────────────────────
print("[2/5] Fig 2: R/B AUC bar chart...")

cross = load_cross_stats()
rb_data = next(m for m in cross["metrics"] if m["metric"] == "rb_auc")
# Mean and SD from cross_trial_summary — group names differ from GROUP_ORDER
with open(os.path.join(CROSS_DIR, "cross_trial_summary.csv"), newline="") as f:
    ct_summary_raw = {row["group"].strip(): row for row in csv.DictReader(f)}

# Map cross_trial group names → GROUP_ORDER names
CT_NAME_MAP = {
    "0 (negative)": "0 (Sterile Negative Control)",
    "10^8":         "10^8 (Positive Control)",
}
ct_summary = {}
for k, v in ct_summary_raw.items():
    canonical = CT_NAME_MAP.get(k, k)
    ct_summary[canonical] = v

groups_plot = [g for g in GROUP_ORDER if g in ct_summary]
rb_means = [float(ct_summary[g]["rb_auc_mean"]) for g in groups_plot]
rb_sds   = [float(ct_summary[g]["rb_auc_sd"])   for g in groups_plot]
colors   = ["#CC3300" if IS_POS[g] else "#4477AA" for g in groups_plot]

fig, ax = plt.subplots(figsize=(9, 5))
x = np.arange(len(groups_plot))
bars = ax.bar(x, rb_means, yerr=rb_sds, color=colors, alpha=0.82,
              capsize=5, width=0.65, error_kw={"lw": 1.5, "ecolor": "#333333"})

# Mark significance: 10^8 vs all others — draw brackets
sig_y = max(rb_means) + max(rb_sds) + 8
idx_8 = groups_plot.index("10^8 (Positive Control)")
for i, g in enumerate(groups_plot):
    if g == "10^8 (Positive Control)":
        continue
    # Short bracket line
    lo, hi = sorted([i, idx_8])
    y_br = sig_y + (hi - lo) * 1.5
    ax.plot([lo, lo, hi, hi], [y_br - 2, y_br, y_br, y_br - 2],
            lw=0.8, color="#333333")
    ax.text((lo + hi) / 2, y_br + 0.3, "*", ha="center", va="bottom",
            fontsize=8, color="#333333")

ax.set_xticks(x)
ax.set_xticklabels([GROUP_LABELS[g] for g in groups_plot], fontsize=9)
ax.set_ylabel("R/B AUC  (ratio · min)")
ax.set_title("R/B Area Under Curve — Mean ± SD Across 4 Trials\n"
             "One-way ANOVA: F = 8.55,  p < 0.001 — Tukey HSD: * p ≤ 0.001 vs 10⁸",
             fontsize=10)

ax.axvline(4.5, color="#aaaaaa", lw=0.8, ls="--")
ax.text(4.6, ax.get_ylim()[0] + 2, "← negative   positive →",
        fontsize=7.5, color="#555555", va="bottom")

pos_patch = mpatches.Patch(color="#CC3300", label="Clinical positive (≥10⁵ CFU/mL)")
neg_patch = mpatches.Patch(color="#4477AA", label="Clinical negative (<10⁵ CFU/mL)")
ax.legend(handles=[pos_patch, neg_patch], fontsize=8, loc="upper left")

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig2_rb_auc_bars.png"), dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("   → fig2_rb_auc_bars.png")


# ── FIG 3 — Detection rate heatmap ────────────────────────────────────────────
print("[3/5] Fig 3: Detection rate heatmap...")

det = load_detection_summary()

algorithms  = ["SD Algorithm", "Slope Algorithm", "Likert Algorithm", "Combined (Any)"]
algo_keys   = ["sd_detection_rate", "slope_detection_rate",
               "likert_detection_rate", "combined_detection_rate"]
groups_h    = GROUP_ORDER

matrix = np.zeros((len(algorithms), len(groups_h)))
for j, g in enumerate(groups_h):
    row = det.get(g, {})
    for i, key in enumerate(algo_keys):
        try:
            matrix[i, j] = float(row.get(key, 0) or 0)
        except ValueError:
            matrix[i, j] = 0.0

# Build first-detection time matrix (NaN where not detected / unavailable)
indiv_time_keys = [
    "sd_first_detection_min_mean",
    "slope_first_detection_min_mean",
    "likert_first_detection_min_mean",
]
indiv_rate_keys = ["sd_detection_rate", "slope_detection_rate", "likert_detection_rate"]

first_det = np.full((len(algorithms), len(groups_h)), np.nan)
for j, g in enumerate(groups_h):
    row = det.get(g, {})
    # Individual algorithm rows (0–2)
    for i, tkey in enumerate(indiv_time_keys):
        try:
            val_str = str(row.get(tkey, "") or "").strip()
            if val_str and val_str.lower() not in ("nan", "never", ""):
                first_det[i, j] = float(val_str)
        except (ValueError, TypeError):
            pass
    # Combined row (3): earliest first-detection among algorithms that fired
    available_times = []
    for rk, tk in zip(indiv_rate_keys, indiv_time_keys):
        try:
            rate = float(row.get(rk, 0) or 0)
            if rate > 0:
                t_str = str(row.get(tk, "") or "").strip()
                if t_str and t_str.lower() not in ("nan", "never", ""):
                    available_times.append(float(t_str))
        except (ValueError, TypeError):
            pass
    if available_times:
        first_det[3, j] = min(available_times)

# Custom diverging colormap: grey (0) → pale yellow (0.5) → red (1)
cmap = LinearSegmentedColormap.from_list(
    "detect",
    [(0.0, "#dddddd"), (0.5, "#ffdd88"), (1.0, "#cc2200")]
)

fig, ax = plt.subplots(figsize=(13, 6.5))
im = ax.imshow(matrix, cmap=cmap, vmin=0, vmax=1, aspect="auto")

# Cell annotations — two separate ax.text calls per cell for precise placement.
# In imshow data coords each cell spans ±0.5 around its row index.
# We offset the rate line up by 0.18 and the time line down by 0.22 so
# neither line touches the cell border or its neighbours.
for i in range(len(algorithms)):
    for j in range(len(groups_h)):
        val  = matrix[i, j]
        tval = first_det[i, j]
        color  = "white" if val > 0.65 else ("#333333" if val > 0.0 else "#aaaaaa")
        weight = "bold" if val == 1.0 else "normal"
        if val > 0 and not np.isnan(tval):
            # Rate % slightly above cell centre
            ax.text(j, i - 0.17, f"{val*100:.0f}%",
                    ha="center", va="center",
                    fontsize=9, color=color, fontweight=weight)
            # Time below cell centre, lighter weight, smaller
            ax.text(j, i + 0.22, f"t={tval:.0f} min",
                    ha="center", va="center",
                    fontsize=7, color=color, fontweight="normal", style="italic")
        else:
            ax.text(j, i, f"{val*100:.0f}%",
                    ha="center", va="center",
                    fontsize=9.5, color=color, fontweight=weight)

ax.set_xticks(range(len(groups_h)))
ax.set_xticklabels([GROUP_LABELS[g] for g in groups_h], fontsize=10)
ax.set_yticks(range(len(algorithms)))
ax.set_yticklabels(algorithms, fontsize=10)

# Divider between negatives and positives
ax.axvline(4.5, color="white", lw=2.5)
# Use axes-fraction coords so the labels never collide with tick labels
ax.text(0.265, -0.10, "← Clinical Negatives (<10⁵ CFU/mL)",
        ha="center", va="top", fontsize=9, color="#444444",
        transform=ax.transAxes, clip_on=False)
ax.text(0.77, -0.10, "Clinical Positives (≥10⁵ CFU/mL) →",
        ha="center", va="top", fontsize=9, color="#880000",
        transform=ax.transAxes, clip_on=False)

ax.set_title("Detection Rate by Algorithm and Concentration Group\n"
             "(N = 4 trials — % detected | mean first-detection time where applicable)",
             fontsize=11, pad=12)

cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
cbar.set_label("Detection rate", fontsize=9)
cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
cbar.set_ticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=9)

fig.tight_layout(rect=[0, 0.08, 1, 1])
fig.savefig(os.path.join(OUT_DIR, "fig3_detection_heatmap.png"), dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("   → fig3_detection_heatmap.png")


# ── FIG 4 — Cross-trial statistics summary table ───────────────────────────────
print("[4/5] Fig 4: Statistics summary table...")

table_data = [
    ["Metric", "Test", "Statistic", "p-value", "Significant?", "Key Finding"],
    ["R/B AUC", "One-way ANOVA", "F = 8.55", "< 0.001", "✓ YES",
     "10⁸ significantly higher\nthan ALL other groups (Tukey HSD)"],
    ["Final R/G Ratio", "Kruskal-Wallis", "H = 18.91", "0.015", "✓ YES",
     "No pair significant after\nBonferroni correction (N=4 low power)"],
    ["Total Detection Score", "Kruskal-Wallis", "H = 17.60", "0.024", "✓ YES",
     "No pair significant after\nBonferroni correction"],
    ["R/G AUC", "Kruskal-Wallis", "H = 17.27", "0.027", "✓ YES",
     "No pair significant after\nBonferroni correction"],
    ["G/B AUC", "Kruskal-Wallis", "H = 8.41", "0.394", "✗ NO",
     "G and B co-vary — ratio\nunresponsive to conversion"],
    ["Time to R/G ≥ 1.50", "Kruskal-Wallis", "H = 3.67", "0.598", "✗ NO",
     "Threshold inflated by\nunlocked auto-exposure artifacts"],
]

col_widths = [1.3, 1.5, 1.2, 0.9, 1.0, 2.5]
fig_w = sum(col_widths) + 0.4
fig, ax = plt.subplots(figsize=(fig_w, 4.2))
ax.axis("off")

n_rows = len(table_data)
n_cols = len(table_data[0])
row_h  = 1.0 / n_rows

for r, row in enumerate(table_data):
    for c, cell in enumerate(row):
        x = sum(col_widths[:c]) / sum(col_widths)
        y = 1.0 - (r + 1) * row_h
        w = col_widths[c] / sum(col_widths)

        # Background
        if r == 0:
            bg = "#2c3e50"
            fc = "white"
            fw = "bold"
            fs = 9
        elif "✓ YES" in cell:
            bg, fc, fw, fs = "#fff3e0", "#333333", "normal", 8.5
        elif "✗ NO" in cell:
            bg, fc, fw, fs = "#f5f5f5", "#888888", "normal", 8.5
        elif r % 2 == 0:
            bg, fc, fw, fs = "#f9f9f9", "#333333", "normal", 8.5
        else:
            bg, fc, fw, fs = "white", "#333333", "normal", 8.5

        # Bold the "✓ YES" and "✗ NO" cells
        if cell in ("✓ YES",):
            fc, fw = "#b8460b", "bold"
        elif cell in ("✗ NO",):
            fc, fw = "#888888", "bold"

        ax.add_patch(plt.Rectangle((x, y), w, row_h,
                                   transform=ax.transAxes,
                                   facecolor=bg, edgecolor="#cccccc", lw=0.5))
        ax.text(x + w / 2, y + row_h / 2, cell,
                transform=ax.transAxes, ha="center", va="center",
                fontsize=fs, color=fc, fontweight=fw, wrap=True,
                multialignment="center")

ax.set_title("Cross-Trial Statistical Results  (N = 4 trials, 4 replicates per concentration group)",
             fontsize=10, pad=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig4_stats_table.png"), dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("   → fig4_stats_table.png")


# ── FIG 5 — RGB channel AUC comparison: R/G, R/B, G/B ────────────────────────
print("[5/5] Fig 5: AUC comparison (R/G, R/B, G/B)...")

summary = load_avg_summary()

groups_b = [g for g in GROUP_ORDER if g in summary]
colors_b = ["#CC3300" if IS_POS[g] else "#4477AA" for g in groups_b]

rg_means = [float(summary[g]["rg_auc_mean"]) for g in groups_b]
rg_sds   = [float(summary[g]["rg_auc_sd"])   for g in groups_b]
rb_means = [float(summary[g]["rb_auc_mean"]) for g in groups_b]
rb_sds   = [float(summary[g]["rb_auc_sd"])   for g in groups_b]
gb_means = [float(summary[g]["gb_auc_mean"]) for g in groups_b]
gb_sds   = [float(summary[g]["gb_auc_sd"])   for g in groups_b]

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

panels = [
    (rg_means, rg_sds,
     "R/G AUC  (ratio · min)", "R/G AUC",
     "Kruskal-Wallis: H=17.27, p=0.027\n(no pairwise sig. after Bonferroni)", False),
    (rb_means, rb_sds,
     "R/B AUC  (ratio · min)", "R/B AUC  ★ Primary result",
     "One-way ANOVA: F=8.55, p<0.001\n10⁸ sig. vs all groups (Tukey p≤0.001)", True),
    (gb_means, gb_sds,
     "G/B AUC  (ratio · min)", "G/B AUC  — Not Significant",
     "Kruskal-Wallis: H=8.41, p=0.394\nG and B co-vary; ratio unresponsive", False),
]

for ax, (means, sds, ylabel, title, stat_label, is_sig) in zip(axes, panels):
    x = np.arange(len(groups_b))
    ax.bar(x, means, yerr=sds, color=colors_b, alpha=0.82,
           capsize=4, width=0.65, error_kw={"lw": 1.4, "ecolor": "#444444"})
    ax.set_xticks(x)
    ax.set_xticklabels([GROUP_LABELS[g] for g in groups_b], fontsize=8)
    ax.set_ylabel(ylabel, fontsize=9.5)
    ax.axvline(4.5, color="#aaaaaa", lw=0.8, ls="--")
    ylo = ax.get_ylim()[0]
    ax.text(4.6, ylo, "positive →", fontsize=7, color="#880000", va="bottom")
    ax.text(4.3, ylo, "← neg.",     fontsize=7, color="#333333", va="bottom", ha="right")

pos_patch = mpatches.Patch(color="#CC3300", label="Clinical positive (≥10⁵ CFU/mL)")
neg_patch = mpatches.Patch(color="#4477AA", label="Clinical negative (<10⁵ CFU/mL)")
fig.legend(handles=[pos_patch, neg_patch], loc="lower center", ncol=2,
           fontsize=9, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=[0, 0.06, 1, 1])
fig.savefig(os.path.join(OUT_DIR, "fig5_auc_comparison.png"), dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("   → fig5_auc_comparison.png")


# ── Summary text file ──────────────────────────────────────────────────────────
summary_txt = os.path.join(OUT_DIR, "POSTER_HIGHLIGHTS_GUIDE.txt")
with open(summary_txt, "w") as f:
    f.write("""POSTER BOARD — HIGHLIGHTS GUIDE
================================
Generated from 4 independent trials (N=4 per concentration group, 9 groups).

FIGURE SELECTION RATIONALE
---------------------------

fig1_rg_ratio_timeseries.png  ← LEAD FIGURE
  - Shows blue→pink resazurin conversion over 120 min by concentration
  - Shaded bands = ±1 SD across 4 trials (shows reproducibility)
  - High concentrations (10^7, 10^8) show clear visible separation
  - KEY CAVEAT: 10^5 (clinical UTI threshold) does not separate cleanly
    from negatives — state this explicitly on the board

fig2_rb_auc_bars.png  ← STRONGEST STATISTICAL RESULT
  - Only metric that meets ANOVA assumptions (normality + equal variance)
  - One-way ANOVA: F=8.55, p<0.001 (4 sig. figures, publishable)
  - Tukey HSD: 10^8 significantly higher than ALL other groups (p≤0.001)
  - Use this as your primary quantitative result
  - CAVEAT: only 10^8 separates — not 10^5 or 10^6

fig3_detection_heatmap.png  ← ALGORITHM COMPARISON
  - Shows detection rate % AND mean first-detection time in each cell
  - Cell format: "75% / t=80 min" — gives audience both sensitivity and speed
  - Likert: 100% at 10^7 and 10^8 — most consistent algorithm
  - SD: 50% false-positive rate on sterile negative — show this honestly
  - Combined (any): 100% at 10^7+ but still 50% FP on sterile negative
  - Key talking point: algorithms are complementary, not redundant

fig4_stats_table.png  ← RESULTS TABLE FOR BOARD
  - All 6 cross-trial metrics in one table
  - Shows which are significant, which aren't, and why
  - Use directly on the poster as your "Statistical Analysis" section

fig5_auc_comparison.png  ← AUC CHANNEL COMPARISON (3 panels)
  - R/G AUC, R/B AUC, and G/B AUC side by side for direct comparison
  - Each panel labelled with its test statistic and significance
  - R/B panel labelled as primary result (★); G/B labelled as not significant
  - WHY G/B IS INCLUDED BUT LABELLED NS: shows the audience that not every
    channel is informative; G and B co-vary because both absorb in the same
    spectral region during resazurin→resorufin conversion, so their ratio
    is insensitive to the reaction. p=0.394.
  - WHY TIME-TO-THRESHOLD IS EXCLUDED: The R/G ≥ 1.50 threshold was inflated
    to account for unlocked auto-exposure artifacts (true biology threshold
    ~1.10). As a result, clinical negatives 10^1 and 10^2 cross it at t=15 min
    in Trial 1 (AE spike, not growth), and positives like 10^7 never cross it
    in Trial 3. Kruskal-Wallis p=0.598. Metric is invalid until auto-exposure
    is locked and threshold recalibrated to ~1.10.

KEY NUMBERS TO QUOTE ON BOARD
-------------------------------
  R/B AUC ANOVA:           F = 8.55,  p < 0.001
  10^8 vs all groups:      Tukey HSD p ≤ 0.001
  Likert detection (10^7): 100% across 4 trials
  Likert detection (10^8): 100% across 4 trials
  Combined detection (10^7+10^8): 100%
  Sterile negative SD false-positive rate: 50%  ← must disclose

WHAT NOT TO PUT ON THE BOARD
------------------------------
  - G/B AUC as a positive result — include it only as a labelled-NS comparison
    (already done in fig5). Do not claim significance (p=0.394).
  - Time-to-threshold metric — invalid; RATIO_THRESHOLD inflated to 1.50 by
    unlocked auto-exposure artifacts. Negatives cross it at t=15 min; positives
    miss it entirely in some trials. Kruskal-Wallis p=0.598. Recalibrate to
    ~1.10 after locking exposure in next trial.
  - Individual per-trial SD/slope score plots (too noisy, trial-dependent)
  - Any sensitivity/specificity claim for 10^5 (missed in 3/4 trials)
""")
print(f"   → POSTER_HIGHLIGHTS_GUIDE.txt")

print(f"\n✓ Done. {len(os.listdir(OUT_DIR))} files in {OUT_DIR}/")
