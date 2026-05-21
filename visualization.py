# =============================================================================
# visualization.py — All plots for the assay analyzer pipeline.
# =============================================================================

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config


def _save(fig, path: str):
    fig.savefig(path, dpi=config.PLOT_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_rgb_traces(results: list[dict], out_dir: str):
    """RGB channel values over time, one subplot per sample."""
    n = len(results)
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), squeeze=False)
    fig.suptitle("RGB Traces Over Time", fontsize=14, fontweight="bold")

    for idx, r in enumerate(results):
        ax = axes[idx // cols][idx % cols]
        t = r["rgb_data"]["times"]
        ax.plot(t, r["rgb_data"]["R"], "r-o", ms=3, label="R")
        ax.plot(t, r["rgb_data"]["G"], "g-o", ms=3, label="G")
        ax.plot(t, r["rgb_data"]["B"], "b-o", ms=3, label="B")
        ax.set_title(r["name"], fontsize=9)
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Mean pixel value")
        ax.legend(fontsize=7)
        ax.set_ylim(0, 270)

    for idx in range(n, rows * cols):
        axes[idx // cols][idx % cols].axis("off")

    plt.tight_layout()
    _save(fig, os.path.join(out_dir, "rgb_traces.png"))


def plot_ratios(results: list[dict], out_dir: str):
    """R/G and R/B ratios over time."""
    n = len(results)
    cols = min(3, n)
    rows = (n + cols - 1) // cols

    for ratio_key, title, fname in [
        ("RG_ratio", "R/G Ratio", "rg_ratio.png"),
        ("RB_ratio", "R/B Ratio", "rb_ratio.png"),
        ("GB_ratio", "G/B Ratio", "gb_ratio.png"),
    ]:
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), squeeze=False)
        fig.suptitle(title, fontsize=14, fontweight="bold")
        for idx, r in enumerate(results):
            ax = axes[idx // cols][idx % cols]
            t = r["rgb_data"]["times"]
            vals = r["rgb_data"].get(ratio_key, [])
            if vals and vals[0] is not None:
                ax.plot(t, vals, "m-o", ms=3)
            ax.set_title(r["name"], fontsize=9)
            ax.set_xlabel("Time (min)")
            ax.set_ylabel(ratio_key)
        for idx in range(n, rows * cols):
            axes[idx // cols][idx % cols].axis("off")
        plt.tight_layout()
        _save(fig, os.path.join(out_dir, fname))


def plot_sd_ratios(results: list[dict], out_dir: str):
    """SD-algorithm weighted ratios over time."""
    n = len(results)
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), squeeze=False)
    fig.suptitle("SD Algorithm — CV (SD/mean) per Channel", fontsize=14, fontweight="bold")

    for idx, r in enumerate(results):
        ax = axes[idx // cols][idx % cols]
        sd = r["sd_result"]
        t = sd["times"]
        ax.plot(t, sd["ratio_R"],  "r-",  label="R",   lw=1.5)
        ax.plot(t, sd["ratio_G"],  "g-",  label="G",   lw=1.5)
        ax.plot(t, sd["ratio_B"],  "b-",  label="B",   lw=1.5)
        ax.plot(t, sd["ratio_RG"], "r--", label="R/G", lw=1.0)
        ax.plot(t, sd["ratio_RB"], "c--", label="R/B", lw=1.0)
        ax.plot(t, sd["ratio_GB"], color="orange", ls="--", label="G/B", lw=1.0)
        ax.axhline(config.SD_THRESHOLD_TWO, color="gray", ls="--", lw=1, label="Thr2")
        ax.axhline(config.SD_THRESHOLD_ALL, color="black", ls=":", lw=1, label="ThrAll")
        detected = r["sd_detected"]
        ax.set_title(f"{r['name']} [{'DET' if detected else 'neg'}]", fontsize=9)
        ax.set_xlabel("Time (min)")
        ax.legend(fontsize=5)

    for idx in range(n, rows * cols):
        axes[idx // cols][idx % cols].axis("off")
    plt.tight_layout()
    _save(fig, os.path.join(out_dir, "sd_ratios.png"))


def plot_slope_scores(results: list[dict], out_dir: str):
    """Slope algorithm weighted scores over time."""
    n = len(results)
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), squeeze=False)
    fig.suptitle("Slope Algorithm — Weighted Scores", fontsize=14, fontweight="bold")

    for idx, r in enumerate(results):
        ax = axes[idx // cols][idx % cols]
        sl = r["slope_result"]
        t = sl["times"]
        ax.plot(t, sl["weighted_scores"], "purple", lw=1.5)
        ax.axhline(config.SLOPE_WEIGHTED_SCORE_THRESHOLD, color="red", ls="--", lw=1)
        detected = r["slope_detected"]
        ax.set_title(f"{r['name']} [{'DET' if detected else 'neg'}]", fontsize=9)
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Weighted score")

    for idx in range(n, rows * cols):
        axes[idx // cols][idx % cols].axis("off")
    plt.tight_layout()
    _save(fig, os.path.join(out_dir, "slope_scores.png"))



def plot_performance_summary(metrics: dict, out_dir: str):
    """
    Grouped bar chart of Sensitivity, Specificity, PPV, NPV for all algorithms
    (SD, Slope, Likert).  Confusion matrix counts shown below each bar.
    """
    algos        = list(metrics.keys())
    metric_names = ["sensitivity_pct", "specificity_pct", "PPV_pct", "NPV_pct"]
    labels       = ["Sensitivity", "Specificity", "PPV", "NPV"]
    n_algos      = len(algos)
    x            = np.arange(len(labels))
    width        = 0.22
    # Centre the group of bars on each x tick
    offsets      = [(i - (n_algos - 1) / 2) * width for i in range(n_algos)]

    algo_colors = {
        "SD":     "#4C9BE8",
        "Slope":  "#E87B4C",
        "Likert": "#6EC46E",
    }

    fig, ax = plt.subplots(figsize=(10, 5))

    for i, algo in enumerate(algos):
        vals = [metrics[algo].get(m, 0) for m in metric_names]
        vals = [0 if (v != v) else v for v in vals]   # replace NaN with 0
        color = algo_colors.get(algo, f"C{i}")
        bars = ax.bar(x + offsets[i], vals, width, label=algo, color=color,
                      edgecolor="white", linewidth=0.5)

        # Annotate each bar: percentage on top, confusion matrix fraction below
        m = metrics[algo]
        cm_labels = {
            "sensitivity_pct": f"TP={m['TP']} FN={m['FN']}",
            "specificity_pct": f"TN={m['TN']} FP={m['FP']}",
            "PPV_pct":         f"TP={m['TP']} FP={m['FP']}",
            "NPV_pct":         f"TN={m['TN']} FN={m['FN']}",
        }
        for bar, val, mname in zip(bars, vals, metric_names):
            if val > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 1.2,
                    f"{val:.0f}%",
                    ha="center", va="bottom", fontsize=7.5, fontweight="bold",
                )
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                2,
                cm_labels[mname],
                ha="center", va="bottom", fontsize=6, color="white",
                rotation=90,
            )

    ax.set_title("Algorithm Performance — Sensitivity / Specificity / PPV / NPV",
                 fontsize=12, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("Percentage (%)")
    ax.set_ylim(0, 118)
    ax.axhline(100, color="green", ls="--", lw=0.8, alpha=0.5)
    ax.legend(title="Algorithm", fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    _save(fig, os.path.join(out_dir, "performance_summary.png"))


def plot_per_group(results: list[dict], out_dir: str):
    """
    One figure per group — shows every chart for that group on a single page.
    Saves as  results/groups/<safe_name>.png
    """
    groups_dir = os.path.join(out_dir, "groups")
    os.makedirs(groups_dir, exist_ok=True)

    for r in results:
        name = r["name"]
        t = r["rgb_data"]["times"]

        fig, axes = plt.subplots(4, 2, figsize=(12, 17))
        conc = r["concentration_cfu_ml"]
        conc_str = "0 (control)" if conc == 0 else f"{conc:.0e} CFU/mL"
        sd_det = "DETECTED" if r["sd_detected"] else "not detected"
        sl_det = "DETECTED" if r["slope_detected"] else "not detected"
        is_pos = r.get("ground_truth_positive", False)
        clinical_str = "CLINICAL POSITIVE (≥10⁵ CFU/mL)" if is_pos else "CLINICAL NEGATIVE (<10⁵ CFU/mL)"
        title_color = "#cc0000" if is_pos else "#003399"
        fig.suptitle(
            f"Group: {name}   |   Concentration: {conc_str}   |   {clinical_str}\n"
            f"SD algorithm: {sd_det}   |   Slope algorithm: {sl_det}",
            fontsize=12, fontweight="bold", y=0.98, color=title_color,
        )

        # ── (0,0) RGB traces ──────────────────────────────────────────────────
        ax = axes[0][0]
        ax.plot(t, r["rgb_data"]["R"], "r-o", ms=4, label="R")
        ax.plot(t, r["rgb_data"]["G"], "g-o", ms=4, label="G")
        ax.plot(t, r["rgb_data"]["B"], "b-o", ms=4, label="B")
        ax.set_title("RGB Channel Traces")
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Mean pixel value (0–255)")
        ax.set_ylim(0, 270)
        ax.legend()
        ax.grid(True, alpha=0.3)

        # ── (0,1) All three ratios (R/G, R/B, G/B) ───────────────────────────
        ax = axes[0][1]
        for ratio_key, color, label in [
            ("RG_ratio", "m",      "R/G"),
            ("RB_ratio", "c",      "R/B"),
            ("GB_ratio", "orange", "G/B"),
        ]:
            vals = r["rgb_data"].get(ratio_key, [])
            if vals and vals[0] is not None:
                ax.plot(t, vals, color=color, marker="o", ms=3, label=label)
        ax.set_title("Color Ratios (R/G, R/B, G/B)")
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Ratio value")
        ax.legend()
        ax.grid(True, alpha=0.3)

        # ── (1,0) SD algorithm — raw channels (diagnostic only) ──────────────
        ax = axes[1][0]
        sd = r["sd_result"]
        ax.plot(sd["times"], sd["ratio_R"], "r-",  lw=1.5, label="R")
        ax.plot(sd["times"], sd["ratio_G"], "g-",  lw=1.5, label="G")
        ax.plot(sd["times"], sd["ratio_B"], "b-",  lw=1.5, label="B")
        ax.axhline(config.SD_THRESHOLD_TWO, color="gray",  ls="--", lw=1,
                   label=f"Thr2={config.SD_THRESHOLD_TWO}")
        ax.axhline(config.SD_THRESHOLD_ALL, color="black", ls=":",  lw=1,
                   label=f"ThrAll={config.SD_THRESHOLD_ALL}")
        ax.set_title("SD — Raw Channels (diagnostic, not used for detection)")
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("SD / mean (CV)")
        ax.legend(fontsize=6)
        ax.grid(True, alpha=0.3)

        # ── (1,1) SD algorithm — ratio channels (detection) ──────────────────
        ax = axes[1][1]
        ax.plot(sd["times"], sd["ratio_RG"], "r--",    lw=1.5, label="R/G")
        ax.plot(sd["times"], sd["ratio_RB"], "c--",    lw=1.5, label="R/B")
        ax.plot(sd["times"], sd["ratio_GB"], color="orange", ls="--", lw=1.5, label="G/B")
        ax.axhline(config.SD_THRESHOLD_TWO, color="gray",  ls="--", lw=1,
                   label=f"Thr2={config.SD_THRESHOLD_TWO}")
        ax.axhline(config.SD_THRESHOLD_ALL, color="black", ls=":",  lw=1,
                   label=f"ThrAll={config.SD_THRESHOLD_ALL}")
        if r["sd_detected"] and r.get("sd_first_detection_min") is not None:
            ax.axvline(r["sd_first_detection_min"], color="black", ls="-.", lw=1.2,
                       label=f"Det @ {r['sd_first_detection_min']:.0f} min")
        ax.set_title(f"SD — Ratio Channels (detection) [{sd_det}]")
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("SD / mean (CV)")
        ax.legend(fontsize=6)
        ax.grid(True, alpha=0.3)

        # ── (2,0) Slope algorithm total score ────────────────────────────────
        ax = axes[2][0]
        sl = r["slope_result"]
        ax.plot(sl["times"], sl["weighted_scores"], color="purple", lw=1.5,
                label="Total score")
        ax.axhline(config.SLOPE_WEIGHTED_SCORE_THRESHOLD, color="red", ls="--", lw=1,
                   label=f"Threshold={config.SLOPE_WEIGHTED_SCORE_THRESHOLD}")
        if r["slope_detected"] and r.get("slope_first_detection_min") is not None:
            ax.axvline(r["slope_first_detection_min"], color="black", ls="-.", lw=1.2,
                       label=f"Det @ {r['slope_first_detection_min']:.0f} min")
        ax.set_title(f"Slope Algorithm  [{sl_det}]")
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Total score (sum across channels)")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        # ── (2,1) Cumulative total detection score over time ─────────────────
        ax = axes[2][1]
        sd_flags     = sd.get("detection_flags", [])
        slope_scores = sl.get("weighted_scores", [])
        n_tp = len(t)
        sd_flags     = list(sd_flags)     + [0]   * (n_tp - len(sd_flags))
        slope_scores = list(slope_scores) + [0.0] * (n_tp - len(slope_scores))
        cumulative, running = [], 0.0
        for sd_f, sl_s in zip(sd_flags, slope_scores):
            running += int(bool(sd_f)) + float(sl_s)
            cumulative.append(running)
        ax.fill_between(t, cumulative, alpha=0.25, color="purple")
        ax.plot(t, cumulative, color="purple", lw=1.5, label="Cumulative total")
        ax.set_title("Cumulative Total Detection Score")
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Cumulative score")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        # ── (3,0) + (3,1) Summary text box ───────────────────────────────────
        ax = axes[3][0]
        ax.axis("off")

        def _fmt(v):
            if v is None:
                return "N/A"
            try:
                return f"{float(v):.4f}"
            except (TypeError, ValueError):
                return str(v)

        ttr = r.get("time_to_ratio_threshold")
        ttr_str = f"{int(ttr)} min" if ttr is not None else "never"
        rg_auc = r.get("rg_auc")
        rb_auc = r.get("rb_auc")
        gb_auc = r.get("gb_auc")
        auc_rg_str = f"{rg_auc:.1f}" if rg_auc is not None else "N/A"
        auc_rb_str = f"{rb_auc:.1f}" if rb_auc is not None else "N/A"
        auc_gb_str = f"{gb_auc:.1f}" if gb_auc is not None else "N/A"
        tds_str = f"{r['total_detection_score']:.2f}" if r.get("total_detection_score") is not None else "N/A"
        roc = r.get("rate_of_change_rg")
        roc_str = f"{roc:.5f} ratio/min" if roc is not None else "N/A"
        concordance = r.get("concordance_status", "N/A")
        disc_min = r.get("concordance_discordance_min")
        disc_str = f" ({disc_min:.0f} min gap)" if disc_min is not None else ""
        clinical_label = "CLINICAL POSITIVE (≥10⁵ CFU/mL)" if r.get("ground_truth_positive") else "CLINICAL NEGATIVE (<10⁵ CFU/mL)"
        box_color = "#ffe6e6" if r.get("ground_truth_positive") else "#e6f0ff"

        summary_lines = [
            f"Group:                {name}",
            f"Concentration:        {conc_str}",
            f"Clinical status:      {clinical_label}",
            "",
            f"SD algorithm:         {sd_det}",
        ]
        if r["sd_detected"] and r.get("sd_first_detection_min") is not None:
            summary_lines.append(f"  First detect:       {r['sd_first_detection_min']:.0f} min")
        summary_lines += [
            f"  Event count:        {r.get('sd_event_count', 0)}",
            "",
            f"Slope algorithm:      {sl_det}",
        ]
        if r["slope_detected"] and r.get("slope_first_detection_min") is not None:
            summary_lines.append(f"  First detect:       {r['slope_first_detection_min']:.0f} min")
        summary_lines += [
            f"  Total score:        {_fmt(r.get('slope_total_score'))}",
            "",
            f"Total detection score:{tds_str}",
            f"Load tier:            {r.get('load_tier', '—')}",
            f"Concordance:          {concordance}{disc_str}",
            "",
            f"AUC  R/G:             {auc_rg_str}",
            f"AUC  R/B:             {auc_rb_str}",
            f"AUC  G/B:             {auc_gb_str}",
            f"Time to threshold:    {ttr_str}",
            f"Rate of change (R/G): {roc_str}",
            "",
            f"Final R/G ratio:      {_fmt(r.get('final_RG_ratio'))}",
            f"Final R/B ratio:      {_fmt(r.get('final_RB_ratio'))}",
            f"Ref-normalized:       {r['rgb_data'].get('normalized', False)}",
        ]
        ax.text(
            0.05, 0.95, "\n".join(summary_lines),
            transform=ax.transAxes,
            fontsize=10, verticalalignment="top",
            fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor=box_color, alpha=0.8),
        )
        axes[3][1].axis("off")

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        safe_name = name.replace("/", "_").replace(" ", "_")
        out_path = os.path.join(groups_dir, f"{safe_name}.png")
        _save(fig, out_path)


def plot_total_detection_scores(results: list[dict], out_dir: str):
    """
    Stacked bar chart: SD event count (bottom) + slope total score (top) per group.
    Total bar height = total_detection_score.
    Groups are sorted by concentration so dose-response is visible left-to-right.
    """
    sorted_results = sorted(results, key=lambda r: r["concentration_cfu_ml"])
    labels   = [r["name"] for r in sorted_results]
    sd_vals  = [r.get("sd_event_count", 0) for r in sorted_results]
    sl_vals  = [r.get("slope_total_score", 0.0) for r in sorted_results]
    lk_vals  = [r.get("likert_event_count", 0) for r in sorted_results]
    totals   = [r.get("total_detection_score", 0.0) for r in sorted_results]
    gt       = [r.get("ground_truth_positive", False) for r in sorted_results]

    # Bottom of each stacked layer
    sd_bottoms = [0] * len(labels)
    sl_bottoms = [s for s in sd_vals]
    lk_bottoms = [s + l for s, l in zip(sd_vals, sl_vals)]

    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 1.4), 5))

    # Shade background: blue for clinical negatives, red for clinical positives
    for i, is_pos in enumerate(gt):
        ax.axvspan(i - 0.5, i + 0.5,
                   color="#ffe6e6" if is_pos else "#e6f0ff",
                   alpha=0.4, zorder=0)

    # Draw a dividing line between the last negative and first positive group
    neg_indices = [i for i, p in enumerate(gt) if not p]
    pos_indices = [i for i, p in enumerate(gt) if p]
    if neg_indices and pos_indices:
        boundary = (max(neg_indices) + min(pos_indices)) / 2
        ax.axvline(boundary, color="black", lw=1.5, ls="--", alpha=0.6,
                   label="Clinical threshold (10⁵ CFU/mL)")

    ax.bar(x, sd_vals, bottom=sd_bottoms, label="SD event count",    color="#4C9BE8", zorder=2)
    ax.bar(x, sl_vals, bottom=sl_bottoms, label="Slope total score",  color="#E87B4C", zorder=2)
    ax.bar(x, lk_vals, bottom=lk_bottoms, label="Likert event count", color="#6EC46E", zorder=2)

    # Annotate total above each bar
    for i, total in enumerate(totals):
        ax.text(x[i], total + 0.3, f"{total:.1f}", ha="center", va="bottom", fontsize=8)

    # Axis labels: mark positives/negatives
    tick_labels = []
    for lbl, is_pos in zip(labels, gt):
        marker = "▲" if is_pos else "▽"
        tick_labels.append(f"{marker} {lbl}")

    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Score")
    ax.set_title(
        "Total Detection Score per Group\n"
        "(SD event count + slope weighted score + Likert event count)   "
        "▽ = clinical negative   ▲ = clinical positive",
        fontweight="bold"
    )
    ax.legend()
    ax.grid(axis="y", alpha=0.3, zorder=1)
    plt.tight_layout()
    _save(fig, os.path.join(out_dir, "total_detection_scores.png"))


def plot_likert_scores(results: list[dict], out_dir: str):
    """
    Step plot of Likert colour-match level (1–5) over time for all groups.

    Each group is one line.  Clinical positives (ground_truth_positive=True)
    are drawn in warm colours; negatives in cool colours.  A horizontal dashed
    line marks the positive threshold.  The right-hand edge shows the reference
    colour swatch for each Likert level.
    """
    ref_colors = getattr(config, "LIKERT_REFERENCE_COLORS",
                         [(94,86,198),(80,68,175),(92,30,154),(179,46,166),(188,2,152)])
    threshold  = int(getattr(config, "LIKERT_POSITIVE_THRESHOLD", 4))
    n_levels   = len(ref_colors)

    pos_results = [r for r in results if r.get("ground_truth_positive")]
    neg_results = [r for r in results if not r.get("ground_truth_positive")]

    # Colour palettes — warm for positives, cool for negatives.
    warm = plt.cm.Reds(np.linspace(0.4, 0.9, max(len(pos_results), 1)))
    cool = plt.cm.Blues(np.linspace(0.4, 0.9, max(len(neg_results), 1)))
    colour_map = {}
    for i, r in enumerate(neg_results):
        colour_map[r["name"]] = cool[i]
    for i, r in enumerate(pos_results):
        colour_map[r["name"]] = warm[i]

    fig, ax = plt.subplots(figsize=(13, 5))

    for r in results:
        lk = r["likert_result"]
        times  = lk["times"]
        scores = lk["likert_scores"]
        lw   = 2.0 if r.get("ground_truth_positive") else 1.2
        ls   = "-" if r.get("ground_truth_positive") else "--"
        col  = colour_map[r["name"]]
        ax.step(times, scores, where="post", color=col, lw=lw, ls=ls,
                label=r["name"], zorder=3)

        # Mark first detection
        fd = r.get("likert_first_detection_min")
        if fd is not None:
            ax.axvline(fd, color=col, alpha=0.4, lw=0.8, zorder=2)

    # Positive threshold line
    ax.axhline(threshold - 0.5, color="#c0392b", lw=1.5, ls=":", zorder=4,
               label=f"Positive threshold (≥{threshold})")

    # Right-side colour swatches showing what each level looks like
    ax_right = ax.twinx()
    ax_right.set_ylim(ax.get_ylim())
    ax_right.set_yticks(list(range(1, n_levels + 1)))
    ax_right.set_yticklabels([])
    # Draw small coloured rectangles as tick labels
    for level_idx, rgb in enumerate(ref_colors):
        level = level_idx + 1
        color_norm = tuple(c / 255 for c in rgb)
        ax_right.axhspan(level - 0.4, level + 0.4,
                         xmin=0.98, xmax=1.0,
                         color=color_norm, alpha=0.9, zorder=5,
                         transform=ax_right.get_xaxis_transform())

    ax.set_xlabel("Time (min)", fontsize=11)
    ax.set_ylabel("Likert level", fontsize=11)
    ax.set_yticks(list(range(1, n_levels + 1)))
    ax.set_yticklabels([f"Level {i}" for i in range(1, n_levels + 1)], fontsize=9)
    ax.set_ylim(0.5, n_levels + 0.5)
    ax.set_title(
        "Likert Algorithm — Colour-Match Level Over Time\n"
        "(warm/solid = clinical positive, cool/dashed = clinical negative,  "
        f"threshold = ≥{threshold})",
        fontweight="bold"
    )
    ax.legend(bbox_to_anchor=(1.06, 1), loc="upper left", fontsize=8,
              borderaxespad=0, framealpha=0.9)
    ax.grid(axis="y", alpha=0.25, zorder=1)
    plt.tight_layout()
    _save(fig, os.path.join(out_dir, "likert_scores.png"))


def generate_all_plots(results: list[dict], metrics: dict, out_dir: str):
    """Generate all plots and save to out_dir."""
    plot_rgb_traces(results, out_dir)
    plot_ratios(results, out_dir)
    plot_sd_ratios(results, out_dir)
    plot_slope_scores(results, out_dir)
    plot_likert_scores(results, out_dir)
    plot_total_detection_scores(results, out_dir)
    plot_performance_summary(metrics, out_dir)
    plot_per_group(results, out_dir)
