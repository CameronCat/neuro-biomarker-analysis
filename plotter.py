"""
plotter.py
----------
Generates publication-ready figures and results tables.
You do not need to edit this file.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import os


def plot_volcano(corr_df: pd.DataFrame, config: dict):
    """
    Volcano plot: effect size vs -log10(p-value).
    Highlights significant biomarkers.
    """
    out = config["output_directory"]
    alpha = config.get("alpha", 0.05)
    label = config.get("biomarker_label", "Biomarker")
    top_n = config.get("top_n_biomarkers", 20)

    df = corr_df.copy()
    df["neg_log10_p"] = -np.log10(df["P_value_corrected"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(10, 7))

    # Non-significant points
    non_sig = df[~df["Significant"]]
    ax.scatter(non_sig["Spearman_rho"], non_sig["neg_log10_p"],
               color="lightgrey", alpha=0.5, s=15, label="Not significant")

    # Significant points
    sig = df[df["Significant"]]
    ax.scatter(sig["Spearman_rho"], sig["neg_log10_p"],
               color="steelblue", alpha=0.8, s=25, label=f"Significant (p<{alpha})")

    # Label top features
    top = df.head(min(top_n, 15))
    for _, row in top.iterrows():
        ax.annotate(
            row["Feature"],
            (row["Spearman_rho"], row["neg_log10_p"]),
            fontsize=7, alpha=0.8,
            xytext=(4, 4), textcoords="offset points"
        )

    # Reference lines
    ax.axhline(-np.log10(alpha), color="red", linestyle="--",
               linewidth=0.8, alpha=0.7, label=f"p={alpha}")
    ax.axvline(0, color="black", linestyle="-", linewidth=0.5, alpha=0.3)

    ax.set_xlabel("Spearman Correlation (rho) with Progression", fontsize=12)
    ax.set_ylabel("-log10(corrected p-value)", fontsize=12)
    ax.set_title(f"{label} Association with Neurological Progression\nVolcano Plot", fontsize=13)
    ax.legend(fontsize=9)
    plt.tight_layout()

    path = f"{out}/volcano_plot.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")


def plot_top_correlations(corr_df: pd.DataFrame, config: dict):
    """
    Horizontal bar chart of top biomarkers by Spearman correlation.
    """
    out = config["output_directory"]
    top_n = min(config.get("top_n_biomarkers", 20), len(corr_df))
    label = config.get("biomarker_label", "Biomarker")

    top = corr_df.head(top_n).copy()
    colors = ["steelblue" if r > 0 else "salmon" for r in top["Spearman_rho"]]

    fig, ax = plt.subplots(figsize=(9, max(5, top_n * 0.4)))
    bars = ax.barh(range(top_n), top["Spearman_rho"].values, color=colors, alpha=0.8)
    ax.set_yticks(range(top_n))
    ax.set_yticklabels(top["Feature"].values, fontsize=9)
    ax.invert_yaxis()
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Spearman Correlation with Progression Score", fontsize=11)
    ax.set_title(f"Top {top_n} {label}s by Correlation with Progression", fontsize=12)

    # Mark significant ones with asterisk
    for i, (_, row) in enumerate(top.iterrows()):
        if row["Significant"]:
            ax.text(row["Spearman_rho"] + 0.002, i, "*", fontsize=10, va="center")

    ax.text(0.98, 0.02, "* = significant after correction",
            transform=ax.transAxes, fontsize=8, ha="right", alpha=0.7)

    plt.tight_layout()
    path = f"{out}/top_correlations.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")


def plot_boxplots(df: pd.DataFrame, group_df: pd.DataFrame, config: dict):
    """
    Box plots for top biomarkers comparing progressors vs non-progressors.
    """
    out = config["output_directory"]
    top_n = min(config.get("top_n_biomarkers", 20), 12, len(group_df))
    neg_label = config.get("negative_class_label", "Non-progressor")
    pos_label = config.get("positive_class_label", "Progressor")

    top_features = group_df.head(top_n)["Feature"].tolist()

    cols = 3
    rows = int(np.ceil(top_n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 3.5))
    axes = axes.flatten()

    progressors = df[df["PROGRESSOR"] == 1]
    non_progressors = df[df["PROGRESSOR"] == 0]

    for i, feat in enumerate(top_features):
        ax = axes[i]
        data = [non_progressors[feat].dropna().values,
                progressors[feat].dropna().values]
        bp = ax.boxplot(data, patch_artist=True, widths=0.5,
                        medianprops=dict(color="black", linewidth=2))
        bp["boxes"][0].set_facecolor("lightblue")
        bp["boxes"][1].set_facecolor("salmon")
        ax.set_xticklabels([neg_label, pos_label], fontsize=8)
        ax.set_title(feat, fontsize=9, fontweight="normal")
        ax.set_ylabel("Value", fontsize=8)

        # Add p-value if available
        row = group_df[group_df["Feature"] == feat]
        if len(row) > 0:
            p = row["P_value_corrected"].values[0]
            sig = row["Significant"].values[0]
            p_str = f"p={p:.3f}" if p >= 0.001 else "p<0.001"
            if sig:
                p_str += " *"
            ax.set_title(f"{feat}\n{p_str}", fontsize=8)

    # Hide unused subplots
    for j in range(top_n, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle(f"Top {top_n} Biomarkers: Progressors vs Non-Progressors",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    path = f"{out}/boxplots.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")


def plot_effect_sizes(group_df: pd.DataFrame, config: dict):
    """
    Bar chart of effect sizes (rank-biserial correlation) for top biomarkers.
    """
    out = config["output_directory"]
    top_n = min(config.get("top_n_biomarkers", 20), len(group_df))
    label = config.get("biomarker_label", "Biomarker")

    top = group_df.head(top_n).copy()
    colors = ["steelblue" if r > 0 else "salmon" for r in top["Effect_Size_RBC"]]

    fig, ax = plt.subplots(figsize=(9, max(5, top_n * 0.4)))
    ax.barh(range(top_n), top["Effect_Size_RBC"].values, color=colors, alpha=0.8)
    ax.set_yticks(range(top_n))
    ax.set_yticklabels(top["Feature"].values, fontsize=9)
    ax.invert_yaxis()
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Effect Size (Rank-Biserial Correlation)", fontsize=11)
    ax.set_title(f"Top {top_n} {label}s by Group Difference Effect Size", fontsize=12)

    # Reference lines for effect size interpretation
    for val, label_text in [(0.1, "small"), (0.3, "medium"), (0.5, "large")]:
        ax.axvline(val, color="grey", linestyle=":", linewidth=0.7, alpha=0.5)
        ax.axvline(-val, color="grey", linestyle=":", linewidth=0.7, alpha=0.5)

    plt.tight_layout()
    path = f"{out}/effect_sizes.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")


def save_tables(corr_df: pd.DataFrame, group_df: pd.DataFrame,
                combined_df: pd.DataFrame, config: dict):
    """Save all results as CSV files."""
    out = config["output_directory"]

    corr_path = f"{out}/correlation_results.csv"
    corr_df.to_csv(corr_path, index=False)
    print(f"  Saved: {corr_path}")

    group_path = f"{out}/group_comparison_results.csv"
    group_df.to_csv(group_path, index=False)
    print(f"  Saved: {group_path}")

    combined_path = f"{out}/combined_results.csv"
    combined_df.to_csv(combined_path, index=False)
    print(f"  Saved: {combined_path}")

    # Plain text summary
    summary_path = f"{out}/results_summary.txt"
    alpha = config.get("alpha", 0.05)
    method = config.get("correction_method", "fdr_bh")
    with open(summary_path, "w") as f:
        f.write("Neuro-Biomarker-Analysis Results Summary\n")
        f.write("=" * 45 + "\n\n")
        f.write(f"Multiple testing correction: {method}\n")
        f.write(f"Significance threshold: p < {alpha}\n\n")
        f.write(f"Correlation analysis:\n")
        f.write(f"  Features tested: {len(corr_df)}\n")
        f.write(f"  Significant: {corr_df['Significant'].sum()}\n")
        f.write(f"  Top feature: {corr_df.iloc[0]['Feature']} "
                f"(rho={corr_df.iloc[0]['Spearman_rho']:.3f})\n\n")
        f.write(f"Group comparison:\n")
        f.write(f"  Features tested: {len(group_df)}\n")
        f.write(f"  Significant: {group_df['Significant'].sum()}\n")
        f.write(f"  Top feature: {group_df.iloc[0]['Feature']} "
                f"(effect size={group_df.iloc[0]['Effect_Size_RBC']:.3f})\n\n")
        f.write("See outputs/ folder for full results tables and figures.\n")
    print(f"  Saved: {summary_path}")
