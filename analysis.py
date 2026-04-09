"""
analysis.py
-----------
Performs correlation analysis and group comparison statistics.
You do not need to edit this file.
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests


def spearman_correlations(df: pd.DataFrame, feature_cols: list, config: dict) -> pd.DataFrame:
    """
    Compute Spearman correlation between each biomarker and
    the continuous progression score.
    Returns a dataframe sorted by absolute correlation strength.
    """
    print(f"  Computing Spearman correlations with progression score...")

    results = []
    progression = df["PROGRESSION"]

    for col in feature_cols:
        values = df[col].dropna()
        common_idx = values.index.intersection(progression.dropna().index)
        if len(common_idx) < 10:
            continue
        rho, pval = stats.spearmanr(values[common_idx], progression[common_idx])
        results.append({
            "Feature": col,
            "Spearman_rho": rho,
            "Abs_rho": abs(rho),
            "P_value_corr": pval,
            "N": len(common_idx),
        })

    results_df = pd.DataFrame(results)

    # Multiple testing correction
    method = config.get("correction_method", "fdr_bh")
    if method != "none" and len(results_df) > 0:
        _, p_corrected, _, _ = multipletests(
            results_df["P_value_corr"], method=method
        )
        results_df["P_value_corrected"] = p_corrected
        results_df["Significant"] = p_corrected < config.get("alpha", 0.05)
    else:
        results_df["P_value_corrected"] = results_df["P_value_corr"]
        results_df["Significant"] = results_df["P_value_corr"] < config.get("alpha", 0.05)

    results_df = results_df.sort_values("Abs_rho", ascending=False).reset_index(drop=True)
    results_df["Rank"] = results_df.index + 1

    n_sig = results_df["Significant"].sum()
    print(f"  Significant correlations (corrected p<{config.get('alpha', 0.05)}): {n_sig}")

    return results_df


def group_comparison(df: pd.DataFrame, feature_cols: list, config: dict) -> pd.DataFrame:
    """
    Compare each biomarker between progressors and non-progressors
    using Mann-Whitney U test with rank-biserial correlation effect size.
    Returns a dataframe sorted by effect size.
    """
    print(f"  Computing group comparisons (Mann-Whitney U test)...")

    progressors = df[df["PROGRESSOR"] == 1]
    non_progressors = df[df["PROGRESSOR"] == 0]

    results = []
    for col in feature_cols:
        prog_vals = progressors[col].dropna().values
        non_prog_vals = non_progressors[col].dropna().values

        if len(prog_vals) < 5 or len(non_prog_vals) < 5:
            continue

        stat, pval = stats.mannwhitneyu(
            prog_vals, non_prog_vals, alternative="two-sided"
        )

        # Rank-biserial correlation as effect size
        n1, n2 = len(prog_vals), len(non_prog_vals)
        rbc = 1 - (2 * stat) / (n1 * n2)

        results.append({
            "Feature": col,
            "Mean_Progressors": prog_vals.mean(),
            "Mean_NonProgressors": non_prog_vals.mean(),
            "Mean_Difference": prog_vals.mean() - non_prog_vals.mean(),
            "Effect_Size_RBC": rbc,
            "Abs_Effect_Size": abs(rbc),
            "U_statistic": stat,
            "P_value_group": pval,
            "N_Progressors": n1,
            "N_NonProgressors": n2,
        })

    results_df = pd.DataFrame(results)

    # Multiple testing correction
    method = config.get("correction_method", "fdr_bh")
    if method != "none" and len(results_df) > 0:
        _, p_corrected, _, _ = multipletests(
            results_df["P_value_group"], method=method
        )
        results_df["P_value_corrected"] = p_corrected
        results_df["Significant"] = p_corrected < config.get("alpha", 0.05)
    else:
        results_df["P_value_corrected"] = results_df["P_value_group"]
        results_df["Significant"] = results_df["P_value_group"] < config.get("alpha", 0.05)

    results_df = results_df.sort_values("Abs_Effect_Size", ascending=False).reset_index(drop=True)
    results_df["Rank"] = results_df.index + 1

    n_sig = results_df["Significant"].sum()
    print(f"  Significant group differences (corrected p<{config.get('alpha', 0.05)}): {n_sig}")

    return results_df


def combine_results(corr_df: pd.DataFrame, group_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge correlation and group comparison results into one summary table.
    """
    corr_sub = corr_df[["Feature", "Spearman_rho", "Abs_rho",
                          "P_value_corrected", "Significant"]].rename(columns={
        "P_value_corrected": "Corr_P_corrected",
        "Significant": "Corr_Significant"
    })
    group_sub = group_df[["Feature", "Effect_Size_RBC", "Abs_Effect_Size",
                            "Mean_Progressors", "Mean_NonProgressors",
                            "P_value_corrected", "Significant"]].rename(columns={
        "P_value_corrected": "Group_P_corrected",
        "Significant": "Group_Significant"
    })

    combined = pd.merge(corr_sub, group_sub, on="Feature", how="outer")
    combined["Overall_Score"] = (
        combined["Abs_rho"].fillna(0) + combined["Abs_Effect_Size"].fillna(0)
    ) / 2
    combined = combined.sort_values("Overall_Score", ascending=False).reset_index(drop=True)
    combined["Overall_Rank"] = combined.index + 1

    return combined
