"""
run.py
------
Main entry point for neuro-biomarker-analysis.

Usage:
    python run.py
    python run.py --config my_config.yaml

This script runs the full statistical analysis pipeline:
  1. Load your data
  2. Compute Spearman correlations with progression
  3. Compare biomarker levels between progressors and non-progressors
  4. Generate volcano plot, box plots, and bar charts
  5. Save full results tables as CSV files

All settings are controlled by config.yaml.
You do not need to edit this file.
"""

import os
import argparse
import time

from loader import load_config, load_and_merge
from analysis import spearman_correlations, group_comparison, combine_results
from plotter import (plot_volcano, plot_top_correlations, plot_boxplots,
                     plot_effect_sizes, save_tables)


def main(config_path: str = "config.yaml"):
    start = time.time()

    print("=" * 55)
    print("  Neuro-Biomarker-Analysis")
    print("  Statistical biomarker analysis for neurological")
    print("  disease progression")
    print("=" * 55)

    # Load config
    print(f"\nLoading config from: {config_path}")
    config = load_config(config_path)
    os.makedirs(config["output_directory"], exist_ok=True)

    # Step 1 - Load data
    print("\n[1/4] Loading data...")
    df, feature_cols = load_and_merge(config)

    # Step 2 - Correlation analysis
    print("\n[2/4] Running correlation analysis...")
    corr_results = spearman_correlations(df, feature_cols, config)

    # Step 3 - Group comparison
    print("\n[3/4] Running group comparison...")
    group_results = group_comparison(df, feature_cols, config)
    combined_results = combine_results(corr_results, group_results)

    # Step 4 - Generate outputs
    print("\n[4/4] Generating figures and tables...")
    plot_volcano(corr_results, config)
    plot_top_correlations(corr_results, config)
    plot_boxplots(df, group_results, config)
    plot_effect_sizes(group_results, config)
    save_tables(corr_results, group_results, combined_results, config)

    # Done
    elapsed = time.time() - start
    print(f"\n{'=' * 55}")
    print(f"  Analysis complete in {elapsed:.1f} seconds")
    print(f"  Results saved to: {config['output_directory']}/")
    print(f"\n  Output files:")
    print(f"    volcano_plot.png       - Effect size vs significance")
    print(f"    top_correlations.png   - Top biomarkers by correlation")
    print(f"    boxplots.png           - Group distributions")
    print(f"    effect_sizes.png       - Group difference effect sizes")
    print(f"    correlation_results.csv  - Full correlation table")
    print(f"    group_comparison_results.csv - Full group comparison table")
    print(f"    combined_results.csv   - Combined ranked summary")
    print(f"    results_summary.txt    - Plain text summary")
    print(f"{'=' * 55}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Neuro-Biomarker-Analysis Pipeline")
    parser.add_argument("--config", type=str, default="config.yaml")
    args = parser.parse_args()
    main(config_path=args.config)
