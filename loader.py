"""
loader.py
---------
Loads and merges user data based on config.yaml settings.
You do not need to edit this file.
"""

import pandas as pd
import yaml
import sys


def load_config(config_path: str = "config.yaml") -> dict:
    try:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"\nERROR: Config file not found at '{config_path}'")
        print("Make sure config.yaml is in the same folder as run.py")
        sys.exit(1)
    return config


def load_and_merge(config: dict) -> tuple:
    """
    Load features and outcomes, compute progression.
    Returns (df, feature_cols) where feature_cols is the list
    of biomarker column names to analyse.
    """
    pid = config["patient_id_column"]
    bl_col = config["baseline_score_column"]
    fu_col = config["followup_score_column"]
    threshold = config["progression_threshold"]

    # Load features
    print(f"Loading features from: {config['features_file']}")
    try:
        features = pd.read_csv(config["features_file"], low_memory=False)
    except FileNotFoundError:
        print(f"\nERROR: Features file not found: {config['features_file']}")
        sys.exit(1)
    features.columns = features.columns.str.strip()

    if pid not in features.columns:
        print(f"\nERROR: Patient ID column '{pid}' not found in features file.")
        print(f"Available columns (first 10): {list(features.columns[:10])}")
        sys.exit(1)

    print(f"  {len(features)} patients, {features.shape[1]-1} features")

    # Load outcomes
    print(f"Loading outcomes from: {config['outcome_file']}")
    try:
        outcomes = pd.read_csv(config["outcome_file"], low_memory=False)
    except FileNotFoundError:
        print(f"\nERROR: Outcome file not found: {config['outcome_file']}")
        sys.exit(1)
    outcomes.columns = outcomes.columns.str.strip()

    for col, key in [(bl_col, "baseline_score_column"), (fu_col, "followup_score_column")]:
        if col not in outcomes.columns:
            print(f"\nERROR: Column '{col}' not found in outcome file.")
            print(f"Available columns: {list(outcomes.columns)}")
            sys.exit(1)

    # Compute progression
    outcomes["PROGRESSION"] = pd.to_numeric(outcomes[fu_col], errors="coerce") - \
                               pd.to_numeric(outcomes[bl_col], errors="coerce")
    outcomes["PROGRESSOR"] = (outcomes["PROGRESSION"] > threshold).astype(int)
    outcomes = outcomes.dropna(subset=["PROGRESSION"])

    progressors = outcomes["PROGRESSOR"].sum()
    non_progressors = (outcomes["PROGRESSOR"] == 0).sum()
    print(f"  Progression threshold: >{threshold} points")
    print(f"  Progressors: {progressors}, Non-progressors: {non_progressors}")

    # Merge
    df = pd.merge(features, outcomes[[pid, bl_col, "PROGRESSION", "PROGRESSOR"]],
                  on=pid, how="inner")
    print(f"  Patients after merge: {len(df)}")

    if len(df) == 0:
        print(f"\nERROR: No patients matched between features and outcome files.")
        sys.exit(1)

    # Identify feature columns (exclude ID, outcome, and meta columns)
    exclude = {pid, bl_col, "PROGRESSION", "PROGRESSOR"}
    feature_cols = [c for c in df.columns
                    if c not in exclude
                    and pd.api.types.is_numeric_dtype(df[c])]

    print(f"  Biomarker features to analyse: {len(feature_cols)}")

    return df, feature_cols
