[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19478292.svg)](https://doi.org/10.5281/zenodo.19478292)
# neuro-biomarker-analysis

A transparent statistical analysis pipeline for identifying biomarkers associated with neurological disease progression.

No machine learning. No black boxes. Just clear, reproducible statistics that any reviewer can understand and verify.

---

## What this does

Given a dataset with baseline patient measurements and a follow-up score, this pipeline will:

1. Compute Spearman correlations between every biomarker and the progression score
2. Compare biomarker levels between progressors and non-progressors using Mann-Whitney U tests
3. Apply multiple testing correction (Benjamini-Hochberg FDR by default)
4. Generate a volcano plot, box plots, correlation bar charts, and effect size charts
5. Save full results tables as CSV files ready for a paper

---

## Why statistics instead of machine learning

Machine learning models can be hard to interpret, require large sample sizes, and produce results that are difficult to explain to reviewers or clinicians. This pipeline uses established statistical methods that are transparent, widely understood, and appropriate for clinical biomarker discovery:

- Spearman correlation is non-parametric and robust to outliers
- Mann-Whitney U test does not assume normally distributed data
- Rank-biserial correlation gives an interpretable effect size
- Benjamini-Hochberg correction controls the false discovery rate

---

## Who this is for

- Clinical researchers with longitudinal patient data
- Neurology studies with baseline biomarkers and follow-up assessments
- Studies of Parkinson's disease, multiple sclerosis, ALS, Alzheimer's disease, or any other progressive neurological condition
- Researchers who want a transparent, reviewer-friendly analysis

---

## Requirements

- A Windows, Mac, or Linux computer
- Python 3.9 or later (free download from python.org)
- Your data in CSV format

---

## Installation

### Step 1 - Install Python

Download from https://www.python.org/downloads/. On Windows, check "Add Python to PATH" during installation.

### Step 2 - Download this tool

Click the green "Code" button on GitHub and select "Download ZIP". Extract to your Desktop.

### Step 3 - Open a terminal

Windows: Search for PowerShell in the Start menu.
Mac: Open Terminal from Applications > Utilities.

### Step 4 - Navigate to the folder

    cd C:\Users\YourName\Desktop\neuro-biomarker-analysis

### Step 5 - Install required packages

    pip install -r requirements.txt

---

## Quick start with example data

The repository includes synthetic example data so you can test immediately:

    python run.py

This analyses 400 synthetic patients with 500 biomarker features and saves all results to the outputs/ folder in about 4 seconds.

---

## Using your own data

### Step 1 - Prepare your CSV files

Features file - one row per patient, one column per biomarker:

    PATIENT_ID, AGE, SEX, BIOMARKER_001, BIOMARKER_002, ...
    PD_0001, 65, 1, 1.23, 0.45, ...

Outcome file - baseline and follow-up scores:

    PATIENT_ID, SCORE_BASELINE, SCORE_FOLLOWUP
    PD_0001, 22.0, 25.0

### Step 2 - Edit config.yaml

Open config.yaml in any text editor and update:

    features_file: data/your_features.csv
    outcome_file: data/your_outcome.csv
    patient_id_column: YOUR_PATIENT_ID_COLUMN
    baseline_score_column: YOUR_BASELINE_COLUMN
    followup_score_column: YOUR_FOLLOWUP_COLUMN

### Step 3 - Run

    python run.py

---

## Output files

All results are saved in the outputs/ folder:

- volcano_plot.png - Spearman rho vs -log10(p-value) for all biomarkers
- top_correlations.png - Bar chart of top biomarkers by correlation strength
- boxplots.png - Distributions of top biomarkers in progressors vs non-progressors
- effect_sizes.png - Effect size bar chart for top biomarkers
- correlation_results.csv - Full table with rho, p-value, corrected p-value for all features
- group_comparison_results.csv - Full table with effect sizes and p-values for all features
- combined_results.csv - Combined ranked summary of both analyses
- results_summary.txt - Plain text summary of key findings

---

## Understanding your results

### Spearman correlation (rho)

Ranges from -1 to +1. A positive rho means higher biomarker levels associate with more progression. A negative rho means higher levels associate with less progression or improvement. Values above 0.2 or below -0.2 are generally considered meaningful in clinical biomarker research.

### Rank-biserial correlation (effect size)

Ranges from -1 to +1. Positive values mean the biomarker is higher in progressors. Rough guidelines: 0.1 = small effect, 0.3 = medium effect, 0.5 = large effect.

### Multiple testing correction

With hundreds of biomarkers tested simultaneously, some will appear significant by chance. Benjamini-Hochberg FDR correction controls the expected proportion of false positives among significant results. Features marked as significant after correction are more reliable findings.

### Volcano plot

Each dot is one biomarker. Position on the x-axis shows the correlation direction and strength. Height on the y-axis shows statistical significance. Biomarkers in the upper left or upper right corners are the most interesting — strong association AND statistically significant.

---

## Companion tool

For machine learning-based progression prediction using the same data format, see the companion repository: https://github.com/CameronPiepkorn/neuro-progression-ml

---

## Limitations

Results are observational and hypothesis-generating. Statistical significance does not imply causation. All findings should be validated in independent cohorts before drawing clinical conclusions.

---

## License

MIT License. Free to use and modify for research purposes.

---

## Questions or issues

Please open an issue on GitHub.
