"""
Plot Bayesian model selection results for driving inputs (C-matrix).

This script:
  - Loads model- and region-level posterior probabilities from CSV files
  - Renames models/regions for clearer plotting
  - Produces a bar plot of model posterior probabilities
  - Saves the figure as a 300 dpi PNG

Adjust the PATH CONFIGURATION block below to match your environment.
"""

import pandas as pd
import numpy as np  # kept for consistency with original script
import matplotlib.pyplot as plt

# -------------------------------------------------------------------------
# PATH CONFIGURATION
# -------------------------------------------------------------------------
SPEC = "GCM_spec"  # main run spec
BASE_DIR = (
    "/.../peb_results"
)
ANALYSIS_SUBDIR = "results_subdir"

path = f"{BASE_DIR}/{ANALYSIS_SUBDIR}/{SPEC}"
models_csv = f"{path}/driving_inputs_models.csv"
regions_csv = f"{path}/driving_inputs_regions.csv"  # loaded for completeness

# -------------------------------------------------------------------------
# LOAD AND PREPARE DATA
# -------------------------------------------------------------------------
df_models = pd.read_csv(models_csv)
df_regions = pd.read_csv(regions_csv)

# Human-readable labels for models / regions
name_dict = {
    "hipp": "Hipp",
    "postvmpfc": "p-vmPFC",
    "amygdala": "Amygdala",
    "hipp+postvmpfc": "Hipp + p-vmPFC",
    "hipp+amygdala": "Hipp + Amygdala",
    "postvmpfc+amygdala": "p-vmPFC + Amygdala",
    "all": "Hipp + p-vmPFC + Amygdala",
}

# Move the first row ("all" model) to the end for plotting order
df_models = pd.concat([df_models.iloc[1:], df_models.iloc[[0]]], ignore_index=True)

# Apply readable labels
df_models["region"] = df_models["region"].replace(name_dict)

# Round probabilities to 2 decimal places
df_models["model_probability"] = df_models["model_probability"].round(2)

# -------------------------------------------------------------------------
# PLOT: MODEL POSTERIOR PROBABILITIES
# -------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 7))

bars = ax.bar(df_models["region"], df_models["model_probability"], color="#E29338")

# Add value labels above each bar
for bar in bars:
    h = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        h + 0.02,  # small vertical offset
        f"{h:.2f}",
        ha="center",
        va="bottom",
        fontsize=10,
    )

# Y-axis limits and labels
ax.set_ylim(0, df_models["model_probability"].max() + 0.1)
ax.set_ylabel("Posterior Probability", fontsize=12)

ax.set_title(
    "Driving Input Bayesian Model Selection",
    fontsize=16,
    pad=25,
)

# Clean up axes
for side in ("top", "right"):
    ax.spines[side].set_visible(False)

plt.xticks(rotation=45, ha="right")

plt.tight_layout()
plt.savefig(f"{path}/models_bms_plot.png", dpi=300, bbox_inches="tight")
plt.show()