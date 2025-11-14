"""
LOO cross-validation plots for DCM B-matrix connections (fear task).

This script:
  - Loads LOO results (qE) and the group design matrix (M)
  - Computes correlation between predicted effects and a STAT-T covariate
  - Performs SPM-style one-tailed testing
  - Plots group-wise violin + strip plots
  - Saves each figure as a 300 dpi PNG

Adjust RESULTS_DIR and LOO_FILES to match your environment.
"""

import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy import stats
from scipy.sparse import issparse
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# -------------------------------------------------------------------------
# PATHS AND FILE LIST
# -------------------------------------------------------------------------
RESULTS_DIR = Path(
    ".../peb_results"
    "/results_subdir/GCM_spec"
)

LOO_FILES = [
    "LOO_hipp_to_postvmpfc_CS+_earlyconditioning.mat",
    "LOO_amygdala_to_postvmpfc_CS+_earlyconditioning.mat",
    "LOO_hipp_to_postvmpfc_CS-_lateextinction.mat",
    "LOO_amygdala_to_postvmpfc_CS-_lateextinction.mat",
]
LOO_FILEPATHS = [RESULTS_DIR / f for f in LOO_FILES]

# Mapping from short codes to nicer names in titles
NAME_DICT = {
    "hipp": "Hipp",
    "postvmpfc": "p-vmPFC",
    "amygdala": "Amygdala",
    "earlyconditioning": "Early conditioning",
    "earlyextinction": "Early extinction",
    "lateconditioning": "Late conditioning",
    "lateextinction": "Late extinction",
}


# -------------------------------------------------------------------------
# SMALL HELPERS
# -------------------------------------------------------------------------
def _to_dense_1d(x) -> np.ndarray:
    """Convert sparse/array-like to dense 1D float array."""
    if issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=float).squeeze().ravel()


def _to_dense_2d(x) -> np.ndarray:
    """Convert sparse/array-like to dense 2D float array."""
    if issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=float)


def _nice_name(fp: Path) -> str:
    """Generate a readable label from a LOO filename."""
    name = fp.stem[4:] if fp.stem.startswith("LOO_") else fp.stem
    for old, new in sorted(NAME_DICT.items(), key=lambda kv: -len(kv[0])):
        name = name.replace(old, new)
    return name.replace("_to_", " → ").replace("_", " ")


def _p_from_r(r: float, n: int) -> float:
    """SPM-style one-tailed p-value from Pearson r and sample size n."""
    df = n - 2
    r = np.clip(r, -0.999999999, 0.999999999)
    t = r * np.sqrt(df / (1 - r**2))
    return stats.t.sf(t, df)  # one-tailed


# -------------------------------------------------------------------------
# MAIN LOOP: LOO PLOTS
# -------------------------------------------------------------------------
if __name__ == "__main__":
    # Plot style
    sns.set_theme(style="whitegrid", context="paper", font_scale=1)
    pastels = sns.color_palette("Set2", 2)

    for filepath in LOO_FILEPATHS:
        # Load LOO output and design matrix
        mat = loadmat(filepath, squeeze_me=True, struct_as_record=False)
        m_mat = loadmat(RESULTS_DIR / "M.mat", squeeze_me=True, struct_as_record=False)

        qE = _to_dense_1d(mat["qE"])          # predicted effects
        X = _to_dense_2d(m_mat["M"].X)        # design matrix
        y = _to_dense_1d(X[:, 1])             # STAT-T covariate (second column)

        # Correlation (SPM-style, mean-centred and normed)
        def norm(z: np.ndarray) -> np.ndarray:
            z0 = z - z.mean()
            return z0 / np.linalg.norm(z0)

        r, _ = stats.pearsonr(norm(y), norm(qE))
        p_one = _p_from_r(r, qE.size)

        # Two-group split based on STAT-T column
        vals = np.unique(y)
        if vals.size != 2:
            raise ValueError("Design column must have exactly two groups.")
        groups = np.where(y == vals.min(), "Lower STAT-T", "Higher STAT-T")

        df_plot = pd.DataFrame({"STAT-T group": groups, "qE": qE})

        # Violin + strip plot
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        sns.violinplot(
            data=df_plot,
            x="STAT-T group",
            y="qE",
            order=["Lower STAT-T", "Higher STAT-T"],
            cut=0,
            inner=None,
            linewidth=1.0,
            palette=pastels,
            saturation=0.6,
            ax=ax,
        )
        sns.stripplot(
            data=df_plot,
            x="STAT-T group",
            y="qE",
            order=["Lower STAT-T", "Higher STAT-T"],
            color="k",
            size=4,
            alpha=0.45,
            jitter=0.06,
            ax=ax,
        )

        ax.set_ylim(-3, 4.5)  # consistent across plots
        ax.set_xlabel("")
        ax.set_ylabel("Predicted effect (qE)")
        ax.set_title(
            f"{_nice_name(filepath)}\n"
            f"LOO prediction · r={r:.2f} · p(one-tailed)={p_one:.3g} · N={qE.size}",
            pad=8,
        )
        sns.despine()
        fig.tight_layout()

        fname = f"{_nice_name(filepath)}_LOO_violin.png"
        fig.savefig(
            RESULTS_DIR / fname,
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.6,
        )
        plt.show()
        plt.close(fig)