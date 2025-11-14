"""
Plot intrinsic and modulatory DCM connections from group DCM / PEB results.

This script:
  - Reads connection summaries from CSV files
  - Draws directed graphs for A- and B-matrix parameters
  - Saves figures as 300 dpi PNGs in a subfolder called 'TNR'

Adjust the PATH CONFIGURATION block below to match your environment.
"""

import itertools as it
from typing import Dict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import pandas as pd
import numpy as np

# -------------------------------------------------------------------------
# Matplotlib defaults
# -------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
})

# -------------------------------------------------------------------------
# Node appearance and naming
# -------------------------------------------------------------------------
NODE_COLOURS: Dict[str, str] = {
    "Hipp":     "#b1e6c2",
    "p-vmPFC":  "#d8c3ff",
    "Amygdala": "#fff6b3",
}

name_dict = {
    "hipp":      "Hipp",
    "postvmpfc": "p-vmPFC",
    "amygdala":  "Amygdala",
}

anchor_xy = {
    "Hipp":     (-0.5, -0.5),   # left
    "p-vmPFC":  (0.5,  -0.5),   # right
    "Amygdala": (0.0,   0.5),   # top
}


# -------------------------------------------------------------------------
# Graph plotting utility
# -------------------------------------------------------------------------
def plot_ep_graph(
    df: pd.DataFrame,
    path: str | Path,
    shell_layout: bool = True,
    scale_labels: float = 1.25,
    width_scale: float = 20.0,
    node_size: int = 500,
    arrowsize: int = 25,
    title: str | None = None,
    pad: float = 0.25,
    *,
    node_colours: Dict[str, str] = NODE_COLOURS,
    fixed_pos: Dict[str, tuple[float, float]] | None = anchor_xy,
    ep_attr: str = "Ep",
    pp_attr: str = "Pp",
    source_col: str = "source",
    target_col: str = "target",
    legend: str = "main",
) -> None:
    """
    Plot a directed connectivity graph from a DataFrame of edges.

    Parameters
    ----------
    df : DataFrame
        Must contain columns [source_col, target_col, ep_attr, pp_attr].
    path : str or Path
        Base directory where a 'TNR' subfolder and PNGs will be saved.
    title : str, optional
        Figure title and base filename ('.png' is added automatically).
    legend : {'main', ''}, optional
        Legend labels: 'main' for excitatory/inhibitory, else increased/decreased.
    """

    # Build directed multigraph from DataFrame
    G = nx.from_pandas_edgelist(
        df,
        source=source_col,
        target=target_col,
        edge_attr=[ep_attr, pp_attr],
        create_using=nx.MultiDiGraph()
    )

    # Node positions
    if fixed_pos:
        G.add_nodes_from(fixed_pos)          # ensure all anchored nodes exist
        pos = dict(fixed_pos)
    else:
        pos = nx.shell_layout(G) if shell_layout else nx.spring_layout(G)

    pos_lbl = {n: (x * scale_labels, y * scale_labels) for n, (x, y) in pos.items()}

    # Nodes
    nx.draw_networkx_nodes(
        G, pos,
        node_size=node_size,
        node_color=[node_colours.get(n, "#d3d3d3") for n in G.nodes()]
    )
    nx.draw_networkx_labels(G, pos_lbl, font_size=20, font_family="Times New Roman")

    # Edge colours & widths (red = Ep ≥ 0, blue = Ep < 0)
    edges = G.edges(keys=True, data=True)
    colours = ['#ff9999' if d[ep_attr] >= 0 else '#9999ff' for *_, d in edges]
    widths  = [abs(d[ep_attr]) * width_scale for *_, d in edges]

    nx.draw_networkx_edges(
        G, pos,
        edge_color=colours,
        width=widths,
        connectionstyle="arc3,rad=0.15",
        arrowsize=arrowsize
    )

    # Edge labels: Ep value and rounded Pp (capped at 0.99)
    connectionstyle = [f"arc3,rad={r}" for r in it.accumulate([0.15] * 4)]
    labels = {
        tuple(e): (
            f"{ep_attr}={d[ep_attr]}, "
            f"{pp_attr}>{min(np.floor(d[pp_attr] * 100) / 100, 0.99)}"
        )
        for *e, d in edges
    }
    nx.draw_networkx_edge_labels(
        G, pos,
        labels,
        connectionstyle=connectionstyle,
        label_pos=0.5,
        font_color="black"
    )

    # Layout
    plt.margins(pad)
    plt.axis("off")
    plt.gca().set_frame_on(False)

    # Legend
    if legend == "main":
        legend_elements = [
            Line2D([0], [0], color="#ff9999", lw=3, label="Excitatory (Ep ≥ 0)"),
            Line2D([0], [0], color="#9999ff", lw=3, label="Inhibitory (Ep < 0)")
        ]
    else:
        legend_elements = [
            Line2D([0], [0], color="#ff9999", lw=3, label="Increased Ep"),
            Line2D([0], [0], color="#9999ff", lw=3, label="Decreased Ep")
        ]

    plt.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.95, 0.8),
        frameon=False,
        fontsize=10
    )

    if title:
        plt.title(title, fontsize=16)

    plt.gcf().set_constrained_layout(True)

    # Save figure
    fig = plt.gcf()
    out_dir = Path(path, "TNR")
    out_dir.mkdir(parents=True, exist_ok=True)
    fname = f"{title or 'plot'}.png"
    fig.savefig(
        out_dir / fname,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.25,
    )

    plt.show()
    plt.close(fig)


# -------------------------------------------------------------------------
# Path configuration and plotting pipeline
# -------------------------------------------------------------------------
if __name__ == "__main__":
    # Base configuration
    SPEC_NAME = "GCM_spec"
    ROOT_RESULTS = Path(
        ".../peb_results"
    )
    ANALYSIS_SUBDIR = "results_subdir"

    base_path = ROOT_RESULTS / ANALYSIS_SUBDIR / SPEC_NAME
    connections_csv = base_path / "connections.csv"
    class_connections_csv = base_path / "class_connections.csv"

    # Load main connection table
    df = pd.read_csv(connections_csv)
    df["source"] = df["source"].replace(name_dict)
    df["target"] = df["target"].replace(name_dict)
    df["Ep"] = df["Ep"].round(2)

    df_A = df[df["matrix"] == "A"][["source", "target", "Ep", "Pp"]].reset_index(drop=True)
    df_B = df[df["matrix"] == "B"][["condition", "source", "target", "Ep", "Pp"]].reset_index(drop=True)

    # Plot A-matrix (intrinsic connections)
    plot_ep_graph(
        df_A,
        path=base_path,
        title="A Matrix: Intrinsic connections",
        width_scale=5,
    )

    # Condition labels for B-matrix plots
    condition_dict = {
        1.0: "Task",
        2.0: "Early conditioning CS+",
        3.0: "Early conditioning CS-",
        4.0: "Early extinction CS+",
        5.0: "Early extinction CS-",
        6.0: "Late conditioning CS+",
        7.0: "Late conditioning CS-",
        8.0: "Late extinction CS+",
        9.0: "Late extinction CS-",
    }

    # Plot B-matrix (condition-specific modulation)
    for con in [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]:
        df_con = df_B[df_B["condition"] == con]
        df_con = df_con[df_con["Pp"] > 0.95].reset_index(drop=True)
        plot_ep_graph(
            df_con,
            path=base_path,
            title=f"B Matrix: {condition_dict[con]}",
            width_scale=3.5,
        )

    # ---------------------------------------------------------------------
    # Class-specific effects (e.g. predictive of trait anxiety)
    # ---------------------------------------------------------------------
    df_class = pd.read_csv(class_connections_csv)
    df_class["source"] = df_class["source"].replace(name_dict)
    df_class["target"] = df_class["target"].replace(name_dict)
    df_class["Ep"] = df_class["Ep"].round(2)
    df_class = df_class[["condition", "source", "target", "Ep", "Pp"]].reset_index(drop=True)

    for con in [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]:
        df_class_con = df_class[df_class["condition"] == con]
        df_class_con = df_class_con[df_class_con["Pp"] > 0.95].reset_index(drop=True)
        if df_class_con.shape[0] > 0:
            plot_ep_graph(
                df_class_con,
                path=base_path,
                title=f"{condition_dict[con]}: Predictive of Trait Anxiety",
                width_scale=5,
                legend="",
            )