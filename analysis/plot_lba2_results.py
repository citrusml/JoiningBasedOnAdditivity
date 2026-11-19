"""Analysis and visualization script for NELSI LBA2 benchmark results.

Generates publication-quality (300 DPI) multi-panel figures comparing
phylogenetic reconstruction methods (MSA+ML, NJ variants, and Greedy-Split/GS
variants) under Long Branch Attraction (LBA ratio = 3.0 and 5.0) on NELSI trees.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

METHOD_ORDER: List[str] = [
    "MSA_ML",
    "MSA_NJ",
    "Poisson_NJ",
    "SSS_Linear_NJ",
    "SSS_Log_NJ",
    "GS_Poisson",
    "GS_SSS_Linear",
    "GS_SSS_Log",
]

METHOD_LABELS: Dict[str, str] = {
    "MSA_ML": "MSA + ML (IQ-TREE)",
    "MSA_NJ": "MSA + NJ",
    "Poisson_NJ": "Poisson + NJ",
    "SSS_Linear_NJ": "SSS-Linear + NJ",
    "SSS_Log_NJ": "SSS-Log + NJ",
    "GS_Poisson": "GS (Poisson)",
    "GS_SSS_Linear": "GS (SSS-Linear)",
    "GS_SSS_Log": "GS (SSS-Log)",
}

METHOD_STYLES: Dict[str, Dict[str, object]] = {
    "MSA_ML": {"color": "#111111", "marker": "o", "linestyle": "-", "linewidth": 2.2, "zorder": 10},
    "MSA_NJ": {"color": "#1f77b4", "marker": "s", "linestyle": "-", "linewidth": 1.8, "zorder": 8},
    "Poisson_NJ": {"color": "#2ca02c", "marker": "^", "linestyle": "-", "linewidth": 1.8, "zorder": 7},
    "SSS_Linear_NJ": {"color": "#ff7f0e", "marker": "D", "linestyle": "-", "linewidth": 1.8, "zorder": 7},
    "SSS_Log_NJ": {"color": "#d62728", "marker": "v", "linestyle": "-", "linewidth": 1.8, "zorder": 6},
    "GS_Poisson": {"color": "#2ca02c", "marker": "^", "linestyle": "--", "linewidth": 1.7, "zorder": 5},
    "GS_SSS_Linear": {"color": "#ff7f0e", "marker": "D", "linestyle": "--", "linewidth": 1.7, "zorder": 5},
    "GS_SSS_Log": {"color": "#d62728", "marker": "v", "linestyle": "--", "linewidth": 1.7, "zorder": 5},
}


def configure_matplotlib() -> None:
    """Configure global matplotlib typography and layout settings.

    Parameters
    ----------
    None
        Configures ``plt.rcParams`` for 300 DPI publication-quality output.
    """
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.titlesize": 13,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linestyle": ":",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": 300,
        }
    )


def compute_mean_sem(
    df: pd.DataFrame,
    group_cols: Sequence[str],
    metric: str,
    scale: float = 1.0,
) -> pd.DataFrame:
    """Compute mean and standard error of the mean (SEM) for a metric.

    Parameters
    ----------
    df : pd.DataFrame
        Raw benchmark results dataframe.
    group_cols : Sequence[str]
        Grouping columns (e.g., ``['lba_ratio', 'n_taxa', 'divergence', 'method']``).
    metric : str
        Column name of the target metric (unit depends on metric, e.g.,
        ``unitless`` for nRF, ``substitutions/site`` for RMSE, ``seconds`` for runtime).
    scale : float, optional
        Multiplicative scaling factor (e.g., ``100.0`` for percentages), by default ``1.0``.

    Returns
    -------
    pd.DataFrame
        Aggregated dataframe with ``mean`` and ``sem`` columns.
    """
    sub = df.dropna(subset=[metric]).copy()
    sub[metric] = sub[metric] * scale
    agg = (
        sub.groupby(list(group_cols))[metric]
        .agg(mean="mean", sem="sem")
        .reset_index()
    )
    return agg


def _plot_method_curves(
    ax: plt.Axes,
    sub_df: pd.DataFrame,
    x_positions: Dict[float, int],
    x_col: str = "divergence",
    methods: Sequence[str] = METHOD_ORDER,
) -> None:
    """Plot Mean +/- SEM curves for each method on a given Axes.

    Parameters
    ----------
    ax : plt.Axes
        Target matplotlib Axes.
    sub_df : pd.DataFrame
        Aggregated dataframe containing ``x_col``, ``method``, ``mean``, and ``sem``.
    x_positions : Dict[float, int]
        Mapping from divergence values (``substitutions/site``) to discrete x-axis indices.
    x_col : str, optional
        X-axis column name (``divergence`` in ``substitutions/site``), by default ``'divergence'``.
    methods : Sequence[str], optional
        Ordered list of method identifiers to plot.
    """
    for method in methods:
        mdf = sub_df[sub_df["method"] == method].sort_values(x_col)
        if mdf.empty:
            continue
        style = METHOD_STYLES[method]
        xs = np.array([x_positions[float(v)] for v in mdf[x_col].to_numpy()])
        ax.errorbar(
            xs,
            mdf["mean"].to_numpy(),
            yerr=mdf["sem"].to_numpy(),
            label=METHOD_LABELS[method],
            color=str(style["color"]),
            marker=str(style["marker"]),
            linestyle=str(style["linestyle"]),
            linewidth=float(style["linewidth"]),
            markersize=4.5,
            capsize=2.5,
            elinewidth=1.0,
            alpha=0.92,
            zorder=int(style["zorder"]),
        )


def _add_shared_legend(fig: plt.Figure, ncol: int = 4, y_pos: float = 0.01) -> None:
    """Add a shared figure-level legend collecting unique handles across all Axes.

    Parameters
    ----------
    fig : plt.Figure
        Target figure.
    ncol : int, optional
        Number of legend columns, by default ``4``.
    y_pos : float, optional
        Vertical anchor position in figure coordinates, by default ``0.01``.
    """
    seen: Dict[str, object] = {}
    for ax in fig.axes:
        handles, labels = ax.get_legend_handles_labels()
        for h, l in zip(handles, labels):
            if l not in seen:
                seen[l] = h
    ordered_labels = [METHOD_LABELS[m] for m in METHOD_ORDER if METHOD_LABELS[m] in seen]
    for l in seen:
        if l not in ordered_labels:
            ordered_labels.append(l)
    ordered_handles = [seen[l] for l in ordered_labels]
    fig.legend(
        ordered_handles,
        ordered_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, y_pos),
        ncol=ncol,
        frameon=True,
        framealpha=0.95,
        edgecolor="#cccccc",
    )


def _build_x_positions(divergences: Sequence[float]) -> Tuple[Dict[float, int], List[int], List[str]]:
    """Build evenly spaced discrete x-axis positions and tick labels for divergence.

    Parameters
    ----------
    divergences : Sequence[float]
        Sorted divergence values in ``substitutions/site``.

    Returns
    -------
    Tuple[Dict[float, int], List[int], List[str]]
        Mapping from divergence to index, tick positions, and formatted tick labels.
    """
    pos_map = {float(d): i for i, d in enumerate(divergences)}
    ticks = list(range(len(divergences)))
    labels = [f"{d:g}" for d in divergences]
    return pos_map, ticks, labels


def plot_metric_grid_by_lba_and_taxa(
    df: pd.DataFrame,
    metric: str,
    ylabel: str,
    title: str,
    out_path: Path,
    ylim: Optional[Tuple[float, float]] = None,
    hline: Optional[float] = None,
    scale: float = 1.0,
) -> None:
    """Create a 2x3 grid (LBA ratio x n_taxa) of metric vs divergence curves.

    Parameters
    ----------
    df : pd.DataFrame
        Raw benchmark dataframe.
    metric : str
        Target metric column name (units vary: ``unitless`` for nRF/quartet attraction).
    ylabel : str
        Y-axis label including unit description.
    title : str
        Overall figure title.
    out_path : Path
        Destination PNG file path.
    ylim : Optional[Tuple[float, float]], optional
        Explicit Y-axis limits, by default ``None``.
    hline : Optional[float], optional
        Optional horizontal reference line value, by default ``None``.
    scale : float, optional
        Multiplicative scaling factor for the metric, by default ``1.0``.
    """
    agg = compute_mean_sem(df, ["lba_ratio", "n_taxa", "divergence", "method"], metric, scale=scale)
    lba_ratios = sorted(df["lba_ratio"].unique())
    taxa_list = sorted(df["n_taxa"].unique())
    divergences = sorted(df["divergence"].unique())
    x_pos, x_ticks, x_labels = _build_x_positions(divergences)

    fig, axes = plt.subplots(
        len(lba_ratios),
        len(taxa_list),
        figsize=(13.5, 7.6),
        sharex=True,
        sharey=True,
    )

    for r_idx, lba in enumerate(lba_ratios):
        for c_idx, n_taxa in enumerate(taxa_list):
            ax = axes[r_idx, c_idx]
            sub = agg[(agg["lba_ratio"] == lba) & (agg["n_taxa"] == n_taxa)]
            _plot_method_curves(ax, sub, x_positions=x_pos)
            if hline is not None:
                ax.axhline(hline, color="#666666", linestyle=":", linewidth=1.0, zorder=1)
            ax.set_title(f"LBA Ratio = {int(lba)}x, N = {int(n_taxa)} taxa")
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_labels)
            if ylim is not None:
                ax.set_ylim(*ylim)
            if r_idx == len(lba_ratios) - 1:
                ax.set_xlabel("Divergence D (substitutions/site)")
            if c_idx == 0:
                ax.set_ylabel(ylabel)

    fig.suptitle(title, fontweight="bold", y=0.98)
    _add_shared_legend(fig, ncol=4, y_pos=0.01)
    fig.tight_layout(rect=[0.0, 0.08, 1.0, 0.94])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_nta_and_false_sister(df: pd.DataFrame, out_path: Path) -> None:
    """Plot NTA (unitless) and False Sister Rate (%) across LBA ratios and divergence.

    Parameters
    ----------
    df : pd.DataFrame
        Raw benchmark dataframe.
    out_path : Path
        Destination PNG file path.
    """
    agg_nta = compute_mean_sem(df, ["lba_ratio", "divergence", "method"], "nta", scale=1.0)
    agg_sis = compute_mean_sem(df, ["lba_ratio", "divergence", "method"], "is_sister", scale=100.0)
    lba_ratios = sorted(df["lba_ratio"].unique())
    divergences = sorted(df["divergence"].unique())
    x_pos, x_ticks, x_labels = _build_x_positions(divergences)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.8), sharex=True)

    for r_idx, lba in enumerate(lba_ratios):
        ax_nta = axes[r_idx, 0]
        _plot_method_curves(ax_nta, agg_nta[agg_nta["lba_ratio"] == lba], x_positions=x_pos)
        ax_nta.axhline(0.0, color="#666666", linestyle=":", linewidth=1.1, zorder=1)
        ax_nta.set_title(f"LBA Ratio = {int(lba)}x: Normalized Topological Attraction (NTA)")
        ax_nta.set_ylabel("NTA (unitless, >0 = attraction)")
        ax_nta.set_ylim(-0.38, 0.42)

        ax_sis = axes[r_idx, 1]
        _plot_method_curves(ax_sis, agg_sis[agg_sis["lba_ratio"] == lba], x_positions=x_pos)
        ax_sis.set_title(f"LBA Ratio = {int(lba)}x: Long-Branch False Sister Rate")
        ax_sis.set_ylabel("False Sister Rate (%)")
        ax_sis.set_ylim(-0.5, 14.0)

        for ax in (ax_nta, ax_sis):
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_labels)
            if r_idx == len(lba_ratios) - 1:
                ax.set_xlabel("Divergence D (substitutions/site)")

    fig.suptitle(
        "LBA Bias Metrics Across Divergence (Pooled over N = 20, 32, 64; Mean +/- SEM)",
        fontweight="bold",
        y=0.98,
    )
    _add_shared_legend(fig, ncol=4, y_pos=0.01)
    fig.tight_layout(rect=[0.0, 0.08, 1.0, 0.94])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_rmse_comparison(df: pd.DataFrame, out_path: Path) -> None:
    """Plot Matrix RMSE and Patristic RMSE (substitutions/site) across LBA ratios.

    Parameters
    ----------
    df : pd.DataFrame
        Raw benchmark dataframe.
    out_path : Path
        Destination PNG file path.
    """
    agg_mat = compute_mean_sem(df, ["lba_ratio", "divergence", "method"], "matrix_rmse")
    agg_pat = compute_mean_sem(df, ["lba_ratio", "divergence", "method"], "patristic_rmse")
    lba_ratios = sorted(df["lba_ratio"].unique())
    divergences = sorted(df["divergence"].unique())
    x_pos, x_ticks, x_labels = _build_x_positions(divergences)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.8), sharex=True, sharey=True)

    for r_idx, lba in enumerate(lba_ratios):
        ax_mat = axes[r_idx, 0]
        _plot_method_curves(ax_mat, agg_mat[agg_mat["lba_ratio"] == lba], x_positions=x_pos)
        ax_mat.set_title(f"LBA Ratio = {int(lba)}x: Distance Matrix RMSE")
        ax_mat.set_ylabel("Matrix RMSE (substitutions/site)")

        ax_pat = axes[r_idx, 1]
        _plot_method_curves(ax_pat, agg_pat[agg_pat["lba_ratio"] == lba], x_positions=x_pos)
        ax_pat.set_title(f"LBA Ratio = {int(lba)}x: Inferred Tree Patristic RMSE")
        ax_pat.set_ylabel("Patristic RMSE (substitutions/site)")

        for ax in (ax_mat, ax_pat):
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_labels)
            if r_idx == len(lba_ratios) - 1:
                ax.set_xlabel("Divergence D (substitutions/site)")

    fig.suptitle(
        "Distance Matrix vs. Patristic Distance RMSE (substitutions/site; Mean +/- SEM)",
        fontweight="bold",
        y=0.98,
    )
    _add_shared_legend(fig, ncol=4, y_pos=0.01)
    fig.tight_layout(rect=[0.0, 0.08, 1.0, 0.94])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_root_dist_margin_and_runtime(df: pd.DataFrame, out_path: Path) -> None:
    """Plot Matrix LBA Margin, Root-to-Tip Max/Min Ratio, and Runtime (seconds).

    Parameters
    ----------
    df : pd.DataFrame
        Raw benchmark dataframe.
    out_path : Path
        Destination PNG file path.
    """
    agg_margin = compute_mean_sem(df, ["lba_ratio", "divergence", "method"], "matrix_lba_margin")
    agg_ratio = compute_mean_sem(df, ["lba_ratio", "divergence", "method"], "root_dist_max_min_ratio")
    agg_time = compute_mean_sem(df, ["lba_ratio", "divergence", "method"], "runtime_sec")
    true_ratio = (
        df.groupby(["lba_ratio", "divergence"])["true_root_dist_max_min_ratio"]
        .mean()
        .reset_index()
    )
    lba_ratios = sorted(df["lba_ratio"].unique())
    divergences = sorted(df["divergence"].unique())
    x_pos, x_ticks, x_labels = _build_x_positions(divergences)

    fig, axes = plt.subplots(2, 3, figsize=(14.0, 7.8), sharex=True)

    for r_idx, lba in enumerate(lba_ratios):
        ax_m = axes[r_idx, 0]
        _plot_method_curves(ax_m, agg_margin[agg_margin["lba_ratio"] == lba], x_positions=x_pos)
        ax_m.axhline(0.0, color="#666666", linestyle=":", linewidth=1.0, zorder=1)
        ax_m.set_title(f"LBA {int(lba)}x: Matrix LBA Margin")
        ax_m.set_ylabel("Matrix LBA Margin (substitutions/site)")

        ax_r = axes[r_idx, 1]
        _plot_method_curves(ax_r, agg_ratio[agg_ratio["lba_ratio"] == lba], x_positions=x_pos)
        tr_sub = true_ratio[true_ratio["lba_ratio"] == lba].sort_values("divergence")
        tr_xs = np.array([x_pos[float(v)] for v in tr_sub["divergence"].to_numpy()])
        ax_r.plot(
            tr_xs,
            tr_sub["true_root_dist_max_min_ratio"].to_numpy(),
            color="#555555",
            linestyle=":",
            linewidth=2.0,
            label="True Tree Ratio",
            zorder=2,
        )
        ax_r.set_title(f"LBA {int(lba)}x: Root-to-Tip Max/Min Ratio")
        ax_r.set_ylabel("Max/Min Ratio (unitless, log scale)")
        ax_r.set_yscale("log")

        ax_t = axes[r_idx, 2]
        _plot_method_curves(ax_t, agg_time[agg_time["lba_ratio"] == lba], x_positions=x_pos)
        ax_t.set_title(f"LBA {int(lba)}x: Runtime")
        ax_t.set_ylabel("Runtime (seconds, log scale)")
        ax_t.set_yscale("log")

        for ax in (ax_m, ax_r, ax_t):
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_labels)
            if r_idx == len(lba_ratios) - 1:
                ax.set_xlabel("Divergence D (substitutions/site)")

    fig.suptitle(
        "Matrix LBA Margin, Root-to-Tip Ratio, and Runtime (Mean +/- SEM)",
        fontweight="bold",
        y=0.98,
    )
    _add_shared_legend(fig, ncol=5, y_pos=0.01)
    fig.tight_layout(rect=[0.0, 0.08, 1.0, 0.94])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    """Parse CLI arguments, load benchmark CSV, and generate all figures.

    Parameters
    ----------
    None
        Reads ``--csv`` and ``--outdir`` from the command line.
    """
    parser = argparse.ArgumentParser(description="Analyze and plot NELSI LBA2 benchmark results.")
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("results_lba2/benchmark_results.csv"),
        help="Path to benchmark_results.csv",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results_lba2/plots"),
        help="Directory to save output PNG figures",
    )
    args = parser.parse_args()

    configure_matplotlib()
    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.csv)

    plot_metric_grid_by_lba_and_taxa(
        df=df,
        metric="nrf",
        ylabel="Normalized RF Distance (unitless)",
        title="Overall Topology Error (nRF) on NELSI LBA Trees (L=500 aa, Indel=0.1, 50 reps)",
        out_path=args.outdir / "lba2_nrf_comparison.png",
        ylim=(0.0, 0.65),
    )

    plot_metric_grid_by_lba_and_taxa(
        df=df,
        metric="lba_quartet_attraction",
        ylabel="LBA Quartet Attraction (unitless)",
        title="LBA Quartet Attraction Rate Across Methods (Mean +/- SEM)",
        out_path=args.outdir / "lba2_quartet_attraction.png",
        ylim=(-0.02, 0.50),
        hline=0.0,
    )

    plot_metric_grid_by_lba_and_taxa(
        df=df,
        metric="pair_nrf",
        ylabel="Pair-Separating Split nRF (unitless)",
        title="LBA Pair-Separating Split Error (pair_nrf) Across Methods (Mean +/- SEM)",
        out_path=args.outdir / "lba2_pair_nrf.png",
        ylim=(0.0, 0.95),
    )

    plot_nta_and_false_sister(df, args.outdir / "lba2_nta_and_false_sister.png")
    plot_rmse_comparison(df, args.outdir / "lba2_rmse_comparison.png")
    plot_root_dist_margin_and_runtime(df, args.outdir / "lba2_root_dist_and_margin.png")
    print(f"Saved 6 figures to {args.outdir.resolve()}")


if __name__ == "__main__":
    main()
