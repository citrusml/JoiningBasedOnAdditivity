"""
Publication-Quality Visualization Script for Large-Scale Phylometrics Benchmark.

Incorporates:
- Taxa = 20, 32, 64, 128
- Standard Deviation (SD) error bars
- Independent y-axis scaling per divergence level for maximum clarity
- Disconnected plots (no connecting lines between taxa) with horizontal dodging
- Displays all evaluated reconstruction methods
"""

from typing import List, Dict, Tuple, Optional
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_benchmark_results(
    df: pd.DataFrame,
    outdir: str = "results_large_scale",
    filename_prefix: str = "nrf_vs_taxa",
) -> None:
    """
    Generates high-resolution publication charts according to user specifications.

    Parameters
    ----------
    df : pd.DataFrame
        Benchmark records with columns like 'n_taxa', 'divergence', and '<method>_nrf'.
    outdir : str, default='results_large_scale'
        Output directory path.
    filename_prefix : str, default='nrf_vs_taxa'
        Prefix for saved figure files.
    """
    os.makedirs(outdir, exist_ok=True)

    taxa_list = sorted(df["n_taxa"].unique())
    divergences = sorted(df["divergence"].unique())

    # Complete list of methods to plot
    method_specs = [
        ("MSA_ML", "MSA + ML (FastTree)", "#d62728", "o"),
        ("Standard_NJ", "Standard NJ", "#1f77b4", "s"),
        ("BioNJ", "BioNJ", "#ff7f0e", "^"),
        ("FastME_SPR", "FastME-SPR", "#2ca02c", "D"),
        ("AddNJ2", "AddNJ2", "#9467bd", "P"),
        ("BioNJ_Add_NNI", "BioNJ + Add-NNI", "#8c564b", "v"),
        ("AddNJ2_Add_NNI", "AddNJ2 + Add-NNI", "#e377c2", "X"),
    ]

    # Filter to methods actually present in the dataframe
    methods = [spec for spec in method_specs if f"{spec[0]}_nrf" in df.columns]

    n_methods = len(methods)
    # Horizontal dodge offsets within each taxa position to prevent marker overlap
    dodge_offsets = np.linspace(-0.28, 0.28, n_methods)

    # 1. nRF Plot: 1 row, 3 columns (one per divergence)
    fig, axes = plt.subplots(1, len(divergences), figsize=(18, 5.5), dpi=300, sharey=False)
    if len(divergences) == 1:
        axes = [axes]

    for ax, div in zip(axes, divergences):
        sub_df = df[df["divergence"] == div]
        max_y_val = 0.05

        for idx, (m_col, m_label, m_color, m_marker) in enumerate(methods):
            x_pts = []
            y_means = []
            y_sds = []

            for t_idx, n in enumerate(taxa_list):
                taxa_df = sub_df[sub_df["n_taxa"] == n]
                if f"{m_col}_nrf" in taxa_df and taxa_df[f"{m_col}_nrf"].notna().any():
                    series = taxa_df[f"{m_col}_nrf"].dropna()
                    if len(series) > 0:
                        x_pts.append(t_idx + dodge_offsets[idx])
                        mean_val = series.mean()
                        sd_val = series.std(ddof=1) if len(series) > 1 else 0.0
                        y_means.append(mean_val)
                        y_sds.append(sd_val)
                        if mean_val + sd_val > max_y_val:
                            max_y_val = mean_val + sd_val

            if x_pts:
                # Plot with error bars, NO connecting lines
                ax.errorbar(
                    x_pts,
                    y_means,
                    yerr=y_sds,
                    fmt=m_marker,
                    color=m_color,
                    label=m_label,
                    markersize=7,
                    capsize=3.5,
                    capthick=1.2,
                    elinewidth=1.2,
                    linestyle="none",
                    alpha=0.9,
                    markeredgecolor="black",
                    markeredgewidth=0.5,
                )

        ax.set_title(f"Divergence Scale D = {div}", fontsize=13, fontweight="bold", pad=10)
        ax.set_xlabel("Number of Taxa (N)", fontsize=11, fontweight="bold")
        ax.set_xticks(range(len(taxa_list)))
        ax.set_xticklabels([f"N = {n}" for n in taxa_list], fontsize=10, fontweight="bold")
        
        # Set individual y-limits per divergence scale for optimal visibility
        y_margin = max_y_val * 0.12
        ax.set_ylim(-max(0.005, y_margin * 0.3), max_y_val + y_margin)
        ax.set_ylabel("Normalized Robinson-Foulds Distance (nRF)", fontsize=11, fontweight="bold")

        # Subtle vertical divider between taxa groupings
        for sep in np.arange(0.5, len(taxa_list) - 0.5, 1.0):
            ax.axvline(sep, color="gray", linestyle=":", alpha=0.35, zorder=0)

        ax.grid(True, axis="y", ls="--", alpha=0.5, zorder=0)

    # Put legend outside on the right of the last axis
    axes[-1].legend(
        bbox_to_anchor=(1.04, 1.0),
        loc="upper left",
        fontsize=10,
        frameon=True,
        edgecolor="gray",
        title="Reconstruction Method",
        title_fontsize=11,
    )

    plt.tight_layout()
    fig_path = os.path.join(outdir, f"{filename_prefix}.png")
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved figure: {fig_path}")


def plot_runtime_results(
    df: pd.DataFrame,
    outdir: str = "results_large_scale",
    filename_prefix: str = "runtime_vs_taxa",
) -> None:
    """
    Generates computational runtime scaling chart with SD error bars.
    """
    os.makedirs(outdir, exist_ok=True)
    taxa_list = sorted(df["n_taxa"].unique())

    method_specs = [
        ("MSA_ML", "MSA + ML (FastTree)", "#d62728", "o"),
        ("Standard_NJ", "Standard NJ", "#1f77b4", "s"),
        ("BioNJ", "BioNJ", "#ff7f0e", "^"),
        ("FastME_SPR", "FastME-SPR", "#2ca02c", "D"),
        ("AddNJ2", "AddNJ2 (Proposed)", "#9467bd", "P"),
        ("BioNJ_Add_NNI", "BioNJ + Add-NNI", "#8c564b", "v"),
        ("AddNJ2_Add_NNI", "AddNJ2 + Add-NNI", "#e377c2", "X"),
    ]
    methods = [spec for spec in method_specs if f"{spec[0]}_time" in df.columns]

    n_methods = len(methods)
    dodge_offsets = np.linspace(-0.28, 0.28, n_methods)

    plt.figure(figsize=(10, 5.5), dpi=300)

    for idx, (m_col, m_label, m_color, m_marker) in enumerate(methods):
        x_pts = []
        time_means = []
        time_sds = []

        for t_idx, n in enumerate(taxa_list):
            taxa_df = df[df["n_taxa"] == n]
            if f"{m_col}_time" in taxa_df and taxa_df[f"{m_col}_time"].notna().any():
                series = taxa_df[f"{m_col}_time"].dropna()
                if len(series) > 0:
                    x_pts.append(t_idx + dodge_offsets[idx])
                    mean_val = series.mean()
                    sd_val = series.std(ddof=1) if len(series) > 1 else 0.0
                    time_means.append(mean_val)
                    time_sds.append(sd_val)

        if x_pts:
            plt.errorbar(
                x_pts,
                time_means,
                yerr=time_sds,
                fmt=m_marker,
                color=m_color,
                label=m_label,
                markersize=7,
                capsize=3.5,
                capthick=1.2,
                elinewidth=1.2,
                linestyle="none",
                alpha=0.9,
                markeredgecolor="black",
                markeredgewidth=0.5,
            )

    plt.yscale("log")
    plt.xticks(range(len(taxa_list)), [f"N = {n}" for n in taxa_list], fontsize=11, fontweight="bold")
    plt.xlabel("Number of Taxa (N)", fontsize=12, fontweight="bold")
    plt.ylabel("Execution Time per Replicate (seconds, log scale)", fontsize=12, fontweight="bold")
    plt.title("Computational Runtime Scaling (L=500 aa)", fontsize=13, fontweight="bold")

    for sep in np.arange(0.5, len(taxa_list) - 0.5, 1.0):
        plt.axvline(sep, color="gray", linestyle=":", alpha=0.35, zorder=0)

    plt.grid(True, which="both", axis="y", ls="--", alpha=0.5, zorder=0)
    plt.legend(bbox_to_anchor=(1.04, 1.0), loc="upper left", fontsize=10, frameon=True, edgecolor="gray")
    plt.tight_layout()

    fig_time = os.path.join(outdir, f"{filename_prefix}.png")
    plt.savefig(fig_time, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved figure: {fig_time}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--large-csv", default="results_large_scale/large_scale_summary.csv")
    parser.add_argument("--taxa20-csv", default=None)
    parser.add_argument("--outdir", default="results_large_scale")
    args = parser.parse_args()

    df_large = pd.read_csv(args.large_csv)
    if args.taxa20_csv and os.path.exists(args.taxa20_csv):
        df_20 = pd.read_csv(args.taxa20_csv)
        df_combined = pd.concat([df_20, df_large], ignore_index=True)
    else:
        df_combined = df_large

    plot_benchmark_results(df_combined, outdir=args.outdir)
    plot_runtime_results(df_combined, outdir=args.outdir)
