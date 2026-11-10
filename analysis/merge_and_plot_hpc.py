"""
Aggregation and Visualization Script for HPC (LSF) Benchmark Results.

Merges all condition shard CSV files from `results_hpc/shards/summary_*.csv` into
`results_hpc/hpc_benchmark_summary.csv` and generates 300 DPI publication-quality
charts following all user visualization specifications:
- Includes Taxa N = 20, 32, 64, 128
- Standard Deviation (SD) error bars
- Independent y-axis scaling per divergence level (sharey=False)
- Disconnected plots (linestyle='none') with horizontal dodge
- Displays all evaluated methods across BLOSUM62, WAG+Gamma, and True distance matrices
- Plots 4-point condition violations (quartet additivity penalty) and SSS relationships
"""

from typing import List, Tuple, Optional
import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def merge_shards(outdir: str = "results_hpc") -> pd.DataFrame:
    """
    Merges all shard CSV files in `<outdir>/shards/summary_*.csv` into a single DataFrame.

    Parameters
    ----------
    outdir : str, default='results_hpc'
        Directory containing `shards/` subdirectory.

    Returns
    -------
    pd.DataFrame
        Combined and sorted benchmark summary DataFrame.
    """
    pattern = os.path.join(outdir, "shards", "summary_*.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        merged_csv = os.path.join(outdir, "hpc_benchmark_summary.csv")
        if os.path.exists(merged_csv):
            return pd.read_csv(merged_csv)
        raise FileNotFoundError(f"No shard CSV files found matching {pattern}")

    dfs = [pd.read_csv(f) for f in files]
    df = pd.concat(dfs, ignore_index=True)
    df = df.drop_duplicates(subset=["n_taxa", "divergence", "replicate"], keep="last")
    df = df.sort_values(["n_taxa", "divergence", "replicate"]).reset_index(drop=True)

    merged_path = os.path.join(outdir, "hpc_benchmark_summary.csv")
    df.to_csv(merged_path, index=False)
    print(f"[OK] Merged {len(files)} shard files ({len(df)} total rows) -> {merged_path}")
    return df


def _plot_metric_by_taxa_and_divergence(
    df: pd.DataFrame,
    col_specs: List[Tuple[str, str, str, str]],
    ylabel: str,
    suptitle: str,
    out_path: str,
) -> None:
    """
    Helper to plot mean +/- SD across Taxa N (dodged, no connecting lines)
    with independent y-axis scaling per Divergence D.
    """
    taxa_list = sorted(df["n_taxa"].unique())
    divergences = sorted(df["divergence"].unique())

    valid_specs = [spec for spec in col_specs if spec[0] in df.columns]
    if not valid_specs:
        return

    n_methods = len(valid_specs)
    dodge_offsets = np.linspace(-0.30, 0.30, n_methods)

    fig, axes = plt.subplots(1, len(divergences), figsize=(18, 5.5), dpi=300, sharey=False)
    if len(divergences) == 1:
        axes = [axes]

    for ax, div in zip(axes, divergences):
        sub_df = df[df["divergence"] == div]
        max_y_val = 0.02

        for idx, (col_name, label, color, marker) in enumerate(valid_specs):
            x_pts = []
            y_means = []
            y_sds = []

            for t_idx, n in enumerate(taxa_list):
                t_df = sub_df[sub_df["n_taxa"] == n]
                if col_name in t_df and t_df[col_name].notna().any():
                    series = t_df[col_name].dropna()
                    if len(series) > 0:
                        x_pts.append(t_idx + dodge_offsets[idx])
                        m_val = float(series.mean())
                        s_val = float(series.std(ddof=1)) if len(series) > 1 else 0.0
                        y_means.append(m_val)
                        y_sds.append(s_val)
                        if m_val + s_val > max_y_val:
                            max_y_val = m_val + s_val

            if x_pts:
                ax.errorbar(
                    x_pts,
                    y_means,
                    yerr=y_sds,
                    fmt=marker,
                    color=color,
                    label=label,
                    markersize=7,
                    capsize=3.5,
                    capthick=1.2,
                    elinewidth=1.2,
                    linestyle="none",
                    alpha=0.9,
                    markeredgecolor="black",
                    markeredgewidth=0.5,
                )

        # Include mean SSS in subplot title if available
        sss_str = ""
        if "sss_mean" in sub_df.columns and sub_df["sss_mean"].notna().any():
            sss_str = f" (Mean SSS = {sub_df['sss_mean'].mean():.3f})"

        ax.set_title(f"Divergence D = {div}{sss_str}", fontsize=12.5, fontweight="bold", pad=10)
        ax.set_xlabel("Number of Taxa (N)", fontsize=11, fontweight="bold")
        ax.set_xticks(range(len(taxa_list)))
        ax.set_xticklabels([f"N = {n}" for n in taxa_list], fontsize=10, fontweight="bold")

        y_margin = max_y_val * 0.14
        ax.set_ylim(-max(0.003, y_margin * 0.25), max_y_val + y_margin)
        ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")

        for sep in np.arange(0.5, len(taxa_list) - 0.5, 1.0):
            ax.axvline(sep, color="gray", linestyle=":", alpha=0.35, zorder=0)

        ax.grid(True, axis="y", ls="--", alpha=0.5, zorder=0)

    axes[-1].legend(
        bbox_to_anchor=(1.03, 1.0),
        loc="upper left",
        fontsize=9.5,
        frameon=True,
        edgecolor="gray",
        title="Method",
        title_fontsize=10.5,
    )
    fig.suptitle(suptitle, fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved figure: {out_path}")


def generate_all_hpc_plots(df: pd.DataFrame, outdir: str = "results_hpc") -> None:
    """
    Generates all nRF, 4-point condition violation, and SSS plots from the merged HPC dataset.
    """
    os.makedirs(outdir, exist_ok=True)

    matrix_labels = {
        "blosum62": "BLOSUM62 + Poisson Distance",
        "wag_gamma": "WAG + Gamma ML Distance",
        "true": "True Patristic Distance",
    }

    # 1. nRF plots for each of the 3 distance matrix types
    for m_type, m_title in matrix_labels.items():
        nrf_specs = [
            ("MSA_ML_nrf", "MSA + ML (FastTree)", "#d62728", "o"),
            (f"{m_type}_Standard_NJ_nrf", "Standard NJ", "#1f77b4", "s"),
            (f"{m_type}_BioNJ_nrf", "BioNJ", "#ff7f0e", "^"),
            (f"{m_type}_FastME_SPR_nrf", "FastME-SPR", "#2ca02c", "D"),
            (f"{m_type}_AddNJ2_nrf", "AddNJ2 (Proposed)", "#9467bd", "P"),
            (f"{m_type}_BioNJ_Add_NNI_nrf", "BioNJ + Add-NNI", "#8c564b", "v"),
            (f"{m_type}_AddNJ2_Add_NNI_nrf", "AddNJ2 + Add-NNI", "#e377c2", "X"),
        ]
        _plot_metric_by_taxa_and_divergence(
            df=df,
            col_specs=nrf_specs,
            ylabel="Normalized Robinson-Foulds Distance (nRF)",
            suptitle=f"Topological Error (nRF ± SD) — {m_title}",
            out_path=os.path.join(outdir, f"nrf_{m_type}_vs_taxa.png"),
        )

    # 2. Direct BLOSUM62 vs WAG+Gamma comparison on key methods
    comp_specs = [
        ("MSA_ML_nrf", "MSA + ML (FastTree)", "#d62728", "o"),
        ("blosum62_Standard_NJ_nrf", "BLOSUM62: Standard NJ", "#1f77b4", "s"),
        ("wag_gamma_Standard_NJ_nrf", "WAG+Γ: Standard NJ", "#aec7e8", "s"),
        ("blosum62_BioNJ_nrf", "BLOSUM62: BioNJ", "#ff7f0e", "^"),
        ("wag_gamma_BioNJ_nrf", "WAG+Γ: BioNJ", "#ffbb78", "^"),
        ("blosum62_AddNJ2_nrf", "BLOSUM62: AddNJ2", "#9467bd", "P"),
        ("wag_gamma_AddNJ2_nrf", "WAG+Γ: AddNJ2", "#c5b0d5", "P"),
        ("blosum62_BioNJ_Add_NNI_nrf", "BLOSUM62: BioNJ+Add-NNI", "#8c564b", "v"),
        ("wag_gamma_BioNJ_Add_NNI_nrf", "WAG+Γ: BioNJ+Add-NNI", "#c49c94", "v"),
    ]
    _plot_metric_by_taxa_and_divergence(
        df=df,
        col_specs=comp_specs,
        ylabel="Normalized Robinson-Foulds Distance (nRF)",
        suptitle="BLOSUM62+Poisson vs WAG+Γ ML Distance Matrix Comparison (nRF ± SD)",
        out_path=os.path.join(outdir, "nrf_blosum62_vs_wag_gamma.png"),
    )

    # 3. 4-Point Condition Violation (Quartet Additivity Penalty) plots
    for m_type, m_title in matrix_labels.items():
        pen_specs = [
            (f"true_tree_penalty_{m_type}", "True Tree", "#222222", "*"),
            (f"MSA_ML_penalty_{m_type}", "MSA + ML (FastTree)", "#d62728", "o"),
            (f"{m_type}_Standard_NJ_penalty_obs", "Standard NJ", "#1f77b4", "s"),
            (f"{m_type}_BioNJ_penalty_obs", "BioNJ", "#ff7f0e", "^"),
            (f"{m_type}_FastME_SPR_penalty_obs", "FastME-SPR", "#2ca02c", "D"),
            (f"{m_type}_AddNJ2_penalty_obs", "AddNJ2 (Proposed)", "#9467bd", "P"),
            (f"{m_type}_BioNJ_Add_NNI_penalty_obs", "BioNJ + Add-NNI", "#8c564b", "v"),
            (f"{m_type}_AddNJ2_Add_NNI_penalty_obs", "AddNJ2 + Add-NNI", "#e377c2", "X"),
        ]
        _plot_metric_by_taxa_and_divergence(
            df=df,
            col_specs=pen_specs,
            ylabel="4-Point Condition Violation L(T; D)",
            suptitle=f"Tree 4-Point Condition Violation (± SD) — {m_title}",
            out_path=os.path.join(outdir, f"penalty_obs_{m_type}_vs_taxa.png"),
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge HPC shard CSVs and generate plots.")
    parser.add_argument("--outdir", type=str, default="results_hpc", help="HPC results directory")
    args = parser.parse_args()

    df_merged = merge_shards(outdir=args.outdir)
    generate_all_hpc_plots(df_merged, outdir=args.outdir)
