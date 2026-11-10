"""
Publication-Quality Visualization Script for LBA (Long Branch Attraction) Benchmark.

Plots:
1. Shortest path edge distance DS(T1, TN) across methods (Control vs LBA).
2. False Sister Rate (% of trials where T1 and TN are erroneously inferred as sister taxa, DS=2).
3. Overall topological error (nRF) under Control vs LBA.
4. Computational runtime comparison across methods.
"""

from typing import List, Dict, Tuple, Optional
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


METHOD_SPECS = [
    ("MSA_ML", "MSA + ML (FastTree)", "#d62728", "o"),
    ("Standard_NJ", "Standard NJ", "#1f77b4", "s"),
    ("BioNJ", "BioNJ", "#ff7f0e", "^"),
    ("FastME_NNI", "FastME-NNI", "#17becf", "v"),
    ("FastME_SPR", "FastME-SPR", "#2ca02c", "D"),
    ("AddNJ", "AddNJ (Original)", "#7f7f7f", "x"),
    ("AddNJ2", "AddNJ2 (Proposed)", "#9467bd", "P"),
    ("NJ_Add_NNI", "NJ + Add-NNI", "#bcbd22", "<"),
    ("NJ_Add_SPR", "NJ + Add-SPR", "#8c564b", ">"),
    ("NJ_Add_Hybrid", "NJ + Add-Hybrid", "#e377c2", "h"),
    ("BioNJ_Add_NNI", "BioNJ + Add-NNI", "#aec7e8", "p"),
    ("BioNJ_Add_SPR", "BioNJ + Add-SPR", "#ffbb78", "*"),
    ("BioNJ_Add_Hybrid", "BioNJ + Add-Hybrid", "#2ca02c", "H"),
    ("AddNJ2_Add_NNI", "AddNJ2 + Add-NNI", "#c5b0d5", "+"),
    ("AddNJ2_Add_Hybrid", "AddNJ2 + Add-Hybrid", "#800080", "X"),
]


def plot_lba_ds_comparison(df: pd.DataFrame, outdir: str = "results_lba") -> None:
    """
    Plots the topological shortest path edge distance DS(T1, TN) between long branches.
    Shows the dramatic drop toward DS=2 when LBA occurs, and how well each method resists.
    """
    os.makedirs(outdir, exist_ok=True)
    taxa_list = sorted(df["n_taxa"].unique())
    divergences = sorted(df["divergence"].unique())
    lba_ratios = sorted(df["lba_ratio"].unique())

    methods = [spec for spec in METHOD_SPECS if f"{spec[0]}_ds_est" in df.columns]

    for n_taxa in taxa_list:
        sub_df = df[df["n_taxa"] == n_taxa]
        fig, axes = plt.subplots(
            len(lba_ratios),
            len(divergences),
            figsize=(6 * len(divergences), 4.5 * len(lba_ratios)),
            dpi=300,
            sharey="row",
        )

        if len(lba_ratios) == 1 and len(divergences) == 1:
            axes = np.array([[axes]])
        elif len(lba_ratios) == 1:
            axes = np.expand_dims(axes, 0)
        elif len(divergences) == 1:
            axes = np.expand_dims(axes, 1)

        for r_idx, lba_r in enumerate(lba_ratios):
            for d_idx, div in enumerate(divergences):
                ax = axes[r_idx, d_idx]
                cell_df = sub_df[(sub_df["divergence"] == div) & (sub_df["lba_ratio"] == lba_r)]

                if cell_df.empty:
                    continue

                ds_true_mean = cell_df["ds_true"].mean()
                ax.axhline(
                    ds_true_mean,
                    color="gray",
                    linestyle="--",
                    linewidth=1.2,
                    label=f"True DS ({ds_true_mean:.1f} edges)" if d_idx == 0 and r_idx == 0 else "",
                )
                ax.axhline(
                    2.0,
                    color="red",
                    linestyle=":",
                    linewidth=1.5,
                    label="Sister Artifact (DS=2)" if d_idx == 0 and r_idx == 0 else "",
                )

                x_positions = np.arange(len(methods))
                means = []
                sems = []
                for m_id, _, _, _ in methods:
                    vals = cell_df[f"{m_id}_ds_est"].values
                    means.append(np.mean(vals))
                    sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0.0)

                colors = [m[2] for m in methods]
                bars = ax.bar(x_positions, means, yerr=sems, capsize=3, color=colors, alpha=0.85, edgecolor="black", linewidth=0.5)

                title_lba = "Control (No LBA, ratio=1.0)" if lba_r == 1.0 else f"LBA Condition (ratio={lba_r})"
                ax.set_title(f"N={n_taxa}, D={div} | {title_lba}", fontsize=11, fontweight="bold")
                ax.set_xticks(x_positions)
                ax.set_xticklabels([m[1] for m in methods], rotation=45, ha="right", fontsize=8)
                if d_idx == 0:
                    ax.set_ylabel("Shortest Path Edges DS(T1, TN)", fontsize=10)
                ax.grid(axis="y", linestyle=":", alpha=0.6)
                ax.set_ylim(0, max(ds_true_mean + 3, 10))

        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.98, 0.98), fontsize=9)

        plt.tight_layout()
        out_path = os.path.join(outdir, f"lba_ds_comparison_N{n_taxa}.png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"[OK] Saved DS comparison plot: {out_path}")


def plot_lba_sister_rate(df: pd.DataFrame, outdir: str = "results_lba") -> None:
    """
    Plots the False Sister Rate (%) indicating direct Long Branch Attraction artifact occurrence.
    """
    os.makedirs(outdir, exist_ok=True)
    lba_df = df[df["lba_ratio"] > 1.0]
    if lba_df.empty:
        print("[WARN] No LBA ratio > 1.0 found for sister rate plotting.")
        return

    taxa_list = sorted(lba_df["n_taxa"].unique())
    divergences = sorted(lba_df["divergence"].unique())
    methods = [spec for spec in METHOD_SPECS if f"{spec[0]}_is_sister" in df.columns]

    fig, axes = plt.subplots(
        len(taxa_list),
        len(divergences),
        figsize=(6 * len(divergences), 4.5 * len(taxa_list)),
        dpi=300,
        sharey=True,
    )

    if len(taxa_list) == 1 and len(divergences) == 1:
        axes = np.array([[axes]])
    elif len(taxa_list) == 1:
        axes = np.expand_dims(axes, 0)
    elif len(divergences) == 1:
        axes = np.expand_dims(axes, 1)

    for t_idx, n_taxa in enumerate(taxa_list):
        for d_idx, div in enumerate(divergences):
            ax = axes[t_idx, d_idx]
            cell_df = lba_df[(lba_df["n_taxa"] == n_taxa) & (lba_df["divergence"] == div)]

            if cell_df.empty:
                continue

            rates = []
            for m_id, _, _, _ in methods:
                sister_vals = cell_df[f"{m_id}_is_sister"].values
                rate = float(np.mean(sister_vals) * 100.0) if len(sister_vals) > 0 else 0.0
                rates.append(rate)

            x_pos = np.arange(len(methods))
            colors = [m[2] for m in methods]
            bars = ax.bar(x_pos, rates, color=colors, alpha=0.85, edgecolor="black", linewidth=0.5)

            for bar in bars:
                h = bar.get_height()
                if h > 0:
                    ax.text(bar.get_x() + bar.get_width() / 2.0, h + 1.5, f"{h:.0f}%", ha="center", va="bottom", fontsize=7)

            ax.set_title(f"False Sister Rate (%): N={n_taxa}, D={div} (LBA Ratio=3.0)", fontsize=10, fontweight="bold")
            ax.set_xticks(x_pos)
            ax.set_xticklabels([m[1] for m in methods], rotation=45, ha="right", fontsize=8)
            ax.set_ylim(0, 105)
            if d_idx == 0:
                ax.set_ylabel("False Sister Rate (%)", fontsize=10)
            ax.grid(axis="y", linestyle=":", alpha=0.6)

    plt.tight_layout()
    out_path = os.path.join(outdir, "lba_false_sister_rate.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved False Sister Rate plot: {out_path}")


def plot_lba_nrf_comparison(df: pd.DataFrame, outdir: str = "results_lba") -> None:
    """
    Plots overall topological error (nRF) under Control vs LBA across all methods.
    """
    os.makedirs(outdir, exist_ok=True)
    taxa_list = sorted(df["n_taxa"].unique())
    divergences = sorted(df["divergence"].unique())
    methods = [spec for spec in METHOD_SPECS if f"{spec[0]}_nrf" in df.columns]

    for n_taxa in taxa_list:
        sub_df = df[df["n_taxa"] == n_taxa]
        fig, axes = plt.subplots(1, len(divergences), figsize=(6 * len(divergences), 5.5), dpi=300, sharey=True)

        if len(divergences) == 1:
            axes = [axes]

        for d_idx, div in enumerate(divergences):
            ax = axes[d_idx]
            cell_ctrl = sub_df[(sub_df["divergence"] == div) & (sub_df["lba_ratio"] == 1.0)]
            cell_lba = sub_df[(sub_df["divergence"] == div) & (sub_df["lba_ratio"] > 1.0)]

            x_pos = np.arange(len(methods))
            width = 0.38

            ctrl_means = [cell_ctrl[f"{m[0]}_nrf"].mean() if not cell_ctrl.empty else 0.0 for m in methods]
            lba_means = [cell_lba[f"{m[0]}_nrf"].mean() if not cell_lba.empty else 0.0 for m in methods]

            ax.bar(x_pos - width / 2, ctrl_means, width, label="Control (ratio=1.0)", color="#4682b4", alpha=0.85, edgecolor="black", linewidth=0.5)
            ax.bar(x_pos + width / 2, lba_means, width, label="LBA (ratio=3.0)", color="#d95f02", alpha=0.85, edgecolor="black", linewidth=0.5)

            ax.set_title(f"nRF Error: N={n_taxa}, D={div}", fontsize=11, fontweight="bold")
            ax.set_xticks(x_pos)
            ax.set_xticklabels([m[1] for m in methods], rotation=45, ha="right", fontsize=8)
            ax.set_ylim(0, max(max(ctrl_means + lba_means) * 1.25, 0.2))
            if d_idx == 0:
                ax.set_ylabel("Normalized Robinson-Foulds (nRF)", fontsize=10)
                ax.legend(loc="upper left", fontsize=9)
            ax.grid(axis="y", linestyle=":", alpha=0.6)

        plt.tight_layout()
        out_path = os.path.join(outdir, f"lba_nrf_comparison_N{n_taxa}.png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"[OK] Saved nRF comparison plot: {out_path}")


def plot_lba_runtime(df: pd.DataFrame, outdir: str = "results_lba") -> None:
    """
    Plots the computational runtime (seconds) across all methods.
    """
    os.makedirs(outdir, exist_ok=True)
    methods = [spec for spec in METHOD_SPECS if f"{spec[0]}_time" in df.columns]

    mean_times = []
    labels = []
    colors = []
    for m in methods:
        val = df[f"{m[0]}_time"].mean()
        mean_times.append(val)
        labels.append(m[1])
        colors.append(m[2])

    fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
    y_pos = np.arange(len(methods))
    bars = ax.barh(y_pos, mean_times, color=colors, alpha=0.85, edgecolor="black", linewidth=0.5)

    for bar in bars:
        w = bar.get_width()
        ax.text(w * 1.05 + 1e-4, bar.get_y() + bar.get_height() / 2.0, f"{w:.4f}s", va="center", fontsize=8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlabel("Mean Execution Time per Tree [seconds] (log scale)", fontsize=10)
    ax.set_title("Computational Runtime Comparison across All Methods", fontsize=11, fontweight="bold")
    ax.grid(axis="x", linestyle=":", alpha=0.6)

    plt.tight_layout()
    out_path = os.path.join(outdir, "lba_runtime_comparison.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved runtime comparison plot: {out_path}")


def generate_all_lba_plots(df: pd.DataFrame, outdir: str = "results_lba") -> None:
    """Generates all 4 publication charts."""
    plot_lba_ds_comparison(df, outdir)
    plot_lba_sister_rate(df, outdir)
    plot_lba_nrf_comparison(df, outdir)
    plot_lba_runtime(df, outdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate publication charts for LBA benchmark.")
    parser.add_argument("--csv", type=str, default="results_lba/lba_summary.csv", help="Path to summary CSV")
    parser.add_argument("--outdir", type=str, default="results_lba", help="Output directory")
    args = parser.parse_args()

    if os.path.exists(args.csv):
        df_summary = pd.read_csv(args.csv)
        generate_all_lba_plots(df_summary, outdir=args.outdir)
    else:
        print(f"[ERROR] CSV not found at {args.csv}")
