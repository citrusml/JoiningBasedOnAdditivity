"""
Publication-Quality Visualization Script for LBA (Long Branch Attraction) Benchmark
across Multiple Distance Matrices (BLOSUM62+Poisson, WAG+Gamma ML, True Patristic)
with SSS and 4-Point Condition Violation Analysis.

Plots:
1. nRF comparison across the 3 distance matrices (BLOSUM62, WAG+Gamma, True) for Control vs LBA.
2. Shortest path edge distance DS(T1, TN) across all 15 methods (Control vs LBA).
3. Sequence Similarity Score (SSS) and 4-Point Condition Violation (Delta_ADD & Tree Penalty).
4. False Sister Rate (%) and Computational Runtime Comparison.
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
            figsize=(6.2 * len(divergences), 4.6 * len(lba_ratios)),
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
                ax.bar(x_positions, means, yerr=sems, capsize=3, color=colors, alpha=0.85, edgecolor="black", linewidth=0.5)

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


def plot_lba_nrf_by_matrix(df: pd.DataFrame, outdir: str = "results_lba") -> None:
    """
    Plots nRF comparison across the 3 distance matrices (BLOSUM62+Poisson, WAG+Gamma, True)
    and MSA+ML under LBA (ratio=3.0) vs Control (ratio=1.0).
    """
    os.makedirs(outdir, exist_ok=True)
    taxa_list = sorted(df["n_taxa"].unique())
    divergences = sorted(df["divergence"].unique())
    dist_methods = [spec for spec in METHOD_SPECS if spec[0] != "MSA_ML" and f"blosum62_{spec[0]}_nrf" in df.columns]

    for n_taxa in taxa_list:
        sub_df = df[df["n_taxa"] == n_taxa]
        fig, axes = plt.subplots(2, len(divergences), figsize=(6.5 * len(divergences), 9.5), dpi=300, sharey="row")

        for r_idx, lba_r in enumerate([1.0, 3.0]):
            for d_idx, div in enumerate(divergences):
                ax = axes[r_idx, d_idx]
                cell_df = sub_df[(sub_df["divergence"] == div) & (sub_df["lba_ratio"] == lba_r)]
                if cell_df.empty:
                    continue

                x_pos = np.arange(len(dist_methods))
                width = 0.27

                b_means = [cell_df[f"blosum62_{m[0]}_nrf"].mean() for m in dist_methods]
                w_means = [cell_df[f"wag_gamma_{m[0]}_nrf"].mean() for m in dist_methods]
                t_means = [cell_df[f"true_{m[0]}_nrf"].mean() for m in dist_methods]

                ax.bar(x_pos - width, b_means, width, label="BLOSUM62 + Poisson", color="#1f77b4", alpha=0.88, edgecolor="black", linewidth=0.4)
                ax.bar(x_pos, w_means, width, label="WAG + Gamma (ML Dist)", color="#ff7f0e", alpha=0.88, edgecolor="black", linewidth=0.4)
                ax.bar(x_pos + width, t_means, width, label="True Patristic Dist", color="#2ca02c", alpha=0.88, edgecolor="black", linewidth=0.4)

                if "MSA_ML_nrf" in cell_df.columns:
                    ml_mean = cell_df["MSA_ML_nrf"].mean()
                    ax.axhline(ml_mean, color="#d62728", linestyle="--", linewidth=1.3, label=f"MSA+ML (FastTree: {ml_mean:.3f})")

                lba_label = "Control (ratio=1.0)" if lba_r == 1.0 else "LBA (ratio=3.0)"
                ax.set_title(f"N={n_taxa}, D={div} | {lba_label}", fontsize=10.5, fontweight="bold")
                ax.set_xticks(x_pos)
                ax.set_xticklabels([m[1] for m in dist_methods], rotation=45, ha="right", fontsize=7.5)
                if d_idx == 0:
                    ax.set_ylabel("Normalized Robinson-Foulds (nRF)", fontsize=10)
                    ax.legend(loc="upper right", fontsize=7.5)
                ax.grid(axis="y", linestyle=":", alpha=0.6)

        plt.tight_layout()
        out_path = os.path.join(outdir, f"lba_nrf_by_matrix_N{n_taxa}.png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"[OK] Saved matrix-comparison nRF plot: {out_path}")


def plot_sss_and_four_point_defect(df: pd.DataFrame, outdir: str = "results_lba") -> None:
    """
    Plots Sequence Similarity Score (SSS) and 4-Point Condition Violation (Delta_ADD)
    across conditions (N, D, LBA ratio).
    """
    os.makedirs(outdir, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), dpi=300)

    cond_labels = []
    sss_all_means, sss_lba_means = [], []
    fp_blosum_means, fp_wag_means, fp_true_means = [], [], []
    pen_nj_b, pen_addnj2_b, pen_hyb_b = [], [], []

    for (n_taxa, div, lba_r), grp in df.groupby(["n_taxa", "divergence", "lba_ratio"]):
        lbl = f"N{n_taxa}\nD{div}\nL{lba_r:.0f}"
        cond_labels.append(lbl)
        sss_all_means.append(grp["sss_mean"].mean())
        sss_lba_means.append(grp["sss_lba_pair"].mean() if "sss_lba_pair" in grp else 0.0)
        fp_blosum_means.append(grp["delta_add_blosum62"].mean())
        fp_wag_means.append(grp["delta_add_wag_gamma"].mean())
        fp_true_means.append(grp["delta_add_true"].mean())
        pen_nj_b.append(grp["blosum62_Standard_NJ_penalty_obs"].mean())
        pen_addnj2_b.append(grp["blosum62_AddNJ2_penalty_obs"].mean())
        pen_hyb_b.append(grp["blosum62_AddNJ2_Add_Hybrid_penalty_obs"].mean())

    x = np.arange(len(cond_labels))
    w = 0.35

    # Panel 1: SSS
    axes[0].bar(x - w / 2, sss_all_means, w, label="Mean SSS (All Pairs)", color="#1f77b4", edgecolor="black", linewidth=0.5)
    axes[0].bar(x + w / 2, sss_lba_means, w, label="LBA Pair SSS w(T1, TN)", color="#d62728", edgecolor="black", linewidth=0.5)
    axes[0].set_title("Sequence Similarity Score (SSS)", fontsize=11, fontweight="bold")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(cond_labels, fontsize=7.5)
    axes[0].set_ylabel("SSS (w_ij in [0, 1])", fontsize=10)
    axes[0].legend(fontsize=8.5)
    axes[0].grid(axis="y", linestyle=":", alpha=0.6)

    # Panel 2: Matrix 4-Point Condition Defect Delta_ADD(D)
    axes[1].bar(x - w, fp_blosum_means, w, label="BLOSUM62+Poisson", color="#1f77b4", edgecolor="black", linewidth=0.5)
    axes[1].bar(x, fp_wag_means, w, label="WAG+Gamma ML", color="#ff7f0e", edgecolor="black", linewidth=0.5)
    axes[1].bar(x + w, fp_true_means, w, label="True Patristic (0.0)", color="#2ca02c", edgecolor="black", linewidth=0.5)
    axes[1].set_title("Matrix 4-Point Condition Defect (Delta_ADD)", fontsize=11, fontweight="bold")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(cond_labels, fontsize=7.5)
    axes[1].set_ylabel("Relative Additivity Defect Delta_ADD(D)", fontsize=10)
    axes[1].legend(fontsize=8.5)
    axes[1].grid(axis="y", linestyle=":", alpha=0.6)

    # Panel 3: Tree Quartet Additivity Penalty L(T | D_obs)
    w3 = 0.27
    axes[2].bar(x - w3, pen_nj_b, w3, label="Standard NJ", color="#1f77b4", edgecolor="black", linewidth=0.5)
    axes[2].bar(x, pen_addnj2_b, w3, label="AddNJ2", color="#9467bd", edgecolor="black", linewidth=0.5)
    axes[2].bar(x + w3, pen_hyb_b, w3, label="AddNJ2 + Add-Hybrid", color="#800080", edgecolor="black", linewidth=0.5)
    axes[2].set_title("Tree Quartet Additivity Penalty L(T | D_obs)", fontsize=11, fontweight="bold")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(cond_labels, fontsize=7.5)
    axes[2].set_ylabel("Quartet Additivity Penalty L(T)", fontsize=10)
    axes[2].legend(fontsize=8.5)
    axes[2].grid(axis="y", linestyle=":", alpha=0.6)

    plt.tight_layout()
    out_path = os.path.join(outdir, "lba_4pc_violation_and_sss.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved SSS and 4PC violation plot: {out_path}")


def plot_lba_sister_rate(df: pd.DataFrame, outdir: str = "results_lba") -> None:
    """
    Plots the False Sister Rate (%) indicating direct Long Branch Attraction artifact occurrence.
    """
    os.makedirs(outdir, exist_ok=True)
    lba_df = df[df["lba_ratio"] > 1.0]
    if lba_df.empty:
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
        fig, axes = plt.subplots(1, len(divergences), figsize=(6.2 * len(divergences), 5.5), dpi=300, sharey=True)

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

    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=300)
    y_pos = np.arange(len(methods))
    bars = ax.barh(y_pos, mean_times, color=colors, alpha=0.85, edgecolor="black", linewidth=0.5)

    for bar in bars:
        w = bar.get_width()
        ax.text(w * 1.05 + 1e-4, bar.get_y() + bar.get_height() / 2.0, f"{w:.4f}s", va="center", fontsize=8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlabel("Mean Execution Time per Tree [seconds] (log scale)", fontsize=10)
    ax.set_title("Computational Runtime Comparison across All 15 Methods", fontsize=11, fontweight="bold")
    ax.grid(axis="x", linestyle=":", alpha=0.6)

    plt.tight_layout()
    out_path = os.path.join(outdir, "lba_runtime_comparison.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved runtime comparison plot: {out_path}")


def generate_all_lba_plots(df: pd.DataFrame, outdir: str = "results_lba") -> None:
    """Generates all publication charts."""
    plot_lba_ds_comparison(df, outdir)
    plot_lba_nrf_by_matrix(df, outdir)
    plot_sss_and_four_point_defect(df, outdir)
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
