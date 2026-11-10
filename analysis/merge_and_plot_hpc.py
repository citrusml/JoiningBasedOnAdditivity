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
        for candidate in [
            os.path.join(outdir, "hpc_benchmark_summary.csv"),
            os.path.join(outdir, "summary", "hpc_benchmark_summary.csv"),
        ]:
            if os.path.exists(candidate):
                return pd.read_csv(candidate)
        raise FileNotFoundError(f"No shard CSV files found matching {pattern}")

    dfs = [pd.read_csv(f) for f in files]
    df = pd.concat(dfs, ignore_index=True)
    dedup_cols = [c for c in ["n_taxa", "divergence", "indel_rate", "replicate"] if c in df.columns]
    df = df.drop_duplicates(subset=dedup_cols, keep="last")
    df = df.sort_values(dedup_cols).reset_index(drop=True)

    merged_path = os.path.join(outdir, "hpc_benchmark_summary.csv")
    df.to_csv(merged_path, index=False)
    print(f"[OK] Merged {len(files)} shard files ({len(df)} total rows) -> {merged_path}")

    # Also merge JSONL tree shards into a single consolidated all_trees.jsonl
    tree_pattern = os.path.join(outdir, "trees", "trees_N*.jsonl")
    tree_files = sorted(glob.glob(tree_pattern))
    if tree_files:
        all_trees_path = os.path.join(outdir, "trees", "all_trees.jsonl")
        with open(all_trees_path, "w") as out_f:
            for tf in tree_files:
                with open(tf, "r") as in_f:
                    for line in in_f:
                        if line.strip():
                            out_f.write(line)
        print(f"[OK] Merged {len(tree_files)} tree JSONL shards -> {all_trees_path}")

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


def _compute_lowess_curve(
    x: np.ndarray,
    y: np.ndarray,
    frac: float = 0.35,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fits a LOWESS smoother over x (SSS, unitless) and y (nRF or 1-nRF, unitless),
    falling back to non-robust local regression (it=0) if MAD collapses when many y=0.

    Parameters
    ----------
    x : np.ndarray
        1D array of Sequence Similarity Scores (unitless, in [0, 1]).
    y : np.ndarray
        1D array of topological error nRF or accuracy 1 - nRF (unitless, in [0, 1]).
    frac : float, default=0.35
        Fraction of data used when estimating each y-value in LOWESS.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Sorted smoothed (x_smooth, y_smooth) arrays (unitless).
    """
    valid = np.isfinite(x) & np.isfinite(y)
    x_v, y_v = x[valid], y[valid]
    if len(x_v) < 5:
        return np.array([]), np.array([])

    order = np.argsort(x_v)
    x_s, y_s = x_v[order], y_v[order]

    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess

        sm = lowess(y_s, x_s, frac=frac, it=3, return_sorted=True)
        # If robust iterations collapse due to many zeros (MAD=0), fallback to it=0
        if np.allclose(sm[:, 1], y_s) or np.max(np.abs(np.diff(sm[:, 1]))) > 0.08:
            sm = lowess(y_s, x_s, frac=frac, it=0, return_sorted=True)
        return sm[:, 0], np.clip(sm[:, 1], 0.0, 1.0)
    except Exception:
        s_ser = pd.Series(y_s).rolling(window=max(5, len(y_s) // 10), min_periods=1, center=True).mean()
        return x_s, np.clip(s_ser.to_numpy(), 0.0, 1.0)


def _plot_nrf_vs_sss_curves(
    df: pd.DataFrame,
    col_specs: List[Tuple[str, str, str, str]],
    subtitle: str,
    out_path: str,
    frac: float = 0.35,
) -> None:
    """
    Generates a 2-panel scatter + LOWESS smoothed curve figure matching
    `nrf_vs_sss_curves.png` (Matsui & Iwasaki 2020 Fig. 3a replication):
    - Left panel: Topological Error (nRF) vs Sequence Similarity Score (SSS)
    - Right panel: Topological Accuracy (1 - nRF) vs Sequence Similarity Score (SSS)

    Parameters
    ----------
    df : pd.DataFrame
        Benchmark summary DataFrame containing `sss_mean` (unitless) and `*_nrf` columns (unitless).
    col_specs : List[Tuple[str, str, str, str]]
        List of (column_name, label, color, linestyle) specifications.
    subtitle : str
        Subtitle describing the distance matrix type or comparison set.
    out_path : str
        Output PNG file path.
    frac : float, default=0.35
        LOWESS smoothing span parameter (unitless).
    """
    if "sss_mean" not in df.columns or df["sss_mean"].dropna().empty:
        return

    valid_specs = [spec for spec in col_specs if spec[0] in df.columns]
    if not valid_specs:
        return

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 7), dpi=300, sharex=True)

    for col_name, label, color, lstyle in valid_specs:
        sub = df[df["sss_mean"].notna() & df[col_name].notna()].copy()
        if len(sub) < 5:
            continue
        sub = sub.sort_values("sss_mean")
        x = sub["sss_mean"].to_numpy(dtype=float)
        y_err = sub[col_name].to_numpy(dtype=float)
        y_acc = 1.0 - y_err

        # Semitransparent scatter points
        ax1.scatter(x, y_err, color=color, alpha=0.12, s=14, edgecolors="none")
        ax2.scatter(x, y_acc, color=color, alpha=0.12, s=14, edgecolors="none")

        # LOWESS smooth curve
        lx_err, ly_err = _compute_lowess_curve(x, y_err, frac=frac)
        lx_acc, ly_acc = _compute_lowess_curve(x, y_acc, frac=frac)

        if len(lx_err) > 0:
            ax1.plot(lx_err, ly_err, color=color, linestyle=lstyle, linewidth=2.5, label=label)
        if len(lx_acc) > 0:
            ax2.plot(lx_acc, ly_acc, color=color, linestyle=lstyle, linewidth=2.5, label=label)

    max_sss = float(df["sss_mean"].max())
    xlim_upper = min(0.65, max(0.20, max_sss + 0.02))

    panel_configs = [
        (
            ax1,
            f"Topological Error (nRF) vs Sequence Similarity Score (SSS)\n[{subtitle}]",
            "Normalized RF Distance (Lower is Better)",
            "upper right",
        ),
        (
            ax2,
            f"Topological Accuracy (1 - nRF) vs Sequence Similarity Score (SSS)\n[Matsui & Iwasaki (2020) Fig. 3a Replication — {subtitle}]",
            "Correct Topology Ratio (1 - nRF, Higher is Better)",
            "lower right",
        ),
    ]

    for ax, title, ylabel, leg_loc in panel_configs:
        ax.set_title(title, fontsize=11.5, fontweight="bold", pad=12)
        ax.set_xlabel(r"Average Sequence Similarity Score ($\bar{w}$)", fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_xlim(-0.01, xlim_upper)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, linestyle="--", alpha=0.5)
        ax.legend(fontsize=9, loc=leg_loc, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved figure: {out_path}")


def _plot_nrf_vs_sss_by_taxa(
    df: pd.DataFrame,
    col_specs: List[Tuple[str, str, str, str]],
    subtitle: str,
    out_path: str,
    frac: float = 0.42,
) -> None:
    """
    Generates a 4-panel figure stratified by Taxon Count N (N = 20, 32, 64, 128),
    plotting Topological Error (nRF) vs SSS with scatter points and LOWESS curves.

    Parameters
    ----------
    df : pd.DataFrame
        Benchmark summary DataFrame containing `n_taxa`, `sss_mean`, and `*_nrf` columns.
    col_specs : List[Tuple[str, str, str, str]]
        List of (column_name, label, color, linestyle) specifications.
    subtitle : str
        Subtitle describing the distance matrix type.
    out_path : str
        Output PNG file path.
    frac : float, default=0.42
        LOWESS smoothing span parameter (unitless).
    """
    if "sss_mean" not in df.columns or "n_taxa" not in df.columns:
        return

    taxa_list = sorted(df["n_taxa"].unique())
    valid_specs = [spec for spec in col_specs if spec[0] in df.columns]
    if not valid_specs or not taxa_list:
        return

    fig, axes = plt.subplots(1, len(taxa_list), figsize=(4.6 * len(taxa_list), 5.2), dpi=300, sharey=True)
    if len(taxa_list) == 1:
        axes = [axes]

    for ax, n_val in zip(axes, taxa_list):
        sub_n = df[df["n_taxa"] == n_val]
        for col_name, label, color, lstyle in valid_specs:
            sub = sub_n[sub_n["sss_mean"].notna() & sub_n[col_name].notna()].copy()
            if len(sub) < 5:
                continue
            sub = sub.sort_values("sss_mean")
            x = sub["sss_mean"].to_numpy(dtype=float)
            y_err = sub[col_name].to_numpy(dtype=float)

            ax.scatter(x, y_err, color=color, alpha=0.16, s=14, edgecolors="none")
            lx_err, ly_err = _compute_lowess_curve(x, y_err, frac=frac)
            if len(lx_err) > 0:
                ax.plot(lx_err, ly_err, color=color, linestyle=lstyle, linewidth=2.2, label=label)

        ax.set_title(f"Taxa N = {n_val}", fontsize=12, fontweight="bold", pad=10)
        ax.set_xlabel(r"Average SSS ($\bar{w}$)", fontsize=10.5)
        ax.set_ylim(-0.02, 0.52)
        ax.grid(True, linestyle="--", alpha=0.5)

    axes[0].set_ylabel("Normalized RF Distance (nRF, Lower is Better)", fontsize=11)
    axes[-1].legend(bbox_to_anchor=(1.02, 1.0), loc="upper left", fontsize=8.5, framealpha=0.9)
    fig.suptitle(f"Topological Error (nRF) vs SSS Stratified by Taxa N — {subtitle}", fontsize=13.5, fontweight="bold", y=1.02)
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
        "sss_linear": "MMseqs2 SSS Linear Distance (d = 1 - w, Exp 18)",
        "sss_log": "MMseqs2 SSS Log Distance (d = -ln(w), Exp 18)",
        "true": "True Patristic Distance",
    }

    # 1. nRF plots for each of the 3 distance matrix types
    for m_type, m_title in matrix_labels.items():
        nrf_specs = [
            ("MSA_ML_nrf", "MSA + ML (FastTree)", "#d62728", "o"),
            (f"{m_type}_Standard_NJ_nrf", "Standard NJ", "#1f77b4", "s"),
            (f"{m_type}_BioNJ_nrf", "BioNJ", "#ff7f0e", "^"),
            (f"{m_type}_FastME_NNI_nrf", "FastME-NNI", "#2ca02c", "v"),
            (f"{m_type}_FastME_SPR_nrf", "FastME-SPR", "#17becf", "D"),
            (f"{m_type}_AddNJ2_nrf", "AddNJ2", "#9467bd", "P"),
            (f"{m_type}_BioNJ_Add_NNI_nrf", "BioNJ + Add-NNI", "#8c564b", "<"),
            (f"{m_type}_AddNJ2_Add_NNI_nrf", "AddNJ2 + Add-NNI", "#e377c2", "X"),
        ]
        _plot_metric_by_taxa_and_divergence(
            df=df,
            col_specs=nrf_specs,
            ylabel="Normalized Robinson-Foulds Distance (nRF)",
            suptitle=f"Topological Error (nRF ± SD) — {m_title}",
            out_path=os.path.join(outdir, f"nrf_{m_type}_vs_taxa.png"),
        )
        if m_type == "blosum62":
            # Also save default nrf_vs_taxa.png matching results_large_scale/nrf_vs_taxa.png
            _plot_metric_by_taxa_and_divergence(
                df=df,
                col_specs=nrf_specs,
                ylabel="Normalized Robinson-Foulds Distance (nRF)",
                suptitle=f"Topological Error (nRF ± SD) — {m_title}",
                out_path=os.path.join(outdir, "nrf_vs_taxa.png"),
            )

        # 1b. KF / Branch Score Distance (BSD) plots if computed
        if "MSA_ML_kf_bsd" in df.columns:
            kf_specs = [
                ("MSA_ML_kf_bsd", "MSA + ML (FastTree)", "#d62728", "o"),
                (f"{m_type}_Standard_NJ_kf_bsd", "Standard NJ", "#1f77b4", "s"),
                (f"{m_type}_BioNJ_kf_bsd", "BioNJ", "#ff7f0e", "^"),
                (f"{m_type}_FastME_NNI_kf_bsd", "FastME-NNI", "#2ca02c", "v"),
                (f"{m_type}_FastME_SPR_kf_bsd", "FastME-SPR", "#17becf", "D"),
                (f"{m_type}_AddNJ2_kf_bsd", "AddNJ2", "#9467bd", "P"),
                (f"{m_type}_BioNJ_Add_NNI_kf_bsd", "BioNJ + Add-NNI", "#8c564b", "<"),
                (f"{m_type}_AddNJ2_Add_NNI_kf_bsd", "AddNJ2 + Add-NNI", "#e377c2", "X"),
            ]
            _plot_metric_by_taxa_and_divergence(
                df=df,
                col_specs=kf_specs,
                ylabel="Kuhner-Felsenstein Distance (KF / BSD)",
                suptitle=f"Branch Score Distance (KF/BSD ± SD) — {m_title}",
                out_path=os.path.join(outdir, f"kf_bsd_{m_type}_vs_taxa.png"),
            )

        # 1c. Patristic Distance RMSE plots if computed
        if "MSA_ML_patristic_rmse" in df.columns:
            prmse_specs = [
                ("MSA_ML_patristic_rmse", "MSA + ML (FastTree)", "#d62728", "o"),
                (f"{m_type}_Standard_NJ_patristic_rmse", "Standard NJ", "#1f77b4", "s"),
                (f"{m_type}_BioNJ_patristic_rmse", "BioNJ", "#ff7f0e", "^"),
                (f"{m_type}_FastME_NNI_patristic_rmse", "FastME-NNI", "#2ca02c", "v"),
                (f"{m_type}_FastME_SPR_patristic_rmse", "FastME-SPR", "#17becf", "D"),
                (f"{m_type}_AddNJ2_patristic_rmse", "AddNJ2", "#9467bd", "P"),
                (f"{m_type}_BioNJ_Add_NNI_patristic_rmse", "BioNJ + Add-NNI", "#8c564b", "<"),
                (f"{m_type}_AddNJ2_Add_NNI_patristic_rmse", "AddNJ2 + Add-NNI", "#e377c2", "X"),
            ]
            _plot_metric_by_taxa_and_divergence(
                df=df,
                col_specs=prmse_specs,
                ylabel="Patristic Distance RMSE (substitutions/site)",
                suptitle=f"Tree Patristic Distance RMSE (± SD) — {m_title}",
                out_path=os.path.join(outdir, f"patristic_rmse_{m_type}_vs_taxa.png"),
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
        ("blosum62_AddNJ2_Add_NNI_nrf", "BLOSUM62: AddNJ2+Add-NNI", "#e377c2", "X"),
        ("wag_gamma_AddNJ2_Add_NNI_nrf", "WAG+Γ: AddNJ2+Add-NNI", "#f7b6d2", "X"),
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
            (f"{m_type}_FastME_NNI_penalty_obs", "FastME-NNI", "#2ca02c", "v"),
            (f"{m_type}_FastME_SPR_penalty_obs", "FastME-SPR", "#17becf", "D"),
            (f"{m_type}_AddNJ2_penalty_obs", "AddNJ2", "#9467bd", "P"),
            (f"{m_type}_BioNJ_Add_NNI_penalty_obs", "BioNJ + Add-NNI", "#8c564b", "<"),
            (f"{m_type}_AddNJ2_Add_NNI_penalty_obs", "AddNJ2 + Add-NNI", "#e377c2", "X"),
        ]
        _plot_metric_by_taxa_and_divergence(
            df=df,
            col_specs=pen_specs,
            ylabel="4-Point Condition Violation L(T; D)",
            suptitle=f"Tree 4-Point Condition Violation (± SD) — {m_title}",
            out_path=os.path.join(outdir, f"penalty_obs_{m_type}_vs_taxa.png"),
        )

    # 4. nRF & Accuracy (1 - nRF) vs SSS Scatter + LOWESS Curves (Matsui & Iwasaki 2020 Fig. 3a style)
    for m_type, m_title in matrix_labels.items():
        sss_curve_specs = [
            ("MSA_ML_nrf", "MSA + ML (FastTree)", "#d62728", "-"),
            (f"{m_type}_Standard_NJ_nrf", "Standard NJ", "#1f77b4", "-"),
            (f"{m_type}_BioNJ_nrf", "BioNJ", "#ff7f0e", "-"),
            (f"{m_type}_FastME_NNI_nrf", "FastME-NNI", "#2ca02c", "--"),
            (f"{m_type}_FastME_SPR_nrf", "FastME-SPR", "#17becf", "-"),
            (f"{m_type}_AddNJ2_nrf", "AddNJ2", "#9467bd", "-"),
            (f"{m_type}_BioNJ_Add_NNI_nrf", "BioNJ + Add-NNI", "#8c564b", "--"),
            (f"{m_type}_AddNJ2_Add_NNI_nrf", "AddNJ2 + Add-NNI", "#e377c2", "-."),
        ]
        _plot_nrf_vs_sss_curves(
            df=df,
            col_specs=sss_curve_specs,
            subtitle=m_title,
            out_path=os.path.join(outdir, f"nrf_vs_sss_curves_{m_type}.png"),
            frac=0.35,
        )
        _plot_nrf_vs_sss_by_taxa(
            df=df,
            col_specs=sss_curve_specs,
            subtitle=m_title,
            out_path=os.path.join(outdir, f"nrf_vs_sss_by_taxa_{m_type}.png"),
            frac=0.45,
        )

    # Main combined nrf_vs_sss_curves.png (MSA+ML vs BLOSUM62 & WAG+Gamma all methods)
    main_sss_specs = [
        ("MSA_ML_nrf", "MSA + ML (FastTree)", "#d62728", "-"),
        ("blosum62_Standard_NJ_nrf", "BLOSUM62: Standard NJ", "#1f77b4", "-"),
        ("blosum62_BioNJ_nrf", "BLOSUM62: BioNJ", "#ff7f0e", "-"),
        ("blosum62_FastME_NNI_nrf", "BLOSUM62: FastME-NNI", "#2ca02c", ":"),
        ("blosum62_FastME_SPR_nrf", "BLOSUM62: FastME-SPR", "#17becf", "-"),
        ("blosum62_AddNJ2_nrf", "BLOSUM62: AddNJ2", "#9467bd", "-"),
        ("blosum62_BioNJ_Add_NNI_nrf", "BLOSUM62: BioNJ + Add-NNI", "#8c564b", "-"),
        ("blosum62_AddNJ2_Add_NNI_nrf", "BLOSUM62: AddNJ2 + Add-NNI", "#e377c2", "-"),
        ("wag_gamma_Standard_NJ_nrf", "WAG+Γ: Standard NJ", "#1f77b4", "--"),
        ("wag_gamma_BioNJ_nrf", "WAG+Γ: BioNJ", "#ff7f0e", "--"),
        ("wag_gamma_FastME_NNI_nrf", "WAG+Γ: FastME-NNI", "#2ca02c", "--"),
        ("wag_gamma_FastME_SPR_nrf", "WAG+Γ: FastME-SPR", "#17becf", "--"),
        ("wag_gamma_AddNJ2_nrf", "WAG+Γ: AddNJ2", "#9467bd", "--"),
        ("wag_gamma_BioNJ_Add_NNI_nrf", "WAG+Γ: BioNJ + Add-NNI", "#8c564b", "--"),
        ("wag_gamma_AddNJ2_Add_NNI_nrf", "WAG+Γ: AddNJ2 + Add-NNI", "#e377c2", "--"),
    ]
    _plot_nrf_vs_sss_curves(
        df=df,
        col_specs=main_sss_specs,
        subtitle="BLOSUM62 & WAG+Γ vs MSA+ML",
        out_path=os.path.join(outdir, "nrf_vs_sss_curves.png"),
        frac=0.35,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge HPC shard CSVs and generate plots.")
    parser.add_argument("--outdir", type=str, default="results_hpc", help="HPC results directory")
    args = parser.parse_args()

    df_merged = merge_shards(outdir=args.outdir)
    generate_all_hpc_plots(df_merged, outdir=args.outdir)

