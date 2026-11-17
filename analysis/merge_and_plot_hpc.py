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
    is_log_x: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fits a LOWESS smoother over x (SSS, unitless) and y (nRF or 1-nRF, unitless).
    When is_log_x is True, smoothing is performed in log10(x) space to preserve
    proper local density and curvature across multiple decades of sequence similarity.

    Parameters
    ----------
    x : np.ndarray
        1D array of Sequence Similarity Scores (unitless, in [0, 1]).
    y : np.ndarray
        1D array of topological error nRF or accuracy 1 - nRF (unitless, in [0, 1]).
    frac : float, default=0.35
        Fraction of data used when estimating each y-value in LOWESS.
    is_log_x : bool, default=True
        Whether to perform LOWESS smoothing in log10(x) space.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Sorted smoothed (x_smooth, y_smooth) arrays (unitless).
    """
    valid = np.isfinite(x) & np.isfinite(y)
    if is_log_x:
        valid = valid & (x > 0)
    x_v, y_v = x[valid], y[valid]
    if len(x_v) < 5:
        return np.array([]), np.array([])

    order = np.argsort(x_v)
    x_s, y_s = x_v[order], y_v[order]

    fit_x = np.log10(x_s) if is_log_x else x_s

    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess

        sm = lowess(y_s, fit_x, frac=frac, it=3, return_sorted=True)
        # If robust iterations collapse due to many zeros (MAD=0), fallback to it=0
        if np.allclose(sm[:, 1], y_s) or np.max(np.abs(np.diff(sm[:, 1]))) > 0.08:
            sm = lowess(y_s, fit_x, frac=frac, it=0, return_sorted=True)
        res_x = 10 ** sm[:, 0] if is_log_x else sm[:, 0]
        return res_x, np.clip(sm[:, 1], 0.0, 1.0)
    except Exception:
        s_ser = pd.Series(y_s).rolling(window=max(5, len(y_s) // 10), min_periods=1, center=True).mean()
        return x_s, np.clip(s_ser.to_numpy(), 0.0, 1.0)


def _plot_nrf_vs_sss_curves(
    df: pd.DataFrame,
    col_specs: List[Tuple[str, str, str, str]],
    subtitle: str,
    out_path: str,
    frac: float = 0.35,
    n_taxa: Optional[int] = None,
) -> None:
    """
    Generates a 2-panel scatter + LOWESS smoothed curve figure matching
    `nrf_vs_sss_curves.png` (Matsui & Iwasaki 2020 Fig. 3a replication) with log-scale SSS:
    - Left panel: Topological Error (nRF) vs Sequence Similarity Score (SSS, log scale)
    - Right panel: Topological Accuracy (1 - nRF) vs Sequence Similarity Score (SSS, log scale)

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
    n_taxa : Optional[int], default=None
        If provided, filters data to this specific taxon count N.
    """
    if "sss_mean" not in df.columns or df["sss_mean"].dropna().empty:
        return

    sub_df = df[df["n_taxa"] == n_taxa].copy() if n_taxa is not None else df.copy()
    valid_specs = [spec for spec in col_specs if spec[0] in sub_df.columns]
    if not valid_specs or sub_df.empty:
        return

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 7), dpi=300)

    all_pos_x = []

    for col_name, label, color, lstyle in valid_specs:
        sub = sub_df[sub_df["sss_mean"].notna() & sub_df[col_name].notna()].copy()
        sub = sub[sub["sss_mean"] > 0]
        if len(sub) < 5:
            continue
        sub = sub.sort_values("sss_mean")
        x = sub["sss_mean"].to_numpy(dtype=float)
        y_err = sub[col_name].to_numpy(dtype=float)
        y_acc = 1.0 - y_err

        all_pos_x.extend(x.tolist())

        # Semitransparent scatter points
        ax1.scatter(x, y_err, color=color, alpha=0.18, s=16, edgecolors="none")
        ax2.scatter(x, y_acc, color=color, alpha=0.18, s=16, edgecolors="none")

        # LOWESS smooth curve in log space
        lx_err, ly_err = _compute_lowess_curve(x, y_err, frac=frac, is_log_x=True)
        lx_acc, ly_acc = _compute_lowess_curve(x, y_acc, frac=frac, is_log_x=True)

        if len(lx_err) > 0:
            ax1.plot(lx_err, ly_err, color=color, linestyle=lstyle, linewidth=2.5, label=label)
        if len(lx_acc) > 0:
            ax2.plot(lx_acc, ly_acc, color=color, linestyle=lstyle, linewidth=2.5, label=label)

    if not all_pos_x:
        plt.close()
        return

    min_x = max(1e-3, float(np.min(all_pos_x)) * 0.8)
    max_x = min(1.05, max(0.04, float(np.max(all_pos_x)) * 1.15))

    taxa_suffix = f" (Taxa N = {n_taxa})" if n_taxa is not None else ""
    panel_configs = [
        (
            ax1,
            f"Topological Error (nRF) vs Sequence Similarity Score (SSS){taxa_suffix}\n[{subtitle}]",
            "Normalized RF Distance (Lower is Better)",
            "upper right",
        ),
        (
            ax2,
            f"Topological Accuracy (1 - nRF) vs Sequence Similarity Score (SSS){taxa_suffix}\n[Matsui & Iwasaki (2020) Fig. 3a Replication — {subtitle}]",
            "Correct Topology Ratio (1 - nRF, Higher is Better)",
            "lower right",
        ),
    ]

    for ax, title, ylabel, leg_loc in panel_configs:
        ax.set_xscale("log")
        ax.set_title(title, fontsize=11.5, fontweight="bold", pad=12)
        ax.set_xlabel(r"Sequence Similarity Score: $\mathrm{SSS}$ ($\bar{w}$, log scale)", fontsize=11, fontweight="bold")
        ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, which="both", linestyle="--", alpha=0.4)
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
    frac: float = 0.35,
) -> None:
    """
    Generates a 4-panel figure stratified by Taxon Count N (N = 20, 32, 64, 128),
    plotting Topological Error (nRF) vs SSS on a logarithmic horizontal axis
    with scatter points and LOWESS curves.

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
    frac : float, default=0.35
        LOWESS smoothing span parameter (unitless).
    """
    if "sss_mean" not in df.columns or "n_taxa" not in df.columns:
        return

    taxa_list = sorted(df["n_taxa"].unique())
    valid_specs = [spec for spec in col_specs if spec[0] in df.columns]
    if not valid_specs or not taxa_list:
        return

    fig, axes = plt.subplots(1, len(taxa_list), figsize=(4.8 * len(taxa_list), 5.5), dpi=300, sharey=True)
    if len(taxa_list) == 1:
        axes = [axes]

    for ax, n_val in zip(axes, taxa_list):
        sub_n = df[df["n_taxa"] == n_val]
        ax.set_xscale("log")
        n_pos_x = []

        for col_name, label, color, lstyle in valid_specs:
            sub = sub_n[sub_n["sss_mean"].notna() & sub_n[col_name].notna()].copy()
            sub = sub[sub["sss_mean"] > 0]
            if len(sub) < 5:
                continue
            sub = sub.sort_values("sss_mean")
            x = sub["sss_mean"].to_numpy(dtype=float)
            y_err = sub[col_name].to_numpy(dtype=float)
            n_pos_x.extend(x.tolist())

            ax.scatter(x, y_err, color=color, alpha=0.18, s=16, edgecolors="none")
            lx_err, ly_err = _compute_lowess_curve(x, y_err, frac=frac, is_log_x=True)
            if len(lx_err) > 0:
                ax.plot(lx_err, ly_err, color=color, linestyle=lstyle, linewidth=2.3, label=label)

        if n_pos_x:
            min_x = max(1e-3, float(np.min(n_pos_x)) * 0.8)
            max_x = min(1.05, max(0.02, float(np.max(n_pos_x)) * 1.15))
            ax.set_xlim(min_x, max_x)

        ax.set_title(f"Taxa N = {n_val}", fontsize=12.5, fontweight="bold", pad=10)
        ax.set_xlabel(r"SSS ($\bar{w}$, log scale)", fontsize=11, fontweight="bold")
        ax.set_ylim(-0.02, 1.0)
        ax.grid(True, which="both", linestyle="--", alpha=0.4)

    axes[0].set_ylabel("Normalized RF Distance (nRF, Lower is Better)", fontsize=11, fontweight="bold")
    axes[-1].legend(bbox_to_anchor=(1.02, 1.0), loc="upper left", fontsize=8.5, framealpha=0.9)
    fig.suptitle(f"Topological Error (nRF) vs SSS Stratified by Taxa N (Log Scale) — {subtitle}", fontsize=13.5, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved figure: {out_path}")


def _plot_nrf_vs_sss_single_taxa(
    df: pd.DataFrame,
    col_specs: List[Tuple[str, str, str, str]],
    subtitle: str,
    n_taxa: int,
    out_path: str,
    frac: float = 0.35,
) -> None:
    """
    Generates a dedicated single-panel figure for a specific Taxon Count N,
    plotting Topological Error (nRF) vs SSS on a logarithmic horizontal axis
    with scatter points and LOWESS curves for all methods.

    Parameters
    ----------
    df : pd.DataFrame
        Benchmark summary DataFrame containing `n_taxa`, `sss_mean`, and `*_nrf` columns.
    col_specs : List[Tuple[str, str, str, str]]
        List of (column_name, label, color, linestyle) specifications.
    subtitle : str
        Subtitle describing the distance matrix type.
    n_taxa : int
        The taxon count N to plot.
    out_path : str
        Output PNG file path.
    frac : float, default=0.35
        LOWESS smoothing span parameter.
    """
    if "sss_mean" not in df.columns or "n_taxa" not in df.columns:
        return

    sub_n = df[df["n_taxa"] == n_taxa].copy()
    valid_specs = [spec for spec in col_specs if spec[0] in sub_n.columns]
    if not valid_specs or sub_n.empty:
        return

    fig, ax = plt.subplots(figsize=(9, 6.2), dpi=300)
    ax.set_xscale("log")
    n_pos_x = []

    for col_name, label, color, lstyle in valid_specs:
        sub = sub_n[sub_n["sss_mean"].notna() & sub_n[col_name].notna()].copy()
        sub = sub[sub["sss_mean"] > 0]
        if len(sub) < 5:
            continue
        sub = sub.sort_values("sss_mean")
        x = sub["sss_mean"].to_numpy(dtype=float)
        y_err = sub[col_name].to_numpy(dtype=float)
        n_pos_x.extend(x.tolist())

        ax.scatter(x, y_err, color=color, alpha=0.20, s=18, edgecolors="none")
        lx_err, ly_err = _compute_lowess_curve(x, y_err, frac=frac, is_log_x=True)
        if len(lx_err) > 0:
            ax.plot(lx_err, ly_err, color=color, linestyle=lstyle, linewidth=2.5, label=label)

    if not n_pos_x:
        plt.close()
        return

    min_x = max(1e-3, float(np.min(n_pos_x)) * 0.8)
    max_x = min(1.05, max(0.02, float(np.max(n_pos_x)) * 1.15))
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(-0.02, 1.0)
    ax.grid(True, which="both", linestyle="--", alpha=0.4)

    ax.set_title(f"Topological Error (nRF) vs SSS (Log Scale) [Taxa N = {n_taxa}]\n{subtitle}", fontsize=12.5, fontweight="bold", pad=12)
    ax.set_xlabel(r"Sequence Similarity Score: $\mathrm{SSS}$ ($\bar{w}$, log scale)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Normalized RF Distance (nRF, Lower is Better)", fontsize=11, fontweight="bold")
    ax.legend(bbox_to_anchor=(1.02, 1.0), loc="upper left", fontsize=9, framealpha=0.9)

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

    # 2b. Standard NJ: Impact of Distance Matrix Construction Methods on nRF
    std_nj_matrix_specs = [
        ("true_Standard_NJ_nrf", "Standard NJ: True Patristic Distance", "#2ca02c", "*"),
        ("blosum62_Standard_NJ_nrf", "Standard NJ: BLOSUM62 + Poisson", "#1f77b4", "s"),
        ("wag_gamma_Standard_NJ_nrf", "Standard NJ: WAG + Γ ML", "#ff7f0e", "^"),
        ("sss_linear_Standard_NJ_nrf", "Standard NJ: MMseqs2 SSS Linear (d = 1 - w)", "#9467bd", "P"),
        ("sss_log_Standard_NJ_nrf", "Standard NJ: MMseqs2 SSS Log (d = -ln(w))", "#17becf", "D"),
        ("MSA_ML_nrf", "MSA + ML (FastTree Ref)", "#d62728", "o"),
    ]
    _plot_metric_by_taxa_and_divergence(
        df=df,
        col_specs=std_nj_matrix_specs,
        ylabel="Normalized Robinson-Foulds Distance (nRF)",
        suptitle="Standard NJ: Impact of Distance Matrix Construction on Topological Error (nRF ± SD)",
        out_path=os.path.join(outdir, "standard_nj_distance_matrix_comparison.png"),
    )
    _plot_metric_by_taxa_and_divergence(
        df=df,
        col_specs=std_nj_matrix_specs,
        ylabel="Normalized Robinson-Foulds Distance (nRF)",
        suptitle="Standard NJ: Impact of Distance Matrix Construction on Topological Error (nRF ± SD)",
        out_path=os.path.join(outdir, "nrf_standard_nj_by_matrix_vs_taxa.png"),
    )

    # Pure distance matrix specs without FastTree ref
    std_nj_pure_matrix_specs = [
        ("true_Standard_NJ_nrf", "True Patristic Distance", "#2ca02c", "*"),
        ("blosum62_Standard_NJ_nrf", "BLOSUM62 + Poisson", "#1f77b4", "s"),
        ("wag_gamma_Standard_NJ_nrf", "WAG + Γ ML", "#ff7f0e", "^"),
        ("sss_linear_Standard_NJ_nrf", "MMseqs2 SSS Linear (d = 1 - w)", "#9467bd", "P"),
        ("sss_log_Standard_NJ_nrf", "MMseqs2 SSS Log (d = -ln(w))", "#17becf", "D"),
    ]
    _plot_metric_by_taxa_and_divergence(
        df=df,
        col_specs=std_nj_pure_matrix_specs,
        ylabel="Normalized Robinson-Foulds Distance (nRF)",
        suptitle="Standard NJ: Distance Matrix Comparison (Pure Distance, nRF ± SD)",
        out_path=os.path.join(outdir, "standard_nj_pure_matrix_comparison.png"),
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

    taxa_list = sorted(df["n_taxa"].unique()) if "n_taxa" in df.columns else []

    # 4. nRF & Accuracy (1 - nRF) vs SSS Scatter + LOWESS Curves (Matsui & Iwasaki 2020 Fig. 3a style, Log-Scale)
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
        # Overall 2-panel curves (all Taxa combined, log x-scale)
        _plot_nrf_vs_sss_curves(
            df=df,
            col_specs=sss_curve_specs,
            subtitle=m_title,
            out_path=os.path.join(outdir, f"nrf_vs_sss_curves_{m_type}.png"),
            frac=0.35,
        )
        # 4-panel subplots stratified by Taxa N (log x-scale)
        _plot_nrf_vs_sss_by_taxa(
            df=df,
            col_specs=sss_curve_specs,
            subtitle=m_title,
            out_path=os.path.join(outdir, f"nrf_vs_sss_by_taxa_{m_type}.png"),
            frac=0.35,
        )
        # Dedicated per-N single-panel and 2-panel curve plots
        for n_val in taxa_list:
            _plot_nrf_vs_sss_single_taxa(
                df=df,
                col_specs=sss_curve_specs,
                subtitle=m_title,
                n_taxa=n_val,
                out_path=os.path.join(outdir, f"nrf_vs_sss_N{n_val}_{m_type}.png"),
                frac=0.35,
            )
            _plot_nrf_vs_sss_curves(
                df=df,
                col_specs=sss_curve_specs,
                subtitle=m_title,
                out_path=os.path.join(outdir, f"nrf_vs_sss_curves_N{n_val}_{m_type}.png"),
                frac=0.35,
                n_taxa=n_val,
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
    for n_val in taxa_list:
        _plot_nrf_vs_sss_curves(
            df=df,
            col_specs=main_sss_specs,
            subtitle="BLOSUM62 & WAG+Γ vs MSA+ML",
            out_path=os.path.join(outdir, f"nrf_vs_sss_curves_N{n_val}.png"),
            frac=0.35,
            n_taxa=n_val,
        )
        _plot_nrf_vs_sss_single_taxa(
            df=df,
            col_specs=main_sss_specs,
            subtitle="BLOSUM62 & WAG+Γ vs MSA+ML",
            n_taxa=n_val,
            out_path=os.path.join(outdir, f"nrf_vs_sss_N{n_val}.png"),
            frac=0.35,
        )


def _plot_standard_nj_cross_indel_comparison(
    df_005: pd.DataFrame,
    df_010: pd.DataFrame,
    out_paths: List[str],
) -> None:
    """
    Generates a 2-row x 3-column error bar plot comparing how distance matrix construction
    impacts Standard NJ performance across two INDELible indel rates (0.05 and 0.10).
    """
    col_specs = [
        ("true_Standard_NJ_nrf", "Standard NJ: True Patristic Distance", "#2ca02c", "*"),
        ("blosum62_Standard_NJ_nrf", "Standard NJ: BLOSUM62 + Poisson", "#1f77b4", "s"),
        ("wag_gamma_Standard_NJ_nrf", "Standard NJ: WAG + Γ ML", "#ff7f0e", "^"),
        ("sss_linear_Standard_NJ_nrf", "Standard NJ: MMseqs2 SSS Linear (d = 1 - w)", "#9467bd", "P"),
        ("sss_log_Standard_NJ_nrf", "Standard NJ: MMseqs2 SSS Log (d = -ln(w))", "#17becf", "D"),
        ("MSA_ML_nrf", "MSA + ML (FastTree Ref)", "#d62728", "o"),
    ]
    divergences = [0.5, 1.0, 2.0]
    taxa_list = [20, 32, 64, 128]
    n_methods = len(col_specs)
    dodge_offsets = np.linspace(-0.30, 0.30, n_methods)

    fig, axes = plt.subplots(2, 3, figsize=(18, 10.5), dpi=300, sharey=False)

    row_data = [
        ("Indel Rate = 0.05", df_005),
        ("Indel Rate = 0.10", df_010),
    ]

    for r_idx, (indel_label, df) in enumerate(row_data):
        for c_idx, div in enumerate(divergences):
            ax = axes[r_idx, c_idx]
            sub_df = df[df["divergence"] == div]
            max_y_val = 0.02

            for idx, (col_name, label, color, marker) in enumerate(col_specs):
                x_pts, y_means, y_sds = [], [], []
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
                        markersize=7.5,
                        capsize=3.5,
                        capthick=1.2,
                        elinewidth=1.2,
                        linestyle="none",
                        alpha=0.9,
                        markeredgecolor="black",
                        markeredgewidth=0.5,
                    )

            sss_str = ""
            if "sss_mean" in sub_df.columns and sub_df["sss_mean"].notna().any():
                sss_str = f" (Mean SSS = {sub_df['sss_mean'].mean():.3f})"

            ax.set_title(f"{indel_label} | Divergence D = {div}{sss_str}", fontsize=12, fontweight="bold", pad=8)
            ax.set_xlabel("Number of Taxa (N)", fontsize=11, fontweight="bold")
            ax.set_xticks(range(len(taxa_list)))
            ax.set_xticklabels([f"N = {n}" for n in taxa_list], fontsize=10, fontweight="bold")

            y_margin = max_y_val * 0.14
            ax.set_ylim(-max(0.003, y_margin * 0.25), max_y_val + y_margin)
            if c_idx == 0:
                ax.set_ylabel("Normalized RF Distance (nRF)", fontsize=11, fontweight="bold")

            for sep in np.arange(0.5, len(taxa_list) - 0.5, 1.0):
                ax.axvline(sep, color="gray", linestyle=":", alpha=0.35, zorder=0)

            ax.grid(True, axis="y", ls="--", alpha=0.5, zorder=0)

    axes[0, -1].legend(
        bbox_to_anchor=(1.03, 1.0),
        loc="upper left",
        fontsize=9.5,
        frameon=True,
        edgecolor="gray",
        title="Distance Matrix Construction",
        title_fontsize=10.5,
    )
    fig.suptitle(
        "Standard NJ: Impact of Distance Matrix Construction on Topological Error (nRF ± SD)\n"
        "[INDELible Zipfian Indel Rates 0.05 vs 0.10]",
        fontsize=14,
        fontweight="bold",
        y=1.01,
    )
    plt.tight_layout()
    for op in out_paths:
        os.makedirs(os.path.dirname(op), exist_ok=True)
        plt.savefig(op, dpi=300, bbox_inches="tight")
        print(f"[OK] Saved cross-indel figure: {op}")
    plt.close()


def generate_cross_indel_comparison(
    dir_005: str = "results_hpc_indel0.05/main",
    dir_010: str = "results_hpc_indel0.1/main",
) -> None:
    """Checks for both indel summary CSVs and produces the unified 2-row comparison plot."""
    path_005 = os.path.join(dir_005, "hpc_benchmark_summary.csv")
    path_010 = os.path.join(dir_010, "hpc_benchmark_summary.csv")
    if os.path.exists(path_005) and os.path.exists(path_010):
        df_005 = pd.read_csv(path_005)
        df_010 = pd.read_csv(path_010)
        out_paths = [
            os.path.join(dir_005, "standard_nj_distance_matrix_cross_indel.png"),
            os.path.join(dir_010, "standard_nj_distance_matrix_cross_indel.png"),
            os.path.join("results", "standard_nj_distance_matrix_cross_indel.png"),
        ]
        _plot_standard_nj_cross_indel_comparison(df_005, df_010, out_paths)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge HPC shard CSVs and generate plots.")
    parser.add_argument("--outdir", type=str, default="results_hpc", help="HPC results directory")
    args = parser.parse_args()

    df_merged = merge_shards(outdir=args.outdir)
    generate_all_hpc_plots(df_merged, outdir=args.outdir)
    generate_cross_indel_comparison()

