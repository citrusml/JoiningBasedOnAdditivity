"""
Main Benchmark Script for Additivity-Based Phylogenetic Tree Reconstruction.

Evaluates AddNJ (Method 1), Add-NNI/SPR/Hybrid (Method 2), Standard NJ, BioNJ,
and FastME BME-NNI/SPR across synthetic WAG+Gamma protein sequences under varying
divergence scales (0.5, 1.0, 2.0) with global alignment.
"""

from typing import List, Dict, Any, Tuple
import os
import sys
import io
import time
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Bio.Phylo as Phylo

from src.tree_simulator import generate_backward_yule_tree, compute_patristic_distance_matrix
from src.sequence_simulator import WAGModel, simulate_wag_gamma_sequences, compute_pairwise_poisson_distance_matrix
from src.tree_inference import matrix_to_nj_newick
from src.evaluation import evaluate_tree_rf, FourPointCore
from src.add_nj import add_nj, add_nj2, add_nj3
from src.fastme_wrapper import run_fastme
from src.topology_search import TreeScorer, TreeDistanceScorer, hill_climb_search


def run_benchmark(
    n_taxa: int = 20,
    seq_length: int = 500,
    divergences: List[float] = [0.5, 1.0, 2.0],
    replicates: int = 10,
    alpha: float = 1.0,
    outdir: str = "results",
    seed: int = 42,
) -> pd.DataFrame:
    """
    Executes the comprehensive phylogenetic reconstruction benchmark.

    Parameters
    ----------
    n_taxa : int, default=20
        Number of taxa (unitless).
    seq_length : int, default=500
        Amino acid sequence length (amino acids).
    divergences : List[float], default=[0.5, 1.0, 2.0]
        Tree divergence scale levels (substitutions/site).
    replicates : int, default=10
        Number of replicates per divergence level (unitless).
    alpha : float, default=1.0
        Gamma shape parameter (unitless).
    outdir : str, default='results'
        Output directory for results (path).
    seed : int, default=42
        Base pseudo-random seed (unitless).

    Returns
    -------
    pd.DataFrame
        DataFrame containing all experimental records.
    """
    os.makedirs(outdir, exist_ok=True)
    wag = WAGModel(wag_path="src/wag.dat", alpha=alpha)
    core = FourPointCore(n_taxa=n_taxa)

    all_records: List[Dict[str, Any]] = []

    print("=" * 80)
    print("PHYLOMETRICS BENCHMARK: ADDITIVITY-BASED METHODS VS BASELINES")
    print(f"Taxa: {n_taxa} | Length: {seq_length} aa | Divergences: {divergences} | Replicates: {replicates}")
    print("=" * 80)

    total_runs = len(divergences) * replicates
    run_idx = 0

    for div in divergences:
        print(f"\n>>> Running Divergence Scale: {div} substitutions/site <<<")
        for rep in range(1, replicates + 1):
            run_idx += 1
            rep_seed = seed + int(div * 1000) + rep * 100

            # 1. Simulate true tree and patristic distances
            tree, taxa_names = generate_backward_yule_tree(
                n_taxa=n_taxa, divergence=div, seed=rep_seed
            )
            true_newick = tree.to_newick() + ";"
            D_true = compute_patristic_distance_matrix(tree, taxa_names)

            # 2. Simulate WAG+Gamma sequences
            seqs = simulate_wag_gamma_sequences(
                root=tree,
                seq_length=seq_length,
                alpha=alpha,
                model=wag,
                seed=rep_seed + 1,
            )

            # 3. Global alignment & Poisson distance matrix
            t_align = time.time()
            D_obs = compute_pairwise_poisson_distance_matrix(seqs, taxa_names, mode="global")
            align_time = time.time() - t_align

            # Initial matrix defect
            delta_add_init, _ = core.compute_additivity_defect(D_obs)
            scorer = TreeScorer(D_obs, taxa_names)
            dist_scorer = TreeDistanceScorer(D_obs, taxa_names, mode="ols_rmse")

            base_info = {
                "run_idx": run_idx,
                "replicate": rep,
                "divergence": div,
                "n_taxa": n_taxa,
                "seq_length": seq_length,
                "delta_add_obs": float(delta_add_init),
                "align_time_sec": float(align_time),
            }

            # Methods dictionary: method_name -> (inferred_newick, runtime_sec, penalty_score)
            method_outputs: Dict[str, Tuple[str, float, float]] = {}

            # 1. Standard NJ
            t0 = time.time()
            nj_nwk = matrix_to_nj_newick(D_obs, taxa_names)
            t_nj = time.time() - t0
            s_nj = scorer.score_tree(Phylo.read(io.StringIO(nj_nwk), "newick"))
            method_outputs["Standard_NJ"] = (nj_nwk, t_nj, s_nj)

            # 2. BioNJ
            t0 = time.time()
            bionj_nwk = run_fastme(D_obs, taxa_names, method="bionj")
            t_bionj = time.time() - t0
            s_bionj = scorer.score_tree(Phylo.read(io.StringIO(bionj_nwk), "newick"))
            method_outputs["BioNJ"] = (bionj_nwk, t_bionj, s_bionj)

            # 3. FastME-NNI
            t0 = time.time()
            fnni_nwk = run_fastme(D_obs, taxa_names, method="nni")
            t_fnni = time.time() - t0
            s_fnni = scorer.score_tree(Phylo.read(io.StringIO(fnni_nwk), "newick"))
            method_outputs["FastME_NNI"] = (fnni_nwk, t_fnni, s_fnni)

            # 4. FastME-SPR
            t0 = time.time()
            fspr_nwk = run_fastme(D_obs, taxa_names, method="spr")
            t_fspr = time.time() - t0
            s_fspr = scorer.score_tree(Phylo.read(io.StringIO(fspr_nwk), "newick"))
            method_outputs["FastME_SPR"] = (fspr_nwk, t_fspr, s_fspr)

            # 5. AddNJ (Method 1)
            t0 = time.time()
            addnj_nwk = add_nj(D_obs, taxa_names)
            t_addnj = time.time() - t0
            s_addnj = scorer.score_tree(Phylo.read(io.StringIO(addnj_nwk), "newick"))
            method_outputs["AddNJ"] = (addnj_nwk, t_addnj, s_addnj)

            # 6. AddNJ2 (Method 1 Improved: Sister-Pair Quartet Consistency)
            t0 = time.time()
            addnj2_nwk = add_nj2(D_obs, taxa_names)
            t_addnj2 = time.time() - t0
            s_addnj2 = scorer.score_tree(Phylo.read(io.StringIO(addnj2_nwk), "newick"))
            method_outputs["AddNJ2"] = (addnj2_nwk, t_addnj2, s_addnj2)

            # 6. NJ + Add-NNI (Method 2)
            t0 = time.time()
            nj_nni_nwk, s_nj_nni = hill_climb_search(nj_nwk, scorer, mode="nni")
            t_nj_nni = time.time() - t0
            method_outputs["NJ_Add_NNI"] = (nj_nni_nwk, t_nj_nni, s_nj_nni)

            # 7. NJ + Add-SPR (Method 2)
            t0 = time.time()
            nj_spr_nwk, s_nj_spr = hill_climb_search(nj_nwk, scorer, mode="spr")
            t_nj_spr = time.time() - t0
            method_outputs["NJ_Add_SPR"] = (nj_spr_nwk, t_nj_spr, s_nj_spr)

            # 8. NJ + Add-Hybrid (Method 2)
            t0 = time.time()
            nj_hyb_nwk, s_nj_hyb = hill_climb_search(nj_nwk, scorer, mode="hybrid")
            t_nj_hyb = time.time() - t0
            method_outputs["NJ_Add_Hybrid"] = (nj_hyb_nwk, t_nj_hyb, s_nj_hyb)

            # 9. BioNJ + Add-NNI
            t0 = time.time()
            bionj_nni_nwk, s_bionj_nni = hill_climb_search(bionj_nwk, scorer, mode="nni")
            t_bionj_nni = time.time() - t0
            method_outputs["BioNJ_Add_NNI"] = (bionj_nni_nwk, t_bionj_nni, s_bionj_nni)

            # 10. BioNJ + Add-SPR
            t0 = time.time()
            bionj_spr_nwk, s_bionj_spr = hill_climb_search(bionj_nwk, scorer, mode="spr")
            t_bionj_spr = time.time() - t0
            method_outputs["BioNJ_Add_SPR"] = (bionj_spr_nwk, t_bionj_spr, s_bionj_spr)

            # 11. BioNJ + Add-Hybrid
            t0 = time.time()
            bionj_hyb_nwk, s_bionj_hyb = hill_climb_search(bionj_nwk, scorer, mode="hybrid")
            t_bionj_hyb = time.time() - t0
            method_outputs["BioNJ_Add_Hybrid"] = (bionj_hyb_nwk, t_bionj_hyb, s_bionj_hyb)

            # 12. AddNJ + Add-Hybrid
            t0 = time.time()
            addnj_hyb_nwk, s_addnj_hyb = hill_climb_search(addnj_nwk, scorer, mode="hybrid")
            t_addnj_hyb = time.time() - t0
            method_outputs["AddNJ_Add_Hybrid"] = (addnj_hyb_nwk, t_addnj_hyb, s_addnj_hyb)

            # 13. AddNJ2 + Add-Hybrid
            t0 = time.time()
            addnj2_hyb_nwk, s_addnj2_hyb = hill_climb_search(addnj2_nwk, scorer, mode="hybrid")
            t_addnj2_hyb = time.time() - t0
            method_outputs["AddNJ2_Add_Hybrid"] = (addnj2_hyb_nwk, t_addnj2_hyb, s_addnj2_hyb)

            # 14. AddNJ3 (Method 1 Improved 2: Patristic-Matrix Distance Discrepancy)
            t0 = time.time()
            addnj3_nwk = add_nj3(D_obs, taxa_names)
            t_addnj3 = time.time() - t0
            s_addnj3 = scorer.score_newick(addnj3_nwk)
            method_outputs["AddNJ3"] = (addnj3_nwk, t_addnj3, s_addnj3)

            # 15. NJ + Add-NNI2 (Method 2 Improved: Distance Discrepancy)
            t0 = time.time()
            nj_nni2_nwk, _ = hill_climb_search(nj_nwk, dist_scorer, mode="nni")
            t_nj_nni2 = time.time() - t0
            s_nj_nni2 = scorer.score_newick(nj_nni2_nwk)
            method_outputs["NJ_Add_NNI2"] = (nj_nni2_nwk, t_nj_nni2, s_nj_nni2)

            # 16. NJ + Add-SPR2 (Method 2 Improved: Distance Discrepancy)
            t0 = time.time()
            nj_spr2_nwk, _ = hill_climb_search(nj_nwk, dist_scorer, mode="spr")
            t_nj_spr2 = time.time() - t0
            s_nj_spr2 = scorer.score_newick(nj_spr2_nwk)
            method_outputs["NJ_Add_SPR2"] = (nj_spr2_nwk, t_nj_spr2, s_nj_spr2)

            # 17. NJ + Add-Hybrid2 (Method 2 Improved: Distance Discrepancy)
            t0 = time.time()
            nj_hyb2_nwk, _ = hill_climb_search(nj_nwk, dist_scorer, mode="hybrid")
            t_nj_hyb2 = time.time() - t0
            s_nj_hyb2 = scorer.score_newick(nj_hyb2_nwk)
            method_outputs["NJ_Add_Hybrid2"] = (nj_hyb2_nwk, t_nj_hyb2, s_nj_hyb2)

            # 18. BioNJ + Add-Hybrid2
            t0 = time.time()
            bionj_hyb2_nwk, _ = hill_climb_search(bionj_nwk, dist_scorer, mode="hybrid")
            t_bionj_hyb2 = time.time() - t0
            s_bionj_hyb2 = scorer.score_newick(bionj_hyb2_nwk)
            method_outputs["BioNJ_Add_Hybrid2"] = (bionj_hyb2_nwk, t_bionj_hyb2, s_bionj_hyb2)

            # 19. AddNJ3 + Add-Hybrid2
            t0 = time.time()
            addnj3_hyb2_nwk, _ = hill_climb_search(addnj3_nwk, dist_scorer, mode="hybrid")
            t_addnj3_hyb2 = time.time() - t0
            s_addnj3_hyb2 = scorer.score_newick(addnj3_hyb2_nwk)
            method_outputs["AddNJ3_Add_Hybrid2"] = (addnj3_hyb2_nwk, t_addnj3_hyb2, s_addnj3_hyb2)

            # Evaluate nRF for each method
            row = dict(base_info)
            for m_name, (inferred_str, runtime, penalty) in method_outputs.items():
                eval_res = evaluate_tree_rf(inferred_str, true_newick, taxa_names)
                row[f"{m_name}_nrf"] = float(eval_res["nrf"])
                row[f"{m_name}_rf"] = float(eval_res["rf"])
                row[f"{m_name}_time"] = float(runtime)
                row[f"{m_name}_penalty"] = float(penalty)

            all_records.append(row)
            print(
                f"  [D={div}] Rep {rep:2d}/{replicates}: "
                f"NJ={row['Standard_NJ_nrf']:.3f} | "
                f"BioNJ={row['BioNJ_nrf']:.3f} | "
                f"FastME-SPR={row['FastME_SPR_nrf']:.3f} | "
                f"AddNJ2={row['AddNJ2_nrf']:.3f} | "
                f"AddNJ3={row['AddNJ3_nrf']:.3f} | "
                f"NJ+Add-Hyb={row['NJ_Add_Hybrid_nrf']:.3f} | "
                f"NJ+Add-Hyb2={row['NJ_Add_Hybrid2_nrf']:.3f}"
            )

    df = pd.DataFrame(all_records)
    csv_path = os.path.join(outdir, "benchmark_summary.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n[OK] Benchmark completed. Saved raw data to: {csv_path}")
    return df


def generate_plots(df: pd.DataFrame, outdir: str = "results") -> None:
    """
    Generates publication-quality comparison charts at 300 DPI.
    """
    methods = [
        "Standard_NJ",
        "BioNJ",
        "FastME_NNI",
        "FastME_SPR",
        "AddNJ",
        "AddNJ2",
        "AddNJ3",
        "NJ_Add_NNI",
        "NJ_Add_SPR",
        "NJ_Add_Hybrid",
        "NJ_Add_NNI2",
        "NJ_Add_SPR2",
        "NJ_Add_Hybrid2",
        "BioNJ_Add_NNI",
        "BioNJ_Add_SPR",
        "BioNJ_Add_Hybrid",
        "BioNJ_Add_Hybrid2",
        "AddNJ_Add_Hybrid",
        "AddNJ2_Add_Hybrid",
        "AddNJ3_Add_Hybrid2",
    ]

    labels = [
        "Standard NJ",
        "BioNJ",
        "FastME (NNI)",
        "FastME (SPR)",
        "AddNJ (Method 1)",
        "AddNJ2 (Improved M1)",
        "AddNJ3 (Dist RMSE)",
        "NJ + Add-NNI",
        "NJ + Add-SPR",
        "NJ + Add-Hybrid",
        "NJ + Add-NNI2",
        "NJ + Add-SPR2",
        "NJ + Add-Hybrid2",
        "BioNJ + Add-NNI",
        "BioNJ + Add-SPR",
        "BioNJ + Add-Hybrid",
        "BioNJ + Add-Hybrid2",
        "AddNJ + Add-Hybrid",
        "AddNJ2 + Add-Hybrid",
        "AddNJ3 + Add-Hybrid2",
    ]

    divergences = sorted(df["divergence"].unique())

    # 1. Figure 1: Mean nRF vs Divergence (Error Bar Plot: Full Scale & Zoomed Scale)
    markers = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "p", "*", "h", "H", "d", "8", "o", "s", "^", "D", "v"]
    colors = plt.cm.tab20(np.linspace(0, 1, len(methods)))

    x = np.arange(len(divergences))
    n_methods = len(methods)
    dodge_offsets = np.linspace(-0.38, 0.38, n_methods)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17.0, 6.6), dpi=300)
    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")

    # Panel 1: Full scale including all methods
    ax1.grid(True, linestyle="--", alpha=0.5)
    lines = []
    for idx, (m, lbl) in enumerate(zip(methods, labels)):
        means = [df[df["divergence"] == d][f"{m}_nrf"].mean() for d in divergences]
        sems = [df[df["divergence"] == d][f"{m}_nrf"].sem() for d in divergences]
        eb = ax1.errorbar(
            x + dodge_offsets[idx],
            means,
            yerr=sems,
            fmt=markers[idx % len(markers)],
            color=colors[idx],
            label=lbl,
            markersize=8,
            capsize=3.5,
            capthick=1.2,
            elinewidth=1.2,
            markeredgecolor="black",
            markeredgewidth=0.6,
            linestyle="none",
            alpha=0.95,
        )
        lines.append(eb)

    ax1.set_xlabel("Tree Divergence Scale (substitutions/site)", fontsize=12, fontweight="bold", labelpad=8)
    ax1.set_ylabel("Normalized Robinson-Foulds Distance (nRF)", fontsize=12, fontweight="bold", labelpad=8)
    ax1.set_title("(a) All Evaluated Methods (Full Scale)", fontsize=13, fontweight="bold", pad=10)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"D = {d}" for d in divergences], fontsize=11, fontweight="bold")
    ax1.set_ylim(-0.02, 0.85)
    for sep in [0.5, 1.5]:
        ax1.axvline(sep, color="gray", linestyle=":", alpha=0.3)

    # Panel 2: Zoomed-in scale for competitive methods (excluding naive AddNJ outlier, keeping same dodge positions)
    ax2.grid(True, linestyle="--", alpha=0.5)
    for idx, (m, lbl) in enumerate(zip(methods, labels)):
        if m == "AddNJ":
            continue
        means = [df[df["divergence"] == d][f"{m}_nrf"].mean() for d in divergences]
        sems = [df[df["divergence"] == d][f"{m}_nrf"].sem() for d in divergences]
        ax2.errorbar(
            x + dodge_offsets[idx],
            means,
            yerr=sems,
            fmt=markers[idx % len(markers)],
            color=colors[idx],
            label=lbl,
            markersize=8,
            capsize=3.5,
            capthick=1.2,
            elinewidth=1.2,
            markeredgecolor="black",
            markeredgewidth=0.6,
            linestyle="none",
            alpha=0.95,
        )

    ax2.set_xlabel("Tree Divergence Scale (substitutions/site)", fontsize=12, fontweight="bold", labelpad=8)
    ax2.set_ylabel("Normalized Robinson-Foulds Distance (nRF)", fontsize=12, fontweight="bold", labelpad=8)
    ax2.set_title("(b) Competitive Methods (Zoomed Scale)", fontsize=13, fontweight="bold", pad=10)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"D = {d}" for d in divergences], fontsize=11, fontweight="bold")
    ax2.set_ylim(-0.005, 0.17)
    for sep in [0.5, 1.5]:
        ax2.axvline(sep, color="gray", linestyle=":", alpha=0.3)

    # Place legend in the dedicated right margin with explicit Left-to-Right cluster ordering
    fig.legend(
        handles=lines,
        labels=labels,
        title="Reconstruction Methods\n(Left-to-Right in cluster)",
        title_fontsize=10.0,
        loc="center left",
        bbox_to_anchor=(0.825, 0.5),
        fontsize=8.5,
        frameon=True,
        framealpha=0.95,
        edgecolor="#cccccc",
    )
    fig.suptitle("Phylogenetic Reconstruction Accuracy Across Divergence Levels (N=20, L=500 aa, 10 reps, Mean ± SEM)", fontsize=13.5, fontweight="bold", y=0.97)
    plt.subplots_adjust(left=0.07, right=0.81, top=0.88, bottom=0.12, wspace=0.22)

    fig1_path = os.path.join(outdir, "nrf_comparison.png")
    plt.savefig(fig1_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved figure: {fig1_path}")

    # 2. Figure 2: Execution Time Comparison
    plt.figure(figsize=(10, 6.5), dpi=300)
    mean_times = [df[f"{m}_time"].mean() for m in methods]
    sem_times = [df[f"{m}_time"].sem() for m in methods]

    y_pos = np.arange(len(methods))
    plt.barh(y_pos, mean_times, xerr=sem_times, color="#2b5c8f", alpha=0.85, capsize=3)
    plt.yticks(y_pos, labels, fontsize=10)
    plt.xlabel("Mean Execution Time per Replicate (seconds)", fontsize=12, fontweight="bold")
    plt.title("Computational Efficiency Benchmark (N=20 taxa)", fontsize=13, fontweight="bold")
    plt.xscale("log")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()

    fig2_path = os.path.join(outdir, "runtime_comparison.png")
    plt.savefig(fig2_path, dpi=300)
    plt.close()
    print(f"[OK] Saved figure: {fig2_path}")

    # 3. Figure 3: Four-Point Condition Defect / Tree Additivity Penalty vs Divergence
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17.0, 6.6), dpi=300)
    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")

    # Panel 1: Full scale including naive AddNJ
    ax1.grid(True, linestyle="--", alpha=0.5)
    lines_pen = []
    for idx, (m, lbl) in enumerate(zip(methods, labels)):
        means = [df[df["divergence"] == d][f"{m}_penalty"].mean() for d in divergences]
        sems = [df[df["divergence"] == d][f"{m}_penalty"].sem() for d in divergences]
        eb = ax1.errorbar(
            x + dodge_offsets[idx],
            means,
            yerr=sems,
            fmt=markers[idx % len(markers)],
            color=colors[idx],
            label=lbl,
            markersize=8,
            capsize=3.5,
            capthick=1.2,
            elinewidth=1.2,
            markeredgecolor="black",
            markeredgewidth=0.6,
            linestyle="none",
            alpha=0.95,
        )
        lines_pen.append(eb)

    ax1.set_xlabel("Tree Divergence Scale (substitutions/site)", fontsize=12, fontweight="bold", labelpad=8)
    ax1.set_ylabel(r"Four-Point Condition Defect $\mathcal{L}(T; D_{\mathrm{obs}})$", fontsize=12, fontweight="bold", labelpad=8)
    ax1.set_title("(a) All Evaluated Methods (Full Scale)", fontsize=13, fontweight="bold", pad=10)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"D = {d}" for d in divergences], fontsize=11, fontweight="bold")
    ax1.set_ylim(-0.002, max(0.06, df[[f"{m}_penalty" for m in methods]].values.max() * 1.15))
    for sep in [0.5, 1.5]:
        ax1.axvline(sep, color="gray", linestyle=":", alpha=0.3)

    # Panel 2: Zoomed-in scale for competitive methods (excluding naive AddNJ)
    ax2.grid(True, linestyle="--", alpha=0.5)
    for idx, (m, lbl) in enumerate(zip(methods, labels)):
        if m == "AddNJ":
            continue
        means = [df[df["divergence"] == d][f"{m}_penalty"].mean() for d in divergences]
        sems = [df[df["divergence"] == d][f"{m}_penalty"].sem() for d in divergences]
        ax2.errorbar(
            x + dodge_offsets[idx],
            means,
            yerr=sems,
            fmt=markers[idx % len(markers)],
            color=colors[idx],
            label=lbl,
            markersize=8,
            capsize=3.5,
            capthick=1.2,
            elinewidth=1.2,
            markeredgecolor="black",
            markeredgewidth=0.6,
            linestyle="none",
            alpha=0.95,
        )

    comp_pen_max = max(0.009, max(df[df["divergence"] == d][[f"{m}_penalty" for m in methods if m != "AddNJ"]].values.max() for d in divergences) * 1.15)
    ax2.set_xlabel("Tree Divergence Scale (substitutions/site)", fontsize=12, fontweight="bold", labelpad=8)
    ax2.set_ylabel(r"Four-Point Condition Defect $\mathcal{L}(T; D_{\mathrm{obs}})$", fontsize=12, fontweight="bold", labelpad=8)
    ax2.set_title("(b) Competitive Methods (Zoomed Scale)", fontsize=13, fontweight="bold", pad=10)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"D = {d}" for d in divergences], fontsize=11, fontweight="bold")
    ax2.set_ylim(-0.0003, comp_pen_max)
    for sep in [0.5, 1.5]:
        ax2.axvline(sep, color="gray", linestyle=":", alpha=0.3)

    fig.legend(
        handles=lines_pen,
        labels=labels,
        title="Reconstruction Methods\n(Left-to-Right in cluster)",
        title_fontsize=10.0,
        loc="center left",
        bbox_to_anchor=(0.825, 0.5),
        fontsize=8.5,
        frameon=True,
        framealpha=0.95,
        edgecolor="#cccccc",
    )
    fig.suptitle(r"Four-Point Condition Defect $\mathcal{L}(T; D_{\mathrm{obs}})$ Across Divergence Levels (N=20, L=500 aa, 10 reps, Mean $\pm$ SEM)", fontsize=13.5, fontweight="bold", y=0.97)
    plt.subplots_adjust(left=0.07, right=0.81, top=0.88, bottom=0.12, wspace=0.22)

    fig3_path = os.path.join(outdir, "penalty_comparison.png")
    plt.savefig(fig3_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved figure: {fig3_path}")

    # Also save with alternative intuitive name
    fig3_alias = os.path.join(outdir, "four_point_defect_comparison.png")
    import shutil
    shutil.copyfile(fig3_path, fig3_alias)
    print(f"[OK] Saved figure: {fig3_alias}")


def update_raw_results(df: pd.DataFrame) -> None:
    """
    Appends the summary table and raw trial records to RawResults.md.
    """
    methods = [
        ("Standard_NJ", "Standard NJ"),
        ("BioNJ", "BioNJ"),
        ("FastME_NNI", "FastME-NNI"),
        ("FastME_SPR", "FastME-SPR"),
        ("AddNJ", "AddNJ (手法1: 縮約後欠損)"),
        ("AddNJ2", "AddNJ2 (手法1改善: 姉妹群整合度)"),
        ("AddNJ3", "AddNJ3 (手法1改善2: 系統樹距離残差)"),
        ("NJ_Add_NNI", "NJ + Add-NNI (手法2)"),
        ("NJ_Add_SPR", "NJ + Add-SPR (手法2)"),
        ("NJ_Add_Hybrid", "NJ + Add-Hybrid (手法2)"),
        ("NJ_Add_NNI2", "NJ + Add-NNI2 (手法2改善: 距離残差)"),
        ("NJ_Add_SPR2", "NJ + Add-SPR2 (手法2改善: 距離残差)"),
        ("NJ_Add_Hybrid2", "NJ + Add-Hybrid2 (手法2改善: 距離残差)"),
        ("BioNJ_Add_NNI", "BioNJ + Add-NNI"),
        ("BioNJ_Add_SPR", "BioNJ + Add-SPR"),
        ("BioNJ_Add_Hybrid", "BioNJ + Add-Hybrid"),
        ("BioNJ_Add_Hybrid2", "BioNJ + Add-Hybrid2"),
        ("AddNJ_Add_Hybrid", "AddNJ + Add-Hybrid"),
        ("AddNJ2_Add_Hybrid", "AddNJ2 + Add-Hybrid"),
        ("AddNJ3_Add_Hybrid2", "AddNJ3 + Add-Hybrid2"),
    ]

    divergences = sorted(df["divergence"].unique())

    content = "\n\n## ベンチマーク集計結果（平均 ± 標準誤差）\n\n"
    content += "| 手法 | D = 0.5 (nRF) | D = 1.0 (nRF) | D = 2.0 (nRF) | 平均実行時間 (秒) |\n"
    content += "| :--- | :---: | :---: | :---: | :---: |\n"

    for m_col, m_label in methods:
        nrf_strs = []
        for d in divergences:
            sub = df[df["divergence"] == d][f"{m_col}_nrf"]
            m_val = sub.mean()
            s_val = sub.sem()
            nrf_strs.append(f"{m_val:.4f} ± {s_val:.4f}")
        avg_time = df[f"{m_col}_time"].mean()
        content += f"| **{m_label}** | {nrf_strs[0]} | {nrf_strs[1]} | {nrf_strs[2]} | {avg_time:.4f}s |\n"

    content += "\n### 4点制約違反度 $\mathcal{L}(T; D_{\mathrm{obs}})$ 集計結果（平均 ± 標準誤差）\n\n"
    content += "| 手法 | D = 0.5 (Penalty) | D = 1.0 (Penalty) | D = 2.0 (Penalty) |\n"
    content += "| :--- | :---: | :---: | :---: |\n"

    for m_col, m_label in methods:
        pen_strs = []
        for d in divergences:
            sub = df[df["divergence"] == d][f"{m_col}_penalty"]
            m_val = sub.mean()
            s_val = sub.sem()
            pen_strs.append(f"{m_val:.6f} ± {s_val:.6f}")
        content += f"| **{m_label}** | {pen_strs[0]} | {pen_strs[1]} | {pen_strs[2]} |\n"

    content += "\n### 全レプリケート生データ (抜粋: 各条件最初の3反復)\n\n"
    content += "| Rep | Div | Standard NJ | BioNJ | FastME-SPR | AddNJ2 | AddNJ3 | NJ+Add-Hyb | NJ+Add-Hyb2 | BioNJ+Add-Hyb2 |\n"
    content += "| :-: | :-: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n"

    for d in divergences:
        sub_df = df[df["divergence"] == d].head(3)
        for _, r in sub_df.iterrows():
            content += (
                f"| {int(r['replicate'])} | {r['divergence']:.1f} | "
                f"{r['Standard_NJ_nrf']:.3f} | {r['BioNJ_nrf']:.3f} | {r['FastME_SPR_nrf']:.3f} | "
                f"{r['AddNJ2_nrf']:.3f} | {r['AddNJ3_nrf']:.3f} | {r['NJ_Add_Hybrid_nrf']:.3f} | "
                f"{r['NJ_Add_Hybrid2_nrf']:.3f} | {r['BioNJ_Add_Hybrid2_nrf']:.3f} |\n"
            )

    with open("RawResults.md", "a") as f:
        f.write(content)
    print("[OK] Appended summary tables to RawResults.md")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Additivity-Based Phylometrics Benchmark.")
    parser.add_argument("--taxa", type=int, default=20, help="Number of taxa (default: 20)")
    parser.add_argument("--length", type=int, default=500, help="Sequence length (default: 500)")
    parser.add_argument("--divergences", type=float, nargs="+", default=[0.5, 1.0, 2.0], help="Divergences")
    parser.add_argument("--replicates", type=int, default=10, help="Replicates per divergence (default: 10)")
    parser.add_argument("--alpha", type=float, default=1.0, help="Gamma alpha (default: 1.0)")
    parser.add_argument("--outdir", type=str, default="results", help="Output directory")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")

    args = parser.parse_args()

    df = run_benchmark(
        n_taxa=args.taxa,
        seq_length=args.length,
        divergences=args.divergences,
        replicates=args.replicates,
        alpha=args.alpha,
        outdir=args.outdir,
        seed=args.seed,
    )

    generate_plots(df, outdir=args.outdir)
    update_raw_results(df)
    print("\nBenchmark and reporting completed successfully!")
