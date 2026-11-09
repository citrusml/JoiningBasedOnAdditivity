"""
Large-Scale Benchmark Script for Phylogenetic Tree Reconstruction (N=32, 64, 128).

Compares MSA+ML (MAFFT + FastTree), Standard NJ, BioNJ, FastME-SPR,
AddNJ2 (Sister-Quartet Consistency), and Add-NNI across sequence sets of N=32, 64, 128.
Designed for execution on high-performance compute nodes or dedicated workstations.
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
from src.evaluation import evaluate_tree_rf
from src.add_nj import add_nj2
from src.fastme_wrapper import run_fastme
from src.ml_wrapper import run_msa_ml
from src.topology_search import TreeScorer, hill_climb_search

# Ensure stdout and stderr flush line-by-line even when redirected with nohup / file redirection
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(line_buffering=True)


def run_large_scale_benchmark(
    taxa_list: List[int] = [32, 64, 128],
    seq_length: int = 500,
    divergences: List[float] = [0.5, 1.0, 2.0],
    replicates: int = 10,
    alpha: float = 1.0,
    outdir: str = "results_large_scale",
    seed: int = 42,
) -> pd.DataFrame:
    """
    Executes large-scale benchmark on varying taxa counts.
    """
    os.makedirs(outdir, exist_ok=True)
    wag = WAGModel(wag_path="src/wag.dat", alpha=alpha)

    all_records: List[Dict[str, Any]] = []

    print("=" * 80)
    print("LARGE-SCALE PHYLOMETRICS BENCHMARK (N=32, 64, 128) WITH MSA + ML")
    print(f"Taxa: {taxa_list} | Length: {seq_length} aa | Divergences: {divergences} | Replicates: {replicates}")
    print("=" * 80)

    for n_taxa in taxa_list:
        print(f"\n=======================================================")
        print(f">>> BENCHMARK FOR TAXA SIZE: N = {n_taxa} <<<")
        print(f"=======================================================")

        for div in divergences:
            print(f"\n--- Divergence Scale: {div} substitutions/site (N={n_taxa}) ---")
            for rep in range(1, replicates + 1):
                rep_seed = seed + n_taxa * 10000 + int(div * 1000) + rep * 100

                # 1. True tree & patristic distance
                tree, taxa_names = generate_backward_yule_tree(
                    n_taxa=n_taxa, divergence=div, seed=rep_seed
                )
                true_newick = tree.to_newick() + ";"

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

                base_info = {
                    "replicate": rep,
                    "n_taxa": n_taxa,
                    "divergence": div,
                    "seq_length": seq_length,
                    "align_time_sec": float(align_time),
                }

                method_outputs: Dict[str, Tuple[str, float]] = {}

                # Method 1: MSA + ML (MAFFT + FastTree)
                t0 = time.time()
                ml_nwk = run_msa_ml(seqs)
                method_outputs["MSA_ML"] = (ml_nwk, time.time() - t0)

                # Method 2: Standard NJ
                t0 = time.time()
                nj_nwk = matrix_to_nj_newick(D_obs, taxa_names)
                method_outputs["Standard_NJ"] = (nj_nwk, time.time() - t0)

                # Method 3: BioNJ
                t0 = time.time()
                bionj_nwk = run_fastme(D_obs, taxa_names, method="bionj")
                method_outputs["BioNJ"] = (bionj_nwk, time.time() - t0)

                # Method 4: FastME-SPR
                t0 = time.time()
                fspr_nwk = run_fastme(D_obs, taxa_names, method="spr")
                method_outputs["FastME_SPR"] = (fspr_nwk, time.time() - t0)

                # Method 5: AddNJ2 (Proposed Method 1 Improved)
                t0 = time.time()
                addnj2_nwk = add_nj2(D_obs, taxa_names)
                method_outputs["AddNJ2"] = (addnj2_nwk, time.time() - t0)

                # Method 6: Add-NNI topology search (if N <= 64 for speed, or bounded)
                if n_taxa <= 64:
                    scorer = TreeScorer(D_obs, taxa_names)
                    t0 = time.time()
                    bionj_nni_nwk, _ = hill_climb_search(bionj_nwk, scorer, mode="nni", max_iter=20)
                    method_outputs["BioNJ_Add_NNI"] = (bionj_nni_nwk, time.time() - t0)

                    t0 = time.time()
                    addnj2_nni_nwk, _ = hill_climb_search(addnj2_nwk, scorer, mode="nni", max_iter=20)
                    method_outputs["AddNJ2_Add_NNI"] = (addnj2_nni_nwk, time.time() - t0)

                # Record results
                row = dict(base_info)
                for m_name, (inferred_str, runtime) in method_outputs.items():
                    eval_res = evaluate_tree_rf(inferred_str, true_newick, taxa_names)
                    row[f"{m_name}_nrf"] = float(eval_res["nrf"])
                    row[f"{m_name}_rf"] = float(eval_res["rf"])
                    row[f"{m_name}_time"] = float(runtime)

                all_records.append(row)

                out_str = f"  [N={n_taxa:3d}, D={div}] Rep {rep:2d}/{replicates}: "
                out_str += f"MSA_ML={row.get('MSA_ML_nrf', -1):.3f} | "
                out_str += f"NJ={row.get('Standard_NJ_nrf', -1):.3f} | "
                out_str += f"BioNJ={row.get('BioNJ_nrf', -1):.3f} | "
                out_str += f"AddNJ2={row.get('AddNJ2_nrf', -1):.3f}"
                if "BioNJ_Add_NNI_nrf" in row:
                    out_str += f" | BioNJ+NNI={row['BioNJ_Add_NNI_nrf']:.3f}"
                print(out_str)

    df = pd.DataFrame(all_records)
    csv_path = os.path.join(outdir, "large_scale_summary.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n[OK] Saved large-scale raw data to: {csv_path}")
    return df


def generate_large_scale_plots(df: pd.DataFrame, outdir: str = "results_large_scale") -> None:
    """
    Generates publication-quality scaling charts for varying Taxa N (300 DPI).
    """
    taxa_list = sorted(df["n_taxa"].unique())
    divergences = sorted(df["divergence"].unique())

    core_methods = ["MSA_ML", "Standard_NJ", "BioNJ", "FastME_SPR", "AddNJ2"]
    labels = ["MSA + ML (FastTree)", "Standard NJ", "BioNJ", "FastME-SPR", "AddNJ2 (Proposed)"]
    colors = ["#d62728", "#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]

    # 1. nRF vs Taxa count for each divergence scale
    fig, axes = plt.subplots(1, len(divergences), figsize=(16, 5), dpi=300, sharey=True)
    if len(divergences) == 1:
        axes = [axes]

    for ax, div in zip(axes, divergences):
        sub_df = df[df["divergence"] == div]
        for m, lbl, col in zip(core_methods, labels, colors):
            means = [sub_df[sub_df["n_taxa"] == n][f"{m}_nrf"].mean() for n in taxa_list]
            sems = [sub_df[sub_df["n_taxa"] == n][f"{m}_nrf"].sem() for n in taxa_list]
            ax.errorbar(taxa_list, means, yerr=sems, label=lbl, color=col, marker="o", capsize=3, lw=2)

        ax.set_title(f"Divergence Scale D = {div}", fontsize=13, fontweight="bold")
        ax.set_xlabel("Number of Taxa (N)", fontsize=11, fontweight="bold")
        ax.set_xticks(taxa_list)
        ax.set_ylim(0, max(0.3, df[[f"{m}_nrf" for m in core_methods if f"{m}_nrf" in df]].values.max() * 1.15))
        ax.grid(True, ls="--", alpha=0.5)

    axes[0].set_ylabel("Normalized RF Distance (nRF)", fontsize=12, fontweight="bold")
    axes[-1].legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=10)
    plt.tight_layout()

    fig_path = os.path.join(outdir, "nrf_vs_taxa.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"[OK] Saved figure: {fig_path}")

    # 2. Runtime vs Taxa count (log-log scale)
    plt.figure(figsize=(9, 5), dpi=300)
    for m, lbl, col in zip(core_methods, labels, colors):
        time_means = [df[df["n_taxa"] == n][f"{m}_time"].mean() for n in taxa_list]
        time_sems = [df[df["n_taxa"] == n][f"{m}_time"].sem() for n in taxa_list]
        plt.errorbar(taxa_list, time_means, yerr=time_sems, label=lbl, color=col, marker="s", capsize=3, lw=2)

    plt.xscale("log", base=2)
    plt.yscale("log")
    plt.xticks(taxa_list, [f"N={n}" for n in taxa_list], fontsize=11)
    plt.xlabel("Number of Taxa N (log scale)", fontsize=12, fontweight="bold")
    plt.ylabel("Execution Time per Replicate (seconds, log scale)", fontsize=12, fontweight="bold")
    plt.title("Computational Runtime Scaling (L=500 aa)", fontsize=13, fontweight="bold")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(fontsize=10)
    plt.tight_layout()

    fig_time = os.path.join(outdir, "runtime_vs_taxa.png")
    plt.savefig(fig_time, dpi=300)
    plt.close()
    print(f"[OK] Saved figure: {fig_time}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Large-Scale Phylometrics Benchmark (N=32, 64, 128).")
    parser.add_argument("--taxa", type=int, nargs="+", default=[32, 64, 128], help="List of taxa sizes")
    parser.add_argument("--length", type=int, default=500, help="Sequence length (default: 500)")
    parser.add_argument("--divergences", type=float, nargs="+", default=[0.5, 1.0, 2.0], help="Divergences")
    parser.add_argument("--replicates", type=int, default=10, help="Replicates per condition (default: 10)")
    parser.add_argument("--outdir", type=str, default="results_large_scale", help="Output directory")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")

    args = parser.parse_args()

    df = run_large_scale_benchmark(
        taxa_list=args.taxa,
        seq_length=args.length,
        divergences=args.divergences,
        replicates=args.replicates,
        outdir=args.outdir,
        seed=args.seed,
    )
    generate_large_scale_plots(df, outdir=args.outdir)
    print("\n[COMPLETE] Large-scale benchmark and visualization finished!")

