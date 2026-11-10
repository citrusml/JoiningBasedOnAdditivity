"""
HPC (LSF) Comprehensive Benchmark Script for Additivity-Based Phylogenetic Reconstruction.

Evaluates:
1. Taxa N in {20, 32, 64, 128} and Divergence D in {0.5, 1.0, 2.0} (default 30 replicates).
2. Three Distance Matrix Conditions:
   - 'blosum62': Pairwise Global Alignment (BLOSUM62) + Poisson distance
   - 'wag_gamma': Pairwise Global Alignment + WAG+Gamma Maximum Likelihood distance
   - 'true': Ground-truth patristic distance matrix from the true tree
3. All Reconstruction Methods across all N (including N=128):
   - MSA_ML (MAFFT + FastTree WAG+Gamma)
   - Standard_NJ, BioNJ, FastME_SPR, AddNJ2, BioNJ_Add_NNI, AddNJ2_Add_NNI
4. Sequence Similarity Score (SSS, w_ij) statistics per replicate.
5. Saves all inferred trees in Newick format (.jsonl and .nwk).
6. Computes 4-point condition violation (quartet additivity penalty) for every tree.
"""

from typing import List, Dict, Any, Tuple
import os
import sys
import json
import time
import argparse
import numpy as np
import pandas as pd

from src.tree_simulator import generate_backward_yule_tree, compute_patristic_distance_matrix
from src.sequence_simulator import (
    WAGModel,
    simulate_wag_gamma_sequences,
    compute_all_pairwise_matrices_and_sss,
)
from src.tree_inference import matrix_to_nj_newick
from src.evaluation import evaluate_tree_rf, evaluate_matrix_accuracy
from src.add_nj import add_nj2
from src.fastme_wrapper import run_fastme
from src.ml_wrapper import run_msa_ml
from src.topology_search import TreeScorer, hill_climb_search

# Ensure line-buffered stdout/stderr for real-time LSF log monitoring
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(line_buffering=True)


DISTANCE_METHODS = [
    "Standard_NJ",
    "BioNJ",
    "FastME_NNI",
    "FastME_SPR",
    "AddNJ2",
    "BioNJ_Add_NNI",
    "AddNJ2_Add_NNI",
]

MATRIX_TYPES = ["blosum62", "wag_gamma", "sss_linear", "sss_log", "true"]


def run_distance_methods_on_matrix(
    D_mat: np.ndarray,
    taxa_names: List[str],
    scorer_search: TreeScorer,
    max_nni_iter: int = 20,
) -> Dict[str, Tuple[str, float]]:
    """
    Runs all 7 distance-based tree reconstruction methods on a given distance matrix.

    Parameters
    ----------
    D_mat : np.ndarray
        Input symmetric distance matrix of shape (n, n) in substitutions/site.
    taxa_names : List[str]
        Ordered list of taxon names (unitless).
    scorer_search : TreeScorer
        Pre-initialized TreeScorer for D_mat used by Add-NNI search.
    max_nni_iter : int, default=20
        Maximum hill-climbing iterations for NNI topology search (unitless).

    Returns
    -------
    Dict[str, Tuple[str, float]]
        Mapping from method name to (inferred_newick_str, runtime_seconds).
    """
    outputs: Dict[str, Tuple[str, float]] = {}

    # 1. Standard NJ
    t0 = time.time()
    nj_nwk = matrix_to_nj_newick(D_mat, taxa_names)
    outputs["Standard_NJ"] = (nj_nwk, time.time() - t0)

    # 2. BioNJ
    t0 = time.time()
    bionj_nwk = run_fastme(D_mat, taxa_names, method="bionj")
    t_bionj = time.time() - t0
    outputs["BioNJ"] = (bionj_nwk, t_bionj)

    # 3. FastME-NNI
    t0 = time.time()
    fnni_nwk = run_fastme(D_mat, taxa_names, method="nni")
    outputs["FastME_NNI"] = (fnni_nwk, time.time() - t0)

    # 4. FastME-SPR
    t0 = time.time()
    fspr_nwk = run_fastme(D_mat, taxa_names, method="spr")
    outputs["FastME_SPR"] = (fspr_nwk, time.time() - t0)

    # 5. AddNJ2 (Sister-Quartet Consistency)
    t0 = time.time()
    addnj2_nwk = add_nj2(D_mat, taxa_names)
    t_addnj2 = time.time() - t0
    outputs["AddNJ2"] = (addnj2_nwk, t_addnj2)

    # 6. BioNJ + Add-NNI
    t0 = time.time()
    bionj_nni_nwk, _ = hill_climb_search(bionj_nwk, scorer_search, mode="nni", max_iter=max_nni_iter)
    outputs["BioNJ_Add_NNI"] = (bionj_nni_nwk, t_bionj + (time.time() - t0))

    # 7. AddNJ2 + Add-NNI
    t0 = time.time()
    addnj2_nni_nwk, _ = hill_climb_search(addnj2_nwk, scorer_search, mode="nni", max_iter=max_nni_iter)
    outputs["AddNJ2_Add_NNI"] = (addnj2_nni_nwk, t_addnj2 + (time.time() - t0))

    return outputs


def run_hpc_benchmark_condition(
    n_taxa: int,
    divergence: float,
    seq_length: int = 500,
    rep_start: int = 1,
    rep_end: int = 30,
    alpha: float = 1.0,
    indel_rate: float = 0.0,
    outdir: str = "results_hpc",
    seed: int = 42,
    nni_max_quartets: int = 25000,
    eval_max_quartets: int = 100000,
) -> pd.DataFrame:
    """
    Executes the full HPC benchmark for a specific (n_taxa, divergence, indel_rate) condition
    across replicates [rep_start, rep_end].

    Parameters
    ----------
    n_taxa : int
        Number of taxa N (unitless).
    divergence : float
        Branch length divergence scale D (substitutions/site).
    seq_length : int, default=500
        Sequence length L (amino acids).
    rep_start : int, default=1
        Starting replicate index (1-based, inclusive).
    rep_end : int, default=30
        Ending replicate index (1-based, inclusive).
    alpha : float, default=1.0
        Gamma rate heterogeneity shape parameter (unitless).
    indel_rate : float, default=0.0
        Zipfian (POW 1.7 50) insertion/deletion rate via INDELible (unitless, e.g. 0.0, 0.05, 0.1).
    outdir : str, default='results_hpc'
        Output root directory path.
    seed : int, default=42
        Base random seed (unitless).
    nni_max_quartets : int, default=25000
        Maximum quartets evaluated per step during NNI search (for speed at N=128).
    eval_max_quartets : int, default=100000
        Maximum quartets evaluated for final 4-point condition violation scoring.

    Returns
    -------
    pd.DataFrame
        DataFrame of benchmark records for this condition.
    """
    shards_dir = os.path.join(outdir, "shards")
    trees_dir = os.path.join(outdir, "trees")
    os.makedirs(shards_dir, exist_ok=True)
    os.makedirs(trees_dir, exist_ok=True)

    indel_suffix = f"_I{indel_rate}" if indel_rate > 0.0 else ""
    tag = f"N{n_taxa}_D{divergence}{indel_suffix}_R{rep_start}_{rep_end}"
    csv_path = os.path.join(shards_dir, f"summary_{tag}.csv")
    jsonl_tree_path = os.path.join(trees_dir, f"trees_{tag}.jsonl")
    if os.path.exists(jsonl_tree_path):
        os.remove(jsonl_tree_path)

    wag = WAGModel(wag_path="src/wag.dat", alpha=alpha)
    records: List[Dict[str, Any]] = []

    print("=" * 88)
    print(
        f"HPC PHYLOMETRICS BENCHMARK | N={n_taxa} | D={divergence} | Indel={indel_rate} | Reps={rep_start}..{rep_end}"
    )
    print("=" * 88)

    for rep in range(rep_start, rep_end + 1):
        t_rep_start = time.time()
        rep_seed = seed + n_taxa * 10000 + int(divergence * 1000) + rep * 100

        # 1. Generate True Tree & True Patristic Distance Matrix
        tree, taxa_names = generate_backward_yule_tree(
            n_taxa=n_taxa, divergence=divergence, seed=rep_seed
        )
        true_newick = tree.to_newick() + ";"
        D_true = compute_patristic_distance_matrix(tree, taxa_names)

        # 2. Simulate WAG+Gamma Sequences (with Zipfian Indels via INDELible when indel_rate > 0)
        seqs = simulate_wag_gamma_sequences(
            root=tree,
            seq_length=seq_length,
            alpha=alpha,
            model=wag,
            seed=rep_seed + 1,
            indel_rate=indel_rate,
        )

        # 3. Compute BLOSUM62+Poisson matrix, WAG+Gamma ML matrix, and Exp 18 SSS distance matrices
        t_mat0 = time.time()
        (
            D_blosum62,
            D_wag_gamma,
            sss_stats,
            D_sss_linear,
            D_sss_log,
        ) = compute_all_pairwise_matrices_and_sss(
            seqs, taxa_names, model=wag, alpha=alpha, mode="global", return_sss_matrices=True
        )
        dist_calc_time = time.time() - t_mat0

        matrices: Dict[str, np.ndarray] = {
            "blosum62": D_blosum62,
            "wag_gamma": D_wag_gamma,
            "sss_linear": D_sss_linear,
            "sss_log": D_sss_log,
            "true": D_true,
        }

        # 4. Initialize evaluation scorers (shared quartet index set across all matrices)
        eval_scorers: Dict[str, TreeScorer] = {
            m_type: TreeScorer(D_mat, taxa_names, max_quartets=eval_max_quartets)
            for m_type, D_mat in matrices.items()
        }
        search_scorers: Dict[str, TreeScorer] = {
            m_type: (
                eval_scorers[m_type]
                if n_taxa <= 32
                else TreeScorer(D_mat, taxa_names, max_quartets=nni_max_quartets)
            )
            for m_type, D_mat in matrices.items()
        }

        row: Dict[str, Any] = {
            "replicate": rep,
            "n_taxa": n_taxa,
            "divergence": divergence,
            "indel_rate": float(indel_rate),
            "seq_length": seq_length,
            "dist_calc_time_sec": float(dist_calc_time),
        }
        row.update(sss_stats)

        # Distance matrix accuracy & intrinsic 4-point condition defect Delta_ADD(D)
        for m_type in [m for m in MATRIX_TYPES if m != "true"]:
            acc = evaluate_matrix_accuracy(matrices[m_type], D_true)
            row[f"rmse_{m_type}"] = acc["rmse"]
            row[f"mae_{m_type}"] = acc["mae"]
            row[f"mre_{m_type}"] = acc["mre"]

        for m_type in MATRIX_TYPES:
            row[f"delta_add_{m_type}"] = eval_scorers[m_type].compute_matrix_defect()

        # Score true tree's 4-point condition violation against each distance matrix
        true_choices = eval_scorers["true"].extract_tree_quartet_choices(true_newick)
        for m_type in MATRIX_TYPES:
            row[f"true_tree_penalty_{m_type}"] = eval_scorers[m_type].score_choices(true_choices)

        saved_trees_rep: Dict[str, str] = {"true_tree": true_newick}

        # 5. Run MSA + ML (MAFFT + FastTree WAG+Gamma)
        t0 = time.time()
        ml_nwk = run_msa_ml(seqs)
        ml_time = time.time() - t0
        saved_trees_rep["MSA_ML"] = ml_nwk

        ml_rf_res = evaluate_tree_rf(ml_nwk, true_newick, taxa_names, D_true=D_true)
        row["MSA_ML_nrf"] = float(ml_rf_res["nrf"])
        row["MSA_ML_rf"] = float(ml_rf_res["rf"])
        row["MSA_ML_kf_bsd"] = float(ml_rf_res["kf_bsd"])
        row["MSA_ML_patristic_rmse"] = float(ml_rf_res["patristic_rmse"])
        row["MSA_ML_time"] = float(ml_time)

        ml_choices = eval_scorers["true"].extract_tree_quartet_choices(ml_nwk)
        for m_type in MATRIX_TYPES:
            row[f"MSA_ML_penalty_{m_type}"] = eval_scorers[m_type].score_choices(ml_choices)

        # 6. Run all distance-based methods on each of the 3 distance matrices
        for m_type in MATRIX_TYPES:
            method_outs = run_distance_methods_on_matrix(
                matrices[m_type],
                taxa_names,
                scorer_search=search_scorers[m_type],
                max_nni_iter=15 if n_taxa >= 128 else 20,
            )

            for m_name, (inferred_nwk, runtime) in method_outs.items():
                col_prefix = f"{m_type}_{m_name}"
                saved_trees_rep[col_prefix] = inferred_nwk

                eval_res = evaluate_tree_rf(inferred_nwk, true_newick, taxa_names, D_true=D_true)
                row[f"{col_prefix}_nrf"] = float(eval_res["nrf"])
                row[f"{col_prefix}_rf"] = float(eval_res["rf"])
                row[f"{col_prefix}_kf_bsd"] = float(eval_res["kf_bsd"])
                row[f"{col_prefix}_patristic_rmse"] = float(eval_res["patristic_rmse"])
                row[f"{col_prefix}_time"] = float(runtime)

                # Extract quartet splits induced by inferred tree once and score wrt D_obs and D_true
                t_choices = eval_scorers[m_type].extract_tree_quartet_choices(inferred_nwk)
                row[f"{col_prefix}_penalty_obs"] = eval_scorers[m_type].score_choices(t_choices)
                row[f"{col_prefix}_penalty_true"] = eval_scorers["true"].score_choices(t_choices)

        records.append(row)

        # 7. Save Newick trees (both consolidated JSONL and individual .nwk files)
        with open(jsonl_tree_path, "a") as jf:
            jf.write(
                json.dumps(
                    {
                        "n_taxa": n_taxa,
                        "divergence": divergence,
                        "indel_rate": float(indel_rate),
                        "replicate": rep,
                        "trees": saved_trees_rep,
                    }
                )
                + "\n"
            )

        rep_tree_dir = os.path.join(trees_dir, f"N{n_taxa}_D{divergence}{indel_suffix}", f"rep_{rep:03d}")
        os.makedirs(rep_tree_dir, exist_ok=True)
        for t_label, nwk_str in saved_trees_rep.items():
            with open(os.path.join(rep_tree_dir, f"{t_label}.nwk"), "w") as tf:
                tf.write(nwk_str.strip() + "\n")

        # Save incremental CSV checkpoint after each replicate
        pd.DataFrame(records).to_csv(csv_path, index=False)

        elapsed = time.time() - t_rep_start
        print(
            f"  [N={n_taxa:3d}, D={divergence}, Indel={indel_rate}, Rep {rep:2d}/{rep_end}] "
            f"SSS={row['sss_mean']:.3f} | "
            f"ML={row['MSA_ML_nrf']:.3f} | "
            f"BLOSUM(NJ={row['blosum62_Standard_NJ_nrf']:.3f}, AddNJ2={row['blosum62_AddNJ2_nrf']:.3f}) | "
            f"WAG(NJ={row['wag_gamma_Standard_NJ_nrf']:.3f}, AddNJ2={row['wag_gamma_AddNJ2_nrf']:.3f}) | "
            f"SSS_lin(NJ={row['sss_linear_Standard_NJ_nrf']:.3f}, AddNJ2={row['sss_linear_AddNJ2_nrf']:.3f}) | "
            f"SSS_log(NJ={row['sss_log_Standard_NJ_nrf']:.3f}, AddNJ2={row['sss_log_AddNJ2_nrf']:.3f}) | "
            f"TRUE(AddNJ2={row['true_AddNJ2_nrf']:.3f}) "
            f"({elapsed:.1f}s)"
        )

    print(f"\n[OK] Saved shard summary CSV: {csv_path}")
    print(f"[OK] Saved Newick trees JSONL: {jsonl_tree_path}")
    return pd.DataFrame(records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run HPC Phylogenetic Benchmark across BLOSUM62, WAG+Gamma, and True distance matrices."
    )
    parser.add_argument("--taxa", type=int, nargs="+", default=[20, 32, 64, 128], help="Taxa sizes N")
    parser.add_argument("--divergence", type=float, nargs="+", default=[0.5, 1.0, 2.0], help="Divergence levels D")
    parser.add_argument(
        "--indel-rate",
        type=float,
        nargs="+",
        default=[0.0],
        help="Zipfian (POW 1.7 50) insertion/deletion rates via INDELible (e.g. 0.05 0.1)",
    )
    parser.add_argument("--length", type=int, default=500, help="Amino acid sequence length L")
    parser.add_argument("--replicates", type=int, default=30, help="Total replicates per condition")
    parser.add_argument("--rep-start", type=int, default=1, help="Start replicate index (1-based)")
    parser.add_argument("--rep-end", type=int, default=None, help="End replicate index (inclusive)")
    parser.add_argument("--outdir", type=str, default="results_hpc", help="Output directory")
    parser.add_argument("--seed", type=int, default=42, help="Base random seed")
    args = parser.parse_args()

    rep_end = args.rep_end if args.rep_end is not None else args.replicates

    for indel_val in args.indel_rate:
        for n_val in args.taxa:
            for d_val in args.divergence:
                run_hpc_benchmark_condition(
                    n_taxa=n_val,
                    divergence=d_val,
                    seq_length=args.length,
                    rep_start=args.rep_start,
                    rep_end=rep_end,
                    indel_rate=indel_val,
                    outdir=args.outdir,
                    seed=args.seed,
                )

