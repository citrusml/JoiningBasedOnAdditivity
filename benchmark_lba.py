"""
Comprehensive Benchmark Script for Long Branch Attraction (LBA) Effects
across Multiple Distance Matrices (BLOSUM62+Poisson, WAG+Gamma ML, True Patristic).

Evaluates:
1. Taxa N in {20, 32}, Divergence D in {0.5, 1.0, 2.0}, LBA Ratio in {1.0, 3.0}, Replicates = 10.
2. Sequence Similarity Score (SSS, w_ij) via MMseqs2 and Smith-Waterman local alignment.
3. Three Distance Matrix Types:
   - 'blosum62': Pairwise Global Alignment (BLOSUM62) + Poisson distance
   - 'wag_gamma': Pairwise Global Alignment + WAG+Gamma Maximum Likelihood distance
   - 'true': Ground-truth patristic distance matrix from the true tree
4. Matrix-level and Tree-level 4-Point Condition Violations (Additivity Defect & Penalty).
5. All 15 Reconstruction Methods (including FastME-NNI and MSA+ML).
6. Saves all generated phylogenetic trees (.nwk, .jsonl, and .csv) and simulated sequences (.fasta).
"""

from typing import List, Dict, Any, Tuple, Optional
import os
import sys
import io
import json
import time
import argparse
import concurrent.futures
import numpy as np
import pandas as pd

from src.tree_simulator import generate_backward_yule_tree, compute_patristic_distance_matrix
from src.sequence_simulator import (
    WAGModel,
    simulate_wag_gamma_sequences,
    compute_all_pairwise_matrices_and_sss,
)
from src.tree_inference import matrix_to_nj_newick
from src.evaluation import evaluate_tree_lba, evaluate_matrix_accuracy, evaluate_four_point_metrics
from src.four_point_core import FourPointCore
from src.add_nj import add_nj, add_nj2
from src.fastme_wrapper import run_fastme
from src.ml_wrapper import run_msa_ml
from src.topology_search import TreeScorer, hill_climb_search


MATRIX_TYPES = ["blosum62", "wag_gamma", "true"]

DISTANCE_METHODS = [
    "Standard_NJ",
    "BioNJ",
    "FastME_NNI",
    "FastME_SPR",
    "AddNJ",
    "AddNJ2",
    "NJ_Add_NNI",
    "NJ_Add_SPR",
    "NJ_Add_Hybrid",
    "BioNJ_Add_NNI",
    "BioNJ_Add_SPR",
    "BioNJ_Add_Hybrid",
    "AddNJ2_Add_NNI",
    "AddNJ2_Add_Hybrid",
]


def run_distance_methods_on_matrix(
    D_mat: np.ndarray,
    taxa_names: List[str],
    scorer: TreeScorer,
    max_iter_search: int = 20,
) -> Dict[str, Tuple[str, float]]:
    """
    Runs all 14 distance-based phylogenetic reconstruction methods on a single distance matrix.

    Parameters
    ----------
    D_mat : np.ndarray
        Symmetric distance matrix of shape (n_taxa, n_taxa) [substitutions/site].
    taxa_names : List[str]
        Ordered list of taxon names [unitless].
    scorer : TreeScorer
        Configured TreeScorer for D_mat [unitless].
    max_iter_search : int, default=20
        Maximum iterations per hill-climbing search phase [unitless].

    Returns
    -------
    Dict[str, Tuple[str, float]]
        Mapping from method_name to (inferred_newick, runtime_sec).
    """
    outputs: Dict[str, Tuple[str, float]] = {}

    # 1. Standard NJ
    t0 = time.time()
    nj_nwk = matrix_to_nj_newick(D_mat, taxa_names)
    t_nj = time.time() - t0
    outputs["Standard_NJ"] = (nj_nwk, t_nj)

    # 2. BioNJ
    t0 = time.time()
    bionj_nwk = run_fastme(D_mat, taxa_names, method="bionj")
    t_bionj = time.time() - t0
    outputs["BioNJ"] = (bionj_nwk, t_bionj)

    # 3. FastME-NNI
    t0 = time.time()
    fnni_nwk = run_fastme(D_mat, taxa_names, method="nni")
    t_fnni = time.time() - t0
    outputs["FastME_NNI"] = (fnni_nwk, t_fnni)

    # 4. FastME-SPR
    t0 = time.time()
    fspr_nwk = run_fastme(D_mat, taxa_names, method="spr")
    t_fspr = time.time() - t0
    outputs["FastME_SPR"] = (fspr_nwk, t_fspr)

    # 5. AddNJ (Method 1 original)
    t0 = time.time()
    addnj_nwk = add_nj(D_mat, taxa_names)
    t_addnj = time.time() - t0
    outputs["AddNJ"] = (addnj_nwk, t_addnj)

    # 6. AddNJ2 (Method 1 sister-pair consistency)
    t0 = time.time()
    addnj2_nwk = add_nj2(D_mat, taxa_names)
    t_addnj2 = time.time() - t0
    outputs["AddNJ2"] = (addnj2_nwk, t_addnj2)

    # 7. NJ + Add-NNI
    t0 = time.time()
    nj_nni_nwk, _ = hill_climb_search(nj_nwk, scorer, mode="nni", max_iter=max_iter_search)
    t_nj_nni = time.time() - t0
    outputs["NJ_Add_NNI"] = (nj_nni_nwk, t_nj_nni)

    # 8. NJ + Add-SPR
    t0 = time.time()
    nj_spr_nwk, _ = hill_climb_search(nj_nwk, scorer, mode="spr", max_iter=max_iter_search)
    t_nj_spr = time.time() - t0
    outputs["NJ_Add_SPR"] = (nj_spr_nwk, t_nj_spr)

    # 9. NJ + Add-Hybrid (NNI -> SPR)
    t0 = time.time()
    nj_hyb_nwk, _ = hill_climb_search(nj_nni_nwk, scorer, mode="spr", max_iter=max_iter_search)
    t_nj_hyb = t_nj_nni + (time.time() - t0)
    outputs["NJ_Add_Hybrid"] = (nj_hyb_nwk, t_nj_hyb)

    # 10. BioNJ + Add-NNI
    t0 = time.time()
    bionj_nni_nwk, _ = hill_climb_search(bionj_nwk, scorer, mode="nni", max_iter=max_iter_search)
    t_bionj_nni = time.time() - t0
    outputs["BioNJ_Add_NNI"] = (bionj_nni_nwk, t_bionj_nni)

    # 11. BioNJ + Add-SPR
    t0 = time.time()
    bionj_spr_nwk, _ = hill_climb_search(bionj_nwk, scorer, mode="spr", max_iter=max_iter_search)
    t_bionj_spr = time.time() - t0
    outputs["BioNJ_Add_SPR"] = (bionj_spr_nwk, t_bionj_spr)

    # 12. BioNJ + Add-Hybrid (NNI -> SPR)
    t0 = time.time()
    bionj_hyb_nwk, _ = hill_climb_search(bionj_nni_nwk, scorer, mode="spr", max_iter=max_iter_search)
    t_bionj_hyb = t_bionj_nni + (time.time() - t0)
    outputs["BioNJ_Add_Hybrid"] = (bionj_hyb_nwk, t_bionj_hyb)

    # 13. AddNJ2 + Add-NNI
    t0 = time.time()
    addnj2_nni_nwk, _ = hill_climb_search(addnj2_nwk, scorer, mode="nni", max_iter=max_iter_search)
    t_addnj2_nni = time.time() - t0
    outputs["AddNJ2_Add_NNI"] = (addnj2_nni_nwk, t_addnj2_nni)

    # 14. AddNJ2 + Add-Hybrid (NNI -> SPR)
    t0 = time.time()
    addnj2_hyb_nwk, _ = hill_climb_search(addnj2_nni_nwk, scorer, mode="spr", max_iter=max_iter_search)
    t_addnj2_hyb = t_addnj2_nni + (time.time() - t0)
    outputs["AddNJ2_Add_Hybrid"] = (addnj2_hyb_nwk, t_addnj2_hyb)

    return outputs


def run_single_lba_trial(
    n_taxa: int,
    divergence: float,
    lba_ratio: float,
    rep: int,
    seq_length: int = 500,
    alpha: float = 1.0,
    seed: int = 42,
    max_iter_search: int = 20,
    outdir: str = "results_lba",
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """
    Executes a single simulation trial for a specific (n_taxa, divergence, lba_ratio, rep)
    across all 3 distance matrices and MSA+ML, saving all generated trees and sequences.

    Returns
    -------
    Tuple[Dict[str, Any], List[Dict[str, Any]]]
        - wide_record : Single dictionary with all metrics for wide summary CSV.
        - tree_rows : List of long-format dictionaries (one per inferred tree) with Newick strings.
    """
    rep_seed = seed + n_taxa * 100000 + int(divergence * 10000) + int(lba_ratio * 1000) + rep * 100

    wag = WAGModel(wag_path="src/wag.dat", alpha=alpha)
    fp_core = FourPointCore(n_taxa=n_taxa)

    # 1. Generate true tree with designated lba_ratio
    tree, taxa_names = generate_backward_yule_tree(
        n_taxa=n_taxa,
        divergence=divergence,
        seed=rep_seed,
        lba_ratio=lba_ratio,
    )
    true_newick = tree.to_newick() + ";"
    D_true = compute_patristic_distance_matrix(tree, taxa_names)
    lba_pair = (taxa_names[0], taxa_names[-1])

    # 2. Simulate WAG+Gamma protein sequences
    seqs = simulate_wag_gamma_sequences(
        root=tree,
        seq_length=seq_length,
        alpha=alpha,
        model=wag,
        seed=rep_seed + 1,
    )

    # 3. Compute BLOSUM62+Poisson matrix, WAG+Gamma ML matrix, and SSS statistics
    t_align_start = time.time()
    D_blosum62, D_wag_gamma, sss_stats = compute_all_pairwise_matrices_and_sss(
        seqs, taxa_names, model=wag, alpha=alpha, mode="global"
    )
    align_time = time.time() - t_align_start

    matrices: Dict[str, np.ndarray] = {
        "blosum62": D_blosum62,
        "wag_gamma": D_wag_gamma,
        "true": D_true,
    }

    scorers: Dict[str, TreeScorer] = {
        m_type: TreeScorer(D_mat, taxa_names)
        for m_type, D_mat in matrices.items()
    }

    true_lba = evaluate_tree_lba(true_newick, true_newick, taxa_names, lba_pair=lba_pair)
    ds_true_val = float(true_lba["ds_true"])

    record: Dict[str, Any] = {
        "replicate": rep,
        "n_taxa": n_taxa,
        "divergence": divergence,
        "lba_ratio": lba_ratio,
        "seq_length": seq_length,
        "align_time_sec": float(align_time),
        "ds_true": ds_true_val,
    }
    record.update(sss_stats)

    # 4. Evaluate distance matrix accuracy and 4-point condition violation Delta_ADD(D)
    for m_type in MATRIX_TYPES:
        fp_metrics = evaluate_four_point_metrics(matrices[m_type], fp_core)
        record[f"delta_add_{m_type}"] = float(fp_metrics["delta_add"])
        record[f"add_abs_{m_type}"] = float(fp_metrics["add_abs"])
        record[f"H_tri_{m_type}"] = float(fp_metrics["H_tri"])
        record[f"n_tri_viol_{m_type}"] = float(fp_metrics["n_tri_viol"])
        if m_type != "true":
            acc = evaluate_matrix_accuracy(matrices[m_type], D_true)
            record[f"rmse_{m_type}"] = float(acc["rmse"])
            record[f"mae_{m_type}"] = float(acc["mae"])
            record[f"mre_{m_type}"] = float(acc["mre"])

    # Score true tree's 4-point condition penalty on each matrix
    true_choices = scorers["true"].extract_tree_quartet_choices(true_newick)
    for m_type in MATRIX_TYPES:
        record[f"true_tree_penalty_{m_type}"] = float(scorers[m_type].score_choices(true_choices))

    # Prepare directory for saving trees and sequences of this replicate
    cond_dir = os.path.join(
        outdir,
        "trees",
        f"N{n_taxa}_D{divergence}_LBA{lba_ratio}",
        f"rep_{rep:03d}",
    )
    os.makedirs(cond_dir, exist_ok=True)

    # Save true tree and simulated FASTA sequences
    with open(os.path.join(cond_dir, "true_tree.nwk"), "w") as f:
        f.write(true_newick.strip() + "\n")
    with open(os.path.join(cond_dir, "sequences.fasta"), "w") as f:
        for t_name in taxa_names:
            f.write(f">{t_name}\n{seqs[t_name]}\n")

    tree_rows: List[Dict[str, Any]] = []
    tree_rows.append(
        {
            "n_taxa": n_taxa,
            "divergence": divergence,
            "lba_ratio": lba_ratio,
            "replicate": rep,
            "sss_mean": record["sss_mean"],
            "sss_lba_pair": record.get("sss_lba_pair", 0.0),
            "matrix_type": "ground_truth",
            "matrix_delta_add": 0.0,
            "method": "True_Tree",
            "nrf": 0.0,
            "rf": 0.0,
            "ds_est": ds_true_val,
            "ds_true": ds_true_val,
            "delta_ds": 0.0,
            "is_sister": 1.0 if ds_true_val == 2.0 else 0.0,
            "penalty_input": 0.0,
            "penalty_blosum62": record["true_tree_penalty_blosum62"],
            "penalty_wag_gamma": record["true_tree_penalty_wag_gamma"],
            "penalty_true": 0.0,
            "time_sec": 0.0,
            "newick": true_newick.strip(),
        }
    )

    # 5. Run MSA + ML (MAFFT + FastTree WAG+Gamma)
    t0 = time.time()
    try:
        ml_nwk = run_msa_ml(seqs)
        t_ml = time.time() - t0
    except Exception:
        ml_nwk = ""
        t_ml = 0.0

    if ml_nwk:
        with open(os.path.join(cond_dir, "MSA_ML.nwk"), "w") as f:
            f.write(ml_nwk.strip() + "\n")
        ml_eval = evaluate_tree_lba(ml_nwk, true_newick, taxa_names, lba_pair=lba_pair)
        ml_choices = scorers["true"].extract_tree_quartet_choices(ml_nwk)
        p_ml_blosum = float(scorers["blosum62"].score_choices(ml_choices))
        p_ml_wag = float(scorers["wag_gamma"].score_choices(ml_choices))
        p_ml_true = float(scorers["true"].score_choices(ml_choices))

        record["MSA_ML_nrf"] = float(ml_eval["nrf"])
        record["MSA_ML_rf"] = float(ml_eval["rf"])
        record["MSA_ML_kf_bsd"] = float(ml_eval["kf_bsd"])
        record["MSA_ML_patristic_rmse"] = float(ml_eval["patristic_rmse"])
        record["MSA_ML_ds_est"] = float(ml_eval["ds_est"])
        record["MSA_ML_delta_ds"] = float(ml_eval["delta_ds"])
        record["MSA_ML_is_sister"] = float(ml_eval["is_sister"])
        record["MSA_ML_penalty"] = p_ml_blosum
        record["MSA_ML_penalty_blosum62"] = p_ml_blosum
        record["MSA_ML_penalty_wag_gamma"] = p_ml_wag
        record["MSA_ML_penalty_true"] = p_ml_true
        record["MSA_ML_time"] = float(t_ml)

        tree_rows.append(
            {
                "n_taxa": n_taxa,
                "divergence": divergence,
                "lba_ratio": lba_ratio,
                "replicate": rep,
                "sss_mean": record["sss_mean"],
                "sss_lba_pair": record.get("sss_lba_pair", 0.0),
                "matrix_type": "msa_ml",
                "matrix_delta_add": record["delta_add_wag_gamma"],
                "method": "MSA_ML",
                "nrf": float(ml_eval["nrf"]),
                "rf": float(ml_eval["rf"]),
                "ds_est": float(ml_eval["ds_est"]),
                "ds_true": ds_true_val,
                "delta_ds": float(ml_eval["delta_ds"]),
                "is_sister": float(ml_eval["is_sister"]),
                "penalty_input": p_ml_wag,
                "penalty_blosum62": p_ml_blosum,
                "penalty_wag_gamma": p_ml_wag,
                "penalty_true": p_ml_true,
                "time_sec": float(t_ml),
                "newick": ml_nwk.strip(),
            }
        )

    # 6. Run all 14 distance-based methods on each of the 3 distance matrices
    for m_type in MATRIX_TYPES:
        m_dir = os.path.join(cond_dir, m_type)
        os.makedirs(m_dir, exist_ok=True)

        method_outs = run_distance_methods_on_matrix(
            matrices[m_type],
            taxa_names,
            scorer=scorers[m_type],
            max_iter_search=max_iter_search,
        )

        for m_name, (inferred_nwk, runtime) in method_outs.items():
            # Save individual .nwk file
            with open(os.path.join(m_dir, f"{m_name}.nwk"), "w") as f:
                f.write(inferred_nwk.strip() + "\n")

            eval_res = evaluate_tree_lba(inferred_nwk, true_newick, taxa_names, lba_pair=lba_pair)
            t_choices = scorers[m_type].extract_tree_quartet_choices(inferred_nwk)

            p_input = float(scorers[m_type].score_choices(t_choices))
            p_blosum = float(scorers["blosum62"].score_choices(t_choices))
            p_wag = float(scorers["wag_gamma"].score_choices(t_choices))
            p_true = float(scorers["true"].score_choices(t_choices))

            col_prefix = f"{m_type}_{m_name}"
            record[f"{col_prefix}_nrf"] = float(eval_res["nrf"])
            record[f"{col_prefix}_rf"] = float(eval_res["rf"])
            record[f"{col_prefix}_kf_bsd"] = float(eval_res["kf_bsd"])
            record[f"{col_prefix}_patristic_rmse"] = float(eval_res["patristic_rmse"])
            record[f"{col_prefix}_ds_est"] = float(eval_res["ds_est"])
            record[f"{col_prefix}_delta_ds"] = float(eval_res["delta_ds"])
            record[f"{col_prefix}_is_sister"] = float(eval_res["is_sister"])
            record[f"{col_prefix}_penalty_obs"] = p_input
            record[f"{col_prefix}_penalty_true"] = p_true
            record[f"{col_prefix}_time"] = float(runtime)

            # Alias default <m_name>_* to blosum62_<m_name>_* for backwards compatibility
            if m_type == "blosum62":
                record[f"{m_name}_nrf"] = float(eval_res["nrf"])
                record[f"{m_name}_rf"] = float(eval_res["rf"])
                record[f"{m_name}_ds_est"] = float(eval_res["ds_est"])
                record[f"{m_name}_delta_ds"] = float(eval_res["delta_ds"])
                record[f"{m_name}_is_sister"] = float(eval_res["is_sister"])
                record[f"{m_name}_penalty"] = p_input
                record[f"{m_name}_time"] = float(runtime)

            tree_rows.append(
                {
                    "n_taxa": n_taxa,
                    "divergence": divergence,
                    "lba_ratio": lba_ratio,
                    "replicate": rep,
                    "sss_mean": record["sss_mean"],
                    "sss_lba_pair": record.get("sss_lba_pair", 0.0),
                    "matrix_type": m_type,
                    "matrix_delta_add": record[f"delta_add_{m_type}"],
                    "method": m_name,
                    "nrf": float(eval_res["nrf"]),
                    "rf": float(eval_res["rf"]),
                    "ds_est": float(eval_res["ds_est"]),
                    "ds_true": ds_true_val,
                    "delta_ds": float(eval_res["delta_ds"]),
                    "is_sister": float(eval_res["is_sister"]),
                    "penalty_input": p_input,
                    "penalty_blosum62": p_blosum,
                    "penalty_wag_gamma": p_wag,
                    "penalty_true": p_true,
                    "time_sec": float(runtime),
                    "newick": inferred_nwk.strip(),
                }
            )

    return record, tree_rows


def run_lba_benchmark(
    taxa_list: List[int] = [20, 32],
    divergences: List[float] = [0.5, 1.0, 2.0],
    lba_ratios: List[float] = [1.0, 3.0],
    replicates: int = 10,
    seq_length: int = 500,
    alpha: float = 1.0,
    outdir: str = "results_lba",
    seed: int = 42,
    threads: int = 6,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Executes the comprehensive LBA benchmark suite across the parameter grid.
    """
    os.makedirs(outdir, exist_ok=True)
    csv_path = os.path.join(outdir, "lba_summary.csv")
    long_trees_csv_path = os.path.join(outdir, "all_inferred_trees.csv")

    tasks = []
    for n_taxa in taxa_list:
        for div in divergences:
            for lba_r in lba_ratios:
                for rep in range(1, replicates + 1):
                    tasks.append((n_taxa, div, lba_r, rep))

    total_tasks = len(tasks)
    print("=" * 95)
    print("LBA BENCHMARK SUITE: 3 DISTANCE MATRICES (BLOSUM62, WAG+GAMMA, TRUE) + SSS + 4PC VIOLATION")
    print(f"Taxa: {taxa_list} | Divergences: {divergences} | LBA Ratios: {lba_ratios} | Replicates: {replicates}")
    print(f"Total Datasets: {total_tasks} (x 3 matrices x 14 dist methods + MSA_ML = {total_tasks * 43} trees)")
    print(f"Worker Processes: {threads}")
    print("=" * 95)

    records: List[Dict[str, Any]] = []
    all_tree_rows: List[Dict[str, Any]] = []

    if threads > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            future_to_task = {
                executor.submit(
                    run_single_lba_trial,
                    n_taxa=task[0],
                    divergence=task[1],
                    lba_ratio=task[2],
                    rep=task[3],
                    seq_length=seq_length,
                    alpha=alpha,
                    seed=seed,
                    outdir=outdir,
                ): task
                for task in tasks
            }

            completed = 0
            for future in concurrent.futures.as_completed(future_to_task):
                task = future_to_task[future]
                completed += 1
                try:
                    res, t_rows = future.result()
                    records.append(res)
                    all_tree_rows.extend(t_rows)
                    print(
                        f"[{completed:3d}/{total_tasks:3d}] "
                        f"N={task[0]:2d}, D={task[1]:3.1f}, LBA={task[2]:3.1f}, Rep={task[3]:2d} | "
                        f"SSS={res['sss_mean']:.3f} | "
                        f"4PC(B/W)={res['delta_add_blosum62']:.3f}/{res['delta_add_wag_gamma']:.3f} | "
                        f"BLOSUM(NJ={res['blosum62_Standard_NJ_nrf']:.3f}, AddNJ2={res['blosum62_AddNJ2_nrf']:.3f}) | "
                        f"WAG(NJ={res['wag_gamma_Standard_NJ_nrf']:.3f}, AddNJ2={res['wag_gamma_AddNJ2_nrf']:.3f}) | "
                        f"TRUE(AddNJ2={res['true_AddNJ2_nrf']:.3f})"
                    )
                except Exception as exc:
                    print(f"[ERROR] Task {task} failed with exception: {exc}")
    else:
        for completed, task in enumerate(tasks, 1):
            res, t_rows = run_single_lba_trial(
                n_taxa=task[0],
                divergence=task[1],
                lba_ratio=task[2],
                rep=task[3],
                seq_length=seq_length,
                alpha=alpha,
                seed=seed,
                outdir=outdir,
            )
            records.append(res)
            all_tree_rows.extend(t_rows)
            print(
                f"[{completed:3d}/{total_tasks:3d}] "
                f"N={task[0]:2d}, D={task[1]:3.1f}, LBA={task[2]:3.1f}, Rep={task[3]:2d} | "
                f"SSS={res['sss_mean']:.3f} | "
                f"4PC(B/W)={res['delta_add_blosum62']:.3f}/{res['delta_add_wag_gamma']:.3f} | "
                f"BLOSUM(NJ={res['blosum62_Standard_NJ_nrf']:.3f}, AddNJ2={res['blosum62_AddNJ2_nrf']:.3f}) | "
                f"WAG(NJ={res['wag_gamma_Standard_NJ_nrf']:.3f}, AddNJ2={res['wag_gamma_AddNJ2_nrf']:.3f}) | "
                f"TRUE(AddNJ2={res['true_AddNJ2_nrf']:.3f})"
            )

    df = pd.DataFrame(records)
    df.sort_values(by=["n_taxa", "divergence", "lba_ratio", "replicate"], inplace=True)
    df.to_csv(csv_path, index=False)

    df_trees = pd.DataFrame(all_tree_rows)
    df_trees.sort_values(by=["n_taxa", "divergence", "lba_ratio", "replicate", "matrix_type", "method"], inplace=True)
    df_trees.to_csv(long_trees_csv_path, index=False)

    print(f"\n[OK] Saved wide summary CSV to: {csv_path}")
    print(f"[OK] Saved all {len(df_trees)} inferred trees CSV to: {long_trees_csv_path}")
    print(f"[OK] Saved individual Newick files (.nwk) and FASTA sequences under: {os.path.join(outdir, 'trees')}")
    return df, df_trees


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Comprehensive LBA Benchmark Suite.")
    parser.add_argument("--taxa", type=int, nargs="+", default=[20, 32], help="Taxa counts (default: 20 32)")
    parser.add_argument("--divergences", type=float, nargs="+", default=[0.5, 1.0, 2.0], help="Divergence scales")
    parser.add_argument("--lba_ratios", type=float, nargs="+", default=[1.0, 3.0], help="LBA branch ratios (default: 1.0 3.0)")
    parser.add_argument("--replicates", type=int, default=10, help="Replicates per condition (default: 10)")
    parser.add_argument("--length", type=int, default=500, help="Sequence length (default: 500)")
    parser.add_argument("--alpha", type=float, default=1.0, help="Gamma shape parameter (default: 1.0)")
    parser.add_argument("--outdir", type=str, default="results_lba", help="Output directory")
    parser.add_argument("--seed", type=int, default=42, help="Base random seed")
    parser.add_argument("--threads", type=int, default=6, help="Number of worker processes (default: 6)")

    args = parser.parse_args()

    run_lba_benchmark(
        taxa_list=args.taxa,
        divergences=args.divergences,
        lba_ratios=args.lba_ratios,
        replicates=args.replicates,
        seq_length=args.length,
        alpha=args.alpha,
        outdir=args.outdir,
        seed=args.seed,
        threads=args.threads,
    )
