"""
Comprehensive Benchmark Script for Long Branch Attraction (LBA) Effects.

Evaluates AddNJ, AddNJ2, Add-NNI/SPR/Hybrid, Standard NJ, BioNJ, FastME, and MSA+ML
across parameter grid of Taxa (N=20, 32), Divergences (D=0.5, 1.0, 2.0),
and LBA Ratios (1.0 vs 3.0).
"""

from typing import List, Dict, Any, Tuple, Optional
import os
import sys
import io
import time
import argparse
import concurrent.futures
import numpy as np
import pandas as pd
import Bio.Phylo as Phylo

from src.tree_simulator import generate_backward_yule_tree, compute_patristic_distance_matrix
from src.sequence_simulator import WAGModel, simulate_wag_gamma_sequences, compute_pairwise_poisson_distance_matrix
from src.tree_inference import matrix_to_nj_newick
from src.evaluation import evaluate_tree_lba
from src.add_nj import add_nj, add_nj2
from src.fastme_wrapper import run_fastme
from src.ml_wrapper import run_msa_ml
from src.topology_search import TreeScorer, hill_climb_search


def run_single_lba_trial(
    n_taxa: int,
    divergence: float,
    lba_ratio: float,
    rep: int,
    seq_length: int = 500,
    alpha: float = 1.0,
    seed: int = 42,
    max_iter_search: int = 20,
) -> Dict[str, Any]:
    """
    Executes a single simulation trial for a specific condition.

    Parameters
    ----------
    n_taxa : int
        Number of taxa [unitless].
    divergence : float
        Evolutionary distance scaling factor D [substitutions/site].
    lba_ratio : float
        LBA branch length ratio for T1 and TN [unitless].
    rep : int
        Replicate index (1-based) [unitless].
    seq_length : int, default=500
        Length of simulated amino acid sequences [amino acids].
    alpha : float, default=1.0
        Gamma shape parameter [unitless].
    seed : int, default=42
        Base random seed [unitless].
    max_iter_search : int, default=20
        Maximum iterations for topology search [unitless].

    Returns
    -------
    Dict[str, Any]
        Record dictionary containing metrics for all evaluated methods.
    """
    # Deterministic seed computation
    rep_seed = seed + n_taxa * 100000 + int(divergence * 10000) + int(lba_ratio * 1000) + rep * 100

    wag = WAGModel(wag_path="src/wag.dat", alpha=alpha)

    # 1. Generate true tree with designated lba_ratio
    tree, taxa_names = generate_backward_yule_tree(
        n_taxa=n_taxa,
        divergence=divergence,
        seed=rep_seed,
        lba_ratio=lba_ratio,
    )
    true_newick = tree.to_newick() + ";"
    lba_pair = (taxa_names[0], taxa_names[-1])

    # 2. Simulate WAG+Gamma protein sequences
    seqs = simulate_wag_gamma_sequences(
        root=tree,
        seq_length=seq_length,
        alpha=alpha,
        model=wag,
        seed=rep_seed + 1,
    )

    # 3. Global alignment & Poisson distance matrix
    t_align_start = time.time()
    D_obs = compute_pairwise_poisson_distance_matrix(seqs, taxa_names, mode="global")
    align_time = time.time() - t_align_start

    scorer = TreeScorer(D_obs, taxa_names)

    record: Dict[str, Any] = {
        "replicate": rep,
        "n_taxa": n_taxa,
        "divergence": divergence,
        "lba_ratio": lba_ratio,
        "seq_length": seq_length,
        "align_time_sec": float(align_time),
    }

    method_outputs: Dict[str, Tuple[str, float, float]] = {}

    # Method 1: MSA + ML (MAFFT + FastTree)
    t0 = time.time()
    try:
        ml_nwk = run_msa_ml(seqs)
        t_ml = time.time() - t0
        s_ml = scorer.score_tree(Phylo.read(io.StringIO(ml_nwk), "newick"))
        method_outputs["MSA_ML"] = (ml_nwk, t_ml, s_ml)
    except Exception as e:
        method_outputs["MSA_ML"] = ("", 0.0, 1.0)

    # Method 2: Standard NJ
    t0 = time.time()
    nj_nwk = matrix_to_nj_newick(D_obs, taxa_names)
    t_nj = time.time() - t0
    s_nj = scorer.score_tree(Phylo.read(io.StringIO(nj_nwk), "newick"))
    method_outputs["Standard_NJ"] = (nj_nwk, t_nj, s_nj)

    # Method 3: BioNJ
    t0 = time.time()
    bionj_nwk = run_fastme(D_obs, taxa_names, method="bionj")
    t_bionj = time.time() - t0
    s_bionj = scorer.score_tree(Phylo.read(io.StringIO(bionj_nwk), "newick"))
    method_outputs["BioNJ"] = (bionj_nwk, t_bionj, s_bionj)

    # Method 4: FastME-NNI
    t0 = time.time()
    fnni_nwk = run_fastme(D_obs, taxa_names, method="nni")
    t_fnni = time.time() - t0
    s_fnni = scorer.score_tree(Phylo.read(io.StringIO(fnni_nwk), "newick"))
    method_outputs["FastME_NNI"] = (fnni_nwk, t_fnni, s_fnni)

    # Method 5: FastME-SPR
    t0 = time.time()
    fspr_nwk = run_fastme(D_obs, taxa_names, method="spr")
    t_fspr = time.time() - t0
    s_fspr = scorer.score_tree(Phylo.read(io.StringIO(fspr_nwk), "newick"))
    method_outputs["FastME_SPR"] = (fspr_nwk, t_fspr, s_fspr)

    # Method 6: AddNJ (Method 1 original)
    t0 = time.time()
    addnj_nwk = add_nj(D_obs, taxa_names)
    t_addnj = time.time() - t0
    s_addnj = scorer.score_tree(Phylo.read(io.StringIO(addnj_nwk), "newick"))
    method_outputs["AddNJ"] = (addnj_nwk, t_addnj, s_addnj)

    # Method 7: AddNJ2 (Method 1 sister-pair consistency)
    t0 = time.time()
    addnj2_nwk = add_nj2(D_obs, taxa_names)
    t_addnj2 = time.time() - t0
    s_addnj2 = scorer.score_tree(Phylo.read(io.StringIO(addnj2_nwk), "newick"))
    method_outputs["AddNJ2"] = (addnj2_nwk, t_addnj2, s_addnj2)

    # Method 8: NJ + Add-NNI
    t0 = time.time()
    nj_nni_nwk, s_nj_nni = hill_climb_search(nj_nwk, scorer, mode="nni", max_iter=max_iter_search)
    t_nj_nni = time.time() - t0
    method_outputs["NJ_Add_NNI"] = (nj_nni_nwk, t_nj_nni, s_nj_nni)

    # Method 9: NJ + Add-SPR
    t0 = time.time()
    nj_spr_nwk, s_nj_spr = hill_climb_search(nj_nwk, scorer, mode="spr", max_iter=max_iter_search)
    t_nj_spr = time.time() - t0
    method_outputs["NJ_Add_SPR"] = (nj_spr_nwk, t_nj_spr, s_nj_spr)

    # Method 10: NJ + Add-Hybrid
    t0 = time.time()
    nj_hyb_nwk, s_nj_hyb = hill_climb_search(nj_nwk, scorer, mode="hybrid", max_iter=max_iter_search)
    t_nj_hyb = time.time() - t0
    method_outputs["NJ_Add_Hybrid"] = (nj_hyb_nwk, t_nj_hyb, s_nj_hyb)

    # Method 11: BioNJ + Add-NNI
    t0 = time.time()
    bionj_nni_nwk, s_bionj_nni = hill_climb_search(bionj_nwk, scorer, mode="nni", max_iter=max_iter_search)
    t_bionj_nni = time.time() - t0
    method_outputs["BioNJ_Add_NNI"] = (bionj_nni_nwk, t_bionj_nni, s_bionj_nni)

    # Method 12: BioNJ + Add-SPR
    t0 = time.time()
    bionj_spr_nwk, s_bionj_spr = hill_climb_search(bionj_nwk, scorer, mode="spr", max_iter=max_iter_search)
    t_bionj_spr = time.time() - t0
    method_outputs["BioNJ_Add_SPR"] = (bionj_spr_nwk, t_bionj_spr, s_bionj_spr)

    # Method 13: BioNJ + Add-Hybrid
    t0 = time.time()
    bionj_hyb_nwk, s_bionj_hyb = hill_climb_search(bionj_nwk, scorer, mode="hybrid", max_iter=max_iter_search)
    t_bionj_hyb = time.time() - t0
    method_outputs["BioNJ_Add_Hybrid"] = (bionj_hyb_nwk, t_bionj_hyb, s_bionj_hyb)

    # Method 14: AddNJ2 + Add-NNI
    t0 = time.time()
    addnj2_nni_nwk, s_addnj2_nni = hill_climb_search(addnj2_nwk, scorer, mode="nni", max_iter=max_iter_search)
    t_addnj2_nni = time.time() - t0
    method_outputs["AddNJ2_Add_NNI"] = (addnj2_nni_nwk, t_addnj2_nni, s_addnj2_nni)

    # Method 15: AddNJ2 + Add-Hybrid
    t0 = time.time()
    addnj2_hyb_nwk, s_addnj2_hyb = hill_climb_search(addnj2_nwk, scorer, mode="hybrid", max_iter=max_iter_search)
    t_addnj2_hyb = time.time() - t0
    method_outputs["AddNJ2_Add_Hybrid"] = (addnj2_hyb_nwk, t_addnj2_hyb, s_addnj2_hyb)

    # Evaluate LBA metrics and RF against true tree
    true_lba = evaluate_tree_lba(true_newick, true_newick, taxa_names, lba_pair=lba_pair)
    record["ds_true"] = float(true_lba["ds_true"])

    for m_name, (inferred_str, runtime, penalty) in method_outputs.items():
        if inferred_str:
            eval_res = evaluate_tree_lba(inferred_str, true_newick, taxa_names, lba_pair=lba_pair)
            record[f"{m_name}_nrf"] = float(eval_res["nrf"])
            record[f"{m_name}_rf"] = float(eval_res["rf"])
            record[f"{m_name}_ds_est"] = float(eval_res["ds_est"])
            record[f"{m_name}_delta_ds"] = float(eval_res["delta_ds"])
            record[f"{m_name}_is_sister"] = float(eval_res["is_sister"])
            record[f"{m_name}_penalty"] = float(penalty)
            record[f"{m_name}_time"] = float(runtime)
        else:
            record[f"{m_name}_nrf"] = 1.0
            record[f"{m_name}_rf"] = float(2 * (n_taxa - 3))
            record[f"{m_name}_ds_est"] = -1.0
            record[f"{m_name}_delta_ds"] = 0.0
            record[f"{m_name}_is_sister"] = 0.0
            record[f"{m_name}_penalty"] = 1.0
            record[f"{m_name}_time"] = 0.0

    return record


def run_lba_benchmark(
    taxa_list: List[int] = [20, 32],
    divergences: List[float] = [0.5, 1.0, 2.0],
    lba_ratios: List[float] = [1.0, 3.0],
    replicates: int = 10,
    seq_length: int = 500,
    alpha: float = 1.0,
    outdir: str = "results_lba",
    seed: int = 42,
    threads: int = 4,
) -> pd.DataFrame:
    """
    Executes the comprehensive LBA benchmark suite across the parameter grid.
    """
    os.makedirs(outdir, exist_ok=True)
    csv_path = os.path.join(outdir, "lba_summary.csv")

    tasks = []
    for n_taxa in taxa_list:
        for div in divergences:
            for lba_r in lba_ratios:
                for rep in range(1, replicates + 1):
                    tasks.append((n_taxa, div, lba_r, rep))

    total_tasks = len(tasks)
    print("=" * 80)
    print("LBA BENCHMARK SUITE: ADDITIVITY-BASED METHODS VS BASELINES & ML")
    print(f"Taxa: {taxa_list} | Divergences: {divergences} | LBA Ratios: {lba_ratios} | Replicates: {replicates}")
    print(f"Total Conditions: {total_tasks} | Worker Threads: {threads}")
    print("=" * 80)

    records: List[Dict[str, Any]] = []

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
                ): task
                for task in tasks
            }

            completed = 0
            for future in concurrent.futures.as_completed(future_to_task):
                task = future_to_task[future]
                completed += 1
                try:
                    res = future.result()
                    records.append(res)
                    print(
                        f"[{completed:3d}/{total_tasks:3d}] "
                        f"N={task[0]:2d}, D={task[1]:3.1f}, LBA={task[2]:3.1f}, Rep={task[3]:2d} | "
                        f"NJ nRF={res.get('Standard_NJ_nrf', -1):.3f} (DS={res.get('Standard_NJ_ds_est', -1):.0f}) | "
                        f"BioNJ nRF={res.get('BioNJ_nrf', -1):.3f} | "
                        f"AddNJ2 nRF={res.get('AddNJ2_nrf', -1):.3f} (DS={res.get('AddNJ2_ds_est', -1):.0f}) | "
                        f"BioNJ+AddHyb nRF={res.get('BioNJ_Add_Hybrid_nrf', -1):.3f}"
                    )
                except Exception as exc:
                    print(f"[ERROR] Task {task} failed with exception: {exc}")
    else:
        completed = 0
        for task in tasks:
            completed += 1
            res = run_single_lba_trial(
                n_taxa=task[0],
                divergence=task[1],
                lba_ratio=task[2],
                rep=task[3],
                seq_length=seq_length,
                alpha=alpha,
                seed=seed,
            )
            records.append(res)
            print(
                f"[{completed:3d}/{total_tasks:3d}] "
                f"N={task[0]:2d}, D={task[1]:3.1f}, LBA={task[2]:3.1f}, Rep={task[3]:2d} | "
                f"NJ nRF={res.get('Standard_NJ_nrf', -1):.3f} (DS={res.get('Standard_NJ_ds_est', -1):.0f}) | "
                f"BioNJ nRF={res.get('BioNJ_nrf', -1):.3f} | "
                f"AddNJ2 nRF={res.get('AddNJ2_nrf', -1):.3f} (DS={res.get('AddNJ2_ds_est', -1):.0f}) | "
                f"BioNJ+AddHyb nRF={res.get('BioNJ_Add_Hybrid_nrf', -1):.3f}"
            )

    df = pd.DataFrame(records)
    # Sort deterministically
    df.sort_values(by=["n_taxa", "divergence", "lba_ratio", "replicate"], inplace=True)
    df.to_csv(csv_path, index=False)
    print(f"\n[OK] Benchmark completed! Saved results to {csv_path}")
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run LBA Benchmark Suite for Additivity-Based Methods.")
    parser.add_argument("--taxa", type=int, nargs="+", default=[20, 32], help="Taxa counts (default: 20 32)")
    parser.add_argument("--divergences", type=float, nargs="+", default=[0.5, 1.0, 2.0], help="Divergence scales")
    parser.add_argument("--lba_ratios", type=float, nargs="+", default=[1.0, 3.0], help="LBA branch ratios (default: 1.0 3.0)")
    parser.add_argument("--replicates", type=int, default=10, help="Replicates per condition (default: 10)")
    parser.add_argument("--length", type=int, default=500, help="Sequence length (default: 500)")
    parser.add_argument("--alpha", type=float, default=1.0, help="Gamma shape parameter (default: 1.0)")
    parser.add_argument("--outdir", type=str, default="results_lba", help="Output directory")
    parser.add_argument("--seed", type=int, default=42, help="Base random seed")
    parser.add_argument("--threads", type=int, default=4, help="Number of worker processes (default: 4)")

    args = parser.parse_args()

    df_res = run_lba_benchmark(
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
