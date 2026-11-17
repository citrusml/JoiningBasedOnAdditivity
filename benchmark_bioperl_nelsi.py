"""
Comprehensive Benchmark Script: BioPerl Backward Yule vs. NELSI Autocorrelated Trees.

Evaluates:
1. Tree Models: BioPerl (Backward Yule) and NELSI (Autocorrelated Kishino).
2. Taxa N in {20, 32, 64}, Divergence D in {0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0},
   Sequence Length L = 500 aa, Zipfian Indel Rate = 0.1, Replicates = 50.
3. General Methods:
   - MSA_ML (MAFFT + FastTree WAG+Gamma)
   - SSS_Linear_NJ, SSS_Log_NJ, Poisson_NJ, MSA_NJ
   - GS_Poisson, GS_SSS_Linear, GS_SSS_Log
4. Additional SSS Linear Methods:
   - SSS_Linear_BioNJ, SSS_Linear_FastME_NNI, SSS_Linear_AddNJ3
   - SSS_Linear_NJ_AddHybrid2, SSS_Linear_BioNJ_AddHybrid2, SSS_Linear_AddNJ3_AddHybrid2
5. Evaluation Metrics:
   - Input distance matrix RMSE (`matrix_rmse`) and tree patristic RMSE (`patristic_rmse`)
   - Topology accuracy (`rf`, `nrf`, `kf_bsd`)
   - Root-to-tip distance statistics (mean, min, max, max/min ratio, variance)
     for both true and estimated trees
   - LBA metrics (`ds_est`, `ds_true`, `delta_ds`, `is_sister`, `nta`, `pair_nrf`,
     `pair_split_recall`, `lba_quartet_attraction`, `lba_quartet_error`, `matrix_lba_margin`)
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

from src.tree_simulator import (
    generate_backward_yule_tree,
    generate_nelsi_tree,
    compute_patristic_distance_matrix,
)
from src.sequence_simulator import (
    WAGModel,
    simulate_wag_gamma_sequences,
    compute_pairwise_poisson_distance_matrix,
    compute_mafft_msa_poisson_distance_matrix,
    compute_mmseqs_similarity_matrix,
    similarity_to_distance_matrix,
)
from src.tree_inference import matrix_to_nj_newick
from src.evaluation import (
    evaluate_matrix_accuracy,
    evaluate_tree_lba,
    compute_root_to_tip_metrics,
    find_max_distance_pair,
    compute_matrix_lba_margin,
)
from src.add_nj import add_nj3
from src.fastme_wrapper import run_fastme
from src.ml_wrapper import run_msa_ml
from src.gs_wrapper import run_gs_from_distance_matrix
from src.topology_search import (
    TreeDistanceScorer,
    hill_climb_search,
    fit_ols_branch_lengths,
    tree_to_newick,
)

# Ensure line-buffered stdout/stderr for real-time monitoring
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(line_buffering=True)


def _fit_gs_branch_lengths(
    gs_newick: str,
    D_mat: np.ndarray,
    taxa_names: List[str],
) -> str:
    """
    Fits OLS branch lengths onto a topology-only Graph Splitting Newick tree
    so patristic RMSE and root-to-tip metrics can be evaluated.

    Parameters
    ----------
    gs_newick : str
        Topology-only Newick string returned by `gs2` [unitless].
    D_mat : np.ndarray
        Distance matrix of shape (n_taxa, n_taxa) [substitutions/site].
    taxa_names : List[str]
        Ordered list of taxon labels [unitless].

    Returns
    -------
    str
        Newick string with fitted OLS branch lengths [substitutions/site].
    """
    tree = Phylo.read(io.StringIO(gs_newick), "newick")
    tree.root_with_outgroup(taxa_names[0])
    fit_ols_branch_lengths(tree, D_mat, taxa_names)
    return tree_to_newick(tree)


def _evaluate_and_build_records(
    meta: Dict[str, Any],
    method_entries: List[Tuple[str, str, Optional[str], Optional[np.ndarray], float]],
    true_newick: str,
    D_true: np.ndarray,
    taxa_names: List[str],
    lba_pair: Tuple[str, str],
    true_rtt: Dict[str, float],
) -> List[Dict[str, Any]]:
    """
    Computes all matrix, topology, root-to-tip, and LBA metrics for each method
    in a single replicate.

    Parameters
    ----------
    meta : Dict[str, Any]
        Replicate metadata dictionary [unitless].
    method_entries : List[Tuple[str, str, Optional[str], Optional[np.ndarray], float]]
        List of `(method_name, est_newick, matrix_type, D_mat, runtime_sec)`.
    true_newick : str
        Ground-truth Newick string [unitless].
    D_true : np.ndarray
        Ground-truth patristic distance matrix [substitutions/site].
    taxa_names : List[str]
        Ordered list of taxon names [unitless].
    lba_pair : Tuple[str, str]
        Maximum patristic distance taxon pair `(u, v)` on the true tree [unitless].
    true_rtt : Dict[str, float]
        Root-to-tip metrics of the true tree [substitutions/site].

    Returns
    -------
    List[Dict[str, Any]]
        List of flat record dictionaries (one per method).
    """
    name_to_idx = {name: i for i, name in enumerate(taxa_names)}
    idx_u = name_to_idx[lba_pair[0]]
    idx_v = name_to_idx[lba_pair[1]]

    # Precompute matrix-level metrics once per matrix type
    matrix_cache: Dict[str, Tuple[float, float]] = {}
    records: List[Dict[str, Any]] = []

    for method_name, est_nwk, mat_type, D_mat, runtime_sec in method_entries:
        if mat_type is not None and D_mat is not None:
            if mat_type not in matrix_cache:
                mat_acc = evaluate_matrix_accuracy(D_mat, D_true)
                mat_margin = compute_matrix_lba_margin(
                    D_mat, idx_u, idx_v, true_newick, taxa_names
                )
                matrix_cache[mat_type] = (mat_acc["rmse"], mat_margin)
            mat_rmse, mat_margin = matrix_cache[mat_type]
        else:
            mat_rmse = float("nan")
            mat_margin = float("nan")

        lba_res = evaluate_tree_lba(
            est_nwk, true_newick, taxa_names, lba_pair=lba_pair
        )
        est_rtt = compute_root_to_tip_metrics(
            est_nwk, taxa_names, use_midpoint_rooting=True
        )

        rec: Dict[str, Any] = dict(meta)
        rec.update(
            {
                "method": method_name,
                "matrix_type": mat_type if mat_type is not None else "msa_ml",
                "runtime_sec": float(runtime_sec),
                "matrix_rmse": float(mat_rmse),
                "matrix_lba_margin": float(mat_margin),
                "rf": float(lba_res["rf"]),
                "nrf": float(lba_res["nrf"]),
                "kf_bsd": float(lba_res["kf_bsd"]),
                "patristic_rmse": float(lba_res["patristic_rmse"]),
                "ds_est": float(lba_res["ds_est"]),
                "ds_true": float(lba_res["ds_true"]),
                "delta_ds": float(lba_res["delta_ds"]),
                "is_sister": float(lba_res["is_sister"]),
                "nta": float(lba_res["nta"]),
                "pair_nrf": float(lba_res["pair_nrf"]),
                "pair_split_recall": float(lba_res["pair_split_recall"]),
                "pair_split_fnr": float(lba_res["pair_split_fnr"]),
                "lba_quartet_attraction": float(lba_res["lba_quartet_attraction"]),
                "lba_quartet_error": float(lba_res["lba_quartet_error"]),
                "root_dist_mean": float(est_rtt["root_dist_mean"]),
                "root_dist_var": float(est_rtt["root_dist_var"]),
                "root_dist_min": float(est_rtt["root_dist_min"]),
                "root_dist_max": float(est_rtt["root_dist_max"]),
                "root_dist_max_min_ratio": float(est_rtt["root_dist_max_min_ratio"]),
                "true_root_dist_mean": float(true_rtt["root_dist_mean"]),
                "true_root_dist_var": float(true_rtt["root_dist_var"]),
                "true_root_dist_min": float(true_rtt["root_dist_min"]),
                "true_root_dist_max": float(true_rtt["root_dist_max"]),
                "true_root_dist_max_min_ratio": float(true_rtt["root_dist_max_min_ratio"]),
            }
        )
        records.append(rec)

    return records


def run_single_replicate(task: Tuple[str, int, float, float, int, int, float, int, int, str]) -> List[Dict[str, Any]]:
    """
    Executes one complete simulation and reconstruction replicate across all methods.

    Parameters
    ----------
    task : Tuple[str, int, float, float, int, int, float, int, int, str]
        Tuple `(tree_model, n_taxa, divergence, lba_ratio, rep, seq_length, indel_rate, base_seed, max_iter_search, methods_set)`.

    Returns
    -------
    List[Dict[str, Any]]
        List of evaluation records (one per method).
    """
    tree_model, n_taxa, divergence, lba_ratio, rep, seq_length, indel_rate, base_seed, max_iter_search, methods_set = task
    model_offset = 0 if tree_model == "bioperl" else 500000
    lba_offset = int(round(lba_ratio * 100000))
    rep_seed = base_seed + model_offset + lba_offset + n_taxa * 10000 + int(round(divergence * 1000)) + rep

    # 1. Generate synthetic true tree
    if tree_model == "bioperl":
        root, taxa_names = generate_backward_yule_tree(
            n_taxa=n_taxa, divergence=divergence, seed=rep_seed, lba_ratio=lba_ratio
        )
    elif tree_model == "nelsi":
        root, taxa_names = generate_nelsi_tree(
            n_taxa=n_taxa, divergence=divergence, seed=rep_seed, lba_ratio=lba_ratio
        )
    else:
        raise ValueError(f"Unsupported tree_model '{tree_model}'.")

    true_newick = root.to_newick() + ";"
    D_true = compute_patristic_distance_matrix(root, taxa_names)
    if lba_ratio > 1.0:
        lba_pair = (taxa_names[0], taxa_names[-1])
    else:
        lba_pair = find_max_distance_pair(true_newick, taxa_names, D_true=D_true)
    true_rtt = compute_root_to_tip_metrics(true_newick, taxa_names, use_midpoint_rooting=False)

    # 2. Simulate WAG+Gamma sequences with Zipfian indels (INDELible)
    wag = WAGModel(wag_path="src/wag.dat", alpha=1.0)
    sequences = simulate_wag_gamma_sequences(
        root=root,
        seq_length=seq_length,
        alpha=1.0,
        model=wag,
        seed=rep_seed + 1,
        indel_rate=indel_rate,
    )

    # 3. Compute distance matrices (Poisson PSA, MSA Poisson, SSS Linear, SSS Log)
    D_poisson = compute_pairwise_poisson_distance_matrix(sequences, taxa_names, mode="global")
    D_msa_poisson = compute_mafft_msa_poisson_distance_matrix(sequences, taxa_names)
    mmseqs_res = compute_mmseqs_similarity_matrix(sequences, taxa_names, threads=1)
    if mmseqs_res is None:
        raise RuntimeError("MMseqs2 execution failed or mmseqs binary is unavailable.")
    W_mmseqs, _ = mmseqs_res
    D_sss_linear = similarity_to_distance_matrix(W_mmseqs, dist_model="linear")
    D_sss_log = similarity_to_distance_matrix(W_mmseqs, dist_model="log")

    method_entries: List[Tuple[str, str, Optional[str], Optional[np.ndarray], float]] = []

    # --- General Methods ---
    # 1) MSA + ML (MAFFT + FastTree WAG+Gamma)
    t0 = time.time()
    msa_ml_nwk = run_msa_ml(sequences)
    method_entries.append(("MSA_ML", msa_ml_nwk, None, None, time.time() - t0))

    # 2) SSS Linear NJ
    t0 = time.time()
    sss_lin_nj_nwk = matrix_to_nj_newick(D_sss_linear, taxa_names)
    method_entries.append(("SSS_Linear_NJ", sss_lin_nj_nwk, "sss_linear", D_sss_linear, time.time() - t0))

    # 3) SSS Log NJ
    t0 = time.time()
    sss_log_nj_nwk = matrix_to_nj_newick(D_sss_log, taxa_names)
    method_entries.append(("SSS_Log_NJ", sss_log_nj_nwk, "sss_log", D_sss_log, time.time() - t0))

    # 4) Poisson NJ (BLOSUM62 PSA + Poisson)
    t0 = time.time()
    poisson_nj_nwk = matrix_to_nj_newick(D_poisson, taxa_names)
    method_entries.append(("Poisson_NJ", poisson_nj_nwk, "poisson", D_poisson, time.time() - t0))

    # 5) MSA NJ (MAFFT MSA + Poisson)
    t0 = time.time()
    msa_nj_nwk = matrix_to_nj_newick(D_msa_poisson, taxa_names)
    method_entries.append(("MSA_NJ", msa_nj_nwk, "msa_poisson", D_msa_poisson, time.time() - t0))

    # 6) GS (Poisson)
    t0 = time.time()
    gs_pois_raw = run_gs_from_distance_matrix(D_poisson, taxa_names)
    gs_pois_nwk = _fit_gs_branch_lengths(gs_pois_raw, D_poisson, taxa_names)
    method_entries.append(("GS_Poisson", gs_pois_nwk, "poisson", D_poisson, time.time() - t0))

    # 7) GS (SSS Linear)
    t0 = time.time()
    gs_lin_raw = run_gs_from_distance_matrix(D_sss_linear, taxa_names)
    gs_lin_nwk = _fit_gs_branch_lengths(gs_lin_raw, D_sss_linear, taxa_names)
    method_entries.append(("GS_SSS_Linear", gs_lin_nwk, "sss_linear", D_sss_linear, time.time() - t0))

    # 8) GS (SSS Log)
    t0 = time.time()
    gs_log_raw = run_gs_from_distance_matrix(D_sss_log, taxa_names)
    gs_log_nwk = _fit_gs_branch_lengths(gs_log_raw, D_sss_log, taxa_names)
    method_entries.append(("GS_SSS_Log", gs_log_nwk, "sss_log", D_sss_log, time.time() - t0))

    # --- Additional SSS Linear Methods (only when methods_set == "all") ---
    if methods_set == "all":
        # 9) SSS Linear BioNJ
        t0 = time.time()
        sss_lin_bionj_nwk = run_fastme(D_sss_linear, taxa_names, method="bionj")
        method_entries.append(("SSS_Linear_BioNJ", sss_lin_bionj_nwk, "sss_linear", D_sss_linear, time.time() - t0))

        # 10) SSS Linear FastME (NNI)
        t0 = time.time()
        sss_lin_fnni_nwk = run_fastme(D_sss_linear, taxa_names, method="nni")
        method_entries.append(("SSS_Linear_FastME_NNI", sss_lin_fnni_nwk, "sss_linear", D_sss_linear, time.time() - t0))

        # 11) SSS Linear AddNJ3
        t0 = time.time()
        sss_lin_addnj3_nwk = add_nj3(D_sss_linear, taxa_names)
        method_entries.append(("SSS_Linear_AddNJ3", sss_lin_addnj3_nwk, "sss_linear", D_sss_linear, time.time() - t0))

        dist_scorer = TreeDistanceScorer(D_sss_linear, taxa_names, mode="ols_rmse")

        # 12) SSS Linear Standard NJ + Add-Hybrid2
        t0 = time.time()
        sss_lin_nj_hyb2_nwk, _ = hill_climb_search(
            sss_lin_nj_nwk, dist_scorer, mode="hybrid", max_iter=max_iter_search
        )
        method_entries.append(("SSS_Linear_NJ_AddHybrid2", sss_lin_nj_hyb2_nwk, "sss_linear", D_sss_linear, time.time() - t0))

        # 13) SSS Linear BioNJ + Add-Hybrid2
        t0 = time.time()
        sss_lin_bionj_hyb2_nwk, _ = hill_climb_search(
            sss_lin_bionj_nwk, dist_scorer, mode="hybrid", max_iter=max_iter_search
        )
        method_entries.append(("SSS_Linear_BioNJ_AddHybrid2", sss_lin_bionj_hyb2_nwk, "sss_linear", D_sss_linear, time.time() - t0))

        # 14) SSS Linear AddNJ3 + Add-Hybrid2
        t0 = time.time()
        sss_lin_addnj3_hyb2_nwk, _ = hill_climb_search(
            sss_lin_addnj3_nwk, dist_scorer, mode="hybrid", max_iter=max_iter_search
        )
        method_entries.append(("SSS_Linear_AddNJ3_AddHybrid2", sss_lin_addnj3_hyb2_nwk, "sss_linear", D_sss_linear, time.time() - t0))

    meta = {
        "tree_model": tree_model,
        "lba_ratio": lba_ratio,
        "n_taxa": n_taxa,
        "divergence": divergence,
        "seq_length": seq_length,
        "indel_rate": indel_rate,
        "replicate": rep,
        "seed": rep_seed,
        "lba_pair_u": lba_pair[0],
        "lba_pair_v": lba_pair[1],
    }

    return _evaluate_and_build_records(
        meta=meta,
        method_entries=method_entries,
        true_newick=true_newick,
        D_true=D_true,
        taxa_names=taxa_names,
        lba_pair=lba_pair,
        true_rtt=true_rtt,
    )


def build_task_list(
    tree_models: List[str],
    taxa_list: List[int],
    divergences: List[float],
    lba_ratios: List[float],
    replicates: int,
    seq_length: int,
    indel_rate: float,
    base_seed: int,
    max_iter_search: int,
    methods_set: str,
) -> List[Tuple[str, int, float, float, int, int, float, int, int, str]]:
    """Builds the ordered list of all simulation replicate tasks."""
    tasks: List[Tuple[str, int, float, float, int, int, float, int, int, str]] = []
    for tm in tree_models:
        for lba_r in lba_ratios:
            for n_taxa in taxa_list:
                for div in divergences:
                    for rep in range(1, replicates + 1):
                        tasks.append(
                            (tm, n_taxa, div, lba_r, rep, seq_length, indel_rate, base_seed, max_iter_search, methods_set)
                        )
    return tasks


def main() -> None:
    parser = argparse.ArgumentParser(
        description="BioPerl vs. NELSI Comprehensive Phylogenetic Benchmark."
    )
    parser.add_argument(
        "--tree-models",
        type=str,
        nargs="+",
        default=["bioperl", "nelsi"],
        choices=["bioperl", "nelsi"],
        help="Tree simulation models (default: bioperl nelsi)",
    )
    parser.add_argument(
        "--taxa",
        type=int,
        nargs="+",
        default=[20, 32, 64],
        help="Taxa counts (default: 20 32 64)",
    )
    parser.add_argument(
        "--divergences",
        type=float,
        nargs="+",
        default=[0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0],
        help="Divergence scales D (default: 0.1 0.2 0.5 1.0 1.5 2.0 3.0)",
    )
    parser.add_argument(
        "--lba-ratios",
        type=float,
        nargs="+",
        default=[1.0],
        help="LBA terminal branch multipliers on T1 and TN (default: 1.0)",
    )
    parser.add_argument(
        "--methods-set",
        type=str,
        default="all",
        choices=["all", "core"],
        help="Method set: 'core' (8 general methods) or 'all' (adds 6 sss_linear methods) (default: all)",
    )
    parser.add_argument(
        "--length",
        type=int,
        default=500,
        help="Sequence length in amino acids (default: 500)",
    )
    parser.add_argument(
        "--indel-rate",
        type=float,
        default=0.1,
        help="Zipfian indel rate relative to substitution rate (default: 0.1)",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=50,
        help="Number of replicates per condition (default: 50)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of parallel worker processes (default: 4)",
    )
    parser.add_argument(
        "--max-iter-search",
        type=int,
        default=20,
        help="Max iterations per hill-climb phase for Add-Hybrid2 (default: 20)",
    )
    parser.add_argument(
        "--shard-idx",
        type=int,
        default=0,
        help="0-based shard index for HPC array jobs (default: 0)",
    )
    parser.add_argument(
        "--num-shards",
        type=int,
        default=1,
        help="Total number of shards for HPC array jobs (default: 1)",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="results_bioperl_nelsi",
        help="Output directory (default: results_bioperl_nelsi)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Base pseudo-random seed (default: 42)",
    )
    parser.add_argument(
        "--merge-only",
        action="store_true",
        help="Merge shard CSV files under <outdir>/shards/ into benchmark_results.csv and summary_by_condition.csv",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    shards_dir = os.path.join(args.outdir, "shards")

    if args.merge_only:
        if not os.path.isdir(shards_dir):
            raise FileNotFoundError(f"Shards directory not found: {shards_dir}")
        shard_files = sorted(
            os.path.join(shards_dir, f)
            for f in os.listdir(shards_dir)
            if f.endswith(".csv")
        )
        if not shard_files:
            raise FileNotFoundError(f"No shard CSV files found in: {shards_dir}")
        dfs = [pd.read_csv(f) for f in shard_files]
        df_merged = pd.concat(dfs, ignore_index=True)
        out_csv = os.path.join(args.outdir, "benchmark_results.csv")
        df_merged.to_csv(out_csv, index=False)
        group_cols = ["tree_model", "lba_ratio", "n_taxa", "divergence", "method"] if "lba_ratio" in df_merged.columns else ["tree_model", "n_taxa", "divergence", "method"]
        summary = (
            df_merged.groupby(group_cols, as_index=False)
            .mean(numeric_only=True)
        )
        summary_csv = os.path.join(args.outdir, "summary_by_condition.csv")
        summary.to_csv(summary_csv, index=False)
        print(f"[OK] Merged {len(shard_files)} shards -> {out_csv} ({len(df_merged)} rows)")
        print(f"[OK] Saved condition-level summary -> {summary_csv}")
        return

    all_tasks = build_task_list(
        tree_models=args.tree_models,
        taxa_list=args.taxa,
        divergences=args.divergences,
        lba_ratios=args.lba_ratios,
        replicates=args.replicates,
        seq_length=args.length,
        indel_rate=args.indel_rate,
        base_seed=args.seed,
        max_iter_search=args.max_iter_search,
        methods_set=args.methods_set,
    )

    if args.num_shards > 1:
        os.makedirs(shards_dir, exist_ok=True)
        tasks = [t for i, t in enumerate(all_tasks) if i % args.num_shards == args.shard_idx]
        m_tag = "_".join(args.tree_models)
        lba_tag = "_".join(str(x) for x in args.lba_ratios)
        t_tag = "_".join(str(x) for x in args.taxa)
        d_tag = "_".join(str(x) for x in args.divergences)
        out_csv = os.path.join(
            shards_dir,
            f"shard_{m_tag}_LBA{lba_tag}_N{t_tag}_D{d_tag}_s{args.shard_idx:03d}.csv",
        )
    else:
        tasks = all_tasks
        out_csv = os.path.join(args.outdir, "benchmark_results.csv")

    print("=" * 80)
    print("BIOPERL VS. NELSI COMPREHENSIVE BENCHMARK")
    print(f"Models: {args.tree_models} | LBA: {args.lba_ratios} | Taxa: {args.taxa} | D: {args.divergences}")
    print(f"L={args.length} aa | Indel={args.indel_rate} | Reps={args.replicates} | Methods={args.methods_set} | Threads={args.threads}")
    print(f"Shard {args.shard_idx}/{args.num_shards}: executing {len(tasks)}/{len(all_tasks)} tasks")
    print("=" * 80)

    all_records: List[Dict[str, Any]] = []
    if args.threads > 1 and len(tasks) > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
            future_to_task = {ex.submit(run_single_replicate, t): t for t in tasks}
            done_cnt = 0
            for fut in concurrent.futures.as_completed(future_to_task):
                done_cnt += 1
                t = future_to_task[fut]
                recs = fut.result()
                all_records.extend(recs)
                print(
                    f"[{done_cnt}/{len(tasks)}] model={t[0]} N={t[1]} D={t[2]} LBA={t[3]} rep={t[4]} completed"
                )
    else:
        for idx, t in enumerate(tasks, start=1):
            recs = run_single_replicate(t)
            all_records.extend(recs)
            print(f"[{idx}/{len(tasks)}] model={t[0]} N={t[1]} D={t[2]} LBA={t[3]} rep={t[4]} completed")

    df = pd.DataFrame(all_records)
    df.to_csv(out_csv, index=False)
    print(f"\n[OK] Saved results to: {out_csv} ({len(df)} rows)")

    if args.num_shards == 1 and not df.empty:
        summary = (
            df.groupby(["tree_model", "lba_ratio", "n_taxa", "divergence", "method"], as_index=False)
            .mean(numeric_only=True)
        )
        summary_csv = os.path.join(args.outdir, "summary_by_condition.csv")
        summary.to_csv(summary_csv, index=False)
        print(f"[OK] Saved condition-level summary to: {summary_csv}")


if __name__ == "__main__":
    main()
