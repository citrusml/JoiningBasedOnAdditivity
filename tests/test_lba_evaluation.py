"""
Unit tests for LBA evaluation and topological shortest path distance.
"""

import numpy as np
import pytest
from src.evaluation import (
    compute_topological_distance_edges,
    evaluate_tree_lba,
    find_max_distance_pair,
    compute_matrix_lba_margin,
    compute_root_to_tip_metrics,
)
from src.gs_wrapper import run_gs_from_distance_matrix


def test_compute_topological_distance_edges():
    taxa = ["T1", "T2", "T3", "T4"]
    # Tree 1: ((T1, T2), (T3, T4));
    # Path T1 -> internal1 -> internal2 -> T4 has length 3 edges.
    tree1 = "((T1:1.0, T2:0.1):0.1, (T3:0.1, T4:1.0):0.1);"
    d1_4 = compute_topological_distance_edges(tree1, "T1", "T4", taxa)
    d1_2 = compute_topological_distance_edges(tree1, "T1", "T2", taxa)
    assert d1_4 == 3
    assert d1_2 == 2

    # Tree 2: ((T1, T4), (T2, T3));
    # Here T1 and T4 are sister taxa -> 2 edges.
    tree2 = "((T1:1.0, T4:1.0):0.1, (T2:0.1, T3:0.1):0.1);"
    d2_1_4 = compute_topological_distance_edges(tree2, "T1", "T4", taxa)
    assert d2_1_4 == 2


def test_evaluate_tree_lba():
    taxa = ["T1", "T2", "T3", "T4"]
    true_tree = "((T1:1.0, T2:0.1):0.1, (T3:0.1, T4:1.0):0.1);"
    lba_tree = "((T1:1.0, T4:1.0):0.1, (T2:0.1, T3:0.1):0.1);"

    # Evaluation when estimated tree is identical to true tree
    res_correct = evaluate_tree_lba(true_tree, true_tree, taxa, lba_pair=("T1", "T4"))
    assert res_correct["rf"] == 0.0
    assert res_correct["nrf"] == 0.0
    assert res_correct["ds_est"] == 3.0
    assert res_correct["ds_true"] == 3.0
    assert res_correct["delta_ds"] == 0.0
    assert res_correct["is_sister"] == 0.0
    assert res_correct["nta"] == 0.0
    assert res_correct["pair_nrf"] == 0.0
    assert res_correct["pair_split_recall"] == 1.0
    assert res_correct["pair_split_fnr"] == 0.0
    assert res_correct["lba_quartet_attraction"] == 0.0
    assert res_correct["lba_quartet_error"] == 0.0

    # Evaluation when estimated tree suffers from LBA (T1 and T4 sister)
    res_lba = evaluate_tree_lba(lba_tree, true_tree, taxa, lba_pair=("T1", "T4"))
    assert res_lba["rf"] > 0.0
    assert res_lba["ds_est"] == 2.0
    assert res_lba["ds_true"] == 3.0
    assert res_lba["delta_ds"] == 1.0
    assert res_lba["is_sister"] == 1.0
    assert res_lba["nta"] == 1.0
    assert res_lba["pair_nrf"] == 1.0
    assert res_lba["pair_split_recall"] == 0.0
    assert res_lba["pair_split_fnr"] == 1.0
    assert res_lba["lba_quartet_attraction"] == 1.0
    assert res_lba["lba_quartet_error"] == 1.0


def test_lba_metrics_local_noise_invariance():
    """
    Verifies that when a non-LBA side-branch rearrangement occurs without attracting
    the long-branch pair (T1, T6), pair_nrf reflects the split change while
    lba_quartet_attraction and nta remain 0.0.
    """
    taxa = ["T1", "T2", "T3", "T4", "T5", "T6"]
    true_tree = "(((T1:1.0, T2:0.1):0.2, T3:0.1):0.2, (T4:0.1, (T5:0.1, T6:1.0):0.2):0.2);"
    # Swap T3 and T4 across the middle split, keeping T1 and T6 separated at both ends
    side_noise_tree = "(((T1:1.0, T2:0.1):0.2, T4:0.1):0.2, (T3:0.1, (T5:0.1, T6:1.0):0.2):0.2);"

    res_noise = evaluate_tree_lba(side_noise_tree, true_tree, taxa, auto_max_pair=True)
    assert res_noise["ds_est"] == res_noise["ds_true"]
    assert res_noise["nta"] == 0.0
    assert res_noise["pair_nrf"] > 0.0
    assert res_noise["lba_quartet_attraction"] == 0.0


def test_find_max_distance_pair_and_matrix_margin():
    taxa = ["T1", "T2", "T3", "T4", "T5", "T6"]
    true_tree = "(((T1:1.0, T2:0.1):0.2, T3:0.1):0.2, (T4:0.1, (T5:0.1, T6:1.0):0.2):0.2);"
    u, v = find_max_distance_pair(true_tree, taxa)
    assert {u, v} == {"T1", "T6"}

    # Construct matrix where T1 and T6 are strongly attracted (short D_16)
    D_attracted = np.full((6, 6), 2.0, dtype=np.float64)
    np.fill_diagonal(D_attracted, 0.0)
    D_attracted[0, 5] = D_attracted[5, 0] = 0.2
    margin = compute_matrix_lba_margin(D_attracted, 0, 5, true_tree, taxa)
    assert margin < 0.0


def test_root_to_tip_metrics_and_gs_wrapper():
    taxa = ["T1", "T2", "T3", "T4"]
    # Ultrametric tree: all root-to-tip distances == 0.5
    ultrametric = "((T1:0.2, T2:0.2):0.3, (T3:0.1, T4:0.1):0.4);"
    rtt = compute_root_to_tip_metrics(ultrametric, taxa, use_midpoint_rooting=False)
    assert abs(rtt["root_dist_mean"] - 0.5) < 1e-6
    assert abs(rtt["root_dist_var"]) < 1e-12
    assert abs(rtt["root_dist_max_min_ratio"] - 1.0) < 1e-6

    # GS wrapper on simple distance matrix
    D = np.array([
        [0.0, 0.2, 0.8, 0.9],
        [0.2, 0.0, 0.8, 0.9],
        [0.8, 0.8, 0.0, 0.3],
        [0.9, 0.9, 0.3, 0.0],
    ], dtype=np.float64)
    gs_nwk = run_gs_from_distance_matrix(D, taxa)
    assert gs_nwk.endswith(";") and "T1" in gs_nwk and "T4" in gs_nwk
