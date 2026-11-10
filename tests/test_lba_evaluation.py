"""
Unit tests for LBA evaluation and topological shortest path distance.
"""

import pytest
from src.evaluation import compute_topological_distance_edges, evaluate_tree_lba


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

    # Evaluation when estimated tree suffers from LBA (T1 and T4 sister)
    res_lba = evaluate_tree_lba(lba_tree, true_tree, taxa, lba_pair=("T1", "T4"))
    assert res_lba["rf"] > 0.0
    assert res_lba["ds_est"] == 2.0
    assert res_lba["ds_true"] == 3.0
    assert res_lba["delta_ds"] == 1.0
    assert res_lba["is_sister"] == 1.0
