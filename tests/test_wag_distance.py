"""
Unit tests for WAG+Gamma pairwise Maximum Likelihood distance estimation
and multi-matrix tree quartet additivity violation scoring.
"""

import math
import numpy as np
from src.tree_simulator import generate_backward_yule_tree, compute_patristic_distance_matrix
from src.sequence_simulator import (
    WAGModel,
    simulate_wag_gamma_sequences,
    compute_pairwise_wag_gamma_distance_matrix,
    compute_both_pairwise_distance_matrices,
    compute_all_pairwise_matrices_and_sss,
)

from src.topology_search import TreeScorer
from src.evaluation import evaluate_tree_rf
from src.add_nj import add_nj2


def test_wag_gamma_distance_identical_sequences() -> None:
    """Identical sequences must have estimated WAG+Gamma ML distance of 0.0."""
    wag = WAGModel(wag_path="src/wag.dat", alpha=1.0)
    seqs = {
        "T1": "ACDEFGHIKLMNPQRSTVWY" * 10,
        "T2": "ACDEFGHIKLMNPQRSTVWY" * 10,
    }
    D_wag = compute_pairwise_wag_gamma_distance_matrix(
        seqs, ["T1", "T2"], model=wag, align=False
    )
    assert D_wag.shape == (2, 2)
    assert np.allclose(D_wag, 0.0, atol=1e-6)


def test_wag_gamma_distance_accuracy_and_both_matrices() -> None:
    """
    Simulated sequences under WAG+Gamma should yield valid, symmetric
    BLOSUM62+Poisson and WAG+Gamma ML distance matrices, and true tree
    additivity violation on true patristic matrix must be 0.0.
    """
    tree, taxa_names = generate_backward_yule_tree(n_taxa=6, divergence=0.8, seed=123)
    true_nwk = tree.to_newick() + ";"
    D_true = compute_patristic_distance_matrix(tree, taxa_names)

    wag = WAGModel(wag_path="src/wag.dat", alpha=1.0)
    seqs = simulate_wag_gamma_sequences(
        root=tree, seq_length=1000, alpha=1.0, model=wag, seed=456
    )

    D_blosum, D_wag = compute_both_pairwise_distance_matrices(
        seqs, taxa_names, model=wag, alpha=1.0
    )

    assert D_blosum.shape == (6, 6)
    assert D_wag.shape == (6, 6)
    assert np.allclose(D_blosum, D_blosum.T)
    assert np.allclose(D_wag, D_wag.T)
    assert np.all(np.diag(D_wag) == 0.0)
    assert np.all(D_wag >= 0.0)

    # Check correlation with true patristic distances
    triu = np.triu_indices(6, k=1)
    corr_wag = np.corrcoef(D_true[triu], D_wag[triu])[0, 1]
    assert corr_wag > 0.90

    # Check TreeScorer on D_true and D_wag
    scorer_true = TreeScorer(D_true, taxa_names)
    scorer_wag = TreeScorer(D_wag, taxa_names)

    assert scorer_true.compute_matrix_defect() < 1e-7
    choices_true = scorer_true.extract_tree_quartet_choices(true_nwk)
    assert scorer_true.score_choices(choices_true) < 1e-7
    assert 0.0 <= scorer_wag.score_choices(choices_true) <= 1.0

    # Check SSS calculation
    D_b2, D_w2, sss_stats = compute_all_pairwise_matrices_and_sss(
        seqs, taxa_names, model=wag, alpha=1.0
    )
    assert np.allclose(D_blosum, D_b2)
    assert np.allclose(D_wag, D_w2)
    assert "sss_mean" in sss_stats
    assert 0.0 <= sss_stats["sss_mean"] <= 1.0
    assert 0.0 <= sss_stats["seq_identity_mean"] <= 1.0


def test_indelible_zipfian_indel_and_kf_patristic_rmse() -> None:
    """
    Verifies that INDELible Zipfian indel simulation (POW 1.7 50) works for
    indel_rate in {0.05, 0.1}, and that evaluate_tree_rf computes kf_bsd and patristic_rmse.
    """
    tree, taxa_names = generate_backward_yule_tree(n_taxa=8, divergence=1.0, seed=77)
    true_nwk = tree.to_newick() + ";"
    D_true = compute_patristic_distance_matrix(tree, taxa_names)

    # Exact recovery on D_true must yield kf_bsd ~ 0 and patristic_rmse ~ 0
    est_true_nwk = add_nj2(D_true, taxa_names)
    res_true = evaluate_tree_rf(est_true_nwk, true_nwk, taxa_names, D_true=D_true)
    assert res_true["nrf"] == 0.0
    assert res_true["kf_bsd"] < 1e-4
    assert res_true["patristic_rmse"] < 1e-4

    # Simulate with Zipfian Indel (indel_rate = 0.05 and 0.10)
    for irate in [0.05, 0.10]:
        seqs_indel = simulate_wag_gamma_sequences(
            root=tree,
            seq_length=300,
            alpha=1.0,
            seed=101,
            indel_rate=irate,
        )
        assert len(seqs_indel) == 8
        lengths = [len(s) for s in seqs_indel.values()]
        assert all(l > 0 for l in lengths)
        # Because of insertions and deletions, sequence lengths vary across taxa
        assert len(set(lengths)) > 1


def test_exp18_sss_distance_matrices() -> None:
    """
    Verifies Exp 18 SSS-derived distance matrix calculation:
    1. similarity_to_distance_matrix transforms w_ij to linear (1 - w) and log (-ln(w)).
    2. compute_all_pairwise_matrices_and_sss with return_sss_matrices=True returns 5 items.
    """
    from src.sequence_simulator import (
        similarity_to_distance_matrix,
        compute_all_pairwise_matrices_and_sss,
    )

    # 1. Test mathematical transformations of w_ij
    w_test = np.array([
        [1.0, 0.8, 0.0],
        [0.8, 1.0, 0.25],
        [0.0, 0.25, 1.0],
    ], dtype=np.float64)

    d_lin = similarity_to_distance_matrix(w_test, dist_model="linear")
    assert np.allclose(np.diag(d_lin), 0.0)
    assert np.allclose(d_lin, d_lin.T)
    assert np.isclose(d_lin[0, 1], 0.2)
    assert np.isclose(d_lin[0, 2], 1.0)

    d_log = similarity_to_distance_matrix(w_test, dist_model="log", epsilon=1e-4)
    assert np.allclose(np.diag(d_log), 0.0)
    assert np.allclose(d_log, d_log.T)
    assert np.isclose(d_log[0, 1], -math.log(0.8))
    assert np.isclose(d_log[0, 2], -math.log(1e-4))

    # 2. Test return_sss_matrices=True in compute_all_pairwise_matrices_and_sss
    tree, taxa_names = generate_backward_yule_tree(n_taxa=5, divergence=0.5, seed=123)
    seqs = simulate_wag_gamma_sequences(tree, seq_length=200, seed=123)

    D_b, D_w, sss_stats, D_sss_lin, D_sss_log = compute_all_pairwise_matrices_and_sss(
        seqs, taxa_names, return_sss_matrices=True
    )
    assert D_sss_lin.shape == (5, 5)
    assert D_sss_log.shape == (5, 5)
    assert np.allclose(np.diag(D_sss_lin), 0.0)
    assert np.allclose(np.diag(D_sss_log), 0.0)
    assert np.allclose(D_sss_lin, D_sss_lin.T)
    assert np.allclose(D_sss_log, D_sss_log.T)
    assert np.all(D_sss_lin >= 0.0)
    assert np.all(D_sss_log >= 0.0)



