"""
Unit tests for WAG+Gamma pairwise Maximum Likelihood distance estimation
and multi-matrix tree quartet additivity violation scoring.
"""

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

