"""
Unit tests for topology search and TreeScorer.
"""

from typing import List
import pytest
import numpy as np
import dendropy
from dendropy.calculate import treecompare
import Bio.Phylo as Phylo
import io
from src.topology_search import TreeScorer, hill_climb_search


def parse_tree(newick: str, tns: dendropy.TaxonNamespace) -> dendropy.Tree:
    t = dendropy.Tree.get(
        data=newick,
        schema="newick",
        taxon_namespace=tns,
        rooting="force-unrooted",
        preserve_underscores=True,
    )
    t.encode_bipartitions()
    return t


def test_tree_scorer_additive_tree() -> None:
    """Verifies that an additive tree receives penalty = 0.0."""
    true_newick = "((T0:0.1,T1:0.2):0.3,T2:0.4,(T3:0.25,T4:0.15):0.35);"
    taxa = ["T0", "T1", "T2", "T3", "T4"]
    tns = dendropy.TaxonNamespace(taxa)
    tree = parse_tree(true_newick, tns)
    pdm = tree.phylogenetic_distance_matrix()

    n = len(taxa)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = pdm(tns.get_taxon(taxa[i]), tns.get_taxon(taxa[j]))

    scorer = TreeScorer(D, taxa)
    bio_tree = Phylo.read(io.StringIO(true_newick), "newick")
    bio_tree.root_with_outgroup(taxa[0])

    score = scorer.score_tree(bio_tree)
    assert abs(score) < 1e-6, f"Expected perfect score 0.0, got {score}"


def test_hill_climb_improves_topology() -> None:
    """
    Verifies that starting from a suboptimal/perturbed topology,
    hill_climb_search successfully finds the optimal topology with penalty 0.0.
    """
    true_newick = "((T0:0.1,T1:0.2):0.3,T2:0.4,(T3:0.25,T4:0.15):0.35);"
    taxa = ["T0", "T1", "T2", "T3", "T4"]
    tns = dendropy.TaxonNamespace(taxa)
    true_tree = parse_tree(true_newick, tns)
    pdm = true_tree.phylogenetic_distance_matrix()

    n = len(taxa)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = pdm(tns.get_taxon(taxa[i]), tns.get_taxon(taxa[j]))

    scorer = TreeScorer(D, taxa)

    # Suboptimal starting tree (T1 and T3 swapped)
    wrong_newick = "((T0:0.1,T3:0.2):0.3,T2:0.4,(T1:0.25,T4:0.15):0.35);"
    wrong_tree = Phylo.read(io.StringIO(wrong_newick), "newick")
    wrong_tree.root_with_outgroup(taxa[0])
    initial_score = scorer.score_tree(wrong_tree)
    assert initial_score > 0.0, "Suboptimal tree must have positive penalty."

    # Test NNI
    best_nni_newick, best_nni_score = hill_climb_search(wrong_newick, scorer, mode="nni")
    assert best_nni_score < initial_score
    t_opt_nni = parse_tree(best_nni_newick, tns)
    rf_nni = treecompare.symmetric_difference(true_tree, t_opt_nni)
    assert rf_nni == 0, f"NNI should recover true topology (rf=0), got {rf_nni}"

    # Test SPR
    best_spr_newick, best_spr_score = hill_climb_search(wrong_newick, scorer, mode="spr")
    assert best_spr_score < initial_score
    t_opt_spr = parse_tree(best_spr_newick, tns)
    rf_spr = treecompare.symmetric_difference(true_tree, t_opt_spr)
    assert rf_spr == 0, f"SPR should recover true topology (rf=0), got {rf_spr}"

    # Test Hybrid
    best_hyb_newick, best_hyb_score = hill_climb_search(wrong_newick, scorer, mode="hybrid")
    assert best_hyb_score < initial_score
    t_opt_hyb = parse_tree(best_hyb_newick, tns)
    rf_hyb = treecompare.symmetric_difference(true_tree, t_opt_hyb)
    assert rf_hyb == 0, f"Hybrid should recover true topology (rf=0), got {rf_hyb}"
