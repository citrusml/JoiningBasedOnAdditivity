"""
Unit tests for AddNJ phylogenetic tree reconstruction.
"""

from typing import List
import pytest
import numpy as np
import dendropy
from dendropy.calculate import treecompare
from src.add_nj import add_nj, add_nj2, compute_additivity_defect_fast, QuartetCache


def parse_tree(newick: str, tns: dendropy.TaxonNamespace) -> dendropy.Tree:
    """Helper to parse unrooted tree preserving underscores."""
    t = dendropy.Tree.get(
        data=newick,
        schema="newick",
        taxon_namespace=tns,
        rooting="force-unrooted",
        preserve_underscores=True,
    )
    t.encode_bipartitions()
    return t


def test_quartet_cache() -> None:
    """Verifies QuartetCache generation and shape."""
    cache = QuartetCache(max_taxa=8)
    q5 = cache.get_quartets(5)
    assert q5.shape == (5, 4)  # 5 choose 4 = 5
    q6 = cache.get_quartets(6)
    assert q6.shape == (15, 4)  # 6 choose 4 = 15


def test_additive_matrix_recovery() -> None:
    """
    Verifies that AddNJ perfectly recovers the true tree topology (nRF = 0.0)
    when given an exactly additive distance matrix.
    """
    true_newick = "((T1:0.1,T2:0.2):0.3,T3:0.4,(T4:0.25,T5:0.15):0.35);"
    taxa_names: List[str] = ["T1", "T2", "T3", "T4", "T5"]
    tns = dendropy.TaxonNamespace(taxa_names)

    true_tree = parse_tree(true_newick, tns)
    pdm = true_tree.phylogenetic_distance_matrix()

    n = len(taxa_names)
    D = np.zeros((n, n), dtype=np.float64)
    for i, t_i in enumerate(taxa_names):
        tax_i = tns.get_taxon(t_i)
        for j, t_j in enumerate(taxa_names):
            tax_j = tns.get_taxon(t_j)
            D[i, j] = pdm(tax_i, tax_j)

    # Reconstruct with AddNJ
    inferred_newick = add_nj(D, taxa_names)
    tns_eval = dendropy.TaxonNamespace(taxa_names)
    t_est = parse_tree(inferred_newick, tns_eval)
    t_true = parse_tree(true_newick, tns_eval)

    rf = treecompare.symmetric_difference(t_true, t_est)
    assert rf == 0, f"Expected perfect topology recovery (rf=0), got rf={rf}"


def test_random_additive_tree_recovery() -> None:
    """Tests AddNJ on an 8-taxon randomly generated tree to verify scalability."""
    taxa_names = [f"Taxon_{i}" for i in range(8)]
    tns = dendropy.TaxonNamespace(taxa_names)
    sim_tree = dendropy.simulate.treesim.birth_death_tree(
        birth_rate=1.0, death_rate=0.0, num_extant_tips=8, taxon_namespace=tns
    )
    sim_tree.is_rooted = False
    true_newick = sim_tree.as_string(schema="newick").strip()
    pdm = sim_tree.phylogenetic_distance_matrix()

    n = len(taxa_names)
    D = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            D[i, j] = pdm(tns.get_taxon(taxa_names[i]), tns.get_taxon(taxa_names[j]))

    inferred_newick = add_nj(D, taxa_names)

    tns_eval = dendropy.TaxonNamespace(taxa_names)
    t_est = parse_tree(inferred_newick, tns_eval)
    t_true = parse_tree(true_newick, tns_eval)

    rf = treecompare.symmetric_difference(t_est, t_true)
    assert rf == 0, f"Expected perfect topology recovery on additive tree (rf=0), got {rf}"


def test_add_nj2_additive_recovery() -> None:
    """Verifies that AddNJ2 perfectly recovers true additive tree topologies."""
    true_newick = "((T1:0.1,T2:0.2):0.3,T3:0.4,(T4:0.25,T5:0.15):0.35);"
    taxa_names: List[str] = ["T1", "T2", "T3", "T4", "T5"]
    tns = dendropy.TaxonNamespace(taxa_names)

    true_tree = parse_tree(true_newick, tns)
    pdm = true_tree.phylogenetic_distance_matrix()

    n = len(taxa_names)
    D = np.zeros((n, n), dtype=np.float64)
    for i, t_i in enumerate(taxa_names):
        tax_i = tns.get_taxon(t_i)
        for j, t_j in enumerate(taxa_names):
            tax_j = tns.get_taxon(t_j)
            D[i, j] = pdm(tax_i, tax_j)

    inferred_newick = add_nj2(D, taxa_names)
    tns_eval = dendropy.TaxonNamespace(taxa_names)
    t_est = parse_tree(inferred_newick, tns_eval)
    t_true = parse_tree(true_newick, tns_eval)

    rf = treecompare.symmetric_difference(t_true, t_est)
    assert rf == 0, f"Expected perfect AddNJ2 recovery (rf=0), got rf={rf}"

