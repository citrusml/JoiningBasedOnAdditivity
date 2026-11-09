"""
Unit tests for FastME wrapper module.
"""

from typing import List
import pytest
import numpy as np
import dendropy
from dendropy.calculate import treecompare
from src.fastme_wrapper import run_fastme


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


def test_fastme_bionj_additive() -> None:
    """Verifies that FastME BioNJ reconstructs an additive tree perfectly."""
    true_newick = "((T1:0.1,T2:0.2):0.3,T3:0.4,(T4:0.25,T5:0.15):0.35);"
    taxa = ["T1", "T2", "T3", "T4", "T5"]
    tns = dendropy.TaxonNamespace(taxa)
    true_tree = parse_tree(true_newick, tns)
    pdm = true_tree.phylogenetic_distance_matrix()

    n = len(taxa)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = pdm(tns.get_taxon(taxa[i]), tns.get_taxon(taxa[j]))

    for m in ["bionj", "nni", "spr"]:
        nwk = run_fastme(D, taxa, method=m)
        t_est = parse_tree(nwk, tns)
        rf = treecompare.symmetric_difference(true_tree, t_est)
        assert rf == 0, f"FastME method {m} failed to recover additive tree (rf={rf})"

