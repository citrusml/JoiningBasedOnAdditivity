"""
Phylogenetic tree inference module.

Reconstructs phylogenetic trees from distance matrices using Neighbor-Joining (NJ)
and provides conversions to Newick format and DendroPy trees.
"""

from typing import List, Optional
import io
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import Bio.Phylo as Phylo
import dendropy


def matrix_to_nj_newick(D: np.ndarray, taxa_names: List[str]) -> str:
    """
    Constructs a Neighbor-Joining (NJ) phylogenetic tree from a distance matrix.

    Parameters
    ----------
    D : np.ndarray
        Distance matrix of shape (n_taxa, n_taxa).
    taxa_names : List[str]
        List of taxa identifiers corresponding to the matrix rows/columns.

    Returns
    -------
    str
        Newick formatted tree string.
    """
    n = len(taxa_names)
    # Build lower triangular matrix format required by Biopython
    lower_matrix: List[List[float]] = []
    for i in range(n):
        row = [float(D[i, j]) for j in range(i + 1)]
        lower_matrix.append(row)

    dm = DistanceMatrix(taxa_names, lower_matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    stream = io.StringIO()
    Phylo.write(tree, stream, "newick")
    newick_str = stream.getvalue().strip()
    return newick_str


def parse_unrooted_tree(
    newick_str: str,
    taxon_namespace: dendropy.TaxonNamespace,
) -> dendropy.Tree:
    """
    Parses a Newick string into a DendroPy Tree with forced unrooting.

    Parameters
    ----------
    newick_str : str
        Newick tree string.
    taxon_namespace : dendropy.TaxonNamespace
        Shared DendroPy taxon namespace.

    Returns
    -------
    dendropy.Tree
        Unrooted DendroPy tree object with encoded bipartitions.
    """
    tree = dendropy.Tree.get(
        data=newick_str,
        schema="newick",
        taxon_namespace=taxon_namespace,
        rooting="force-unrooted",
    )
    tree.encode_bipartitions()
    return tree
