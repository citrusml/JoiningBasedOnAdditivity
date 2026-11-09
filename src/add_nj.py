"""
AddNJ: Additivity-based Neighbor-Joining Phylogenetic Tree Reconstruction.

Reconstructs phylogenetic trees from distance matrices by greedily minimizing
the four-point condition relative additivity defect (Delta_ADD) during
matrix reduction steps.
"""

from typing import List, Tuple, Dict, Optional, Any
import itertools
import numpy as np


class TreeNode:
    """
    Representation of a node in a phylogenetic tree.

    Parameters
    ----------
    name : str
        Taxon name for leaf nodes, or internal node identifier (unitless).
    children : Optional[List[Tuple['TreeNode', float]]], default=None
        List of child node and branch length pairs (branch length in substitutions/site).
    """

    def __init__(
        self,
        name: str,
        children: Optional[List[Tuple["TreeNode", float]]] = None,
    ) -> None:
        self.name: str = name
        self.children: List[Tuple["TreeNode", float]] = children if children is not None else []

    def is_leaf(self) -> bool:
        """
        Checks whether the node is a leaf (terminal taxon).

        Returns
        -------
        bool
            True if leaf, False otherwise (unitless).
        """
        return len(self.children) == 0

    def to_newick(self) -> str:
        """
        Recursively converts the subtree rooted at this node to Newick format.

        Returns
        -------
        str
            Newick formatted representation of the subtree.
        """
        if self.is_leaf():
            return self.name

        sub_strs: List[str] = []
        for child, length in self.children:
            sub_strs.append(f"{child.to_newick()}:{max(0.0, float(length)):.6f}")
        return f"({','.join(sub_strs)})"


class QuartetCache:
    """
    Precomputed quartet combinations for matrix dimensions m = 4 to max_taxa.

    Parameters
    ----------
    max_taxa : int
        Maximum number of taxa (unitless).
    """

    def __init__(self, max_taxa: int) -> None:
        self.cache: Dict[int, np.ndarray] = {}
        for m in range(4, max_taxa + 1):
            combos = list(itertools.combinations(range(m), 4))
            self.cache[m] = np.array(combos, dtype=np.int64)

    def get_quartets(self, m: int) -> np.ndarray:
        """
        Retrieves quartet index array of shape (n_quartets, 4) for dimension m.

        Parameters
        ----------
        m : int
            Dimension of square distance matrix (unitless).

        Returns
        -------
        np.ndarray
            Array of shape (n_quartets, 4) containing taxon indices (unitless).
        """
        if m not in self.cache:
            combos = list(itertools.combinations(range(m), 4))
            self.cache[m] = np.array(combos, dtype=np.int64)
        return self.cache[m]


def compute_additivity_defect_fast(
    D: np.ndarray,
    quartets: np.ndarray,
    eps: float = 1e-7,
) -> float:
    """
    Computes mean relative additivity defect Delta_ADD for a distance matrix.

    Parameters
    ----------
    D : np.ndarray
        Distance matrix of shape (m, m) in substitutions/site.
    quartets : np.ndarray
        Quartet index array of shape (Q, 4) (unitless).
    eps : float, default=1e-7
        Numerical epsilon to avoid division by zero (unitless).

    Returns
    -------
    float
        Mean relative additivity defect Delta_ADD (unitless).
    """
    if len(quartets) == 0:
        return 0.0

    i = quartets[:, 0]
    j = quartets[:, 1]
    k = quartets[:, 2]
    l = quartets[:, 3]

    s1 = D[i, j] + D[k, l]
    s2 = D[i, k] + D[j, l]
    s3 = D[i, l] + D[j, k]

    sums = np.stack([s1, s2, s3], axis=1)
    sorted_sums = np.sort(sums, axis=1)  # Ascending: C <= B <= A
    A = sorted_sums[:, 2]
    B = sorted_sums[:, 1]

    delta = (A - B) / (A + eps)
    return float(np.mean(delta))


def add_nj(
    D: np.ndarray,
    taxa_names: List[str],
    eps: float = 1e-7,
    tie_tolerance: float = 1e-6,
) -> str:
    """
    Constructs a phylogenetic tree Newick string using the AddNJ algorithm.

    Parameters
    ----------
    D : np.ndarray
        Symmetric evolutionary distance matrix of shape (n, n) in substitutions/site.
    taxa_names : List[str]
        Identifiers for each taxon (unitless).
    eps : float, default=1e-7
        Numerical stability constant (unitless).
    tie_tolerance : float, default=1e-6
        Tolerance threshold for tie-breaking using the NJ Q-criterion (unitless).

    Returns
    -------
    str
        Unrooted Newick tree string with branch lengths in substitutions/site.
    """
    n = len(taxa_names)
    if n < 3:
        raise ValueError(f"At least 3 taxa required, got {n}.")

    curr_D = np.copy(D).astype(np.float64)
    nodes: List[TreeNode] = [TreeNode(name=name) for name in taxa_names]
    internal_id = 1
    quartet_cache = QuartetCache(max_taxa=n)

    # Iterative reduction while m > 4
    while len(nodes) > 4:
        m = len(nodes)
        q_sub = quartet_cache.get_quartets(m - 1)
        r = np.sum(curr_D, axis=1)

        best_score = float("inf")
        best_q = float("inf")
        best_pair: Tuple[int, int] = (0, 1)

        # Evaluate all candidate pairs (a, b)
        for a in range(m):
            for b in range(a + 1, m):
                rem_idx = [x for x in range(m) if x != a and x != b]
                d_u = 0.5 * (curr_D[a, rem_idx] + curr_D[b, rem_idx] - curr_D[a, b])
                D_sub = curr_D[np.ix_(rem_idx, rem_idx)]

                D_prime = np.zeros((m - 1, m - 1), dtype=np.float64)
                D_prime[: m - 2, : m - 2] = D_sub
                D_prime[m - 2, : m - 2] = d_u
                D_prime[: m - 2, m - 2] = d_u

                score = compute_additivity_defect_fast(D_prime, q_sub, eps=eps)
                q_val = (m - 2) * curr_D[a, b] - r[a] - r[b]

                if score < best_score - tie_tolerance:
                    best_score = score
                    best_q = q_val
                    best_pair = (a, b)
                elif abs(score - best_score) <= tie_tolerance:
                    if q_val < best_q:
                        best_score = score
                        best_q = q_val
                        best_pair = (a, b)

        a, b = best_pair
        # Calculate NJ branch lengths from nodes a and b to parent u
        r = np.sum(curr_D, axis=1)
        r_a = r[a]
        r_b = r[b]
        d_ab = curr_D[a, b]

        v_a = 0.5 * d_ab + (r_a - r_b) / (2.0 * (m - 2))
        v_b = d_ab - v_a

        # Create new parent node
        node_u = TreeNode(
            name=f"Node_{internal_id}",
            children=[(nodes[a], float(v_a)), (nodes[b], float(v_b))],
        )
        internal_id += 1

        # Form permanent reduced distance matrix and update node list
        rem_idx = [x for x in range(m) if x != a and x != b]
        d_u = 0.5 * (curr_D[a, rem_idx] + curr_D[b, rem_idx] - curr_D[a, b])
        D_sub = curr_D[np.ix_(rem_idx, rem_idx)]

        new_D = np.zeros((m - 1, m - 1), dtype=np.float64)
        new_D[: m - 2, : m - 2] = D_sub
        new_D[m - 2, : m - 2] = d_u
        new_D[: m - 2, m - 2] = d_u

        new_nodes: List[TreeNode] = [nodes[x] for x in rem_idx]
        new_nodes.append(node_u)

        curr_D = new_D
        nodes = new_nodes

    # When m == 4, choose pair minimizing quartet sum
    if len(nodes) == 4:
        s1 = curr_D[0, 1] + curr_D[2, 3]
        s2 = curr_D[0, 2] + curr_D[1, 3]
        s3 = curr_D[0, 3] + curr_D[1, 2]

        if s1 <= s2 and s1 <= s3:
            a, b = 0, 1
        elif s2 <= s1 and s2 <= s3:
            a, b = 0, 2
        else:
            a, b = 0, 3

        r = np.sum(curr_D, axis=1)
        r_a = r[a]
        r_b = r[b]
        d_ab = curr_D[a, b]
        v_a = 0.5 * d_ab + (r_a - r_b) / 4.0
        v_b = d_ab - v_a

        node_u = TreeNode(
            name=f"Node_{internal_id}",
            children=[(nodes[a], float(v_a)), (nodes[b], float(v_b))],
        )
        internal_id += 1

        rem_idx = [x for x in range(4) if x != a and x != b]
        d_u = 0.5 * (curr_D[a, rem_idx] + curr_D[b, rem_idx] - curr_D[a, b])
        D_sub = curr_D[np.ix_(rem_idx, rem_idx)]

        new_D = np.zeros((3, 3), dtype=np.float64)
        new_D[:2, :2] = D_sub
        new_D[2, :2] = d_u
        new_D[:2, 2] = d_u

        new_nodes = [nodes[x] for x in rem_idx]
        new_nodes.append(node_u)

        curr_D = new_D
        nodes = new_nodes

    # When m == 3, resolve central 3-taxon star root
    if len(nodes) == 3:
        d01 = curr_D[0, 1]
        d02 = curr_D[0, 2]
        d12 = curr_D[1, 2]

        v0 = 0.5 * (d01 + d02 - d12)
        v1 = 0.5 * (d01 + d12 - d02)
        v2 = 0.5 * (d02 + d12 - d01)

        root = TreeNode(
            name="Root",
            children=[
                (nodes[0], float(v0)),
                (nodes[1], float(v1)),
                (nodes[2], float(v2)),
            ],
        )
        return root.to_newick() + ";"

    raise RuntimeError(f"Unexpected termination state with {len(nodes)} nodes.")


def add_nj2(
    D: np.ndarray,
    taxa_names: List[str],
    eps: float = 1e-7,
    tie_tolerance: float = 1e-6,
) -> str:
    """
    Constructs a phylogenetic tree Newick string using the AddNJ2 algorithm.

    Selects joining pairs by minimizing Sister-Pair Quartet Inconsistency:
    for candidate pair (a, b), evaluates all other pairs (k, l) to verify
    whether D_ab + D_kl <= min(D_ak + D_bl, D_al + D_bk).

    Parameters
    ----------
    D : np.ndarray
        Symmetric evolutionary distance matrix of shape (n, n) in substitutions/site.
    taxa_names : List[str]
        Identifiers for each taxon (unitless).
    eps : float, default=1e-7
        Numerical stability constant (unitless).
    tie_tolerance : float, default=1e-6
        Tolerance threshold for tie-breaking using the NJ Q-criterion (unitless).

    Returns
    -------
    str
        Unrooted Newick tree string with branch lengths in substitutions/site.
    """
    n = len(taxa_names)
    if n < 3:
        raise ValueError(f"At least 3 taxa required, got {n}.")

    curr_D = np.copy(D).astype(np.float64)
    nodes: List[TreeNode] = [TreeNode(name=name) for name in taxa_names]
    internal_id = 1

    while len(nodes) > 4:
        m = len(nodes)
        r = np.sum(curr_D, axis=1)
        all_indices = list(range(m))

        best_score = float("inf")
        best_q = float("inf")
        best_pair: Tuple[int, int] = (0, 1)

        for a in range(m):
            for b in range(a + 1, m):
                others = [x for x in all_indices if x != a and x != b]
                pairs_other = list(itertools.combinations(others, 2))
                k = np.array([p[0] for p in pairs_other], dtype=np.int64)
                l = np.array([p[1] for p in pairs_other], dtype=np.int64)

                s1 = curr_D[a, b] + curr_D[k, l]
                s2 = curr_D[a, k] + curr_D[b, l]
                s3 = curr_D[a, l] + curr_D[b, k]

                s_min = np.minimum(s1, np.minimum(s2, s3))
                s_max = np.maximum(s1, np.maximum(s2, s3))

                penalty = (s1 - s_min) / (s_max + eps)
                score = float(np.mean(penalty))
                q_val = (m - 2) * curr_D[a, b] - r[a] - r[b]

                if score < best_score - tie_tolerance:
                    best_score = score
                    best_q = q_val
                    best_pair = (a, b)
                elif abs(score - best_score) <= tie_tolerance:
                    if q_val < best_q:
                        best_score = score
                        best_q = q_val
                        best_pair = (a, b)

        a, b = best_pair
        r = np.sum(curr_D, axis=1)
        r_a = r[a]
        r_b = r[b]
        d_ab = curr_D[a, b]

        v_a = 0.5 * d_ab + (r_a - r_b) / (2.0 * (m - 2))
        v_b = d_ab - v_a

        node_u = TreeNode(
            name=f"Node_{internal_id}",
            children=[(nodes[a], float(v_a)), (nodes[b], float(v_b))],
        )
        internal_id += 1

        rem_idx = [x for x in range(m) if x != a and x != b]
        d_u = 0.5 * (curr_D[a, rem_idx] + curr_D[b, rem_idx] - curr_D[a, b])
        D_sub = curr_D[np.ix_(rem_idx, rem_idx)]

        new_D = np.zeros((m - 1, m - 1), dtype=np.float64)
        new_D[: m - 2, : m - 2] = D_sub
        new_D[m - 2, : m - 2] = d_u
        new_D[: m - 2, m - 2] = d_u

        new_nodes: List[TreeNode] = [nodes[x] for x in rem_idx]
        new_nodes.append(node_u)

        curr_D = new_D
        nodes = new_nodes

    if len(nodes) == 4:
        s1 = curr_D[0, 1] + curr_D[2, 3]
        s2 = curr_D[0, 2] + curr_D[1, 3]
        s3 = curr_D[0, 3] + curr_D[1, 2]

        if s1 <= s2 and s1 <= s3:
            a, b = 0, 1
        elif s2 <= s1 and s2 <= s3:
            a, b = 0, 2
        else:
            a, b = 0, 3

        r = np.sum(curr_D, axis=1)
        r_a = r[a]
        r_b = r[b]
        d_ab = curr_D[a, b]
        v_a = 0.5 * d_ab + (r_a - r_b) / 4.0
        v_b = d_ab - v_a

        node_u = TreeNode(
            name=f"Node_{internal_id}",
            children=[(nodes[a], float(v_a)), (nodes[b], float(v_b))],
        )
        internal_id += 1

        rem_idx = [x for x in range(4) if x != a and x != b]
        d_u = 0.5 * (curr_D[a, rem_idx] + curr_D[b, rem_idx] - curr_D[a, b])
        D_sub = curr_D[np.ix_(rem_idx, rem_idx)]

        new_D = np.zeros((3, 3), dtype=np.float64)
        new_D[:2, :2] = D_sub
        new_D[2, :2] = d_u
        new_D[:2, 2] = d_u

        new_nodes = [nodes[x] for x in rem_idx]
        new_nodes.append(node_u)

        curr_D = new_D
        nodes = new_nodes

    if len(nodes) == 3:
        d01 = curr_D[0, 1]
        d02 = curr_D[0, 2]
        d12 = curr_D[1, 2]

        v0 = 0.5 * (d01 + d02 - d12)
        v1 = 0.5 * (d01 + d12 - d02)
        v2 = 0.5 * (d02 + d12 - d01)

        root = TreeNode(
            name="Root",
            children=[
                (nodes[0], float(v0)),
                (nodes[1], float(v1)),
                (nodes[2], float(v2)),
            ],
        )
        return root.to_newick() + ";"

    raise RuntimeError(f"Unexpected termination state with {len(nodes)} nodes.")

