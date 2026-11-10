"""
Evaluation metrics for distance matrix accuracy, four-point condition defect, and tree topology.
"""

from typing import Dict, List, Tuple, Optional
import numpy as np
import dendropy
from dendropy.calculate import treecompare
from src.four_point_core import FourPointCore
from src.tree_inference import parse_unrooted_tree


def evaluate_matrix_accuracy(
    D_est: np.ndarray,
    D_true: np.ndarray,
    eps: float = 1e-7,
) -> Dict[str, float]:
    """
    Computes distance matrix estimation errors against ground truth.

    Parameters
    ----------
    D_est : np.ndarray
        Estimated/refined distance matrix of shape (n, n) in substitutions/site.
    D_true : np.ndarray
        Ground-truth patristic distance matrix of shape (n, n) in substitutions/site.
    eps : float, default=1e-7
        Numerical epsilon.

    Returns
    -------
    Dict[str, float]
        Dictionary of 'rmse', 'mae', 'mre'.
    """
    n = D_est.shape[0]
    triu_indices = np.triu_indices(n, k=1)

    est_vals = D_est[triu_indices]
    true_vals = D_true[triu_indices]

    diff = est_vals - true_vals
    rmse = float(np.sqrt(np.mean(diff**2)))
    mae = float(np.mean(np.abs(diff)))
    mre = float(np.mean(np.abs(diff) / (true_vals + eps)))

    return {
        "rmse": rmse,
        "mae": mae,
        "mre": mre,
    }


def evaluate_four_point_metrics(
    D: np.ndarray,
    core: FourPointCore,
) -> Dict[str, float]:
    """
    Computes four-point condition defects and metric violations.

    Parameters
    ----------
    D : np.ndarray
        Distance matrix of shape (n, n).
    core : FourPointCore
        Preconfigured FourPointCore instance.

    Returns
    -------
    Dict[str, float]
        Dictionary containing 'delta_add' (relative), 'add_abs' (absolute),
        'H_tri' (penalty), and 'n_tri_viol' (count of violating triplets).
    """
    Delta_ADD, _ = core.compute_additivity_defect(D)
    H, _ = core.compute_triangle_penalty(D)

    # Compute absolute defect
    i = core.quartets[:, 0]
    j = core.quartets[:, 1]
    k = core.quartets[:, 2]
    l = core.quartets[:, 3]

    s1 = D[i, j] + D[k, l]
    s2 = D[i, k] + D[j, l]
    s3 = D[i, l] + D[j, k]

    sums = np.stack([s1, s2, s3], axis=1)
    sorted_sums = np.sort(sums, axis=1)[:, ::-1]
    add_abs = float(np.mean(sorted_sums[:, 0] - sorted_sums[:, 1]))

    # Triangle violations count
    ti = core.triplets[:, 0]
    tj = core.triplets[:, 1]
    tk = core.triplets[:, 2]
    dists = np.stack([D[ti, tj], D[ti, tk], D[tj, tk]], axis=1)
    sorted_dists = np.sort(dists, axis=1)[:, ::-1]
    tri_viol_count = int(np.sum((sorted_dists[:, 0] - sorted_dists[:, 1] - sorted_dists[:, 2]) > 1e-7))

    return {
        "delta_add": Delta_ADD,
        "add_abs": add_abs,
        "H_tri": H,
        "delta_um": core.compute_ultrametric_defect(D)[0],
        "n_tri_viol": float(tri_viol_count),
    }


def compute_tree_patristic_matrix(
    tree: dendropy.Tree,
    n_taxa: int,
) -> np.ndarray:
    """
    Computes the pairwise leaf-to-leaf patristic distance matrix from an unrooted
    DendroPy tree with encoded bipartitions using vectorized split-incidence multiplication.

    Parameters
    ----------
    tree : dendropy.Tree
        Unrooted DendroPy tree with `encode_bipartitions()` already called.
    n_taxa : int
        Number of terminal taxa in the tree namespace (unitless).

    Returns
    -------
    np.ndarray
        Symmetric patristic distance matrix of shape (n_taxa, n_taxa) in substitutions/site.
    """
    edges_b: List[float] = []
    edges_m: List[np.ndarray] = []
    all_ones = (1 << n_taxa) - 1

    for edge in tree.postorder_edge_iter():
        if edge.bipartition is None or edge.length is None:
            continue
        b_val = float(edge.length)
        if b_val == 0.0:
            continue
        bm = edge.bipartition.split_bitmask
        if bm == 0 or bm == all_ones:
            continue
        row = np.array([(bm >> i) & 1 for i in range(n_taxa)], dtype=np.float64)
        edges_b.append(b_val)
        edges_m.append(row)

    if not edges_b:
        return np.zeros((n_taxa, n_taxa), dtype=np.float64)

    M = np.array(edges_m, dtype=np.float64)  # (E, N)
    b = np.array(edges_b, dtype=np.float64)  # (E,)
    wM = M * b[:, None]
    one_minus_M = 1.0 - M
    D_pat = wM.T @ one_minus_M + one_minus_M.T @ wM
    np.fill_diagonal(D_pat, 0.0)
    return D_pat


def evaluate_tree_rf(
    est_tree_newick: str,
    true_tree_newick: str,
    taxa_names: List[str],
    D_true: Optional[np.ndarray] = None,
) -> Dict[str, float]:
    """
    Computes unrooted Robinson-Foulds (RF) distance, normalized RF (nRF),
    Kuhner-Felsenstein Branch Score Distance (KF / BSD), and Patristic Distance RMSE.

    Parameters
    ----------
    est_tree_newick : str
        Newick string of inferred tree.
    true_tree_newick : str
        Newick string of ground truth tree.
    taxa_names : List[str]
        List of taxa names (unitless).
    D_true : Optional[np.ndarray], default=None
        Precomputed ground-truth patristic distance matrix of shape (n_taxa, n_taxa)
        in substitutions/site. If None, computed directly from `true_tree_newick`.

    Returns
    -------
    Dict[str, float]
        Dictionary with:
        - 'rf': Robinson-Foulds symmetric difference (unitless).
        - 'nrf': Normalized RF distance in [0, 1] (unitless).
        - 'kf_bsd': Kuhner-Felsenstein Branch Score Distance (substitutions/site).
        - 'patristic_rmse': Root Mean Square Error of tree patristic distances vs true tree (substitutions/site).
    """
    n_taxa = len(taxa_names)
    tns = dendropy.TaxonNamespace(taxa_names)

    t_est = parse_unrooted_tree(est_tree_newick, tns)
    t_true = parse_unrooted_tree(true_tree_newick, tns)

    rf = float(treecompare.symmetric_difference(t_est, t_true))
    max_rf = 2 * (n_taxa - 3) if n_taxa > 3 else 1
    nrf = rf / max_rf if max_rf > 0 else 0.0

    # Kuhner-Felsenstein (1994) Branch Score Distance (BSD)
    kf_bsd = float(treecompare.euclidean_distance(t_est, t_true))

    # Patristic Distance RMSE between inferred tree and true tree
    D_pat_est = compute_tree_patristic_matrix(t_est, n_taxa)
    D_pat_true = D_true if D_true is not None else compute_tree_patristic_matrix(t_true, n_taxa)
    triu_idx = np.triu_indices(n_taxa, k=1)
    diff = D_pat_est[triu_idx] - D_pat_true[triu_idx]
    patristic_rmse = float(np.sqrt(np.mean(diff ** 2))) if len(diff) > 0 else 0.0

    return {
        "rf": rf,
        "nrf": float(nrf),
        "kf_bsd": kf_bsd,
        "patristic_rmse": patristic_rmse,
    }


def compute_topological_distance_edges(
    tree_newick: str,
    taxon1: str,
    taxon2: str,
    taxa_names: List[str],
) -> int:
    """
    Computes topological shortest path distance (number of edges) between two taxa
    in an unrooted phylogenetic tree via Breadth-First Search (BFS).

    Parameters
    ----------
    tree_newick : str
        Newick string of the phylogenetic tree.
    taxon1 : str
        First terminal taxon label.
    taxon2 : str
        Second terminal taxon label.
    taxa_names : List[str]
        List of all valid taxa labels in the tree namespace.

    Returns
    -------
    int
        Shortest path distance in number of edges (branches).
        Returns -1 if either taxon is not found.
    """
    tns = dendropy.TaxonNamespace(taxa_names)
    tree = parse_unrooted_tree(tree_newick, tns)
    tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    node1 = tree.find_node_with_taxon_label(taxon1)
    node2 = tree.find_node_with_taxon_label(taxon2)
    if not node1 or not node2:
        return -1

    visited = {node1: 0}
    queue = [node1]
    while queue:
        cur = queue.pop(0)
        dist = visited[cur]
        if cur == node2:
            return dist
        neighbors = []
        if cur.parent_node:
            neighbors.append(cur.parent_node)
        for child in cur.child_nodes():
            neighbors.append(child)
        for nbr in neighbors:
            if nbr not in visited:
                visited[nbr] = dist + 1
                queue.append(nbr)
    return -1


def evaluate_tree_lba(
    est_tree_newick: str,
    true_tree_newick: str,
    taxa_names: List[str],
    lba_pair: Optional[Tuple[str, str]] = None,
) -> Dict[str, float]:
    """
    Evaluates both general tree topology error (RF/nRF) and specific Long Branch
    Attraction (LBA) metrics between designated long-branch taxa.

    Parameters
    ----------
    est_tree_newick : str
        Newick string of inferred tree.
    true_tree_newick : str
        Newick string of ground-truth tree.
    taxa_names : List[str]
        Complete list of leaf taxa labels.
    lba_pair : Optional[Tuple[str, str]], default=None
        Tuple of the two long-branch taxa. If None, defaults to (taxa_names[0], taxa_names[-1]).

    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'rf': Robinson-Foulds symmetric distance [unitless].
        - 'nrf': Normalized Robinson-Foulds distance in [0, 1] [unitless].
        - 'ds_est': Shortest path edge distance between LBA taxa in estimated tree [edges].
        - 'ds_true': Shortest path edge distance between LBA taxa in true tree [edges].
        - 'delta_ds': Path reduction ds_true - ds_est (positive indicates attraction) [edges].
        - 'is_sister': 1.0 if estimated as immediate sister taxa (ds_est == 2), 0.0 otherwise [unitless].
    """
    rf_res = evaluate_tree_rf(est_tree_newick, true_tree_newick, taxa_names)
    t1 = lba_pair[0] if lba_pair else taxa_names[0]
    t2 = lba_pair[1] if lba_pair else taxa_names[-1]

    ds_est = compute_topological_distance_edges(est_tree_newick, t1, t2, taxa_names)
    ds_true = compute_topological_distance_edges(true_tree_newick, t1, t2, taxa_names)
    delta_ds = (ds_true - ds_est) if (ds_est >= 0 and ds_true >= 0) else 0

    is_sister = 1.0 if ds_est == 2 else 0.0

    return {
        "rf": rf_res["rf"],
        "nrf": rf_res["nrf"],
        "kf_bsd": rf_res["kf_bsd"],
        "patristic_rmse": rf_res["patristic_rmse"],
        "ds_est": float(ds_est),
        "ds_true": float(ds_true),
        "delta_ds": float(delta_ds),
        "is_sister": is_sister,
    }

