"""
Evaluation metrics for distance matrix accuracy, four-point condition defect, and tree topology.
"""

from typing import Dict, List, Tuple, Optional, Set
import itertools
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
    unit_lengths: bool = False,
) -> np.ndarray:
    """
    Computes the pairwise leaf-to-leaf patristic (or topological unit-edge) distance matrix
    from an unrooted DendroPy tree with encoded bipartitions using vectorized split-incidence
    multiplication.

    Parameters
    ----------
    tree : dendropy.Tree
        Unrooted DendroPy tree with `encode_bipartitions()` already called.
    n_taxa : int
        Number of terminal taxa in the tree namespace (unitless).
    unit_lengths : bool, default=False
        If True, treats every unique unrooted edge as length 1.0 to compute exact
        topological edge-count distances (edges).

    Returns
    -------
    np.ndarray
        Symmetric distance matrix of shape (n_taxa, n_taxa) in substitutions/site
        (or edges if `unit_lengths=True`).
    """
    edges_b: List[float] = []
    edges_m: List[np.ndarray] = []
    seen_canon: set = set()
    all_ones = (1 << n_taxa) - 1

    for edge in tree.postorder_edge_iter():
        if edge.bipartition is None:
            continue
        bm = edge.bipartition.split_bitmask
        if bm == 0 or bm == all_ones:
            continue

        if unit_lengths:
            canon = min(bm, all_ones ^ bm)
            if canon in seen_canon:
                continue
            seen_canon.add(canon)
            b_val = 1.0
        else:
            if edge.length is None:
                continue
            b_val = float(edge.length)
            if b_val == 0.0:
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


def find_max_distance_pair(
    tree_newick: str,
    taxa_names: List[str],
    D_true: Optional[np.ndarray] = None,
) -> Tuple[str, str]:
    """
    Identifies the pair of terminal taxa with the maximum patristic distance
    in the reference tree (or precomputed distance matrix).

    Parameters
    ----------
    tree_newick : str
        Newick string of the reference tree (unitless).
    taxa_names : List[str]
        Ordered list of taxon names (unitless).
    D_true : Optional[np.ndarray], default=None
        Precomputed pairwise patristic distance matrix of shape (n_taxa, n_taxa)
        in substitutions/site. If None, computed from `tree_newick`.

    Returns
    -------
    Tuple[str, str]
        Pair of taxon labels (u, v) maximizing D_true[u, v].
    """
    n_taxa = len(taxa_names)
    if n_taxa < 2:
        raise ValueError(f"At least 2 taxa required, got {n_taxa}.")

    if D_true is None:
        tns = dendropy.TaxonNamespace(taxa_names)
        tree = parse_unrooted_tree(tree_newick, tns)
        D_mat = compute_tree_patristic_matrix(tree, n_taxa)
    else:
        D_mat = D_true

    triu_i, triu_j = np.triu_indices(n_taxa, k=1)
    max_flat = int(np.argmax(D_mat[triu_i, triu_j]))
    u_idx = int(triu_i[max_flat])
    v_idx = int(triu_j[max_flat])
    return taxa_names[u_idx], taxa_names[v_idx]


def _extract_pair_separating_splits(
    tree: dendropy.Tree,
    idx_u: int,
    idx_v: int,
    n_taxa: int,
) -> Set[int]:
    """
    Extracts canonical bitmasks of nontrivial internal splits that separate
    taxon `idx_u` and taxon `idx_v` (i.e., internal edges on the path between u and v).
    """
    all_ones = (1 << n_taxa) - 1
    splits: Set[int] = set()
    for edge in tree.postorder_edge_iter():
        if edge.bipartition is None:
            continue
        bm = edge.bipartition.split_bitmask
        if bm == 0 or bm == all_ones:
            continue
        pop = bin(bm).count("1")
        if pop <= 1 or pop >= n_taxa - 1:
            continue
        if ((bm >> idx_u) & 1) != ((bm >> idx_v) & 1):
            splits.add(min(bm, all_ones ^ bm))
    return splits


def compute_pair_separating_rf(
    est_tree: dendropy.Tree,
    true_tree: dendropy.Tree,
    idx_u: int,
    idx_v: int,
    n_taxa: int,
) -> Dict[str, float]:
    """
    Computes Robinson-Foulds metrics restricted exclusively to nontrivial internal
    splits that separate the target pair (u, v).

    Parameters
    ----------
    est_tree : dendropy.Tree
        Unrooted estimated DendroPy tree with encoded bipartitions.
    true_tree : dendropy.Tree
        Unrooted ground-truth DendroPy tree with encoded bipartitions.
    idx_u : int
        Index of first target taxon in the shared TaxonNamespace (unitless).
    idx_v : int
        Index of second target taxon in the shared TaxonNamespace (unitless).
    n_taxa : int
        Total number of terminal taxa (unitless).

    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'pair_nrf': Normalized RF over pair-separating splits in [0, 1] (unitless).
        - 'pair_split_recall': Fraction of true pair-separating splits recovered in [0, 1] (unitless).
        - 'pair_split_fnr': 1 - pair_split_recall in [0, 1] (unitless).
    """
    s_true = _extract_pair_separating_splits(true_tree, idx_u, idx_v, n_taxa)
    s_est = _extract_pair_separating_splits(est_tree, idx_u, idx_v, n_taxa)

    n_true = len(s_true)
    n_est = len(s_est)
    n_shared = len(s_true & s_est)

    recall = float(n_shared) / float(n_true) if n_true > 0 else 1.0
    fnr = 1.0 - recall
    denom = n_true + n_est
    pair_nrf = float(len(s_true ^ s_est)) / float(denom) if denom > 0 else 0.0

    return {
        "pair_nrf": pair_nrf,
        "pair_split_recall": recall,
        "pair_split_fnr": fnr,
    }


def compute_lba_quartet_attraction(
    est_tree: dendropy.Tree,
    true_tree: dendropy.Tree,
    idx_u: int,
    idx_v: int,
    n_taxa: int,
) -> Dict[str, float]:
    """
    Computes noise-invariant LBA quartet metrics over all quartets (u, v, k, l)
    containing the target pair (u, v).

    For quartets where u and v are separated in `true_tree` (topology uk|vl or ul|vk),
    `lba_quartet_attraction` measures the fraction falsely resolved as uv|kl in `est_tree`.

    Parameters
    ----------
    est_tree : dendropy.Tree
        Unrooted estimated DendroPy tree with encoded bipartitions.
    true_tree : dendropy.Tree
        Unrooted ground-truth DendroPy tree with encoded bipartitions.
    idx_u : int
        Index of first target taxon (unitless).
    idx_v : int
        Index of second target taxon (unitless).
    n_taxa : int
        Total number of terminal taxa (unitless).

    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'lba_quartet_attraction': Fraction of separated quartets attracted into uv|kl [unitless].
        - 'lba_quartet_error': Fraction of all (n-2 choose 2) quartets (u,v,k,l) mismatched [unitless].
    """
    others = [x for x in range(n_taxa) if x != idx_u and x != idx_v]
    if len(others) < 2:
        return {"lba_quartet_attraction": 0.0, "lba_quartet_error": 0.0}

    D_topo_true = compute_tree_patristic_matrix(true_tree, n_taxa, unit_lengths=True)
    D_topo_est = compute_tree_patristic_matrix(est_tree, n_taxa, unit_lengths=True)

    pairs_kl = np.array(list(itertools.combinations(others, 2)), dtype=np.int64)
    k = pairs_kl[:, 0]
    l = pairs_kl[:, 1]

    # True quartet topologies for (u, v, k, l)
    st1 = D_topo_true[idx_u, idx_v] + D_topo_true[k, l]
    st2 = D_topo_true[idx_u, k] + D_topo_true[idx_v, l]
    st3 = D_topo_true[idx_u, l] + D_topo_true[idx_v, k]
    St = np.stack([st1, st2, st3], axis=1)
    true_choice = np.argmin(St, axis=1)
    separated_mask = st1 > np.minimum(st2, st3)

    # Estimated quartet topologies for (u, v, k, l)
    se1 = D_topo_est[idx_u, idx_v] + D_topo_est[k, l]
    se2 = D_topo_est[idx_u, k] + D_topo_est[idx_v, l]
    se3 = D_topo_est[idx_u, l] + D_topo_est[idx_v, k]
    Se = np.stack([se1, se2, se3], axis=1)
    est_choice = np.argmin(Se, axis=1)
    est_attracted = se1 < np.minimum(se2, se3)

    n_sep = int(np.sum(separated_mask))
    attraction = float(np.mean(est_attracted[separated_mask])) if n_sep > 0 else 0.0
    q_error = float(np.mean(est_choice != true_choice))

    return {
        "lba_quartet_attraction": attraction,
        "lba_quartet_error": q_error,
    }


def compute_matrix_lba_margin(
    D: np.ndarray,
    idx_u: int,
    idx_v: int,
    true_tree_newick: str,
    taxa_names: List[str],
    eps: float = 1e-7,
) -> float:
    """
    Computes the mean normalized four-point condition LBA margin M_LBA(u, v; D)
    on an observed distance matrix over quartets where (u, v) are separated
    in the true tree.

    Positive values indicate the matrix favors separating (u, v); negative values
    indicate the matrix has entered the Felsenstein attraction zone.

    Parameters
    ----------
    D : np.ndarray
        Symmetric distance matrix of shape (n_taxa, n_taxa) in substitutions/site.
    idx_u : int
        Index of first target taxon (unitless).
    idx_v : int
        Index of second target taxon (unitless).
    true_tree_newick : str
        Newick string of ground-truth tree (unitless).
    taxa_names : List[str]
        Ordered list of taxon labels (unitless).
    eps : float, default=1e-7
        Numerical epsilon to avoid division by zero (unitless).

    Returns
    -------
    float
        Mean normalized LBA quartet margin in [-1, 1] (unitless).
    """
    n_taxa = len(taxa_names)
    others = [x for x in range(n_taxa) if x != idx_u and x != idx_v]
    if len(others) < 2:
        return 0.0

    tns = dendropy.TaxonNamespace(taxa_names)
    true_tree = parse_unrooted_tree(true_tree_newick, tns)
    D_topo_true = compute_tree_patristic_matrix(true_tree, n_taxa, unit_lengths=True)

    pairs_kl = np.array(list(itertools.combinations(others, 2)), dtype=np.int64)
    k = pairs_kl[:, 0]
    l = pairs_kl[:, 1]

    st1 = D_topo_true[idx_u, idx_v] + D_topo_true[k, l]
    st2 = D_topo_true[idx_u, k] + D_topo_true[idx_v, l]
    st3 = D_topo_true[idx_u, l] + D_topo_true[idx_v, k]
    separated_mask = st1 > np.minimum(st2, st3)

    if int(np.sum(separated_mask)) == 0:
        return 0.0

    ks = k[separated_mask]
    ls = l[separated_mask]

    s1 = D[idx_u, idx_v] + D[ks, ls]
    s2 = D[idx_u, ks] + D[idx_v, ls]
    s3 = D[idx_u, ls] + D[idx_v, ks]
    s_max = np.maximum(s1, np.maximum(s2, s3))
    margin = (s1 - np.minimum(s2, s3)) / (s_max + eps)
    return float(np.mean(margin))


def evaluate_tree_lba(
    est_tree_newick: str,
    true_tree_newick: str,
    taxa_names: List[str],
    lba_pair: Optional[Tuple[str, str]] = None,
    auto_max_pair: bool = False,
) -> Dict[str, float]:
    """
    Evaluates both general tree topology error (RF/nRF) and specific Long Branch
    Attraction (LBA) metrics between designated (or max-distance) long-branch taxa.

    Parameters
    ----------
    est_tree_newick : str
        Newick string of inferred tree.
    true_tree_newick : str
        Newick string of ground-truth tree.
    taxa_names : List[str]
        Complete list of leaf taxa labels.
    lba_pair : Optional[Tuple[str, str]], default=None
        Tuple of the two long-branch taxa. If None and `auto_max_pair=True`, selects
        the pair maximizing true patristic distance; otherwise defaults to (taxa_names[0], taxa_names[-1]).
    auto_max_pair : bool, default=False
        If True and `lba_pair` is None, automatically identifies the maximum patristic
        distance pair in `true_tree_newick`.

    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'rf': Robinson-Foulds symmetric distance [unitless].
        - 'nrf': Normalized Robinson-Foulds distance in [0, 1] [unitless].
        - 'kf_bsd': Kuhner-Felsenstein Branch Score Distance [substitutions/site].
        - 'patristic_rmse': Patristic distance RMSE vs true tree [substitutions/site].
        - 'ds_est': Shortest path edge distance between LBA taxa in estimated tree [edges].
        - 'ds_true': Shortest path edge distance between LBA taxa in true tree [edges].
        - 'delta_ds': Path reduction ds_true - ds_est (positive indicates attraction) [edges].
        - 'is_sister': 1.0 if estimated as immediate sister taxa (ds_est == 2), 0.0 otherwise [unitless].
        - 'nta': Normalized Topological Attraction in (-inf, 1] [unitless].
        - 'pair_nrf': Normalized RF restricted to pair-separating internal splits [unitless].
        - 'pair_split_recall': Recall of true pair-separating internal splits [unitless].
        - 'pair_split_fnr': False negative rate of pair-separating internal splits [unitless].
        - 'lba_quartet_attraction': Fraction of separated quartets (u,v,k,l) attracted to uv|kl [unitless].
        - 'lba_quartet_error': Fraction of all (n-2 choose 2) quartets (u,v,k,l) mismatched [unitless].
    """
    rf_res = evaluate_tree_rf(est_tree_newick, true_tree_newick, taxa_names)
    if lba_pair is not None:
        t1, t2 = lba_pair
    elif auto_max_pair:
        t1, t2 = find_max_distance_pair(true_tree_newick, taxa_names)
    else:
        t1, t2 = taxa_names[0], taxa_names[-1]

    ds_est = compute_topological_distance_edges(est_tree_newick, t1, t2, taxa_names)
    ds_true = compute_topological_distance_edges(true_tree_newick, t1, t2, taxa_names)
    delta_ds = (ds_true - ds_est) if (ds_est >= 0 and ds_true >= 0) else 0

    is_sister = 1.0 if ds_est == 2 else 0.0
    nta = float(delta_ds) / float(ds_true - 2) if ds_true > 2 else 0.0

    n_taxa = len(taxa_names)
    name_to_idx = {name: i for i, name in enumerate(taxa_names)}
    idx_u = name_to_idx[t1]
    idx_v = name_to_idx[t2]

    tns = dendropy.TaxonNamespace(taxa_names)
    t_est = parse_unrooted_tree(est_tree_newick, tns)
    t_true = parse_unrooted_tree(true_tree_newick, tns)

    pair_rf_res = compute_pair_separating_rf(t_est, t_true, idx_u, idx_v, n_taxa)
    quartet_res = compute_lba_quartet_attraction(t_est, t_true, idx_u, idx_v, n_taxa)

    return {
        "rf": rf_res["rf"],
        "nrf": rf_res["nrf"],
        "kf_bsd": rf_res["kf_bsd"],
        "patristic_rmse": rf_res["patristic_rmse"],
        "ds_est": float(ds_est),
        "ds_true": float(ds_true),
        "delta_ds": float(delta_ds),
        "is_sister": is_sister,
        "nta": nta,
        "pair_nrf": pair_rf_res["pair_nrf"],
        "pair_split_recall": pair_rf_res["pair_split_recall"],
        "pair_split_fnr": pair_rf_res["pair_split_fnr"],
        "lba_quartet_attraction": quartet_res["lba_quartet_attraction"],
        "lba_quartet_error": quartet_res["lba_quartet_error"],
    }


def compute_root_to_tip_metrics(
    tree_newick: str,
    taxa_names: List[str],
    use_midpoint_rooting: bool = False,
) -> Dict[str, float]:
    """
    Computes root-to-tip distance summary statistics (mean, variance, min, max,
    and max/min ratio) across all terminal taxa of a phylogenetic tree.

    Parameters
    ----------
    tree_newick : str
        Newick string of the tree (unitless).
    taxa_names : List[str]
        List of terminal taxon labels (unitless).
    use_midpoint_rooting : bool, default=False
        If True (recommended for unrooted inferred trees such as NJ/BioNJ/AddNJ),
        places the root at the midpoint of the longest leaf-to-leaf path before
        computing root-to-tip distances.

    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'root_dist_mean': Mean root-to-tip distance [substitutions/site].
        - 'root_dist_var': Population variance of root-to-tip distances [(substitutions/site)^2].
        - 'root_dist_min': Minimum root-to-tip distance [substitutions/site].
        - 'root_dist_max': Maximum root-to-tip distance [substitutions/site].
        - 'root_dist_max_min_ratio': Ratio max / max(min, 1e-12) [unitless].
    """
    tns = dendropy.TaxonNamespace(taxa_names)
    tree = dendropy.Tree.get(
        data=tree_newick,
        schema="newick",
        taxon_namespace=tns,
        rooting="force-unrooted" if use_midpoint_rooting else "force-rooted",
        preserve_underscores=True,
    )
    if use_midpoint_rooting:
        tree.reroot_at_midpoint(update_bipartitions=False)

    dists: List[float] = []
    for leaf in tree.leaf_node_iter():
        d = 0.0
        curr = leaf
        while curr is not None and curr != tree.seed_node:
            if curr.edge_length is not None:
                d += max(0.0, float(curr.edge_length))
            curr = curr.parent_node
        dists.append(d)

    if not dists:
        return {
            "root_dist_mean": 0.0,
            "root_dist_var": 0.0,
            "root_dist_min": 0.0,
            "root_dist_max": 0.0,
            "root_dist_max_min_ratio": 1.0,
        }

    arr = np.array(dists, dtype=np.float64)
    d_min = float(np.min(arr))
    d_max = float(np.max(arr))
    return {
        "root_dist_mean": float(np.mean(arr)),
        "root_dist_var": float(np.var(arr)),
        "root_dist_min": d_min,
        "root_dist_max": d_max,
        "root_dist_max_min_ratio": d_max / max(d_min, 1e-12),
    }

