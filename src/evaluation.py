"""
Evaluation metrics for distance matrix accuracy, four-point condition defect, and tree topology.
"""

from typing import Dict, List, Tuple
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


def evaluate_tree_rf(
    est_tree_newick: str,
    true_tree_newick: str,
    taxa_names: List[str],
) -> Dict[str, float]:
    """
    Computes unrooted Robinson-Foulds (RF) distance and normalized RF (nRF).

    Parameters
    ----------
    est_tree_newick : str
        Newick string of inferred tree.
    true_tree_newick : str
        Newick string of ground truth tree.
    taxa_names : List[str]
        List of taxa names.

    Returns
    -------
    Dict[str, float]
        Dictionary with 'rf' (symmetric difference) and 'nrf' (normalized RF, [0, 1]).
    """
    n_taxa = len(taxa_names)
    tns = dendropy.TaxonNamespace(taxa_names)

    t_est = parse_unrooted_tree(est_tree_newick, tns)
    t_true = parse_unrooted_tree(true_tree_newick, tns)

    rf = float(treecompare.symmetric_difference(t_est, t_true))
    max_rf = 2 * (n_taxa - 3) if n_taxa > 3 else 1
    nrf = rf / max_rf if max_rf > 0 else 0.0

    return {
        "rf": rf,
        "nrf": float(nrf),
    }
