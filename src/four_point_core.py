"""
Core mathematical functions for four-point condition and triangle inequality optimization.

This module implements the additivity defect loss, triangle inequality penalty,
and analytic subgradients as described in Hyb-Adam-UA (Chaikovskii et al., 2026)
and method.md.
"""

from typing import Tuple, List, Dict, Optional
import itertools
import numpy as np


class FourPointCore:
    """
    Evaluator for Buneman four-point additivity and metric constraints.

    Parameters
    ----------
    n_taxa : int
        Number of taxa (dimension of distance matrix n x n).
    eps : float, default=1e-7
        Small numerical stability constant to prevent division by zero (unitless).
    """

    def __init__(self, n_taxa: int, eps: float = 1e-7) -> None:
        self.n_taxa: int = n_taxa
        self.eps: float = eps

        # Precompute quartet indices (n choose 4)
        quartets: List[Tuple[int, int, int, int]] = list(
            itertools.combinations(range(n_taxa), 4)
        )
        self.quartets: np.ndarray = np.array(quartets, dtype=np.int64)  # Shape: (Q, 4)
        self.n_quartets: int = len(quartets)

        # Precompute triplet indices (n choose 3)
        triplets: List[Tuple[int, int, int]] = list(
            itertools.combinations(range(n_taxa), 3)
        )
        self.triplets: np.ndarray = np.array(triplets, dtype=np.int64)  # Shape: (T, 3)
        self.n_triplets: int = len(triplets)

    def compute_additivity_defect(self, D: np.ndarray) -> Tuple[float, np.ndarray]:
        """
        Calculates relative additivity defect delta_ADD and its mean Delta_ADD.

        Parameters
        ----------
        D : np.ndarray
            Symmetric distance matrix of shape (n_taxa, n_taxa) in substitution units.

        Returns
        -------
        Tuple[float, np.ndarray]
            Delta_ADD : float
                Mean relative additivity defect across all quartets (unitless).
            delta_all : np.ndarray
                Relative additivity defects for each quartet of shape (n_quartets,).
        """
        if self.n_quartets == 0:
            return 0.0, np.zeros(0, dtype=np.float64)

        i = self.quartets[:, 0]
        j = self.quartets[:, 1]
        k = self.quartets[:, 2]
        l = self.quartets[:, 3]

        s1 = D[i, j] + D[k, l]
        s2 = D[i, k] + D[j, l]
        s3 = D[i, l] + D[j, k]

        sums = np.stack([s1, s2, s3], axis=1)  # (Q, 3)
        sorted_sums = np.sort(sums, axis=1)[:, ::-1]  # A >= B >= C

        A = sorted_sums[:, 0]
        B = sorted_sums[:, 1]

        delta = (A - B) / (A + self.eps)
        Delta_ADD = float(np.mean(delta))
        return Delta_ADD, delta

    def compute_triangle_penalty(self, D: np.ndarray) -> Tuple[float, np.ndarray]:
        """
        Calculates metric violation penalty h and its mean H across triplets.

        Parameters
        ----------
        D : np.ndarray
            Symmetric distance matrix of shape (n_taxa, n_taxa) in substitution units.

        Returns
        -------
        Tuple[float, np.ndarray]
            H : float
                Mean triangle penalty across all triplets (unitless).
            h_all : np.ndarray
                Penalties for each triplet of shape (n_triplets,).
        """
        if self.n_triplets == 0:
            return 0.0, np.zeros(0, dtype=np.float64)

        i = self.triplets[:, 0]
        j = self.triplets[:, 1]
        k = self.triplets[:, 2]

        d1 = D[i, j]
        d2 = D[i, k]
        d3 = D[j, k]

        dists = np.stack([d1, d2, d3], axis=1)  # (T, 3)
        sorted_dists = np.sort(dists, axis=1)[:, ::-1]  # a >= b >= c

        a = sorted_dists[:, 0]
        b = sorted_dists[:, 1]
        c = sorted_dists[:, 2]

        diff = a - b - c
        viol = np.maximum(0.0, diff) / (a + self.eps)
        h = viol**2
        H = float(np.mean(h))
        return H, h

    def compute_ultrametric_defect(self, D: np.ndarray) -> Tuple[float, np.ndarray]:
        """
        Calculates relative ultrametricity defect delta_UM(i,j,k) = (a - b) / (a + eps)
        for sorted triplet distances a >= b >= c, and its mean Delta_UM.

        Parameters
        ----------
        D : np.ndarray
            Symmetric distance matrix of shape (n_taxa, n_taxa) in substitution units.

        Returns
        -------
        Tuple[float, np.ndarray]
            Delta_UM : float
                Mean relative ultrametricity defect across all triplets (unitless).
            delta_all : np.ndarray
                Defects for each triplet of shape (n_triplets,).
        """
        if self.n_triplets == 0:
            return 0.0, np.zeros(0, dtype=np.float64)

        i = self.triplets[:, 0]
        j = self.triplets[:, 1]
        k = self.triplets[:, 2]
        dists = np.stack([D[i, j], D[i, k], D[j, k]], axis=1)
        sorted_dists = np.sort(dists, axis=1)[:, ::-1]
        a = sorted_dists[:, 0]
        b = sorted_dists[:, 1]
        delta = (a - b) / (a + self.eps)
        return float(np.mean(delta)), delta

    def compute_loss(
        self,
        D: np.ndarray,
        D_init: Optional[np.ndarray] = None,
        lambda_fid: float = 0.0,
        w_tri: float = 1.0,
    ) -> Dict[str, float]:
        """
        Computes overall scalar loss L(D) and individual component losses.

        Parameters
        ----------
        D : np.ndarray
            Current distance matrix of shape (n_taxa, n_taxa).
        D_init : Optional[np.ndarray], default=None
            Initial noisy distance matrix used for fidelity anchor.
        lambda_fid : float, default=0.0
            Regularization weight for fidelity penalty (unitless).
        w_tri : float, default=1.0
            Weight multiplier for triangle inequality penalty term (unitless).

        Returns
        -------
        Dict[str, float]
            Dictionary containing 'loss', 'delta_add', 'H_tri', and 'loss_fid'.
        """
        Delta_ADD, _ = self.compute_additivity_defect(D)
        H, _ = self.compute_triangle_penalty(D)

        loss_add = 0.5 * (Delta_ADD**2)
        loss_tri = 0.5 * w_tri * H
        loss_fid = 0.0

        if lambda_fid > 0.0 and D_init is not None:
            norm_init_sq = np.sum(D_init**2) + self.eps
            loss_fid = 0.5 * lambda_fid * (np.sum((D - D_init) ** 2) / norm_init_sq)

        total_loss = loss_add + loss_tri + loss_fid

        return {
            "loss": float(total_loss),
            "delta_add": float(Delta_ADD),
            "H_tri": float(H),
            "loss_add": float(loss_add),
            "loss_tri": float(loss_tri),
            "loss_fid": float(loss_fid),
        }

    def compute_gradients(
        self,
        D: np.ndarray,
        D_init: Optional[np.ndarray] = None,
        lambda_fid: float = 0.0,
        w_tri: float = 1.0,
    ) -> np.ndarray:
        """
        Computes exact analytic subgradient matrix dL/dD of shape (n_taxa, n_taxa).

        Gradient is computed with respect to symmetric matrix entries D_ij (i < j).
        The returned gradient matrix G is symmetric with zero diagonal.

        Parameters
        ----------
        D : np.ndarray
            Current symmetric distance matrix of shape (n_taxa, n_taxa).
        D_init : Optional[np.ndarray], default=None
            Initial distance matrix for fidelity loss.
        lambda_fid : float, default=0.0
            Fidelity penalty weight.

        Returns
        -------
        np.ndarray
            Symmetric gradient matrix G of shape (n_taxa, n_taxa).
        """
        n = self.n_taxa
        grad = np.zeros((n, n), dtype=np.float64)

        # 1. Four-point additivity gradient
        if self.n_quartets > 0:
            i = self.quartets[:, 0]
            j = self.quartets[:, 1]
            k = self.quartets[:, 2]
            l = self.quartets[:, 3]

            s1 = D[i, j] + D[k, l]
            s2 = D[i, k] + D[j, l]
            s3 = D[i, l] + D[j, k]

            sums = np.stack([s1, s2, s3], axis=1)  # (Q, 3)
            # argsort descending
            order = np.argsort(-sums, axis=1)  # (Q, 3)

            idx_A = order[:, 0]
            idx_B = order[:, 1]

            A = np.take_along_axis(sums, idx_A[:, None], axis=1).squeeze(1)
            B = np.take_along_axis(sums, idx_B[:, None], axis=1).squeeze(1)

            delta = (A - B) / (A + self.eps)
            Delta_ADD = np.mean(delta)

            # Derivatives d delta / dA and d delta / dB
            d_delta_dA = (B + self.eps) / ((A + self.eps) ** 2)
            d_delta_dB = -1.0 / (A + self.eps)

            # Factor from d(0.5 * Delta_ADD^2) / d delta = Delta_ADD * (1 / |Q|)
            factor = (Delta_ADD / self.n_quartets)
            grad_A = factor * d_delta_dA
            grad_B = factor * d_delta_dB

            # Pair mappings for sum 0, 1, 2
            # s1 (0): (i, j) and (k, l)
            # s2 (1): (i, k) and (j, l)
            # s3 (2): (i, l) and (j, k)
            pairs_map = [
                ((i, j), (k, l)),
                ((i, k), (j, l)),
                ((i, l), (j, k)),
            ]

            for s_idx in range(3):
                mask_A = (idx_A == s_idx)
                mask_B = (idx_B == s_idx)

                # Contribution when s_idx is A
                if np.any(mask_A):
                    g_vals = grad_A[mask_A]
                    (u1, v1), (u2, v2) = pairs_map[s_idx]
                    np.add.at(grad, (u1[mask_A], v1[mask_A]), g_vals)
                    np.add.at(grad, (u2[mask_A], v2[mask_A]), g_vals)

                # Contribution when s_idx is B
                if np.any(mask_B):
                    g_vals = grad_B[mask_B]
                    (u1, v1), (u2, v2) = pairs_map[s_idx]
                    np.add.at(grad, (u1[mask_B], v1[mask_B]), g_vals)
                    np.add.at(grad, (u2[mask_B], v2[mask_B]), g_vals)

        # 2. Triangle inequality penalty gradient
        if self.n_triplets > 0:
            i = self.triplets[:, 0]
            j = self.triplets[:, 1]
            k = self.triplets[:, 2]

            d1 = D[i, j]
            d2 = D[i, k]
            d3 = D[j, k]

            dists = np.stack([d1, d2, d3], axis=1)  # (T, 3)
            order_tri = np.argsort(-dists, axis=1)

            idx_a = order_tri[:, 0]
            idx_b = order_tri[:, 1]
            idx_c = order_tri[:, 2]

            a = np.take_along_axis(dists, idx_a[:, None], axis=1).squeeze(1)
            b = np.take_along_axis(dists, idx_b[:, None], axis=1).squeeze(1)
            c = np.take_along_axis(dists, idx_c[:, None], axis=1).squeeze(1)

            diff = a - b - c
            viol_mask = diff > 0.0

            if np.any(viol_mask):
                a_v = a[viol_mask]
                b_v = b[viol_mask]
                c_v = c[viol_mask]
                diff_v = diff[viol_mask]

                rho = diff_v / (a_v + self.eps)
                # d(0.5 * H) / d h = 0.5 * (1 / |T|)
                # dh/da = 2 * rho * (b + c + eps) / (a + eps)^2
                # dh/db = - 2 * rho / (a + eps)
                # dh/dc = - 2 * rho / (a + eps)
                # Multiplied by 0.5 * w_tri / |T|:
                factor_tri = w_tri / self.n_triplets
                grad_a = factor_tri * rho * (b_v + c_v + self.eps) / ((a_v + self.eps) ** 2)
                grad_b = -factor_tri * rho / (a_v + self.eps)
                grad_c = -factor_tri * rho / (a_v + self.eps)

                tri_pairs = [(i, j), (i, k), (j, k)]
                for d_idx in range(3):
                    u, v = tri_pairs[d_idx]
                    u_v = u[viol_mask]
                    v_v = v[viol_mask]

                    mask_a = (idx_a[viol_mask] == d_idx)
                    mask_b = (idx_b[viol_mask] == d_idx)
                    mask_c = (idx_c[viol_mask] == d_idx)

                    if np.any(mask_a):
                        np.add.at(grad, (u_v[mask_a], v_v[mask_a]), grad_a[mask_a])
                    if np.any(mask_b):
                        np.add.at(grad, (u_v[mask_b], v_v[mask_b]), grad_b[mask_b])
                    if np.any(mask_c):
                        np.add.at(grad, (u_v[mask_c], v_v[mask_c]), grad_c[mask_c])

        # 3. Optional fidelity gradient
        if lambda_fid > 0.0 and D_init is not None:
            norm_init_sq = np.sum(D_init**2) + self.eps
            grad_fid = lambda_fid * (D - D_init) / norm_init_sq
            # Zero out diagonal
            np.fill_diagonal(grad_fid, 0.0)
            grad += grad_fid

        # Symmetrize gradient: since D is symmetric, dL/dD_ij = dL/dD_ij + dL/dD_ji
        sym_grad = grad + grad.T
        np.fill_diagonal(sym_grad, 0.0)
        return sym_grad
