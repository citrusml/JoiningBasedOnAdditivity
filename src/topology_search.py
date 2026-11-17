"""
Topology Search Module for Additivity-based Phylogenetic Inference.

Performs Nearest Neighbor Interchange (NNI) and Subtree Pruning and Regrafting (SPR)
hill-climbing heuristic search using the tree-induced quartet additivity penalty.
"""

from typing import List, Tuple, Dict, Optional, Set, Union
import io
import math
import copy
import itertools
import numpy as np
import Bio.Phylo as Phylo
from Bio.Phylo.BaseTree import Tree, Clade
from Bio.Phylo.TreeConstruction import NNITreeSearcher, ParsimonyScorer


class TreeScorer:
    """
    Evaluator for the tree-induced quartet additivity penalty.

    Parameters
    ----------
    D_obs : np.ndarray
        Observed symmetric evolutionary distance matrix of shape (n, n) in substitutions/site.
    taxa_names : List[str]
        List of taxon labels corresponding to matrix rows/columns (unitless).
    eps : float, default=1e-7
        Numerical epsilon to avoid zero division (unitless).
    max_quartets : Optional[int], default=100000
        Maximum number of quartets to evaluate when (n choose 4) is very large (e.g., N=128).
        If None, evaluates all (n choose 4) quartets exhaustively.
    """

    def __init__(
        self,
        D_obs: np.ndarray,
        taxa_names: List[str],
        eps: float = 1e-7,
        max_quartets: Optional[int] = 100000,
    ) -> None:
        self.D_obs: np.ndarray = D_obs
        self.taxa_names: List[str] = taxa_names
        self.n: int = len(taxa_names)
        self.eps: float = eps
        self.taxa_to_idx: Dict[str, int] = {name: i for i, name in enumerate(taxa_names)}

        total_q = math.comb(self.n, 4) if self.n >= 4 else 0
        if total_q == 0:
            self.quartets = np.zeros((0, 4), dtype=np.int64)
            self.n_quartets = 0
        elif max_quartets is not None and total_q > max_quartets:
            step = max(1, total_q // max_quartets)
            combos = list(
                itertools.islice(
                    itertools.combinations(range(self.n), 4),
                    0,
                    total_q,
                    step,
                )
            )[:max_quartets]
            self.quartets = np.array(combos, dtype=np.int64)
            self.n_quartets = len(combos)
        else:
            combos = list(itertools.combinations(range(self.n), 4))
            self.quartets = np.array(combos, dtype=np.int64)  # (Q, 4)
            self.n_quartets = len(combos)

        # Precompute observed quartet pairwise sums
        if self.n_quartets > 0:
            i = self.quartets[:, 0]
            j = self.quartets[:, 1]
            k = self.quartets[:, 2]
            l = self.quartets[:, 3]

            s0 = D_obs[i, j] + D_obs[k, l]
            s1 = D_obs[i, k] + D_obs[j, l]
            s2 = D_obs[i, l] + D_obs[j, k]

            self.S_obs: np.ndarray = np.stack([s0, s1, s2], axis=1)  # (Q, 3)
            self.S_min: np.ndarray = np.min(self.S_obs, axis=1)
            self.S_max: np.ndarray = np.max(self.S_obs, axis=1)
        else:
            self.S_obs = np.zeros((0, 3), dtype=np.float64)
            self.S_min = np.zeros(0, dtype=np.float64)
            self.S_max = np.zeros(0, dtype=np.float64)


    def score_tree(self, tree: Tree) -> float:
        """
        Computes the mean tree-induced quartet additivity penalty L(T).

        Parameters
        ----------
        tree : Tree
            Biopython Tree object to evaluate (unitless).

        Returns
        -------
        float
            Mean quartet additivity penalty in [0, 1] (unitless).
        """
        if self.n_quartets == 0:
            return 0.0

        # Compute topological distance (edge counts) between leaves
        terminals = [t.name for t in tree.get_terminals()]
        # Map each taxon to its terminal clade
        name_to_clade = {t.name: t for t in tree.get_terminals()}

        # Build adjacency graph for leaf-to-leaf topological distances
        adj: Dict[int, List[int]] = {}
        all_clades = list(tree.find_clades())
        clade_to_id = {c: i for i, c in enumerate(all_clades)}

        for c in all_clades:
            c_id = clade_to_id[c]
            adj.setdefault(c_id, [])
            for child in c.clades:
                ch_id = clade_to_id[child]
                adj[c_id].append(ch_id)
                adj.setdefault(ch_id, []).append(c_id)

        # BFS distances between all leaf pairs
        D_topo = np.zeros((self.n, self.n), dtype=np.float64)
        for idx, t_name in enumerate(self.taxa_names):
            if t_name not in name_to_clade:
                continue
            start_id = clade_to_id[name_to_clade[t_name]]
            dist: Dict[int, int] = {start_id: 0}
            queue = [start_id]
            while queue:
                curr = queue.pop(0)
                d = dist[curr]
                for nxt in adj[curr]:
                    if nxt not in dist:
                        dist[nxt] = d + 1
                        queue.append(nxt)

            for jdx, target_name in enumerate(self.taxa_names):
                if target_name in name_to_clade:
                    t_id = clade_to_id[name_to_clade[target_name]]
                    D_topo[idx, jdx] = dist.get(t_id, 0)

        # Determine tree's induced topology on each quartet via topological distance
        i = self.quartets[:, 0]
        j = self.quartets[:, 1]
        k = self.quartets[:, 2]
        l = self.quartets[:, 3]

        st0 = D_topo[i, j] + D_topo[k, l]
        st1 = D_topo[i, k] + D_topo[j, l]
        st2 = D_topo[i, l] + D_topo[j, k]
        St = np.stack([st0, st1, st2], axis=1)

        # Tree quartet partition choice (smallest sum)
        choice = np.argmin(St, axis=1)

        chosen_s = self.S_obs[np.arange(self.n_quartets), choice]
        penalty = (chosen_s - self.S_min) / (self.S_max + self.eps)
        return float(np.mean(penalty))

    def score_newick(self, newick_str: str) -> float:

        """
        Parses a Newick string and computes its mean quartet additivity penalty L(T).

        Parameters
        ----------
        newick_str : str
            Newick string representing the phylogenetic tree (unitless).

        Returns
        -------
        float
            Mean quartet additivity penalty in [0, 1] (unitless).
        """
        tree: Tree = Phylo.read(io.StringIO(newick_str), "newick")
        return self.score_tree(tree)

    def extract_tree_quartet_choices(self, tree_or_newick: Union[Tree, str]) -> np.ndarray:
        """
        Computes the quartet partition choices (0, 1, or 2) induced by a tree topology.

        Parameters
        ----------
        tree_or_newick : Union[Tree, str]
            Biopython Tree object or Newick string (unitless).

        Returns
        -------
        np.ndarray
            1D integer array of shape (n_quartets,) with values in {0, 1, 2}.
        """
        if self.n_quartets == 0:
            return np.zeros(0, dtype=np.int64)

        if isinstance(tree_or_newick, str):
            tree: Tree = Phylo.read(io.StringIO(tree_or_newick), "newick")
        else:
            tree = tree_or_newick

        name_to_clade = {t.name: t for t in tree.get_terminals()}
        adj: Dict[int, List[int]] = {}
        all_clades = list(tree.find_clades())
        clade_to_id = {c: i for i, c in enumerate(all_clades)}

        for c in all_clades:
            c_id = clade_to_id[c]
            adj.setdefault(c_id, [])
            for child in c.clades:
                ch_id = clade_to_id[child]
                adj[c_id].append(ch_id)
                adj.setdefault(ch_id, []).append(c_id)

        D_topo = np.zeros((self.n, self.n), dtype=np.float64)
        for idx, t_name in enumerate(self.taxa_names):
            if t_name not in name_to_clade:
                continue
            start_id = clade_to_id[name_to_clade[t_name]]
            dist: Dict[int, int] = {start_id: 0}
            queue = [start_id]
            while queue:
                curr = queue.pop(0)
                d = dist[curr]
                for nxt in adj[curr]:
                    if nxt not in dist:
                        dist[nxt] = d + 1
                        queue.append(nxt)

            for jdx, target_name in enumerate(self.taxa_names):
                if target_name in name_to_clade:
                    t_id = clade_to_id[name_to_clade[target_name]]
                    D_topo[idx, jdx] = dist.get(t_id, 0)

        i = self.quartets[:, 0]
        j = self.quartets[:, 1]
        k = self.quartets[:, 2]
        l = self.quartets[:, 3]

        st0 = D_topo[i, j] + D_topo[k, l]
        st1 = D_topo[i, k] + D_topo[j, l]
        st2 = D_topo[i, l] + D_topo[j, k]
        St = np.stack([st0, st1, st2], axis=1)
        return np.argmin(St, axis=1)

    def score_choices(self, choice: np.ndarray) -> float:
        """
        Scores pre-extracted tree quartet choices against this scorer's distance matrix.

        Parameters
        ----------
        choice : np.ndarray
            1D integer array of shape (n_quartets,) with values in {0, 1, 2}.

        Returns
        -------
        float
            Mean quartet additivity penalty in [0, 1] (unitless).
        """
        if self.n_quartets == 0 or len(choice) == 0:
            return 0.0
        chosen_s = self.S_obs[np.arange(self.n_quartets), choice]
        penalty = (chosen_s - self.S_min) / (self.S_max + self.eps)
        return float(np.mean(penalty))

    def compute_matrix_defect(self) -> float:
        """
        Computes the intrinsic 4-point condition additivity defect Delta_ADD(D)
        of the distance matrix itself over the evaluated quartets.

        Returns
        -------
        float
            Mean relative additivity defect (S_max - S_mid) / (S_max + eps) in [0, 1] (unitless).
        """
        if self.n_quartets == 0:
            return 0.0
        s_mid = np.sum(self.S_obs, axis=1) - self.S_max - self.S_min
        defect = np.maximum(0.0, self.S_max - s_mid) / (self.S_max + self.eps)
        return float(np.mean(defect))


def _build_split_incidence_matrix(
    tree: Tree,
    name_to_idx: Dict[str, int],
    triu_i: np.ndarray,
    triu_j: np.ndarray,
    reset_root_and_empty: bool = False,
) -> Tuple[np.ndarray, List[Tuple[bool, ...]], Dict[Tuple[bool, ...], List[Clade]]]:
    """
    Constructs the split-path incidence matrix A_T for all nontrivial edges of a tree.

    Parameters
    ----------
    tree : Tree
        Biopython Tree object (unitless).
    name_to_idx : Dict[str, int]
        Mapping from terminal taxon names to row/column indices (unitless).
    triu_i : np.ndarray
        Upper-triangular row indices of shape (n*(n-1)//2,) (unitless).
    triu_j : np.ndarray
        Upper-triangular column indices of shape (n*(n-1)//2,) (unitless).
    reset_root_and_empty : bool, default=False
        If True, resets root and empty/full clade branch lengths to 0.0 in-place.

    Returns
    -------
    Tuple[np.ndarray, List[Tuple[bool, ...]], Dict[Tuple[bool, ...], List[Clade]]]
        A : np.ndarray
            Incidence matrix of shape (n*(n-1)//2, n_unique_splits) (unitless).
        canons : List[Tuple[bool, ...]]
            Ordered list of canonical split keys corresponding to columns of A (unitless).
        seen_splits : Dict[Tuple[bool, ...], List[Clade]]
            Mapping from canonical split key to corresponding Biopython Clade objects.
    """
    n = len(name_to_idx)
    canons: List[Tuple[bool, ...]] = []
    masks: List[np.ndarray] = []
    seen_splits: Dict[Tuple[bool, ...], List[Clade]] = {}

    if reset_root_and_empty:
        tree.root.branch_length = 0.0

    # Postorder memoization avoids repeated get_terminals() traversals across internal clades
    clade_masks: Dict[Clade, np.ndarray] = {}
    for c in tree.find_clades(order="postorder"):
        if c.is_terminal():
            mask = np.zeros(n, dtype=bool)
            if c.name in name_to_idx:
                mask[name_to_idx[c.name]] = True
            clade_masks[c] = mask
        else:
            mask = np.zeros(n, dtype=bool)
            for ch in c.clades:
                mask |= clade_masks[ch]
            clade_masks[c] = mask

        if c == tree.root:
            continue
        n_in = int(np.sum(mask))
        if n_in == 0 or n_in == n:
            if reset_root_and_empty:
                c.branch_length = 0.0
            continue
        canon = tuple(~mask) if mask[0] else tuple(mask)
        if canon in seen_splits:
            seen_splits[canon].append(c)
        else:
            seen_splits[canon] = [c]
            canons.append(canon)
            masks.append(mask)

    if not masks:
        return np.zeros((len(triu_i), 0), dtype=np.float64), canons, seen_splits

    A = np.column_stack([(m[triu_i] ^ m[triu_j]).astype(np.float64) for m in masks])
    return A, canons, seen_splits


class TreeDistanceScorer:
    """
    Evaluator for the discrepancy between tree-induced patristic distances D_T
    (fitted via least squares on the candidate topology) and observed matrix
    distances D_obs.

    Supports objective modes:
    - 'ols_rmse': Unweighted OLS fit with RMSE discrepancy [substitutions/site].
    - 'ols_mae': Unweighted OLS fit with MAE discrepancy [substitutions/site].
    - 'wls_fm': Fitch-Margoliash WLS (w = 1 / (d + eps)) with relative squared error [unitless].
    - 'wls_mre': WLS fit with Mean Relative Error [unitless].

    Parameters
    ----------
    D_obs : np.ndarray
        Observed symmetric evolutionary distance matrix of shape (n, n) in substitutions/site.
    taxa_names : List[str]
        List of taxon labels corresponding to matrix rows/columns (unitless).
    mode : str, default='ols_rmse'
        Distance discrepancy objective mode (unitless).
    eps : float, default=1e-7
        Numerical epsilon to avoid division by zero (unitless).
    reg : float, default=1e-12
        Ridge regularization constant for solving normal equations (unitless).
    """

    SUPPORTED_MODES = ("ols_rmse", "ols_mae", "wls_fm", "wls_mre")

    def __init__(
        self,
        D_obs: np.ndarray,
        taxa_names: List[str],
        mode: str = "ols_rmse",
        eps: float = 1e-7,
        reg: float = 1e-12,
    ) -> None:
        if mode not in self.SUPPORTED_MODES:
            raise ValueError(f"Unsupported mode '{mode}'. Choose from {self.SUPPORTED_MODES}.")
        self.D_obs: np.ndarray = D_obs
        self.taxa_names: List[str] = taxa_names
        self.n: int = len(taxa_names)
        self.mode: str = mode
        self.eps: float = eps
        self.reg: float = reg
        self.name_to_idx: Dict[str, int] = {name: i for i, name in enumerate(taxa_names)}

        self.triu_i, self.triu_j = np.triu_indices(self.n, k=1)
        self.d_vec: np.ndarray = D_obs[self.triu_i, self.triu_j].astype(np.float64)
        self.weights: Optional[np.ndarray] = (
            1.0 / (self.d_vec + self.eps) if mode.startswith("wls") else None
        )

    def score_tree(self, tree: Tree) -> float:
        """
        Computes the distance discrepancy between tree-induced patristic distances
        and the observed distance matrix D_obs.

        Parameters
        ----------
        tree : Tree
            Biopython Tree object to evaluate (unitless).

        Returns
        -------
        float
            Discrepancy score (RMSE/MAE in substitutions/site or dimensionless WLS error).
        """
        if len(self.d_vec) == 0:
            return 0.0

        A, _, _ = _build_split_incidence_matrix(
            tree, self.name_to_idx, self.triu_i, self.triu_j, reset_root_and_empty=False
        )
        if A.shape[1] == 0:
            return float(np.sqrt(np.mean(self.d_vec ** 2)))

        if self.weights is not None:
            A_fit = A * self.weights[:, None]
            d_fit = self.d_vec * self.weights
        else:
            A_fit = A
            d_fit = self.d_vec

        ata = A_fit.T @ A_fit
        ata[np.diag_indices_from(ata)] += self.reg
        atd = A_fit.T @ d_fit
        try:
            b = np.linalg.solve(ata, atd)
        except np.linalg.LinAlgError:
            b, _, _, _ = np.linalg.lstsq(A_fit, d_fit, rcond=None)

        b = np.maximum(0.0, b)
        d_tree = A @ b
        diff = d_tree - self.d_vec

        if self.mode == "ols_rmse":
            return float(np.sqrt(np.mean(diff ** 2)))
        if self.mode == "ols_mae":
            return float(np.mean(np.abs(diff)))
        if self.mode == "wls_fm":
            return float(np.mean((diff / (self.d_vec + self.eps)) ** 2))
        return float(np.mean(np.abs(diff) / (self.d_vec + self.eps)))

    def score_newick(self, newick_str: str) -> float:
        """
        Parses a Newick string and computes its tree-vs-matrix distance discrepancy.

        Parameters
        ----------
        newick_str : str
            Newick string representing the phylogenetic tree (unitless).

        Returns
        -------
        float
            Discrepancy score (substitutions/site or unitless).
        """
        tree: Tree = Phylo.read(io.StringIO(newick_str), "newick")
        return self.score_tree(tree)


def get_nni_neighbors(tree: Tree) -> List[Tree]:
    """
    Generates all Nearest Neighbor Interchange (NNI) neighbors for a binary tree.

    Parameters
    ----------
    tree : Tree
        Biopython Tree object (unitless).

    Returns
    -------
    List[Tree]
        List of neighboring Tree objects generated by NNI (unitless).
    """
    class _DummyScorer(ParsimonyScorer):
        def __init__(self) -> None:
            pass

    searcher = NNITreeSearcher(_DummyScorer())
    return searcher._get_neighbors(tree)


def get_spr_neighbors(tree: Tree) -> List[Tree]:
    """
    Generates all Subtree Pruning and Regrafting (SPR) neighbors for a binary tree.

    Parameters
    ----------
    tree : Tree
        Biopython Tree object (unitless).

    Returns
    -------
    List[Tree]
        List of neighboring Tree objects generated by SPR (unitless).
    """
    neighbors: List[Tree] = []
    parents: Dict[Clade, Clade] = {}
    for c in tree.find_clades():
        for ch in c.clades:
            parents[ch] = c

    all_clades = list(tree.find_clades())
    for c in all_clades:
        if c == tree.root:
            continue
        p = parents[c]
        if len(p.clades) != 2:
            continue
        s = p.clades[1] if p.clades[0] == c else p.clades[0]
        g = parents.get(p)

        sub_c: Set[Clade] = set(c.find_clades())

        for target in all_clades:
            if target in sub_c or target == p or target == s:
                continue
            if g is not None and target == g and len(g.clades) == 2:
                continue

            # Clone and apply SPR
            tree_copy = copy.deepcopy(tree)
            parents_copy: Dict[Clade, Clade] = {}
            for node in tree_copy.find_clades():
                for ch in node.clades:
                    parents_copy[ch] = node

            # Locate c in copy via index path
            curr = tree_copy.root
            for step in tree.get_path(c):
                orig_parent = parents[step]
                idx = orig_parent.clades.index(step)
                curr = curr.clades[idx]
            c_cp = curr
            p_cp = parents_copy[c_cp]
            s_cp = p_cp.clades[1] if p_cp.clades[0] == c_cp else p_cp.clades[0]
            g_cp = parents_copy.get(p_cp)

            # Locate target in copy
            curr_t = tree_copy.root
            if target != tree.root:
                for step in tree.get_path(target):
                    orig_parent = parents[step]
                    idx = orig_parent.clades.index(step)
                    curr_t = curr_t.clades[idx]
            target_cp = curr_t
            p_target_cp = parents_copy.get(target_cp)

            # 1. Prune
            if g_cp is None:
                tree_copy.root = s_cp
            else:
                idx = g_cp.clades.index(p_cp)
                g_cp.clades[idx] = s_cp

            # 2. Regraft
            p_cp.clades = [target_cp, c_cp]
            if p_target_cp is None or target_cp == tree_copy.root:
                tree_copy.root = p_cp
            else:
                idx = p_target_cp.clades.index(target_cp)
                p_target_cp.clades[idx] = p_cp

            neighbors.append(tree_copy)

    return neighbors


def fit_ols_branch_lengths(
    tree: Tree,
    D_obs: np.ndarray,
    taxa_names: List[str],
) -> None:
    """
    Fits Ordinary Least Squares (OLS) branch lengths (clamped to >= 0.0) on a fixed
    Biopython Tree topology against the observed distance matrix D_obs.

    Parameters
    ----------
    tree : Tree
        Biopython Tree object whose clade branch lengths are updated in-place.
    D_obs : np.ndarray
        Observed pairwise distance matrix of shape (n, n) in substitutions/site.
    taxa_names : List[str]
        Ordered list of taxon names corresponding to D_obs rows/columns (unitless).
    """
    n = len(taxa_names)
    if n < 3:
        return

    name_to_idx = {name: i for i, name in enumerate(taxa_names)}
    triu_i, triu_j = np.triu_indices(n, k=1)
    d_vec = D_obs[triu_i, triu_j]

    A, canons, seen_splits = _build_split_incidence_matrix(
        tree, name_to_idx, triu_i, triu_j, reset_root_and_empty=True
    )
    if A.shape[1] == 0:
        return

    b_ols, _, _, _ = np.linalg.lstsq(A, d_vec, rcond=None)
    b_ols = np.maximum(0.0, b_ols)

    for canon, val in zip(canons, b_ols):
        clades = seen_splits[canon]
        for cl in clades:
            cl.branch_length = float(val) / len(clades)


def tree_to_newick(tree: Tree) -> str:
    """Converts a Biopython Tree to a Newick string."""
    stream = io.StringIO()
    Phylo.write(tree, stream, "newick")
    return stream.getvalue().strip()


def hill_climb_search(
    initial_newick: str,
    scorer: Union[TreeScorer, TreeDistanceScorer],
    mode: str = "hybrid",
    max_iter: int = 50,
    verbose: bool = False,
) -> Tuple[str, float]:
    """
    Executes hill-climbing topology search to minimize the scorer objective
    (quartet additivity penalty or tree-vs-matrix distance discrepancy).

    Parameters
    ----------
    initial_newick : str
        Starting tree in Newick format (unitless).
    scorer : Union[TreeScorer, TreeDistanceScorer]
        Scorer instance configured with observed distance matrix (unitless).
    mode : str, default='hybrid'
        Search strategy: 'nni', 'spr', or 'hybrid' (NNI followed by SPR) (unitless).
    max_iter : int, default=50
        Maximum iterations per search phase (unitless).
    verbose : bool, default=False
        Whether to print search progression (unitless).

    Returns
    -------
    Tuple[str, float]
        best_newick : str
            Optimized tree in Newick format.
        best_score : float
            Final objective score (penalty or distance discrepancy).
    """
    tree = Phylo.read(io.StringIO(initial_newick), "newick")
    outgroup = scorer.taxa_names[0]
    tree.root_with_outgroup(outgroup)

    current_tree = tree
    best_score = scorer.score_tree(current_tree)
    topology_changed = False
    tol = 1e-9 if isinstance(scorer, TreeDistanceScorer) else 1e-7

    def _step(current_t: Tree, current_s: float, search_type: str) -> Tuple[Tree, float, bool]:
        get_nb = get_nni_neighbors if search_type == "nni" else get_spr_neighbors
        neighbors = get_nb(current_t)
        step_best_t = current_t
        step_best_s = current_s
        improved = False

        for nb in neighbors:
            s = scorer.score_tree(nb)
            if s < step_best_s - tol:
                step_best_s = s
                step_best_t = nb
                improved = True

        return step_best_t, step_best_s, improved

    phases = ["nni"] if mode == "nni" else (["spr"] if mode == "spr" else ["nni", "spr"])

    for phase in phases:
        for it in range(max_iter):
            if best_score <= tol:
                break
            next_tree, next_score, improved = _step(current_tree, best_score, phase)
            if not improved:
                break
            current_tree = next_tree
            best_score = next_score
            topology_changed = True
            if verbose:
                print(f"[{phase.upper()} iter {it+1}] Improved penalty: {best_score:.6f}")

    if topology_changed:
        fit_ols_branch_lengths(current_tree, scorer.D_obs, scorer.taxa_names)

    return tree_to_newick(current_tree), best_score

