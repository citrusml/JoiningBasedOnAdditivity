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


def tree_to_newick(tree: Tree) -> str:
    """Converts a Biopython Tree to a Newick string."""
    stream = io.StringIO()
    Phylo.write(tree, stream, "newick")
    return stream.getvalue().strip()


def hill_climb_search(
    initial_newick: str,
    scorer: TreeScorer,
    mode: str = "hybrid",
    max_iter: int = 50,
    verbose: bool = False,
) -> Tuple[str, float]:
    """
    Executes hill-climbing topology search to minimize quartet additivity penalty.

    Parameters
    ----------
    initial_newick : str
        Starting tree in Newick format (unitless).
    scorer : TreeScorer
        TreeScorer instance configured with observed distance matrix (unitless).
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
            Final quartet additivity penalty score in [0, 1].
    """
    tree = Phylo.read(io.StringIO(initial_newick), "newick")
    outgroup = scorer.taxa_names[0]
    tree.root_with_outgroup(outgroup)

    current_tree = tree
    best_score = scorer.score_tree(current_tree)

    def _step(current_t: Tree, current_s: float, search_type: str) -> Tuple[Tree, float, bool]:
        get_nb = get_nni_neighbors if search_type == "nni" else get_spr_neighbors
        neighbors = get_nb(current_t)
        step_best_t = current_t
        step_best_s = current_s
        improved = False

        for nb in neighbors:
            s = scorer.score_tree(nb)
            if s < step_best_s - 1e-7:
                step_best_s = s
                step_best_t = nb
                improved = True

        return step_best_t, step_best_s, improved

    phases = ["nni"] if mode == "nni" else (["spr"] if mode == "spr" else ["nni", "spr"])

    for phase in phases:
        for it in range(max_iter):
            if best_score <= 1e-7:
                break
            next_tree, next_score, improved = _step(current_tree, best_score, phase)
            if not improved:
                break
            current_tree = next_tree
            best_score = next_score
            if verbose:
                print(f"[{phase.upper()} iter {it+1}] Improved penalty: {best_score:.6f}")

    return tree_to_newick(current_tree), best_score

