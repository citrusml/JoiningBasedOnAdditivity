"""
Phylogenetic tree simulator adhering to Matsui & Iwasaki (2020, Systematic Biology).

Generates synthetic phylogenetic trees using BioPerl (Bio::Tree::RandomFactory,
rand_yule_c_tree) via src/generate_tree.pl (matching PhyloTree_Methods_Comparison),
with a pure-Python coalescent fallback, and computes exact patristic distance matrices.
"""

from typing import List, Tuple, Optional, Dict
import os
import sys
import math
import random
import shutil
import tempfile
import subprocess
import numpy as np


class SimTreeNode:
    """
    Binary tree node representation for synthetic phylogenetic trees.

    Parameters
    ----------
    name : str, default=""
        Taxon name for leaf nodes (e.g., 'T1', 'T2'), or empty string for internal nodes [unitless].
    length : float, default=0.0
        Branch length leading from this node to its parent [substitutions/site].
    left : Optional['SimTreeNode'], default=None
        Left child node.
    right : Optional['SimTreeNode'], default=None
        Right child node.
    """

    def __init__(
        self,
        name: str = "",
        length: float = 0.0,
        left: Optional["SimTreeNode"] = None,
        right: Optional["SimTreeNode"] = None,
    ) -> None:
        self.name: str = name
        self.length: float = length
        self.left: Optional["SimTreeNode"] = left
        self.right: Optional["SimTreeNode"] = right

    def is_leaf(self) -> bool:
        """
        Checks if the node is a leaf.

        Returns
        -------
        bool
            True if leaf node, False otherwise [unitless].
        """
        return self.left is None and self.right is None

    def to_newick(self) -> str:
        """
        Converts subtree to Newick string representation (without trailing semicolon).

        Returns
        -------
        str
            Newick format string [unitless].
        """
        if self.is_leaf():
            return f"{self.name}:{self.length:.6f}"
        left_str = self.left.to_newick() if self.left else ""
        right_str = self.right.to_newick() if self.right else ""
        if self.length > 0.0:
            return f"({left_str},{right_str}):{self.length:.6f}"
        return f"({left_str},{right_str})"


def parse_newick_to_tree(newick_str: str) -> SimTreeNode:
    """
    Parses a binary Newick tree string into a SimTreeNode hierarchy.

    Parameters
    ----------
    newick_str : str
        Newick formatted tree string [unitless].

    Returns
    -------
    SimTreeNode
        Root node of the parsed binary tree.
    """
    s = newick_str.strip()
    if s.endswith(";"):
        s = s[:-1].strip()

    def _parse_subtree(sub: str) -> SimTreeNode:
        sub = sub.strip()
        if not sub.startswith("("):
            # Leaf node: Name:Length or Name
            if ":" in sub:
                name_part, len_part = sub.rsplit(":", 1)
                return SimTreeNode(name=name_part.strip(), length=float(len_part.strip()))
            return SimTreeNode(name=sub.strip(), length=0.0)

        # Find matching closing parenthesis for the opening '(' at index 0
        depth = 0
        close_idx = -1
        top_comma_idx = -1
        for i, ch in enumerate(sub):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    close_idx = i
                    break
            elif ch == "," and depth == 1 and top_comma_idx == -1:
                top_comma_idx = i

        if close_idx == -1 or top_comma_idx == -1:
            raise ValueError(f"Malformed binary Newick substring: {sub}")

        left_sub = sub[1:top_comma_idx]
        right_sub = sub[top_comma_idx + 1 : close_idx]
        remainder = sub[close_idx + 1 :].strip()

        branch_len = 0.0
        name_part = ""
        if remainder:
            if ":" in remainder:
                name_str, len_str = remainder.rsplit(":", 1)
                name_part = name_str.strip()
                branch_len = float(len_str.strip())
            else:
                name_part = remainder.strip()

        left_node = _parse_subtree(left_sub)
        right_node = _parse_subtree(right_sub)
        return SimTreeNode(name=name_part, length=branch_len, left=left_node, right=right_node)

    root = _parse_subtree(s)
    # Root branch length is conventionally 0.0 for unrooted/rooted patristic evaluation
    root.length = 0.0
    return root


def get_perl_executable() -> Optional[str]:
    """
    Finds a Perl executable capable of running BioPerl (Bio::Tree::RandomFactory).

    Returns
    -------
    Optional[str]
        Path to Perl executable with BioPerl installed, or None if unavailable.
    """
    candidates: List[str] = []

    # 1. Same directory as current Python interpreter (conda/micromamba env)
    env_perl = os.path.join(os.path.dirname(sys.executable), "perl")
    if os.path.isfile(env_perl) and os.access(env_perl, os.X_OK):
        candidates.append(env_perl)

    # 2. CONDA_PREFIX if set
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        cp_perl = os.path.join(conda_prefix, "bin", "perl")
        if os.path.isfile(cp_perl) and os.access(cp_perl, os.X_OK) and cp_perl not in candidates:
            candidates.append(cp_perl)

    # 3. Standard phylomethod_env location
    home_env_perl = os.path.expanduser("~/.micromamba/envs/phylomethod_env/bin/perl")
    if os.path.isfile(home_env_perl) and os.access(home_env_perl, os.X_OK) and home_env_perl not in candidates:
        candidates.append(home_env_perl)

    # 4. perl from PATH
    which_perl = shutil.which("perl")
    if which_perl and which_perl not in candidates:
        candidates.append(which_perl)

    for p in candidates:
        try:
            res = subprocess.run(
                [p, "-MBio::Tree::RandomFactory", "-e", "1"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if res.returncode == 0:
                return p
        except Exception:
            continue

    return None


def generate_tree_with_bioperl(
    n_taxa: int = 15,
    divergence: float = 0.5,
    seed: Optional[int] = None,
    rate_sd: float = 0.0,
    lba_ratio: float = 1.0,
    outtree: Optional[str] = None,
) -> str:
    """
    Generates a Newick tree strictly using BioPerl Bio::Tree::RandomFactory via src/generate_tree.pl.

    Parameters
    ----------
    n_taxa : int, default=15
        Number of terminal taxa N [unitless].
    divergence : float, default=0.5
        Evolutionary distance scaling factor D [substitutions/site].
    seed : Optional[int], default=None
        Random seed for reproducibility [unitless].
    rate_sd : float, default=0.0
        Rate heterogeneity standard deviation across branches [unitless].
    lba_ratio : float, default=1.0
        Long-branch attraction ratio on T1 and TN [unitless].
    outtree : Optional[str], default=None
        Optional file path to save the Newick tree.

    Returns
    -------
    str
        Generated Newick tree string ending with ';' [unitless].
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    pl_script = os.path.join(script_dir, "generate_tree.pl")
    if not os.path.exists(pl_script):
        raise FileNotFoundError(f"BioPerl tree generator script not found at {pl_script}")

    perl_bin = get_perl_executable()
    if perl_bin is None:
        raise RuntimeError(
            "BioPerl (Bio::Tree::RandomFactory) is not available in any found perl executable."
        )

    cleanup_temp = False
    if outtree is None:
        fd, target_path = tempfile.mkstemp(suffix=".nwk")
        os.close(fd)
        cleanup_temp = True
    else:
        target_path = outtree

    cmd = [
        perl_bin,
        pl_script,
        "--taxa", str(n_taxa),
        "--scale", str(divergence),
        "--rate_sd", str(rate_sd),
        "--lba_ratio", str(lba_ratio),
        "--outtree", target_path,
    ]
    if seed is not None:
        cmd.extend(["--seed", str(seed)])

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            err_msg = proc.stderr.strip() or proc.stdout.strip()
            raise RuntimeError(
                f"BioPerl execution failed (code {proc.returncode}) using {perl_bin}:\n{err_msg}"
            )
        with open(target_path, "r", encoding="utf-8") as f:
            newick = f.read().strip()
        return newick
    finally:
        if cleanup_temp and os.path.exists(target_path):
            try:
                os.remove(target_path)
            except OSError:
                pass


def sample_branch_length(
    rng: random.Random,
    divergence: float = 0.5,
    min_length: float = 1e-6,
) -> float:
    """
    Samples branch length according to Matsui & Iwasaki (2020) log-distribution.

    Formula:
        l = (1.0 - ln(u * (e - 1.0) + 1.0)) * divergence

    Parameters
    ----------
    rng : random.Random
        Random number generator instance.
    divergence : float, default=0.5
        Divergence scaling parameter D [substitutions/site].
    min_length : float, default=1e-6
        Minimum branch length guard [substitutions/site].

    Returns
    -------
    float
        Sampled branch length [substitutions/site].
    """
    u = rng.random()
    raw_len = (1.0 - math.log(u * (math.e - 1.0) + 1.0)) * divergence
    return max(min_length, raw_len)


def _generate_paper_yule_tree_python(
    n_taxa: int = 15,
    divergence: float = 0.5,
    seed: Optional[int] = None,
    rate_sd: float = 0.0,
    lba_ratio: float = 1.0,
    min_length: float = 1e-6,
) -> SimTreeNode:
    """
    Pure-Python fallback implementing BioPerl's rand_yule_c_tree coalescent pairing algorithm.

    Parameters
    ----------
    n_taxa : int, default=15
        Number of terminal taxa N [unitless].
    divergence : float, default=0.5
        Divergence scale parameter D [substitutions/site].
    seed : Optional[int], default=None
        Random seed [unitless].
    rate_sd : float, default=0.0
        Rate heterogeneity standard deviation [unitless].
    lba_ratio : float, default=1.0
        Long-branch attraction ratio [unitless].
    min_length : float, default=1e-6
        Minimum allowable branch length [substitutions/site].

    Returns
    -------
    SimTreeNode
        Root node of the generated binary tree.
    """
    rng = random.Random(seed)

    # 1. Sample N - 1 coalescence times and sort ascending (matching BioPerl rand_yule_c_tree)
    times: List[float] = []
    for _ in range(n_taxa - 1):
        u = rng.random()
        base_t = 1.0 - math.log(u * (math.e - 1.0) + 1.0)
        times.append(base_t * divergence)
    times.sort()

    # 2. Initialize leaves: T1, T2, ..., TN
    nodes: List[SimTreeNode] = [SimTreeNode(name=f"T{i + 1}", length=0.0) for i in range(n_taxa)]

    # 3. Backward Yule coalescence
    while len(nodes) > 1:
        t = times.pop(0)

        idx1, idx2 = sorted(rng.sample(range(len(nodes)), 2), reverse=True)
        v_node = nodes.pop(idx1)
        u_node = nodes.pop(idx2)

        lu = t
        lv = t

        if rate_sd > 0.0:
            lu = max(min_length, lu + rng.gauss(0.0, rate_sd * divergence))
            lv = max(min_length, lv + rng.gauss(0.0, rate_sd * divergence))

        if lba_ratio > 1.0:
            if u_node.name in ("T1", f"T{n_taxa}"):
                lu *= lba_ratio
            if v_node.name in ("T1", f"T{n_taxa}"):
                lv *= lba_ratio

        u_node.length = max(min_length, lu)
        v_node.length = max(min_length, lv)

        parent = SimTreeNode(name="", length=0.0, left=u_node, right=v_node)
        nodes.append(parent)

    return nodes[0]


def generate_backward_yule_tree(
    n_taxa: int,
    divergence: float = 0.5,
    seed: Optional[int] = None,
    rate_sd: float = 0.0,
    lba_ratio: float = 1.0,
    engine: str = "bioperl",
) -> Tuple[SimTreeNode, List[str]]:
    """
    Generates a synthetic phylogenetic tree via backward Yule coalescence using BioPerl
    (Bio::Tree::RandomFactory, matching PhyloTree_Methods_Comparison).

    Parameters
    ----------
    n_taxa : int
        Number of taxa N [unitless].
    divergence : float, default=0.5
        Divergence scale parameter D [substitutions/site].
    seed : Optional[int], default=None
        Random seed for reproducibility [unitless].
    rate_sd : float, default=0.0
        Rate heterogeneity standard deviation across branches [unitless].
    lba_ratio : float, default=1.0
        Long-branch attraction ratio [unitless].
    engine : str, default="bioperl"
        Simulation backend ('bioperl' or 'python'). If 'bioperl' is selected and BioPerl
        is installed, calls src/generate_tree.pl; otherwise falls back to 'python'.

    Returns
    -------
    Tuple[SimTreeNode, List[str]]
        root : SimTreeNode
            Root node of generated binary tree.
        taxa_names : List[str]
            List of leaf taxon names in index order ['T1', 'T2', ..., 'TN'].
    """
    if n_taxa < 2:
        raise ValueError(f"n_taxa must be at least 2, got {n_taxa}")
    if divergence <= 0.0:
        raise ValueError(f"divergence must be strictly positive, got {divergence}")

    taxa_names: List[str] = [f"T{i + 1}" for i in range(n_taxa)]

    if engine == "bioperl" and get_perl_executable() is not None:
        newick_str = generate_tree_with_bioperl(
            n_taxa=n_taxa,
            divergence=divergence,
            seed=seed,
            rate_sd=rate_sd,
            lba_ratio=lba_ratio,
        )
        root = parse_newick_to_tree(newick_str)
        return root, taxa_names

    if engine not in ("bioperl", "python"):
        raise ValueError(f"Unsupported engine '{engine}'. Choose 'bioperl' or 'python'.")

    root = _generate_paper_yule_tree_python(
        n_taxa=n_taxa,
        divergence=divergence,
        seed=seed,
        rate_sd=rate_sd,
        lba_ratio=lba_ratio,
    )
    return root, taxa_names


def compute_patristic_distance_matrix(
    root: SimTreeNode,
    taxa_names: List[str],
) -> np.ndarray:
    """
    Computes true tree-metric (patristic distance matrix) from tree.

    Parameters
    ----------
    root : SimTreeNode
        Root node of tree.
    taxa_names : List[str]
        Ordered list of taxa names [unitless].

    Returns
    -------
    np.ndarray
        True patristic distance matrix of shape (n_taxa, n_taxa) [substitutions/site].
    """
    n = len(taxa_names)
    taxa_map: Dict[str, int] = {name: idx for idx, name in enumerate(taxa_names)}
    D_true = np.zeros((n, n), dtype=np.float64)

    def get_leaf_distances(node: SimTreeNode) -> List[Tuple[int, float]]:
        if node.is_leaf():
            idx = taxa_map[node.name]
            return [(idx, 0.0)]

        left_leaves: List[Tuple[int, float]] = []
        right_leaves: List[Tuple[int, float]] = []

        if node.left:
            for l_idx, d in get_leaf_distances(node.left):
                left_leaves.append((l_idx, d + node.left.length))
        if node.right:
            for r_idx, d in get_leaf_distances(node.right):
                right_leaves.append((r_idx, d + node.right.length))

        for l_idx, d_l in left_leaves:
            for r_idx, d_r in right_leaves:
                dist = d_l + d_r
                D_true[l_idx, r_idx] = dist
                D_true[r_idx, l_idx] = dist

        return left_leaves + right_leaves

    get_leaf_distances(root)
    return D_true
