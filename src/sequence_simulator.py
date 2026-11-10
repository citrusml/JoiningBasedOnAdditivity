"""
Amino acid sequence evolution simulator under WAG + Gamma and pairwise alignment distance estimator.

Implements:
1. WAG + Gamma (discrete 4-category Gamma rate heterogeneity) amino acid sequence evolution.
2. Pairwise sequence alignment (PSA) using Biopython PairwiseAligner with BLOSUM62 substitution matrix.
3. Global alignment (Needleman-Wunsch) and Local alignment (Smith-Waterman).
4. Amino acid Poisson correction distance calculation: d = -ln(1 - p).
"""

from typing import Dict, List, Tuple, Optional, Union
import os
import sys
import math
import shutil
import tempfile
import subprocess
import numpy as np

import scipy.linalg
import scipy.stats
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from src.tree_simulator import SimTreeNode


AMINO_ACIDS = [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
]
AA_TO_INT = {aa: idx for idx, aa in enumerate(AMINO_ACIDS)}


class WAGModel:
    """
    Whelan and Goldman (2001) amino acid substitution model with discrete Gamma rate heterogeneity.

    Parameters
    ----------
    wag_path : Optional[str], default=None
        Path to wag.dat file. If None, looks in current directory or src/wag.dat.
    alpha : float, default=1.0
        Shape parameter for discrete Gamma distribution across sites (unitless).
    n_categories : int, default=4
        Number of discrete rate categories for Gamma distribution.
    """

    def __init__(
        self,
        wag_path: Optional[str] = None,
        alpha: float = 1.0,
        n_categories: int = 4,
    ) -> None:
        self.alpha: float = alpha
        self.n_categories: int = n_categories

        if wag_path is None:
            base_dir = os.path.dirname(__file__)
            wag_path = os.path.join(base_dir, "wag.dat")
            if not os.path.exists(wag_path):
                wag_path = "src/wag.dat"

        self._load_wag_matrix(wag_path)
        self._setup_gamma_rates()
        self._decompose_rate_matrix()

    def _load_wag_matrix(self, path: str) -> None:
        """Parses exchangeabilities S and frequencies pi from wag.dat."""
        with open(path, "r") as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith("//")]

        S = np.zeros((20, 20), dtype=np.float64)
        for row_idx in range(1, 20):
            vals = [float(x) for x in lines[row_idx - 1].split()]
            for col_idx, v in enumerate(vals):
                S[row_idx, col_idx] = v
                S[col_idx, row_idx] = v

        pi = np.array([float(x) for x in lines[19].split()], dtype=np.float64)
        self.pi: np.ndarray = pi / np.sum(pi)

        # Rate matrix Q_ij = S_ij * pi_j
        Q = S * self.pi[None, :]
        for i in range(20):
            Q[i, i] = -np.sum(Q[i, :]) + Q[i, i]

        # Normalize rate matrix so mean substitution rate is 1.0
        mu = -np.sum(self.pi * np.diag(Q))
        self.Q: np.ndarray = Q / mu

    def _setup_gamma_rates(self) -> None:
        """Computes discrete Gamma rate categories using Yang (1994) quantile medians."""
        k = self.n_categories
        probs = (2 * np.arange(k) + 1) / (2.0 * k)
        raw_rates = scipy.stats.gamma.ppf(probs, a=self.alpha, scale=1.0 / self.alpha)
        self.gamma_rates: np.ndarray = raw_rates / np.mean(raw_rates)

    def _decompose_rate_matrix(self) -> None:
        """Eigendecomposes Q for fast continuous-time matrix exponentiation."""
        w, V = scipy.linalg.eig(self.Q)
        self.eigenvals: np.ndarray = w
        self.V: np.ndarray = V
        self.V_inv: np.ndarray = np.linalg.inv(V)

    def get_transition_matrices(self, branch_length: float) -> np.ndarray:
        """
        Computes transition probability matrices for each Gamma rate category along a branch.

        Parameters
        ----------
        branch_length : float
            Branch length t (expected substitutions per site).

        Returns
        -------
        np.ndarray
            Array of shape (n_categories, 20, 20) with transition matrices P_k.
        """
        k = self.n_categories
        P = np.zeros((k, 20, 20), dtype=np.float64)
        if branch_length <= 0.0:
            for c in range(k):
                P[c] = np.eye(20)
            return P

        for c in range(k):
            effective_t = branch_length * self.gamma_rates[c]
            exp_diag = np.exp(self.eigenvals * effective_t)
            p_mat = np.real(self.V @ np.diag(exp_diag) @ self.V_inv)
            p_mat = np.maximum(0.0, p_mat)
            p_mat = p_mat / np.sum(p_mat, axis=1, keepdims=True)
            P[c] = p_mat

        return P

    def get_marginal_log_prob_grid(
        self,
        max_dist: float = 5.0,
        n_grid: int = 500,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Precomputes marginal transition log-probability matrices across a distance grid
        under the discrete Gamma WAG model.

        Parameters
        ----------
        max_dist : float, default=5.0
            Maximum evolutionary distance on the grid (substitutions/site).
        n_grid : int, default=500
            Number of grid points (unitless).

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            - grid_d : 1D array of shape (n_grid,) in substitutions/site.
            - log_P_flat : 2D array of shape (n_grid, 400) containing ln P_ab(d) (unitless).
        """
        cache_key = (max_dist, n_grid)
        if hasattr(self, "_grid_cache") and cache_key in self._grid_cache:
            return self._grid_cache[cache_key]

        if not hasattr(self, "_grid_cache"):
            self._grid_cache: Dict[Tuple[float, int], Tuple[np.ndarray, np.ndarray]] = {}

        # Combine dense small distances and linear grid up to max_dist
        grid_d = np.concatenate([
            np.linspace(1e-5, 0.2, n_grid // 3, endpoint=False),
            np.linspace(0.2, max_dist, n_grid - (n_grid // 3)),
        ])

        log_P_flat = np.zeros((len(grid_d), 400), dtype=np.float64)
        for idx, d in enumerate(grid_d):
            P_cats = self.get_transition_matrices(float(d))  # (k, 20, 20)
            P_marg = np.mean(P_cats, axis=0)  # (20, 20)
            P_marg = np.maximum(P_marg, 1e-15)
            log_P_flat[idx] = np.log(P_marg).ravel()

        self._grid_cache[cache_key] = (grid_d, log_P_flat)
        return grid_d, log_P_flat

    def marginal_log_likelihood(self, counts_flat: np.ndarray, dist: float) -> float:
        """
        Computes exact WAG+Gamma log-likelihood for a single flattened 20x20 transition count vector.

        Parameters
        ----------
        counts_flat : np.ndarray
            Flattened 20x20 amino acid transition counts of shape (400,) (sites).
        dist : float
            Candidate evolutionary distance (substitutions/site).

        Returns
        -------
        float
            Log-likelihood value (unitless).
        """
        P_cats = self.get_transition_matrices(dist)
        P_marg = np.maximum(np.mean(P_cats, axis=0), 1e-15)
        return float(np.dot(counts_flat, np.log(P_marg).ravel()))


_INDELIBLE_BIN_CACHE: Optional[str] = None


def find_or_build_indelible_binary() -> str:
    """
    Locates a working INDELible v1.03 binary for the current OS/architecture,
    or builds one via `scripts/install_indelible.sh` if necessary.

    Returns
    -------
    str
        Absolute path to the executable INDELible binary.
    """
    global _INDELIBLE_BIN_CACHE
    if _INDELIBLE_BIN_CACHE and os.path.isfile(_INDELIBLE_BIN_CACHE) and os.access(_INDELIBLE_BIN_CACHE, os.X_OK):
        return _INDELIBLE_BIN_CACHE

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    candidates = [
        shutil.which("indelible"),
        os.path.join(repo_root, "bin", "indelible"),
        os.path.abspath(os.path.join(repo_root, "..", "PhyloTree_Methods_Comparison", "bin", "indelible")),
    ]

    for cand in candidates:
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            try:
                res = subprocess.run(
                    [cand],
                    input="\n",
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
                out_str = (res.stdout or "") + (res.stderr or "")
                if "INDELible" in out_str:
                    _INDELIBLE_BIN_CACHE = cand
                    return cand
            except Exception:
                continue

    # Build via scripts/install_indelible.sh
    install_script = os.path.join(repo_root, "scripts", "install_indelible.sh")
    if os.path.isfile(install_script):
        subprocess.run(["bash", install_script], check=True)
        built_bin = os.path.join(repo_root, "bin", "indelible")
        if os.path.isfile(built_bin) and os.access(built_bin, os.X_OK):
            _INDELIBLE_BIN_CACHE = built_bin
            return built_bin

    raise RuntimeError("Could not find or build a working INDELible binary.")


def simulate_indelible_wag_gamma_sequences(
    root: SimTreeNode,
    seq_length: int = 500,
    alpha: float = 1.0,
    indel_rate: float = 0.05,
    indel_power: float = 1.7,
    indel_max_len: int = 50,
    seed: Optional[int] = None,
) -> Dict[str, str]:
    """
    Simulates amino acid sequence evolution with Zipfian (power-law) insertions and
    deletions using INDELible v1.03 under WAG + discrete 4-category Gamma, matching
    Matsui & Iwasaki (2020) and `PhyloTree_Methods_Comparison`.

    Parameters
    ----------
    root : SimTreeNode
        Root node of the phylogenetic tree.
    seq_length : int, default=500
        Root amino acid sequence length L (residues).
    alpha : float, default=1.0
        Shape parameter of discrete 4-category Gamma rate heterogeneity (unitless).
    indel_rate : float, default=0.05
        Insertion and deletion rate relative to substitution rate (e.g., 0.05 or 0.10, unitless).
    indel_power : float, default=1.7
        Zipfian power-law exponent a for `POW a M` indel length distribution (unitless).
    indel_max_len : int, default=50
        Maximum indel length M for `POW a M` distribution (residues).
    seed : Optional[int], default=None
        Random seed for deterministic simulation (unitless).

    Returns
    -------
    Dict[str, str]
        Mapping from taxon leaf name to unaligned evolved amino acid sequence string.
    """
    indelible_bin = find_or_build_indelible_binary()
    tree_str = root.to_newick().strip()
    if not tree_str.endswith(";"):
        tree_str += ";"

    expected_leaves: List[str] = []

    def _collect_leaves(node: SimTreeNode) -> None:
        if node.is_leaf():
            expected_leaves.append(node.name)
            return
        if node.left:
            _collect_leaves(node.left)
        if node.right:
            _collect_leaves(node.right)

    _collect_leaves(root)
    n_expected = len(expected_leaves)
    base_seed = int(seed) if seed is not None else 12345

    for attempt in range(5):
        cur_seed = abs(base_seed + attempt * 100003) % 2147483647
        if cur_seed == 0:
            cur_seed = 42

        with tempfile.TemporaryDirectory(prefix="indelible_sim_") as tmpdir:
            control_content = (
                "[TYPE] AMINOACID 1\n"
                "[SETTINGS]\n"
                f"  [randomseed] {cur_seed}\n"
                "  [output] FASTA\n"
                "[MODEL] mymodel\n"
                "  [submodel] WAG\n"
                f"  [rates] 0 {alpha} 4\n"
                f"  [insertmodel] POW {indel_power} {indel_max_len}\n"
                f"  [deletemodel] POW {indel_power} {indel_max_len}\n"
                f"  [insertrate] {indel_rate}\n"
                f"  [deleterate] {indel_rate}\n"
                f"[TREE] mytree {tree_str}\n"
                f"[PARTITIONS] mypart [mytree mymodel {seq_length}]\n"
                "[EVOLVE] mypart 1 sim\n"
            )
            ctrl_path = os.path.join(tmpdir, "control.txt")
            with open(ctrl_path, "w") as f:
                f.write(control_content)

            proc = subprocess.run(
                [indelible_bin],
                cwd=tmpdir,
                capture_output=True,
                text=True,
            )
            fas_path = os.path.join(tmpdir, "sim.fas")
            if proc.returncode == 0 and os.path.isfile(fas_path):
                seqs: Dict[str, str] = {}
                for rec in SeqIO.parse(fas_path, "fasta"):
                    s_clean = str(rec.seq).replace("-", "").replace(" ", "").strip().upper()
                    if len(s_clean) > 0:
                        seqs[rec.id] = s_clean
                if len(seqs) == n_expected and all(k in seqs for k in expected_leaves):
                    return seqs

    raise RuntimeError(
        f"INDELible simulation failed to produce valid non-empty sequences for all {n_expected} taxa after 5 attempts."
    )


def simulate_wag_gamma_sequences(
    root: SimTreeNode,
    seq_length: int = 300,
    alpha: float = 1.0,
    model: Optional[WAGModel] = None,
    seed: Optional[int] = None,
    indel_rate: float = 0.0,
    indel_power: float = 1.7,
    indel_max_len: int = 50,
) -> Dict[str, str]:
    """
    Simulates amino acid sequences along a tree under WAG + Gamma, optionally with
    Zipfian insertions and deletions via INDELible when `indel_rate > 0.0`.

    Parameters
    ----------
    root : SimTreeNode
        Root node of the phylogenetic tree.
    seq_length : int, default=300
        Length of root amino acid sequence L (residues).
    alpha : float, default=1.0
        Shape parameter of Gamma rate distribution (unitless).
    model : Optional[WAGModel], default=None
        Pre-initialized WAGModel. If None, initialized with alpha.
    seed : Optional[int], default=None
        Random seed (unitless).
    indel_rate : float, default=0.0
        Insertion/deletion rate relative to substitution rate (0.0 for substitution-only,
        or e.g. 0.05 / 0.10 for INDELible Zipfian indels, unitless).
    indel_power : float, default=1.7
        Zipfian power-law exponent a for `POW a M` indel length model (unitless).
    indel_max_len : int, default=50
        Maximum indel length M for `POW a M` indel length model (residues).

    Returns
    -------
    Dict[str, str]
        Mapping from taxon leaf name to amino acid sequence string.
    """
    if indel_rate > 0.0:
        return simulate_indelible_wag_gamma_sequences(
            root=root,
            seq_length=seq_length,
            alpha=alpha,
            indel_rate=indel_rate,
            indel_power=indel_power,
            indel_max_len=indel_max_len,
            seed=seed,
        )

    if model is None or model.alpha != alpha:
        model = WAGModel(alpha=alpha)

    rng = np.random.default_rng(seed)

    # 1. Assign each site to a discrete Gamma rate category (0..3) uniformly
    site_categories = rng.integers(0, model.n_categories, size=seq_length)

    # 2. Sample root sequence from equilibrium distribution pi
    root_seq = rng.choice(20, size=seq_length, p=model.pi)

    leaf_sequences_int: Dict[str, np.ndarray] = {}

    def evolve_branch(parent_seq: np.ndarray, branch_length: float) -> np.ndarray:
        if branch_length <= 0.0:
            return parent_seq.copy()

        trans_mats = model.get_transition_matrices(branch_length)  # (4, 20, 20)
        child_seq = np.empty(seq_length, dtype=np.int64)

        for c in range(model.n_categories):
            mask = (site_categories == c)
            n_sites = int(np.sum(mask))
            if n_sites == 0:
                continue

            parent_states = parent_seq[mask]
            # Transition probabilities for these sites: shape (n_sites, 20)
            probs = trans_mats[c][parent_states, :]

            # Cumulative sum for fast inverse CDF sampling
            cumprobs = np.cumsum(probs, axis=1)
            u = rng.random(size=(n_sites, 1))
            new_states = np.argmax(cumprobs > u, axis=1)
            child_seq[mask] = new_states

        return child_seq

    def traverse(node: SimTreeNode, current_seq: np.ndarray) -> None:
        if node.is_leaf():
            leaf_sequences_int[node.name] = current_seq
            return

        if node.left:
            left_seq = evolve_branch(current_seq, node.left.length)
            traverse(node.left, left_seq)
        if node.right:
            right_seq = evolve_branch(current_seq, node.right.length)
            traverse(node.right, right_seq)

    traverse(root, root_seq)

    # Convert integer sequences to amino acid letter strings
    leaf_sequences_str: Dict[str, str] = {}
    for name, seq_arr in leaf_sequences_int.items():
        leaf_sequences_str[name] = "".join(AMINO_ACIDS[idx] for idx in seq_arr)

    return leaf_sequences_str


def create_pairwise_aligner(mode: str = "global") -> PairwiseAligner:
    """
    Creates Biopython PairwiseAligner configured with BLOSUM62.

    Parameters
    ----------
    mode : str, default='global'
        'global' for Needleman-Wunsch, 'local' for Smith-Waterman.

    Returns
    -------
    PairwiseAligner
        Configured aligner.
    """
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -1.0
    aligner.mode = mode
    return aligner


def compute_pairwise_poisson_distance_matrix(
    sequences: Dict[str, str],
    taxa_names: List[str],
    mode: str = "global",
    max_dist: float = 3.0,
) -> np.ndarray:
    """
    Computes pairwise amino acid Poisson correction distance matrix using BLOSUM62 alignment.

    Formula:
        d = -ln(1 - p)
    where p = mismatches / aligned_length (excluding gap positions).

    Parameters
    ----------
    sequences : Dict[str, str]
        Taxon amino acid sequences.
    taxa_names : List[str]
        Ordered list of taxon names.
    mode : str, default='global'
        Alignment mode: 'global' or 'local'.
    max_dist : float, default=3.0
        Maximum distance cap when p >= 0.95 (substitutions/site).

    Returns
    -------
    np.ndarray
        Estimated Poisson distance matrix D_est of shape (n_taxa, n_taxa).
    """
    n = len(taxa_names)
    D_est = np.zeros((n, n), dtype=np.float64)
    aligner = create_pairwise_aligner(mode=mode)

    for i in range(n):
        s_i = sequences[taxa_names[i]]
        for j in range(i + 1, n):
            s_j = sequences[taxa_names[j]]

            alignments = aligner.align(s_i, s_j)
            top_aln = alignments[0]

            # Extract aligned string representations with gaps
            aln_i, aln_j = top_aln[0], top_aln[1]

            mismatches = 0
            aligned_sites = 0

            for char_i, char_j in zip(aln_i, aln_j):
                if char_i != "-" and char_j != "-":
                    aligned_sites += 1
                    if char_i != char_j:
                        mismatches += 1

            if aligned_sites == 0:
                p = 0.0
            else:
                p = mismatches / aligned_sites

            if p >= 0.95:
                d = max_dist
            elif p <= 0.0:
                d = 0.0
            else:
                d = -math.log(1.0 - p)

            D_est[i, j] = d
            D_est[j, i] = d

    return D_est


def compute_mafft_msa_poisson_distance_matrix(
    sequences: Dict[str, str],
    taxa_names: List[str],
    max_dist: float = 3.0,
) -> np.ndarray:
    """
    Computes pairwise amino acid Poisson correction distances from MAFFT Multiple Sequence Alignment.

    Parameters
    ----------
    sequences : Dict[str, str]
        Dictionary mapping taxon name to unaligned amino acid sequence.
    taxa_names : List[str]
        Ordered list of taxon names.
    max_dist : float, default=3.0
        Cap distance when p >= 0.95 (substitutions/site).

    Returns
    -------
    np.ndarray
        Estimated Poisson distance matrix D_msa of shape (n_taxa, n_taxa).
    """
    n = len(taxa_names)
    D_msa = np.zeros((n, n), dtype=np.float64)

    with tempfile.TemporaryDirectory() as tmpdir:
        in_fasta = os.path.join(tmpdir, "in_seqs.fasta")
        out_fasta = os.path.join(tmpdir, "out_msa.fasta")

        with open(in_fasta, "w") as f:
            for name in taxa_names:
                f.write(f">{name}\n{sequences[name]}\n")

        # Run MAFFT MSA
        subprocess.run(
            ["mafft", "--auto", "--quiet", in_fasta],
            stdout=open(out_fasta, "w"),
            stderr=subprocess.DEVNULL,
            check=True,
        )

        msa_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(out_fasta, "fasta")}

    for i in range(n):
        s_i = msa_dict[taxa_names[i]]
        for j in range(i + 1, n):
            s_j = msa_dict[taxa_names[j]]

            mismatches = 0
            valid_sites = 0

            for c_i, c_j in zip(s_i, s_j):
                if c_i != "-" and c_j != "-":
                    valid_sites += 1
                    if c_i != c_j:
                        mismatches += 1

            if valid_sites == 0:
                p = 0.0
            else:
                p = mismatches / valid_sites

            if p >= 0.95:
                d = max_dist
            elif p <= 0.0:
                d = 0.0
            else:
                d = -math.log(1.0 - p)

            D_msa[i, j] = d
            D_msa[j, i] = d

    return D_msa


_ASCII_TO_AA = np.full(256, -1, dtype=np.int64)
for _aa_char, _aa_idx in AA_TO_INT.items():
    _ASCII_TO_AA[ord(_aa_char)] = _aa_idx


def _extract_alignment_counts_and_poisson(
    aln_i: str,
    aln_j: str,
    max_poisson_dist: float = 3.0,
) -> Tuple[np.ndarray, float]:
    """
    Extracts flattened 20x20 amino acid transition count vector and Poisson distance
    from two aligned sequence strings using vectorized byte lookup.

    Parameters
    ----------
    aln_i : str
        First aligned sequence string (amino acids and '-').
    aln_j : str
        Second aligned sequence string (amino acids and '-').
    max_poisson_dist : float, default=3.0
        Maximum cap for Poisson distance (substitutions/site).

    Returns
    -------
    Tuple[np.ndarray, float]
        - counts_flat : Flattened 20x20 transition count array of shape (400,) (sites).
        - d_poisson : Estimated Poisson distance (substitutions/site).
    """
    b_i = np.frombuffer(aln_i.encode("ascii"), dtype=np.uint8)
    b_j = np.frombuffer(aln_j.encode("ascii"), dtype=np.uint8)
    aa_i = _ASCII_TO_AA[b_i]
    aa_j = _ASCII_TO_AA[b_j]

    valid = (aa_i >= 0) & (aa_j >= 0)
    v_i = aa_i[valid]
    v_j = aa_j[valid]

    n_valid = len(v_i)
    if n_valid == 0:
        return np.zeros(400, dtype=np.float64), 0.0

    counts_flat = np.bincount(v_i * 20 + v_j, minlength=400).astype(np.float64)
    mismatches = int(np.sum(v_i != v_j))
    p = mismatches / n_valid

    if p >= 0.95:
        d_poisson = max_poisson_dist
    elif p <= 0.0:
        d_poisson = 0.0
    else:
        d_poisson = -math.log(1.0 - p)

    return counts_flat, d_poisson


def _estimate_ml_distances_from_pair_counts(
    pair_counts: np.ndarray,
    model: WAGModel,
    max_dist: float = 5.0,
) -> np.ndarray:
    """
    Vectorized Maximum Likelihood distance estimator under WAG+Gamma for multiple pairs.

    Parameters
    ----------
    pair_counts : np.ndarray
        2D array of shape (n_pairs, 400) containing flattened 20x20 transition counts (sites).
    model : WAGModel
        Configured WAG+Gamma substitution model.
    max_dist : float, default=5.0
        Maximum allowed ML distance (substitutions/site).

    Returns
    -------
    np.ndarray
        1D array of shape (n_pairs,) with estimated ML distances (substitutions/site).
    """
    grid_d, log_P_flat = model.get_marginal_log_prob_grid(max_dist=max_dist, n_grid=500)
    # Matrix multiplication: (n_pairs, 400) @ (400, n_grid) -> (n_pairs, n_grid)
    log_L_mat = pair_counts @ log_P_flat.T
    best_indices = np.argmax(log_L_mat, axis=1)

    n_pairs = pair_counts.shape[0]
    n_grid = len(grid_d)
    diag_indices = np.arange(20) * 21  # indices of (a, a) in flattened 20x20
    total_counts = np.sum(pair_counts, axis=1)
    ident_counts = np.sum(pair_counts[:, diag_indices], axis=1)

    ml_dists = np.zeros(n_pairs, dtype=np.float64)

    for p_idx in range(n_pairs):
        if total_counts[p_idx] <= 0 or ident_counts[p_idx] >= total_counts[p_idx]:
            ml_dists[p_idx] = 0.0
            continue

        m = int(best_indices[p_idx])
        if 0 < m < n_grid - 1:
            # 3-point parabolic refinement around grid peak
            x0, x1, x2 = grid_d[m - 1], grid_d[m], grid_d[m + 1]
            y0, y1, y2 = log_L_mat[p_idx, m - 1], log_L_mat[p_idx, m], log_L_mat[p_idx, m + 1]
            denom = (x0 - x1) * (x0 - x2) * (x1 - x2)
            if abs(denom) > 1e-15:
                A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
                B = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom
                if A < 0.0:
                    x_vertex = -B / (2.0 * A)
                    if x0 <= x_vertex <= x2:
                        ml_dists[p_idx] = float(x_vertex)
                        continue
        ml_dists[p_idx] = float(grid_d[m])

    return ml_dists


def compute_pairwise_wag_gamma_distance_matrix(
    sequences: Dict[str, str],
    taxa_names: List[str],
    model: Optional[WAGModel] = None,
    alpha: float = 1.0,
    align: bool = True,
    mode: str = "global",
    max_dist: float = 5.0,
) -> np.ndarray:
    """
    Computes pairwise Maximum Likelihood evolutionary distance matrix under the WAG+Gamma model.

    Parameters
    ----------
    sequences : Dict[str, str]
        Mapping from taxon name to amino acid sequence.
    taxa_names : List[str]
        Ordered list of taxon names (unitless).
    model : Optional[WAGModel], default=None
        Pre-initialized WAGModel instance.
    alpha : float, default=1.0
        Gamma shape parameter (unitless).
    align : bool, default=True
        Whether to perform pairwise sequence alignment before counting transitions.
    mode : str, default='global'
        Alignment mode ('global' or 'local').
    max_dist : float, default=5.0
        Maximum distance cap (substitutions/site).

    Returns
    -------
    np.ndarray
        Symmetric WAG+Gamma ML distance matrix of shape (n_taxa, n_taxa) in substitutions/site.
    """
    if model is None or model.alpha != alpha:
        model = WAGModel(alpha=alpha)

    n = len(taxa_names)
    n_pairs = n * (n - 1) // 2
    pair_counts = np.zeros((n_pairs, 400), dtype=np.float64)
    pair_ij: List[Tuple[int, int]] = []

    aligner = create_pairwise_aligner(mode=mode) if align else None

    p_idx = 0
    for i in range(n):
        s_i = sequences[taxa_names[i]]
        for j in range(i + 1, n):
            s_j = sequences[taxa_names[j]]
            if aligner is not None:
                top_aln = aligner.align(s_i, s_j)[0]
                aln_i, aln_j = top_aln[0], top_aln[1]
            else:
                aln_i, aln_j = s_i, s_j

            counts_flat, _ = _extract_alignment_counts_and_poisson(aln_i, aln_j)
            pair_counts[p_idx] = counts_flat
            pair_ij.append((i, j))
            p_idx += 1

    ml_dists = _estimate_ml_distances_from_pair_counts(pair_counts, model, max_dist=max_dist)

    D_wag = np.zeros((n, n), dtype=np.float64)
    for (i, j), d_val in zip(pair_ij, ml_dists):
        D_wag[i, j] = d_val
        D_wag[j, i] = d_val

    return D_wag


def compute_both_pairwise_distance_matrices(
    sequences: Dict[str, str],
    taxa_names: List[str],
    model: Optional[WAGModel] = None,
    alpha: float = 1.0,
    mode: str = "global",
    max_poisson_dist: float = 3.0,
    max_wag_dist: float = 5.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes both BLOSUM62+Poisson distance matrix and WAG+Gamma ML distance matrix
    from a single pass of pairwise global alignments.

    Parameters
    ----------
    sequences : Dict[str, str]
        Mapping from taxon name to amino acid sequence.
    taxa_names : List[str]
        Ordered list of taxon names (unitless).
    model : Optional[WAGModel], default=None
        Pre-initialized WAGModel instance.
    alpha : float, default=1.0
        Gamma shape parameter (unitless).
    mode : str, default='global'
        Alignment mode ('global' or 'local').
    max_poisson_dist : float, default=3.0
        Maximum cap for Poisson distance (substitutions/site).
    max_wag_dist : float, default=5.0
        Maximum cap for WAG+Gamma ML distance (substitutions/site).

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        - D_blosum62 : Poisson distance matrix of shape (n_taxa, n_taxa) in substitutions/site.
        - D_wag_gamma : WAG+Gamma ML distance matrix of shape (n_taxa, n_taxa) in substitutions/site.
    """
    if model is None or model.alpha != alpha:
        model = WAGModel(alpha=alpha)

    n = len(taxa_names)
    n_pairs = n * (n - 1) // 2
    pair_counts = np.zeros((n_pairs, 400), dtype=np.float64)
    pair_ij: List[Tuple[int, int]] = []
    D_blosum62 = np.zeros((n, n), dtype=np.float64)

    aligner = create_pairwise_aligner(mode=mode)

    p_idx = 0
    for i in range(n):
        s_i = sequences[taxa_names[i]]
        for j in range(i + 1, n):
            s_j = sequences[taxa_names[j]]
            top_aln = aligner.align(s_i, s_j)[0]
            aln_i, aln_j = top_aln[0], top_aln[1]

            counts_flat, d_poisson = _extract_alignment_counts_and_poisson(
                aln_i, aln_j, max_poisson_dist=max_poisson_dist
            )
            pair_counts[p_idx] = counts_flat
            pair_ij.append((i, j))
            D_blosum62[i, j] = d_poisson
            D_blosum62[j, i] = d_poisson
            p_idx += 1

    ml_dists = _estimate_ml_distances_from_pair_counts(pair_counts, model, max_dist=max_wag_dist)

    D_wag_gamma = np.zeros((n, n), dtype=np.float64)
    for (i, j), d_val in zip(pair_ij, ml_dists):
        D_wag_gamma[i, j] = d_val
        D_wag_gamma[j, i] = d_val

    return D_blosum62, D_wag_gamma


def compute_mmseqs_similarity_matrix(
    sequences: Dict[str, str],
    taxa_names: List[str],
    sensitivity: float = 7.5,
    threads: int = 1,
) -> Optional[Tuple[np.ndarray, Dict[str, float]]]:
    """
    Computes all-to-all MMseqs2 symmetrical relative sequence similarity score w_ij
    matrix and summary statistics following Matsui & Iwasaki (2020) and Exp 18:
        w_ij = (S_bits(i, j) + S_bits(j, i)) / (S_bits(i, i) + S_bits(j, j))

    Parameters
    ----------
    sequences : Dict[str, str]
        Mapping from taxon name to amino acid sequence.
    taxa_names : List[str]
        Ordered list of taxon names (unitless).
    sensitivity : float, default=7.5
        MMseqs2 sensitivity parameter (-s).
    threads : int, default=1
        Number of threads for MMseqs2 search.

    Returns
    -------
    Optional[Tuple[np.ndarray, Dict[str, float]]]
        - w_mat : (N, N) symmetrical relative similarity matrix in [0.0, 1.0] with w_ii = 1.0.
        - summary_stats : Summary statistics dictionary, or None if MMseqs2 is unavailable.
    """
    mmseqs_bin = shutil.which("mmseqs")
    if not mmseqs_bin:
        candidates = [
            os.path.join(os.path.dirname(sys.executable), "mmseqs"),
            os.path.expanduser("~/.micromamba/envs/phylomethod_env/bin/mmseqs"),
        ]
        for c in candidates:
            if os.path.isfile(c) and os.access(c, os.X_OK):
                mmseqs_bin = c
                break
    if not mmseqs_bin:
        return None

    n = len(taxa_names)
    w_mat = np.zeros((n, n), dtype=np.float64)
    np.fill_diagonal(w_mat, 1.0)
    if n < 2:
        return w_mat, {
            "sss_mean": 1.0,
            "sss_median": 1.0,
            "sss_std": 0.0,
            "sss_min": 1.0,
            "sss_max": 1.0,
            "sss_zero_frac": 0.0,
            "sss_lba_pair": 1.0,
        }

    id_map = {name: i for i, name in enumerate(taxa_names)}
    S = np.zeros((n, n), dtype=np.float64)

    with tempfile.TemporaryDirectory(prefix="mmseqs_sss_") as tmpdir:
        fasta_path = os.path.join(tmpdir, "seqs.fasta")
        db = os.path.join(tmpdir, "db")
        res = os.path.join(tmpdir, "res")
        tmp = os.path.join(tmpdir, "tmp")
        out_m8 = os.path.join(tmpdir, "out.m8")
        os.makedirs(tmp, exist_ok=True)

        with open(fasta_path, "w") as f:
            for name in taxa_names:
                f.write(f">{name}\n{sequences[name].replace('-', '').upper()}\n")

        try:
            subprocess.run(
                [mmseqs_bin, "createdb", fasta_path, db],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
            subprocess.run(
                [mmseqs_bin, "search", db, db, res, tmp, "--threads", str(threads), "-e", "10", "-s", str(sensitivity)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
            subprocess.run(
                [mmseqs_bin, "convertalis", db, db, res, out_m8, "--format-mode", "0"],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
        except Exception:
            return None

        if os.path.exists(out_m8):
            with open(out_m8, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 12:
                        q, t = parts[0], parts[1]
                        if q in id_map and t in id_map:
                            try:
                                bits = float(parts[11])
                            except ValueError:
                                bits = 0.0
                            qi, ti = id_map[q], id_map[t]
                            if bits > S[qi, ti]:
                                S[qi, ti] = bits

    scores: List[float] = []
    w_lba_pair = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            self_score = S[i, i] + S[j, j]
            comp_score = S[i, j] + S[j, i]
            w_ij = (comp_score / self_score) if self_score > 0.0 else 0.0
            w_clamped = min(1.0, max(0.0, float(w_ij)))
            w_mat[i, j] = w_clamped
            w_mat[j, i] = w_clamped
            scores.append(w_clamped)
            if i == 0 and j == n - 1:
                w_lba_pair = w_clamped

    arr = np.array(scores, dtype=np.float64)
    stats = {
        "sss_mean": float(np.mean(arr)),
        "sss_median": float(np.median(arr)),
        "sss_std": float(np.std(arr)),
        "sss_min": float(np.min(arr)),
        "sss_max": float(np.max(arr)),
        "sss_zero_frac": float(np.mean(arr <= 1e-6)),
        "sss_lba_pair": float(w_lba_pair),
    }
    return w_mat, stats


def compute_mmseqs_sss_statistics(
    sequences: Dict[str, str],
    taxa_names: List[str],
    sensitivity: float = 7.5,
    threads: int = 1,
) -> Optional[Dict[str, float]]:
    """
    Computes Sequence Similarity Score (SSS, w_ij) summary statistics using MMseqs2
    following Matsui & Iwasaki (2020):
        w_ij = (S_bits(i, j) + S_bits(j, i)) / (S_bits(i, i) + S_bits(j, j))
    Returns None if MMseqs2 binary is not available in PATH.
    """
    res = compute_mmseqs_similarity_matrix(sequences, taxa_names, sensitivity=sensitivity, threads=threads)
    if res is None:
        return None
    return res[1]


def similarity_to_distance_matrix(
    w_mat: np.ndarray,
    dist_model: str = "linear",
    epsilon: float = 1e-4,
) -> np.ndarray:
    """
    Converts a symmetrical relative sequence similarity matrix w_ij in [0, 1]
    into an evolutionary distance matrix d_ij matching Matsui & Iwasaki (2020) and Exp 18.

    Parameters
    ----------
    w_mat : np.ndarray
        (N, N) symmetrical relative sequence similarity matrix in [0.0, 1.0] (unitless).
    dist_model : str, default="linear"
        Transformation formula:
        - 'linear' / 'mmseqs_linear' / 'sss_linear': d_ij = 1.0 - w_ij (unitless in [0, 1])
        - 'log' / 'mmseqs_log' / 'sss_log':          d_ij = -ln(max(w_ij, epsilon)) (unitless)
    epsilon : float, default=1e-4
        Small positive floor for logarithmic transformation when w_ij == 0 (unitless).

    Returns
    -------
    d_mat : np.ndarray
        (N, N) symmetrical pairwise distance matrix with zero diagonal (unitless).
    """
    n = w_mat.shape[0]
    d_mat = np.zeros((n, n), dtype=np.float64)
    mode = dist_model.lower().strip()
    is_linear = mode in ("linear", "mmseqs_linear", "sss_linear", "sss")

    for i in range(n):
        for j in range(i + 1, n):
            w = float(w_mat[i, j])
            if is_linear:
                d_val = max(0.0, 1.0 - w)
            else:
                eff_w = max(w, epsilon)
                d_val = max(0.0, -math.log(eff_w))
            d_mat[i, j] = d_val
            d_mat[j, i] = d_val

    return d_mat


def compute_all_pairwise_matrices_and_sss(
    sequences: Dict[str, str],
    taxa_names: List[str],
    model: Optional[WAGModel] = None,
    alpha: float = 1.0,
    mode: str = "global",
    max_poisson_dist: float = 3.0,
    max_wag_dist: float = 5.0,
    return_sss_matrices: bool = False,
) -> Union[
    Tuple[np.ndarray, np.ndarray, Dict[str, float]],
    Tuple[np.ndarray, np.ndarray, Dict[str, float], np.ndarray, np.ndarray],
]:
    """
    Computes BLOSUM62+Poisson distance matrix, WAG+Gamma ML distance matrix,
    Sequence Similarity Score (SSS) statistics, and optionally SSS distance matrices (Exp 18).

    SSS (w_ij) is computed from MMseqs2 bit-scores matching Matsui & Iwasaki (2020) and Exp 18,
    or falls back to normalized Smith-Waterman local alignment scores if MMseqs2 is unavailable.

    Parameters
    ----------
    sequences : Dict[str, str]
        Mapping from taxon name to amino acid sequence.
    taxa_names : List[str]
        Ordered list of taxon names (unitless).
    model : Optional[WAGModel], default=None
        Pre-initialized WAGModel instance.
    alpha : float, default=1.0
        Gamma shape parameter (unitless).
    mode : str, default='global'
        Alignment mode ('global' or 'local').
    max_poisson_dist : float, default=3.0
        Maximum cap for Poisson distance (substitutions/site).
    max_wag_dist : float, default=5.0
        Maximum cap for WAG+Gamma ML distance (substitutions/site).
    return_sss_matrices : bool, default=False
        If True, also returns (D_sss_linear, D_sss_log) computed via Exp 18 transformations.

    Returns
    -------
    Union[Tuple[3], Tuple[5]]
        - D_blosum62 : Poisson distance matrix of shape (n_taxa, n_taxa) in substitutions/site.
        - D_wag_gamma : WAG+Gamma ML distance matrix of shape (n_taxa, n_taxa) in substitutions/site.
        - sss_stats : Dictionary with 'sss_mean', 'sss_median', 'sss_std', 'sss_min', 'sss_max',
          'sss_zero_frac', 'psa_sss_mean', and 'seq_identity_mean' (unitless).
        - (Optional if return_sss_matrices=True) D_sss_linear : d_ij = 1 - w_ij (Exp 18 linear SSS distance).
        - (Optional if return_sss_matrices=True) D_sss_log : d_ij = -ln(max(w_ij, eps)) (Exp 18 log SSS distance).
    """
    if model is None or model.alpha != alpha:
        model = WAGModel(alpha=alpha)

    n = len(taxa_names)
    n_pairs = n * (n - 1) // 2
    pair_counts = np.zeros((n_pairs, 400), dtype=np.float64)
    pair_ij: List[Tuple[int, int]] = []
    D_blosum62 = np.zeros((n, n), dtype=np.float64)
    w_mat_psa = np.zeros((n, n), dtype=np.float64)
    np.fill_diagonal(w_mat_psa, 1.0)

    aligner = create_pairwise_aligner(mode=mode)
    # Also use local aligner with BLOSUM62 for non-negative Smith-Waterman similarity score S_local(i, j)
    local_aligner = create_pairwise_aligner(mode="local")
    self_scores = np.array(
        [float(local_aligner.score(sequences[name], sequences[name])) for name in taxa_names],
        dtype=np.float64,
    )

    psa_w_scores: List[float] = []
    identities: List[float] = []
    psa_w_lba: float = 0.0

    p_idx = 0
    for i in range(n):
        s_i = sequences[taxa_names[i]]
        for j in range(i + 1, n):
            s_j = sequences[taxa_names[j]]
            top_aln = aligner.align(s_i, s_j)[0]
            aln_i, aln_j = top_aln[0], top_aln[1]

            counts_flat, d_poisson = _extract_alignment_counts_and_poisson(
                aln_i, aln_j, max_poisson_dist=max_poisson_dist
            )
            pair_counts[p_idx] = counts_flat
            pair_ij.append((i, j))
            D_blosum62[i, j] = d_poisson
            D_blosum62[j, i] = d_poisson

            # Sequence identity = exp(-d_poisson) when d_poisson < max_poisson_dist
            tot_sites = float(np.sum(counts_flat))
            ident_sites = float(np.sum(counts_flat[np.arange(20) * 21]))
            identities.append(ident_sites / tot_sites if tot_sites > 0 else 0.0)

            # Normalized Smith-Waterman similarity score w_ij in [0, 1]
            denom = self_scores[i] + self_scores[j]
            if denom > 0.0:
                loc_score = float(local_aligner.score(s_i, s_j))
                w_ij = min(1.0, max(0.0, (2.0 * loc_score) / denom))
            else:
                w_ij = 0.0
            psa_w_scores.append(w_ij)
            w_mat_psa[i, j] = w_ij
            w_mat_psa[j, i] = w_ij
            if i == 0 and j == n - 1:
                psa_w_lba = w_ij
            p_idx += 1

    ml_dists = _estimate_ml_distances_from_pair_counts(pair_counts, model, max_dist=max_wag_dist)

    D_wag_gamma = np.zeros((n, n), dtype=np.float64)
    for (i, j), d_val in zip(pair_ij, ml_dists):
        D_wag_gamma[i, j] = d_val
        D_wag_gamma[j, i] = d_val

    psa_arr = np.array(psa_w_scores, dtype=np.float64) if psa_w_scores else np.array([0.0])
    ident_arr = np.array(identities, dtype=np.float64) if identities else np.array([0.0])

    # Check if MMseqs2 is available for exact MMseqs2 SSS; otherwise use Smith-Waterman PSA SSS
    mmseqs_res = compute_mmseqs_similarity_matrix(sequences, taxa_names)
    if mmseqs_res is not None:
        w_mat, sss_stats = mmseqs_res
    else:
        w_mat = w_mat_psa
        sss_stats = {
            "sss_mean": float(np.mean(psa_arr)),
            "sss_median": float(np.median(psa_arr)),
            "sss_std": float(np.std(psa_arr)),
            "sss_min": float(np.min(psa_arr)),
            "sss_max": float(np.max(psa_arr)),
            "sss_zero_frac": float(np.mean(psa_arr <= 1e-6)),
            "sss_lba_pair": float(psa_w_lba) if n >= 2 else 0.0,
        }

    sss_stats["psa_sss_mean"] = float(np.mean(psa_arr))
    sss_stats["psa_sss_median"] = float(np.median(psa_arr))
    sss_stats["psa_sss_std"] = float(np.std(psa_arr))
    sss_stats["psa_sss_min"] = float(np.min(psa_arr))
    sss_stats["psa_sss_max"] = float(np.max(psa_arr))
    sss_stats["psa_sss_lba_pair"] = float(psa_w_lba) if n >= 2 else 0.0
    sss_stats["seq_identity_mean"] = float(np.mean(ident_arr))

    D_sss_linear = similarity_to_distance_matrix(w_mat, dist_model="linear")
    D_sss_log = similarity_to_distance_matrix(w_mat, dist_model="log", epsilon=1e-4)

    if return_sss_matrices:
        return D_blosum62, D_wag_gamma, sss_stats, D_sss_linear, D_sss_log

    return D_blosum62, D_wag_gamma, sss_stats



