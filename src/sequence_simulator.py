"""
Amino acid sequence evolution simulator under WAG + Gamma and pairwise alignment distance estimator.

Implements:
1. WAG + Gamma (discrete 4-category Gamma rate heterogeneity) amino acid sequence evolution.
2. Pairwise sequence alignment (PSA) using Biopython PairwiseAligner with BLOSUM62 substitution matrix.
3. Global alignment (Needleman-Wunsch) and Local alignment (Smith-Waterman).
4. Amino acid Poisson correction distance calculation: d = -ln(1 - p).
"""

from typing import Dict, List, Tuple, Optional
import os
import math
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


def simulate_wag_gamma_sequences(
    root: SimTreeNode,
    seq_length: int = 300,
    alpha: float = 1.0,
    model: Optional[WAGModel] = None,
    seed: Optional[int] = None,
) -> Dict[str, str]:
    """
    Simulates amino acid sequences along a tree under WAG + Gamma.

    Parameters
    ----------
    root : SimTreeNode
        Root node of the phylogenetic tree.
    seq_length : int, default=300
        Length of amino acid alignment L (residues).
    alpha : float, default=1.0
        Shape parameter of Gamma rate distribution.
    model : Optional[WAGModel], default=None
        Pre-initialized WAGModel. If None, initialized with alpha.
    seed : Optional[int], default=None
        Random seed.

    Returns
    -------
    Dict[str, str]
        Mapping from taxon leaf name to amino acid sequence string.
    """
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

