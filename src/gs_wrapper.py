"""
Wrapper module for Graph Splitting (gs2) phylogenetic tree reconstruction.

Supports running `gs2` directly from a pairwise distance matrix (`-d` mode) or
a sequence similarity matrix (`-w` mode).
"""

from typing import List, Optional
import os
import shutil
import tempfile
import subprocess
import numpy as np


def find_gs_binary() -> str:
    """
    Locates the `gs2` executable on the system.

    Returns
    -------
    str
        Absolute path to the `gs2` executable.
    """
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    candidates = [
        shutil.which("gs2"),
        os.path.join(repo_root, "bin", "gs2"),
        os.path.abspath(os.path.join(repo_root, "..", "PhyloTree_Methods_Comparison", "bin", "gs2")),
    ]
    for path in candidates:
        if path and os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    raise FileNotFoundError(
        "Graph Splitting binary `gs2` not found in PATH, bin/gs2, or ../PhyloTree_Methods_Comparison/bin/gs2."
    )


def _format_phylip_matrix(
    M: np.ndarray,
    taxa_names: List[str],
) -> str:
    """
    Formats a symmetric square matrix into PHYLIP-like text input for `gs2`.

    Parameters
    ----------
    M : np.ndarray
        Symmetric matrix of shape (n_taxa, n_taxa) [substitutions/site or unitless].
    taxa_names : List[str]
        Ordered list of taxon names [unitless].

    Returns
    -------
    str
        Formatted matrix string.
    """
    n = len(taxa_names)
    lines = [f"{n}"]
    for i in range(n):
        vals = " ".join(f"{float(M[i, j]):.8f}" for j in range(n))
        lines.append(f"{taxa_names[i]} {vals}")
    return "\n".join(lines) + "\n"


def run_gs_from_distance_matrix(
    D: np.ndarray,
    taxa_names: List[str],
    gs_bin: Optional[str] = None,
    threads: int = 1,
) -> str:
    """
    Reconstructs a phylogenetic tree topology using Graph Splitting (`gs2 -d`)
    from an observed pairwise distance matrix.

    Parameters
    ----------
    D : np.ndarray
        Symmetric distance matrix of shape (n_taxa, n_taxa) in substitutions/site.
    taxa_names : List[str]
        Ordered list of taxon labels [unitless].
    gs_bin : Optional[str], default=None
        Path to `gs2` binary. If None, auto-detected via `find_gs_binary()`.
    threads : int, default=1
        Number of CPU threads for `gs2` [unitless].

    Returns
    -------
    str
        Reconstructed Newick tree string (topology only, ending with ';') [unitless].
    """
    if gs_bin is None:
        gs_bin = find_gs_binary()

    mat_str = _format_phylip_matrix(D, taxa_names)
    with tempfile.TemporaryDirectory(prefix="gs_dist_") as tmpdir:
        mat_path = os.path.join(tmpdir, "dist.mat")
        with open(mat_path, "w", encoding="utf-8") as f:
            f.write(mat_str)

        cmd = [gs_bin, "-s", "-l", "-t", str(threads), "-d", mat_path]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"gs2 execution failed (code {res.returncode}):\n{res.stderr}")

        nwk = res.stdout.strip().replace('"', "")
        if not nwk or "(" not in nwk:
            raise ValueError(f"Invalid Newick output from gs2: {nwk}\nStderr: {res.stderr}")
        if not nwk.endswith(";"):
            nwk += ";"
        return nwk
