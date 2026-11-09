"""
FastME wrapper module for BioNJ, FastME-NNI, and FastME-SPR execution.
"""

from typing import List, Optional
import os
import shutil
import subprocess
import tempfile
import numpy as np


def run_fastme(
    D: np.ndarray,
    taxa_names: List[str],
    method: str = "bionj",
    fastme_bin: Optional[str] = None,
) -> str:
    """
    Executes FastME to construct phylogenetic trees using BioNJ, FastME-NNI, or FastME-SPR.

    Parameters
    ----------
    D : np.ndarray
        Distance matrix of shape (n, n) in substitutions/site.
    taxa_names : List[str]
        Identifiers for taxa (unitless).
    method : str, default='bionj'
        FastME tree inference method: 'bionj' (BioNJ), 'nni' (BME-NNI), or 'spr' (BME-SPR) (unitless).
    fastme_bin : Optional[str], default=None
        Path to fastme executable. If None, auto-detected from environment or fallback path.

    Returns
    -------
    str
        Inferred Newick tree string.
    """
    if fastme_bin is None:
        fastme_bin = shutil.which("fastme")
        if not fastme_bin:
            default_path = "/Users/aibarakazuki/.micromamba/envs/phylomethod_env/bin/fastme"
            if os.path.exists(default_path):
                fastme_bin = default_path
            else:
                raise FileNotFoundError("FastME binary could not be located in PATH or default prefix.")

    n = len(taxa_names)

    with tempfile.TemporaryDirectory() as tmpdir:
        mat_path = os.path.join(tmpdir, "matrix.phy")
        with open(mat_path, "w") as f:
            f.write(f" {n}\n")
            for i, name in enumerate(taxa_names):
                row_str = " ".join(f"{D[i, j]:.8f}" for j in range(n))
                f.write(f"{name:<10} {row_str}\n")

        out_nwk = os.path.join(tmpdir, "tree.nwk")
        cmd = [fastme_bin, "-i", mat_path, "-o", out_nwk]

        if method.lower() == "bionj":
            cmd += ["-m", "B"]
        elif method.lower() == "nni":
            cmd += ["-m", "B", "-n"]
        elif method.lower() == "spr":
            cmd += ["-m", "B", "-s"]
        else:
            raise ValueError(f"Unknown method '{method}'. Choose 'bionj', 'nni', or 'spr'.")

        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"FastME execution failed: {res.stderr}")

        with open(out_nwk, "r") as f:
            newick_str = f.read().strip()

        return newick_str

