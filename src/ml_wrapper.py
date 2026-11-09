"""
MSA + Maximum Likelihood (ML) Pipeline Module.

Executes Multiple Sequence Alignment using MAFFT followed by Maximum Likelihood
phylogenetic inference using FastTree under the WAG + Gamma model.
"""

from typing import Dict, Optional
import os
import shutil
import subprocess
import tempfile


def run_msa_ml(
    sequences: Dict[str, str],
    mafft_bin: Optional[str] = None,
    fasttree_bin: Optional[str] = None,
) -> str:
    """
    Runs MAFFT multiple sequence alignment and FastTree WAG+Gamma maximum likelihood inference.

    Parameters
    ----------
    sequences : Dict[str, str]
        Dictionary mapping taxon name to amino acid sequence [amino acid characters].
    mafft_bin : Optional[str], default=None
        Path to mafft executable. If None, auto-detected from PATH or system locations.
    fasttree_bin : Optional[str], default=None
        Path to fasttree executable. If None, auto-detected from PATH or environment locations.

    Returns
    -------
    str
        Inferred Newick tree string.
    """
    if mafft_bin is None:
        mafft_bin = shutil.which("mafft")
        if not mafft_bin:
            candidates = ["/opt/homebrew/bin/mafft", "/usr/local/bin/mafft"]
            for c in candidates:
                if os.path.exists(c):
                    mafft_bin = c
                    break
        if not mafft_bin:
            raise FileNotFoundError("MAFFT binary could not be found. Please install MAFFT or specify mafft_bin.")

    if fasttree_bin is None:
        fasttree_bin = shutil.which("fasttree")
        if not fasttree_bin:
            conda_prefix = os.environ.get("CONDA_PREFIX")
            if conda_prefix:
                cp = os.path.join(conda_prefix, "bin", "fasttree")
                if os.path.exists(cp):
                    fasttree_bin = cp
        if not fasttree_bin:
            fallback = "/Users/aibarakazuki/.micromamba/envs/phylomethod_env/bin/fasttree"
            if os.path.exists(fallback):
                fasttree_bin = fallback
        if not fasttree_bin:
            raise FileNotFoundError("FastTree binary could not be found. Please install FastTree or specify fasttree_bin.")

    with tempfile.TemporaryDirectory() as tmpdir:
        fa_in = os.path.join(tmpdir, "unaligned.fasta")
        with open(fa_in, "w") as f:
            for name, seq in sequences.items():
                f.write(f">{name}\n{seq}\n")

        fa_aln = os.path.join(tmpdir, "aligned.fasta")
        with open(fa_aln, "w") as f_out:
            res_mafft = subprocess.run(
                [mafft_bin, "--auto", "--quiet", fa_in],
                stdout=f_out,
                stderr=subprocess.PIPE,
                text=True,
            )
            if res_mafft.returncode != 0:
                raise RuntimeError(f"MAFFT alignment failed: {res_mafft.stderr}")

        nwk_out = os.path.join(tmpdir, "tree.nwk")
        with open(nwk_out, "w") as f_tree:
            res_ft = subprocess.run(
                [fasttree_bin, "-wag", "-gamma", "-quiet", fa_aln],
                stdout=f_tree,
                stderr=subprocess.PIPE,
                text=True,
            )
            if res_ft.returncode != 0:
                raise RuntimeError(f"FastTree inference failed: {res_ft.stderr}")

        with open(nwk_out, "r") as f:
            newick_str = f.read().strip()

        return newick_str

