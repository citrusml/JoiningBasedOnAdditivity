"""
Unit test for MSA + ML pipeline module.
"""

from typing import Dict
import pytest
import dendropy
from src.ml_wrapper import run_msa_ml


def test_run_msa_ml() -> None:
    """Verifies that MAFFT + FastTree executes cleanly and outputs valid Newick."""
    seqs: Dict[str, str] = {
        "T1": "MKTIIALSYIFCLVFA",
        "T2": "MKTIIALSYIFCLVFA",
        "T3": "MKTIIALSNIFCLVFA",
        "T4": "MKTIIALSNIFCLVFA",
    }
    nwk = run_msa_ml(seqs)
    assert len(nwk) > 0
    tns = dendropy.TaxonNamespace(list(seqs.keys()))
    tree = dendropy.Tree.get(data=nwk, schema="newick", taxon_namespace=tns)
    assert len(tree.leaf_nodes()) == 4

