"""Visualize N=32 NELSI ultrametric, autocorrelated, and uncorrelated trees.

Generates 300 DPI rectangular phylogram figures:
1. The underlying ultrametric Yule chronogram (before rate variation).
2. Autocorrelated relaxed-clock trees (Kishino model) across increasing v.
3. Uncorrelated lognormal relaxed-clock trees across increasing sd_log.

All relaxed-clock trees color each branch by its local evolutionary rate
(branch_length / chronogram_branch_length, in substitutions/site/unit time)
on a log2 diverging color scale centered at rate = 1.0.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

from src.tree_simulator import SimTreeNode, generate_nelsi_tree, parse_newick_to_tree


def _find_rscript() -> str:
    """Locate the Rscript executable with ape/NELSI installed.

    Returns
    -------
    str
        Path to Rscript executable.
    """
    candidates = [
        os.environ.get("NELSI_RSCRIPT"),
        os.path.join(
            os.environ.get("MAMBA_ROOT_PREFIX", os.path.expanduser("~/.micromamba")),
            "envs",
            "nelsi_env",
            "bin",
            "Rscript",
        ),
        shutil.which("Rscript"),
    ]
    rscript_bin = next(
        (p for p in candidates if p and os.path.isfile(p) and os.access(p, os.X_OK)),
        None,
    )
    if rscript_bin is None:
        raise RuntimeError("Rscript executable not found.")
    return rscript_bin


def generate_ultrametric_chronogram(n_taxa: int = 32, seed: int = 42) -> SimTreeNode:
    """Generate the exact normalized ultrametric Yule chronogram (root height = 1.0).

    Parameters
    ----------
    n_taxa : int, default=32
        Number of terminal taxa N (unitless).
    seed : int, default=42
        Random seed matching ``generate_tree_nelsi.R`` (unitless).

    Returns
    -------
    SimTreeNode
        Root node of the ultrametric time tree (unit time, root height = 1.0).
    """
    rscript_bin = _find_rscript()
    fd, target_path = tempfile.mkstemp(suffix=".nwk")
    os.close(fd)
    r_code = (
        "suppressPackageStartupMessages(library(ape));"
        f"set.seed({int(seed)});"
        f"tr <- ape::rphylo(n={int(n_taxa)}, birth=1.0, death=0.0);"
        "tr <- ape::reorder.phylo(tr, order='cladewise');"
        f"tr$tip.label <- paste0('T', seq_len({int(n_taxa)}));"
        "tr$node.label <- NULL;"
        "rh <- max(ape::node.depth.edgelength(tr));"
        "tr$edge.length <- tr$edge.length / rh;"
        f"ape::write.tree(tr, file='{target_path}');"
    )
    try:
        proc = subprocess.run([rscript_bin, "-e", r_code], capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"Failed to generate chronogram: {proc.stderr}")
        with open(target_path, "r", encoding="utf-8") as f:
            newick_str = f.read().strip()
        return parse_newick_to_tree(newick_str)
    finally:
        if os.path.exists(target_path):
            os.remove(target_path)


def collect_clade_lengths(root: SimTreeNode) -> Dict[FrozenSet[str], float]:
    """Map each node's descendant leaf set to its incoming branch length.

    Parameters
    ----------
    root : SimTreeNode
        Root node of the tree.

    Returns
    -------
    Dict[FrozenSet[str], float]
        Mapping from clade leaf frozenset to incoming branch length.
    """
    mapping: Dict[FrozenSet[str], float] = {}

    def _walk(node: SimTreeNode) -> FrozenSet[str]:
        if node.is_leaf():
            clade = frozenset([node.name])
            mapping[clade] = node.length
            return clade
        leaves: set[str] = set()
        if node.left is not None:
            leaves.update(_walk(node.left))
        if node.right is not None:
            leaves.update(_walk(node.right))
        clade = frozenset(leaves)
        mapping[clade] = node.length
        return clade

    _walk(root)
    return mapping


def compute_branch_rates(
    phylogram: SimTreeNode,
    chronogram: SimTreeNode,
) -> Dict[int, float]:
    """Compute local evolutionary rate (phylogram_length / chronogram_length) per node.

    Parameters
    ----------
    phylogram : SimTreeNode
        Rate-modulated phylogenetic tree (substitutions/site).
    chronogram : SimTreeNode
        Underlying ultrametric time tree (unit time).

    Returns
    -------
    Dict[int, float]
        Mapping from ``id(node)`` in ``phylogram`` to evolutionary rate
        (substitutions/site/unit time).
    """
    chrono_lengths = collect_clade_lengths(chronogram)
    rates: Dict[int, float] = {}

    def _walk(node: SimTreeNode) -> FrozenSet[str]:
        if node.is_leaf():
            clade = frozenset([node.name])
        else:
            leaves: set[str] = set()
            if node.left is not None:
                leaves.update(_walk(node.left))
            if node.right is not None:
                leaves.update(_walk(node.right))
            clade = frozenset(leaves)
        t_len = chrono_lengths.get(clade, 1.0)
        rates[id(node)] = node.length / max(t_len, 1e-12)
        return clade

    _walk(phylogram)
    return rates


def compute_layout_coordinates(
    root: SimTreeNode,
) -> Tuple[Dict[int, float], Dict[int, float], List[SimTreeNode]]:
    """Compute horizontal (x) and vertical (y) layout coordinates.

    Parameters
    ----------
    root : SimTreeNode
        Root node of the tree.

    Returns
    -------
    Tuple[Dict[int, float], Dict[int, float], List[SimTreeNode]]
        - ``x_coords``: Mapping from ``id(node)`` to root-to-node distance.
        - ``y_coords``: Mapping from ``id(node)`` to vertical coordinate.
        - ``leaves``: Ordered list of terminal leaf nodes.
    """
    x_coords: Dict[int, float] = {}
    y_coords: Dict[int, float] = {}
    leaves: List[SimTreeNode] = []

    def _walk(node: SimTreeNode, current_x: float) -> None:
        nid = id(node)
        x_coords[nid] = current_x
        if node.is_leaf():
            y_coords[nid] = float(len(leaves))
            leaves.append(node)
            return
        if node.left is not None:
            _walk(node.left, current_x + node.left.length)
        if node.right is not None:
            _walk(node.right, current_x + node.right.length)
        child_ys = [y_coords[id(c)] for c in (node.left, node.right) if c is not None]
        y_coords[nid] = float(np.mean(child_ys))

    _walk(root, 0.0)
    return x_coords, y_coords, leaves


def plot_phylogram_by_rate(
    root: SimTreeNode,
    branch_rates: Optional[Dict[int, float]],
    title: str,
    xlabel: str,
    out_path: Path,
    ref_x: Optional[float] = 1.0,
    log2_rate_lim: float = 3.0,
) -> Dict[str, float]:
    """Render a 300 DPI rectangular phylogram with branches colored by evolutionary rate.

    Parameters
    ----------
    root : SimTreeNode
        Root node of the tree to plot.
    branch_rates : Optional[Dict[int, float]]
        Mapping from ``id(node)`` to branch evolutionary rate (substitutions/site/unit time).
        If ``None``, draws uniform dark-slate branches (used for the ultrametric chronogram).
    title : str
        Figure title.
    xlabel : str
        X-axis label with units.
    out_path : Path
        Output PNG file path.
    ref_x : Optional[float], default=1.0
        Optional vertical reference line (e.g., ``D = 1.0``).
    log2_rate_lim : float, default=3.0
        Symmetric ``log2(rate)`` color scale limit (e.g., ``3.0`` spans ``1/8x`` to ``8x``).

    Returns
    -------
    Dict[str, float]
        Root-to-tip distance summary statistics (mean, std, min, max, max_min_ratio).
    """
    x_coords, y_coords, leaves = compute_layout_coordinates(root)
    leaf_dists = np.array([x_coords[id(lf)] for lf in leaves], dtype=float)
    min_d = float(np.min(leaf_dists))
    max_d = float(np.max(leaf_dists))
    mean_d = float(np.mean(leaf_dists))
    std_d = float(np.std(leaf_dists))
    ratio_d = max_d / min_d if min_d > 0 else float("inf")

    fig, ax = plt.subplots(figsize=(10.2, 7.4), dpi=300)
    cmap = plt.get_cmap("coolwarm")
    norm = mcolors.TwoSlopeNorm(vmin=-log2_rate_lim, vcenter=0.0, vmax=log2_rate_lim)

    def _color_for_node(node: SimTreeNode) -> object:
        if branch_rates is None:
            return "#2c3e50"
        r = max(branch_rates.get(id(node), 1.0), 1e-6)
        val = float(np.clip(np.log2(r), -log2_rate_lim, log2_rate_lim))
        return cmap(norm(val))

    def _draw_node(node: SimTreeNode) -> None:
        nid = id(node)
        x0 = x_coords[nid]
        children = [c for c in (node.left, node.right) if c is not None]
        if children:
            ys = [y_coords[id(c)] for c in children]
            ax.plot([x0, x0], [min(ys), max(ys)], color="#444444", linewidth=1.4, zorder=2)
            for child in children:
                cid = id(child)
                x1 = x_coords[cid]
                yc = y_coords[cid]
                ax.plot(
                    [x0, x1],
                    [yc, yc],
                    color=_color_for_node(child),
                    linewidth=2.2,
                    zorder=3,
                )
                _draw_node(child)

    _draw_node(root)

    x_offset = max_d * 0.012
    for lf in leaves:
        lid = id(lf)
        lx = x_coords[lid]
        ly = y_coords[lid]
        ax.scatter(
            [lx],
            [ly],
            color=_color_for_node(lf),
            edgecolor="#222222",
            s=28,
            linewidths=0.6,
            zorder=4,
        )
        ax.text(
            lx + x_offset,
            ly,
            f"{lf.name} ({lx:.2f})",
            va="center",
            ha="left",
            fontsize=8.0,
            color="#222222",
        )

    if ref_x is not None:
        ax.axvline(
            ref_x,
            color="#2ca02c",
            linestyle="--",
            linewidth=1.3,
            label=f"Strict Clock / Mean D = {ref_x:.1f}",
            zorder=1,
        )
        ax.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=8.8)

    if branch_rates is not None:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.032, aspect=28)
        ticks = [-3, -2, -1, 0, 1, 2, 3]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(["1/8x", "1/4x", "1/2x", "1x", "2x", "4x", "8x"])
        cbar.set_label("Branch Evolutionary Rate (relative to 1x)", fontsize=9.5, fontweight="bold")

    ax.set_xlim(-max_d * 0.02, max_d * 1.18)
    ax.set_ylim(-2.2, len(leaves) + 0.6)
    ax.set_yticks([])
    ax.set_xlabel(xlabel, fontsize=10.5, fontweight="bold")
    ax.set_title(title, fontsize=11.0, fontweight="bold", pad=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", linestyle=":", alpha=0.35)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)

    return {
        "mean": mean_d,
        "std": std_d,
        "min": min_d,
        "max": max_d,
        "max_min_ratio": ratio_d,
    }


def main() -> None:
    """Generate ultrametric, autocorrelated, and uncorrelated N=32 tree plots.

    Parameters
    ----------
    None
        Parses CLI options and writes PNG figures.
    """
    parser = argparse.ArgumentParser(description="Plot N=32 ultrametric, correlated, and uncorrelated NELSI trees.")
    parser.add_argument("--taxa", type=int, default=32, help="Number of taxa (default: 32)")
    parser.add_argument("--divergence", type=float, default=1.0, help="Mean root-to-tip distance D (default: 1.0)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for fixed topology (default: 42)")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/nelsi_rate_trees"),
        help="Output directory for tree PNG files",
    )
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        default=None,
        help="Optional artifact directory to copy PNGs for inline display",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.artifact_dir is not None:
        args.artifact_dir.mkdir(parents=True, exist_ok=True)

    # 1. Ultrametric chronogram before rate variation
    chronogram = generate_ultrametric_chronogram(n_taxa=args.taxa, seed=args.seed)
    chrono_fname = f"nelsi_ultrametric_N{args.taxa}.png"
    chrono_path = args.outdir / chrono_fname
    plot_phylogram_by_rate(
        root=chronogram,
        branch_rates=None,
        title=(
            f"Ultrametric Yule Chronogram Before Rate Variation (N = {args.taxa}, Seed = {args.seed})\n"
            "All Root-to-Tip Time Heights = 1.00 (Strict Molecular Clock)"
        ),
        xlabel="Relative Time from Root (unit time, height = 1.0)",
        out_path=chrono_path,
        ref_x=1.0,
    )
    if args.artifact_dir is not None:
        shutil.copy2(chrono_path, args.artifact_dir / chrono_fname)
    print(f"[OK] Ultrametric -> {chrono_path}")

    # 2. Autocorrelated Kishino trees colored by branch evolutionary rate
    for v_val in [0.3, 0.6, 1.0, 1.5, 2.0, 3.0]:
        root, _ = generate_nelsi_tree(
            n_taxa=args.taxa,
            divergence=args.divergence,
            seed=args.seed,
            model="nelsi_kishino",
            v=v_val,
            scale_mode="root_to_tip",
        )
        rates = compute_branch_rates(root, chronogram)
        fname = f"nelsi_autocor_rate_N{args.taxa}_v{v_val:.1f}.png"
        out_path = args.outdir / fname
        title = (
            f"NELSI Autocorrelated Clock (Kishino, N = {args.taxa}, v = {v_val:g}) — Colored by Branch Rate\n"
            f"Root-to-Tip Mean = {args.divergence:.2f}"
        )
        stats = plot_phylogram_by_rate(
            root=root,
            branch_rates=rates,
            title=title,
            xlabel="Evolutionary Distance from Root (substitutions/site)",
            out_path=out_path,
            ref_x=args.divergence,
        )
        if args.artifact_dir is not None:
            shutil.copy2(out_path, args.artifact_dir / fname)
        print(
            f"[OK] Autocor v={v_val:.1f} -> {out_path} | "
            f"Max/Min={stats['max_min_ratio']:.2f}x, SD={stats['std']:.4f}"
        )

    # 3. Uncorrelated Lognormal trees colored by branch evolutionary rate
    for sd_val in [0.3, 0.5, 1.0, 1.5, 2.0]:
        root, _ = generate_nelsi_tree(
            n_taxa=args.taxa,
            divergence=args.divergence,
            seed=args.seed,
            model="nelsi_uncor_lnorm",
            sd_log=sd_val,
            scale_mode="root_to_tip",
        )
        rates = compute_branch_rates(root, chronogram)
        fname = f"nelsi_uncor_rate_N{args.taxa}_sd_log{sd_val:.1f}.png"
        out_path = args.outdir / fname
        title = (
            f"NELSI Uncorrelated Lognormal Clock (N = {args.taxa}, sd_log = {sd_val:g}) — Colored by Branch Rate\n"
            f"Root-to-Tip Mean = {args.divergence:.2f}"
        )
        stats = plot_phylogram_by_rate(
            root=root,
            branch_rates=rates,
            title=title,
            xlabel="Evolutionary Distance from Root (substitutions/site)",
            out_path=out_path,
            ref_x=args.divergence,
        )
        if args.artifact_dir is not None:
            shutil.copy2(out_path, args.artifact_dir / fname)
        print(
            f"[OK] Uncor sd_log={sd_val:.1f} -> {out_path} | "
            f"Max/Min={stats['max_min_ratio']:.2f}x, SD={stats['std']:.4f}"
        )


if __name__ == "__main__":
    main()
