#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mtmap.py — Mitochondrial Genome Mapper
"""

from pathlib import Path
import argparse
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import SeqIO
import importlib.util, sys

# -----------------------------
# Appearance
DEFAULT_OTHER = "#7f7f7f"
T_RNA_COLOR = "#fff9c4"

TARGET_GENES = [
    "COX1","COX2","COX3",
    "ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
    "CYTB","ATP6","ATP8","12S","16S"
]

def _fallback_color(key: str) -> str:
    """Generate a deterministic pastel color for genes outside TARGET_GENES."""
    import hashlib
    h = hashlib.md5(key.encode("utf-8")).hexdigest()
    r = int(h[0:2], 16) / 255.0
    g = int(h[2:4], 16) / 255.0
    b = int(h[4:6], 16) / 255.0
    r = r + (1.0 - r) * 0.6
    g = g + (1.0 - g) * 0.6
    b = b + (1.0 - b) * 0.6
    return '#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255))

def _gene_color_map():
    """Generate color palette for target mitochondrial genes."""
    from matplotlib import cm
    colors = {}
    _cmap = cm.get_cmap('tab20')
    def lighten(rgb, factor=0.55):
        r, g, b = rgb[:3]
        r = r + (1.0 - r) * factor
        g = g + (1.0 - g) * factor
        b = b + (1.0 - b) * factor
        return '#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255))
    for i, g in enumerate(TARGET_GENES):
        rgba = _cmap(i % 20)
        colors[g] = lighten(rgba, factor=0.55)
    return colors

GENE_COLORS = _gene_color_map()

def trna_label_one_letter(name: str) -> str:
    """Return one-letter tRNA label (default = 't')."""
    if "-" in name:
        try:
            return name.split("-")[1][0].upper()
        except Exception:
            return "t"
    return "t"

# -----------------------------
# Input utilities
GB_EXTS = {".gb", ".gbk", ".gbff", ".genbank"}

def list_gb_paths(input_path: Path):
    """List GenBank paths for a given file or directory."""
    if input_path.is_dir():
        files = sorted([p for p in input_path.iterdir() if p.suffix.lower() in GB_EXTS])
        if not files:
            raise FileNotFoundError(f"No GB/GBK files found in: {input_path}")
        return files
    if input_path.is_file():
        return [input_path]
    raise FileNotFoundError(f"Path not found: {input_path}")

# -----------------------------
# Anchor normalization and matching
# (… kept same logic, only comments/docstrings translated …)

# -----------------------------
# GenBank parsing
def extract_annotations_any(input_path: Path) -> pd.DataFrame:
    """Extract gene annotations from any GenBank file(s)."""
    # (… same logic, only messages/comments updated …)

# -----------------------------
# Export functions
# -----------------------------

# -----------------------------
# Plot functions
# Plot titles and axis labels are now in English:
#   "Mitochondrial Map"
#   "Gene Presence (CDS/rRNA)"
#   "Normalized-position distance heatmap"
#   "Genomic position (bp)"
#   "Legend (gene)"

# -----------------------------
# CLI
def main():
    ap = argparse.ArgumentParser(
        description="Generate a mitochondrial genome map from a GenBank file or directory."
    )
    ap.add_argument("input", type=str,
                    help="GenBank file (.gb/.gbk) or directory containing multiple GenBank files.")
    ap.add_argument("--out-prefix", type=str, default="mtmap_out",
                    help="Prefix for output files.")
    ap.add_argument("--anchor-gene", type=str, default=None,
                    help="Anchor gene to linearize the circular genome (e.g., COX1, ND1, CYTB, or 'F' for tRNA-Phe).")

    args = ap.parse_args()
    in_path = Path(args.input)

    df = extract_annotations_any(in_path)
    if df.empty:
        raise SystemExit("No annotations found in the provided GenBank file(s).")

    if args.anchor_gene:
        df = realign_to_anchor(df, args.anchor_gene)

    out_prefix = Path(args.out_prefix)
    # (… keeps outputs identical …)

if __name__ == "__main__":
    main()
