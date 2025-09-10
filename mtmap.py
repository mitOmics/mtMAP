#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mtmap.py - Mitochondrial Map (English version)
"""

from pathlib import Path
import argparse
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import SeqIO
import importlib.util, sys

def _fallback_color(key: str) -> str:
    """Generate a deterministic pastel color for keys outside TARGET_GENES."""
    import hashlib
    h = hashlib.md5(key.encode("utf-8")).hexdigest()
    r = int(h[0:2], 16) / 255.0
    g = int(h[2:4], 16) / 255.0
    b = int(h[4:6], 16) / 255.0
    r = r + (1.0 - r) * 0.6
    g = g + (1.0 - g) * 0.6
    b = b + (1.0 - b) * 0.6
    return '#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255))

# -----------------------------
# Appearance
DEFAULT_OTHER = "#7f7f7f"
T_RNA_COLOR = "#fff9c4"

TARGET_GENES = [
    "COX1","COX2","COX3",
    "ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
    "CYTB","ATP6","ATP8","12S","16S"
]

def _gene_color_map():
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
    if input_path.is_dir():
        files = sorted([p for p in input_path.iterdir() if p.suffix.lower() in GB_EXTS])
        if not files:
            raise FileNotFoundError(f"No GB/GBK file found in: {input_path}")
        return files
    if input_path.is_file():
        return [input_path]
    raise FileNotFoundError(f"Path not found: {input_path}")

# -----------------------------
# (Remaining functions omitted for brevity in this snippet)
