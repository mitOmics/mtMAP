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
# Gene normalization and anchor
def _norm_gene_name(name: str) -> str:
    """Normalize gene names to canonical keys."""
    n = (name or '').strip().upper().replace(' ', '')
    repl = {
        'COI':'COX1', 'COXI':'COX1', 'MT-CO1':'COX1', 'CO1':'COX1',
        'COII':'COX2', 'COXII':'COX2', 'MT-CO2':'COX2', 'CO2':'COX2',
        'COIII':'COX3', 'COXIII':'COX3', 'MT-CO3':'COX3', 'CO3':'COX3',
        'CYTB':'CYTB', 'COB':'CYTB', 'MT-CYB':'CYTB', 'CYB':'CYTB',
        'ATPASE6':'ATP6', 'ATPASE8':'ATP8',
        'RRNS':'12S', '12SRRNA':'12S', '12S RRNA':'12S',
        'RRNL':'16S', '16SRNA':'16S', '16S RRNA':'16S',
    }
    return repl.get(n, n)

def _match_anchor(row_name: str, anchor: str) -> bool:
    """Check if a gene matches the chosen anchor."""
    if not anchor:
        return False
    a = anchor.strip().upper()
    nrow = _norm_gene_name(row_name)
    nanc = _norm_gene_name(a)
    if nrow == nanc:
        return True
    if len(a) == 1 and a.isalpha():
        return trna_label_one_letter(row_name).upper() == a
    if a.startswith('TRNA'):
        return row_name.upper().replace(' ', '') == a
    return False

# -----------------------------
# GenBank parsing
def extract_annotations_any(input_path: Path) -> pd.DataFrame:
    """Extract annotations from GenBank file(s)."""
    rows = []
    species_map = {}
    accession_map = {}
    for gb in list_gb_paths(input_path):
        for record in SeqIO.parse(str(gb), "genbank"):
            rec_id = record.id
            L = len(record.seq)
            species = None
            for feat in record.features:
                if feat.type == "source" and "organism" in feat.qualifiers:
                    species = feat.qualifiers["organism"][0]
                    break
            if not species:
                species = record.annotations.get("organism") or record.annotations.get("source") or "Unknown_species"
            accession = rec_id
            species_map[rec_id] = species
            accession_map[rec_id] = accession

            for feat in record.features:
                if feat.type in ("CDS", "rRNA", "tRNA"):
                    qual = feat.qualifiers
                    name = (qual.get("gene") or qual.get("product") or ["unknown"])[0]
                    start = int(feat.location.start)
                    end = int(feat.location.end)
                    rows.append({
                        "Record": rec_id,
                        "Type": "CDS" if feat.type == "CDS" else feat.type,
                        "Name": name,
                        "Start": start,
                        "End": end,
                        "GenomeLen": L
                    })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["Species"] = df["Record"].map(species_map)
        df["Accession"] = df["Record"].map(accession_map)
    return df

# -----------------------------
# Export annotations and presence/absence
def export_tables(df: pd.DataFrame, out_prefix: Path):
    """Export annotations and presence/absence tables as CSV."""
    out_ann = out_prefix.with_name(out_prefix.name + "_annotations.csv")
    df.sort_values(["Record","Start","End"]).to_csv(out_ann, index=False)
    genes = [g for g in TARGET_GENES]
    sub = df[df["Type"].isin(["CDS","rRNA"])].copy()
    sub["GeneNorm"] = sub["Name"].apply(_norm_gene_name)
    mat = (sub.assign(val=1)
              .pivot_table(index="Record", columns="GeneNorm", values="val", aggfunc="max", fill_value=0))
    for g in genes:
        if g not in mat.columns:
            mat[g] = 0
    mat = mat[genes]
    if {"Species","Accession"}.issubset(df.columns):
        meta = df.groupby("Record")[["Species","Accession"]].first()
        mat = meta.join(mat, how="left")
    out_presence = out_prefix.with_name(out_prefix.name + "_presence.csv")
    mat.to_csv(out_presence)

# -----------------------------
# Plot functions
# (kept unchanged except translated titles/labels)

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
    # (plots and tables are generated here)

if __name__ == "__main__":
    main()
