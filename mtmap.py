#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mtmap.py — Mitochondrial Map (English version)
This file is a faithful translation of comments, docstrings, and user-facing messages.
Functionality, function/variable names, and logic are preserved.
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
    """Return a one-letter label for a tRNA feature name (e.g., 'tRNA-Phe' -> 'F')."""
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
    """Return a list of GenBank file paths from a directory or a single file path."""
    if input_path.is_dir():
        files = sorted([p for p in input_path.iterdir() if p.suffix.lower() in GB_EXTS])
        if not files:
            raise FileNotFoundError(f"No GB/GBK files found in: {input_path}")
        return files
    if input_path.is_file():
        return [input_path]
    raise FileNotFoundError(f"Path not found: {input_path}")

# -----------------------------
# Normalization and anchor
def _norm_gene_name(name: str) -> str:
    """Normalize common gene synonyms to canonical keys (e.g., COI -> COX1)."""
    n = (name or '').strip().upper().replace(' ', '')
    repl = {
        'COI':'COX1', 'COXI':'COX1', 'MT-CO1':'COX1', 'CO1':'COX1',
        'COII':'COX2', 'COXII':'COX2', 'MT-CO2':'COX2', 'CO2':'COX2',
        'COIII':'COX3', 'COXIII':'COX3', 'MT-CO3':'COX3', 'CO3':'COX3',
        'CYTB':'CYTB', 'COB':'CYTB', 'MT-CYB':'CYTB', 'CYB':'CYTB',
        'ATPASE6':'ATP6', 'ATPASE8':'ATP8',
        'RRNS':'12S', '12SRRNA':'12S', '12S RRNA':'12S',
        'RRNL':'16S', '16SRRNA':'16S', '16S RRNA':'16S',
    }
    return repl.get(n, n)

def _match_anchor(row_name: str, anchor: str) -> bool:
    """Return True if row_name matches the user-provided anchor gene spec."""
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


def _normalize_gene_key(name: str) -> str:
    """
    Normalize gene names to canonical keys:
    - Case-insensitive, remove spaces and 'MT-' prefix
    - Common synonyms: COI->COX1, COB->CYTB, etc.
    - Convert NADx/nadx/nadhx -> NDx (e.g., NAD1 -> ND1)
    - Normalize rRNAs: 12S, 16S
    - Heuristics from 'product' (e.g., 'NADH dehydrogenase subunit 1' -> ND1)
    """
    if not name:
        return ""

    s = name.upper().strip()
    s = s.replace(" ", "").replace("MT-", "")

    # rRNA patterns
    if "12S" in s or "RRNS" in s:
        return "12S"
    if "16S" in s or "RRNL" in s or "16SRRNA" in s or "16SRNA" in s:
        return "16S"

    # COX synonyms
    s_cox = s.replace("COI", "COX1").replace("CO1", "COX1") \
             .replace("COII", "COX2").replace("CO2", "COX2") \
             .replace("COIII", "COX3").replace("CO3", "COX3")
    if s_cox.startswith("COX"):
        s = s_cox

    # CYTB synonyms
    if s in {"CYB","COB","CYTOCHROMEB","MTCYB","MT-CYB"}:
        return "CYTB"

    # ATP6/8 synonyms
    if s in {"ATPASE6", "ATP6"}:
        return "ATP6"
    if s in {"ATPASE8", "ATP8"}:
        return "ATP8"

    # NADx -> NDx
    m = re.match(r"^(NADH?|ND)(\d)(L?)$", s)
    if m:
        base = "ND" + m.group(2)
        if m.group(3) == "L":
            base += "L"
        return base

    # Heuristics based on 'product' (when Name comes from 'product')
    sp = name.upper()
    if "CYTOCHROME B" in sp:
        return "CYTB"
    if "CYTOCHROME C OXIDASE" in sp:
        if " SUBUNIT I" in sp or " SUBUNIT 1" in sp:
            return "COX1"
        if " SUBUNIT II" in sp or " SUBUNIT 2" in sp:
            return "COX2"
        if " SUBUNIT III" in sp or " SUBUNIT 3" in sp:
            return "COX3"
    if "ATP SYNTHASE" in sp and " SUBUNIT 6" in sp:
        return "ATP6"
    if "ATP SYNTHASE" in sp and " SUBUNIT 8" in sp:
        return "ATP8"
    if "NADH DEHYDROGENASE SUBUNIT" in sp:
        mm = re.search(r"SUBUNIT\s+(\d)"," "+sp)
        if mm:
            nd = "ND" + mm.group(1)
            if " ND4L" in sp or " SUBUNIT 4L" in sp:
                nd = "ND4L"
            return nd
    if "SMALL SUBUNIT RIBOSOMAL RNA" in sp or "12S RIBOSOMAL RNA" in sp:
        return "12S"
    if "LARGE SUBUNIT RIBOSOMAL RNA" in sp or "16S RIBOSOMAL RNA" in sp:
        return "16S"

    return s

def realign_to_anchor(df: pd.DataFrame, anchor: str) -> pd.DataFrame:
    """Realign features so the anchor gene starts at position 0 for each record."""
    if not anchor or df.empty:
        return df
    out = []
    for rec_id, sub in df.groupby('Record', sort=False):
        L = int(sub['GenomeLen'].iloc[0])
        candidates = sub[sub['Name'].apply(lambda n: _match_anchor(n, anchor))]
        if candidates.empty:
            out.append(sub)
            continue
        offset = int(candidates['Start'].min())
        s = sub.copy()
        ns = (s['Start'] - offset) % L
        ne = (s['End'] - offset) % L
        wrap = ne <= ns
        s.loc[:, 'Start'] = ns.astype(int)
        s.loc[:, 'End'] = (ne + wrap * L).astype(int)
        out.append(s)
    return pd.concat(out, ignore_index=True)

# -----------------------------
# GenBank parsing
def extract_annotations_any(input_path: Path) -> pd.DataFrame:
    """Parse GenBank entries (CDS, rRNA, tRNA) and build a tidy annotations table."""
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
# Table exports
def export_tables(df: pd.DataFrame, out_prefix: Path):
    """Export full annotations and presence/absence matrices as CSV files."""
    out_ann = out_prefix.with_name(out_prefix.name + "_annotations.csv")
    df.sort_values(["Record","Start","End"]).to_csv(out_ann, index=False)
    genes = [g for g in TARGET_GENES]
    sub = df[df["Type"].isin(["CDS","rRNA"])].copy()
    sub["GeneNorm"] = sub["Name"].apply(_normalize_gene_key)
    syn = {
        "COI":"COX1","COXI":"COX1","MT-CO1":"COX1","CO1":"COX1",
        "COII":"COX2","COXII":"COX2","MT-CO2":"COX2","CO2":"COX2",
        "COIII":"COX3","COXIII":"COX3","MT-CO3":"COX3","CO3":"COX3",
        "CYTB":"CYTB","COB":"CYTB","MT-CYB":"CYTB","CYB":"CYTB",
        "ATPASE6":"ATP6","ATPASE8":"ATP8",
        "RRNS":"12S","12SRRNA":"12S","12S RRNA":"12S",
        "RRNL":"16S","16SRRNA":"16S","16S RRNA":"16S",
    }
    sub["GeneNorm"] = sub["GeneNorm"].replace(syn)
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
# Main plot
def plot_all_annotations(df: pd.DataFrame, out_png: Path, out_svg: Path, out_pdf: Path):
    """Draw linear mitochondrial gene maps for each record and export PNG/SVG/PDF."""
    import matplotlib as mpl
    mpl.rcParams['svg.fonttype'] = 'none'

    if df.empty:
        raise ValueError("Empty DataFrame.")

    len_by_rec = df.groupby("Record")["GenomeLen"].max().to_dict()
    max_len = max(len_by_rec.values()) if len_by_rec else 16000

    ordered_records = list(dict.fromkeys(df['Record'].tolist()))

    rec2label = {}
    if "Species" in df.columns and "Accession" in df.columns:
        g = df.groupby("Record").agg({"Species":"first","Accession":"first"})
        for rid, row in g.iterrows():
            rec2label[rid] = f"{row['Species']}\n({row['Accession']})"
    else:
        for rid in ordered_records:
            rec2label[rid] = rid

    W = 32
    H = max(6, 1.2 * len(ordered_records)) if ordered_records else 6
    fig = plt.figure(figsize=(W, H))
    ax = plt.gca()

    y = 0.0
    yticks, yticklabels = [], []
    y_step = 1.6
    bar_h = 0.6

    for rec_id in ordered_records:
        sub = df[df["Record"] == rec_id].sort_values("Start")
        L = int(sub["GenomeLen"].iloc[0]) if not sub.empty else max_len
        ax.hlines(y, 0, L, color="gray", linewidth=0.6, zorder=0)

        if not sub.empty:
            for _, r in sub.iterrows():
                name = r["Name"]
                width = r["End"] - r["Start"]
                cx = (r["Start"] + r["End"]) / 2.0
                if r["Type"] == "tRNA":
                    color = T_RNA_COLOR
                    ax.barh(y, width, left=r["Start"], color=color, edgecolor="black",
                            linewidth=1.0, height=bar_h, zorder=1)
                    aa_label = trna_label_one_letter(name)
                    ax.text(cx, y + 0.55, aa_label, ha='center', va='center',
                            fontsize=5, color='black', zorder=2)
                else:
                    color = GENE_COLORS.get(_normalize_gene_key(name), _fallback_color(_normalize_gene_key(name)))
                    ax.barh(y, width, left=r["Start"], color=color, edgecolor="black",
                            linewidth=1.0, height=bar_h, zorder=1)
                    ax.text(cx, y - 0.60, name, ha='center', va='center',
                            fontsize=5, color='#222222', zorder=2)

        yticks.append(y)
        yticklabels.append(rec2label.get(rec_id, rec_id))
        y += y_step

    ax.set_ylim(-2.0, y - 0.1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=8)

    ax.set_xlim(-max_len * 0.02, max_len * 1.03)
    ax.margins(x=0.03)

    ax.set_xlabel("Genomic position (bp)", fontsize=13)
    ax.tick_params(axis='x', labelsize=11, pad=25)
    ax.tick_params(axis='y', length=0)
    plt.title("Mitochondrial Map", fontsize=16, pad=10)
    ax.grid(axis='x', color='#dddddd', linewidth=0.6, alpha=0.8)

    legend_patches = [Patch(facecolor=GENE_COLORS[g], edgecolor="black", label=g) for g in TARGET_GENES]
    legend_patches.append(Patch(facecolor=T_RNA_COLOR, edgecolor="black", label="tRNA"))
    leg = ax.legend(handles=legend_patches, title="Legend (gene)", loc="center left",
                    bbox_to_anchor=(1.05, 0.5), frameon=True)
    leg.get_frame().set_linewidth(1.0)
    leg.get_frame().set_edgecolor('#444444')

    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.35)
    plt.tight_layout(rect=[0, 0.08, 0.83, 1])

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

# -----------------------------
# Presence heatmap
def plot_presence_heatmap(df: pd.DataFrame, out_png: Path, out_svg: Path, out_pdf: Path):
    """Draw presence/absence heatmap for canonical CDS/rRNA genes and export."""
    import matplotlib as mpl, numpy as np, matplotlib.colors as mcolors
    mpl.rcParams['svg.fonttype'] = 'none'

    genes = [g for g in TARGET_GENES]
    sub = df[df["Type"].isin(["CDS","rRNA"])].copy()
    if sub.empty:
        raise ValueError("No CDS/rRNA for presence heatmap.")

    syn = {
        "COI":"COX1","COXI":"COX1","MT-CO1":"COX1","CO1":"COX1",
        "COII":"COX2","COXII":"COX2","MT-CO2":"COX2","CO2":"COX2",
        "COIII":"COX3","COXIII":"COX3","MT-CO3":"COX3","CO3":"COX3",
        "CYTB":"CYTB","COB":"CYTB","MT-CYB":"CYTB","CYB":"CYTB",
        "ATPASE6":"ATP6","ATPASE8":"ATP8",
        "RRNS":"12S","12SRRNA":"12S","12S RRNA":"12S",
        "RRNL":"16S","16SRRNA":"16S","16S RRNA":"16S",
    }
    sub["GeneNorm"] = sub["Name"].apply(_normalize_gene_key)
    sub["GeneNorm"] = sub["GeneNorm"].replace(syn)
    mat = (sub.assign(val=1)
              .pivot_table(index="Record", columns="GeneNorm", values="val", aggfunc="max", fill_value=0))
    for g in genes:
        if g not in mat.columns:
            mat[g] = 0
    mat = mat[genes]
    ordered_records = list(dict.fromkeys(df['Record'].tolist()))
    mat = mat.reindex(index=ordered_records)

    if {"Species","Accession"}.issubset(df.columns):
        meta = df.groupby("Record")[["Species","Accession"]].first()
        ylabels = [f"{meta.loc[r,'Species']}\n({meta.loc[r,'Accession']})" if r in meta.index else r
                   for r in mat.index]
    else:
        ylabels = list(mat.index)

    data = mat.to_numpy()
    H = max(4, 0.4 * len(ordered_records))
    W = max(8, 0.6 * len(genes))
    fig, ax = plt.subplots(figsize=(W, H))
    cmap = mcolors.LinearSegmentedColormap.from_list("bin", ["#f0f0f0", "#6fbf73"])
    im = ax.imshow(data, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)

    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, fontsize=10, rotation=45, ha="right")
    ax.set_title("Gene Presence (CDS/rRNA)", fontsize=14, pad=10)
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

# -----------------------------
# CLI
# -----------------------------
# Normalized-position distance heatmap
def plot_normpos_distance_heatmap(df: pd.DataFrame, out_png: Path, out_svg: Path, out_pdf: Path, out_csv: Path):
    """
    Distance matrix between records based on the mean normalized positions of shared genes (CDS/rRNA).
    Distance = mean |Δ position| in [0,1]. Also exports the matrix to CSV.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    sub = df[df["Type"].isin(["CDS","rRNA"])].copy()
    if sub.empty:
        raise ValueError("No CDS/rRNA for normalized-position distance heatmap.")

    sub["mid_norm"] = ((sub["Start"] + sub["End"]) / 2.0) / sub["GenomeLen"]

    def _n(name: str) -> str:
        s = (name or "").upper().replace(" ", "")
        s = re.sub(r'^MT[-_]?','', s)
        if s in {"COX1","COI","CO1","COXI"}: return "COX1"
        if s in {"COX2","COII","CO2","COXII"}: return "COX2"
        if s in {"COX3","COIII","CO3","COXIII"}: return "COX3"
        if s in {"CYTB","COB","MTCYB","MT-CYB","CYB"}: return "CYTB"
        if s in {"ATP6","ATPASE6"}: return "ATP6"
        if s in {"ATP8","ATPASE8"}: return "ATP8"
        m = re.fullmatch(r"(?:NADH?|ND)(\d)(L?)", s)
        if m:
            return "ND" + m.group(1) + ("L" if m.group(2) == "L" else "")
        if "12S" in s or "RRNS" in s: return "12S"
        if "16S" in s or "RRNL" in s: return "16S"
        return s

    sub["GeneNorm"] = sub["Name"].map(_n)
    sub = sub.sort_values(["Record","GeneNorm","Start"]).drop_duplicates(["Record","GeneNorm"], keep="first")
    pivot = sub.pivot(index="Record", columns="GeneNorm", values="mid_norm")

    records = list(pivot.index)
    n = len(records)
    import numpy as np
    dist = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(n):
            pi = pivot.iloc[i]
            pj = pivot.iloc[j]
            common = pi.notna() & pj.notna()
            if not common.any():
                dist[i, j] = np.nan
            else:
                dist[i, j] = np.abs(pi[common].values - pj[common].values).mean()

    # export CSV
    dist_df = pd.DataFrame(dist, index=records, columns=records)
    dist_df.to_csv(out_csv)

    # plot
    vmin = np.nanmin(dist)
    vmax = np.nanmax(dist)
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(dist, interpolation="nearest", aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(records, rotation=45, ha="right")
    ax.set_yticklabels(records)
    ax.set_xlabel("Record")
    ax.set_ylabel("Record")
    ax.set_title("Normalized-position distance heatmap")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Mean |Δ position|")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description="Generate mitochondrial map (single file or directory with GB/GBK).")
    ap.add_argument("input", type=str, help="A .gb/.gbk (multi-entry) file OR a directory with GB/GBK files.")
    ap.add_argument("--out-prefix", type=str, default="mtmap_out", help="Prefix of output files.")
    ap.add_argument("--anchor-gene", type=str, default=None,
                    help="Anchoring gene to linearize the circular genome (e.g., COX1, ND1, CYTB, 'F' for tRNA-Phe).")

    args = ap.parse_args()
    in_path = Path(args.input)

    df = extract_annotations_any(in_path)
    if df.empty:
        raise SystemExit("No annotations found in the provided GB/GBK files.")

    if args.anchor_gene:
        df = realign_to_anchor(df, args.anchor_gene)

    out_prefix = Path(args.out_prefix)
    out_png = out_prefix.with_name(out_prefix.name + "_map.png")
    out_svg = out_prefix.with_name(out_prefix.name + "_map.svg")
    out_pdf = out_prefix.with_name(out_prefix.name + "_map.pdf")
    plot_all_annotations(df, out_png, out_svg, out_pdf)

    export_tables(df, out_prefix)
    plot_presence_heatmap(df,
                          out_prefix.with_name(out_prefix.name + "_presence.png"),
                          out_prefix.with_name(out_prefix.name + "_presence.svg"),
                          out_prefix.with_name(out_prefix.name + "_presence.pdf"))

    # new: normalized-position distance heatmap
    plot_normpos_distance_heatmap(
        df,
        out_prefix.with_name(out_prefix.name + "_normposdist.png"),
        out_prefix.with_name(out_prefix.name + "_normposdist.svg"),
        out_prefix.with_name(out_prefix.name + "_normposdist.pdf"),
        out_prefix.with_name(out_prefix.name + "_normposdist.csv")
    )

if __name__ == "__main__":
    main()
