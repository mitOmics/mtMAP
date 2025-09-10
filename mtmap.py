#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mtmap_plus_refactored.py — Mitochondrial Map (refactor)
- Centralizes gene-name normalization (normalize_gene)
- Presence heatmap (CDS/rRNA) with synonyms applied
- Normalized-position distance heatmap (CDS/rRNA)
"""
from pathlib import Path
import argparse
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import SeqIO

# -----------------------------
# Fallback color for non-target genes
def _fallback_color(key: str) -> str:
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
TARGET_GENES = ["COX1","COX2","COX3","CYTB","ND1","ND2","ND3","ND4","ND4L","ND5","ND6","ATP6","ATP8","12S","16S"]
T_RNA_COLOR = "#ffd27a"

def _gene_color_map():
    import matplotlib as mpl
    colors = {}
    tab = list(mpl.colormaps.get_cmap('tab20').colors)
    for i,g in enumerate(TARGET_GENES):
        rgba = tab[i % 20]
        # lighten
        r,gc,b = rgba[:3]
        r = r + (1.0 - r) * 0.35
        gc = gc + (1.0 - gc) * 0.35
        b = b + (1.0 - b) * 0.35
        colors[g] = (r, gc, b)
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
# Centralized normalization
def normalize_gene(name: str) -> str:
    import re as _re
    if not name:
        return ""
    s = name.upper().strip().replace(" ", "")
    s = _re.sub(r'^MT[-_]?','', s)

    if "12S" in s or "RRNS" in s or "12SRRNA" in s or "12SRNA" in s: return "12S"
    if "16S" in s or "RRNL" in s or "16SRRNA" in s or "16SRNA" in s: return "16S"

    if s in {"COX1","COI","CO1","COXI"}: return "COX1"
    if s in {"COX2","COII","CO2","COXII"}: return "COX2"
    if s in {"COX3","COIII","CO3","COXIII"}: return "COX3"

    if s in {"CYTB","CYB","COB","CYTOCHROMEB","MTCYB","MT-CYB"}: return "CYTB"

    if s in {"ATP6","ATPASE6"}: return "ATP6"
    if s in {"ATP8","ATPASE8"}: return "ATP8"

    m = _re.fullmatch(r"(?:NADH?|ND)(\d)(L?)", s)
    if m: return "ND" + m.group(1) + ("L" if m.group(2) == "L" else "")

    sp = name.upper()
    if "CYTOCHROME B" in sp: return "CYTB"
    if "CYTOCHROME C OXIDASE" in sp:
        if _re.search(r"\bSUBUNIT\s*(I|1)\b", sp): return "COX1"
        if _re.search(r"\bSUBUNIT\s*(II|2)\b", sp): return "COX2"
        if _re.search(r"\bSUBUNIT\s*(III|3)\b", sp): return "COX3"
    if "ATP SYNTHASE" in sp:
        if _re.search(r"\bSUBUNIT\s*6\b", sp): return "ATP6"
        if _re.search(r"\bSUBUNIT\s*8\b", sp): return "ATP8"
    if "NADH DEHYDROGENASE SUBUNIT" in sp:
        mm = _re.search(r"SUBUNIT\s+(\d)", " "+sp)
        if mm:
            nd = "ND" + mm.group(1)
            if " ND4L" in sp or " SUBUNIT 4L" in sp: nd = "ND4L"
            return nd
    if "SMALL SUBUNIT RIBOSOMAL RNA" in sp or "12S RIBOSOMAL RNA" in sp: return "12S"
    if "LARGE SUBUNIT RIBOSOMAL RNA" in sp or "16S RIBOSOMAL RNA" in sp: return "16S"
    return s

# -----------------------------
def list_gb_paths(input_path: Path):
    GB_EXTS = {".gb",".gbk"}
    if input_path.is_dir():
        files = sorted([p for p in input_path.iterdir() if p.suffix.lower() in GB_EXTS])
        if not files:
            raise FileNotFoundError(f"Nenhum arquivo GB/GBK encontrado em: {input_path}")
        return files
    if input_path.is_file():
        return [input_path]
    raise FileNotFoundError(f"Caminho não encontrado: {input_path}")

# -----------------------------
# GenBank parsing
def extract_annotations_any(input_path: Path) -> pd.DataFrame:
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
                if feat.type not in {"CDS","rRNA","tRNA"}:
                    continue
                name = None
                if "gene" in feat.qualifiers:
                    name = feat.qualifiers.get("gene",[None])[0]
                elif "product" in feat.qualifiers:
                    name = feat.qualifiers.get("product",[None])[0]
                start = int(feat.location.start)
                end = int(feat.location.end)
                rows.append({
                    "Record": rec_id,
                    "Type": "CDS" if feat.type=="CDS" else feat.type,
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
def realign_to_anchor(df: pd.DataFrame, anchor: str) -> pd.DataFrame:
    if not anchor or df.empty:
        return df
    out = []
    for rec_id, sub in df.groupby("Record"):
        s = sub.copy()
        L = int(s["GenomeLen"].iloc[0])
        # search for the anchor gene midpoint among normalized CDS/rRNA/tRNA
        ss = s.copy()
        ss["GeneNorm"] = ss["Name"].apply(normalize_gene)
        hit = ss[ss["GeneNorm"] == anchor]
        offset = 0
        if not hit.empty:
            offset = int((hit["Start"].iloc[0] + hit["End"].iloc[0])//2)
        ns = (s['Start'] - offset) % L
        ne = (s['End'] - offset) % L
        wrap = ne <= ns
        s.loc[:, 'Start'] = ns.astype(int)
        s.loc[:, 'End'] = (ne + wrap * L).astype(int)
        out.append(s)
    return pd.concat(out, ignore_index=True)

# -----------------------------
# Main plot (simple linear map)
def plot_all_annotations(df: pd.DataFrame, out_png: Path, out_svg: Path, out_pdf: Path):
    import matplotlib as mpl
    mpl.rcParams['svg.fonttype'] = 'none'

    fig, ax = plt.subplots(figsize=(11, max(4, 0.5 * df["Record"].nunique())))
    bar_h = 8
    ordered_records = list(dict.fromkeys(df['Record'].tolist()))
    max_len = df.groupby("Record")["GenomeLen"].max().max()

    y = 0
    for rec_id in ordered_records:
        sub = df[df["Record"] == rec_id].sort_values("Start")
        L = int(sub["GenomeLen"].iloc[0]) if not sub.empty else max_len
        ax.hlines(y, 0, L, color="gray", linewidth=0.6, zorder=0)
        for _, r in sub.iterrows():
            name = r["Name"]
            width = r["End"] - r["Start"]
            if r["Type"] == "tRNA":
                color = T_RNA_COLOR
                ax.barh(y, width, left=r["Start"], color=color, edgecolor="black",
                        linewidth=1.0, height=bar_h, zorder=1)
                aa_label = trna_label_one_letter(name or "")
                ax.text(r["Start"]+width/2, y, aa_label, ha="center", va="center", fontsize=7)
            else:
                gn = normalize_gene(name or "")
                color = GENE_COLORS.get(gn, _fallback_color(gn))
                ax.barh(y, width, left=r["Start"], color=color, edgecolor="black",
                        linewidth=1.0, height=bar_h, zorder=1)
                if gn in TARGET_GENES:
                    ax.text(r["Start"]+width/2, y, gn, ha="center", va="center", fontsize=7, rotation=0)
        y += bar_h + 10

    ax.set_xlim(0, max_len)
    ax.set_ylim(-bar_h, y + bar_h)
    ax.set_xlabel("Genomic position (bp)")
    ax.set_yticks([])
    legend_patches = [Patch(facecolor=GENE_COLORS[g], edgecolor="black", label=g) for g in TARGET_GENES]
    legend_patches.append(Patch(facecolor=T_RNA_COLOR, edgecolor="black", label="tRNA"))
    leg = ax.legend(handles=legend_patches, title="Legend (gene)", loc="center left",
                    bbox_to_anchor=(1.02, 0.5), frameon=True)
    leg.get_frame().set_linewidth(1.0)
    leg.get_frame().set_edgecolor('#444444')
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

# -----------------------------
# Presence heatmap (CDS/rRNA)
def plot_presence_heatmap(df: pd.DataFrame, out_png: Path, out_svg: Path, out_pdf: Path, out_csv: Path):
    import matplotlib as mpl, numpy as np, matplotlib.colors as mcolors
    mpl.rcParams['svg.fonttype'] = 'none'
    genes = [g for g in TARGET_GENES]
    sub = df[df["Type"].isin(["CDS","rRNA"])].copy()
    if sub.empty:
        raise ValueError("No CDS/rRNA for the presence heatmap.")
    syn = {
        "COI":"COX1","COXI":"COX1","MT-CO1":"COX1","CO1":"COX1",
        "COII":"COX2","COXII":"COX2","MT-CO2":"COX2","CO2":"COX2",
        "COIII":"COX3","COXIII":"COX3","MT-CO3":"COX3","CO3":"COX3",
        "CYTB":"CYTB","COB":"CYTB","MT-CYB":"CYTB","CYB":"CYTB",
        "ATPASE6":"ATP6","ATPASE8":"ATP8",
        "RRNS":"12S","12SRRNA":"12S","12S RRNA":"12S",
        "RRNL":"16S","16SRRNA":"16S","16S RRNA":"16S",
    }
    sub["GeneNorm"] = sub["Name"].apply(normalize_gene)
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
    mat.to_csv(out_csv)

    ordered_records = list(dict.fromkeys(df['Record'].tolist()))
    mat = mat.set_index(mat.columns[0] if "Species" in mat.columns else "Record")
    mat = mat.reindex(index=ordered_records)

    data = mat[genes].to_numpy()
    H = max(4, 0.4 * len(ordered_records))
    W = max(8, 0.6 * len(genes))
    fig, ax = plt.subplots(figsize=(W, H))
    cmap = mcolors.LinearSegmentedColormap.from_list("bin", ["#f0f0f0", "#6fbf73"])
    im = ax.imshow(data, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    ax.set_yticks(range(len(ordered_records)))
    ax.set_yticklabels(ordered_records, fontsize=8)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_title("Presence Heatmap (CDS/rRNA)")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_ticks([0,1]); cbar.set_ticklabels(["Absent","Present"])
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

# -----------------------------
# Normalized-position distance heatmap
def plot_normpos_distance_heatmap(df: pd.DataFrame, out_png: Path, out_svg: Path, out_pdf: Path, out_csv: Path):
    import numpy as np
    import matplotlib.pyplot as plt
    sub = df[df["Type"].isin(["CDS","rRNA"])].copy()
    if sub.empty:
        raise ValueError("No CDS/rRNA for the normalized-position distance heatmap.")
    sub["mid_norm"] = ((sub["Start"] + sub["End"]) / 2.0) / sub["GenomeLen"]
    sub["GeneNorm"] = sub["Name"].apply(normalize_gene)
    sub = sub.sort_values(["Record","GeneNorm","Start"]).drop_duplicates(["Record","GeneNorm"], keep="first")
    pivot = sub.pivot(index="Record", columns="GeneNorm", values="mid_norm")
    records = list(pivot.index)
    n = len(records)
    dist = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            pi = pivot.iloc[i]
            pj = pivot.iloc[j]
            common = pi.notna() & pj.notna()
            dist[i, j] = np.nan if not common.any() else abs(pi[common].values - pj[common].values).mean()
    dist_df = pd.DataFrame(dist, index=records, columns=records)
    dist_df.to_csv(out_csv)

    vmin = float(np.nanmin(dist)); vmax = float(np.nanmax(dist))
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(dist, interpolation="nearest", aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(records, rotation=45, ha="right"); ax.set_yticklabels(records)
    ax.set_xlabel("Record"); ax.set_ylabel("Record")
    ax.set_title("Normalized-position distance heatmap")
    cbar = plt.colorbar(im, ax=ax); cbar.set_label("Mean |Δ position|")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

# -----------------------------
def export_tables(df: pd.DataFrame, out_prefix: Path):
    out_ann = out_prefix.with_name(out_prefix.name + "_annotations.csv")
    df.sort_values(["Record","Start","End"]).to_csv(out_ann, index=False)

# -----------------------------
def main():
    ap = argparse.ArgumentParser(description="Generate maps and heatmaps from GB/GBK files.")
    ap.add_argument("input", type=str, help="A .gb/.gbk file (multi-entry) OR a directory with GB/GBK files.")
    ap.add_argument("--out-prefix", type=str, default="mtmap_out", help="Output file prefix.")
    ap.add_argument("--anchor-gene", type=str, default=None,
                    help="Anchor gene to linearize the circular genome (e.g., COX1, ND1, CYTB, 'F' for tRNA-Phe).")
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

    # presence
    plot_presence_heatmap(
        df,
        out_prefix.with_name(out_prefix.name + "_presence.png"),
        out_prefix.with_name(out_prefix.name + "_presence.svg"),
        out_prefix.with_name(out_prefix.name + "_presence.pdf"),
        out_prefix.with_name(out_prefix.name + "_presence.csv"),
    )

    # normalized-position distance
    plot_normpos_distance_heatmap(
        df,
        out_prefix.with_name(out_prefix.name + "_normposdist.png"),
        out_prefix.with_name(out_prefix.name + "_normposdist.svg"),
        out_prefix.with_name(out_prefix.name + "_normposdist.pdf"),
        out_prefix.with_name(out_prefix.name + "_normposdist.csv"),
    )

    print("\nGenerated outputs:")
    outs = [
        "_map.png","_map.svg","_map.pdf",
        "_annotations.csv",
        "_presence.csv","_presence.png","_presence.svg","_presence.pdf",
        "_normposdist.csv","_normposdist.png","_normposdist.svg","_normposdist.pdf"
    ]
    for suf in outs:
        print(" -", (out_prefix.with_name(out_prefix.name + suf)))

if __name__ == "__main__":
    main()
