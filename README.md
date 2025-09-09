# MTMap — Mitochondrial Genome Mapper (`mtmap.py`)

## Overview  

`mtmap.py` is a Python-based tool designed for the **systematic extraction, normalization, and visualization of mitochondrial genome annotations** from GenBank files. The pipeline integrates **data parsing, canonical gene normalization, and multi-format visualization** to support comparative mitochondrial genomics, evolutionary studies, and molecular systematics.  

The tool automatically generates:  

- **Linear mitochondrial genome maps**, displaying coding sequences (CDS), rRNAs, and tRNAs.  
- **Annotation tables** (`CSV`) with standardized gene nomenclature.  
- **Binary gene presence/absence matrices** across multiple records.  
- **Heatmaps** representing gene presence across species.  
- **Multi-format outputs** (`PNG`, `SVG`, `PDF`) suitable for publications or presentations.  

---

## Scientific Rationale  

Mitochondrial genomes (mtDNA) are a cornerstone in phylogenetics, population genetics, and molecular systematics. However, annotations from GenBank are often inconsistent due to synonym usage (`COI`, `COX1`, `MT-CO1`) and variation in feature naming. This tool addresses these issues by:  

1. **Normalization of gene nomenclature**: Standardizing synonyms to canonical forms (e.g., `COI` → `COX1`, `COB` → `CYTB`).  
2. **Anchored genome linearization**: Re-aligning the circular mitochondrial genome based on a selected reference gene (e.g., `COX1`, `ND1`, `CYTB`, or a specific tRNA).  
3. **Comparative visualization**: Producing presence/absence matrices and annotated maps for comparative mitochondrial genomics.  
4. **Multi-format outputs**: Allowing downstream use in both automated pipelines and human-readable figures.  

---

## Key Features  

- **Flexible Input**: Accepts single GenBank files or directories containing multiple entries.  
- **Robust Parsing**: Extracts `CDS`, `rRNA`, and `tRNA` features, including metadata such as species and accession.  
- **Gene Anchoring**: Linearization of circular genomes based on user-defined anchor gene.  
- **Normalization Layer**: Maps synonyms and heuristic product descriptions to canonical keys.  
- **Data Export**:  
  - Complete annotations (`*_annotations.csv`)  
  - Gene presence/absence matrix (`*_presence.csv`)  
- **Visualization**:  
  - Linear mitochondrial genome maps  
  - Heatmap of presence/absence for canonical mitochondrial genes  

---

## Target Genes  

The presence/absence matrix evaluates canonical mitochondrial genes:  

```
COX1, COX2, COX3,
ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
CYTB, ATP6, ATP8, 12S, 16S
```  

In addition, all **tRNAs** are included and represented with one-letter amino acid codes.  

---

## Installation  

### Requirements  

- Python **≥ 3.8**  
- Required libraries:  
  ```bash
  pip install biopython pandas matplotlib
  ```

### Repository Structure  

```
mtmap.py             # Main script
examples/            # Example GenBank files (user-provided)
html/                # MTMap — Mitochondrial Genome Mapper (Pyodide)
README.md            # Documentation
```

---

## 🚀 Usage

### Basic command
```bash
python mtmap.py my_genome.gbk --out-prefix results/mtmap
```

### Multiple genomes
```bash
python mtmap.py genomes_directory/ --out-prefix results/comparison
```

### With anchor gene
```bash
python mtmap.py genomes_directory/ --out-prefix results/anchored --anchor-gene COX1
```

---

## 📊 Output Files

Given `--out-prefix results/mtmap`, the tool will generate:

- **Tables**
  - `mtmap_annotations.csv` → Full annotation table
  - `mtmap_presence.csv` → Gene presence/absence matrix
  - `mtmap_normposdist.csv` → Distance matrix (normalized positions)

- **Plots**
  - `mtmap_map.(png|svg|pdf)` → Linear mitochondrial genome map
  - `mtmap_presence.(png|svg|pdf)` → Presence/absence heatmap
  - `mtmap_normposdist.(png|svg|pdf)` → Distance heatmap

---

## Example  

### Command  

```bash
python mtmap.py examples/ --out-prefix results/mtDNA --anchor-gene COX1
```

### Generated Outputs  

- `results/mtDNA_map.png` → Annotated mitochondrial genome maps.  
- `results/mtDNA_presence.csv` → Presence/absence matrix of canonical genes.  
- `results/mtDNA_presence.png` → Heatmap summarizing gene distribution.  

---

## Technical Notes  

- **Normalization Strategy**  
  - Synonym resolution: e.g., `COI`, `COXI`, `MT-CO1` → `COX1`.  
  - Heuristic recognition based on product names (`NADH dehydrogenase subunit 4L` → `ND4L`).  
- **Visualization**  
  - Coding genes: pastel color scheme (deterministic mapping).  
  - rRNAs: distinct category (`12S`, `16S`).  
  - tRNAs: pale yellow with one-letter amino acid labels.  
- **Robustness**  
  - Handles multi-entry GenBank files.  
  - Fallback colors generated deterministically for unknown/non-target genes.  

---

## Applications  

- Comparative genomics of mitochondrial gene content.  
- Phylogenetic and evolutionary studies.  
- Validation of mitochondrial genome assemblies.  
- Educational visualization of mtDNA architecture.  

---

## Citation  

If you use this tool in scientific work, please cite as:  

> **LaBiOmicS (2025)**. *Mitochondrial Map (mtmap.py): a Python tool for comparative mitochondrial genome visualization and annotation normalization*. GitHub Repository.  

---

## License  

Distributed under the **MIT License**.  
You are free to use, modify, and distribute this tool with appropriate credit.  


---

## Reproducible Simulation (Demo)

This section reproduces a complete run of `mtmap.py` using a panel of publicly available vertebrate mitochondrial GenBank files (uploaded to this repository's `/examples` directory). We anchor the circular mtDNA on **COX1** for consistent linearization across species.

### Command

```bash
python mtmap.py /path/to/examples --out-prefix mtmap_demo --anchor-gene COX1
```

### Generated Outputs (this demo)

- **Genome Map**: [mtmap_demo_map.png](outputs/mtmap_demo_map.png) · [SVG](output/mtmap_demo_map.svg) · [PDF](output/mtmap_demo_map.pdf)
- **Presence/Absence Heatmap**: [mtmap_demo_presence.png](output/mtmap_demo_presence.png) · [SVG](output/mtmap_demo_presence.svg) · [PDF](output/mtmap_demo_presence.pdf)
- **Annotation Table (CSV)**: [mtmap_demo_annotations.csv](output/mtmap_demo_annotations.csv)
- **Gene Presence Matrix (CSV)**: [mtmap_demo_presence.csv](output/mtmap_demo_presence.csv)

> Notes
> - The demo uses the following GenBank files:
>   NC_000861_Salvelinus_alpinus.gbk, NC_000890_Mustelus_manazo.gbk, NC_000893_Amblyraja_radiata.gbk, NC_001131_Lampetra_fluviatilis.gbk, NC_001606_Cyprinus_carpio.gbk, NC_001626_Petromyzon_marinus.gbk, NC_001708_Protopterus_dolloi.gbk, NC_001717_Oncorhynchus_mykiss.gbk, NC_001727_Formosania_lacustris.gbk, NC_001778_Polypterus_ornatipinnis.gbk
> - Anchoring to `COX1` aligns the 0 bp position to the start of `COX1` (if present); records lacking the anchor remain unshifted.
> - Figures are exported in publication-ready formats (PNG 300 dpi, SVG with text preserved, and PDF).

