# MTMap — Mitochondrial Genome Mapper (Pyodide)

Client-side mitochondrial genome mapper built with Pyodide + Biopython.
- Upload up to **20** GenBank files (.gb/.gbk/.gbff)
- Generates: linear genome map (PNG), CDS/rRNA presence heatmap (PNG),
  annotations CSV, presence CSV, and a ZIP with all outputs.
- Runs **fully in the browser** (no server).

## Use
Open: https://<seu-usuario>.github.io/mtmap-pyodide/

## Develop
Just edit `index.html` and push. No build required.

## Limits
GitHub repo limit: ~100 MB/file, ~1 GB/repo. Browser memory/CPU bound.
