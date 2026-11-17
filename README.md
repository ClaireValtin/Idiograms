# Idiograms

Interactive human chromosome idiograms, a microscope simulator, structural abnormalities and FISH probes.

**Open it:** https://clairevaltin.github.io/Idiograms/

## What it does

A cytogenetics teaching and reference site with five pages:

- **Idiograms** – every human chromosome in G- or R-banding at 400, 550 and ~700 band resolution, with genes and GRCh38 coordinates
- **Microscope** – a simulated banded metaphase (46,XX or 46,XY) with adjustable staining, focus and chromosome bending, to practise recognising chromosomes
- **Abnormalities** – type an ISCN karyotype (e.g. `46,XY,t(9;22)(q34.1;q11.2)`) to see it drawn, or build rearrangements by cutting chromosomes: translocations, deletions, duplications, inversions, isochromosomes, rings
- **FISH probes** – interphase and metaphase signal patterns for common probes
- **Probe catalogue** – a searchable list of FISH probes with their target genes

Band names, lengths and coordinates come from NCBI's ideogram tables, which match the UCSC hg38 cytoband table. Each page has an FAQ covering its sources and limits.

## How to use

Open the link above. Each page is a self-contained HTML file, so it also works offline.

## Repository contents

- `src/` – page templates and shared code
- `data/` – band, gene, probe and catalogue data (JSON)
- `tools/` – scripts that generate and check the data
- `build.py` – rebuilds the HTML pages from `src/` and `data/`. Edit those, then run `python build.py`; never edit the pages by hand.

## Disclaimer

**This tool is NOT validated for clinical or diagnostic use.** It is provided for research and educational purposes only, without any warranty.

## License & citation

© 2026 Rehioui Soufiyan. All rights reserved.

You're welcome to explore this code. If you'd like to use, adapt, or share it, please get in touch first: rehiouisoufiyan@gmail.com

If this work helps you, please cite it or link back to this repository.
