# Idiograms

Five self-contained HTML pages for human cytogenetics, with no dependencies:

- R- and G-banded idiograms of all 24 chromosomes at the 400, 550 and ~700-band levels.
- A microscope simulator that draws the same chromosomes as they look on a slide.
- An abnormality builder that writes ISCN nomenclature.
- FISH probe signal patterns (interphase and metaphase).
- A searchable catalogue of FISH probes from six vendors.

Notes, sources and limits are behind the **FAQ** button on every page (`src/_faq.html`).

## Pages

| Page | What it does |
|---|---|
| `index.html` (Idiograms) | R/G switch, a one-row karyogram, and the detail view of one chromosome at 400/550/~700 bands. Hover a band to outline it at every level. Search by band (`9q34`) or landmark gene (`ABL1`). |
| `simulator.html` (Microscope) | Simulated R- (default) or G-banded chromosomes: a homologue pair or a full karyogram. Sliders: bands, focus, staining, chromatids, bending, field (debris, cytoplasm, nuclei). Giemsa or grey, heterochromatin tone. |
| `abnormalities.html` (Abnormalities) | Type a karyotype (a trailing `[n]` cell count is kept) or pick an example. Or build one by hand on the bench: click bands, drag the breakpoint markers. Shows the short and detailed ISCN, the genes at the breakpoints, and the result as idiograms or under the microscope. `#k=<karyotype>` opens a karyotype. |
| `fish.html` (FISH probes) | 148 probes: every Vysis design, plus CytoCell and MetaSystems probes for the haematology targets Vysis lacks. Each probe has a map on GRCh38 with genes, interphase nuclei for its normal and abnormal patterns (signal code and `nuc ish` string), and a metaphase FISH view. `#probe-id` links to a probe. |
| `catalogue.html` (Probe catalogue) | 2,964 human FISH probes from Vysis, CytoCell, ZytoVision, MetaSystems, Leica (Kreatech) and Empire Genomics. Search, vendor, chromosome and status filters, with links to each vendor page and to the FISH page. |

Pages remember the last chromosome across tabs; `#chr9` etc. opens one directly.

The site is published with GitHub Pages from `main` at `https://clairevaltin.github.io/Idiograms/`. Links between pages carry a `?v=…` tag that changes with every build, so browsers never show a stale page.

## How this repo is organised

| Path | What it is |
|---|---|
| `*.html` (root) | The website. **Generated, do not edit by hand.** |
| `src/*.template.html` | One template per page. |
| `src/_shared.css`, `src/_shared.js`, `src/_faq.html` | App bar, styles, FAQ and helpers shared by every page. |
| `src/_scope.js` | Microscope image engine (bright-field R/G and fluorescence). Used by the Microscope, Abnormalities and FISH pages. |
| `src/_fish.js` | FISH engine: applies a rearrangement to a diploid cell and derives the visible signals, the signal code and the ISCN `nuc ish` string. Also runs under node. |
| `data/bands.json`, `data/landmark-genes.json`, `data/breakpoint-genes.json` | Band and gene data. |
| `data/probes.json`, `data/probe-genes.json` | FISH probe designs placed on GRCh38, and the gene coordinates they use. Generated. |
| `data/catalogue.json` | The vendor catalogue. Generated. |
| `build.py` | Combines templates and data into the pages. Runs the probe check first. |
| `tools/make_data.py` | Rebuilds the band and gene data from NCBI/Ensembl. |
| `tools/probe_specs.py`, `tools/make_probes.py` | Hand-curated probe designs (from vendor maps and inserts) and the script that places them on GRCh38 using Ensembl gene coordinates (`--refresh` re-fetches them). |
| `tools/check_probes.mjs` | For every pattern, checks that the signal code derived from the probe positions equals the vendor's (or textbook) code. `-v` prints each pattern with its ISCN. |
| `tools/make_catalogue.py` | Builds the catalogue from saved vendor listings (`python3 tools/make_catalogue.py RAW_DIR`). |

After changing anything in `src/` or `data/`:

```
python3 tools/make_probes.py      # only if tools/probe_specs.py changed
python3 build.py                  # checks the probes, then writes the pages
```

## Data

| File | Contents | Source |
|---|---|---|
| `data/bands.json` | Every band at 400, 550 and GRCh38 (~700, key `850`): `[band, iscn_start, iscn_stop, bp_start, bp_stop, ncbi_stain]` | NCBI human ideogram data, GRCh38, via [ideogram.js](https://github.com/eweitz/ideogram) 1.53.0. The ~700 level matches UCSC hg38 cytoBand exactly. The 400 level matches NCBI's file exactly. |
| `data/landmark-genes.json`, `data/breakpoint-genes.json` | Symbol, GRCh38 start/end, band, note | Ensembl release 110 GRCh38 |
| `data/probes.json` | Probe components (colour, size, GRCh38 segments) and signal patterns (events, expected code, source) | Abbott/Vysis product pages, probe maps and package inserts; OGT CytoCell and MetaSystems product pages; Ensembl REST canonical transcripts |
| `data/catalogue.json` | Vendor, name, catalogue numbers, targets, colours, status, URL | Vendor websites, October 2026 |

- **R pattern:** drawn as the reciprocal of NCBI's G-stain grades, as ISCN treats R-banded idiograms.
- **550-level gaps:** NCBI's file leaves gaps at 7q11.2, 8q11.2, 11p11.1, 12q24.3, 18p11.3, Xp11.2 and Yq11.2. The sub-bands there were rescaled to fill their parent band.

## Main limitations

The FAQ has the full list. The main ones:

- **The microscope is a model, not a photograph.** Band intensities start from NCBI stain grades, and condensation is simulated by smoothing. Heterochromatin tones and the acrocentric short-arm lengths are judgements.
- **The finest level is GRCh38 (~700), not ISCN 850.** ISCN 300- and 850-band idiograms are not in the source data.
- **FISH probe positions are approximate**, placed from vendor maps. MetaSystems pages give no sizes, so those positions are illustrative. Signal patterns are derived and checked against expected codes; the `nuc ish` strings have not yet been checked against the ISCN 2024 text.
- **Not for diagnostic use.** Vendor names are trademarks of their owners; this site is not affiliated with any of them.

## Landmark genes

| Chr | Genes (band) |
|---|---|
| 1 | MUTYH (1p34.1), NRAS (1p13.2), CKS1B (1q21.3) |
| 2 | MYCN (2p24.3), ALK (2p23.2–p23.1), MSH2 (2p21–p16.3) |
| 3 | VHL (3p25.3), MLH1 (3p22.2), MECOM (3q26.2) |
| 4 | FGFR3 (4p16.3), KIT (4q12), TET2 (4q24) |
| 5 | APC (5q22.2), NPM1 (5q35.1), NSD1 (5q35.3) |
| 6 | DEK (6p22.3), HLA-A (6p22.1), MYB (6q23.3) |
| 7 | EGFR (7p11.2), ELN (7q11.23), CFTR (7q31.2–q31.31) |
| 8 | FGFR1 (8p11.23), RUNX1T1 (8q21.3), MYC (8q24.21) |
| 9 | JAK2 (9p24.1), CDKN2A (9p21.3), ABL1 (9q34.12) |
| 10 | RET (10q11.21), PTEN (10q23.31) |
| 11 | WT1 (11p13), CCND1 (11q13.3), ATM (11q22.3), KMT2A (11q23.3) |
| 12 | ETV6 (12p13.2), KRAS (12p12.1), PAH (12q23.2) |
| 13 | FLT3 (13q12.2), BRCA2 (13q13.1), RB1 (13q14.2) |
| 14 | SERPINA1 (14q32.13), IGH (14q32.33) |
| 15 | SNRPN (15q11.2), FBN1 (15q21.1), PML (15q24.1) |
| 16 | TSC2 (16p13.3), MYH11 (16p13.11), CBFB (16q22.1) |
| 17 | TP53 (17p13.1), NF1 (17q11.2), RARA (17q21.2), BRCA1 (17q21.31) |
| 18 | SMAD4 (18q21.2), BCL2 (18q21.33) |
| 19 | STK11 (19p13.3), LDLR (19p13.2), CEBPA (19q13.11) |
| 20 | JAG1 (20p12.2), ASXL1 (20q11.21), GNAS (20q13.32) |
| 21 | APP (21q21.3), RUNX1 (21q22.12), ERG (21q22.2) |
| 22 | TBX1 (22q11.21), BCR (22q11.23), NF2 (22q12.2) |
| X | SHOX (Xp22.33), DMD (Xp21.2–p21.1), AR (Xq12), XIST (Xq13.2), FMR1 (Xq27.3), MECP2 (Xq28) |
| Y | SRY (Yp11.2), USP9Y (Yq11.221), DAZ1 (Yq11.223) |
