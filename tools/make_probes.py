#!/usr/bin/env python3
"""Build data/probes.json from tools/probe_specs.py.

    python tools/make_probes.py            # uses the cached gene coordinates in data/probe-genes.json
    python tools/make_probes.py --refresh  # re-fetches gene coordinates from the Ensembl REST API

Each probe region in probe_specs.py is written relative to genes ("200 kb telomeric of the
3' end of BCR"). This script looks up the GRCh38 coordinates of those genes (the Ensembl
canonical transcript), resolves every region and breakpoint to base pairs, and writes the
result for the FISH page. Centromeres, heterochromatin blocks and chromosome ends come from
data/bands.json, so probes and idiograms share one coordinate system.
"""
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import probe_specs as S  # noqa: E402

GENES_FILE = ROOT / "data" / "probe-genes.json"
OUT = ROOT / "data" / "probes.json"

# Genes whose canonical transcript starts downstream of the exons where breakpoints occur:
# use the whole gene (ABL1 from exon 1b, CCND3 and ERG with their upstream first exons).
WHOLE_GENE = {
    "ABL1": {"chr": "9", "start": 130713043, "end": 130887675, "strand": 1},
    "CCND3": {"chr": "6", "start": 41934577, "end": 42050357, "strand": -1},
    "ERG": {"chr": "21", "start": 38380027, "end": 38661963, "strand": -1},
}

# Immunoglobulin and T-cell receptor loci are not single genes in Ensembl; GRCh38 extents of the loci.
PSEUDO = {
    "TRA": {"chr": "14", "start": 21621904, "end": 22552132, "strand": 1},
    "TRB": {"chr": "7", "start": 142299011, "end": 142813287, "strand": 1},
    "IGL": {"chr": "22", "start": 22026076, "end": 22922913, "strand": 1},
}

BANDS = json.loads((ROOT / "data" / "bands.json").read_text())["chromosomes"]
CHR_LEN = {c: v["850"][-1][4] for c, v in BANDS.items()}
CEN = {}
for c, v in BANDS.items():
    ac = [b for b in v["850"] if b[5] == "acen"]
    CEN[c] = (ac[0][3], ac[-1][4], ac[0][4])  # start, end, p/q boundary
PAR1_X = (10001, 2781479)
PAR2_X = (155701383, 156030895)
PAR2_Y0 = 56887903


def walk(obj, found):
    if isinstance(obj, tuple):
        if obj and obj[0] in ("gene", "5", "3", "s", "e", "f"):
            found.add(obj[1])
        for x in obj:
            walk(x, found)
    elif isinstance(obj, (list, dict)):
        for x in (obj.values() if isinstance(obj, dict) else obj):
            walk(x, found)


def fetch_genes(symbols):
    out = {}
    symbols = [s for s in sorted(symbols) if s not in PSEUDO]
    for i in range(0, len(symbols), 50):
        chunk = symbols[i:i + 50]
        req = urllib.request.Request(
            "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
            data=json.dumps({"symbols": chunk, "expand": 1}).encode(),
            headers={"Content-Type": "application/json", "Accept": "application/json"})
        data = json.loads(urllib.request.urlopen(req, timeout=120).read())
        for sym in chunk:
            g = data.get(sym)
            if not g:
                raise SystemExit(f"Ensembl has no gene {sym}")
            tx = [t for t in g.get("Transcript", []) if t.get("is_canonical")]
            t = tx[0] if tx else g
            out[sym] = {"chr": g["seq_region_name"], "start": t["start"], "end": t["end"], "strand": g["strand"], "id": g["id"]}
    return out


def load_genes(refresh):
    need = set()
    walk(S.PROBES, need)
    walk(S.BP, need)
    cache = json.loads(GENES_FILE.read_text()) if GENES_FILE.exists() and not refresh else {}
    missing = [g for g in need if g not in cache and g not in PSEUDO]
    if missing:
        cache.update(fetch_genes(missing))
        GENES_FILE.write_text(json.dumps(dict(sorted(cache.items())), indent=1) + "\n")
    genes = dict(cache)
    genes.update(PSEUDO)
    for g, v in WHOLE_GENE.items():
        genes[g] = dict(genes.get(g, {}), **v)
    for g in genes.values():  # pseudoautosomal genes: Ensembl may report the Y copy; use X coordinates
        if g["chr"] == "Y" and g["end"] <= PAR1_X[1]:
            g["chr"] = "X"
    return genes, need


GENES = {}


def arm_dir(c, pos):
    """+1 if the telomere of this arm lies at higher coordinates (q arm), -1 for the p arm."""
    return 1 if pos >= CEN[c][2] else -1


def anchor(a):
    kind = a[0]
    if kind == "bp":
        return a[1], int(a[2])
    if kind in ("ts", "cs"):
        c, p = anchor(a[1])
        t = arm_dir(c, p) * (1 if kind == "ts" else -1)
        return c, int(p + t * a[2] * 1000)
    g = GENES[a[1]]
    s, e, strand = g["start"], g["end"], g["strand"]
    five, three = (s, e) if strand == 1 else (e, s)
    if kind == "5":
        return g["chr"], five
    if kind == "3":
        return g["chr"], three
    if kind == "s":
        return g["chr"], s
    if kind == "e":
        return g["chr"], e
    if kind == "f":
        return g["chr"], int(five + (three - five) * a[2])
    raise ValueError(a)


def region(r):
    kind = r[0]
    if kind == "gene":
        g = GENES[r[1]]
        return [[g["chr"], g["start"], g["end"]]]
    if kind == "span":
        return [[r[1], int(r[2]), int(r[3])]]
    if kind == "cen":
        c = r[1]
        mid = CEN[c][2]
        return [[c, mid - 600000, mid + 600000]]
    if kind == "het":
        c, band = r[1], r[2]
        b = [x for x in BANDS[c]["850"] if x[0] == band]
        if not b:
            raise ValueError(f"no band {c}{band}")
        s, e = b[0][3], b[0][4]
        mid = (s + e) // 2
        return [[c, max(s, mid - 1500000), min(e, mid + 1500000)]]
    if kind == "tel":
        c, arm = r[1], r[2]
        L = CHR_LEN[c]
        return [[c, 100000, 300000]] if arm == "p" else [[c, L - 300000, L - 100000]]
    if kind in ("ts", "cs"):
        c, p = anchor(r[1])
        t = arm_dir(c, p) * (1 if kind == "ts" else -1)
        a, b = p + t * r[2] * 1000, p + t * (r[2] + r[3]) * 1000
        lo, hi = int(min(a, b)), int(max(a, b))
        return [[c, max(1, lo), min(CHR_LEN[c], hi)]]
    raise ValueError(r)


def with_par(segs):
    """Pseudoautosomal regions of X are present on Y at the same sequence."""
    out = list(segs)
    for c, s, e in segs:
        if c == "X" and e <= PAR1_X[1]:
            out.append(["Y", s, e])
        if c == "X" and s >= PAR2_X[0]:
            out.append(["Y", s - PAR2_X[0] + PAR2_Y0, e - PAR2_X[0] + PAR2_Y0])
    return out


CEP_LOCUS = {"1": "D1Z5", "2": "D2Z1", "3": "D3Z1", "4": "D4Z1", "6": "D6Z1", "7": "D7Z1", "8": "D8Z2", "9": "D9Z1",
             "10": "D10Z1", "11": "D11Z1", "12": "D12Z3", "15": "D15Z4", "16": "D16Z3", "17": "D17Z1", "18": "D18Z1",
             "20": "D20Z1", "X": "DXZ1", "Y": "DYZ3"}

BA_GENE = {"v-rara-ba": "RARA", "v-cbfb-ba": "CBFB", "v-mll-ba": "KMT2A", "v-pdgfrb-ba": "PDGFRB", "v-etv6-ba": "ETV6",
           "v-myc-ba": "MYC", "v-igh-ba": "IGH", "v-bcl2-ba": "BCL2", "v-bcl6-ba": "BCL6", "v-ccnd1-ba": "CCND1",
           "v-malt1-ba": "MALT1", "v-tra-ba": "TRA", "v-alk-ba": "ALK", "v-ros1-ba": "ROS1", "v-ret-ba": "RET",
           "v-fgfr2-ba": "FGFR2", "v-ewsr1-ba": "EWSR1", "v-fus-ba": "FUS", "v-ddit3-ba": "DDIT3", "v-ss18-ba": "SS18",
           "v-foxo1-ba": "FOXO1", "v-ntrk1-ba": "NTRK1", "v-erg-ba": "ERG", "v-etv6-runx1-es": None,
           "c-runx1-ba": "RUNX1", "c-tcf3-ba": "TCF3", "c-crlf2-ba": "CRLF2", "m-jak2-ba": "JAK2", "m-nup98-ba": "NUP98",
           "m-pax5-ba": "PAX5", "m-abl2-ba": "ABL2", "m-abl1-ba": "ABL1", "m-p2ry8": "P2RY8", "c-igk-ba": "IGK",
           "c-igl-ba": "IGL", "c-tcl1-ba": "TCL1A", "c-trb-ba": "TRB", "c-tlx1-ba": "TLX1", "c-tlx3-ba": "TLX3"}
GENE_STRAND_FALLBACK = {"IGH": -1, "IGK": -1}
EXTRA_LOCUS = {"c-mecom-ba": ["LRRC34", "MECOM", "D3S1614"], "v-4q12": ["FIP1L1", "CHIC2", "PDGFRA"]}


def locus_name(p, i, c, segs):
    if c.get("locus"):
        return c["locus"]
    if p["id"] in EXTRA_LOCUS:
        return EXTRA_LOCUS[p["id"]][i]
    gene = BA_GENE.get(p["id"])
    if gene and len(p["comps"]) == 2:
        g = GENES.get(gene)
        if g:
            gmid, strand = (g["start"] + g["end"]) / 2, g["strand"]
        else:  # IG loci
            gmid = {"IGH": 106200000, "IGK": 89300000}[gene]
            strand = GENE_STRAND_FALLBACK[gene]
        cmid = sum((s + e) / 2 for _, s, e in segs) / len(segs)
        low = cmid < gmid
        five = low if strand == 1 else not low
        return ("5'" if five else "3'") + gene
    lab = c["label"]
    if lab.startswith("CEP "):
        return CEP_LOCUS.get(lab.split()[1], lab)
    first = lab.split(" (")[0].split()[0].rstrip(",")
    if "-" in first and first.split("-")[0] in GENES:
        first = first.split("-")[0]
    return {"HER-2/neu": "ERBB2", "TUPLE1": "HIRA", "KAL": "ANOS1", "LIS1": "PAFAH1B1", "SMS": "RAI1", "WHS": "NSD2",
            "N25": "D22S75", "D5S23,": "D5S721", "p58": "CDK11B", "1pTel": "1pTEL", "9q34": "ASS1",
            "D7S486,": "D7S486", "13q34": "LAMP1", "CHIC2": "CHIC2"}.get(first, first)


def ev(e):
    (k, v), = e.items()
    if k == "t":
        return {"t": [list(anchor(v[0])), list(anchor(v[1]))]}
    if k in ("del", "del2"):
        c = v[0]
        a = anchor(v[1]) if isinstance(v[1], tuple) else (c, v[1])
        b = anchor(v[2]) if isinstance(v[2], tuple) else (c, v[2])
        assert a[0] == b[0] == c, (e, a, b)
        return {k: [c, min(a[1], b[1]), max(a[1], b[1])]}
    if k == "delter":
        a = anchor(v[0]) if v[0][0] != "tel" else (v[0][1], region(v[0])[0][1])
        return {"delter": list(a)}
    if k == "inv":
        a, b = anchor(v[0]), anchor(v[1])
        assert a[0] == b[0]
        return {"inv": [a[0], min(a[1], b[1]), max(a[1], b[1])]}
    if k in ("gain", "loss", "armloss", "dup_der"):
        return {k: list(v)}
    if k in ("amp", "gainreg"):
        c, a, b = region(v[0])[0]
        pad = 500000 if k == "amp" else 0  # amplicons are typically several hundred kb to a few Mb
        return {k: [[c, max(1, a - pad), min(CHR_LEN[c], b + pad)], v[1]]}
    if k == "der":
        return {"der": [list(anchor(v[0])), list(anchor(v[1])), v[2]]}
    if k == "ins":
        return {"ins": [v[0], v[1][2], v[2], v[3][2], v[4][2]]}
    raise ValueError(e)


def main():
    global GENES
    genes, need = load_genes("--refresh" in sys.argv)
    GENES = genes
    out = []
    seen = set()
    for p in S.PROBES:
        if p["id"] in seen:
            raise SystemExit(f"duplicate id {p['id']}")
        seen.add(p["id"])
        comps = []
        for i, c in enumerate(p["comps"]):
            segs = []
            for r in c["region"]:
                segs += region(r)
            segs = with_par(segs)
            comps.append({"label": c["label"], "color": c["color"], "size": c["size"], "note": c["note"],
                          "set": c.get("set") or 1, "locus": locus_name(p, i, c, [s for s in segs if s[0] != "Y" or len(segs) == 1]),
                          "segs": segs})
        pats = []
        for pt in p["patterns"]:
            pats.append({"name": pt["name"], "events": [ev(e) for e in pt["events"]], "expect": pt["expect"],
                         "note": pt["note"], "src": pt["src"], "karyo": pt["karyo"]})
        out.append({k: p[k] for k in ("id", "vendor", "name", "group", "design", "codes", "status", "url", "note")} |
                   {"comps": comps, "patterns": pats})
    gene_tab = {g: [v["chr"], v["start"], v["end"], v["strand"]] for g, v in sorted(GENES.items()) if g in need}
    OUT.write_text(json.dumps({"_about": "Generated by tools/make_probes.py from tools/probe_specs.py. Positions are GRCh38 and approximate.",
                               "genes": gene_tab, "probes": out}, separators=(",", ":"), ensure_ascii=False) + "\n")
    print(f"wrote {OUT} ({len(out)} probes, {len(gene_tab)} genes)")


if __name__ == "__main__":
    main()
