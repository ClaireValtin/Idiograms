"""Curated FISH probe designs used by the FISH page (fish.html).

Every probe below was transcribed from the manufacturer's own material:

* Vysis (Abbott): the per-product "Results of hybridization" text and probe tables on
  molecular.abbott/int/en/chromosome-main/<n> (chromosomes 2-12, X, Y), the probe-map images
  served by molecular.abbott/api/assets/add/molecular/chromosome-probe-images/..., the
  Vysis CLL FISH Probe Kit package insert (30-608712/R4) and the 2019 US ordering list
  (AMD_US_Product_List.pdf). Sizes and gaps are the ones printed on the probe maps.
* CytoCell (OGT) and MetaSystems: the "Probe specification" / "Product description" text of
  each product page. They are only used for haematology targets that Vysis does not sell.

Positions are approximate: each probe is placed on GRCh38 from the genes and sizes on the
vendor map (see tools/make_probes.py). They are drawn to show where a probe lies, not to
replace the vendor's map.

Region helpers (all sizes in kb):
  G(gene)                    the whole gene
  R(chr, start, end)         explicit GRCh38 interval (bp)
  CEN(chr)                   centromeric alpha-satellite (or satellite III where stated)
  HET(chr, band)             a heterochromatin block, e.g. HET("9", "q12")
  TEL(chr, arm)              subtelomeric unique sequence (~100-300 kb from the end)
  TS(anchor, gap, size)      `size` kb on the telomeric side of `anchor`, starting `gap` kb away
  CS(anchor, gap, size)      the same on the centromeric side
  anchors: ("5", gene) 5' end, ("3", gene) 3' end, ("s", gene) lowest coordinate,
           ("e", gene) highest coordinate, ("f", gene, x) fraction x of the way 5'->3',
           ("bp", chr, pos)
Colours: O orange, G green, A aqua, Y gold, R red, B blue (CytoCell "blue" is an aqua dye).

Patterns are written as events applied to a normal diploid cell; the page derives the
signals (and the ISCN nuc ish string) from them. `expect` is the signal code the vendor
gives for that pattern; tools/check_probes.mjs fails the build if the derived code differs.
"""

def G(g): return ("gene", g)
def R(c, s, e): return ("span", c, s, e)
def CEN(c, kind="alpha"): return ("cen", c, kind)
def HET(c, band): return ("het", c, band)
def TEL(c, arm): return ("tel", c, arm)
def TS(anchor, gap, size): return ("ts", anchor, gap, size)
def CS(anchor, gap, size): return ("cs", anchor, gap, size)
def a5(g): return ("5", g)
def a3(g): return ("3", g)
def af(g, x): return ("f", g, x)
def bp(c, p): return ("bp", c, p)

# ---------------------------------------------------------------- breakpoints
# Typical breakpoints, as fractions along the gene (5'->3') or explicit positions.
BP = {
    "ABL1": af("ABL1", 0.45),          # intron 1, between exons 1b and a2
    "BCR_M": af("BCR", 0.80),          # M-bcr (exons 12-16)
    "BCR_m": af("BCR", 0.25),          # m-bcr (intron 1)
    "PML": af("PML", 0.70),            # bcr1, intron 6
    "RARA": af("RARA", 0.30),          # intron 2
    "RUNX1_t821": af("RUNX1", 0.68),   # intron 5 (t(8;21))
    "RUNX1T1": af("RUNX1T1", 0.10),    # intron 1
    "RUNX1_t1221": af("RUNX1", 0.30),  # large intron 1 (t(12;21))
    "ETV6": af("ETV6", 0.78),          # intron 5
    "CBFB": af("CBFB", 0.85),          # intron 5
    "MYH11": af("MYH11", 0.45),
    "KMT2A": af("KMT2A", 0.38),        # breakpoint cluster region, exons 9-11
    "AFF1": af("AFF1", 0.25),
    "MLLT3": af("MLLT3", 0.30),
    "MLLT1": af("MLLT1", 0.40),
    "AFDN": af("AFDN", 0.10),
    "MECOM_inv": TS(a5("MECOM"), 60, 0),   # inv(3): upstream (5') of MECOM
    "GATA2": TS(a5("GATA2"), 40, 0),       # 3q21.3, between GATA2 and RPN1 (distal GATA2 enhancer)
    "DEK": af("DEK", 0.70),
    "NUP214": af("NUP214", 0.75),
    "TCF3": af("TCF3", 0.85),
    "PBX1": af("PBX1", 0.25),
    "IGH": ("bp", "14", 105865000),    # IGHJ region
    "IGH_sw": ("bp", "14", 105700000), # switch region (myeloma translocations)
    "MYC": CS(a5("MYC"), 20, 0),       # just 5' of MYC (t(8;14))
    "BCL2": af("BCL2", 0.99),          # MBR, 3' UTR
    "CCND1": CS(a5("CCND1"), 110, 0),  # MTC, ~110 kb centromeric of CCND1
    "CCND3": CS(a5("CCND3"), 10, 0),
    "FGFR3_NSD2": af("NSD2", 0.05),
    "MAF": CS(a3("MAF"), 700, 0),      # in WWOX, centromeric (3') of MAF
    "MAFB": CS(a3("MAFB"), 300, 0),
    "MALT1": af("MALT1", 0.30),
    "BIRC3": af("BIRC3", 0.70),
    "BCL6": af("BCL6", 0.15),          # MBR, 5' of BCL6
    "PDGFRB": af("PDGFRB", 0.40),
    "JAK2": af("JAK2", 0.70),
    "PCM1": af("PCM1", 0.50),
    "NUP98": af("NUP98", 0.40),
    "NSD1": af("NSD1", 0.30),
    "EWSR1": af("EWSR1", 0.55),
    "FLI1": af("FLI1", 0.40),
    "ALK": af("ALK", 0.75),
    "EML4": af("EML4", 0.40),
    "ROS1": af("ROS1", 0.70),
    "CD74": af("CD74", 0.60),
    "RET": af("RET", 0.45),
    "KIF5B": af("KIF5B", 0.60),
    "SS18": af("SS18", 0.90),
    "SSX1": af("SSX1", 0.40),
    "FUS": af("FUS", 0.40),
    "DDIT3": af("DDIT3", 0.10),
    "FOXO1": af("FOXO1", 0.30),
    "PAX3": af("PAX3", 0.80),
    "TRA": af("TRA", 0.50),
    "TCL1A": TS(a5("TCL1A"), 5, 0),
    "CRLF2": af("CRLF2", 0.50),
    "P2RY8": af("P2RY8", 0.30),
    "PAX5": af("PAX5", 0.20),
    "ABL2": af("ABL2", 0.30),
    "TCF3_BA": af("TCF3", 0.85),
    "TLX1": TS(a3("TLX1"), 30, 0),
    "TLX3": TS(a3("TLX3"), 30, 0),
    "TRB": af("TRB", 0.50),
    "IGK": CS(("e", "IGKC"), 10, 0),
    "IGL": ("bp", "22", 22900000),
    "NTRK1": af("NTRK1", 0.55),
    "TPM3": af("TPM3", 0.50),
    "ERG": af("ERG", 0.40),
    "TMPRSS2": af("TMPRSS2", 0.10),
    "FGFR2": af("FGFR2", 0.95),
    "ETV1": af("ETV1", 0.50),
}

# Common rearrangements as event lists (positions resolved through BP).
def T(a, b): return {"t": [a, b]}                 # balanced translocation between two breakpoints
def DEL(c, s, e): return {"del": [c, s, e]}       # interstitial deletion of one homologue (anchors or bp)
def DELTER(anchor, side): return {"delter": [anchor, side]}  # terminal deletion: lose everything on `side` ('tel'/'cen'... 'tel' only)
def INV(a, b): return {"inv": [a, b]}             # paracentric/pericentric inversion of one homologue
def GAIN(c, n=1): return {"gain": [c, n]}
def LOSS(c): return {"loss": [c]}
def AMP(region, n): return {"amp": [region, n]}   # extra copies of a region (clustered)
def ISO(c, arm): return {"iso": [c, arm]}         # isochromosome: lose the other arm, gain this one
def DELBOTH(c, s, e): return {"del2": [c, s, e]}  # biallelic interstitial deletion
def GAINREG(region, n=1): return {"gainreg": [region, n]}  # extra copies of a region as separate small chromosomes
def ARMLOSS(c, arm): return {"armloss": [c, arm]}  # one homologue loses a whole arm

def pat(name, events, expect=None, note=None, src="standard", karyo=None):
    return {"name": name, "events": events, "expect": expect, "note": note, "src": src, "karyo": karyo}

NORMAL = lambda expect=None, note=None, src="vendor": pat("Normal", [], expect, note, src)

V = "Vysis (Abbott)"
C = "CytoCell (OGT)"
M = "MetaSystems"

def P(id, vendor, name, group, design, comps, patterns, codes=None, status=None, url=None, note=None, uses=None):
    return {"id": id, "vendor": vendor, "name": name, "group": group, "design": design,
            "comps": comps, "patterns": patterns, "codes": codes or [], "status": status or [],
            "url": url, "note": note, "uses": uses}

def comp(label, color, region, size=None, note=None, set=None, locus=None):
    return {"label": label, "color": color, "region": region if isinstance(region, list) else [region],
            "size": size, "note": note, "set": set, "locus": locus}

ABB_INT = "https://www.molecular.abbott/int/en/chromosome-main/"

# ------------------------------------------------------------------ components
# Vysis components reused by several kits (from the probe maps).
V_ABL1_ASS1_O = comp("ABL1-ASS1", "O", TS(a5("ASS1"), -60, 671), "~671 kb", "ASS1 to beyond the 3' end of ABL1")
V_ASS1_A = comp("ASS1", "A", TS(a5("ASS1"), -60, 329), "~329 kb", "ASS1 to ABL1 exon 1b")
V_BCR_DF_G = comp("BCR", "G", [CS(a3("BCR"), 0, 579), TS(a3("BCR"), 326, 645)], "~579 + ~645 kb", "IGLV to the 3' end of BCR, ~326 kb gap, then a telomeric segment")
V_BCR_G = comp("BCR", "G", TS(a5("BCR"), -140, 287), "~287 kb", "spans BCR (m-bcr and M-bcr)")
V_ABL1_SF_O = comp("ABL1", "O", TS(af("ABL1", 0.62), 0, 190), "~190 kb", "ABL1 exon 5 to beyond the 3' end")
V_IGH_G = comp("IGH", "G", R("14", 105400000, 107000000), "~1.6 Mb", "whole IGH locus: constant, J, D and variable segments")
C_IGH_G = comp("IGH", "G", [R("14", 105550000, 105726000), R("14", 105726000, 105896000), R("14", 106270000, 106594000), R("14", 106594000, 106883000)],
                "176 + 170 + 324 + 289 kb", "constant, J and D segments, and the variable segments")
V_IGH_G2 = comp("IGH", "G", [R("14", 105580000, 105918000), R("14", 106395000, 106832000)], "~338 + ~437 kb", "two segments with a ~477 kb gap")

def cep(c, color, label=None, size=None, kind="alpha"):
    names = {"1": "D1Z5", "2": "D2Z1", "3": "D3Z1", "4": "D4Z1", "6": "D6Z1", "7": "D7Z1", "8": "D8Z2",
             "9": "D9Z1", "10": "D10Z1", "11": "D11Z1", "12": "D12Z3", "15": "D15Z4", "16": "D16Z3",
             "17": "D17Z1", "18": "D18Z1", "20": "D20Z1", "X": "DXZ1", "Y": "DYZ3"}
    return comp(label or f"CEP {c}", color, CEN(c, kind), size or "alpha satellite", names.get(c, ""))

PROBES = []
add = PROBES.append

# ================================================================== MYELOID
add(P("v-bcr-abl-df", V, "Vysis LSI BCR/ABL Dual Color, Dual Fusion Translocation Probe", "Myeloid", "DF",
      [V_ABL1_ASS1_O, V_BCR_DF_G],
      [NORMAL("2O2G"),
       pat("t(9;22)", [T(BP["ABL1"], BP["BCR_M"])], "1O1G2F", "Simple balanced t(9;22): one fusion on the der(9), one on the der(22).", "vendor", "t(9;22)(q34.1;q11.2)"),
       pat("t(9;22), der(9) deletion", [T(BP["ABL1"], BP["BCR_M"]), DEL("9", TS(a5("ASS1"), -60, 0), BP["ABL1"]), DEL("22", BP["BCR_M"], TS(a3("BCR"), 980, 0))], "1O1G1F", "Deletion of sequences 5' of the ABL1 and 3' of the BCR breakpoint on the der(9) leaves a single fusion.", "vendor"),
       pat("t(9;22) + extra Ph", [T(BP["ABL1"], BP["BCR_M"]), {"dup_der": ["22"]}], "1O1G3F", "Duplicated der(22): a common finding in progression.", "standard")],
      ["05J82-001 (ASR, 20 µL)", "05J82-010 (ASR, 50 µL)", "08L10-001 (CE kit)"], ["ASR", "CE"], ABB_INT + "9"))

add(P("v-bcr-abl-sf", V, "Vysis LSI BCR/ABL Dual Color, Single Fusion Translocation Probe", "Myeloid", "SF",
      [V_ABL1_SF_O, V_BCR_G],
      [NORMAL("2O2G"),
       pat("t(9;22)", [T(BP["ABL1"], BP["BCR_M"])], "1O1G1F", None, "vendor", "t(9;22)(q34.1;q11.2)")],
      ["05J77-001 (ASR)", "08L56-050 (CE kit)"], ["ASR", "CE"], ABB_INT + "9"))

add(P("v-bcr-abl-es", V, "Vysis LSI BCR/ABL ES Dual Color Translocation Probe", "Myeloid", "ES",
      [V_ABL1_ASS1_O, V_BCR_G],
      [NORMAL("2O2G"),
       pat("t(9;22), M-bcr", [T(BP["ABL1"], BP["BCR_M"])], "2O1G1F", "One large orange (native ABL1), one smaller orange (the 'extra signal' on the der(9)), one green and one fusion.", "vendor", "t(9;22)(q34.1;q11.2)"),
       pat("t(9;22), m-bcr", [T(BP["ABL1"], BP["BCR_m"])], "1O1G2F", "Minor breakpoint: the BCR probe is split, giving two fusions.", "vendor")],
      ["08L55-020 (CE kit)"], ["CE"], ABB_INT + "9"))

add(P("v-bcr-abl-ass1", V, "Vysis BCR/ABL1/ASS1 Tri-Color DF FISH Probe", "Myeloid", "TC-DF",
      [V_ABL1_ASS1_O, V_ASS1_A, V_BCR_DF_G],
      [NORMAL("2AO2G", "Two aqua/orange signals and two green."),
       pat("t(9;22)", [T(BP["ABL1"], BP["BCR_M"])], "1AO1G1OG1AOG", "Without deletion the der(9) keeps ASS1 (aqua) next to its orange/green fusion.", "standard", "t(9;22)(q34.1;q11.2)"),
       pat("t(9;22), der(9) deletion", [T(BP["ABL1"], BP["BCR_M"]), DEL("9", TS(a5("ASS1"), -60, 0), BP["ABL1"]), DEL("22", BP["BCR_M"], TS(a3("BCR"), 980, 0))], "1AO1G1OG", "Loss of the aqua ASS1 signal marks the der(9) deletion.", "vendor")],
      ["05N54-020 (CE, 20 µL)", "05N54-050 (CE, 50 µL)"], ["CE"], ABB_INT + "9"))

add(P("v-9q34", V, "Vysis LSI 9q34 SpectrumAqua Probe", "Myeloid", "LSI",
      [comp("9q34 (ASS1)", "A", TS(a5("ASS1"), -60, 329), "~329 kb", "ASS1 to ABL1 exon 1b")],
      [NORMAL("2A"), pat("9q34 deletion", [DEL("9", TS(a5("ASS1"), -100, 0), TS(a5("ASS1"), 300, 0))], "1A", None, "vendor")],
      ["05J79-011 (ASR)", "05N53-020 (CE kit)"], ["ASR", "CE"], ABB_INT + "9"))

add(P("v-pml-rara-df", V, "Vysis LSI PML/RARA Dual Color, Dual Fusion Translocation Probe", "Myeloid", "DF",
      [comp("PML", "O", [CS(a3("PML"), 0, 239), TS(a3("PML"), 15, 278)], "~239 + ~278 kb", "spans PML with a ~15 kb gap"),
       comp("RARA", "G", [CS(a3("RARA"), 0, 410), TS(a3("RARA"), 8, 417)], "~410 + ~417 kb", "THRA to RARA, ~8 kb gap, then TOP2A to CCR7")],
      [NORMAL("2O2G", src="standard"),
       pat("t(15;17)", [T(BP["PML"], BP["RARA"])], "1O1G2F", None, "standard", "t(15;17)(q24.1;q21.2)")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None))

add(P("v-pml-rara-sf", V, "Vysis LSI PML/RARA Dual Color, Single Fusion Probes", "Myeloid", "SF",
      [comp("PML", "O", CS(a3("PML"), 0, 239), "~239 kb", "SHGC-144626 to PML exon 7"),
       comp("RARA", "G", TS(a3("RARA"), 0, 417), "~417 kb", "3' RARA to D17S1228")],
      [NORMAL("2O2G", src="standard"),
       pat("t(15;17)", [T(BP["PML"], BP["RARA"])], "1O1G1F", "Fusion on the der(15) only.", "standard", "t(15;17)(q24.1;q21.2)")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None))

add(P("v-rara-ba", V, "Vysis LSI RARA Dual Color Break Apart Rearrangement Probe", "Myeloid", "BA",
      [comp("5' RARA", "O", CS(af("RARA", 0.25), 0, 150), "~150 kb", "D17S1238 to RARA 5'"),
       comp("3' RARA", "G", TS(af("RARA", 0.25), 0, 417), "~417 kb", "RARA 3', TOP2A, CCR7")],
      [NORMAL("2F", src="standard"),
       pat("RARA rearranged", [T(BP["RARA"], BP["PML"])], "1O1G1F", None, "standard", "t(15;17)(q24.1;q21.2)")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None))

add(P("v-runx1-runx1t1", V, "Vysis LSI RUNX1/RUNX1T1 DF FISH Probe Kit", "Myeloid", "DF",
      [comp("RUNX1T1", "O", TS(a3("RUNX1T1"), -150, 655), "~655 kb", "covers RUNX1T1"),
       comp("RUNX1", "G", TS(a3("RUNX1"), -150, 1400), "~1.4 Mb", "covers RUNX1")],
      [NORMAL("2O2G"),
       pat("t(8;21)", [T(BP["RUNX1T1"], BP["RUNX1_t821"])], "1O1G2F", None, "vendor", "t(8;21)(q21.3;q22.1)")],
      ["08L70-020 (CE kit)"], ["CE"], ABB_INT + "8"))

add(P("v-cbfb-ba", V, "Vysis LSI CBFB Dual Color Break Apart Rearrangement Probe", "Myeloid", "BA",
      [comp("5' CBFB", "R", CS(af("CBFB", 0.70), 22, 130), "~130 kb", "centromeric, to CBFB exon 1"),
       comp("3' CBFB", "G", TS(af("CBFB", 0.70), 0, 204), "~204 kb", "CBFB 3' to D16S318")],
      [NORMAL("2F", src="standard"),
       pat("inv(16)", [INV(BP["MYH11"], BP["CBFB"])], "1R1G1F", None, "standard", "inv(16)(p13.1q22.1)")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None,
      "Partner-independent: a split also occurs with rare CBFB partners other than MYH11."))

add(P("v-rpn1-mecom", V, "Vysis LSI RPN1/MECOM DF FISH Probe Kit", "Myeloid", "DF",
      [comp("MECOM", "O", TS(a3("MECOM"), -400, 1200), "~1.2 Mb", "D3S1243 to WI-3520, spans MECOM"),
       comp("RPN1", "G", TS(a5("GATA2"), -200, 580), "~580 kb", "EEFSEC, GATA2, RPN1")],
      [NORMAL("2O2G", "The two probes are on the same arm and may appear fused in some normal nuclei."),
       pat("inv(3)", [INV(BP["GATA2"], BP["MECOM_inv"])], "1O1G2F", None, "vendor", "inv(3)(q21.3q26.2)")],
      ["06N60-010 (CE kit)"], ["CE"], ABB_INT + "3"))

add(P("v-dek-nup214", V, "Vysis LSI DEK/NUP214 Dual Color, Dual Fusion Translocation FISH Probe Kit", "Myeloid", "DF",
      [comp("NUP214", "O", TS(a5("ABL1"), -20, 871), "~871 kb", "ABL1 to POMT1, spans NUP214"),
       comp("DEK", "G", TS(a5("DEK"), -500, 951), "~951 kb", "KIF13A to DEK")],
      [NORMAL("2O2G", src="standard"),
       pat("t(6;9)", [T(BP["DEK"], BP["NUP214"])], "1O1G2F", None, "standard", "t(6;9)(p22.3;q34.1)")],
      ["RUO"], ["RUO"], ABB_INT + "6"))

add(P("v-mll-ba", V, "Vysis LSI MLL Dual Color, Break Apart Rearrangement Probe", "Myeloid", "BA",
      [comp("5' KMT2A", "G", CS(BP["KMT2A"], 0, 384), "~384 kb", "SHGC-85140, SCN4B, CD3D, into KMT2A"),
       comp("3' KMT2A", "O", TS(BP["KMT2A"], 0, 256), "~256 kb", "KMT2A 3', ARCN1, TREH")],
      [NORMAL("2F", src="standard"),
       pat("KMT2A rearranged, t(9;11)", [T(BP["MLLT3"], BP["KMT2A"])], "1O1G1F", None, "vendor", "t(9;11)(p21.3;q23.3)"),
       pat("KMT2A rearranged, 3' deletion", [T(BP["MLLT3"], BP["KMT2A"]), DEL("11", BP["KMT2A"], TS(BP["KMT2A"], 300, 0))], "1G1F", "Loss of the translocated 3' part: one fusion and one green only.", "standard")],
      ["08L57-020 (CE)"], ["CE"], ABB_INT + "11"))

add(P("v-egr1", V, "Vysis LSI EGR1/D5S23, D5S721 Dual Color Probe Kit", "Myeloid", "LSI+control",
      [comp("EGR1 (5q31)", "O", TS(a5("EGR1"), -100, 209), "~209 kb", None),
       comp("D5S23, D5S721 (5p15.2)", "G", R("5", 9900000, 10461000), "~561 kb", "control")],
      [NORMAL("2O2G"),
       pat("del(5q)", [DEL("5", bp("5", 98000000), bp("5", 158000000))], "1O2G", None, "vendor", "del(5)(q21q33)"),
       pat("monosomy 5", [LOSS("5")], "1O1G", None, "standard")],
      ["08L68-020 (CE kit)", "05J76-001 (ASR, Cri-du-Chat region)"], ["CE", "ASR"], ABB_INT + "5"))

add(P("v-csf1r", V, "Vysis LSI CSF1R/D5S23, D5S721 FISH Probe Kit", "Myeloid", "LSI+control",
      [comp("CSF1R (5q33-34)", "O", TS(a5("SLC26A2"), -20, 322), "~322 kb", "SLC26A2, CSF1R, PDGFRB"),
       comp("D5S23, D5S721 (5p15.2)", "G", R("5", 9900000, 10461000), "~561 kb", "control")],
      [NORMAL("2O2G"), pat("del(5q)", [DEL("5", bp("5", 98000000), bp("5", 158000000))], "1O2G", None, "vendor", "del(5)(q21q33)")],
      ["05N03-020 (CE kit)", "05J60-001 (ASR)"], ["CE", "ASR"], ABB_INT + "5"))

add(P("v-d7s486", V, "Vysis LSI D7S486/CEP 7 FISH Probe Kit", "Myeloid", "LSI+CEP",
      [comp("D7S486 (7q31)", "O", TS(a3("TES"), 20, 308), "~308 kb", None), cep("7", "G")],
      [NORMAL("2O2G"),
       pat("del(7q)", [DEL("7", bp("7", 100000000), bp("7", 140000000))], "1O2G", None, "standard", "del(7)(q22q34)"),
       pat("monosomy 7", [LOSS("7")], "1O1G", None, "standard")],
      ["05N07-020 (CE kit)", "05J61-001 (ASR)"], ["CE", "ASR"], ABB_INT + "7"))

add(P("v-d7s522", V, "Vysis LSI D7S522/CEP 7 FISH Probe Kit", "Myeloid", "LSI+CEP",
      [comp("D7S522 (7q31)", "O", TS(a3("TES"), 230, 224), "~224 kb", None), cep("7", "G")],
      [NORMAL("2O2G"), pat("del(7q)", [DEL("7", bp("7", 100000000), bp("7", 140000000))], "1O2G", None, "vendor", "del(7)(q22q34)")],
      ["05N08-020 (CE kit)", "05J85-001 (ASR)"], ["CE", "ASR"], ABB_INT + "7"))

add(P("v-d20s108", V, "Vysis LSI D20S108 SpectrumOrange Probe", "Myeloid", "LSI",
      [comp("D20S108 (20q12)", "O", TS(af("PTPRT", 0.55), 0, 201), "~201 kb", None)],
      [NORMAL("2O", src="standard"), pat("del(20q)", [DEL("20", bp("20", 38000000), bp("20", 48000000))], "1O", None, "standard", "del(20)(q11.2q13.1)")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None))

add(P("v-4q12", V, "Vysis 4q12 Tri-Color Rearrangement FISH Probe Kit", "Myeloid", "TC-BA",
      [comp("FIP1L1 side", "G", CS(a5("FIP1L1"), 0, 703), "~703 kb", "SCFD2 to the 5' end of FIP1L1"),
       comp("CHIC2 region", "O", CS(("e", "CHIC2"), -10, 448), "~448 kb", "LNX1 to CHIC2"),
       comp("PDGFRA/KIT", "A", TS(("s", "PDGFRA"), -20, 578), "~578 kb", "PDGFRA, KIT")],
      [NORMAL("2GOA", "Two tri-colour (green/orange/aqua) fusions."),
       pat("FIP1L1::PDGFRA (CHIC2 deletion)", [DEL("4", af("FIP1L1", 0.55), af("PDGFRA", 0.45))], "1GOA1GA", "One tri-colour fusion and one green/aqua fusion without orange.", "vendor"),
       pat("PDGFRA translocation", [T(af("PDGFRA", 0.45), af("ETV6", 0.5))], "1GOA1GO1A", "One tri-colour fusion, one orange/green fusion and a separate aqua signal.", "vendor")],
      ["05N52-020 (CE kit)", "01N79-020 (ASR)"], ["CE", "ASR"], ABB_INT + "4"))

add(P("v-pdgfrb-ba", V, "Vysis LSI PDGFRB Break Apart FISH Probe Kit", "Myeloid", "BA",
      [comp("PDGFRB (Cen)", "O", CS(af("PDGFRB", 0.95), 0, 351), "~351 kb", "D5S1907 to CSF1R"),
       comp("PDGFRB (Tel)", "G", TS(a5("PDGFRB"), 10, 368), "~368 kb", "D5S2015, TCOF1, D5S551")],
      [NORMAL("2F"), pat("PDGFRB rearranged", [T(BP["PDGFRB"], af("ETV6", 0.5))], "1O1G1F", None, "vendor", "t(5;12)(q32;p13.2)")],
      ["06N24-010 (CE kit)"], ["CE"], ABB_INT + "5"))

add(P("v-fgfr1-cep8", V, "Vysis 8p12 FGFR1 SpectrumRed/CEP 8 SpectrumAqua FISH", "Solid tumour", "LSI+CEP",
      [comp("FGFR1 (8p12)", "R", TS(a5("FGFR1"), -280, 531), "~531 kb", "SHGC-56788 to WI-8882"), cep("8", "A")],
      [NORMAL("2R2A", src="standard"), pat("FGFR1 amplification", [AMP(G("FGFR1"), 12)], "2R2A amp R", None, "standard")],
      ["08N21-060 (RUO)"], ["RUO"], ABB_INT + "8"))

add(P("v-cep8", V, "Vysis CEP 8 SpectrumOrange DNA Probe Kit", "Myeloid", "CEP",
      [cep("8", "O")],
      [NORMAL("2O", src="standard"), pat("trisomy 8", [GAIN("8")], "3O", None, "standard", "+8")],
      ["CE kit (with/without control slides)", "CEP 8 (D8Z2) also as SpectrumAqua and SpectrumGreen ASR"], ["CE", "ASR"], ABB_INT + "8"))

add(P("v-cep12", V, "Vysis CEP 12 SpectrumOrange DNA Probe Kit", "CLL", "CEP",
      [cep("12", "O")],
      [NORMAL("2O", src="standard"), pat("trisomy 12", [GAIN("12")], "3O", None, "standard", "+12")],
      ["CE kit (with/without control slides)", "CEP 12 (D12Z3) also as SpectrumGreen ASR"], ["CE", "ASR"], ABB_INT + "12"))

add(P("v-cepxy", V, "Vysis CEP X SpectrumOrange/Y SpectrumGreen DNA Probe Kit", "Constitutional", "CEP",
      [cep("X", "O"), comp("CEP Y (Yq12)", "G", HET("Y", "q12"), "satellite III", "DYZ1", locus="DYZ1")],
      [NORMAL("2O", "Female (46,XX): two orange."),
       pat("Male", [LOSS("X"), GAIN("Y")], "1O1G", "Normal male: one orange, one green.", "vendor"),
       pat("47,XXY", [GAIN("Y")], "2O1G", None, "standard")],
      ["07J22-050 (CE kit)"], ["CE"], ABB_INT + "x",
      "Used for sex-mismatched transplant monitoring and sex chromosome enumeration."))

# ================================================================== LYMPHOID / ALL
add(P("v-etv6-ba", V, "Vysis ETV6 Break Apart FISH Probe Kit", "Lymphoid", "BA",
      [comp("ETV6 (Tel)", "O", TS(af("ETV6", 0.70), 0, 486), "~486 kb", "RH71401 to ETV6 (5' side)"),
       comp("ETV6 (Cen)", "G", CS(af("ETV6", 0.70), 140, 632), "~632 kb", "BCL2L14, DUSP16")],
      [NORMAL("2F"), pat("ETV6 rearranged", [T(BP["ETV6"], BP["RUNX1_t1221"])], "1O1G1F", None, "vendor", "t(12;21)(p13.2;q22.1)")],
      ["04N09-020 (CE kit)"], ["CE"], ABB_INT + "12"))

add(P("v-etv6-runx1-es", V, "Vysis LSI ETV6 (TEL)/RUNX1 (AML1) ES Dual Color Translocation Probe", "Lymphoid", "ES",
      [comp("ETV6", "G", TS(af("ETV6", 0.82), 0, 347), "~347 kb", "UT6136 to RH7051: telomeric flank and the 5' part of ETV6"),
       comp("RUNX1", "O", TS(a3("RUNX1"), -180, 732), "~732 kb", "SHGC-52032 to SHGC-110084, whole RUNX1")],
      [NORMAL("2O2G", src="standard"),
       pat("t(12;21)", [T(BP["ETV6"], BP["RUNX1_t1221"])], "2O1G1F", "One smaller orange: the part of RUNX1 left on the der(12).", "standard", "t(12;21)(p13.2;q22.1)"),
       pat("t(12;21) + loss of other ETV6", [T(BP["ETV6"], BP["RUNX1_t1221"]), DEL("12", TS(a5("ETV6"), 300, 0), CS(a3("ETV6"), 300, 0))], "2O1F", "Deletion of the non-translocated ETV6 is common in ETV6::RUNX1 ALL.", "standard"),
       pat("iAMP21", [AMP(G("RUNX1"), 5)], "2G2O amp O", "Five or more RUNX1 signals, usually clustered on one abnormal chromosome 21.", "standard")],
      ["05J62-001 (ASR)", "08L66-020 (CE)"], ["ASR", "CE"], ABB_INT + "12"))

add(P("v-etv6-runx1-df", V, "Vysis ETV6/RUNX1 DF FISH Probe Kit", "Lymphoid", "DF",
      [comp("ETV6", "O", TS(a5("ETV6"), -600, 1300), "~1.3 Mb", "GDB-214813 to SHGC-83859, whole ETV6"),
       comp("RUNX1", "G", TS(a3("RUNX1"), -150, 1400), "~1.4 Mb", "whole RUNX1")],
      [NORMAL("2O2G", src="standard"),
       pat("t(12;21)", [T(BP["ETV6"], BP["RUNX1_t1221"])], "1O1G2F", None, "standard", "t(12;21)(p13.2;q22.1)")],
      ["05N96-010 (CE kit)"], ["CE"], ABB_INT + "12"))

add(P("v-tcf3-pbx1", V, "Vysis LSI TCF3 SpectrumGreen + LSI PBX1 SpectrumOrange Probes", "Lymphoid", "DF",
      [comp("TCF3", "G", TS(a5("TCF3"), -350, 730), "~730 kb", "RH55657 to RH75938, whole TCF3"),
       comp("PBX1", "O", TS(a5("PBX1"), -150, 635), "~635 kb", "SHGC-81164 to SHGC-147938, whole PBX1")],
      [NORMAL("2O2G", src="standard"),
       pat("t(1;19)", [T(BP["PBX1"], BP["TCF3"])], "1O1G2F", "Both probes span the breakpoints, so the balanced t(1;19) gives two fusions.", "standard", "t(1;19)(q23.3;p13.3)"),
       pat("der(19)t(1;19)", [{"der": [BP["PBX1"], BP["TCF3"], "19"]}], "2O1G1F", "Unbalanced form: der(19) with two normal chromosomes 1.", "standard", "der(19)t(1;19)(q23.3;p13.3)")],
      ["Sold as two separate ASR probes"], ["ASR"], None))

add(P("v-d5-cep9-cep15", V, "Vysis LSI D5S23, D5S721/CEP 9/CEP 15 FISH Probe Kit", "Lymphoid", "CEP",
      [comp("D5S23, D5S721", "G", R("5", 9900000, 10461000), "~561 kb", None), cep("9", "A"), cep("15", "O")],
      [NORMAL("2G2A2O"), pat("High hyperdiploid", [GAIN("9"), GAIN("15")], "2G3A3O", None, "vendor")],
      ["05N35-020 (CE kit)"], ["CE"], ABB_INT + "5"))

add(P("v-myc-ba", V, "Vysis LSI MYC Dual Color Break Apart Rearrangement Probe", "Lymphoid", "BA",
      [comp("5' MYC", "O", CS(a5("MYC"), 119, 277), "~277 kb", "WI-1203 to just 5' of MYC"),
       comp("3' MYC", "G", TS(a3("MYC"), 1500, 407), "~407 kb", "AFMA084ZD1, D8S1207, AFMB335XE9")],
      [NORMAL("2F"),
       pat("t(8;14)", [T(TS(a3("MYC"), 400, 0), BP["IGH"])], "1O1G1F", "Breakpoint within the ~1.5 Mb gap.", "vendor", "t(8;14)(q24.2;q32.3)")],
      ["01N63-020 (CE kit)", "05J91-001 (ASR)"], ["CE", "ASR"], ABB_INT + "8",
      "Breakpoints outside the gap, common with IGL/IGK partners, can be missed."))

add(P("v-igh-myc-cep8", V, "Vysis LSI IGH/MYC/CEP 8 Tri-Color Dual Fusion Probe Kit", "Lymphoid", "TC-DF",
      [comp("MYC", "O", TS(a5("MYC"), -450, 821), "~821 kb", "RH48719 to SHGC-9463, MYC and PVT1"), V_IGH_G, cep("8", "A")],
      [NORMAL("2O2G2A"),
       pat("t(8;14)", [T(BP["MYC"], BP["IGH"])], "1O1G2F2A", None, "standard", "t(8;14)(q24.2;q32.3)"),
       pat("t(8;14) + trisomy 8", [T(BP["MYC"], BP["IGH"]), GAIN("8")], "2O1G2F3A", None, "standard")],
      ["04N10-020 (CE kit)", "05J75-001 (ASR)"], ["CE", "ASR"], ABB_INT + "8"))

add(P("v-igh-ba", V, "Vysis LSI IGH Dual Color, Break Apart Rearrangement Probe", "Lymphoid", "BA",
      [comp("3' IGH (constant)", "O", R("14", 105300000, 105582000), "~282 kb", "3' flanking, centromeric"),
       comp("5' IGH (variable)", "G", R("14", 106076000, 106903000), "~827 kb", "variable segments")],
      [NORMAL("2F", src="standard"),
       pat("IGH rearranged, t(14;18)", [T(BP["IGH"], BP["BCL2"])], "1O1G1F", None, "standard", "t(14;18)(q32.3;q21.3)")],
      ["US ordering list 2019 (ASR)", "CE kit"], ["ASR", "CE"], None))

add(P("v-igh-bcl2", V, "Vysis LSI IGH/BCL2 Dual Color Dual Fusion Probes", "Lymphoid", "DF",
      [comp("BCL2", "O", TS(a3("BCL2"), -250, 750), "~750 kb", "D18S91 to D18S87, whole BCL2"), V_IGH_G],
      [NORMAL("2O2G", src="standard"), pat("t(14;18)", [T(BP["IGH"], BP["BCL2"])], "1O1G2F", None, "standard", "t(14;18)(q32.3;q21.3)")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None))

add(P("v-bcl2-ba", V, "Vysis LSI BCL2 Break Apart FISH Probe Kit", "Lymphoid", "BA",
      [comp("BCL2 (Cen)", "G", CS(a3("BCL2"), 0, 855), "~855 kb", "TNFRSF11A, ZCCHC2"),
       comp("BCL2 (Tel)", "O", TS(a3("BCL2"), 260, 608), "~608 kb", "BCL2 5', KDSR, SERPINB2")],
      [NORMAL("2F", src="standard"), pat("BCL2 rearranged", [T(BP["IGH"], af("BCL2", 0.5))], "1O1G1F", None, "standard", "t(14;18)(q32.3;q21.3)")],
      ["US ordering list 2019"], ["ASR"], None))

add(P("v-bcl6-ba", V, "Vysis LSI BCL6 (ABR) Dual Color Break Apart Rearrangement Probe", "Lymphoid", "BA",
      [comp("BCL6 (Cen)", "G", CS(a5("BCL6"), 0, 816), "~816 kb", "ST6GAL1, MASP1, BCL6 to the major breakpoint region"),
       comp("BCL6 (Tel)", "O", TS(a5("BCL6"), 265, 349), "~349 kb", "telomeric of the alternative breakpoint region (ABR), ~265 kb gap")],
      [NORMAL("2F", src="standard"), pat("BCL6 rearranged", [T(BP["BCL6"], BP["IGH"])], "1O1G1F", None, "vendor", "t(3;14)(q27.3;q32.3)")],
      ["01N23-020 (ASR)"], ["ASR"], ABB_INT + "3"))

add(P("v-ccnd1-ba", V, "Vysis LSI CCND1 Break Apart Rearrangement FISH Probe Kit", "Lymphoid", "BA",
      [comp("CCND1 (Cen)", "G", [CS(a5("CCND1"), 450, 169), CS(a5("CCND1"), 649, 496)], "~496 + ~169 kb", "RH46314 to SHGC-154860"),
       comp("CCND1 (Tel)", "O", TS(a5("CCND1"), -5, 532), "~532 kb", "CCND1, FGF19, FGF4, FGF3")],
      [NORMAL("2F"), pat("t(11;14)", [T(BP["CCND1"], BP["IGH"])], "1O1G1F", "Breakpoint within the gap between the two probes.", "vendor", "t(11;14)(q13.3;q32.3)")],
      ["05N38-020 (CE kit)", "05J96-001 (ASR)"], ["CE", "ASR"], ABB_INT + "11"))

add(P("v-igh-ccnd1", V, "Vysis LSI IGH/CCND1 DF FISH Probe Kit", "Lymphoid", "DF",
      [comp("CCND1", "O", TS(a5("CCND1"), -170, 378), "~378 kb", "D11S1076 to FGF4"), V_IGH_G2],
      [NORMAL("2O2G"),
       pat("t(11;14)", [T(BP["CCND1"], BP["IGH"])], "1O1G2F", "Breakpoints at the MTC and the IGHJ region. The IGH probe gap can make normal IGH loci look like two close green signals.", "vendor", "t(11;14)(q13.3;q32.3)")],
      ["08L58-020 (CE kit)"], ["CE"], ABB_INT + "11"))

add(P("v-igh-ccnd1-xt", V, "Vysis IGH/CCND1 XT DF FISH Probe Kit", "Myeloma", "DF",
      [comp("CCND1/MYEOV (XT)", "O", [R("11", 68838000, 69014000), R("11", 69044000, 69263000), R("11", 69325000, 69780000)], "~176 + ~219 + ~455 kb", "TPCN2, MYEOV, CCND1, FGF19, FGF4", locus="CCND1"),
       V_IGH_G],
      [NORMAL("2O2G"), pat("t(11;14)", [T(BP["CCND1"], BP["IGH_sw"])], "1O1G2F", "Extended CCND1 probe covering both MYEOV and CCND1 breakpoints (myeloma).", "vendor", "t(11;14)(q13.3;q32.3)")],
      ["05N33-020 (CE kit)", "05J72-001 (ASR)"], ["CE", "ASR"], ABB_INT + "11"))

add(P("v-igh-malt1", V, "Vysis LSI IGH/MALT1 Dual Color Dual Fusion Probes", "Lymphoid", "DF",
      [comp("MALT1", "O", TS(a5("MALT1"), -400, 670), "~670 kb", "RH15698, ALPK2, MALT1, D18S531"), V_IGH_G],
      [NORMAL("2O2G", src="standard"), pat("t(14;18)(q32;q21) IGH::MALT1", [T(BP["IGH"], BP["MALT1"])], "1O1G2F", None, "standard", "t(14;18)(q32.3;q21.3)")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None))

add(P("v-birc3-malt1", V, "Vysis LSI BIRC3/MALT1 DF FISH Probe Kit", "Lymphoid", "DF",
      [comp("BIRC3", "G", TS(a5("BIRC3"), -650, 932), "~932 kb", "D11S1247, YAP1, BIRC3, BIRC2"),
       comp("MALT1", "O", TS(a5("MALT1"), -400, 670), "~670 kb", None)],
      [NORMAL("2O2G"), pat("t(11;18)", [T(BP["BIRC3"], BP["MALT1"])], "1O1G2F", None, "vendor", "t(11;18)(q22.2;q21.3)")],
      ["05N50-020 (CE kit)"], ["CE"], ABB_INT + "11"))

add(P("v-malt1-ba", V, "Vysis LSI MALT1 Dual Color Break Apart Rearrangement Probe", "Lymphoid", "BA",
      [comp("MALT1 (Cen)", "O", CS(af("MALT1", 0.3), 0, 600), "~600 kb", "NEDD4L, ALPK2, MALT1 5'"),
       comp("MALT1 (Tel)", "G", TS(af("MALT1", 0.3), 10, 765), "~765 kb", "GRP, LMAN1")],
      [NORMAL("2F", src="standard"), pat("MALT1 rearranged", [T(BP["BIRC3"], BP["MALT1"])], "1O1G1F", None, "standard", "t(11;18)(q22.2;q21.3)")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None))

add(P("v-tra-ba", V, "Vysis LSI TRA/D Dual Color Break Apart Rearrangement Probe", "Lymphoid", "BA",
      [comp("TRA/D (Cen)", "O", CS(a5("TRA"), 0, 659), "~659 kb", "NDRG2, RPGRIP1"),
       comp("TRA/D (Tel)", "G", TS(a3("TRA"), 0, 714), "~714 kb", "RH53928 to SHGC-89225")],
      [NORMAL("2F", src="standard"), pat("inv(14)(q11.2q32)", [INV(BP["TRA"], BP["TCL1A"])], "1O1G1F", "T-PLL. The ~978 kb gap covers the TRA/TRD locus.", "standard", "inv(14)(q11.2q32.1)")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None))

# ================================================================== MYELOMA
add(P("v-igh-fgfr3", V, "Vysis IGH/FGFR3 DF FISH Probe Kit", "Myeloma", "DF",
      [comp("FGFR3/NSD2", "O", [R("4", 1280000, 2030000), R("4", 2080000, 2330000)], "~1 Mb (two parts)", "MAEA to D4S43, spans FGFR3 and NSD2", locus="FGFR3"), V_IGH_G],
      [NORMAL("2O2G"), pat("t(4;14)", [T(BP["FGFR3_NSD2"], BP["IGH_sw"])], "1O1G2F", None, "vendor", "t(4;14)(p16.3;q32.3)")],
      ["01N69-020 (CE kit)", "05J74-001 (ASR)"], ["CE", "ASR"], ABB_INT + "4"))

add(P("v-igh-maf", V, "Vysis LSI IGH/MAF Dual Color Dual Fusion Probes", "Myeloma", "DF",
      [comp("MAF", "O", [CS(a3("MAF"), 2200 - 14, 350), TS(a3("MAF"), -14, 346)], "~350 + ~346 kb", "flanks WWOX and MAF, ~2.2 Mb gap"), V_IGH_G],
      [NORMAL("2O2G", src="standard"), pat("t(14;16)", [T(BP["IGH_sw"], BP["MAF"])], "1O1G2F", None, "standard", "t(14;16)(q32.3;q23.2)")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None))

add(P("v-cks1b-cdkn2c", V, "Vysis 1q21 CKS1B SpectrumOrange/1p32 CDKN2C SpectrumGreen FISH Probe Kit", "Myeloma", "LSI",
      [comp("CKS1B (1q21)", "O", TS(a5("CKS1B"), -150, 303), "~303 kb", None),
       comp("CDKN2C (1p32)", "G", TS(a5("CDKN2C"), -60, 228), "~228 kb", None)],
      [NORMAL("2O2G", src="standard"),
       pat("1q gain", [GAINREG(G("CKS1B"), 1)], "3O2G", "Three or more CKS1B signals.", "standard"),
       pat("1q amplification + 1p deletion", [GAINREG(G("CKS1B"), 2), DEL("1", CS(a5("CDKN2C"), 500, 0), TS(a5("CDKN2C"), 500, 0))], "4O1G", None, "standard")],
      ["US ordering list 2019"], ["ASR"], None))

# ================================================================== CLL
add(P("v-cll", V, "Vysis CLL FISH Probe Kit", "CLL", "Panel",
      [comp("TP53 (17p13.1)", "O", TS(a5("TP53"), -60, 172), "~172 kb", "probe set 1", set=1),
       comp("ATM (11q22.3)", "G", TS(a5("ATM"), -250, 731), "~731 kb", "probe set 1", set=1),
       comp("D13S319 (13q14.3)", "O", TS(a5("DLEU1"), -20, 135), "~135 kb", "probe set 2", set=2),
       comp("13q34", "A", TS(a5("CUL4A"), -120, 612), "~612 kb", "probe set 2, control", set=2, locus="LAMP1"),
       comp("CEP 12", "G", CEN("12"), "alpha satellite", "probe set 2", set=2, locus="D12Z3")],
      [NORMAL("2O2G | 2O2A2G", "Set 1: 2 orange (TP53), 2 green (ATM). Set 2: 2 orange (D13S319), 2 aqua (13q34), 2 green (CEP 12)."),
       pat("del(17p)", [DEL("17", bp("17", 1), bp("17", 20000000))], None, "Set 1: one orange (TP53) signal.", "vendor"),
       pat("del(11q)", [DEL("11", bp("11", 100000000), bp("11", 115000000))], None, "Set 1: one green (ATM) signal.", "vendor"),
       pat("del(13q14), monoallelic", [DEL("13", bp("13", 49000000), bp("13", 51500000))], None, "Set 2: one orange, two aqua.", "vendor"),
       pat("del(13q14), biallelic", [DELBOTH("13", bp("13", 49800000), bp("13", 51000000))], None, "Set 2: no orange, two aqua.", "vendor"),
       pat("trisomy 12", [GAIN("12")], None, "Set 2: three green (CEP 12).", "vendor")],
      ["04N02-021 (US-IVD)"], ["US-IVD"], "https://www.molecular.abbott/content/dam/add/molecular/vysis-fish-knowledge-center/us-package-inserts/Vysis-CLL-FISH-Probe-Kit-PI_30-608712R4.pdf",
      "Two hybridisation areas on one slide. A 2021 field notice reported missed 13q deletions with some lots."))

add(P("v-atm-cep11", V, "Vysis LSI ATM/CEP 11 FISH Probe Kit", "CLL", "LSI+CEP",
      [comp("ATM", "O", TS(a5("ATM"), -250, 732), "~732 kb", "RAB39, CUL5, ATM, EXPH5"), cep("11", "G")],
      [NORMAL("2O2G"), pat("del(11q)", [DEL("11", bp("11", 100000000), bp("11", 115000000))], "1O2G", None, "vendor", "del(11)(q22.1q23.3)")],
      ["05N55-020 (CE kit)", "01N18-020 (ASR)", "01N33-020 (ATM only, CE)"], ["CE", "ASR"], ABB_INT + "11"))

add(P("v-tp53", V, "Vysis LSI TP53 (17p13.1) SpectrumOrange Probe", "CLL", "LSI",
      [comp("TP53", "O", TS(a5("TP53"), -60, 172), "~172 kb", "D17S812E, SAT2, TP53, EFNB3")],
      [NORMAL("2O", src="standard"), pat("del(17p)", [DEL("17", bp("17", 1), bp("17", 20000000))], "1O", None, "standard", "del(17)(p12)")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None, "Often paired with a chromosome 17 centromere probe (CEP 17)."))

add(P("v-d13s319", V, "Vysis LSI D13S319 (13q14.3) SpectrumOrange + LSI 13q34 Probes", "CLL", "LSI+control",
      [comp("D13S319", "O", TS(a5("DLEU1"), -20, 135), "~135 kb", "DLEU1 region"),
       comp("13q34", "G", TS(a5("CUL4A"), -120, 612), "~612 kb", "PROZ, CUL4A, LAMP1")],
      [NORMAL("2O2G", src="standard"),
       pat("del(13q14), monoallelic", [DEL("13", bp("13", 49000000), bp("13", 51500000))], "1O2G", None, "standard", "del(13)(q14.2q14.3)"),
       pat("del(13q14), biallelic", [DELBOTH("13", bp("13", 49800000), bp("13", 51000000))], "2G", None, "standard")],
      ["US ordering list 2019 (ASR); 13q34 also as SpectrumAqua"], ["ASR"], None))

# ================================================================== SOLID TUMOURS
add(P("v-alk-ba", V, "Vysis ALK Break Apart FISH Probe Kit", "Solid tumour", "BA",
      [comp("3' ALK", "O", TS(af("ALK", 0.78), 0, 300), "~300 kb", "telomeric, 3' of the breakpoint"),
       comp("5' ALK", "G", CS(af("ALK", 0.78), 0, 442), "~442 kb", "5' ALK and centromeric flank")],
      [NORMAL("2F"),
       pat("ALK rearranged (split)", [INV(BP["EML4"], BP["ALK"])], "1O1G1F", "EML4::ALK inv(2): the probes separate by more than two signal diameters.", "vendor", "inv(2)(p21p23.2)"),
       pat("ALK rearranged (isolated 3')", [INV(BP["EML4"], BP["ALK"]), DEL("2", BP["ALK"], CS(BP["ALK"], 442, 0))], "1O1F", "Isolated 3' (orange) signal, also counted as positive.", "vendor")],
      ["06N38-033 / 06N38-053 (IVD)", "CE kit"], ["IVD", "CE"], ABB_INT + "2"))

add(P("v-ros1-ba", V, "Vysis ROS1 Break Apart FISH Probe Kit", "Solid tumour", "BA",
      [comp("ROS1 (Cen)", "G", CS(af("ROS1", 0.8), 0, 557), "~557 kb", "SHGC-14567 to ROS1 3'"),
       comp("ROS1 (Tel)", "O", TS(af("ROS1", 0.8), 0, 317), "~317 kb", "ROS1 5', GOPC")],
      [NORMAL("2F", src="standard"), pat("ROS1 rearranged", [T(BP["ROS1"], BP["CD74"])], "1O1G1F", None, "standard", "t(5;6)(q32;q22.1)")],
      ["08N29-021 (CE kit)"], ["CE"], ABB_INT + "6"))

add(P("v-ret-ba", V, "Vysis 10q11 RET Break-Apart FISH Probe Kit", "Solid tumour", "BA",
      [comp("RET (Cen)", "G", CS(a5("RET"), 30, 545), "~545 kb", "SHGC-8894 side"),
       comp("RET (Tel)", "R", TS(a3("RET"), 0, 469), "~469 kb", "to RH69007")],
      [NORMAL("2F", src="standard"), pat("RET rearranged", [INV(BP["KIF5B"], BP["RET"])], "1R1G1F", None, "standard", "inv(10)(p11.2q11.2)")],
      ["RUO"], ["RUO"], ABB_INT + "10"))

add(P("v-fgfr2-ba", V, "Vysis 10q26 FGFR2 Break Apart FISH Probe Kit", "Solid tumour", "BA",
      [comp("FGFR2 (Cen)", "O", CS(a3("FGFR2"), 0, 414), "~414 kb", "D10S294 side"),
       comp("FGFR2 (Tel)", "G", TS(a3("FGFR2"), 44, 491), "~491 kb", "FGFR2 5', ATE1")],
      [NORMAL("2F", src="standard"), pat("FGFR2 rearranged", [T(BP["FGFR2"], af("BICC1", 0.8))], "1O1G1F", None, "standard")],
      ["RUO"], ["RUO"], ABB_INT + "10"))

add(P("v-ewsr1-ba", V, "Vysis LSI EWSR1 (22q12) Dual Color Break Apart Rearrangement Probe", "Solid tumour", "BA",
      [comp("5' EWSR1", "O", CS(af("EWSR1", 0.55), 0, 497), "~497 kb", "RH27882 to EWSR1"),
       comp("3' EWSR1", "G", TS(af("EWSR1", 0.55), 0, 1100), "~1.1 Mb", "EWSR1 3', NF2")],
      [NORMAL("2F", src="standard"), pat("EWSR1 rearranged, t(11;22)", [T(BP["FLI1"], BP["EWSR1"])], "1O1G1F", None, "standard", "t(11;22)(q24.3;q12.2)")],
      ["US ordering list 2019"], ["ASR"], None))

add(P("v-fus-ba", V, "Vysis LSI FUS Dual Color Break Apart Rearrangement Probe", "Solid tumour", "BA",
      [comp("FUS (Cen)", "O", CS(af("FUS", 0.4), 0, 275), "~275 kb", "FUS 3', ITGAX to RH65072"),
       comp("FUS (Tel)", "G", TS(af("FUS", 0.4), 30, 505), "~505 kb", "SHGC-84622, RNF40, BCL7C")],
      [NORMAL("2F", src="standard"), pat("FUS rearranged, t(12;16)", [T(BP["DDIT3"], BP["FUS"])], "1O1G1F", None, "standard", "t(12;16)(q13.3;p11.2)")],
      ["Components: FUS (Cen) SpectrumOrange, FUS (Tel) SpectrumGreen (ASR)"], ["ASR"], None))

add(P("v-ddit3-ba", V, "Vysis DDIT3 Break Apart FISH Probe Kit", "Solid tumour", "BA",
      [comp("DDIT3 (Cen)", "O", CS(a3("DDIT3"), 0, 700), "~700 kb", "HSD17B6, STAT6, GLI1"),
       comp("DDIT3 (Tel)", "G", TS(a3("DDIT3"), 134, 663), "~663 kb", "OS9, CDK4, XRCC6BP1")],
      [NORMAL("2F"), pat("DDIT3 rearranged, t(12;16)", [T(BP["DDIT3"], BP["FUS"])], "1O1G1F", None, "vendor", "t(12;16)(q13.3;p11.2)")],
      ["03N57-020 (CE kit)"], ["CE"], ABB_INT + "12"))

add(P("v-ss18-ba", V, "Vysis LSI SS18 Dual Color Break Apart Rearrangement Probe", "Solid tumour", "BA",
      [comp("SS18 (Cen)", "G", CS(af("SS18", 0.9), 0, 1044), "~1044 kb", "ZNF521 to SS18 3'"),
       comp("SS18 (Tel)", "O", TS(af("SS18", 0.9), 0, 650), "~650 kb", "SS18 5', TAF4B")],
      [NORMAL("2F", src="standard"), pat("SS18 rearranged, t(X;18)", [T(BP["SSX1"], BP["SS18"])], "1O1G1F", None, "standard", "t(X;18)(p11.2;q11.2)")],
      ["Components: SS18 (Cen) SpectrumGreen, SS18 (Tel) SpectrumOrange (ASR)"], ["ASR"], None))

add(P("v-foxo1-ba", V, "Vysis LSI FOXO1 Dual Color Break Apart Rearrangement Probe", "Solid tumour", "BA",
      [comp("FOXO1 (Cen)", "G", CS(a3("FOXO1"), 0, 724), "~724 kb", "SHGC-4344 to the 3' end of FOXO1"),
       comp("FOXO1 (Tel)", "O", TS(a5("FOXO1"), 0, 655), "~655 kb", "from the 5' end of FOXO1 past ELF1; the ~109 kb gap covers FOXO1")],
      [NORMAL("2F", src="standard"), pat("FOXO1 rearranged, t(2;13)", [T(BP["PAX3"], BP["FOXO1"])], "1O1G1F", None, "standard", "t(2;13)(q36.1;q14.1)")],
      ["Components: FOXO1 (Cen) SpectrumGreen, FOXO1 (Tel) SpectrumOrange (ASR)"], ["ASR"], None))

add(P("v-ntrk1-ba", V, "Vysis LSI NTRK1 (Cen) SpectrumGreen + NTRK1 (Tel) SpectrumRed Probes", "Solid tumour", "BA",
      [comp("NTRK1 (Cen)", "G", CS(af("NTRK1", 0.55), 0, 440), "~440 kb", "SHGC-103490 to NTRK1"),
       comp("NTRK1 (Tel)", "R", TS(af("NTRK1", 0.55), 0, 332), "~332 kb", "to SHGC-145655")],
      [NORMAL("2F", src="standard"), pat("NTRK1 rearranged", [INV(BP["TPM3"], BP["NTRK1"])], "1R1G1F", None, "standard", "inv(1)(q21.3q23.1)")],
      ["ASR components"], ["ASR"], None))

add(P("v-erg-ba", V, "Vysis LSI ERG (Cen) SpectrumRed + ERG (Tel) SpectrumGreen Probes", "Solid tumour", "BA",
      [comp("ERG (Cen)", "R", CS(a3("ERG"), 0, 523), "~523 kb", "AFM234YD2 side, 3' of ERG"),
       comp("ERG (Tel)", "G", TS(a5("ERG"), 0, 719), "~719 kb", "5' of ERG to STSG622746")],
      [NORMAL("2F", src="standard"), pat("TMPRSS2::ERG by deletion", [DEL("21", BP["ERG"], BP["TMPRSS2"])], "1R1F", "Interstitial deletion removes the green (5') signal.", "standard")],
      ["ASR components"], ["ASR"], None))

add(P("v-pathvysion", V, "PathVysion HER-2 DNA Probe Kit", "Solid tumour", "LSI+CEP",
      [comp("HER-2/neu (ERBB2)", "O", TS(a5("ERBB2"), -100, 226), "~226 kb", None), cep("17", "A")],
      [NORMAL("2O2A", src="standard"), pat("HER2 amplified", [AMP(G("ERBB2"), 14)], "2O2A amp O", "ERBB2/CEP17 ratio ≥ 2.0 (ASCO/CAP algorithm).", "standard")],
      ["FDA-approved IVD"], ["IVD"], None))

add(P("v-egfr-cep7", V, "Vysis EGFR/CEP 7 FISH Probe Kit", "Solid tumour", "LSI+CEP",
      [comp("EGFR", "O", TS(a5("EGFR"), -30, 303), "~303 kb", "RH70813 to SHGC-101379"), cep("7", "G")],
      [NORMAL("2O2G"), pat("EGFR amplified", [AMP(G("EGFR"), 14)], "2O2G amp O", None, "vendor"), pat("Polysomy 7", [GAIN("7", 2)], "4O4G", None, "standard")],
      ["01N35-020 (CE kit)"], ["CE"], ABB_INT + "7"))

add(P("v-met", V, "Vysis MET SpectrumRed FISH Probe Kit", "Solid tumour", "LSI",
      [comp("MET", "R", TS(a5("MET"), -150, 456), "~456 kb", None)],
      [NORMAL("2R"), pat("MET amplified", [AMP(G("MET"), 12)], "2R amp R", None, "vendor")],
      ["06N05-020 (CE kit)"], ["CE"], ABB_INT + "7"))

add(P("v-mdm2-cep12", V, "Vysis MDM2/CEP 12 FISH Probe Kit", "Solid tumour", "LSI+CEP",
      [comp("MDM2", "O", TS(a5("MDM2"), -40, 209), "~209 kb", "NUP107 to CPM"), cep("12", "G")],
      [NORMAL("2O2G"), pat("MDM2 amplified", [AMP(G("MDM2"), 16)], "2O2G amp O", "Typical of well-differentiated/dedifferentiated liposarcoma.", "vendor")],
      ["01N15-010 (CE kit)"], ["CE"], ABB_INT + "12"))

add(P("v-nmyc-cep2", V, "Vysis LSI N-MYC SpectrumGreen/CEP 2 SpectrumOrange Probes", "Solid tumour", "LSI+CEP",
      [comp("MYCN", "G", TS(a5("MYCN"), -90, 200), "~200 kb", None), cep("2", "O")],
      [NORMAL("2G2O"), pat("MYCN amplified", [AMP(G("MYCN"), 20)], "2G2O amp G", None, "vendor")],
      ["07J72-001 (ASR)"], ["ASR"], ABB_INT + "2"))

add(P("v-pten-cep10", V, "Vysis LSI PTEN/CEP 10 FISH Probe Kit", "Solid tumour", "LSI+CEP",
      [comp("PTEN", "O", TS(a5("PTEN"), -100, 368), "~368 kb", "D10S215 to RH93626"), cep("10", "G")],
      [NORMAL("2O2G"), pat("PTEN deletion", [DEL("10", bp("10", 86000000), bp("10", 90500000))], "1O2G", None, "vendor", "del(10)(q23.2q23.31)")],
      ["04N62-020 (CE kit)"], ["CE"], ABB_INT + "10"))

add(P("v-cdkn2a-cep9", V, "Vysis LSI CDKN2A/CEP 9 FISH Probe Kit", "Solid tumour", "LSI+CEP",
      [comp("CDKN2A (p16)", "O", CS(a5("MTAP"), -20, 222), "~222 kb", "MTAP, CDKN2A, CDKN2B"), cep("9", "G")],
      [NORMAL("2O2G"),
       pat("CDKN2A deletion", [DEL("9", bp("9", 21000000), bp("9", 23000000))], "1O2G", None, "vendor", "del(9)(p21.3p21.3)"),
       pat("Biallelic CDKN2A deletion", [DELBOTH("9", bp("9", 21300000), bp("9", 22500000))], "2G", "Common in ALL and many solid tumours.", "standard")],
      ["04N61-020 (CE kit)"], ["CE"], ABB_INT + "9"))

add(P("v-1p19q", V, "Vysis LSI 1p36/1q25 and 19q13/19p13 Probes", "Solid tumour", "LSI+control",
      [comp("1p36", "O", TS(a5("MEGF6"), -100, 435), "~435 kb", "MEGF6, TP73", set=1, locus="TP73"),
       comp("1q25", "G", TS(a5("ANGPTL1"), -100, 618), "~618 kb", "ANGPTL1, ABL2", set=1, locus="ABL2"),
       comp("19q13", "O", TS(a5("BICRA"), -50, 380), "~380 kb", "GLTSCR1 (BICRA), GLTSCR2, CRX region", set=2, locus="BICRA"),
       comp("19p13", "G", TS(a5("ZNF443"), -250, 502), "~502 kb", "ZNF44, ZNF443, MAN2B1", set=2, locus="ZNF443")],
      [NORMAL("2O2G | 2O2G", "Two hybridisation areas: 1p36/1q25 and 19q13/19p13."), pat("1p/19q codeletion", [ARMLOSS("1", "p"), ARMLOSS("19", "q")], "1O2G | 1O2G", "Oligodendroglioma: whole-arm loss of 1p and 19q, usually through an unbalanced t(1;19)(q10;p10).", "standard")],
      ["US ordering list 2019 (ASR)"], ["ASR"], None))

add(P("v-urovysion", V, "UroVysion Bladder Cancer Kit", "Solid tumour", "CEP",
      [cep("3", "R"), cep("7", "G"), cep("17", "A"), comp("9p21", "Y", CS(a5("MTAP"), -20, 222), "~222 kb", "CDKN2A", locus="CDKN2A")],
      [NORMAL("2R2G2A2Y", src="standard"),
       pat("Polysomy (aneuploid urothelial cell)", [GAIN("3", 2), GAIN("7", 1), GAIN("17", 2)], "4R3G4A2Y", "Gains of two or more of chromosomes 3, 7 and 17 in 4 or more of 25 cells.", "standard"),
       pat("9p21 loss", [DELBOTH("9", bp("9", 21300000), bp("9", 22500000))], "2R2G2A", "Homozygous 9p21 loss in 12 or more of 25 cells.", "standard")],
      ["US-IVD", "CE"], ["IVD", "CE"], ABB_INT + "9"))

add(P("v-melanoma", V, "Vysis Melanoma FISH Probe Kit", "Solid tumour", "LSI+CEP",
      [comp("RREB1 (6p25)", "R", CS(a5("RREB1"), -100, 638), "~638 kb", None),
       comp("MYB (6q23)", "Y", TS(a5("MYB"), -500, 740), "~740 kb", None),
       comp("CCND1 (11q13)", "G", TS(a5("CCND1"), -170, 378), "~378 kb", None),
       cep("6", "A")],
      [NORMAL("2R2Y2G2A"), pat("RREB1 gain", [GAINREG(G("RREB1"), 2)], "4R2Y2G2A", None, "vendor")],
      ["01N89-020 (CE kit)"], ["CE"], ABB_INT + "6"))

# ================================================================== CONSTITUTIONAL
add(P("v-aneuvysion", V, "AneuVysion Multicolor DNA Probe Kit", "Constitutional", "CEP",
      [dict(cep("18", "A"), set=1), dict(cep("X", "G"), set=1), comp("CEP Y", "O", CEN("Y"), "alpha satellite", "DYZ3", set=1, locus="DYZ3"),
       comp("LSI 13 (13q14)", "G", TS(a5("RB1"), -100, 422), "~422 kb", "vial 2", set=2, locus="RB1"),
       comp("LSI 21 (21q22.13-q22.2)", "O", TS(a5("DYRK1A"), -50, 200), "~200 kb", "vial 2", set=2, locus="D21S259")],
      [NORMAL("2A2G | 2G2O", "Female. Vial 1 (18/X/Y): 2 aqua, 2 green. Vial 2 (13/21): 2 green, 2 orange."),
       pat("Trisomy 21", [GAIN("21")], None, "Vial 2: three orange (21) signals.", "standard", "47,XX,+21"),
       pat("Trisomy 18, male", [GAIN("18"), LOSS("X"), GAIN("Y")], None, "Vial 1: three aqua, one green, one orange.", "vendor", "47,XY,+18")],
      ["05J38-010/-030/-050 (IVD)"], ["IVD", "CE"], ABB_INT + "x"))

add(P("v-multivysion", V, "Vysis MultiVysion PGT Multi-color Probe", "Constitutional", "CEP",
      [comp("LSI 13", "R", TS(a5("RB1"), -100, 422), "~422 kb", None, locus="RB1"), cep("18", "A"),
       comp("LSI 21", "G", TS(a5("DYRK1A"), -200, 350), "~350 kb", None, locus="D21S259"), cep("X", "B", "CEP X (SpectrumBlue)"), comp("CEP Y", "Y", CEN("Y"), "alpha satellite", None)],
      [NORMAL("2R2A2G2B", src="standard")],
      ["08L69-010 (CE)"], ["CE"], ABB_INT + "x"))

def cons(id, name, lsi_label, lsi_region, size, ctrl, ctrl_label, ctrl_region, ctrl_size, deleted, codes, url=None, ctrl_color="G", lsi_color="O", note=None, del_ev=None):
    return P(id, V, name, "Constitutional", "LSI+control",
             [comp(lsi_label, lsi_color, lsi_region, size), comp(ctrl_label, ctrl_color, ctrl_region, ctrl_size, "control")],
             [NORMAL(f"2{lsi_color}2{ctrl_color}", src="standard"),
              pat(deleted, [del_ev] if del_ev else [], f"1{lsi_color}2{ctrl_color}", "Seen best on metaphases: one chromosome lacks the orange signal.", "vendor")],
             codes, ["ASR"], url, note)

add(cons("v-digeorge-tuple1", "Vysis DiGeorge Region LSI TUPLE 1 SpectrumOrange/LSI ARSA SpectrumGreen Probes", "TUPLE1 (HIRA)", TS(a3("HIRA"), -20, 117), "~117 kb", "G", "ARSA (22q13.3)", TS(a5("SBF1"), -20, 334), "~334 kb",
         "22q11.2 deletion", ["05J21-028 (ASR)"], del_ev=DEL("22", bp("22", 18900000), bp("22", 21500000))))
add(cons("v-digeorge-n25", "Vysis DiGeorge Region LSI N25 SpectrumOrange/LSI ARSA SpectrumGreen Probes", "N25 (D22S75)", TS(a3("CLTCL1"), -50, 126), "~126 kb", "G", "ARSA (22q13.3)", TS(a5("SBF1"), -20, 334), "~334 kb",
         "22q11.2 deletion", ["05N24-010 (ASR)"], del_ev=DEL("22", bp("22", 18900000), bp("22", 21500000))))
add(cons("v-williams", "Vysis Williams Region Probe - LSI ELN SpectrumOrange/LSI D7S486, D7S522 SpectrumGreen Probes", "ELN (7q11.23)", TS(a5("ELN"), -20, 215), "~215 kb", "G", "D7S486, D7S522 (7q31)", TS(a3("TES"), -100, 468), "~468 kb",
         "Williams-Beuren deletion", ["05J30-045 (ASR)", "CE kit"], ABB_INT + "7", del_ev=DEL("7", bp("7", 73300000), bp("7", 75000000))))
add(cons("v-pws-snrpn", "Vysis Prader-Willi/Angelman Region LSI SNRPN SpectrumOrange/CEP 15 SpectrumAqua/PML SpectrumGreen Probes", "SNRPN (15q11.2)", TS(a5("SNRPN"), -20, 125), "~125 kb", "G", "PML (15q24)", CS(a3("PML"), 0, 239), "~239 kb",
         "15q11.2-q13 deletion", ["01N12-010 (ASR)"], del_ev=DEL("15", bp("15", 22800000), bp("15", 28500000)), note="CEP 15 SpectrumAqua is also in the mix."))
add(cons("v-miller-dieker", "Vysis Miller-Dieker Region/Isolated Lissencephaly LSI LIS1 SpectrumOrange/RARA SpectrumGreen Probes", "LIS1 (PAFAH1B1)", TS(a5("PAFAH1B1"), -40, 110), "~110 kb", "G", "RARA (17q21)", TS(a3("RARA"), 0, 417), "~417 kb",
         "17p13.3 deletion", ["05J88-001 (ASR)"], del_ev=DEL("17", bp("17", 1), bp("17", 3000000))))
add(cons("v-smith-magenis", "Vysis Smith-Magenis Region LSI SMS SpectrumOrange/LSI RARA SpectrumGreen Probes", "SMS (17p11.2)", TS(a5("SHMT1"), -20, 146), "~146 kb", "G", "RARA (17q21)", TS(a3("RARA"), 0, 417), "~417 kb",
         "17p11.2 deletion", ["US ordering list 2019 (ASR)"], del_ev=DEL("17", bp("17", 16700000), bp("17", 20500000))))
add(cons("v-sotos", "Vysis Sotos Region LSI NSD1 SpectrumOrange/LSI D5S23, D5S721 SpectrumGreen Probes", "NSD1 (5q35)", TS(a5("NSD1"), -10, 180), "~180 kb", "G", "D5S23, D5S721 (5p15.2)", R("5", 9900000, 10461000), "~561 kb",
         "NSD1 deletion", ["05J48-007 (ASR)"], ABB_INT + "5", del_ev=DEL("5", bp("5", 176800000), bp("5", 178000000))))
add(cons("v-whs", "Vysis Wolf-Hirschhorn Region LSI WHS SpectrumOrange/CEP 4 SpectrumGreen Probes", "WHS (4p16.3)", TS(af("NSD2", 0.7), 0, 98), "~98 kb", "G", "CEP 4", CEN("4"), "alpha satellite",
         "4p16.3 deletion", ["05J29-074 (ASR)"], ABB_INT + "4", del_ev=DELTER(bp("4", 3000000), "tel")))
add(cons("v-kallmann", "Vysis Kallmann Region LSI KAL SpectrumOrange/CEP X SpectrumGreen Probes", "KAL (ANOS1)", R("X", 8540000, 8700000), "~160 kb", "G", "CEP X", CEN("X"), "alpha satellite",
         "Xp22.3 deletion", ["05J23-070 (ASR)"], ABB_INT + "x", del_ev=DEL("X", bp("X", 6000000), bp("X", 9000000))))
add(cons("v-sts", "Vysis Steroid Sulfatase Deficiency LSI STS SpectrumOrange/CEP X SpectrumGreen Probes", "STS", CS(a5("STS"), -10, 244), "~244 kb", "G", "CEP X", CEN("X"), "alpha satellite",
         "STS deletion", ["05J28-004 (ASR)"], ABB_INT + "x", del_ev=DEL("X", bp("X", 6400000), bp("X", 8200000))))
add(cons("v-sry", "Vysis LSI SRY SpectrumOrange/CEP X SpectrumGreen Probes", "SRY (Yp11.3)", CS(a5("SRY"), 0, 122), "~122 kb", "G", "CEP X", CEN("X"), "alpha satellite",
         "SRY-positive female (46,XX testicular DSD)", ["05J27-007 (ASR)", "06N29-020 (CE kit)"], ABB_INT + "y"))
PROBES[-1]["patterns"] = [
    pat("Normal female", [], "2G", "Two CEP X, no SRY.", "standard"),
    pat("Normal male", [LOSS("X"), GAIN("Y")], "1O1G", None, "standard"),
    pat("46,XX with SRY on an X", [{"ins": ["X", bp("X", 2700000), "Y", bp("Y", 2780000), bp("Y", 2900000)]}], "1O2G", "SRY translocated to Xp: testicular DSD.", "standard"),
]
add(P("v-1p36-microdel", V, "Vysis LSI 1p36 Microdeletion Region Probe", "Constitutional", "LSI",
      [comp("1pTel", "G", TS(a5("AGRN"), -30, 101), "~101 kb", "AGRN"),
       comp("p58 (CDK11)", "O", TS(a3("CDK11B"), -60, 111), "~111 kb", "CDK11B/CDK11A, ~577 kb from the green probe")],
      [NORMAL("2F", "The two probes are ~577 kb apart and read as adjacent signals.", "standard"), pat("1p36 deletion", [DELTER(bp("1", 3000000), "tel")], "1F", "Terminal 1p36 deletion removes both signals from one chromosome.", "standard")],
      ["05J21-020 (ASR)"], ["ASR"], None))

# ---- single CEP probes (enumeration)
# CEP 8 and CEP 12 are listed above as kits (v-cep8, v-cep12).
CEP_LIST = [("1", "O"), ("2", "O"), ("3", "O"), ("4", "OGA"), ("6", "OGA"), ("7", "OGA"), ("9", "OGA"), ("10", "OGA"),
            ("11", "OGA"), ("15", "OGA"), ("16", "OGA"), ("17", "OGA"), ("18", "OGA"), ("20", "O"), ("X", "OGA"), ("Y", "OGA")]
for c, cols in CEP_LIST:
    add(P(f"v-cep{c.lower()}-single", V, f"Vysis CEP {c} Probe", "Enumeration (CEP)", "CEP",
          [cep(c, cols[0])],
          [NORMAL(f"2{cols[0]}" if c not in ("Y",) else f"1{cols[0]}", src="standard"),
           pat(f"Gain of {c}", [GAIN(c)], None, None, "standard"),
           pat(f"Loss of {c}", [LOSS(c)], None, None, "standard")],
          [f"Available as Spectrum{', Spectrum'.join({'O':'Orange','G':'Green','A':'Aqua'}[x] for x in cols)}"], ["ASR"], ABB_INT + (c.lower() if c in "XY" else c)))
    if c == "Y":
        PROBES[-1]["patterns"] = [pat("Normal male", [LOSS("X"), GAIN("Y")], "1" + cols[0], None, "standard"),
                                  pat("Female (no Y)", [], "", None, "standard"),
                                  pat("47,XYY", [LOSS("X"), GAIN("Y", 2)], "2" + cols[0], None, "standard")]
    if c == "X":
        PROBES[-1]["patterns"] = [pat("Normal female", [], "2" + cols[0], None, "standard"),
                                  pat("Normal male / 45,X", [LOSS("X")], "1" + cols[0], None, "standard"),
                                  pat("47,XXX", [GAIN("X")], "3" + cols[0], None, "standard")]

# ---- TelVysion subtelomeric probes
for c in [str(i) for i in range(1, 23)] + ["X"]:
    arms = ["q"] if c in ("13", "14", "15", "21", "22") else ["p", "q"]
    comps = [comp(f"{c}{a}ter", "G" if a == "p" else "O", TEL(c, a), "~100-300 kb from the end", None) for a in arms]
    label = "Xp/Yp and Xq/Yq (pseudoautosomal)" if c == "X" else f"{c}{'p and ' + c + 'q' if len(arms) == 2 else 'q'}"
    add(P(f"v-telvysion-{c.lower()}", V, f"TelVysion {label}", "Subtelomere", "TEL", comps,
          [NORMAL("".join(f"2{x['color']}" for x in comps), src="standard"),
           pat("Terminal deletion", [DELTER(TEL(c, arms[-1]), "tel")], None, None, "standard")],
          ["ASR; also combined in ToTelVysion (15 multicolour mixes)"], ["ASR"], ABB_INT + (c.lower() if c in "XY" else c)))

# ================================================================== OTHER VENDORS (targets Vysis does not cover)
def cyto(slug): return "https://www.ogt.com/products/product-search/" + slug + "/"
def mets(code): return "https://metasystems-probes.com/en/probes/xl/" + code + "/"

def kmt2a_df(id, partner, partner_region, partner_size, partner_note, bpk, kt, status, slug):
    return P(id, C, f"CytoCell MLL (KMT2A)/{partner} Translocation, Dual Fusion", "Myeloid", "DF",
             [comp("KMT2A", "G", TS(a5("KMT2A"), -60, 200), "~200 kb", "covers KMT2A"),
              comp(partner, "R", partner_region, partner_size, partner_note)],
             [NORMAL("2R2G", src="standard"), pat(kt.split(")(")[0] + ")", [T(BP[bpk], BP["KMT2A"])], "1R1G2F", None, "standard", kt)],
             [], [status], cyto(slug))

add(kmt2a_df("c-kmt2a-aff1", "AFF1", [CS(a5("AFF1"), 30, 220), TS(a5("AFF1"), -30, 271), TS(a3("AFF1"), 20, 211)], "220 + 61 + 110 + 211 kb", "four clones covering AFF1 and flanks",
             "AFF1", "t(4;11)(q21.3;q23.3)", "IVDD", "cytocell-mll-kmt2a-aff1-translocation-dual-fusion-fish-probe"))
add(kmt2a_df("c-kmt2a-mllt3", "MLLT3", [R("9", 20400000, 20570000), R("9", 20570000, 20710000), R("9", 20710000, 20896000)], "140 + 167 + 186 kb", "5' MLLT3 and flanks",
             "MLLT3", "t(9;11)(p21.3;q23.3)", "RUO", "cytocell-mll-kmt2amllt3-translocation-dual-fusion-fish-probe"))
add(kmt2a_df("c-kmt2a-mllt1", "MLLT1", [R("19", 5984000, 6210000), R("19", 6280000, 6483000)], "203 + 226 kb", "two clones flanking MLLT1",
             "MLLT1", "t(11;19)(q23.3;p13.3)", "RUO", "cytocell-mll-kmt2amllt1-translocation-dual-fusion-fish-probe"))
add(kmt2a_df("c-kmt2a-afdn", "MLLT4", [R("6", 167640000, 167827000), R("6", 167980000, 168116000)], "187 + 136 kb", "two clones flanking MLLT4 (AFDN)",
             "AFDN", "t(6;11)(q27;q23.3)", "RUO", "cytocell-mll-kmt2amllt4-translocation-dual-fusion-fish-probe"))

add(P("c-cbfb-myh11", C, "CytoCell CBFβ (CBFB)/MYH11 Translocation, Dual Fusion", "Myeloid", "DF",
      [comp("CBFB", "R", TS(a5("CBFB"), -270, 617), "~617 kb", "includes CBFB"),
       comp("MYH11", "G", TS(a3("MYH11"), -230, 621), "~621 kb", "includes MYH11")],
      [NORMAL("2R2G", src="standard"), pat("inv(16)", [INV(BP["MYH11"], BP["CBFB"])], "1R1G2F", None, "standard", "inv(16)(p13.1q22.1)")],
      ["CE-LPH 022"], ["IVDR"], cyto("cytocell-cbf-cbfbmyh11-translocation-dual-fusion-fish-probe")))

add(P("c-mecom-ba", C, "CytoCell EVI1 (MECOM) Breakapart", "Myeloid", "TC-BA",
      [comp("Telomeric of MECOM", "R", TS(a5("MECOM"), 140, 158), "~158 kb", "telomeric to D3S4415, LRRC34"),
       comp("MECOM (centromeric part)", "G", TS(a3("MECOM"), -130, 181), "~181 kb", "centromeric (3') part of MECOM, beyond D3S1282"),
       comp("Centromeric of MECOM", "B", CS(a3("MECOM"), 130, 563), "~563 kb", "D3S1614")],
      [NORMAL("2RGB", "Two red/green/blue fusions.", "standard"),
       pat("inv(3)", [INV(BP["GATA2"], BP["MECOM_inv"])], "1RGB1R1GB", "Break 5' (telomeric) of MECOM separates red from green/blue.", "standard", "inv(3)(q21.3q26.2)"),
       pat("t(3;3)", [T(BP["GATA2"], CS(a3("MECOM"), 100, 0))], "1RGB1RG1B", "Break 3' (centromeric) of MECOM separates blue from red/green.", "standard", "t(3;3)(q21.3;q26.2)")],
      ["CE-LPH 036"], ["IVDR"], cyto("cytocell-evi1-mecom-breakapart-fish-probe")))

add(P("c-runx1-ba", C, "CytoCell AML1 (RUNX1) Breakapart", "Myeloid", "BA",
      [comp("Centromeric of RUNX1", "R", CS(a3("RUNX1"), 30, 156), "~156 kb", "spans CLIC6"),
       comp("RUNX1", "G", TS(a3("RUNX1"), 120, 169), "~169 kb", "part of RUNX1 (SHGC-87606, D21S1921)")],
      [NORMAL("2F", src="standard"), pat("RUNX1 rearranged", [T(BP["RUNX1T1"], af("RUNX1", 0.7))], "1R1G1F", None, "standard", "t(8;21)(q21.3;q22.1)")],
      ["CE-LPH 027"], ["IVDR"], cyto("cytocell-aml1-runx1-breakapart-fish-probe")))

add(P("c-tcf3-ba", C, "CytoCell E2A (TCF3) Breakapart", "Lymphoid", "BA",
      [comp("Centromeric of TCF3", "R", CS(a5("TCF3"), 0, 191), "~191 kb", "RH98588"),
       comp("Telomeric of TCF3", "G", TS(a3("TCF3"), 0, 165), "~165 kb", "D19S883")],
      [NORMAL("2F", src="standard"), pat("TCF3 rearranged", [T(BP["PBX1"], BP["TCF3"])], "1R1G1F", None, "standard", "t(1;19)(q23.3;p13.3)")],
      ["LPH 019"], ["IVDD"], cyto("cytocell-e2a-tcf3-breakapart-fish-probe")))

add(P("c-crlf2-ba", C, "CytoCell CRLF2 Breakapart", "Lymphoid", "BA",
      [comp("Centromeric of CRLF2", "R", CS(a5("CRLF2"), 0, 243), "~243 kb", "PAR1 of X and Y"),
       comp("Telomeric of CRLF2", "G", TS(a3("CRLF2"), 0, 202), "71 + 131 kb", None)],
      [NORMAL("2F", "Female: two fusions (X and X). Male: X and Y.", "standard"),
       pat("P2RY8::CRLF2 (PAR1 deletion)", [DEL("X", af("P2RY8", 0.3), af("CRLF2", 0.5))], "1G1F", "Interstitial PAR1 deletion removes the red probe.", "standard"),
       pat("IGH::CRLF2", [T(BP["CRLF2"], BP["IGH"])], "1R1G1F", None, "standard", "t(X;14)(p22.3;q32.3)")],
      ["RU-LPH 085"], ["RUO"], cyto("cytocell-crlf2-breakapart-fish-probe")))

add(P("m-jak2-ba", M, "XL JAK2 BA Break Apart Probe", "Lymphoid", "BA",
      [comp("Proximal JAK2", "O", CS(BP["JAK2"], 0, 400), "approx.", "proximal (centromeric)"),
       comp("Distal JAK2", "G", TS(BP["JAK2"], 0, 400), "approx.", "distal (telomeric)")],
      [NORMAL("2F", src="standard"), pat("JAK2 rearranged, t(8;9)", [T(BP["PCM1"], BP["JAK2"])], "1O1G1F", None, "standard", "t(8;9)(p22;p24.1)")],
      ["D-5098-100-OG"], ["IVDR"], mets("d-5098-100-og"), "Probe sizes are not given on the product page; positions here are illustrative."))

add(P("m-nup98-ba", M, "XL NUP98 Break Apart Probe", "Myeloid", "BA",
      [comp("Proximal NUP98", "O", CS(BP["NUP98"], 0, 400), "approx.", "proximal (centromeric)"),
       comp("Distal NUP98", "G", TS(BP["NUP98"], 0, 400), "approx.", "distal (telomeric)")],
      [NORMAL("2F", src="standard"), pat("NUP98 rearranged, t(5;11)", [T(BP["NSD1"], BP["NUP98"])], "1O1G1F", "NUP98::NSD1 is cryptic on karyotype.", "standard", "t(5;11)(q35.3;p15.4)")],
      ["D-5077-100-OG"], ["IVDR"], mets("d-5077-100-og"), "Probe sizes are not given on the product page; positions here are illustrative."))

add(P("m-pax5-ba", M, "XL PAX5 BA Break Apart Probe", "Lymphoid", "BA",
      [comp("PAX5 and distal", "O", TS(af("PAX5", 0.15), 0, 600), "approx.", "PAX5 and distal (telomeric) region"),
       comp("Proximal PAX5", "G", CS(a5("PAX5"), 20, 400), "approx.", "proximal (centromeric)")],
      [NORMAL("2F", src="standard"), pat("PAX5 rearranged", [T(BP["PAX5"], BP["ETV6"])], "1O1G1F", None, "standard")],
      ["D-5143-100-OG"], ["IVDD"], mets("d-5143-100-og"), "Probe sizes are not given on the product page; positions here are illustrative."))

add(P("m-abl2-ba", M, "XL ABL2 BA Break Apart Probe", "Lymphoid", "BA",
      [comp("Proximal ABL2", "O", CS(BP["ABL2"], 0, 400), "approx.", "proximal (centromeric) of the breakpoint"),
       comp("Distal ABL2", "G", TS(BP["ABL2"], 0, 400), "approx.", "distal (telomeric)")],
      [NORMAL("2F", src="standard"), pat("ABL2 rearranged", [T(BP["ABL2"], BP["ETV6"])], "1O1G1F", "Ph-like ALL (ABL-class fusion).", "standard")],
      ["D-5138-100-OG"], ["IVDD"], mets("d-5138-100-og"), "Probe sizes are not given on the product page; positions here are illustrative."))

add(P("m-abl1-ba", M, "XL ABL1 BA Break Apart Probe", "Lymphoid", "BA",
      [comp("Proximal ABL1", "G", CS(BP["ABL1"], 0, 400), "approx.", "proximal (centromeric)"),
       comp("Distal ABL1", "O", TS(BP["ABL1"], 0, 400), "approx.", "distal (telomeric)")],
      [NORMAL("2F", src="standard"), pat("ABL1 rearranged (e.g. NUP214::ABL1, BCR::ABL1)", [T(BP["ABL1"], BP["BCR_M"])], "1O1G1F", "Detects ABL1 rearrangement regardless of partner.", "standard")],
      ["D-5148-100-OG"], ["IVDD"], mets("d-5148-100-og"), "Probe sizes are not given on the product page; positions here are illustrative."))

add(P("m-p2ry8", M, "XL P2RY8 del Deletion Probe", "Lymphoid", "Deletion",
      [comp("P2RY8 and distal", "O", TS(a5("P2RY8"), 0, 250), "approx.", "P2RY8 and the distal (telomeric) region, PAR1 of X and Y"),
       comp("Proximal P2RY8", "G", CS(a5("P2RY8"), 0, 250), "approx.", "proximal (centromeric)")],
      [NORMAL("2F", "Two green-orange fusions (2GO).", "vendor"),
       pat("P2RY8::CRLF2 deletion", [DEL("X", af("P2RY8", 0.3), af("CRLF2", 0.5))], "1G1F", "One fusion and one green: loss of the orange signal.", "vendor")],
      ["D-5150-100-OG"], ["IVDR"], mets("d-5150-100-og"), "Probe sizes are not given on the product page; positions here are illustrative."))

add(P("c-igk-ba", C, "CytoCell IGK Breakapart", "Lymphoid", "BA",
      [comp("Distal IGK variable", "R", CS(("e", "IGKC"), 900, 301), "~301 kb", "distal variable region"),
       comp("Telomeric of IGKJ/IGKC", "G", TS(("s", "IGKC"), 0, 606), "~606 kb", "D2S2216 to D2S2510")],
      [NORMAL("2F", src="standard"), pat("IGK rearranged, t(2;8)", [T(BP["IGK"], TS(a3("MYC"), 400, 0))], "1R1G1F", None, "standard", "t(2;8)(p11.2;q24.2)")],
      ["LPH 034"], ["IVDD"], cyto("cytocell-igk-breakapart-fish-probe")))

add(P("c-igl-ba", C, "CytoCell IGL Breakapart", "Lymphoid", "BA",
      [comp("Centromeric of IGLV", "R", CS(("s", "IGL"), 0, 278), "~278 kb", "covers MAPK1"),
       comp("Telomeric of IGLC", "G", TS(("e", "IGL"), 30, 307), "~307 kb", "includes the 5' end of BCR")],
      [NORMAL("2F", src="standard"), pat("IGL rearranged, t(8;22)", [T(TS(a3("MYC"), 400, 0), BP["IGL"])], "1R1G1F", None, "standard", "t(8;22)(q24.2;q11.2)")],
      ["LPH 033"], ["IVDD"], cyto("cytocell-igl-breakapart-fish-probe")))

add(P("c-tcl1-ba", C, "CytoCell TCL1 Breakapart", "Lymphoid", "BA",
      [comp("Centromeric of TCL1A/B", "R", CS(a3("TCL1A"), 10, 138), "~138 kb", "GLRX5, D14S62"),
       comp("Telomeric of TCL1A/B", "G", TS(a5("TCL1A"), 10, 124), "~124 kb", "D14S131")],
      [NORMAL("2F", src="standard"), pat("inv(14)(q11.2q32)", [INV(BP["TRA"], BP["TCL1A"])], "1R1G1F", "T-PLL.", "standard", "inv(14)(q11.2q32.1)")],
      ["LPH 046"], ["IVDD"], cyto("cytocell-tcl1-breakapart-fish-probe")))

add(P("c-trb-ba", C, "CytoCell TCRB (TRB) Breakapart", "Lymphoid", "BA",
      [comp("Centromeric TRB", "R", TS(a5("TRB"), -100, 209), "~209 kb", "centromeric end of TRB, D7S2473"),
       comp("Telomeric of TRB", "G", TS(a3("TRB"), 20, 133), "~133 kb", "SHGC-81465")],
      [NORMAL("2F", src="standard"), pat("TRB rearranged", [T(BP["TRB"], TS(a3("TLX1"), 30, 0))], "1R1G1F", "T-ALL.", "standard", "t(7;10)(q34;q24.3)")],
      ["LPH 048"], ["IVDD"], cyto("cytocell-tcrb-trb-breakapart-fish-probe")))

add(P("c-tlx1-ba", C, "CytoCell TLX1 Breakapart", "Lymphoid", "BA",
      [comp("Centromeric of TLX1", "R", CS(a5("TLX1"), 0, 179), "~179 kb", "KAZALD1, D10S1629"),
       comp("Telomeric of TLX1", "G", TS(a3("TLX1"), 0, 124), "~124 kb", "LBX1")],
      [NORMAL("2F", src="standard"), pat("TLX1 rearranged, t(10;14)", [T(BP["TLX1"], BP["TRA"])], "1R1G1F", "T-ALL.", "standard", "t(10;14)(q24.3;q11.2)")],
      ["LPH 049"], ["IVDD"], cyto("cytocell-tlx1-breakapart-fish-probe")))

add(P("c-tlx3-ba", C, "CytoCell TLX3 Breakapart", "Lymphoid", "BA",
      [comp("Centromeric of TLX3", "R", CS(a5("TLX3"), 0, 162), "~162 kb", "RANBP17"),
       comp("Telomeric of TLX3", "G", TS(a3("TLX3"), 0, 173), "~173 kb", "NPM1, FGF18")],
      [NORMAL("2F", src="standard"), pat("TLX3 rearranged, t(5;14)", [T(BP["TLX3"], bp("14", 98500000))], "1R1G1F", "T-ALL (BCL11B partner).", "standard", "t(5;14)(q35.1;q32.2)")],
      ["LPH 050"], ["IVDD"], cyto("cytocell-tlx3-breakapart-fish-probe")))

add(P("c-igh-ccnd3", C, "CytoCell IGH/CCND3 Translocation, Dual Fusion", "Myeloma", "DF",
      [comp("CCND3", "R", [TS(a3("CCND3"), 0, 158), CS(("e", "TAF8"), 0, 244)], "158 + 244 kb", "telomeric of CCND3 (D6S1017) and centromeric of TAF8"), C_IGH_G],
      [NORMAL("2R2G", src="standard"), pat("t(6;14)", [T(BP["CCND3"], BP["IGH_sw"])], "1R1G2F", None, "standard", "t(6;14)(p21.1;q32.3)")],
      ["LPH 040"], ["IVDD"], cyto("cytocell-ighccnd3-translocation-dual-fusion-fish-probe")))

add(P("c-igh-mafb", C, "CytoCell IGH/MAFB Translocation, Dual Fusion", "Myeloma", "DF",
      [comp("MAFB", "R", [CS(a3("MAFB"), 30, 600), TS(a3("MAFB"), -10, 160)], "195 + 161 + 401 kb", "either side of the breakpoint region"), C_IGH_G],
      [NORMAL("2R2G", src="standard"), pat("t(14;20)", [T(BP["IGH_sw"], BP["MAFB"])], "1R1G2F", None, "standard", "t(14;20)(q32.3;q12)")],
      ["LPH 044"], ["IVDD"], cyto("cytocell-ighmafb-translocation-dual-fusion-fish-probe")))

add(P("m-pcm1-jak2", M, "XL t(8;9) PCM1/JAK2 DF Translocation/Dual Fusion Probe", "Myeloid", "DF",
      [comp("PCM1", "O", TS(a5("PCM1"), -250, 600), "approx.", "PCM1 region"), comp("JAK2", "G", TS(a5("JAK2"), -250, 600), "approx.", "JAK2 region")],
      [NORMAL("2O2G", "Two green (2G) and two orange (2O).", "vendor"),
       pat("t(8;9)", [T(BP["PCM1"], BP["JAK2"])], "1O1G2F", "One green, one orange and two green-orange fusions.", "vendor", "t(8;9)(p22;p24.1)")],
      ["D-5120-100-OG"], ["IVDR"], mets("d-5120-100-og"), "Probe sizes are not given on the product page; positions here are illustrative."))

add(P("m-nup98-nsd1", M, "XL t(5;11) NSD1/NUP98 DF Translocation/Dual Fusion Probe", "Myeloid", "DF",
      [comp("NSD1", "G", TS(a5("NSD1"), -250, 600), "approx.", None), comp("NUP98", "O", TS(a5("NUP98"), -250, 600), "approx.", None)],
      [NORMAL("2O2G", "Two green (2G) and two orange (2O).", "vendor"),
       pat("t(5;11)", [T(BP["NSD1"], BP["NUP98"])], "1O1G2F", None, "vendor", "t(5;11)(q35.3;p15.4)")],
      ["D-5141-100-OG"], ["IVDR"], mets("d-5141-100-og"), "Probe sizes are not given on the product page; positions here are illustrative."))
