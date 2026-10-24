#!/usr/bin/env python3
"""Build data/catalogue.json: one searchable list of human FISH probes from several vendors.

    python tools/make_catalogue.py RAW_DIR

RAW_DIR holds the files saved from each vendor's public website (October 2026):

  Vysis (Abbott)    amd_items.txt        rows of AMD_US_Product_List.pdf ("Abbott Molecular Ordering
                                         Information", 10/03/2019), via pdftotext -layout
                    abb_int.json         rendered text of molecular.abbott/int/en/chromosome-main/<n>
                    abb_products.json    per-product blocks parsed from the same pages
  CytoCell (OGT)    ogt_all.json         ogt.com /umbraco/Api/SearchApi/ProductDetailsPageSearchAdvanced
  ZytoVision        zyto.json            titles and intended-use text of zytovision.com/products/zytolight/*
  MetaSystems       ms_details.json      metasystems-probes.com/en/probes/* product pages
  Leica (Kreatech)  leica_probe_slugs.txt  product slugs of shop.leicabiosystems.com/.../ish-probes-molecular-pathology
  Empire Genomics   eg_product-sitemap.xml, eg_product-sitemap2.xml  empiregenomics.com product sitemaps

Names, codes, regions and labels are copied as published. Where a vendor publishes only a URL
slug (Leica, Empire Genomics) the name is rebuilt from the slug and flagged.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "catalogue.json"
ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]
BAND = re.compile(r"\b(?:[1-9]|1\d|2[0-2]|X|Y)[pq]\d{1,2}(?:\.\d{1,3})?(?:-[pq]?\d{1,2}(?:\.\d{1,3})?)?\b")


def chroms_of(*texts):
    found = []
    for t in texts:
        for m in BAND.findall(t or ""):
            c = re.match(r"(\d+|X|Y)", m).group(1)
            if c not in found:
                found.append(c)
    return found


def bands_of(*texts):
    out = []
    for t in texts:
        for m in BAND.findall(t or ""):
            if m not in out:
                out.append(m)
    return out


COLORS = {"orange": "O", "green": "G", "red": "R", "aqua": "A", "blue": "B", "gold": "Y", "yellow": "Y", "spectrumorange": "O",
          "spectrumgreen": "G", "spectrumred": "R", "spectrumaqua": "A", "spectrumgold": "Y", "spectrumblue": "B", "texasred": "R"}


def colors_in(text):
    out = []
    for w in re.findall(r"Spectrum(Orange|Green|Red|Aqua|Gold|Blue)|\b(orange|green|red|aqua|blue|gold|yellow)\b", text or "", re.I):
        k = (w[0] or w[1]).lower()
        c = COLORS.get(k)
        if c and c not in out:
            out.append(c)
    return out


def vysis(raw, probes_by_name):
    rows = {}
    skip = re.compile(r"ProbeChek|Control Slides|COT-1|Control DNA|Target Slides|Hybridization|Buffer|Wash|DAPI|Counterstain|"
                      r"Denaturation|Pretreatment|Paraffin|Protease|Formalin|Kit Accessories|Analyzer|ThermoBrite|Instrument|Software|Imaging", re.I)

    def add(name, code=None, unit=None, status=None, chrom=None, url=None, source=None):
        name = re.sub(r"\s+", " ", name).strip()
        m = re.search(r"\((ASR|RUO|IVD|CE|GPR|US-IVD)\)\s*$", name)
        if m:
            status = status or m.group(1)
            name = name[:m.start()].strip()
        if skip.search(name):
            return
        key = re.sub(r"[^a-z0-9]", "", name.lower())
        r = rows.setdefault(key, {"vendor": "Vysis (Abbott)", "name": name, "codes": [], "status": [], "chroms": [], "url": url, "source": source})
        full = code + (f" ({unit})" if unit else "") if code else None
        if full and full not in r["codes"]:
            r["codes"].append(full)
        if status and status not in r["status"]:
            r["status"].append(status)
        for c in chrom or []:
            if c not in r["chroms"]:
                r["chroms"].append(c)
        if url and not r["url"]:
            r["url"] = url

    for line in (raw / "amd_items.txt").read_text().splitlines():
        _, code, unit, name = line.split("\t")
        add(name, code, unit, None, None, None, "Abbott US ordering list, 2019")
    pages = json.loads((raw / "abb_int.json").read_text())
    for c, v in pages.items():
        t = v["text"]
        if "VYSIS FISH CHROMOSOME SEARCH:" not in t:
            continue
        seg = t.split("VYSIS FISH CHROMOSOME SEARCH:")[1]
        tag = None
        for line in [x.strip() for x in seg.split("\n") if x.strip()]:
            if line in ("CE", "ASR", "RUO", "US-IVD", "IVD", "GPR"):
                tag = line
                continue
            if re.match(r"(Vysis|TelVysion|UroVysion|AneuVysion|PathVysion|MultiVysion)", line):
                add(line, None, None, tag, [c.upper()], f"https://www.molecular.abbott/int/en/chromosome-main/{c}", "molecular.abbott chromosome search")
    for name, v in json.loads((raw / "abb_products.json").read_text()).items():
        m = re.match(r"(CE|ASR|RUO|US-IVD)?\s*(.+)", name)
        if m and re.match(r"(Vysis|AneuVysion|MultiVysion|UroVysion|PathVysion|TelVysion)", m.group(2)):
            add(m.group(2), None, None, m.group(1), [x.upper() for x in v["pages"]], f"https://www.molecular.abbott/int/en/chromosome-main/{v['pages'][0]}", "molecular.abbott chromosome search")
    out = []
    for r in rows.values():
        r["chroms"] = r["chroms"] or chroms_of(r["name"])
        r["targets"] = ", ".join(bands_of(r["name"]))
        r["labels"] = colors_in(r["name"])
        p = probes_by_name.get(norm(r["name"]))
        if p:
            r["fish"] = p
        out.append(r)
    return out


def cytocell(raw):
    out = []
    for p in json.loads((raw / "ogt_all.json").read_text()):
        if not p.get("region"):
            continue
        custom = p["productType"] == "Custom_Probe"
        out.append({"vendor": "CytoCell (OGT)", "name": p["title"], "codes": p.get("productNumber") or [],
                    "status": [p["regulatoryStatus"]] if p.get("regulatoryStatus") else [],
                    "chroms": [re.sub(r"<[^>]+>", "", x) for x in re.findall(r"<span>([^<]+)</span>", p.get("friendlyChromosome") or "")],
                    "targets": ", ".join(p["region"]), "labels": [COLORS.get(x.lower(), x[:1]) for x in p.get("labels") or []],
                    "area": " / ".join(x for x in [p.get("friendlyApplication"), p.get("friendlyAreaOfInterest")] if x),
                    "url": "https://www.ogt.com" + p["url"], "custom": custom,
                    "source": "ogt.com product search" + (" (myProbes custom design)" if custom else "")})
    return out


def zytovision(raw):
    out = []
    for o in json.loads((raw / "zyto.json").read_text()):
        name = o["title"]
        slug = o["f"].split("/")[-1].replace(".html", "")
        if not re.search(r"Probe", name) or "Implementation Kit" in name or not re.match(r"(zytolight|flexish)_", slug):
            continue
        code = slug.split("_", 1)[1].upper()
        if not re.match(r"Z-\d{4}", code):
            continue
        desc = o.get("desc", "")
        out.append({"vendor": "ZytoVision", "name": name, "codes": [code], "status": [], "chroms": chroms_of(desc, name),
                    "targets": ", ".join(bands_of(desc)), "labels": [],
                    "url": "https://zytovision.com/products/" + slug.replace("_", "/"), "source": "zytovision.com"})
    return out


def metasystems(raw):
    out = []
    for k, v in json.loads((raw / "ms_details.json").read_text()).items():
        fam = k.split("_")[0]
        if fam in ("xmp", "xrp", "xod") or "Mouse" in v["title"] or "Rat" in v["title"] or "Customized" in v["title"]:
            continue
        name = v["title"].split(" - ")[0].strip()
        kind = v["title"].split(" - ")[1].strip() if " - " in v["title"] else ""
        desc = re.split(r"Probe maps for selected products|Need assistance", v.get("desc") or "")[0]
        labels = colors_in(" ".join((v.get("labels") or "").split()[:6]))
        reg = (v.get("reg") or "").split()[0] if v.get("reg") else ""
        chrom = chroms_of(desc, name)
        out.append({"vendor": "MetaSystems", "name": name + (f" ({kind})" if kind else ""), "codes": [k.split("_")[1].upper()],
                    "status": [reg] if reg in ("IVDR", "IVDD", "RUO", "CE") else [], "chroms": chrom,
                    "targets": ", ".join(bands_of(desc, name)), "labels": labels, "note": desc[:400],
                    "url": f"https://metasystems-probes.com/en/probes/{fam}/{k.split('_')[1]}/", "source": "metasystems-probes.com"})
    return out


TOK = {"tc": "TC", "dc": "DC", "es": "ES", "se": "SE", "xl": "XL", "ivd": "IVD", "ruo": "RUO", "asr": "ASR"}


def humanize(slug):
    s = slug.replace("pid-", "")
    status = None
    m = re.match(r"(ivd|ruo)-(.*)", s)
    if m:
        status, s = m.group(1).upper(), m.group(2)
    m = re.match(r"(.*)-(asr|ivd)$", s)
    if m:
        s, status = m.group(1), (status or m.group(2).upper())
    s = s.replace("triple-color", "triple_color").replace("dual-color", "dual_color").replace("break-tissue", "break_tissue")
    toks = s.split("-")
    words = []
    for t in toks:
        if t in TOK:
            words.append(TOK[t])
        elif t == "break":
            words.append("Break")
        elif t == "break_tissue":
            words.append("Break (tissue)")
        elif t in ("triple_color", "dual_color"):
            words.append(t.replace("_", " ").title())
        elif re.fullmatch(r"\d{1,2}[pq]", t):
            words.append(t)
        elif t == "dic" and words == [] and "920" in toks:
            words.append("dic(9;20)")
        elif t == "920" and words == ["dic(9;20)"]:
            continue
        elif re.fullmatch(r"\d{1,2}[pq]\d+(\.\d+)?|[xy][pq]\d+", t):
            words.append(t[0].upper() + t[1:] if t[0] in "xy" else t)
        elif re.fullmatch(r"\d+", t) and words and words[-1] == "SE":
            words[-1] = "SE " + t
        else:
            words.append(t.upper())
    out = []
    for w in words:  # join consecutive loci with "/"
        if out and re.fullmatch(r"[A-Z0-9.]+|\d{1,2}[pq][\d.]+|[XY][pq]\d+|SE \d+", w) and re.fullmatch(r"[A-Z0-9.]+|\d{1,2}[pq][\d.]+|[XY][pq]\d+|SE \d+", out[-1]) and w not in ("TC", "DC", "ES", "XL"):
            out[-1] += "/" + w
        else:
            out.append(w)
    return " ".join(out), status


def leica(raw):
    rows = {}
    for slug in (raw / "leica_probe_slugs.txt").read_text().split():
        if not re.search(r"-(ivd|asr)$|^pid-(ivd|ruo)-", slug) or "bond" in slug or "ancillar" in slug or "pretreatment" in slug or "kits" in slug or "media" in slug or "mouse" in slug:
            continue
        name, status = humanize(slug)
        key = name.lower()
        r = rows.setdefault(key, {"vendor": "Leica (Kreatech)", "name": name, "codes": [], "status": [], "chroms": chroms_of(name),
                                  "targets": ", ".join(bands_of(name)), "labels": [], "slug_name": True,
                                  "url": "https://www.leicabiosystems.com/ihc-ish/kreatech-fish-probes/", "source": "shop.leicabiosystems.com (name rebuilt from the product URL)"})
        if status and status not in r["status"]:
            r["status"].append(status)
    return list(rows.values())


def empire(raw):
    out = []
    for f in ("eg_product-sitemap.xml", "eg_product-sitemap2.xml"):
        for u in re.findall(r"<loc>([^<]+)</loc>", (raw / f).read_text()):
            m = re.search(r"/product/([^/]+)/?$", u)
            if not m:
                continue
            slug = m.group(1)
            if not re.search(r"fish-probe|break-apart|enumeration|fusion", slug) or re.search(r"oncore|wash|dapi|buffer|mouse|cish|chromogen|kit-ancillary", slug):
                continue
            words = []
            for t in slug.split("-"):
                if t in ("fish", "probe", "and"):
                    words.append({"fish": "FISH", "probe": "Probe", "and": "and"}[t])
                elif t in ("break", "apart", "green", "orange", "red", "gold", "aqua", "blue", "telomere", "centromere", "extended", "fusion", "dual", "color", "tri", "enumeration", "amplification", "deletion"):
                    words.append(t.title())
                else:
                    words.append(t.upper())
            name = " ".join(words).replace("Break Apart", "Break-Apart")
            out.append({"vendor": "Empire Genomics", "name": name, "codes": [], "status": [], "chroms": [], "targets": "",
                        "labels": colors_in(name), "slug_name": True, "url": u, "source": "empiregenomics.com (name rebuilt from the product URL)"})
    return out


STOP = re.compile(r"\b(vysis|lsi|fish|probes?|kit|translocation|rearrangement|dual|color|colour|spectrum\w*|the|and)\b")


def norm(name):
    return re.sub(r"[^a-z0-9]", "", STOP.sub("", name.lower()))


VYSIS_LINKS = [
    (r"ETV6 \(TEL\)/RUNX1 \(AML1\) ES", "v-etv6-runx1-es"), (r"LSI (?:TCF3|PBX1) ", "v-tcf3-pbx1"),
    (r"D13S319|LSI \(?13q34", "v-d13s319"), (r"LSI FUS \((?:Cen|Tel)\)", "v-fus-ba"), (r"LSI SS18 \((?:Cen|Tel)\)", "v-ss18-ba"),
    (r"LSI FOXO1 \((?:Cen|Tel)\)", "v-foxo1-ba"), (r"LSI NTRK1 \((?:Cen|Tel)\)", "v-ntrk1-ba"), (r"LSI ERG \((?:Cen|Tel)\)", "v-erg-ba"),
    (r"LSI RET \((?:Cen|Tel)\)", "v-ret-ba"), (r"LSI ROS1 \((?:Cen|Tel)\)", "v-ros1-ba"), (r"LSI DDIT3 \((?:Cen|Tel)\)", "v-ddit3-ba"),
    (r"LSI ETV6 \((?:CEN|TEL)\) Spectrum", "v-etv6-ba"), (r"LSI RARA \(Tel\)", "v-rara-ba"),
    (r"MultiVysion", "v-multivysion"), (r"LSI STS ", "v-sts"), (r"LSI SRY", "v-sry"), (r"CEP X .*Y .*Kit", "v-cepxy"),
    (r"Vysis CEP ?(\d+|X|Y)\b", "v-cep{0}-single"), (r"TelVysion (\d+|X)[pq]", "v-telvysion-{0}"),
    (r"PathVysion|LSI HER-?2", "v-pathvysion"), (r"LSI ATM ", "v-atm-cep11"), (r"LSI CDKN2A ", "v-cdkn2a-cep9"),
    (r"LSI MYC ", "v-myc-ba"), (r"LSI EGFR ", "v-egfr-cep7"), (r"LSI MDM2 ", "v-mdm2-cep12"), (r"LSI N-MYC", "v-nmyc-cep2"),
    (r"LSI PTEN ", "v-pten-cep10"), (r"LSI EGR1", "v-egr1"), (r"Cri-du-Chat", "v-egr1"), (r"LSI D5S23, D5S721 Spectrum", "v-egr1"),
    (r"LSI CKS1B|LSI CDKN2C", "v-cks1b-cdkn2c"), (r"LSI BIRC3 ", "v-birc3-malt1"), (r"LSI MALT1 Spectrum", "v-igh-malt1"),
    (r"LSI 4q12 ", "v-4q12"), (r"LSI DEK |LSI NUP214 ", "v-dek-nup214"), (r"LSI RUNX1 Spectrum", "v-runx1-runx1t1"),
    (r"LSI 22 \(BCR\)", "v-bcr-abl-es"), (r"LSI 1p36 Spectrum|1p36 Microdeletion", "v-1p36-microdel"), (r"DiGeorge.*TelVysion", "v-digeorge-tuple1"),
    (r"Prader-Willi", "v-pws-snrpn"), (r"BRAF|MET Spectrum|PIK3CA|KRAS", None),
]
VYSIS_LINKS = [(a, b) for a, b in VYSIS_LINKS if b]


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    raw = Path(sys.argv[1])
    probes = json.loads((ROOT / "data" / "probes.json").read_text())["probes"]
    by_name = {norm(p["name"]): p["id"] for p in probes}
    items = vysis(raw, by_name) + cytocell(raw) + zytovision(raw) + metasystems(raw) + leica(raw) + empire(raw)
    # link catalogue rows to FISH-page probes by name where possible (same vendor only)
    prefix = {"Vysis (Abbott)": "v-", "CytoCell (OGT)": "c-", "MetaSystems": "m-"}
    ids = {p["id"] for p in probes}
    for it in items:
        if it.get("fish") and (not it["fish"].startswith(prefix.get(it["vendor"], "?")) or it["fish"] not in ids):
            del it["fish"]
        if "fish" not in it and it["vendor"] in prefix:
            k = norm(it["name"])
            if k in by_name and by_name[k].startswith(prefix[it["vendor"]]):
                it["fish"] = by_name[k]
        if "fish" not in it and it["vendor"] == "Vysis (Abbott)":
            for pat, pid in VYSIS_LINKS:
                m = re.search(pat, it["name"])
                if m:
                    pid = pid.format(*[g.lower() for g in m.groups()])
                    it["fish"] = pid if pid in ids else pid.replace("-single", "")
                    break
        it["chroms"] = sorted(set(it["chroms"]), key=lambda c: ORDER.index(c) if c in ORDER else 99)
    items.sort(key=lambda r: (r["vendor"], r["name"].lower()))
    counts = {}
    for it in items:
        counts[it["vendor"]] = counts.get(it["vendor"], 0) + 1
    OUT.write_text(json.dumps({"_about": "Generated by tools/make_catalogue.py from vendor websites, October 2026.",
                               "counts": counts, "items": items}, separators=(",", ":"), ensure_ascii=False) + "\n")
    print(f"wrote {OUT}: {counts}")


if __name__ == "__main__":
    main()
