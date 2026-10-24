#!/usr/bin/env python3
"""Rebuild the website pages from src/ templates and data/ JSON.

Usage:
    python build.py                  # writes the pages into the repo root (the GitHub Pages site)
    python build.py --out DIR        # writes the same pages into another folder

Each page is a single self-contained HTML file: the template gets the shared CSS/JS, the
app bar, the FAQ and the data embedded, so the pages work offline and on GitHub Pages
without a server. Edit src/ or data/, then run this script; never edit the pages by hand.
Before building, tools/check_probes.mjs is run (if node is available) so a FISH pattern
that no longer matches its expected signal code stops the build.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
DATA = ROOT / "data"

# page file -> (template, tab label)
PAGES = {
    "index.html": ("idiograms.template.html", "Idiograms"),
    "simulator.html": ("simulator.template.html", "Microscope"),
    "abnormalities.html": ("abnormalities.template.html", "Abnormalities"),
    "fish.html": ("fish.template.html", "FISH probes"),
    "catalogue.html": ("catalogue.template.html", "Probe catalogue"),
}
SHORT = {"index.html": "Ideo", "simulator.html": "Scope", "abnormalities.html": "Abn", "fish.html": "FISH", "catalogue.html": "Catalogue"}  # tab labels on phones
FONTS = ('<link rel="preconnect" href="https://fonts.googleapis.com">\n'
         '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>\n'
         '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Geist:wght@400;500;600;700&family=Geist+Mono:wght@400;500&display=swap">\n')


def compact(obj):
    return json.dumps(obj, separators=(",", ":"), ensure_ascii=False)


def load_data():
    j = lambda name: json.loads((DATA / name).read_text(encoding="utf-8"))
    probes = j("probes.json")
    cat = j("catalogue.json")
    # catalogue rows as arrays to keep the page small
    keys = ["vendor", "name", "codes", "status", "chroms", "targets", "labels", "url", "fish", "custom", "slug_name", "area", "source"]
    rows = [[it.get(k) for k in keys] for it in cat["items"]]
    return {
        "__DATA__": compact(j("bands.json")["chromosomes"]),
        "__GENES__": compact(j("landmark-genes.json")["chromosomes"]),
        "__BPGENES__": compact(j("breakpoint-genes.json")["chromosomes"]),
        "__PROBES__": compact(probes["probes"]),
        "__PROBEGENES__": compact(probes["genes"]),
        "__CATALOGUE__": compact({"keys": keys, "rows": rows, "counts": cat["counts"]}),
    }


def nav(current, links):
    tabs = "".join(
        f'<a href="{links[page]}" data-chr-link{" aria-current=\"page\"" if page == current else ""}><span class="lg">{label}</span><span class="sm">{SHORT[page]}</span></a>'
        for page, (_, label) in PAGES.items())
    return ('<nav class="appbar" aria-label="Site"><div class="appbar-in">'
            f'<a class="brand" href="{links["index.html"]}"><span class="brand-mark" aria-hidden="true"></span><span class="brand-name">Idiograms</span></a>'
            f'<div class="tabs">{tabs}</div>'
            '<button class="faq-btn" id="faqBtn" type="button" aria-haspopup="dialog"><svg viewBox="0 0 24 24" width="15" height="15" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><circle cx="12" cy="12" r="10"/><path d="M9.09 9a3 3 0 0 1 5.83 1c0 2-3 3-3 3"/><path d="M12 17h.01"/></svg><span>FAQ</span></button>'
            '</div></nav>')


def wrap(title, head_extra, body):
    return ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
            f'<title>{title}</title>\n{FONTS}{head_extra}</head>\n<body>\n{body}\n</body>\n</html>\n')


def render(page, template, data, links):
    t = template
    title = t[t.index("<title>") + 7:t.index("</title>")]
    t = t[t.index("</title>") + 8:]
    css_end = t.index("</style>") + len("</style>")
    head, body = t[:css_end], t[css_end:]
    head = head.replace("/*__SHARED_CSS__*/", (SRC / "_shared.css").read_text(encoding="utf-8"))
    body = body.replace("<!--__NAV__-->", nav(page, links))
    body = body.replace("<!--__FAQ__-->", (SRC / "_faq.html").read_text(encoding="utf-8"))
    body = body.replace("/*__SHARED_JS__*/", (SRC / "_shared.js").read_text(encoding="utf-8"))
    body = body.replace("/*__SCOPE_JS__*/", (SRC / "_scope.js").read_text(encoding="utf-8") if (SRC / "_scope.js").exists() else "")
    body = body.replace("/*__FISH_JS__*/", (SRC / "_fish.js").read_text(encoding="utf-8"))
    for k, v in data.items():
        body = body.replace(k, v)
    for k, v in links.items():
        body = body.replace(f"__HREF:{k}__", v)
    left = [tok for tok in ("__DATA__", "__GENES__", "__PROBES__", "__CATALOGUE__", "__HREF:", "/*__", "<!--__") if tok in body]
    if left:
        raise SystemExit(f"{page}: unfilled placeholders {left}")
    return wrap(title, "", head + body)


def check_probes():
    if not shutil.which("node"):
        print("node not found: skipping the FISH pattern check")
        return
    r = subprocess.run(["node", str(ROOT / "tools" / "check_probes.mjs")], capture_output=True, text=True)
    print(r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr)
    if r.returncode:
        raise SystemExit(r.stdout + r.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(ROOT), help="output folder (default: repo root)")
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    check_probes()
    data = load_data()
    # Version tag on links between pages: it changes whenever a source or the data changes,
    # so following a tab always fetches the current page instead of a stale cached copy.
    digest = hashlib.sha1()
    for path in sorted(SRC.glob("*")) + sorted(DATA.glob("*.json")):
        digest.update(path.read_bytes())
    v = "?v=" + digest.hexdigest()[:8]
    links = {page: page + v for page in PAGES}
    for page, (template_name, _) in PAGES.items():
        html = render(page, (SRC / template_name).read_text(encoding="utf-8"), data, links)
        (out_dir / page).write_text(html, encoding="utf-8")
        print(f"wrote {out_dir / page} ({len(html) // 1024} KB)")


if __name__ == "__main__":
    main()
