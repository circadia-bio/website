#!/usr/bin/env python3
"""Pre-render script: resolves Zenodo concept DOIs for tracked packages and
injects version/DOI/release-date/citation/license fields into the products
carousel array in index.qmd, between the `// PRODUCTS:START` /
`// PRODUCTS:END` sentinels inside the <script> block.

Data flow:
  projects/zenodo.yml    -> which packages are tracked, their repo, concept
                            DOI, and SPDX license id
  doi.org redirect + Zenodo REST API -> current version + publication date
  .zenodo_cache.json     -> last-known-good fallback if the network call fails
  index.qmd              -> products array gets extra fields injected in place
"""

import json
import os
import re
import sys
import unicodedata
import urllib.request
import urllib.error
from datetime import datetime, date

ZENODO_YML = os.path.join("projects", "zenodo.yml")
CACHE_PATH = ".zenodo_cache.json"
INDEX_QMD = "index.qmd"
GITHUB_ORG = "circadia-bio"

START_MARKER = "// PRODUCTS:START"
END_MARKER = "// PRODUCTS:END"


def parse_zenodo_yml(path):
    entries = []
    current = None
    with open(path, encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.split("#", 1)[0].rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("- "):
                if current is not None:
                    entries.append(current)
                current = {}
                line = line[2:]
            stripped = line.strip()
            if not stripped or current is None:
                continue
            if ":" not in stripped:
                continue
            key, _, value = stripped.partition(":")
            key = key.strip()
            value = value.strip().strip('"').strip("'")
            if value == "null" or value == "":
                value = None
            elif value == "true":
                value = True
            elif value == "false":
                value = False
            current[key] = value
    if current is not None:
        entries.append(current)
    return entries


RECORD_ID_RE = re.compile(r"/records?/(\d+)")


def resolve_concept_doi(doi):
    doi_url = f"https://doi.org/{doi}"
    try:
        req = urllib.request.Request(doi_url, method="HEAD")
        with urllib.request.urlopen(req, timeout=10) as resp:
            landing_url = resp.geturl()
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        print(f"  ! could not resolve {doi}: {e}", file=sys.stderr)
        return None

    match = RECORD_ID_RE.search(landing_url)
    if not match:
        print(f"  ! could not extract record id from {landing_url}", file=sys.stderr)
        return None
    record_id = match.group(1)

    api_url = f"https://zenodo.org/api/records/{record_id}"
    try:
        with urllib.request.urlopen(api_url, timeout=10) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ValueError) as e:
        print(f"  ! could not fetch Zenodo record {record_id}: {e}", file=sys.stderr)
        return None

    metadata = data.get("metadata", {})
    return {
        "version": metadata.get("version"),
        "doi_url": f"https://doi.org/{data.get('doi', doi)}",
        "released_date": metadata.get("publication_date"),
    }


def relative_date(iso_date):
    if not iso_date:
        return None
    try:
        d = datetime.strptime(iso_date, "%Y-%m-%d").date()
    except ValueError:
        return None
    days = (date.today() - d).days
    if days < 0:
        return "released " + d.strftime("%-d %b %Y")
    if days == 0:
        return "released today"
    if days == 1:
        return "released yesterday"
    if days < 30:
        return f"released {days} days ago"
    if days < 365:
        months = days // 30
        return f"released {months} month{'s' if months != 1 else ''} ago"
    years = days // 365
    return f"released {years} year{'s' if years != 1 else ''} ago"


def load_cache():
    if os.path.exists(CACHE_PATH):
        try:
            with open(CACHE_PATH, encoding="utf-8") as f:
                return json.load(f)
        except (OSError, ValueError):
            return {}
    return {}


def save_cache(cache):
    with open(CACHE_PATH, "w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2, sort_keys=True)
        f.write("\n")


LICENSE_LABELS = {
    "MIT": "MIT",
    "GPL-2.0-or-later": "GPL-2.0+",
    "GPL-3.0-or-later": "GPL-3.0+",
    "Apache-2.0": "Apache 2.0",
    "BSD-2-Clause": "BSD-2",
    "BSD-3-Clause": "BSD-3",
}


def license_label(spdx_id):
    if not spdx_id:
        return None
    return LICENSE_LABELS.get(spdx_id, spdx_id[:12])


def slugify_ascii(s):
    nfkd = unicodedata.normalize("NFKD", s)
    ascii_str = nfkd.encode("ascii", "ignore").decode("ascii")
    return re.sub(r"[^a-z0-9]", "", ascii_str.lower())


def build_bibtex(name, repo, title, authors_str, version, doi_url, released_date):
    if not authors_str:
        return None

    authors = [a.strip() for a in authors_str.split(";") if a.strip()]
    if not authors:
        return None
    author_field = " and ".join(authors)

    first_family = authors[0].split(",")[0].strip()
    key = f"{slugify_ascii(first_family)}_{slugify_ascii(name)}_"
    year = released_date[:4] if released_date else str(date.today().year)
    key += year

    lines = [f"@software{{{key},"]
    lines.append(f"  author  = {{{author_field}}},")
    lines.append(f"  title   = {{{{{title or name}}}}},")
    lines.append(f"  year    = {{{year}}},")
    if version:
        v = version if str(version).lower().startswith("v") else f"v{version}"
        lines.append(f"  version = {{{v}}},")
    if doi_url:
        doi_bare = doi_url.replace("https://doi.org/", "")
        lines.append(f"  doi     = {{{doi_bare}}},")
    if repo:
        lines.append(f"  url     = {{https://github.com/{GITHUB_ORG}/{repo}}}")
    lines.append("}")
    return "\n".join(lines)


def build_product_data(entries, cache):
    result = {}
    for entry in entries:
        name = entry.get("name")
        if not name:
            continue
        if not entry.get("zenodo_tracked"):
            continue

        repo = entry.get("repo")
        citation_url = (
            f"https://github.com/{GITHUB_ORG}/{repo}/blob/main/CITATION.cff"
            if repo else None
        )
        license_url = (
            f"https://github.com/{GITHUB_ORG}/{repo}/blob/main/LICENSE"
            if repo else None
        )
        concept_doi = entry.get("zenodo_concept_doi")

        version = None
        doi_url = None
        released_date = None

        if concept_doi:
            resolved = resolve_concept_doi(concept_doi)
            if resolved is None:
                cached = cache.get(name)
                if cached:
                    print(f"  using cached data for {name}")
                    resolved = cached
            if resolved:
                version = resolved.get("version")
                doi_url = resolved.get("doi_url")
                released_date = resolved.get("released_date")
                cache[name] = resolved

        entry_license = entry.get("license")
        citation_bibtex = build_bibtex(
            name=name,
            repo=repo,
            title=entry.get("title"),
            authors_str=entry.get("authors"),
            version=version,
            doi_url=doi_url,
            released_date=released_date,
        )

        result[name] = {
            "version": version,
            "doi_url": doi_url,
            "released_date": relative_date(released_date),
            "citation_url": citation_url,
            "citation_bibtex": citation_bibtex,
            "license_label": license_label(entry_license),
            "license_url": license_url if entry_license else None,
        }
    return result


OBJECT_RE = re.compile(r"\{[^{}]*\bname:\s*\"([^\"]+)\"[^{}]*\}")


def js_string(value):
    if value is None:
        return "null"
    escaped = (
        value.replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("\r\n", "\\n")
        .replace("\n", "\\n")
    )
    return f'"{escaped}"'


def inject_fields(match, product_data):
    obj_text = match.group(0)
    name = match.group(1)
    data = product_data.get(name)

    obj_text = re.sub(
        r",\s*zenodoTracked:.*?(?=\})", "", obj_text, flags=re.DOTALL
    )

    if data is None:
        return obj_text

    fields = (
        f", zenodoTracked: true"
        f", version: {js_string(data['version'])}"
        f", doiUrl: {js_string(data['doi_url'])}"
        f", releasedDate: {js_string(data['released_date'])}"
        f", citationUrl: {js_string(data['citation_url'])}"
        f", citationBibtex: {js_string(data['citation_bibtex'])}"
        f", licenseLabel: {js_string(data['license_label'])}"
        f", licenseUrl: {js_string(data['license_url'])}"
    )
    return obj_text[:-1].rstrip() + fields + " }"


def update_index_qmd(product_data):
    with open(INDEX_QMD, encoding="utf-8") as f:
        content = f.read()

    start = content.find(START_MARKER)
    end = content.find(END_MARKER)
    if start == -1 or end == -1:
        raise RuntimeError(
            "Could not find // PRODUCTS:START / // PRODUCTS:END sentinels in index.qmd"
        )

    before = content[: start + len(START_MARKER)]
    block = content[start + len(START_MARKER): end]
    after = content[end:]

    new_block = OBJECT_RE.sub(lambda m: inject_fields(m, product_data), block)

    with open(INDEX_QMD, "w", encoding="utf-8") as f:
        f.write(before + new_block + after)


def main():
    if not os.path.exists(ZENODO_YML):
        print(f"! {ZENODO_YML} not found, skipping Zenodo badge injection")
        return

    entries = parse_zenodo_yml(ZENODO_YML)
    cache = load_cache()
    product_data = build_product_data(entries, cache)
    save_cache(cache)
    update_index_qmd(product_data)

    tracked = len(product_data)
    with_doi = sum(1 for d in product_data.values() if d["version"])
    print(f"✓ index.qmd updated ({tracked} tracked packages, {with_doi} with a Zenodo release)")


if __name__ == "__main__":
    main()
