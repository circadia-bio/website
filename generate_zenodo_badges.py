#!/usr/bin/env python3
"""Pre-render script: resolves Zenodo concept DOIs for tracked packages and
injects version/DOI/release-date/citation fields into the products carousel
array in index.qmd, between the `// PRODUCTS:START` / `// PRODUCTS:END`
sentinels inside the <script> block.

Data flow:
  projects/zenodo.yml    -> which packages are tracked, their repo + concept DOI
  doi.org redirect + Zenodo REST API -> current version + publication date
  .zenodo_cache.json     -> last-known-good fallback if the network call fails
  index.qmd              -> products array gets extra fields injected in place
"""

import json
import os
import re
import sys
import urllib.request
import urllib.error
from datetime import datetime, date

ZENODO_YML = os.path.join("projects", "zenodo.yml")
CACHE_PATH = ".zenodo_cache.json"
INDEX_QMD = "index.qmd"
GITHUB_ORG = "circadia-bio"

START_MARKER = "// PRODUCTS:START"
END_MARKER = "// PRODUCTS:END"


# ── Minimal YAML reading (avoids adding a PyYAML dependency to CI) ─────────────
# projects/zenodo.yml is a flat list of simple key: value maps, so a small
# hand-rolled parser is enough and keeps this script dependency-free.

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


# ── Zenodo / DOI resolution ─────────────────────────────────────────────────

RECORD_ID_RE = re.compile(r"/records?/(\d+)")


def resolve_concept_doi(doi):
    """Resolve a Zenodo concept DOI to {version, doi_url, released_date}.

    Two steps, no Zenodo API key needed:
      1. Follow the doi.org redirect for the concept DOI to find the latest
         version's record ID (Zenodo always redirects concept DOIs to the
         current latest record's landing page).
      2. Call the plain Zenodo REST API for that record, which exposes
         metadata.version (the free-text version string, usually matching
         the git tag) and metadata.publication_date.

    Returns None on any failure so the caller falls back to cache.
    """
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


# ── Load cache ───────────────────────────────────────────────────────────────

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


# ── Build per-product data ──────────────────────────────────────────────────

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

        result[name] = {
            "version": version,
            "doi_url": doi_url,
            "released_date": relative_date(released_date),
            "citation_url": citation_url,
        }
    return result


# ── Inject fields into the existing product object literals in index.qmd ────

# Matches one product literal like: { href: "...", ..., tag: "R pkg" }
OBJECT_RE = re.compile(r"\{[^{}]*\bname:\s*\"([^\"]+)\"[^{}]*\}")


def js_string(value):
    if value is None:
        return "null"
    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def inject_fields(match, product_data):
    obj_text = match.group(0)
    name = match.group(1)
    data = product_data.get(name)

    # Strip any previously-injected fields so re-runs are idempotent.
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
    )
    # Insert right before the closing brace.
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


# ── Main ─────────────────────────────────────────────────────────────────────

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
