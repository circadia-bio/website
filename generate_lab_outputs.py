#!/usr/bin/env python3
"""Pre-render script: pulls every record in the lab's Zenodo Community and
renders a card grid into lab-outputs.qmd, between the
`<!-- LAB_OUTPUTS:START -->` / `<!-- LAB_OUTPUTS:END -->` sentinels.

Data flow:
  projects/lab_outputs.yml -> the community slug (null until it exists)
  Zenodo community-records API -> title, creators, date, DOI, resource type
  .lab_outputs_cache.json  -> last-known-good fallback if the API call fails
  lab-outputs.qmd          -> card grid injected in place

If no community slug is configured yet, this script is a no-op: the
fallback message already committed in lab-outputs.qmd is left untouched.
"""

import html
import json
import os
import re
import sys
import urllib.request
import urllib.error
from datetime import datetime

CONFIG_PATH = os.path.join("projects", "lab_outputs.yml")
CACHE_PATH = ".lab_outputs_cache.json"
PAGE_PATH = "lab-outputs.qmd"

START_MARKER = "<!-- LAB_OUTPUTS:START -->"
END_MARKER = "<!-- LAB_OUTPUTS:END -->"

RESOURCE_TYPE_LABELS = {
    "software": "Software",
    "publication": "Publication",
    "dataset": "Dataset",
    "poster": "Poster",
    "presentation": "Presentation",
    "other": "Other",
}


def read_community_slug(path):
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line.startswith("community_slug:"):
                value = line.split(":", 1)[1].strip().strip('"').strip("'")
                if value in ("", "null", "~"):
                    return None
                return value
    return None


def fetch_community_records(slug):
    """Fetch all records in the community, newest first. Returns None on
    any failure so the caller can fall back to cache.

    Unauthenticated Zenodo API requests cap page size at 25, so this pages
    through results (bounded to a handful of pages -- a lab's total output
    is small; raise MAX_PAGES if it ever isn't)."""
    MAX_PAGES = 8
    records = []
    url = f"https://zenodo.org/api/communities/{slug}/records?size=25&sort=newest"

    for _ in range(MAX_PAGES):
        try:
            with urllib.request.urlopen(url, timeout=15) as resp:
                data = json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ValueError) as e:
            print(f"  ! could not fetch community '{slug}': {e}", file=sys.stderr)
            return None if not records else records

        hits = data.get("hits", {}).get("hits", [])
        for h in hits:
            metadata = h.get("metadata", {})
            creators = metadata.get("creators", []) or []
            author_names = [c.get("name", "") for c in creators if c.get("name")]
            resource_type = (metadata.get("resource_type") or {}).get("type", "other")

            records.append({
                "title": metadata.get("title", "Untitled"),
                "authors": author_names,
                "date": metadata.get("publication_date"),
                "doi": h.get("doi"),
                "url": (h.get("links", {}) or {}).get("self_html") or (h.get("links", {}) or {}).get("html"),
                "resource_type": resource_type,
            })

        next_url = data.get("links", {}).get("next")
        if not next_url:
            break
        url = next_url

    return records


def load_cache():
    if os.path.exists(CACHE_PATH):
        try:
            with open(CACHE_PATH, encoding="utf-8") as f:
                return json.load(f)
        except (OSError, ValueError):
            return None
    return None


def save_cache(records):
    with open(CACHE_PATH, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2, sort_keys=True)
        f.write("\n")


def format_authors(authors):
    if not authors:
        return ""
    if len(authors) == 1:
        return authors[0]
    if len(authors) <= 3:
        return ", ".join(authors[:-1]) + " & " + authors[-1]
    return ", ".join(authors[:2]) + f" et al. ({len(authors)} authors)"


def format_date(iso_date):
    if not iso_date:
        return ""
    try:
        d = datetime.strptime(iso_date, "%Y-%m-%d").date()
        return d.strftime("%-d %b %Y")
    except ValueError:
        return iso_date


def esc(s):
    return html.escape(s or "", quote=True)


def card_html(record):
    rtype = record["resource_type"]
    rlabel = RESOURCE_TYPE_LABELS.get(rtype, rtype.title() if rtype else "Other")
    authors = format_authors(record["authors"])
    date_str = format_date(record["date"])
    doi_badge = ""
    if record["doi"]:
        doi_url = f"https://doi.org/{record['doi']}"
        doi_badge = (
            f'<a class="lo-doi" href="{esc(doi_url)}" target="_blank" rel="noopener">'
            f'<svg viewBox="0 0 24 24" width="9" height="9" fill="none" stroke="currentColor" '
            f'stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round">'
            f'<path d="M7 17L17 7M9 7h8v8"/></svg>DOI</a>'
        )

    url = record["url"] or "#"
    return (
        '<a class="lo-card" href="{url}" target="_blank" rel="noopener">'
        '<span class="lo-type lo-type-{rtype}">{rlabel}</span>'
        '<p class="lo-title">{title}</p>'
        '<p class="lo-meta">{authors}{sep}{date}</p>'
        '<span class="lo-footer">{doi_badge}</span>'
        '</a>'
    ).format(
        url=esc(url),
        rtype=esc(rtype),
        rlabel=esc(rlabel),
        title=esc(record["title"]),
        authors=esc(authors),
        sep=" · " if authors and date_str else "",
        date=esc(date_str),
        doi_badge=doi_badge,
    )


def grid_html(records):
    if not records:
        return (
            '<p class="lo-empty">No outputs in the community yet. '
            "Records added to the Circadia Lab Zenodo community will appear here automatically.</p>"
        )
    cards = "\n    ".join(card_html(r) for r in records)
    return f'<div class="lo-grid">\n    {cards}\n  </div>'


def update_page(records):
    with open(PAGE_PATH, encoding="utf-8") as f:
        content = f.read()

    start = content.find(START_MARKER)
    end = content.find(END_MARKER)
    if start == -1 or end == -1:
        raise RuntimeError(
            f"Could not find {START_MARKER} / {END_MARKER} sentinels in {PAGE_PATH}"
        )

    before = content[: start + len(START_MARKER)]
    after = content[end:]

    injected = "\n```{=html}\n" + grid_html(records) + "\n```\n"

    with open(PAGE_PATH, "w", encoding="utf-8") as f:
        f.write(before + injected + after)


def main():
    slug = read_community_slug(CONFIG_PATH)
    if not slug:
        print(f"! no community_slug set in {CONFIG_PATH}, skipping lab outputs (fallback message left as-is)")
        return

    records = fetch_community_records(slug)
    if records is None:
        cached = load_cache()
        if cached is not None:
            print("  using cached community records")
            records = cached
        else:
            print("! could not fetch community and no cache available, leaving page untouched")
            return
    else:
        save_cache(records)

    update_page(records)
    print(f"✓ {PAGE_PATH} updated ({len(records)} community records)")


if __name__ == "__main__":
    main()
