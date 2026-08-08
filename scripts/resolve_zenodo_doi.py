#!/usr/bin/env python3
"""Called from the zenodo-doi-lookup reusable workflow after a package
release. Looks up the newly-archived Zenodo record for a repo (by matching
its GitHub URL in related_identifiers), and if the corresponding entry in
projects/zenodo.yml still has zenodo_concept_doi: null, fills it in.

Exits 0 with no changes if:
  - the record isn't found yet after retries (Zenodo archiving can lag a
    release by a minute or two -- a later manual run or the next release
    will pick it up)
  - the entry already has a concept DOI set (nothing to do)

Prints "CHANGED=1" to stdout (last line) if projects/zenodo.yml was
modified, so the calling workflow knows whether to open a PR.
"""

import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import urllib.error

ZENODO_YML = os.path.join("projects", "zenodo.yml")
GITHUB_ORG = "circadia-bio"

MAX_ATTEMPTS = 10
RETRY_DELAY_SECONDS = 30


def search_concept_doi(package_query, repo):
    """Search Zenodo for a record whose related_identifiers reference this
    repo, and return its concept DOI. Returns None if not found."""
    url = f"https://zenodo.org/api/records?q={urllib.parse.quote(package_query)}&size=25"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ValueError) as e:
        print(f"  ! search failed: {e}", file=sys.stderr)
        return None

    prefix = f"https://github.com/{GITHUB_ORG}/{repo}/"
    exact = f"https://github.com/{GITHUB_ORG}/{repo}"

    for h in data.get("hits", {}).get("hits", []):
        related = h.get("metadata", {}).get("related_identifiers", []) or []
        for r in related:
            ident = r.get("identifier", "")
            if ident == exact or ident.startswith(prefix):
                return h.get("conceptdoi")
    return None


def find_entry_span(content, package_name):
    """Return (start, end) character offsets of the YAML list entry whose
    `name:` field matches package_name, or None if not found."""
    pattern = re.compile(
        r'(?m)^- name:\s*"' + re.escape(package_name) + r'"\s*$'
    )
    m = pattern.search(content)
    if not m:
        return None
    start = m.start()
    next_entry = re.search(r'(?m)^- name:', content[m.end():])
    end = m.end() + next_entry.start() if next_entry else len(content)
    return start, end


def main():
    if len(sys.argv) != 3:
        print("usage: resolve_zenodo_doi.py <package-name> <repo>", file=sys.stderr)
        sys.exit(2)

    package_name, repo = sys.argv[1], sys.argv[2]

    with open(ZENODO_YML, encoding="utf-8") as f:
        content = f.read()

    span = find_entry_span(content, package_name)
    if span is None:
        print(f"! no entry named '{package_name}' found in {ZENODO_YML}", file=sys.stderr)
        sys.exit(1)

    start, end = span
    entry_text = content[start:end]

    if re.search(r'zenodo_concept_doi:\s*"10\.5281', entry_text):
        print(f"'{package_name}' already has a concept DOI set -- nothing to do")
        print("CHANGED=0")
        return

    concept_doi = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        concept_doi = search_concept_doi(package_name, repo)
        if concept_doi:
            break
        print(f"  attempt {attempt}/{MAX_ATTEMPTS}: not found yet, waiting {RETRY_DELAY_SECONDS}s")
        time.sleep(RETRY_DELAY_SECONDS)

    if not concept_doi:
        print(f"! could not find a Zenodo record for '{repo}' after {MAX_ATTEMPTS} attempts")
        print("CHANGED=0")
        return

    new_entry_text = re.sub(
        r'zenodo_concept_doi:\s*null',
        f'zenodo_concept_doi: "{concept_doi}"',
        entry_text,
        count=1,
    )
    if new_entry_text == entry_text:
        print(f"! could not find 'zenodo_concept_doi: null' in the '{package_name}' entry to replace")
        print("CHANGED=0")
        return

    new_content = content[:start] + new_entry_text + content[end:]
    with open(ZENODO_YML, "w", encoding="utf-8") as f:
        f.write(new_content)

    print(f"✓ set {package_name}'s concept DOI to {concept_doi}")
    print("CHANGED=1")


if __name__ == "__main__":
    main()
