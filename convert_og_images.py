#!/usr/bin/env python3
"""
Convert project SVG images to PNG for OG/social media compatibility.
WhatsApp and other platforms don't support SVG in og:image.

Run from the website root:
    pip install cairosvg Pillow
    python3 convert_og_images.py

Outputs 1200x630 PNG (standard OG image size) with the SVG logo
centred on a white background.
"""

import cairosvg
import urllib.request
from pathlib import Path
from PIL import Image
import io

IMAGES_DIR = Path(__file__).parent / "projects" / "images"
OG_WIDTH = 1200
OG_HEIGHT = 630
PADDING = 80  # px padding around the logo

# Local SVGs to convert (filename without extension, in projects/images/)
LOCAL_SVGS = [
    "slumbr_logo",
    "tallier_logo",
    "circadiabase_docker",
    "actigraphy_clustering",
    "actt_validation",
    "cycling_entropy",
    "pinetime_research",
]

# External SVGs: (output_stem, url)
EXTERNAL_SVGS = [
    ("r-itable_og", "https://r-itable.circadia-lab.uk/logo.svg"),
]


def svg_to_og_png(svg_bytes: bytes, output_path: Path):
    """Convert SVG bytes to PNG at natural resolution with transparency."""
    png_bytes = cairosvg.svg2png(bytestring=svg_bytes, output_width=600, output_height=600)

    img = Image.open(io.BytesIO(png_bytes)).convert("RGBA")
    img.save(output_path, "PNG", optimize=True)
    print(f"  ✅ {output_path.name} ({img.width}×{img.height})")


def main():
    # Local SVGs
    for stem in LOCAL_SVGS:
        svg_path = IMAGES_DIR / f"{stem}.svg"
        png_path = IMAGES_DIR / f"{stem}_og.png"

        if not svg_path.exists():
            print(f"  ⚠️  Not found: {svg_path.name} — skipping")
            continue

        print(f"Converting {svg_path.name}...")
        svg_to_og_png(svg_path.read_bytes(), png_path)

    # External SVGs
    for stem, url in EXTERNAL_SVGS:
        png_path = IMAGES_DIR / f"{stem}.png"
        print(f"Fetching {url}...")
        try:
            with urllib.request.urlopen(url, timeout=10) as r:
                svg_bytes = r.read()
            svg_to_og_png(svg_bytes, png_path)
        except Exception as e:
            print(f"  ❌ Failed to fetch {url}: {e}")


if __name__ == "__main__":
    main()
