#!/usr/bin/env python3
"""Pre-render script: generates the news carousel HTML block and injects it
into index.qmd between <!-- NEWS:START --> and <!-- NEWS:END --> sentinels."""

import os
import re
from datetime import datetime

# ── Helpers ────────────────────────────────────────────────────────────────────

def read_frontmatter(path):
    """Parse YAML frontmatter from a .qmd/.md file."""
    try:
        with open(path, encoding="utf-8") as f:
            content = f.read()
    except OSError:
        return None
    m = re.match(r"^---\n(.*?)\n---", content, re.DOTALL)
    if not m:
        return None
    fm = {}
    # Simple key: value parser (handles quoted and unquoted, multiline folded)
    raw = m.group(1)
    # Collapse folded/literal block scalars for description
    # Replace multiline values joined with spaces
    raw = re.sub(r":\s*\|\n((?:  .+\n)+)", lambda m: ": " + " ".join(
        l.strip() for l in m.group(1).splitlines()) + "\n", raw)
    raw = re.sub(r":\s*>\n((?:  .+\n)+)", lambda m: ": " + " ".join(
        l.strip() for l in m.group(1).splitlines()) + "\n", raw)
    raw = re.sub(r":\s*>\-\n((?:  .+\n)+)", lambda m: ": " + " ".join(
        l.strip() for l in m.group(1).splitlines()) + "\n", raw)
    for line in raw.splitlines():
        m2 = re.match(r'^(\w[\w\-]*)\s*:\s*"?(.*?)"?\s*$', line)
        if m2:
            fm[m2.group(1)] = m2.group(2).strip('"').strip("'")
    return fm or None

def find_image(directory, candidates):
    for f in candidates:
        if os.path.exists(os.path.join(directory, f)):
            return f
    return None

def slide_html(href, img_src, img_alt, label, see_all_href,
               title, journal, year, desc, flip, active):
    classes = "nc-slide"
    if active:
        classes += " nc-active"
    if flip:
        classes += " nc-slide-flip"
    journal_line = ""
    if journal and year:
        journal_line = f'          <p class="nc-journal"><em>{journal}</em> &middot; {year}</p>\n'
    return (
        f'      <div class="{classes}">\n'
        f'        <a href="{href}" class="nc-image">'
        f'<img src="{img_src}" alt="{img_alt}"></a>\n'
        f'        <div class="nc-body">\n'
        f'          <div class="nc-label-row">'
        f'<span class="nc-label">{label}</span>'
        f'<a href="{see_all_href}" class="news-see-all">See all</a></div>\n'
        f'          <a href="{href}" class="nc-title">{title}</a>\n'
        f'{journal_line}'
        f'          <p class="nc-desc">{desc}</p>\n'
        f'        </div>\n'
        f'      </div>'
    )

def carousel_html(carousel_id, slides):
    track = "\n\n".join(slides)
    return (
        f'  <div class="nc-carousel" id="{carousel_id}">\n'
        f'    <div class="nc-track">\n\n'
        f'{track}\n\n'
        f'    </div>\n'
        f'    <div class="nc-controls">\n'
        f'      <button class="nc-btn nc-prev" aria-label="Previous"><svg viewBox="0 0 24 24" width="13" height="13" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><polyline points="17 18 11 12 17 6"/></svg></button>\n'
        f'      <div class="nc-dots"></div>\n'
        f'      <button class="nc-btn nc-next" aria-label="Next"><svg viewBox="0 0 24 24" width="13" height="13" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><polyline points="7 18 13 12 7 6"/></svg></button>\n'
        f'    </div>\n'
        f'  </div>'
    )

# ── Tutorials ──────────────────────────────────────────────────────────────────

blog_dirs = sorted([
    d for d in [os.path.join("blog", e) for e in os.listdir("blog")]
    if os.path.isdir(d)
])

tutorials = []
for d in blog_dirs:
    f = next((os.path.join(d, n) for n in ["index.qmd", "index.md"]
               if os.path.exists(os.path.join(d, n))), None)
    if not f:
        continue
    fm = read_frontmatter(f)
    if not fm:
        continue
    img = find_image(d, ["cover.png", "cover.jpg", "cover.svg",
                          "featured.png", "featured.jpg"])
    if not img:
        continue
    try:
        date = datetime.strptime(fm.get("date", "1970-01-01")[:10], "%Y-%m-%d")
    except ValueError:
        continue
    tutorials.append({
        "date":  date,
        "href":  d.lstrip("./") + "/index.html",
        "img":   d.lstrip("./") + "/" + img,
        "title": fm.get("title", ""),
        "desc":  re.sub(r"\s+", " ", fm.get("description", "")).strip(),
    })

tutorials.sort(key=lambda x: x["date"], reverse=True)

tut_slides = [
    slide_html(
        href=t["href"], img_src=t["img"], img_alt=t["title"],
        label="Tutorials", see_all_href="blog.html",
        title=t["title"], journal=None, year=None, desc=t["desc"],
        flip=(i % 2 == 1), active=(i == 0)
    )
    for i, t in enumerate(tutorials)
]

# ── Publications ───────────────────────────────────────────────────────────────

pub_dirs = sorted([
    d for d in [os.path.join("publications", e) for e in os.listdir("publications")]
    if os.path.isdir(d)
])

publications = []
for d in pub_dirs:
    qmd = os.path.join(d, "index.qmd")
    if not os.path.exists(qmd):
        continue
    fm = read_frontmatter(qmd)
    if not fm:
        continue
    img = find_image(d, ["featured.jpg", "featured.png", "featured.svg",
                          "cover.jpg", "cover.png"])
    if not img:
        continue
    try:
        date = datetime.strptime(fm.get("date", "1970-01-01")[:10], "%Y-%m-%d")
    except ValueError:
        continue
    journal = re.sub(r"\*", "", fm.get("subtitle", "")).strip() or None
    publications.append({
        "date":    date,
        "href":    d.lstrip("./") + "/index.html",
        "img":     d.lstrip("./") + "/" + img,
        "title":   fm.get("title", ""),
        "journal": journal,
        "year":    str(date.year),
        "desc":    re.sub(r"\s+", " ", fm.get("description", "")).strip(),
    })

publications.sort(key=lambda x: x["date"], reverse=True)

pub_slides = [
    slide_html(
        href=p["href"], img_src=p["img"], img_alt=p["title"],
        label="Publications", see_all_href="publications.html",
        title=p["title"], journal=p["journal"], year=p["year"], desc=p["desc"],
        flip=(i % 2 == 1), active=(i == 0)
    )
    for i, p in enumerate(publications)
]

# ── JS ─────────────────────────────────────────────────────────────────────────

js = """
<script>
(function(){
  document.querySelectorAll('.nc-carousel').forEach(function(carousel) {
    var slides = carousel.querySelectorAll('.nc-slide');
    var dotsEl = carousel.querySelector('.nc-dots');
    var N = slides.length;
    var current = 0;
    var timer;
    slides.forEach(function(_, i) {
      var d = document.createElement('button');
      d.className = 'nc-dot';
      d.setAttribute('aria-label', 'Slide ' + (i + 1));
      d.onclick = function() { goTo(i); };
      dotsEl.appendChild(d);
    });
    function goTo(idx) {
      slides[current].classList.remove('nc-active');
      dotsEl.children[current].classList.remove('nc-dot-on');
      current = (idx + N) % N;
      slides[current].classList.add('nc-active');
      dotsEl.children[current].classList.add('nc-dot-on');
    }
    carousel.querySelector('.nc-prev').onclick = function() { goTo(current - 1); };
    carousel.querySelector('.nc-next').onclick = function() { goTo(current + 1); };
    function startTimer() { timer = setInterval(function(){ goTo(current + 1); }, 4500); }
    function stopTimer()  { clearInterval(timer); }
    carousel.addEventListener('mouseenter', stopTimer);
    carousel.addEventListener('mouseleave', startTimer);
    goTo(0);
    startTimer();
  });
})();
</script>"""

# ── Assemble ───────────────────────────────────────────────────────────────────

output = "\n".join([
    '<div class="nc-wrap">',
    "",
    carousel_html("nc-tutorials",    tut_slides),
    "",
    carousel_html("nc-publications", pub_slides),
    "",
    "</div>",
    js,
])

# ── Inject into index.qmd ──────────────────────────────────────────────────────

injected = "```{=html}\n" + output + "\n```"

with open("index.qmd", encoding="utf-8") as f:
    lines = f.readlines()

start = next((i for i, l in enumerate(lines) if l.strip() == "<!-- NEWS:START -->"), None)
end   = next((i for i, l in enumerate(lines) if l.strip() == "<!-- NEWS:END -->"),   None)

if start is None or end is None:
    raise RuntimeError("Could not find NEWS:START / NEWS:END sentinels in index.qmd")

new_lines = (
    lines[:start + 1] +
    [l + "\n" for l in injected.splitlines()] +
    lines[end:]
)

with open("index.qmd", "w", encoding="utf-8") as f:
    f.writelines(new_lines)

print(f"✓ index.qmd updated ({len(tut_slides)} tutorials, {len(pub_slides)} publications)")
