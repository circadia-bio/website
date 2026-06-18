#!/usr/bin/env python3
"""Pre-render script: generates the people page carousels and injects them
into people.qmd between sentinel comments.

Usage (from repo root):
    python generate_people.py
"""

import os
import re


# ── Helpers ──────────────────────────────────────────────────────────────────────────────

def read_frontmatter(path):
    try:
        with open(path, encoding="utf-8") as f:
            content = f.read()
    except OSError:
        return None
    m = re.match(r"^---\n(.*?)\n---", content, re.DOTALL)
    if not m:
        return None
    fm = {}
    for line in m.group(1).splitlines():
        m2 = re.match(r'^([\w\-]+)\s*:\s*"?(.*?)"?\s*, line)
        if m2:
            fm[m2.group(1)] = m2.group(2).strip('"').strip("'")
    return fm or None


def inject(qmd_path, sentinel, html):
    with open(qmd_path, encoding="utf-8") as f:
        lines = f.readlines()
    start = next(
        (i for i, l in enumerate(lines) if l.strip() == f"<!-- {sentinel}:START -->"),
        None,
    )
    end = next(
        (i for i, l in enumerate(lines) if l.strip() == f"<!-- {sentinel}:END -->"),
        None,
    )
    if start is None or end is None:
        raise RuntimeError(f"Could not find {sentinel} sentinels in {qmd_path}")
    block = "```{=html}\n" + html + "\n```"
    new_lines = (
        lines[: start + 1]
        + [l + "\n" for l in block.splitlines()]
        + lines[end:]
    )
    with open(qmd_path, "w", encoding="utf-8") as f:
        f.writelines(new_lines)


# ── Staff spotlight carousel ─────────────────────────────────────────────────────

staff_qmds = sorted(
    [
        os.path.join("people/staff", f)
        for f in os.listdir("people/staff")
        if f.endswith(".qmd")
    ]
)

staff = []
for qmd in staff_qmds:
    fm = read_frontmatter(qmd)
    if not fm:
        continue
    slug = os.path.splitext(os.path.basename(qmd))[0]
    img_file = fm.get("image")
    staff.append(
        {
            "sortby": int(fm.get("sortby", 99)),
            "href": f"people/staff/{slug}.html",
            "img": f"people/staff/{img_file}" if img_file else None,
            "alt": fm.get("image-alt", f"Photo of {fm.get('title', '')}"),
            "role": fm.get("subtitle", ""),
            "name": fm.get("title", ""),
            "bio": fm.get("description", ""),
        }
    )

staff.sort(key=lambda x: x["sortby"])


def staff_slide(s, active):
    img_html = f'<img src="{s["img"]}" alt="{s["alt"]}">' if s["img"] else ""
    active_class = " rs-active" if active else ""
    return (
        f'    <div class="rs-slide{active_class}">\n'
        f'      <a class="rs-img-col" href="{s["href"]}">\n'
        f"        {img_html}\n"
        f"      </a>\n"
        f'      <div class="rs-body">\n'
        f'        <p class="rs-role">{s["role"]}</p>\n'
        f'        <h3 class="rs-name">{s["name"]}</h3>\n'
        f'        <p class="rs-bio">{s["bio"]}</p>\n'
        f'        <a href="{s["href"]}" class="rs-link">View profile \u2192</a>\n'
        f"      </div>\n"
        f"    </div>"
    )


staff_slides_html = "\n\n".join(staff_slide(s, i == 0) for i, s in enumerate(staff))

staff_html = f"""<div id="rs-wrap">
  <div id="rs-track">

{staff_slides_html}

  </div>
  <div id="rs-controls">
    <button class="rs-btn" id="rs-prev" aria-label="Previous">&#8249;</button>
    <div id="rs-dots"></div>
    <button class="rs-btn" id="rs-next" aria-label="Next">&#8250;</button>
  </div>
</div>

<script>
(function(){{
  var slides = document.querySelectorAll('#rs-track .rs-slide');
  var dotsEl = document.getElementById('rs-dots');
  var N = slides.length;
  var cur = 0;
  var timer;

  slides.forEach(function(_, i){{
    var d = document.createElement('button');
    d.className = 'rs-dot' + (i === 0 ? ' rs-dot-on' : '');
    d.setAttribute('aria-label', 'Slide ' + (i + 1));
    d.onclick = function(){{ goTo(i); }};
    dotsEl.appendChild(d);
  }});

  function goTo(i){{
    slides[cur].classList.remove('rs-active');
    dotsEl.children[cur].classList.remove('rs-dot-on');
    cur = (i + N) % N;
    slides[cur].classList.add('rs-active');
    dotsEl.children[cur].classList.add('rs-dot-on');
  }}

  function startTimer(){{ timer = setInterval(function(){{ goTo(cur + 1); }}, 4000); }}
  function stopTimer(){{ clearInterval(timer); }}

  document.getElementById('rs-prev').onclick = function(){{ stopTimer(); goTo(cur - 1); startTimer(); }};
  document.getElementById('rs-next').onclick = function(){{ stopTimer(); goTo(cur + 1); startTimer(); }};

  var wrap = document.getElementById('rs-wrap');
  wrap.addEventListener('mouseenter', stopTimer);
  wrap.addEventListener('mouseleave', startTimer);

  startTimer();
}})();
</script>"""


# ── Collaborators strip carousel ────────────────────────────────────────────────────

collab_qmds = sorted(
    [
        os.path.join("people/collaborators", f)
        for f in os.listdir("people/collaborators")
        if f.endswith(".qmd")
    ]
)

collabs = []
for qmd in collab_qmds:
    fm = read_frontmatter(qmd)
    if not fm:
        continue
    slug = os.path.splitext(os.path.basename(qmd))[0]
    img_file = fm.get("image")
    collabs.append(
        {
            "sortby": fm.get("sortby", slug),
            "href": f"people/collaborators/{slug}.html",
            "img": f"people/collaborators/{img_file}" if img_file else None,
            "alt": fm.get("image-alt", f"Photo of {fm.get('title', '')}"),
            "name": fm.get("title", ""),
            "inst": fm.get("subtitle", ""),
        }
    )

collabs.sort(key=lambda x: x["sortby"])


def collab_card(c):
    img_html = f'<img src="{c["img"]}" alt="{c["alt"]}">' if c["img"] else ""
    return (
        f'      <a class="co-card" href="{c["href"]}">\n'
        f"        {img_html}\n"
        f'        <p class="co-name">{c["name"]}</p>\n'
        f'        <p class="co-inst">{c["inst"]}</p>\n'
        f"      </a>"
    )


collab_cards_html = "\n\n".join(collab_card(c) for c in collabs)

collabs_html = f"""<div id="co-wrap">
  <div id="co-viewport">
    <div id="co-track">

{collab_cards_html}

    </div>
  </div>
</div>
<div id="co-dots">
  <button class="co-btn" id="co-prev" aria-label="Previous">&#8249;</button>
  <div id="co-dot-list"></div>
  <button class="co-btn" id="co-next" aria-label="Next">&#8250;</button>
</div>

<script>
(function(){{
  var track = document.getElementById('co-track');
  var dotsEl = document.getElementById('co-dot-list');
  var VISIBLE = window.innerWidth < 600 ? 1 : 3;
  var origCards = Array.from(track.querySelectorAll('.co-card'));
  var N = origCards.length;
  var cur = 0;

  origCards.forEach(function(c){{ track.appendChild(c.cloneNode(true)); }});

  for (var i = 0; i < N; i++) {{
    (function(idx) {{
      var d = document.createElement('button');
      d.className = 'co-dot' + (idx === 0 ? ' co-dot-on' : '');
      d.setAttribute('aria-label', 'Go to collaborator ' + (idx + 1));
      d.onclick = function() {{ goTo(idx); }};
      dotsEl.appendChild(d);
    }})(i);
  }}

  function cardW() {{
    var c = track.querySelector('.co-card');
    return c ? c.getBoundingClientRect().width + 14 : 0;
  }}

  function setPos(animated) {{
    track.style.transition = animated ? 'transform 0.4s cubic-bezier(0.4,0,0.2,1)' : 'none';
    track.style.transform = 'translateX(-' + (cur * cardW()) + 'px)';
    dotsEl.querySelectorAll('.co-dot').forEach(function(d, i) {{
      d.classList.toggle('co-dot-on', i === cur % N);
    }});
  }}

  function goTo(i) {{ cur = i; setPos(true); }}
  function next() {{ cur++; setPos(true); }}
  function prev() {{ cur--; setPos(true); }}

  track.addEventListener('transitionend', function() {{
    if (cur >= N) {{ cur = cur - N; setPos(false); }}
    if (cur < 0)  {{ cur = cur + N; setPos(false); }}
  }});

  document.getElementById('co-prev').onclick = prev;
  document.getElementById('co-next').onclick = next;

  var timer = setInterval(next, 3500);
  var wrap = document.getElementById('co-wrap');
  wrap.addEventListener('mouseenter', function() {{ clearInterval(timer); }});
  wrap.addEventListener('mouseleave', function() {{ timer = setInterval(next, 3500); }});

  setPos(false);
}})();
</script>"""


# ── Inject ──────────────────────────────────────────────────────────────────────────────

inject("people.qmd", "STAFF", staff_html)
inject("people.qmd", "COLLABORATORS", collabs_html)

print(f"\u2713 people.qmd updated ({len(staff)} staff, {len(collabs)} collaborators)")
