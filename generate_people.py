#!/usr/bin/env python3
"""Pre-render script: generates the people page carousels and injects them
into people.qmd between sentinel comments.

Usage (from repo root):
    python generate_people.py
"""

import os
import re


# -- Helpers -------------------------------------------------------------------

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
        m2 = re.match(r"^([\w\-]+)\s*:\s*[\"']?(.*?)[\"']?\s*$", line)
        if m2:
            fm[m2.group(1)] = m2.group(2).strip('"').strip("'")
    return fm or None


def inject(qmd_path, sentinel, html):
    with open(qmd_path, encoding="utf-8") as f:
        lines = f.readlines()
    start = next(
        (i for i, l in enumerate(lines) if l.strip() == "<!-- " + sentinel + ":START -->"),
        None,
    )
    end = next(
        (i for i, l in enumerate(lines) if l.strip() == "<!-- " + sentinel + ":END -->"),
        None,
    )
    if start is None or end is None:
        raise RuntimeError("Could not find " + sentinel + " sentinels in " + qmd_path)
    block = "```{=html}\n" + html + "\n```"
    new_lines = (
        lines[: start + 1]
        + [l + "\n" for l in block.splitlines()]
        + lines[end:]
    )
    with open(qmd_path, "w", encoding="utf-8") as f:
        f.writelines(new_lines)


# -- Staff spotlight carousel --------------------------------------------------

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
            "href": "people/staff/" + slug + ".html",
            "img": "people/staff/" + img_file if img_file else None,
            "alt": fm.get("image-alt", "Photo of " + fm.get("title", "")),
            "role": fm.get("subtitle", ""),
            "name": fm.get("title", ""),
            "bio": fm.get("description", ""),
        }
    )

staff.sort(key=lambda x: x["sortby"])


def staff_slide(s, active):
    img_html = '<img src="' + s["img"] + '" alt="' + s["alt"] + '">' if s["img"] else ""
    active_class = " rs-active" if active else ""
    lines = [
        '    <div class="rs-slide' + active_class + '">',
        '      <a class="rs-img-col" href="' + s["href"] + '">',
        "        " + img_html,
        "      </a>",
        '      <div class="rs-body">',
        '        <p class="rs-role">' + s["role"] + "</p>",
        '        <h3 class="rs-name">' + s["name"] + "</h3>",
        '        <p class="rs-bio">' + s["bio"] + "</p>",
        '        <a href="' + s["href"] + '" class="rs-link">View profile \u2192</a>',
        "      </div>",
        "    </div>",
    ]
    return "\n".join(lines)


staff_slides_html = "\n\n".join(staff_slide(s, i == 0) for i, s in enumerate(staff))

STAFF_JS = """(function(){
  var slides = document.querySelectorAll('#rs-track .rs-slide');
  var dotsEl = document.getElementById('rs-dots');
  var N = slides.length;
  var cur = 0;
  var timer;

  slides.forEach(function(_, i){
    var d = document.createElement('button');
    d.className = 'rs-dot' + (i === 0 ? ' rs-dot-on' : '');
    d.setAttribute('aria-label', 'Slide ' + (i + 1));
    d.onclick = function(){ goTo(i); };
    dotsEl.appendChild(d);
  });

  function goTo(i){
    slides[cur].classList.remove('rs-active');
    dotsEl.children[cur].classList.remove('rs-dot-on');
    cur = (i + N) % N;
    slides[cur].classList.add('rs-active');
    dotsEl.children[cur].classList.add('rs-dot-on');
  }

  function startTimer(){ timer = setInterval(function(){ goTo(cur + 1); }, 4000); }
  function stopTimer(){ clearInterval(timer); }

  document.getElementById('rs-prev').onclick = function(){ stopTimer(); goTo(cur - 1); startTimer(); };
  document.getElementById('rs-next').onclick = function(){ stopTimer(); goTo(cur + 1); startTimer(); };

  var wrap = document.getElementById('rs-wrap');
  wrap.addEventListener('mouseenter', stopTimer);
  wrap.addEventListener('mouseleave', startTimer);

  startTimer();
})();"""

staff_html = (
    '<div id="rs-wrap">\n'
    '  <div id="rs-track">\n\n'
    + staff_slides_html
    + '\n\n  </div>\n'
    '  <div id="rs-controls">\n'
    '    <button class="rs-btn" id="rs-prev" aria-label="Previous"><svg viewBox="0 0 24 24" width="13" height="13" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><polyline points="17 18 11 12 17 6"/></svg></button>\n'
    '    <div id="rs-dots"></div>\n'
    '    <button class="rs-btn" id="rs-next" aria-label="Next"><svg viewBox="0 0 24 24" width="13" height="13" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><polyline points="7 18 13 12 7 6"/></svg></button>\n'
    '  </div>\n'
    '</div>\n\n'
    '<script>\n' + STAFF_JS + '\n</script>'
)


# -- Collaborators strip carousel ----------------------------------------------

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
            "href": "people/collaborators/" + slug + ".html",
            "img": "people/collaborators/" + img_file if img_file else None,
            "alt": fm.get("image-alt", "Photo of " + fm.get("title", "")),
            "name": fm.get("title", ""),
            "inst": fm.get("subtitle", ""),
        }
    )

collabs.sort(key=lambda x: x["sortby"])


def collab_card(c):
    img_html = '<img src="' + c["img"] + '" alt="' + c["alt"] + '">' if c["img"] else ""
    lines = [
        '      <a class="co-card" href="' + c["href"] + '">',
        "        " + img_html,
        '        <p class="co-name">' + c["name"] + "</p>",
        '        <p class="co-inst">' + c["inst"] + "</p>",
        "      </a>",
    ]
    return "\n".join(lines)


collab_cards_html = "\n\n".join(collab_card(c) for c in collabs)

COLLABS_JS = """(function(){
  var track = document.getElementById('co-track');
  var dotsEl = document.getElementById('co-dot-list');
  var VISIBLE = window.innerWidth < 600 ? 1 : 3;
  var origCards = Array.from(track.querySelectorAll('.co-card'));
  var N = origCards.length;
  var cur = 0;

  origCards.forEach(function(c){ track.appendChild(c.cloneNode(true)); });

  for (var i = 0; i < N; i++) {
    (function(idx) {
      var d = document.createElement('button');
      d.className = 'co-dot' + (idx === 0 ? ' co-dot-on' : '');
      d.setAttribute('aria-label', 'Go to collaborator ' + (idx + 1));
      d.onclick = function() { goTo(idx); };
      dotsEl.appendChild(d);
    })(i);
  }

  function cardW() {
    var c = track.querySelector('.co-card');
    return c ? c.getBoundingClientRect().width + 14 : 0;
  }

  function setPos(animated) {
    track.style.transition = animated ? 'transform 0.4s cubic-bezier(0.4,0,0.2,1)' : 'none';
    track.style.transform = 'translateX(-' + (cur * cardW()) + 'px)';
    dotsEl.querySelectorAll('.co-dot').forEach(function(d, i) {
      d.classList.toggle('co-dot-on', i === cur % N);
    });
  }

  function goTo(i) { cur = i; setPos(true); }
  function next() { cur++; setPos(true); }
  function prev() { cur--; setPos(true); }

  track.addEventListener('transitionend', function() {
    if (cur >= N) { cur = cur - N; setPos(false); }
    if (cur < 0)  { cur = cur + N; setPos(false); }
  });

  document.getElementById('co-prev').onclick = prev;
  document.getElementById('co-next').onclick = next;

  var timer = setInterval(next, 3500);
  var wrap = document.getElementById('co-wrap');
  wrap.addEventListener('mouseenter', function() { clearInterval(timer); });
  wrap.addEventListener('mouseleave', function() { timer = setInterval(next, 3500); });

  setPos(false);
})();"""

collabs_html = (
    '<div id="co-wrap">\n'
    '  <div id="co-viewport">\n'
    '    <div id="co-track">\n\n'
    + collab_cards_html
    + '\n\n    </div>\n'
    '  </div>\n'
    '</div>\n'
    '<div id="co-dots">\n'
    '  <button class="co-btn" id="co-prev" aria-label="Previous"><svg viewBox="0 0 24 24" width="13" height="13" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><polyline points="17 18 11 12 17 6"/></svg></button>\n'
    '  <div id="co-dot-list"></div>\n'
    '  <button class="co-btn" id="co-next" aria-label="Next"><svg viewBox="0 0 24 24" width="13" height="13" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><polyline points="7 18 13 12 7 6"/></svg></button>\n'
    '</div>\n\n'
    '<script>\n' + COLLABS_JS + '\n</script>'
)


# -- Inject --------------------------------------------------------------------

inject("people.qmd", "STAFF", staff_html)
inject("people.qmd", "COLLABORATORS", collabs_html)

print("\u2713 people.qmd updated (" + str(len(staff)) + " staff, " + str(len(collabs)) + " collaborators)")
