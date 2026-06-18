library(yaml)
library(glue)

# ── Helpers ────────────────────────────────────────────────────────────────────

read_frontmatter <- function(path) {
  lines <- readLines(path, warn = FALSE)
  delims <- which(lines == "---")
  if (length(delims) < 2) return(NULL)
  fm <- tryCatch(
    yaml.load(paste(lines[(delims[1] + 1):(delims[2] - 1)], collapse = "\n")),
    error = function(e) NULL
  )
  if (is.null(fm)) return(NULL)
  fm$`_path` <- path
  fm
}

find_image <- function(dir, candidates) {
  for (f in candidates) {
    if (file.exists(file.path(dir, f))) return(f)
  }
  NULL
}

slide_html <- function(href, img_src, img_alt, label, see_all_href,
                       title, journal = NULL, year = NULL, desc, flip, active) {
  classes <- paste(c("nc-slide", if (active) "nc-active", if (flip) "nc-slide-flip"), collapse = " ")
  journal_line <- if (!is.null(journal) && !is.null(year)) {
    glue('          <p class="nc-journal"><em>{journal}</em> &middot; {year}</p>\n')
  } else ""
  glue(
    '      <div class="{classes}">\n',
    '        <a href="{href}" class="nc-image">',
    '<img src="{img_src}" alt="{img_alt}"></a>\n',
    '        <div class="nc-body">\n',
    '          <div class="nc-label-row">',
    '<span class="nc-label">{label}</span>',
    '<a href="{see_all_href}" class="news-see-all">See all</a></div>\n',
    '          <a href="{href}" class="nc-title">{title}</a>\n',
    '{journal_line}',
    '          <p class="nc-desc">{desc}</p>\n',
    '        </div>\n',
    '      </div>'
  )
}

carousel_html <- function(id, slides_html) {
  glue(
    '  <div class="nc-carousel" id="{id}">\n',
    '    <div class="nc-track">\n\n',
    '{paste(slides_html, collapse = "\n\n")}\n\n',
    '    </div>\n',
    '    <div class="nc-controls">\n',
    '      <button class="nc-btn nc-prev" aria-label="Previous">&#8249;</button>\n',
    '      <div class="nc-dots"></div>\n',
    '      <button class="nc-btn nc-next" aria-label="Next">&#8250;</button>\n',
    '    </div>\n',
    '  </div>'
  )
}

# ── Tutorials ──────────────────────────────────────────────────────────────────

blog_dirs <- list.dirs("blog", recursive = FALSE)
tutorials <- Filter(Negate(is.null), lapply(blog_dirs, function(d) {
  qmd <- file.path(d, "index.qmd")
  md  <- file.path(d, "index.md")
  f   <- if (file.exists(qmd)) qmd else if (file.exists(md)) md else return(NULL)
  fm  <- read_frontmatter(f)
  if (is.null(fm)) return(NULL)
  img <- find_image(d, c("cover.png", "cover.jpg", "cover.svg",
                          "featured.png", "featured.jpg"))
  if (is.null(img)) return(NULL)
  list(
    date  = as.Date(as.character(fm$date)),
    href  = paste0(gsub("^\\./", "", d), "/index.html"),
    img   = paste0(gsub("^\\./", "", d), "/", img),
    title = fm$title,
    desc  = gsub("\n", " ", as.character(fm$description))
  )
}))
tutorials <- tutorials[order(sapply(tutorials, `[[`, "date"), decreasing = TRUE)]

tut_slides <- lapply(seq_along(tutorials), function(i) {
  t <- tutorials[[i]]
  slide_html(
    href         = t$href,
    img_src      = t$img,
    img_alt      = t$title,
    label        = "Tutorials",
    see_all_href = "blog.html",
    title        = t$title,
    desc         = t$desc,
    flip         = (i %% 2 == 0),
    active       = (i == 1)
  )
})

# ── Publications ───────────────────────────────────────────────────────────────

pub_dirs <- list.dirs("publications", recursive = FALSE)
publications <- Filter(Negate(is.null), lapply(pub_dirs, function(d) {
  qmd <- file.path(d, "index.qmd")
  if (!file.exists(qmd)) return(NULL)
  fm <- read_frontmatter(qmd)
  if (is.null(fm)) return(NULL)
  img <- find_image(d, c("featured.jpg", "featured.png", "featured.svg",
                          "cover.jpg", "cover.png"))
  if (is.null(img)) return(NULL)
  # Parse journal from subtitle (e.g. "*PLOS ONE*" -> "PLOS ONE")
  journal <- if (!is.null(fm$subtitle)) gsub("\\*", "", fm$subtitle) else NULL
  year    <- format(as.Date(as.character(fm$date)), "%Y")
  list(
    date    = as.Date(as.character(fm$date)),
    href    = paste0(gsub("^\\./", "", d), "/index.html"),
    img     = paste0(gsub("^\\./", "", d), "/", img),
    title   = fm$title,
    journal = journal,
    year    = year,
    desc    = gsub("\n", " ", as.character(fm$description))
  )
}))
publications <- publications[order(sapply(publications, `[[`, "date"), decreasing = TRUE)]

pub_slides <- lapply(seq_along(publications), function(i) {
  p <- publications[[i]]
  slide_html(
    href         = p$href,
    img_src      = p$img,
    img_alt      = p$title,
    label        = "Publications",
    see_all_href = "publications.html",
    title        = p$title,
    journal      = p$journal,
    year         = p$year,
    desc         = p$desc,
    flip         = (i %% 2 == 0),
    active       = (i == 1)
  )
})

# ── Assemble & write ───────────────────────────────────────────────────────────

js <- '
<script>
(function(){
  document.querySelectorAll(\'.nc-carousel\').forEach(function(carousel) {
    var slides = carousel.querySelectorAll(\'.nc-slide\');
    var dotsEl = carousel.querySelector(\'.nc-dots\');
    var N = slides.length;
    var current = 0;
    var timer;
    slides.forEach(function(_, i) {
      var d = document.createElement(\'button\');
      d.className = \'nc-dot\';
      d.setAttribute(\'aria-label\', \'Slide \' + (i + 1));
      d.onclick = function() { goTo(i); };
      dotsEl.appendChild(d);
    });
    function goTo(idx) {
      slides[current].classList.remove(\'nc-active\');
      dotsEl.children[current].classList.remove(\'nc-dot-on\');
      current = (idx + N) % N;
      slides[current].classList.add(\'nc-active\');
      dotsEl.children[current].classList.add(\'nc-dot-on\');
    }
    carousel.querySelector(\'.nc-prev\').onclick = function() { goTo(current - 1); };
    carousel.querySelector(\'.nc-next\').onclick = function() { goTo(current + 1); };
    function startTimer() { timer = setInterval(function(){ goTo(current + 1); }, 4500); }
    function stopTimer()  { clearInterval(timer); }
    carousel.addEventListener(\'mouseenter\', stopTimer);
    carousel.addEventListener(\'mouseleave\', startTimer);
    goTo(0);
    startTimer();
  });
})();
</script>'

output <- paste(
  '```{=html}',
  '<div class="nc-wrap">',
  '',
  carousel_html("nc-tutorials",    tut_slides),
  '',
  carousel_html("nc-publications", pub_slides),
  '',
  '</div>',
  js,
  '```',
  sep = "\n"
)

writeLines(output, "_news-carousel.qmd")
message("✓ _news-carousel.qmd generated (",
        length(tut_slides), " tutorials, ",
        length(pub_slides), " publications)")
