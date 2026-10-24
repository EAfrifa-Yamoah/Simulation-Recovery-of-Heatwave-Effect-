# =============================================================================
# scripts/build_latex.R
# Build a LaTeX source (and, if a TeX distribution is present, the PDF) from the
# same resolved manuscript that build_manuscript.R writes to
# results/manuscript_build.md. Run build_manuscript.R first (make doc).
#   Rscript scripts/build_latex.R
#
# What this does beyond a plain pandoc conversion
#   * title, author and abstract go into the document preamble
#   * display equations become numbered amsmath equation environments
#   * text that follows an equation is not indented (\noindent)
#   * the manuscript's own figure and table numbers are kept (caption labels
#     are suppressed, so "Table 5" stays Table 5 wherever it appears)
#   * tables are set in a small font; tables with more than nine columns are
#     rotated to landscape
#   * characters outside the pdflatex input encoding are mapped to LaTeX
# =============================================================================
src_md <- "results/manuscript_build.md"
if (!file.exists(src_md)) stop("run scripts/build_manuscript.R first (it writes ", src_md, ")")
dir.create("results/latex", showWarnings = FALSE, recursive = TRUE)

x <- readLines(src_md, warn = FALSE, encoding = "UTF-8")
Sys.setlocale("LC_CTYPE", "C.UTF-8")

## 1. title, abstract and keywords -> metadata -------------------------------
title <- sub("^# ", "", x[1]); x <- x[-1]
ia <- which(x == "## Abstract"); ib <- which(grepl("^## ", x) & seq_along(x) > ia)[1]
abs_block <- x[(ia + 1):(ib - 1)]; abs_block <- abs_block[nzchar(abs_block)]
x <- x[-(ia:(ib - 1))]
author <- "Eben Afrifa-Yamoah"
affil  <- "School of Science, Edith Cowan University, Joondalup, Western Australia, Australia. e.afrifayamoah@ecu.edu.au"

## 2. section headings: drop the manual numbers, LaTeX numbers them ----------
x <- sub("^(#{2,3}) [0-9]+(\\.[0-9]+)*\\.? ", "\\1 ", x)
unnum <- c("References", "Data and code availability")
for (u in unnum) x <- sub(paste0("^## ", u, "$"), paste0("## ", u, " {-}"), x)

## 3. display equations -> numbered equation environments --------------------
# (all display math in the source is one $$...$$ paragraph on a single line)
is_eq <- grepl("^\\$\\$.*\\$\\$$", x)
x[is_eq] <- paste0("\\begin{equation}\n", sub("^\\$\\$(.*)\\$\\$$", "\\1", x[is_eq]), "\n\\end{equation}")

## 4. table captions: the paragraph after a pipe table becomes its caption ---
# The manuscript numbers tables itself (Table 5 is discussed after Table 8), so
# the label text is kept in the caption and LaTeX's own label is suppressed.
i <- 1
while (i <= length(x)) {
  if (grepl("^Table [0-9]+\\. ", x[i]) && i > 2 && grepl("^\\|", x[i - 2])) {
    x[i] <- paste0(": ", x[i])
  }
  i <- i + 1
}
# figure captions already begin with "Figure n." and are kept as they are;
# reference the PDF version of each figure
x <- gsub("\\(figures/([A-Za-z0-9_]+)\\.png\\)", "(figures/\\1.pdf)", x)

## 5. characters pdflatex cannot take from utf8 ------------------------------
txt <- paste(x, collapse = "\n")
txt <- gsub("2\u2074", "$2^{4}$", txt)            # 2⁴
txt <- gsub("\u2074", "$^{4}$", txt)
txt <- gsub("\u209c", "$_{t}$", txt)              # subscript t
abs_txt <- paste(abs_block, collapse = "\n\n")
abs_txt <- gsub("2\u2074", "$2^{4}$", abs_txt)

## 6. YAML header and pandoc ----------------------------------------------------
yaml_q <- function(s) paste0("'", gsub("'", "''", s), "'")
hdr <- c("---",
         paste0("title: ", yaml_q(title)),
         paste0("author: ", yaml_q(paste0(author, "\\thanks{", affil, "}"))),
         "abstract: |",
         paste0("  ", strsplit(abs_txt, "\n")[[1]]),
         "---", "")
writeLines(enc2utf8(c(hdr, "\\setstretch{1.15}", "", txt)), "results/latex/manuscript_pre.md", useBytes = TRUE)

header_tex <- c(
  "\\usepackage{amsmath,amssymb}",
  "\\usepackage{pdflscape}",
  "\\usepackage[section]{placeins}",
  "\\usepackage{setspace}",
  "\\usepackage{caption}\\captionsetup{labelformat=empty,justification=raggedright,singlelinecheck=false}",
  "\\usepackage{booktabs}",
  "\\providecommand{\\tightlist}{\\setlength{\\itemsep}{0pt}\\setlength{\\parskip}{0pt}}")
writeLines(header_tex, "results/latex/header.tex")

out_tex <- "results/latex/manuscript_raw.tex"
st <- system2("pandoc", c("results/latex/manuscript_pre.md", "-s",
                          "-f", "markdown+tex_math_dollars+raw_tex", "-t", "latex",
                          "--number-sections", "--shift-heading-level-by=-1", "-V", "secnumdepth=2",
                          "-V", "documentclass=article", "-V", "fontsize=11pt",
                          "-V", "papersize=a4", "-V", "geometry:margin=25mm",
                          "-V", "colorlinks=true", "-V", "linkcolor=black",
                          "-V", "urlcolor=blue", "-V", "citecolor=black",
                          "-H", "results/latex/header.tex", "-o", out_tex))
if (st != 0) stop("pandoc failed")

## 7. post processing of the LaTeX -------------------------------------------
tex <- readLines(out_tex, warn = FALSE, encoding = "UTF-8")

# (a) no indentation of the text that follows an equation
starts_cmd <- function(s) grepl("^\\\\(section|subsection|subsubsection|paragraph|begin|end|item|hypertarget|label|noindent)", s)
j <- 1
while (j <= length(tex)) {
  if (tex[j] == "\\end{equation}") {
    k <- j + 1
    while (k <= length(tex) && !nzchar(trimws(tex[k]))) k <- k + 1
    if (k <= length(tex) && !starts_cmd(tex[k])) tex[k] <- paste0("\\noindent ", tex[k])
  }
  j <- j + 1
}

# (b) small font for every table, landscape for wide ones, and column widths
#     that give text columns twice the width of numeric columns
out <- character(0); k <- 1
while (k <= length(tex)) {
  ln <- tex[k]
  if (grepl("^\\\\begin\\{longtable\\}", ln)) {
    m <- k
    while (!grepl("@\\{\\}\\}\\s*$", tex[m])) m <- m + 1
    spec <- tex[(k + 1):m]
    is_text <- grepl("raggedright", spec)
    w <- ifelse(is_text, 2, 1); w <- w / sum(w)
    spec <- mapply(function(s, wi) sub("\\\\real\\{[0-9.]+\\}", sprintf("\\\\real{%.4f}", wi), s), spec, w)
    wide <- length(spec) > 9
    if (wide) out <- c(out, "\\begin{landscape}")
    out <- c(out, "\\FloatBarrier", "\\begingroup", if (wide || length(spec) >= 8) "\\scriptsize" else "\\footnotesize",
             "\\setlength{\\tabcolsep}{3pt}", "\\setlength{\\LTcapwidth}{\\linewidth}", ln, unname(spec))
    k <- m + 1
    # copy the body to \end{longtable}
    while (!grepl("^\\\\end\\{longtable\\}", tex[k])) { out <- c(out, tex[k]); k <- k + 1 }
    out <- c(out, tex[k], "\\endgroup", if (wide) "\\end{landscape}")
  } else out <- c(out, ln)
  k <- k + 1
}
tex <- out

# (c) hanging indent for the reference list
r <- which(grepl("^\\\\section\\*?\\{References\\}", tex))
if (length(r)) {
  e <- which(grepl("^\\\\end\\{document\\}", tex))
  body <- tex[(r + 1):(e - 1)]
  # the first paragraph after a heading is set without \parindent, so its
  # first line would not hang: pull it back by the same amount
  f <- which(nzchar(trimws(body)) & !grepl("^\\\\addcontentsline", body))[1]
  body[f] <- paste0("\\hspace*{-1.5em}", body[f])
  tex <- c(tex[1:r], "\\begin{hangparas}", body, "\\end{hangparas}", tex[e:length(tex)])
  tex <- sub("^\\\\usepackage\\{pdflscape\\}$",
             "\\\\usepackage{pdflscape}\n\\\\newenvironment{hangparas}{\\\\setlength{\\\\parindent}{-1.5em}\\\\setlength{\\\\leftskip}{1.5em}\\\\setlength{\\\\parskip}{4pt}}{\\\\par}",
             tex)
}

final_tex <- "manuscript_sim-rec_revised.tex"
writeLines(enc2utf8(tex), final_tex, useBytes = TRUE)
cat("wrote", final_tex, "\n")

## 8. compile when a TeX distribution is available ------------------------------
if (nzchar(Sys.which("latexmk"))) {
  st <- system2("latexmk", c("-pdf", "-interaction=nonstopmode", "-halt-on-error",
                             "-outdir=results/latex/build", final_tex),
                stdout = "results/latex/latexmk.log", stderr = "results/latex/latexmk.log")
  if (st == 0) {
    file.copy("results/latex/build/manuscript_sim-rec_revised.pdf", "manuscript_sim-rec_revised.pdf", overwrite = TRUE)
    cat("compiled manuscript_sim-rec_revised.pdf\n")
  } else cat("latexmk failed; see results/latex/latexmk.log\n")
} else cat("latexmk not found; compile manuscript_sim-rec_revised.tex with pdflatex (twice)\n")
