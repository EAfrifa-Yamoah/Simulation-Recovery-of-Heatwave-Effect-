# =============================================================================
# scripts/build_manuscript.R
# Build manuscript.docx from the single source manuscript.md with pandoc.
#   Rscript scripts/build_manuscript.R
#
# The source may contain:
#   {{TABLE:name}}   replaced by results/tables/<name>.md
#   {{r expr}}       replaced by the value of the R expression, evaluated with
#                    the results loaded (fs = factorial summary, cs = coverage
#                    summary, fr = frame, tp = temporal, cl = classifier, id =
#                    identifiability, sz = run sizes, mm = matched, an = anova),
#                    so that every number quoted in the text is reproducible.
# LaTeX math ($...$ and $$...$$) is converted by pandoc into native Word
# equations. Figures are referenced by their PNG in figures/.
# =============================================================================
for (f in list.files("R", full.names = TRUE)) source(f)
stopifnot(file.exists("manuscript.md"))
dir.create("config", showWarnings = FALSE)

env <- new.env()
for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = env)
env$fs <- read.csv("results/factorial_summary.csv")
env$cs <- read.csv("results/coverage_summary.csv")
env$fr <- readRDS("results/frame.rds")
env$tp <- readRDS("results/temporal.rds")
env$cl <- readRDS("results/classifier.rds")
env$id <- readRDS("results/identifiability.rds")
env$sz <- readRDS("results/run_sizes.rds")
env$mm <- read.csv("results/matched.csv")
env$an <- read.csv("results/factorial_anova.csv")
env$fac <- readRDS("results/factorial.rds")
env$cfg <- readRDS("results/config.rds")
env$f2 <- function(x) formatC(x, digits = 2, format = "f")
env$f3 <- function(x) formatC(x, digits = 3, format = "f")
env$sc <- function(s, col) env$fs[[col]][env$fs$scenario == s]
env$pct <- function(x) paste0(round(100 * x), " percent")
env$ane <- function(resp, term) { a <- env$an[env$an$response == resp & env$an$term == term, ]; a$estimate[1] }
env$anse <- function(resp, term) { a <- env$an[env$an$response == resp & env$an$term == term, ]; a$se[1] }
env$cov <- function(every, design, coverage, col, cv = env$cfg$cv_surv) {
  d <- env$cs[env$cs$every == every & env$cs$design == design & abs(env$cs$coverage - coverage) < 1e-9 & abs(env$cs$cv_surv - cv) < 1e-9, ]; d[[col]][1] }
env$frm <- function(frame, scen, col) { d <- env$fr[env$fr$frame == frame & env$fr$scenario == scen, ]; mean(d[[col]], na.rm = TRUE) }
env$tpo <- function(col) mean(env$tp$summary_obs[[col]]); env$tpt <- function(col) mean(env$tp$summary_tru[[col]])
env$sens <- read.csv("results/sensitivity_design.csv")
env$srcv <- function(resp, par) { z <- src(env$sens, resp); z$src[z$parameter == par] }
env$srcr2 <- function(resp) src(env$sens, resp)$r2[1]
env$mcse_max <- function() signif(max(env$fs[, grep("^mcse_beta", names(env$fs))], na.rm = TRUE), 1)
env$cells3 <- function(scen, col) { d <- env$fac[env$fac$scenario == scen, ]; mean(d[[col]], na.rm = TRUE) }

txt <- paste(readLines("manuscript.md", warn = FALSE), collapse = "\n")

# tables
tab_tok <- unique(regmatches(txt, gregexpr("\\{\\{TABLE:[A-Za-z0-9_]+\\}\\}", txt))[[1]])
for (tk in tab_tok) {
  nm <- sub("\\{\\{TABLE:(.*)\\}\\}", "\\1", tk)
  path <- file.path("results/tables", paste0(nm, ".md"))
  if (!file.exists(path)) stop("missing table: ", path)
  txt <- gsub(tk, paste(readLines(path, warn = FALSE), collapse = "\n"), txt, fixed = TRUE)
}
# inline R
toks <- regmatches(txt, gregexpr("\\{\\{r [^}]+\\}\\}", txt))[[1]]
for (tk in unique(toks)) {
  expr <- sub("^\\{\\{r (.*)\\}\\}$", "\\1", tk)
  val <- tryCatch(eval(parse(text = expr), envir = env), error = function(e) stop("in {{r ", expr, "}}: ", conditionMessage(e)))
  if (length(val) != 1 || is.na(val)) stop("token gave no value: {{r ", expr, "}} -> ", paste(val, collapse = ","))
  txt <- gsub(tk, as.character(val), txt, fixed = TRUE)
}
left <- regmatches(txt, gregexpr("\\{\\{[^}]*\\}\\}", txt))[[1]]
if (length(left)) stop("unresolved tokens: ", paste(left, collapse = " "))
writeLines(txt, "results/manuscript_build.md")

ref <- "config/reference.docx"
if (!file.exists(ref)) system2("pandoc", c("-o", ref, "--print-default-data-file", "reference.docx"))
out <- "manuscript_sim-rec_revised.docx"
st <- system2("pandoc", c("results/manuscript_build.md", "-f", "markdown+tex_math_dollars",
                          "-t", "docx", "--reference-doc", ref, "-o", out))
if (st != 0) stop("pandoc failed")

# verify: a valid zip, well formed XML, equations present, no unresolved markers
chk <- system2("unzip", c("-p", out, "word/document.xml"), stdout = TRUE)
xml <- paste(chk, collapse = "")
stopifnot(grepl("<m:oMath", xml), !grepl("\\{\\{", xml))
n_eq <- length(gregexpr("<m:oMath>", xml)[[1]])
cat("wrote", out, "with", n_eq, "equations\n")
