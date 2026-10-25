# =============================================================================
# scripts/run_all.R
# One command reproduction of every result in the manuscript.
#   Rscript scripts/run_all.R            full run (stages cached in results/parts)
#   MHW_QUICK=1 Rscript scripts/run_all.R   small sizes, for checking the pipeline
#   MHW_CORES=n                          forked workers for replicates (default: all - 1)
# Each stage is skipped if its RDS already exists; delete results/parts/<stage>.rds
# to rerun one stage.
# =============================================================================
t_start <- Sys.time()
for (f in list.files("R", full.names = TRUE)) source(f)
dir.create("results/parts", showWarnings = FALSE, recursive = TRUE)

cores <- as.integer(Sys.getenv("MHW_CORES", max(1L, parallel::detectCores() - 1L)))
options(mhw.cores = cores)
message("using ", cores, " core(s)")

cfg <- default_config()
saveRDS(cfg, "results/config.rds")

SZ <- list(fac_R = 400, cov_R = 300, mat_R = 80, tmp_R = 300,
           sens_design = 200, sens_R = 20, frame_R = 250)
if (nzchar(Sys.getenv("MHW_QUICK")))
  SZ <- list(fac_R = 12, cov_R = 10, mat_R = 6, tmp_R = 8,
             sens_design = 10, sens_R = 3, frame_R = 8)
saveRDS(SZ, "results/run_sizes.rds")

## cache key: run sizes, configuration and the R source; a cached stage whose
## key differs is ignored and recomputed, so a quick run never reuses full
## stages and a code change never leaves stale results in place
cache_key <- function() {
  src <- unlist(lapply(sort(list.files("R", full.names = TRUE)), readLines, warn = FALSE))
  tf <- tempfile(); writeLines(c(deparse(SZ), deparse(cfg), src, R.version.string), tf)
  on.exit(unlink(tf)); substr(unname(tools::md5sum(tf)), 1, 12)
}
KEY <- cache_key()
writeLines(KEY, "results/parts/cache_key.txt")
stage <- function(nm, expr) {
  path <- file.path("results/parts", paste0(nm, ".rds"))
  if (file.exists(path)) {
    v <- readRDS(path)
    if (identical(attr(v, "cache_key"), KEY)) { message("=== ", nm, " (cached)"); return(v) }
    message("=== ", nm, " (cached copy has a different key; recomputing)")
  }
  message("=== ", nm, " === ", format(Sys.time(), "%H:%M:%S"))
  v <- force(expr)
  attr(v, "cache_key") <- KEY
  saveRDS(v, path)
  message("    done at ", round(as.numeric(Sys.time() - t_start, units = "mins"), 1), " min")
  v
}

ident <- stage("identifiability", identifiability_experiment(cfg))
fac   <- stage("factorial", exp_factorial(cfg, R = SZ$fac_R))
frame <- stage("frame", exp_frame(cfg, R = SZ$frame_R))
cov   <- stage("coverage", exp_coverage(cfg, R = SZ$cov_R, cv_grid = c(cfg$cv_surv, 0.90)))
mat   <- stage("matched", exp_matched(cfg, R = SZ$mat_R))
tmp   <- stage("temporal", exp_temporal(cfg, R = SZ$tmp_R))
sens  <- stage("sensitivity", exp_sensitivity(cfg, n_design = SZ$sens_design, R = SZ$sens_R))

## assemble machine readable outputs -------------------------------------------
saveRDS(ident, "results/identifiability.rds")
saveRDS(fac, "results/factorial.rds")
write.csv(fac, "results/factorial_replicates.csv", row.names = FALSE)
write.csv(factorial_summary(fac), "results/factorial_summary.csv", row.names = FALSE)
anova_tab <- do.call(rbind, lapply(c("beta_nominal", "beta_delta", "beta_survey"),
                                   function(r) factorial_anova(fac, r)))
write.csv(anova_tab, "results/factorial_anova.csv", row.names = FALSE)
saveRDS(exp_classifier(fac), "results/classifier.rds")
saveRDS(frame, "results/frame.rds")
saveRDS(cov, "results/coverage.rds")
write.csv(coverage_summary(cov), "results/coverage_summary.csv", row.names = FALSE)
write.csv(mat, "results/matched.csv", row.names = FALSE)
saveRDS(tmp, "results/temporal.rds")
saveRDS(sens, "results/sensitivity.rds")
write.csv(sens, "results/sensitivity_design.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
message("all stages complete in ", round(as.numeric(Sys.time() - t_start, units = "mins"), 1), " min")
