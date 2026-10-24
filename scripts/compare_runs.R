# =============================================================================
# scripts/compare_runs.R
# Compare the quantities quoted in the manuscript between a saved results
# directory (before the fix of review finding 1) and results/ (after).
#   Rscript scripts/compare_runs.R <path to old results dir>
# Writes results/review_checks/moved_numbers.csv and a LaTeX table body.
# =============================================================================
Sys.setlocale("LC_CTYPE", "C.UTF-8")
args <- commandArgs(trailingOnly = TRUE)
old_dir <- if (length(args)) args[1] else stop("give the old results directory")
for (f in list.files("R", full.names = TRUE)) source(f)
dir.create("results/review_checks", showWarnings = FALSE)

load_set <- function(dir) {
  fs <- read.csv(file.path(dir, "factorial_summary.csv")); fac <- readRDS(file.path(dir, "factorial.rds"))
  cs <- read.csv(file.path(dir, "coverage_summary.csv")); mm <- read.csv(file.path(dir, "matched.csv"))
  tp <- readRDS(file.path(dir, "temporal.rds")); cl <- readRDS(file.path(dir, "classifier.rds"))
  fr <- readRDS(file.path(dir, "frame.rds")); sens <- read.csv(file.path(dir, "sensitivity_design.csv"))
  cfg <- readRDS(file.path(dir, "config.rds")); KM <- cfg$cell_km * sd(make_grid(cfg)$lat)
  sc <- function(s, col) fs[[col]][fs$scenario == s]
  mg <- function(k, d, col) mm[[col]][abs(mm$kappa - k) < 1e-9 & abs(mm$delta_dep - d) < 1e-9]
  cov <- function(every, coverage, col, cv = 0.45) cs[[col]][cs$every == every & cs$design == "srs" &
    abs(cs$coverage - coverage) < 1e-9 & abs(cs$cv_surv - cv) < 1e-9][1]
  cg <- function(s, col) KM * mean(fac[[col]][fac$scenario == s], na.rm = TRUE)
  cm <- as.matrix(cl$confusion)
  c(`True effect, decline (S1001)` = sc("S1001", "truth_total"),
    `Nominal, decline (S1001)` = sc("S1001", "beta_nominal"),
    `Delta, decline (S1001)` = sc("S1001", "beta_delta"),
    `Survey, decline (S1001)` = sc("S1001", "beta_survey"),
    `Nominal, catchability loss (S0011)` = sc("S0011", "beta_nominal"),
    `Delta, catchability loss (S0011)` = sc("S0011", "beta_delta"),
    `Survey, catchability loss (S0011)` = sc("S0011", "beta_survey"),
    `True effect, decline + catchability (S1011)` = sc("S1011", "truth_total"),
    `Nominal, decline + catchability (S1011)` = sc("S1011", "beta_nominal"),
    `Survey, decline + catchability (S1011)` = sc("S1011", "beta_survey"),
    `In-frame effect, redistribution (S0101)` = sc("S0101", "truth_in"),
    `Survey, redistribution (S0101)` = sc("S0101", "beta_survey"),
    `Spatial delta, redistribution (S0101)` = sc("S0101", "beta_deltaint"),
    `Nominal, warm-edge loss (SPQ1)` = sc("SPQ1", "beta_nominal"),
    `Matched grid: nominal, loss only` = mg(0, 0.07, "beta_nominal"),
    `Matched grid: nominal, mixture` = mg(0.06, 0.035, "beta_nominal"),
    `Matched grid: nominal, decline only` = mg(0.12, 0, "beta_nominal"),
    `Matched grid: true effect, decline only` = mg(0.12, 0, "beta_truth"),
    `Fishery CoG shift, redistribution (km)` = cg("S0101", "dG_fishery"),
    `Survey CoG shift, redistribution (km)` = cg("S0101", "dG_survey"),
    `Fishery CoG shift, warm-edge loss (km)` = cg("SPQ1", "dG_fishery"),
    `Survey RMSE, 7 stations monthly` = cov(1, 0.05, "rmse"),
    `Survey RMSE, 14 stations monthly` = cov(1, 0.10, "rmse"),
    `Survey coverage, 14 stations monthly` = cov(1, 0.10, "cover_survey"),
    `Survey power, 14 stations monthly` = cov(1, 0.10, "power_survey"),
    `Nominal RMSE, design scenario` = cov(1, 0.10, "rmse_nominal"),
    `Delta RMSE, design scenario` = cov(1, 0.10, "rmse_delta"),
    `Survey RMSE, 29 stations quarterly` = cov(3, 0.20, "rmse"),
    `Distributed lag R2 (survey target)` = mean(tp$summary_obs$r2_dl),
    `Distributed lag corr. with truth` = mean(tp$summary_obs$cor_obs_truth_dl),
    `Distributed lag OOS RMSE` = mean(tp$summary_obs$oos_dl),
    `Same-month OOS RMSE` = mean(tp$summary_obs$oos_acute),
    `Key accuracy` = cl$accuracy, `Key indeterminate` = cl$indeterminate,
    `Frame: in-frame effect, warming refuge` = mean(fr$beta_cf_in[fr$frame == "open_dynamic" & fr$scenario == "S0101"]),
    `Frame: survey, static refuge` = mean(fr$beta_survey[fr$frame == "open_static" & fr$scenario == "S0101"]),
    `SRC nominal bias on catchability loss` = src(sens, "gap_nominal")$src[src(sens, "gap_nominal")$parameter == "delta_dep"],
    `Max MCSE` = max(fs[, grep("^mcse_beta", names(fs))], na.rm = TRUE))
}
a <- load_set(old_dir); b <- load_set("results")
d <- data.frame(quantity = names(a), before = round(a, 4), after = round(b, 4), change = round(b - a, 4), row.names = NULL)
write.csv(d, "results/review_checks/moved_numbers.csv", row.names = FALSE)
fmt <- function(x) { v <- formatC(x, digits = 3, format = "f"); ifelse(substr(v, 1, 1) == "-", paste0("$-$", substring(v, 2)), v) }
body <- c("\\begin{tabular}{@{}lrrr@{}}", "\\toprule", "Quantity & Before & After & Change \\\\", "\\midrule",
          sprintf("%s & %s & %s & %s \\\\", d$quantity, fmt(d$before), fmt(d$after), fmt(d$change)),
          "\\bottomrule", "\\end{tabular}")
writeLines(body, "results/review_checks/moved_numbers.tex")
print(d, digits = 4)
cat("largest absolute change in a heatwave coefficient or target:",
    max(abs(d$change[!grepl("km|RMSE|R2|corr|accuracy|indeterminate|SRC|MCSE|coverage|power", d$quantity)])), "\n")
