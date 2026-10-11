# =============================================================================
# scripts/make_tables.R
# Build every manuscript table from results/ as markdown and CSV.
#   Rscript scripts/make_tables.R
# =============================================================================
for (f in list.files("R", full.names = TRUE)) source(f)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

md_table <- function(df, digits = 3, align = NULL) {
  fmt <- function(x) {
    if (!is.numeric(x)) return(as.character(x))
    if (all(is.na(x) | abs(x - round(x)) < 1e-9)) return(formatC(x, digits = 0, format = "f"))
    formatC(x, digits = digits, format = "f")
  }
  cells <- as.data.frame(lapply(df, fmt), stringsAsFactors = FALSE)
  cells[is.na(df)] <- "—"
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  al  <- if (is.null(align)) vapply(df, function(x) if (is.numeric(x)) "---:" else ":---", "") else align
  sep <- paste0("| ", paste(al, collapse = " | "), " |")
  body <- apply(cells, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, body), collapse = "\n")
}
save_tab <- function(nm, df, digits = 3) {
  write.csv(df, file.path("results/tables", paste0(nm, ".csv")), row.names = FALSE)
  writeLines(md_table(df, digits), file.path("results/tables", paste0(nm, ".md")))
}

label_scenario <- function(sc) {
  if (grepl("^S[01]{4}$", sc)) {
    b <- strsplit(substr(sc, 2, 5), "")[[1]] == "1"
    parts <- c("Abundance", "Redistribution", "Catchability")[b[1:3]]
    base <- if (!length(parts)) "None" else paste(parts, collapse = " + ")
    if (b[4]) paste0(base, " (fleet contraction)") else base
  } else switch(sc, SPQ0 = "Warm edge catchability", SPQ1 = "Warm edge catchability (fleet contraction)",
                SPQT0 = "Warm edge catch. + redistribution",
                SPQT1 = "Warm edge catch. + redistribution (fleet contraction)", sc)
}

fs  <- read.csv("results/factorial_summary.csv")
fs  <- fs[match(c(scenario_levels, scenario_warm), fs$scenario), ]
fs$Scenario <- vapply(fs$scenario, label_scenario, "")

## Table 1: heatwave coefficients against both estimands ----------------------
t1 <- data.frame(Scenario = fs$Scenario,
                 `Truth (total)` = fs$truth_total, `Truth (in frame)` = fs$truth_in,
                 Nominal = fs$beta_nominal, `Positive only` = fs$beta_naive,
                 Delta = fs$beta_delta, `Delta, fished cells` = fs$beta_delta_fished,
                 `Delta + interaction` = fs$beta_deltaint, Survey = fs$beta_survey,
                 check.names = FALSE)
save_tab("table1_coefficients", t1)
cat("Table 1 max MCSE of any estimator mean:",
    signif(max(fs[, grep("^mcse_beta", names(fs))], na.rm = TRUE), 2), "\n")

## Table 2: 2^4 ANOVA, main effects and two way interactions -----------------
an <- read.csv("results/factorial_anova.csv")
an <- an[!grepl(":.*:", an$term) & an$term != "(Intercept)", ]
an$term <- gsub("A", "Abundance", an$term); an$term <- gsub("B", "Redistribution", an$term)
an$term <- gsub("C", "Catchability", an$term); an$term <- gsub("D", "Contraction", an$term)
an$term <- gsub(":", " × ", an$term)
w <- reshape(an, idvar = "term", timevar = "response", direction = "wide")
t2 <- data.frame(Term = w$term,
                 `Nominal` = w$estimate.beta_nominal, `(SE)` = w$se.beta_nominal,
                 `Delta` = w$estimate.beta_delta, `(SE) ` = w$se.beta_delta,
                 `Survey` = w$estimate.beta_survey, `(SE)  ` = w$se.beta_survey,
                 check.names = FALSE)
save_tab("table2_anova", t2, digits = 4)

## Table 3: centroid diagnostics ---------------------------------------------
fac <- readRDS("results/factorial.rds")
cells3 <- c("S0001", "S1001", "S0101", "S0011", "SPQ1", "SPQT1", "S1111")
t3 <- do.call(rbind, lapply(cells3, function(sc) {
  d <- fac[fac$scenario == sc, ]
  data.frame(Scenario = label_scenario(sc),
             `Latent truth` = mean(d$dG_true, na.rm = TRUE),
             `Fishery catch centroid` = mean(d$dG_fishery, na.rm = TRUE),
             `Delta surface centroid` = mean(d$dG_grid, na.rm = TRUE),
             `Survey centroid` = mean(d$dG_survey, na.rm = TRUE),
             `Frame share change` = mean(d$dp_out, na.rm = TRUE),
             `Footprint change` = mean(d$dfoot, na.rm = TRUE), check.names = FALSE)
}))
save_tab("table3_centroids", t3)

## Table 4: survey design curve ----------------------------------------------
cs <- read.csv("results/coverage_summary.csv")
cs <- cs[order(cs$cv_surv, cs$design, cs$every, cs$coverage), ]
t4 <- data.frame(`Survey CV` = cs$cv_surv, Design = cs$design,
                 `Every (months)` = cs$every, `Stations` = cs$n_station,
                 `Tows per year` = cs$tows_per_year,
                 `Survey RMSE` = cs$rmse, `95% interval` = sprintf("[%.3f, %.3f]", cs$lo, cs$hi),
                 `Survey MAE` = cs$mae, `Empty months` = cs$empty_months,
                 `CI coverage` = cs$cover_survey, `Power` = cs$power_survey,
                 `Nominal RMSE` = cs$rmse_nominal, `Delta RMSE` = cs$rmse_delta,
                 check.names = FALSE)
save_tab("table4_design_curve", t4)

## Table 5: classifier confusion -----------------------------------------------
cl <- readRDS("results/classifier.rds")
cm <- as.data.frame.matrix(cl$confusion); cm <- cbind(Truth = rownames(cm), cm)
save_tab("table5_classifier", cm, digits = 0)
writeLines(sprintf("accuracy %.3f indeterminate %.3f", cl$accuracy, cl$indeterminate),
           "results/tables/table5_accuracy.txt")

## Table 6: frame geometry ------------------------------------------------------
fr <- readRDS("results/frame.rds")
t6 <- do.call(rbind, lapply(split(fr, list(fr$frame, fr$scenario), drop = TRUE), function(d)
  data.frame(Frame = d$frame[1], Scenario = label_scenario(d$scenario[1]),
             `Truth (total)` = mean(d$beta_cf), `Truth (in frame)` = mean(d$beta_cf_in),
             `Frame share change` = mean(d$dp_out, na.rm = TRUE),
             Nominal = mean(d$beta_nominal, na.rm = TRUE), Delta = mean(d$beta_delta, na.rm = TRUE),
             Survey = mean(d$beta_survey, na.rm = TRUE), check.names = FALSE)))
t6 <- t6[order(t6$Scenario, factor(t6$Frame, c("closed", "open_static", "open_dynamic"))), ]
save_tab("table6_frame", t6)

## Table 7: temporal ------------------------------------------------------------
tp <- readRDS("results/temporal.rds")
row7 <- function(s, nm) data.frame(Target = nm,
  `R² acute` = mean(s$r2_acute), `R² dist. lag` = mean(s$r2_dl),
  `Corr. acute` = mean(s$cor_acute), `Corr. dist. lag` = mean(s$cor_dl),
  `OOS RMSE acute` = mean(s$oos_acute), `OOS RMSE dist. lag` = mean(s$oos_dl),
  `Effective df` = mean(s$edf), check.names = FALSE)
t7 <- rbind(row7(tp$summary_obs, "Survey index (observed)"),
            row7(tp$summary_tru, "True event effect (upper bound)"))
save_tab("table7_temporal", t7)
writeLines(sprintf("fit to survey index, correlation of fitted series with the TRUE effect: acute %.3f, distributed lag %.3f",
                   mean(tp$summary_obs$cor_obs_truth_ac), mean(tp$summary_obs$cor_obs_truth_dl)),
           "results/tables/table7_note.txt")

## Table 8: sensitivity, standardised regression coefficients -------------------
sens <- read.csv("results/sensitivity_design.csv")
s8 <- lapply(c("gap_nominal", "rmse_survey", "bias_survey_total", "sign_err_nominal"),
             function(r) src(sens, r))
r2s <- sapply(s8, function(z) z$r2[1])
nm8 <- c("Nominal gap", "Survey RMSE (in frame)", "Survey bias (total)", "Nominal sign error")
t8 <- Reduce(function(a, b) merge(a, b, by = "Parameter", all = TRUE),
             lapply(seq_along(s8), function(i) { z <- s8[[i]][, c("parameter", "src")]
               names(z) <- c("Parameter", nm8[i]); z }))
t8 <- t8[order(-abs(t8[[2]])), ]
save_tab("table8_sensitivity", t8)
writeLines(paste(sprintf("%.2f", r2s), collapse = ", "), "results/tables/table8_r2.txt")

cat("tables written to results/tables/\n")
