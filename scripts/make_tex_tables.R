# =============================================================================
# scripts/make_tex_tables.R
# LaTeX bodies for the manuscript tables, written from results/tables/*.csv to
# results/tables/tex/<name>.tex. The manuscript embeds these fragments, so every
# number in a table can be regenerated from the results with one command:
#   Rscript scripts/make_tex_tables.R
# =============================================================================
Sys.setlocale("LC_CTYPE", "C.UTF-8")
dir.create("results/tables/tex", showWarnings = FALSE, recursive = TRUE)
num <- function(x, d = 3) {
  out <- formatC(round(x, d), digits = d, format = "f")
  out[abs(round(x, d)) < 0.5 * 10^-d] <- formatC(0, digits = d, format = "f")
  out <- ifelse(substr(out, 1, 1) == "-", paste0("$-$", substring(out, 2)), out)
  out[is.na(x)] <- "--"
  out
}
rows <- function(df) paste0(apply(df, 1, paste, collapse = " & "), " \\\\")
write_tab <- function(nm, lines) { writeLines(lines, file.path("results/tables/tex", paste0(nm, ".tex"))); invisible(nm) }

plain <- c("None" = "None", "Abundance" = "Abundance decline", "Redistribution" = "Redistribution",
           "Catchability" = "Uniform catchability loss",
           "Abundance + Redistribution" = "Decline + redistribution",
           "Abundance + Catchability" = "Decline + catchability loss",
           "Redistribution + Catchability" = "Redistribution + catchability loss",
           "Abundance + Redistribution + Catchability" = "All three",
           "Warm edge catchability" = "Warm-edge catchability loss",
           "Warm edge catch. + redistribution" = "Warm-edge loss + redistribution")
split_label <- function(lab) {
  contr <- grepl("\\(fleet contraction\\)", lab)
  base <- trimws(sub("\\(fleet contraction\\)", "", lab))
  list(base = unname(plain[base]), contr = contr)
}

## Table: heatwave coefficients ------------------------------------------------
t1 <- read.csv("results/tables/table1_coefficients.csv", check.names = FALSE)
sl <- split_label(t1$Scenario)
numcols <- c("Truth (total)", "Truth (in frame)", "Nominal", "Positive only", "Delta",
             "Delta, fished cells", "Delta + interaction", "Survey")
body <- function(idx, with_contr) {
  d <- data.frame(Scenario = sl$base[idx], check.names = FALSE)
  if (with_contr) d$Contraction <- ifelse(sl$contr[idx], "yes", "no")
  for (cc in numcols) d[[cc]] <- num(t1[[cc]][idx])
  rows(d)
}
ord <- c("None", "Abundance decline", "Redistribution", "Uniform catchability loss",
         "Decline + redistribution", "Decline + catchability loss", "Redistribution + catchability loss",
         "All three", "Warm-edge catchability loss", "Warm-edge loss + redistribution")
i_main <- unlist(lapply(ord, function(o) which(sl$base == o & sl$contr)))
i_full <- unlist(lapply(ord, function(o) which(sl$base == o)))
write_tab("coef_main", body(i_main, FALSE))
write_tab("coef_full", body(i_full, TRUE))

## Table: centre of gravity -----------------------------------------------------
t3 <- read.csv("results/tables/table3_centroids.csv", check.names = FALSE)
s3 <- split_label(t3$Scenario)
# centre-of-gravity shifts are in units of the standardised northing; report km
for (f in list.files("R", full.names = TRUE)) source(f)
cfg <- default_config(); KM <- cfg$cell_km * sd(make_grid(cfg)$lat)
d3 <- data.frame(Scenario = s3$base, check.names = FALSE)
for (cc in c("Latent truth", "Fishery catch centroid", "Delta surface centroid", "Survey centroid"))
  d3[[cc]] <- num(t3[[cc]] * KM, 1)
for (cc in c("Frame share change", "Footprint change")) d3[[cc]] <- num(t3[[cc]])
write_tab("centroids", rows(d3))

## Table: ANOVA ------------------------------------------------------------------
t2 <- read.csv("results/tables/table2_anova.csv", check.names = FALSE, encoding = "UTF-8")
names(t2) <- c("Term", "Nominal", "seN", "Delta", "seD", "Survey", "seS")
t2$Term <- gsub("Abundance", "Decline", t2$Term)
t2$Term <- gsub("Catchability", "Uniform catchability", t2$Term)
t2$Term <- gsub("\u00d7", "$\\\\times$", t2$Term)
d2 <- data.frame(Term = t2$Term,
                 Nominal = paste0(num(t2$Nominal, 4), " (", num(t2$seN, 4), ")"),
                 Delta = paste0(num(t2$Delta, 4), " (", num(t2$seD, 4), ")"),
                 Survey = paste0(num(t2$Survey, 4), " (", num(t2$seS, 4), ")"))
write_tab("anova", rows(d2))

## Table: survey design curve -----------------------------------------------------
t4 <- read.csv("results/tables/table4_design_curve.csv", check.names = FALSE)
freq <- c(`1` = "monthly", `3` = "quarterly", `6` = "twice-yearly")
d4 <- data.frame(CV = formatC(t4$`Survey CV`, digits = 2, format = "f"),
                 Design = ifelse(t4$Design == "srs", "simple random", "stratified"),
                 Frequency = freq[as.character(t4$`Every (months)`)],
                 Stations = t4$Stations, Tows = t4$`Tows per year`,
                 RMSE = num(t4$`Survey RMSE`),
                 CI = gsub("-", "$-$", gsub("\\[|\\]", "", t4$`95% interval`)),
                 MAE = num(t4$`Survey MAE`), Empty = formatC(t4$`Empty months`, digits = 2, format = "f"),
                 Coverage = formatC(t4$`CI coverage`, digits = 2, format = "f"),
                 Power = formatC(t4$Power, digits = 2, format = "f"), check.names = FALSE)
d4$CI <- paste0("[", d4$CI, "]")
write_tab("design_curve", rows(d4))
writeLines(sprintf("%s %s", num(t4$`Nominal RMSE`[1]), num(t4$`Delta RMSE`[1])),
           "results/tables/tex/design_curve_fishery_rmse.txt")

## Table: frame geometry -----------------------------------------------------------
t6 <- read.csv("results/tables/table6_frame.csv", check.names = FALSE)
s6 <- split_label(t6$Scenario)
fr <- c(closed = "closed", open_static = "open, static refuge", open_dynamic = "open, warming refuge")
d6 <- data.frame(Scenario = s6$base, Contraction = ifelse(s6$contr, "yes", "no"),
                 Frame = unname(fr[t6$Frame]), check.names = FALSE)
for (cc in c("Truth (total)", "Truth (in frame)", "Frame share change", "Nominal", "Delta", "Survey"))
  d6[[cc]] <- num(t6[[cc]])
d6 <- d6[order(d6$Scenario != "Redistribution", d6$Contraction, match(d6$Frame, fr)), ]
write_tab("frame", rows(d6))

## Table: temporal ------------------------------------------------------------------
t7 <- read.csv("results/tables/table7_temporal.csv", check.names = FALSE)
tp <- readRDS("results/temporal.rds")
# correlation of the fitted series with the TRUE effect: for the survey-index
# target it is a separate quantity; for the true-effect target it equals the
# correlation with the target
ctrue_a <- c(mean(tp$summary_obs$cor_obs_truth_ac), t7[[4]][2])
ctrue_d <- c(mean(tp$summary_obs$cor_obs_truth_dl), t7[[5]][2])
d7 <- data.frame(Target = c("Survey index (detrended)", "True heatwave effect"),
                 R2a = num(t7[[2]], 2), R2d = num(t7[[3]], 2), Ca = num(t7[[4]], 2), Cd = num(t7[[5]], 2),
                 Ta = num(ctrue_a, 2), Td = num(ctrue_d, 2),
                 Oa = num(t7[[6]], 2), Od = num(t7[[7]], 2), edf = formatC(t7[[8]], digits = 1, format = "f"))
write_tab("temporal", rows(d7))

## Table: classifier ------------------------------------------------------------------
t5 <- read.csv("results/tables/table5_classifier.csv", check.names = FALSE)
cap1 <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))
lab5 <- c(none = "No effect", abundance = "Abundance decline", redistribution = "Redistribution",
          `uniform catchability` = "Uniform catchability", `warm edge catchability` = "Warm-edge catchability",
          mixed = "Mixed")
d5 <- data.frame(Truth = unname(lab5[t5$Truth]), t5[, -1], check.names = FALSE)
d5$Total <- rowSums(t5[, -1])
d5$Correct <- sprintf("%.0f", 100 * diag(as.matrix(t5[, 2:7])) / d5$Total)
write_tab("classifier", rows(d5))

## Table: sensitivity -----------------------------------------------------------------
t8 <- read.csv("results/tables/table8_sensitivity.csv", check.names = FALSE)
pl8 <- c(beta_J = "Fleet contraction strength ($\\beta_J$)", rho_move = "Movement rate ($\\rho_m$)",
         delta_dep = "Uniform catchability loss ($\\delta_{\\text{dep}}$)", n_events = "Number of events ($K$)",
         sigma_th = "Niche breadth ($\\sigma_{\\text{th}}$)", kappa = "Abundance sensitivity ($\\kappa$)",
         beta_mu = "Targeting strength ($\\beta_\\mu$)", warm_trend = "Warming trend ($\\omega$)",
         belief_lr = "Fleet learning rate ($\\psi$)", sigma_log = "Catch variability ($\\sigma_{\\log}$)",
         F_scale = "Fishing mortality scale ($f$)", cv_surv = "Survey CV",
         a_surv = "Survey encounter scale ($a^{\\text{surv}}$)", a0 = "Fishery encounter scale ($a_0$)",
         surv_cover = "Survey coverage ($c$)")
d8 <- data.frame(Parameter = unname(pl8[t8$Parameter]), check.names = FALSE)
for (k in 2:5) d8[[k]] <- num(t8[[k]])
write_tab("sensitivity", rows(d8))
file.copy("results/tables/table8_r2.txt", "results/tables/tex/sensitivity_r2.txt", overwrite = TRUE)
cat("LaTeX table bodies written to results/tables/tex/\n")
