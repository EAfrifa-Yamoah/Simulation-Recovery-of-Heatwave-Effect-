# =============================================================================
# scripts/review_checks.R
# Three checks made in response to the code review of 23 October 2026, each
# on the manuscript seeds. Results are quoted in Supporting Information S2 and
# S3 and in the response to the review.
#   Rscript scripts/review_checks.R        (about 15 minutes on one core)
# =============================================================================
Sys.setlocale("LC_CTYPE", "C.UTF-8")
for (f in list.files("R", full.names = TRUE)) source(f)
cfg <- default_config(); grid <- make_grid(cfg); S <- nrow(grid); Tn <- cfg$n_months
dir.create("results/review_checks", showWarnings = FALSE)

paired <- function(scenario, seed, cfg. = cfg) {
  set.seed(seed); sw <- scenario_switches(scenario)
  ev <- draw_events(cfg.); field <- sim_field(cfg., grid); eta <- rnorm(Tn, 0, cfg.$sd_proc)
  fleet <- draw_fleet(cfg., grid); parts <- temperature_parts(cfg., grid, ev, field)
  rf_c <- refuge_suitability(cfg., parts, FALSE)
  wf <- run_world(cfg., grid, sw, parts$base + parts$event, refuge_suitability(cfg., parts, TRUE),
                  ev$J, eta, fleet, Tmat_base = parts$base, phi_out_base = rf_c)
  wc <- run_world(cfg., grid, sw, parts$base, rf_c, rep(0, Tn), eta, fleet, E_fixed = wf$E,
                  observe = FALSE, Tmat_base = parts$base, phi_out_base = rf_c)
  list(wf = wf, wc = wc, ev = ev, parts = parts, I = (ev$J - mean(ev$J)) / sd(ev$J))
}

## 1. Paired contrast before the first heatwave anomaly (should be exactly zero)
z <- t(sapply(1:40, function(i) { p <- paired("S1001", 20260724 + i)
  first <- which(colSums(p$parts$event) > 0)[1]; d <- p$wf$logN - p$wc$logN
  c(first_anomaly = first, max_abs_contrast_before = if (first > 2) max(abs(d[1:(first - 1)])) else NA) }))
cat("Check 1: max |paired contrast| before the first anomaly, 40 replicates of S1001:",
    max(z[, 2], na.rm = TRUE), "\n")

## 2. Thermal stress omitted beyond the frame (Supporting Information S2)
z2 <- t(sapply(1:40, function(i) { p <- paired("S1101", 20260724 + i)
  Tref <- p$parts$refuge_base + p$parts$refuge_event
  d_in <- p$wf$dose; d_out <- p$wf$p_out * pmax(Tref - cfg$T_crit, 0)
  c(frac_months_refuge_stressed = mean(Tref > cfg$T_crit), omitted_share = sum(d_out) / (sum(d_in) + sum(d_out)),
    max_p_out = max(p$wf$p_out), cumulative_effect_logN = cfg$kappa * sum(d_out)) }))
cat("Check 2 (mean, max over 40 replicates):\n"); print(rbind(mean = colMeans(z2), max = apply(z2, 2, max)), digits = 3)
write.csv(z2, "results/review_checks/refuge_stress.csv", row.names = FALSE)

## 3. Effort-adjusted delta index (Supporting Information S3)
est_delta_effort <- function(sim) {
  Y <- sim$fish$Y; E <- sim$fish$E
  f <- !is.na(Y); d <- data.frame(y = Y[f], e = E[f], cell = factor(row(Y)[f], levels = 1:S),
                                  tim = factor(col(Y)[f], levels = 1:Tn)); d$enc <- as.numeric(d$y > 0)
  m1 <- suppressWarnings(glm(enc ~ cell + tim, family = binomial("cloglog"), offset = log(e), data = d,
                             control = glm.control(maxit = 100)))
  m2 <- lm(log(y) ~ cell + tim + log(e), data = d[d$enc == 1, ])
  Eref <- mean(d$e)
  nd <- expand.grid(cell = factor(1:S, levels = 1:S), tim = factor(1:Tn, levels = 1:Tn)); nd$e <- Eref
  cf1 <- coef(m1); cf1[is.na(cf1)] <- 0; X1 <- model.matrix(~ cell + tim, nd)
  p <- 1 - exp(-exp(as.numeric(X1 %*% cf1[colnames(X1)]) + log(Eref)))
  X2 <- model.matrix(~ cell + tim + log(e), nd); cf2 <- coef(m2)[match(colnames(X2), names(coef(m2)))]; cf2[is.na(cf2)] <- 0
  mu <- exp(as.numeric(X2 %*% cf2) + 0.5 * summary(m2)$sigma^2)
  idx <- log(colMeans(matrix(p * mu, S, Tn))); slope_se(idx, sim$I)$beta
}
rows <- list()
for (sc in c("S1011", "S1001", "S0011", "S0001")) {
  z3 <- t(sapply(1:40, function(i) { p <- paired(sc, 20260724 + i)
    sim <- list(fish = list(Y = p$wf$Y, E = p$wf$E), grid = grid, J = p$ev$J, I = p$I)
    c(truth = slope_se(p$wf$logN - p$wc$logN, p$I)$beta, delta = est_delta_fe(sim, p$I, grid)$beta,
      delta_effort = est_delta_effort(sim),
      effort_per_fished_cell_event_minus_calm = mean(p$wf$E[p$wf$E > 0 & rep(p$ev$J > 0, each = S)]) -
        mean(p$wf$E[p$wf$E > 0 & rep(p$ev$J == 0, each = S)])) }))
  rows[[sc]] <- data.frame(scenario = sc, truth = mean(z3[, 1]), delta = mean(z3[, 2]), delta_se = sd(z3[, 2]) / sqrt(40),
                           delta_effort = mean(z3[, 3]), delta_effort_se = sd(z3[, 3]) / sqrt(40),
                           difference = mean(z3[, 3] - z3[, 2]), d_effort_per_cell = mean(z3[, 4]))
}
out <- do.call(rbind, rows); print(out, digits = 3)
write.csv(out, "results/review_checks/effort_adjusted_delta.csv", row.names = FALSE)
