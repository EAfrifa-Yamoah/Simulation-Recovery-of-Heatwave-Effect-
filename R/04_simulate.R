# =============================================================================
# 04_simulate.R
# Replicate driver and estimand definition.
#
# The factual world is simulated in ONE sequential loop: effort at t is
# allocated from the fleet's belief at t-1; abundance at t carries the removals
# implied by effort at t; catch at t is drawn from abundance at t. The survey
# observes the same abundance. The counterfactual world deletes the event
# anomaly and re-runs the state process with the SAME realised effort
# allocation, process noise, temperature noise field and fleet, so the pair
# differ only through the event: temperature, suitability, shares, catchability
# and, through them, fishing mortality.
#
# Two estimands are carried, because they answer different questions:
#   beta_cf     slope of  Delta_t = log N_t - log N_t^0          (total abundance)
#   beta_cf_in  slope of  Delta_t^in = log N_t(1-p_t) - log N_t^0(1-p_t^0)
#                                                                (in frame abundance)
# A fishery independent survey of the frame is unbiased for the second, not
# the first, whenever the event moves biomass across the frame boundary.
# =============================================================================

#' Run the state and observation processes for one world.
#' If `E_fixed` is supplied the fleet allocation is taken as given (counterfactual).
run_world <- function(cfg, grid, sw, Tmat, phi_out_t, J, eta, fleet,
                      E_fixed = NULL, observe = TRUE) {
  S <- nrow(grid); Tn <- ncol(Tmat)
  Js  <- scale_severity(J)
  phi <- suitability(Tmat, cfg)
  ds  <- distribution_shares(phi, cfg, sw$tracking, phi_out_t)
  q   <- catchability_field(cfg, grid, J, sw)
  dose <- thermal_dose(Tmat, ds$pi, cfg)

  E   <- if (is.null(E_fixed)) matrix(0, S, Tn) else E_fixed
  Y   <- matrix(NA_real_, S, Tn)
  logN <- numeric(Tn); Fv <- numeric(Tn)
  logN[1] <- log(cfg$N0)

  belief <- rep(0, S)
  for (t in seq_len(Tn)) {
    if (is.null(E_fixed)) {
      if (t == 1) {                              # prior belief: climatological density
        b0 <- log(ds$pi[, 1] * S + 1e-8); belief <- (b0 - mean(b0)) / (sd(b0) + 1e-8)
      }
      E[, t] <- allocate_effort_month(belief, grid, fleet, cfg, Js[t], sw$sampling)
    }
    Fv[t] <- fishing_mortality(q[, t], ds$pi[, t], E[, t], cfg, S)
    if (t > 1) {
      lagd <- if (t - cfg$dose_lag >= 1) dose[t - cfg$dose_lag] else 0
      logN[t] <- abundance_step(logN[t - 1], Fv[t], lagd, eta[t], cfg, sw$decline)
    }
    if (observe) {
      fished <- E[, t] > 0
      if (any(fished)) {
        mu <- q[fished, t] * exp(logN[t]) * ds$pi[fished, t] * S
        dc <- draw_catch_month(mu, E[fished, t], cfg)
        Y[fished, t] <- dc$y
      }
      if (is.null(E_fixed)) belief <- update_belief(belief, Y[, t], cfg)
    }
  }
  list(logN = logN, pi = ds$pi, p_out = ds$p_out, q = q, E = E, Y = Y, F = Fv,
       dose = dose, logN_in = logN + log(1 - ds$p_out))
}

simulate_replicate <- function(cfg, scenario = "S1111", seed = NULL,
                               estimators = c("nominal", "naive", "delta",
                                              "deltaint", "survey"),
                               coverage = NULL, every = NULL, design = NULL,
                               grid = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sw   <- scenario_switches(scenario)
  grid <- if (is.null(grid)) make_grid(cfg) else grid
  S <- nrow(grid); Tn <- cfg$n_months

  ## ---- common random numbers -------------------------------------------
  ev    <- draw_events(cfg)
  field <- sim_field(cfg, grid)
  eta   <- rnorm(Tn, 0, cfg$sd_proc)
  fleet <- draw_fleet(cfg, grid)
  parts <- temperature_parts(cfg, grid, ev, field)
  J     <- ev$J
  Ivec  <- if (sd(J) > 0) (J - mean(J)) / sd(J) else rep(0, Tn)

  ## ---- factual world (sequential) ---------------------------------------
  wf <- run_world(cfg, grid, sw, parts$base + parts$event,
                  refuge_suitability(cfg, parts, TRUE), J, eta, fleet)

  ## ---- counterfactual: no event, same effort, same noise ----------------
  wc <- run_world(cfg, grid, sw, parts$base,
                  refuge_suitability(cfg, parts, FALSE), rep(0, Tn), eta, fleet,
                  E_fixed = wf$E, observe = FALSE)

  dlogN    <- wf$logN - wc$logN            # total abundance effect
  dlogN_in <- wf$logN_in - wc$logN_in      # in frame abundance effect
  sf <- slope_se(dlogN, Ivec); si <- slope_se(dlogN_in, Ivec)

  ## ---- survey of the factual frame ---------------------------------------
  surv <- simulate_survey(cfg, grid, exp(wf$logN), wf$pi, coverage = coverage,
                          every = every, design = design)

  out <- list(scenario = scenario,
              beta_cf = sf$beta, beta_cf_in = si$beta,
              beta_assoc = slope(wf$logN, Ivec),          # check column only
              loss_max = min(dlogN), loss_cum = sum(dlogN),
              n_events = ev$n_events, cum_int = ev$cum_int,
              p_out = mean(wf$p_out),
              dp_out = event_contrast(wf$p_out, J),        # frame share, event contrast
              dlogN_in_minus_tot = si$beta - sf$beta,
              F_mean = mean(wf$F), F_min = min(wf$F),
              frac_fished = mean(wf$E > 0),
              dfoot = event_contrast(colMeans(wf$E > 0), J),  # footprint, event contrast
              frac_zero = mean(wf$Y[!is.na(wf$Y)] == 0),
              surv_n = surv$n_per_month, surv_every = surv$every)

  ## ---- estimators ------------------------------------------------------
  sim <- list(fish = list(Y = wf$Y, E = wf$E), grid = grid, J = J, I = Ivec)
  dfit <- NULL
  if ("nominal" %in% estimators) {
    nf <- est_nominal(sim, Ivec); out$beta_nominal <- nf$beta; out$se_nominal <- nf$se
  }
  if ("naive" %in% estimators) out$beta_naive <- est_naive(sim, Ivec)$beta
  if ("delta" %in% estimators) {
    dfit <- est_delta_fe(sim, Ivec, grid)
    out$beta_delta <- dfit$beta; out$se_delta <- dfit$se
    out$beta_delta_fished <- dfit$beta_fished
    out$delta_conv <- dfit$converged
    out$imputed_event <- dfit$imputed_event; out$imputed_calm <- dfit$imputed_calm
  }
  ifit <- NULL
  if ("deltaint" %in% estimators) {
    ifit <- est_delta_int(sim, Ivec, grid)
    out$beta_deltaint <- ifit$beta
  }
  if ("survey" %in% estimators) {
    sfit <- est_survey(surv, Ivec, grid)
    out$beta_survey <- sfit$beta; out$se_survey <- sfit$se
    out$surv_empty <- sfit$n_empty
    out$dG_survey   <- event_contrast(sfit$cg, J)
    out$dGx_survey  <- event_contrast(sfit$cgx, J)
    out$beta_survey_mb <- est_survey_model(surv, Ivec, grid)$beta
  }

  ## ---- spatial diagnostics ---------------------------------------------
  cg <- centroids(sim, if (is.null(dfit)) list(pgrid = NULL) else dfit, grid, Ivec)
  out$dG_fishery  <- event_contrast(cg$raw, J)
  out$dGx_fishery <- event_contrast(cg$raw_x, J)
  if (!is.null(ifit) && !is.null(ifit$dens)) {
    d <- ifit$dens
    out$dG_grid  <- event_contrast(as.numeric(colSums(grid$y * d) / colSums(d)), J)
    out$dGx_grid <- event_contrast(as.numeric(colSums(grid$x * d) / colSums(d)), J)
  } else { out$dG_grid <- event_contrast(cg$grid, J); out$dGx_grid <- NA_real_ }
  cg_true <- as.numeric(colSums(grid$y * wf$pi) / colSums(wf$pi))
  out$dG_true <- event_contrast(cg_true, J)
  out
}

#' Run R replicates of one scenario. Seeds are a deterministic function of the
#' replicate index only, so scenarios share random numbers. Parallel over
#' replicates with forked workers when `cores > 1`.
run_scenario <- function(cfg, scenario, R, seed0 = 1000, estimators, grid = NULL,
                         coverage = NULL, every = NULL, design = NULL,
                         cores = getOption("mhw.cores", 1L)) {
  grid <- if (is.null(grid)) make_grid(cfg) else grid
  one <- function(r) simulate_replicate(cfg, scenario, seed = seed0 + r,
                                        estimators = estimators, grid = grid,
                                        coverage = coverage, every = every,
                                        design = design)
  res <- if (cores > 1) parallel::mclapply(seq_len(R), one, mc.cores = cores)
         else lapply(seq_len(R), one)
  bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z), logical(1))
  if (any(bad)) warning(sum(bad), " replicates failed in ", scenario)
  do.call(rbind, lapply(res[!bad], function(z) as.data.frame(z, stringsAsFactors = FALSE)))
}

# ---- performance measures (ADEMP style) -------------------------------------
perf <- function(est, truth, se = NULL) {
  ok <- is.finite(est) & is.finite(truth)
  e <- est[ok]; th <- truth[ok]; R <- length(e)
  base <- c(mean = NA, bias = NA, sd = NA, rmse = NA, mcse_mean = NA,
            mcse_bias = NA, mcse_rmse = NA, sign_err = NA, q025 = NA, q975 = NA,
            coverage = NA, power = NA, R = R)
  if (R < 2) return(base)
  d <- e - th
  rmse <- sqrt(mean(d^2))
  out <- c(mean = mean(e), bias = mean(d), sd = sd(e), rmse = rmse,
           mcse_mean = sd(e) / sqrt(R), mcse_bias = sd(d) / sqrt(R),
           mcse_rmse = sqrt(var(d^2) / (4 * R * rmse^2 + 1e-12)),
           sign_err = mean(sign(e) != sign(th) & abs(th) > 1e-8),
           q025 = unname(quantile(e, 0.025)), q975 = unname(quantile(e, 0.975)),
           coverage = NA, power = NA, R = R)
  if (!is.null(se)) {
    s <- se[ok]; oks <- is.finite(s)
    out["coverage"] <- mean(abs(d[oks]) <= 1.96 * s[oks])     # 95% CI covers truth
    out["power"]    <- mean(abs(e[oks]) > 1.96 * s[oks])      # CI excludes zero
  }
  out
}
