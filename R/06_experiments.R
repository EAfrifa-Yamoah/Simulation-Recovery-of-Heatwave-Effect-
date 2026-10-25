# =============================================================================
# 06_experiments.R
# The experiments reported in the manuscript. Every experiment benchmarks
# against the paired counterfactual estimand (beta_cf, total abundance) and,
# where the survey is involved, also against the in frame estimand (beta_cf_in).
# =============================================================================

EST_ALL <- c("nominal", "naive", "delta", "deltaint", "survey")

EST_LABELS <- c(
  nominal      = "Nominal CPUE",
  naive        = "Positive only demeaned",
  delta        = "Delta, two way fixed effects",
  delta_fished = "Delta, fished cells only",
  deltaint     = "Delta with spatial severity interaction",
  survey       = "Survey, design based",
  survey_mb    = "Survey, model based"
)

# ---- Experiment A: 2^4 factorial plus warm edge cells -----------------------
exp_factorial <- function(cfg, R = 200, seed0 = 20260724,
                          scenarios = c(scenario_levels, scenario_warm),
                          estimators = EST_ALL) {
  grid <- make_grid(cfg)
  out <- list()
  for (sc in scenarios) {
    message("  factorial: ", sc, "  ", format(Sys.time(), "%H:%M:%S"))
    d <- run_scenario(cfg, sc, R, seed0 = seed0, estimators = estimators, grid = grid)
    d$scenario <- sc
    out[[sc]] <- d
  }
  do.call(rbind, out)
}

#' Main effects and interactions over the 2^4 core design.
factorial_anova <- function(fac, response = "beta_nominal") {
  d <- fac[fac$scenario %in% scenario_levels, ]
  b <- do.call(rbind, strsplit(substr(d$scenario, 2, 5), ""))
  d$A <- as.numeric(b[, 1] == "1")   # abundance decline
  d$B <- as.numeric(b[, 2] == "1")   # redistribution
  d$C <- as.numeric(b[, 3] == "1")   # uniform catchability
  d$D <- as.numeric(b[, 4] == "1")   # fleet contraction (sampling shock)
  d$yy <- d[[response]]
  d <- d[is.finite(d$yy), ]
  fit <- lm(yy ~ A * B * C * D, data = d)
  cf <- summary(fit)$coefficients
  data.frame(response = response, term = rownames(cf),
             estimate = cf[, 1], se = cf[, 2], row.names = NULL)
}

#' Scenario level summary of every estimator against both estimands.
factorial_summary <- function(fac) {
  ests <- c("beta_nominal", "beta_naive", "beta_delta", "beta_delta_fished",
            "beta_deltaint", "beta_survey", "beta_survey_mb")
  do.call(rbind, lapply(split(fac, fac$scenario), function(d) {
    row <- data.frame(scenario = d$scenario[1], R = nrow(d),
                      truth_total = mean(d$beta_cf), truth_in = mean(d$beta_cf_in),
                      mcse_truth = sd(d$beta_cf) / sqrt(nrow(d)),
                      dp_out = mean(d$dp_out, na.rm = TRUE),
                      beta_assoc = mean(d$beta_assoc),
                      imputed_event = mean(d$imputed_event, na.rm = TRUE),
                      imputed_calm = mean(d$imputed_calm, na.rm = TRUE),
                      delta_conv = mean(d$delta_conv, na.rm = TRUE),
                      F_min = min(d$F_min, na.rm = TRUE))
    for (e in ests) if (e %in% names(d)) {
      v <- d[[e]]; ok <- is.finite(v)
      row[[e]] <- mean(v[ok]); row[[paste0("mcse_", e)]] <- sd(v[ok]) / sqrt(sum(ok))
    }
    sek <- c(beta_nominal = "se_nominal", beta_delta = "se_delta", beta_survey = "se_survey")
    for (e in names(sek)) if (sek[[e]] %in% names(d)) {
      tr <- if (e == "beta_survey") d$beta_cf_in else d$beta_cf
      pf <- perf(d[[e]], tr, d[[sek[[e]]]])
      row[[paste0("cover_", e)]] <- pf["coverage"]; row[[paste0("power_", e)]] <- pf["power"]
    }
    row
  }))
}

# ---- Experiment B: survey design curve --------------------------------------
# Coverage x survey frequency x design, for the decline plus uniform
# catchability cell with the fleet contraction on (S1011), the case in which
# the fishery is most badly misled. Reported against the estimand with RMSE,
# Monte Carlo intervals, CI coverage and power.
exp_coverage <- function(cfg, coverages = c(0.02, 0.05, 0.10, 0.20, 0.40),
                         everys = c(1, 3, 6), designs = c("srs", "strat"),
                         cv_grid = NULL, R = 150, seed0 = 771, scenario = "S1011",
                         estimators = c("nominal", "delta", "survey")) {
  grid <- make_grid(cfg)
  cvs <- if (is.null(cv_grid)) cfg$cv_surv else cv_grid
  res <- list()
  for (cv in cvs) for (ev in everys) for (des in designs) for (co in coverages) {
    # the full grid is run at the default CV; the second CV and the
    # stratified design are run at monthly frequency only
    if (cv != cvs[1] && (ev != 1 || des != "srs")) next
    if (des == "strat" && ev != 1) next
    cc <- cfg; cc$cv_surv <- cv
    message("  coverage ", co, " every ", ev, " ", des, " cv ", cv,
            "  ", format(Sys.time(), "%H:%M:%S"))
    d <- run_scenario(cc, scenario, R, seed0 = seed0, estimators = estimators,
                      grid = grid, coverage = co, every = ev, design = des)
    d$coverage <- co; d$cv_surv <- cv; d$every <- ev; d$design <- des
    d$n_station <- max(1L, round(co * nrow(grid)))
    d$tows_per_year <- d$n_station * (12 / ev)
    res[[paste(cv, ev, des, co)]] <- d
  }
  do.call(rbind, res)
}

coverage_summary <- function(cov, boot = 400) {
  spl <- split(cov, list(cov$cv_surv, cov$every, cov$design, cov$coverage), drop = TRUE)
  do.call(rbind, lapply(spl, function(d) {
    tr <- d$beta_cf
    mk <- function(e) {
      ok <- is.finite(e) & is.finite(tr)
      err <- (e - tr)[ok]
      if (length(err) < 3) return(c(rmse = NA, lo = NA, hi = NA, bias = NA, sd = NA))
      r <- sqrt(mean(err^2))
      bs <- replicate(boot, sqrt(mean(sample(err, length(err), TRUE)^2)))
      c(rmse = r, lo = unname(quantile(bs, 0.025)), hi = unname(quantile(bs, 0.975)),
        bias = mean(err), sd = sd(err), mae = median(abs(err)))
    }
    ps <- perf(d$beta_survey, tr, d$se_survey)
    pn <- perf(d$beta_nominal, tr, d$se_nominal)
    pd <- perf(d$beta_delta, tr, d$se_delta)
    zero <- mk(rep(0, nrow(d)))
    data.frame(cv_surv = d$cv_surv[1], every = d$every[1], design = d$design[1],
               coverage = d$coverage[1], n_station = d$n_station[1],
               tows_per_year = d$tows_per_year[1],
               t(mk(d$beta_survey))[, , drop = FALSE],
               cover_survey = ps["coverage"], power_survey = ps["power"],
               empty_months = mean(d$surv_empty, na.rm = TRUE),
               n_estimable = sum(is.finite(d$beta_survey)),
               rmse_nominal = mk(d$beta_nominal)["rmse"], power_nominal = pn["power"],
               cover_nominal = pn["coverage"],
               rmse_delta = mk(d$beta_delta)["rmse"], power_delta = pd["power"],
               cover_delta = pd["coverage"],
               rmse_nochange = zero["rmse"], row.names = NULL)
  }))
}

# ---- Experiment C: matched observable calibration ---------------------------
exp_matched <- function(cfg, kappas = c(0, 0.06, 0.12, 0.18, 0.24),
                        deltas = c(0, 0.035, 0.07, 0.105, 0.14),
                        R = 60, seed0 = 4242) {
  grid <- make_grid(cfg)
  res <- list()
  for (k in kappas) for (dd in deltas) {
    cc <- cfg; cc$kappa <- k; cc$delta_dep <- dd
    sc <- paste0("S", as.integer(k > 0), "0", as.integer(dd > 0), "1")
    d <- run_scenario(cc, sc, R, seed0 = seed0,
                      estimators = c("nominal", "delta", "survey"), grid = grid)
    res[[paste(k, dd)]] <- data.frame(
      kappa = k, delta_dep = dd,
      beta_truth = mean(d$beta_cf, na.rm = TRUE),
      se_truth = sd(d$beta_cf, na.rm = TRUE) / sqrt(R),
      beta_nominal = mean(d$beta_nominal, na.rm = TRUE),
      se_nominal = sd(d$beta_nominal, na.rm = TRUE) / sqrt(R),
      beta_delta = mean(d$beta_delta, na.rm = TRUE),
      beta_survey = mean(d$beta_survey, na.rm = TRUE),
      loss_max = mean(d$loss_max, na.rm = TRUE))
  }
  do.call(rbind, res)
}

# ---- Experiment D: temporal recovery --------------------------------------
# Operating model given a genuine lag and slow recovery. Two targets:
#   (i)  the survey index (observed): what an analyst can actually fit;
#   (ii) the true event effect Delta_t: an upper bound on what any model of
#        the lag structure could recover.
# For each target: a contemporaneous model and a penalised distributed lag with
# the smoothing parameter chosen by GCV INSIDE each training block of a
# blocked cross validation, so that no information from the held out block
# enters the choice of lambda.
MASS_ginv <- function(A, tol = 1e-10) { e <- eigen(A, symmetric = TRUE); d <- ifelse(e$values > tol * max(e$values), 1 / e$values, 0)
  e$vectors %*% (d * t(e$vectors)) }

exp_temporal <- function(cfg, R = 150, seed0 = 909, L = 24, coverage = 0.20) {
  cfgL <- cfg
  cfgL$dose_lag <- 3; cfgL$kappa <- 0.16; cfgL$r0 <- 0.05; cfgL$sd_proc <- 0.010
  grid <- make_grid(cfgL); Tn <- cfgL$n_months; S <- nrow(grid)
  tt <- seq_len(Tn)
  D <- diff(diag(L + 1), differences = 2)
  P <- matrix(0, L + 2, L + 2); P[2:(L + 2), 2:(L + 2)] <- crossprod(D)
  lam_grid <- 10^seq(-1, 5, length.out = 40)

  gcv_lambda <- function(X, y) {
    n <- length(y)
    sc <- sapply(lam_grid, function(lam) {
      A <- crossprod(X) + lam * P + diag(1e-8, ncol(X))
      H <- tryCatch(X %*% solve(A, t(X)), error = function(e) X %*% MASS_ginv(A) %*% t(X))
      edf <- sum(diag(H))
      n * sum((y - H %*% y)^2) / (n - edf)^2
    })
    lam_grid[which.min(sc)]
  }
  fit_one <- function(target, I) {
    ok <- is.finite(target)
    tc <- target; tc[ok] <- tc[ok] - mean(tc[ok])
    lagm <- sapply(0:L, function(l) c(rep(0, l), I[seq_len(Tn - l)]))
    Xac <- cbind(1, I); Xdl <- cbind(1, lagm)
    # rank robust solvers: a training block with no event months leaves the
    # severity column constant, in which case its coefficient is set to zero
    rsolve <- function(A, b) { z <- tryCatch(qr.coef(qr(A), b), error = function(e) rep(NA, length(b)))
      z[!is.finite(z)] <- 0; as.numeric(z) }
    f_ac <- function(idx) rsolve(crossprod(Xac[idx, ]) + diag(1e-8, 2),
                                 crossprod(Xac[idx, ], tc[idx]))
    f_dl <- function(idx, lam) rsolve(crossprod(Xdl[idx, ]) + lam * P +
                                        diag(1e-8, ncol(Xdl)),
                                      crossprod(Xdl[idx, ], tc[idx]))
    all <- which(ok)
    lam_full <- gcv_lambda(Xdl[all, ], tc[all])
    fa <- as.numeric(Xac %*% f_ac(all)); fd <- as.numeric(Xdl %*% f_dl(all, lam_full))
    r2 <- function(f) 1 - sum((tc[all] - f[all])^2) / sum((tc[all] - mean(tc[all]))^2)
    folds <- split(all, cut(seq_along(all), 4, labels = FALSE))
    oos <- function(pred_fn) {
      e <- unlist(lapply(folds, function(te) {
        tr <- setdiff(all, te); tc[te] - pred_fn(tr, te) }))
      sqrt(mean(e^2))
    }
    oos_ac <- oos(function(tr, te) as.numeric(Xac[te, , drop = FALSE] %*% f_ac(tr)))
    oos_dl <- oos(function(tr, te) {
      lam <- gcv_lambda(Xdl[tr, ], tc[tr])          # nested: lambda from training block only
      as.numeric(Xdl[te, , drop = FALSE] %*% f_dl(tr, lam)) })
    A <- crossprod(Xdl[all, ]) + lam_full * P + diag(1e-8, ncol(Xdl))
    edf <- sum(diag(Xdl[all, ] %*% solve(A, t(Xdl[all, ]))))
    list(r2_acute = r2(fa), r2_dl = r2(fd),
         cor_acute = cor(tc[all], fa[all]), cor_dl = cor(tc[all], fd[all]),
         oos_acute = oos_ac, oos_dl = oos_dl, lambda = lam_full, edf = edf,
         amp_true = diff(range(tc[all])), amp_acute = diff(range(fa[all])),
         amp_dl = diff(range(fd[all])), w_hat = f_dl(all, lam_full)[-1],
         fit_acute = fa, fit_dl = fd, target = tc)
  }

  res <- vector("list", R); keep_example <- NULL
  for (r in seq_len(R)) {
    set.seed(seed0 + r)
    sw <- scenario_switches("S1001")
    ev <- draw_events(cfgL); field <- sim_field(cfgL, grid)
    eta <- rnorm(Tn, 0, cfgL$sd_proc); fleet <- draw_fleet(cfgL, grid)
    parts <- temperature_parts(cfgL, grid, ev, field)
    J <- ev$J; I <- (J - mean(J)) / sd(J)
    rf_c <- refuge_suitability(cfgL, parts, FALSE)
    wf <- run_world(cfgL, grid, sw, parts$base + parts$event,
                    refuge_suitability(cfgL, parts, TRUE), J, eta, fleet,
                    Tmat_base = parts$base, phi_out_base = rf_c)
    wc <- run_world(cfgL, grid, sw, parts$base, rf_c,
                    rep(0, Tn), eta, fleet, E_fixed = wf$E, observe = FALSE,
                    Tmat_base = parts$base, phi_out_base = rf_c)
    truth <- wf$logN - wc$logN
    surv <- simulate_survey(cfgL, grid, exp(wf$logN), wf$pi, coverage = coverage)
    sidx <- est_survey(surv, I, grid)$index
    # remove trend and season from the survey index so both targets are on the
    # same footing (the true effect carries neither)
    Xc <- cbind(1, tt / Tn, sin(2 * pi * tt / 12), cos(2 * pi * tt / 12))
    okk <- is.finite(sidx)
    sres <- sidx; sres[okk] <- sidx[okk] - Xc[okk, ] %*% qr.solve(Xc[okk, ], sidx[okk])
    fo <- fit_one(sres, I); ft <- fit_one(truth, I)
    res[[r]] <- list(
      obs = fo[c("r2_acute","r2_dl","cor_acute","cor_dl","oos_acute","oos_dl",
                 "lambda","edf","amp_true","amp_acute","amp_dl")],
      tru = ft[c("r2_acute","r2_dl","cor_acute","cor_dl","oos_acute","oos_dl",
                 "lambda","edf","amp_true","amp_acute","amp_dl")],
      w_obs = fo$w_hat, w_tru = ft$w_hat,
      cor_obs_truth_dl = cor(fo$fit_dl, truth - mean(truth)),
      cor_obs_truth_ac = cor(fo$fit_acute, truth - mean(truth)))
    if (r == 7) keep_example <- list(month = tt, J = J, event = J > 0,
                                     truth = truth - mean(truth),
                                     survey_resid = sres,
                                     fit_acute_obs = fo$fit_acute, fit_dl_obs = fo$fit_dl,
                                     fit_acute_tru = ft$fit_acute, fit_dl_tru = ft$fit_dl,
                                     w_obs = fo$w_hat, w_tru = ft$w_hat)
  }
  tab <- function(key) do.call(rbind, lapply(res, function(z) as.data.frame(z[[key]])))
  list(summary_obs = cbind(tab("obs"),
                           cor_obs_truth_dl = sapply(res, `[[`, "cor_obs_truth_dl"),
                           cor_obs_truth_ac = sapply(res, `[[`, "cor_obs_truth_ac")),
       summary_tru = tab("tru"),
       lag_weights_obs = do.call(rbind, lapply(res, `[[`, "w_obs")),
       lag_weights_tru = do.call(rbind, lapply(res, `[[`, "w_tru")),
       example = keep_example, L = L)
}

# ---- Experiment E: global sensitivity -----------------------------------------
#' Ranges sampled by the global sensitivity design. Each must contain the
#' default value in default_config() (checked in tests/testthat/test-core.R).
SENS_RANGES <- list(
    beta_mu    = c(0.5, 3.5),  beta_J    = c(0.0, 1.6),   # contraction strength, centred on the default 0.8
    a0         = c(0.4, 1.8),  sigma_log = c(0.35, 1.0),
    cv_surv    = c(0.2, 0.9),  surv_cover = c(0.05, 0.35),
    a_surv     = c(1.0, 5.0),  belief_lr = c(0.1, 0.9),
    rho_move   = c(0.2, 0.95), n_events  = c(3, 10),
    warm_trend = c(0.0, 0.20), kappa     = c(0.04, 0.24),
    delta_dep  = c(0.0, 0.16), sigma_th  = c(1.4, 3.0),
    F_scale    = c(0, 3e-4))

exp_sensitivity <- function(cfg, n_design = 180, R = 12, seed0 = 606,
                            scenario = "S1111") {
  set.seed(seed0)
  rng <- SENS_RANGES
  k <- length(rng)
  lhs <- sapply(seq_len(k), function(j) (sample(n_design) - runif(n_design)) / n_design)
  design <- as.data.frame(mapply(function(u, rg) rg[1] + u * diff(rg),
                                 as.data.frame(lhs), rng))
  names(design) <- names(rng)
  design$n_events <- round(design$n_events)

  grid <- make_grid(cfg)
  out <- vector("list", n_design)
  for (i in seq_len(n_design)) {
    if (i %% 20 == 0) message("  sensitivity design point ", i, "/", n_design,
                              "  ", format(Sys.time(), "%H:%M:%S"))
    cc <- cfg
    for (nm in names(rng)) cc[[nm]] <- design[[nm]][i]
    d <- tryCatch(run_scenario(cc, scenario, R, seed0 = 3000,
                               estimators = c("nominal", "delta", "survey"),
                               grid = grid), error = function(e) NULL)
    if (is.null(d)) next
    out[[i]] <- cbind(design[i, ], data.frame(
      gap_nominal = mean(d$beta_nominal - d$beta_cf, na.rm = TRUE),
      gap_delta   = mean(d$beta_delta - d$beta_cf, na.rm = TRUE),
      rmse_survey = sqrt(mean((d$beta_survey - d$beta_cf_in)^2, na.rm = TRUE)),
      bias_survey_total = mean(d$beta_survey - d$beta_cf, na.rm = TRUE),
      sign_err_nominal = mean(sign(d$beta_nominal) != sign(d$beta_cf) &
                                abs(d$beta_cf) > 1e-8, na.rm = TRUE),
      dG_fishery = mean(d$dG_fishery, na.rm = TRUE),
      dG_true    = mean(d$dG_true, na.rm = TRUE)))
  }
  do.call(rbind, out)
}

SENS_RESPONSES <- c("gap_nominal", "gap_delta", "rmse_survey", "bias_survey_total",
                    "sign_err_nominal", "dG_fishery", "dG_true")

#' Standardised regression coefficients as a global sensitivity measure.
src <- function(sens, response) {
  preds <- setdiff(names(sens), SENS_RESPONSES)
  d <- sens[is.finite(sens[[response]]), ]
  z <- as.data.frame(scale(d[, preds]))
  z$y <- as.numeric(scale(d[[response]]))
  fit <- lm(y ~ ., data = z)
  cf <- summary(fit)$coefficients
  out <- data.frame(response = response, parameter = rownames(cf)[-1],
                    src = cf[-1, 1], se = cf[-1, 2], row.names = NULL)
  out$r2 <- summary(fit)$r.squared
  out[order(-abs(out$src)), ]
}

# ---- Experiment F: frame geometry -------------------------------------------
# The redistribution cells under three frame settings: closed; open with a
# static refuge (legacy thermostat); open with a dynamic refuge (default).
exp_frame <- function(cfg, R = 100, seed0 = 20260724,
                      scenarios = c("S0100", "S0101", "S1100", "S1101")) {
  settings <- list(closed = list(out_open = FALSE),
                   open_static = list(out_open = TRUE, refuge_static = TRUE),
                   open_dynamic = list(out_open = TRUE, refuge_static = FALSE))
  out <- list()
  for (nm in names(settings)) {
    cc <- modifyList(cfg, settings[[nm]]); grid <- make_grid(cc)
    for (sc in scenarios) {
      message("  frame ", nm, " ", sc, "  ", format(Sys.time(), "%H:%M:%S"))
      d <- run_scenario(cc, sc, R, seed0 = seed0,
                        estimators = c("nominal", "delta", "survey"), grid = grid)
      d$scenario <- sc; d$frame <- nm
      out[[paste(nm, sc)]] <- d
    }
  }
  do.call(rbind, out)
}

# ---- Experiment G: diagnostic key as a validated classifier -----------------
exp_classifier <- function(fac, train_frac = 0.5, seed = 99) {
  truth <- c(S0001 = "none", S1001 = "abundance", S0101 = "redistribution",
             S0011 = "uniform catchability", S1111 = "mixed",
             SPQ1 = "warm edge catchability", SPQT1 = "redistribution")
  d <- fac[fac$scenario %in% names(truth), ]
  d <- d[is.finite(d$beta_delta) & is.finite(d$beta_survey) &
           is.finite(d$dG_fishery) & is.finite(d$dG_survey), ]
  d$truth <- truth[d$scenario]
  set.seed(seed)
  itr <- sample(nrow(d), floor(train_frac * nrow(d)))
  tr <- d[itr, ]; te <- d[-itr, ]
  thr <- list(
    idx  = quantile(abs(tr$beta_delta[tr$truth == "none"]), 0.95),
    surv = quantile(abs(tr$beta_survey[tr$truth == "none"]), 0.95),
    cgf  = quantile(abs(tr$dG_fishery[tr$truth == "none"]), 0.95),
    cgs  = quantile(abs(tr$dG_survey[tr$truth == "none"]), 0.95))
  classify <- function(x) {
    idx_dn <- x$beta_delta < -thr$idx
    ab_dn  <- x$beta_survey < -thr$surv
    cgf_mv <- x$dG_fishery > thr$cgf
    cgs_mv <- x$dG_survey  > thr$cgs
    out <- rep("indeterminate", nrow(x))
    out[!idx_dn & !ab_dn & !cgf_mv & !cgs_mv] <- "none"
    out[idx_dn & ab_dn & !cgf_mv & !cgs_mv]   <- "abundance"
    out[!ab_dn & cgf_mv & cgs_mv]             <- "redistribution"
    out[idx_dn & !ab_dn & !cgf_mv & !cgs_mv]  <- "uniform catchability"
    out[idx_dn & !ab_dn & cgf_mv & !cgs_mv]   <- "warm edge catchability"
    out[idx_dn & ab_dn & cgf_mv & cgs_mv]     <- "mixed"
    out
  }
  pred <- classify(te)
  lv <- c("none", "abundance", "redistribution", "uniform catchability",
          "warm edge catchability", "mixed", "indeterminate")
  cm <- table(truth = factor(te$truth, levels = lv[1:6]),
              predicted = factor(pred, levels = lv))
  acc <- sum(diag(as.matrix(cm[, rownames(cm), drop = FALSE]))) / sum(cm)
  list(thresholds = unlist(thr), confusion = cm, accuracy = acc,
       indeterminate = mean(pred == "indeterminate"))
}
