# Run with: Rscript -e 'testthat::test_dir("tests/testthat")'  (from the repo root)
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
cfg  <- default_config()
grid <- make_grid(cfg)

test_that("event severity is non negative, zero scaled severity is zero", {
  set.seed(3); ev <- draw_events(cfg)
  expect_true(all(ev$J >= 0)); expect_true(any(ev$J > 0))
  expect_equal(length(ev$J), cfg$n_months)
  expect_true(ev$n_events <= cfg$n_events)            # merged runs counted once
  expect_equal(scale_severity(rep(0, 10)), rep(0, 10)) # counterfactual world
})

test_that("distribution shares are a proper distribution with an open frame", {
  set.seed(4); ev <- draw_events(cfg)
  parts <- temperature_parts(cfg, grid, ev, sim_field(cfg, grid))
  phi <- suitability(parts$base + parts$event, cfg)
  ds  <- distribution_shares(phi, cfg, TRUE, refuge_suitability(cfg, parts, TRUE))
  tot <- colSums(ds$pi) + ds$p_out
  expect_true(all(abs(tot - 1) < 1e-8))
  expect_true(all(ds$pi >= 0)); expect_true(all(ds$p_out >= 0))
  expect_equal(length(unique(round(refuge_suitability(cfg, parts, TRUE), 6))) > 1, TRUE,
               info = "dynamic refuge suitability varies in time")
})

test_that("tow level delta lognormal observation is mean coherent at any effort", {
  set.seed(5)
  for (E in c(1, 4)) {
    m <- rep(0.8, 40000)
    y <- draw_catch_month(m, rep(E, length(m)), cfg)$y
    expect_lt(abs(mean(y) - 0.8) / 0.8, 0.03)
  }
  # more effort, less variance of the cell catch rate
  v1 <- var(draw_catch_month(rep(0.8, 20000), rep(1, 20000), cfg)$y)
  v4 <- var(draw_catch_month(rep(0.8, 20000), rep(4, 20000), cfg)$y)
  expect_lt(v4, v1)
})

test_that("a constant shift cannot change the allocation; the contraction interaction does", {
  set.seed(6); fleet <- draw_fleet(cfg, grid); belief <- rnorm(nrow(grid))
  set.seed(1); a <- allocate_effort_month(belief, grid, fleet, cfg, 1.5, FALSE)
  set.seed(1); b <- allocate_effort_month(belief, grid, fleet, cfg, 1.5, TRUE)
  expect_false(identical(a, b))
  # contraction concentrates effort: fewer cells fished in expectation
  set.seed(2); nf0 <- mean(replicate(200, sum(allocate_effort_month(belief, grid, fleet, cfg, 2, FALSE) > 0)))
  set.seed(2); nf1 <- mean(replicate(200, sum(allocate_effort_month(belief, grid, fleet, cfg, 2, TRUE) > 0)))
  expect_lt(nf1, nf0)
})

test_that("the world is sequential and consistent: catch comes from the truth state, F >= 0", {
  set.seed(7); sw <- scenario_switches("S1111")
  ev <- draw_events(cfg); field <- sim_field(cfg, grid)
  eta <- rnorm(cfg$n_months, 0, cfg$sd_proc); fleet <- draw_fleet(cfg, grid)
  parts <- temperature_parts(cfg, grid, ev, field)
  w <- run_world(cfg, grid, sw, parts$base + parts$event,
                 refuge_suitability(cfg, parts, TRUE), ev$J, eta, fleet)
  expect_true(all(w$F >= 0))
  expect_true(all(is.finite(w$logN)))
  # the catch rates in a month are mean coherent with the abundance of that month
  t <- which(ev$J > 0)[2]; fished <- w$E[, t] > 0
  mu <- w$q[fished, t] * exp(w$logN[t]) * w$pi[fished, t] * nrow(grid)
  expect_true(cor(w$Y[fished, t], mu) > 0.2)
  # removals are the same functional in both worlds: with no event and the
  # same effort, F differs only through q and pi
  wc <- run_world(cfg, grid, sw, parts$base, refuge_suitability(cfg, parts, FALSE),
                  rep(0, cfg$n_months), eta, fleet, E_fixed = w$E, observe = FALSE)
  expect_true(all(wc$F >= 0)); expect_equal(dim(wc$E), dim(w$E))
})

test_that("the two estimands are distinct objects and agree when the frame share is fixed", {
  z0 <- simulate_replicate(cfg, "S1001", seed = 8, estimators = "survey", grid = grid)
  expect_equal(z0$beta_cf, z0$beta_cf_in)          # no tracking: p_out constant
  z1 <- simulate_replicate(cfg, "S1101", seed = 8, estimators = "survey", grid = grid)
  expect_false(isTRUE(all.equal(z1$beta_cf, z1$beta_cf_in)))
  expect_true(is.finite(z1$dp_out))
})

test_that("fishery only information is singular in the uniform direction; spatial catchability is identified", {
  z <- identifiability_experiment(cfg, seed = 9)
  expect_lt(min(abs(z$eig2_fishery)), 1e-6 * max(z$eig2_fishery))
  expect_gt(min(abs(z$eig2_both)), 1e-3 * max(z$eig2_both))
  expect_lt(diff(range(z$profile_fishery)), 1e-3)
  expect_gt(diff(range(z$profile_both)), 10)
  e3 <- sort(abs(z$eig3_fishery))
  expect_lt(e3[1], 1e-6 * e3[3]); expect_gt(e3[2], 1e-3 * e3[3])   # exactly one flat direction
  expect_gt(diff(range(z$profile_gamma_fishery)), 5)                # gamma curved from fishery alone
})

test_that("survey is unbiased for in frame abundance under redistribution; CPUE is misled by catchability", {
  skip_on_cran()
  grid2 <- make_grid(cfg)
  d <- do.call(rbind, lapply(1:30, function(r)
    as.data.frame(simulate_replicate(cfg, "S0101", seed = 100 + r,
                                     estimators = c("nominal", "survey"), grid = grid2))))
  err_in  <- mean(d$beta_survey - d$beta_cf_in)
  se      <- sd(d$beta_survey - d$beta_cf_in) / sqrt(nrow(d))
  expect_lt(abs(err_in), 3 * se + 0.003)
  d2 <- do.call(rbind, lapply(1:20, function(r)
    as.data.frame(simulate_replicate(cfg, "S0011", seed = 200 + r,
                                     estimators = c("nominal", "survey"), grid = grid2))))
  expect_lt(mean(d2$beta_nominal, na.rm = TRUE), -0.01)
  expect_lt(abs(mean(d2$beta_survey, na.rm = TRUE)), 0.01)
})

test_that("two way solvers agree with glm/lm on a small separable free problem", {
  set.seed(10)
  ni <- 6; nj <- 8; n <- 400
  i <- sample(ni, n, TRUE); j <- sample(nj, n, TRUE)
  y <- 0.5 + rnorm(ni)[i] + rnorm(nj)[j] + rnorm(n)
  f <- twoway_gauss(y, i, j, ni, nj)
  fit <- lm(y ~ factor(i) + factor(j))
  expect_lt(max(abs((f$mu + f$a[i] + f$b[j]) - fitted(fit))), 1e-6)
  expect_true(f$converged)
  eta <- rnorm(ni)[i] + rnorm(nj)[j]; yb <- rbinom(n, 1, plogis(eta))
  g <- twoway_logit(yb, i, j, ni, nj, ridge = 1e-6)
  gfit <- glm(yb ~ factor(i) + factor(j), family = binomial())
  expect_lt(max(abs(plogis(g$a[i] + g$b[j]) - fitted(gfit))), 1e-3)
  expect_true(g$converged)
})

test_that("slope_se returns a finite HAC standard error and handles rank deficiency", {
  set.seed(11); x <- rnorm(96); y <- 0.3 * x + cumsum(rnorm(96, 0, 0.3))
  f <- slope_se(y, x)
  expect_true(is.finite(f$beta) && is.finite(f$se) && f$se > 0)
  expect_equal(slope(y, x), f$beta)
  y6 <- y; y6[!(seq_len(96) %in% seq(1, 96, 6))] <- NA     # semiannual sampling
  f6 <- slope_se(y6, x)
  expect_true(is.finite(f6$beta)); expect_lt(f6$rank, 5)
})

test_that("sensitivity ranges contain the defaults and keep contraction non-negative", {
  cfg <- default_config()
  for (nm in names(SENS_RANGES)) {
    r <- SENS_RANGES[[nm]]
    expect_true(r[1] <= cfg[[nm]] && cfg[[nm]] <= r[2], info = nm)
  }
  # beta_J multiplies targeting strength; a negative value would reverse targeting
  expect_gte(SENS_RANGES$beta_J[1], 0)
})

test_that("without redistribution the paired contrast is zero before the first heatwave anomaly", {
  cfg <- default_config(); grid <- make_grid(cfg); Tn <- cfg$n_months
  for (seed in c(5, 11, 23)) {
    set.seed(seed); sw <- scenario_switches("S1001")
    ev <- draw_events(cfg); field <- sim_field(cfg, grid); eta <- rnorm(Tn, 0, cfg$sd_proc)
    fleet <- draw_fleet(cfg, grid); parts <- temperature_parts(cfg, grid, ev, field)
    rf_c <- refuge_suitability(cfg, parts, FALSE)
    wf <- run_world(cfg, grid, sw, parts$base + parts$event, refuge_suitability(cfg, parts, TRUE),
                    ev$J, eta, fleet, Tmat_base = parts$base, phi_out_base = rf_c)
    wc <- run_world(cfg, grid, sw, parts$base, rf_c, rep(0, Tn), eta, fleet, E_fixed = wf$E,
                    observe = FALSE, Tmat_base = parts$base, phi_out_base = rf_c)
    first <- which(colSums(parts$event) > 0)[1]
    expect_equal(wf$p_out[1], wc$p_out[1])
    if (first > 2) expect_equal(wf$logN[1:(first - 1)], wc$logN[1:(first - 1)])
  }
})
