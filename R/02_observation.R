# =============================================================================
# 02_observation.R
# Observation and sampling process.
#
# The fishery observation is tow level: a cell receiving E_st effort units is
# observed through E_st tows, so the encounter probability and the variance of
# the cell's catch rate both depend on effort. The delta lognormal is
# parameterised so that the unconditional mean equals the latent catch rate,
# E[Y_st] = mu_st, whatever the effort.
#
# Effort allocation and catch drawing are exposed as month level functions and
# called from the sequential driver in 04_simulate.R, so that effort at t uses
# information available before t, abundance at t carries the removals implied
# by effort at t, and catch at t is drawn from the abundance the truth is
# defined on.
# =============================================================================

#' Vessel structure, drawn once per replicate: per cell vessel assignment and
#' vessel level log efficiency.
draw_fleet <- function(cfg, grid) {
  vessel_eff <- rnorm(cfg$n_vessel, 0, cfg$sd_vessel)
  vid <- sample.int(cfg$n_vessel, nrow(grid), replace = TRUE)
  list(eff = vessel_eff[vid])
}

#' Multinomial allocation of a fixed monthly effort total, from the fleet's
#' current belief about cell catch rates (lagged, standardised), a travel cost,
#' vessel efficiency and, when the sampling switch is on, a direct event driven
#' contraction of the footprint.
allocate_effort_month <- function(belief, grid, fleet, cfg, Js_t, sampling) {
  # The contraction enters as an interaction with targeting strength: a
  # constant shift in the score would cancel in the softmax and do nothing.
  targ <- cfg$beta_mu * (1 + if (sampling) cfg$beta_J * Js_t else 0)
  lin <- targ * belief - cfg$beta_dist * grid$dist * 10 + fleet$eff
  w <- exp(lin - max(lin)); w <- w / sum(w)
  as.numeric(rmultinom(1, cfg$E_total, w))
}

#' Tow level delta lognormal catch for the fished cells of one month.
#' m is the latent catch rate per unit effort, E the effort units.
#' Encounter: at least one of E tows encounters, P = 1 - exp(-a0 m E).
#' Positive part: the mean of E tows, lognormal on the log scale with variance
#' sigma_log^2 / E and location chosen so that E[Y] = P * (m / P) = m.
draw_catch_month <- function(m, E, cfg) {
  P <- 1 - exp(-cfg$a0 * m * E)
  P <- pmin(pmax(P, 1e-6), 1 - 1e-9)
  hit <- runif(length(m)) < P
  s2 <- cfg$sigma_log^2 / E
  y <- rep(0, length(m))
  if (any(hit))
    y[hit] <- exp(rnorm(sum(hit), log(m[hit] / P[hit]) - 0.5 * s2[hit], sqrt(s2[hit])))
  list(y = y, hit = hit)
}

#' Belief update: last observed log catch rate where fished, carried forward
#' where not, standardised, blended with the previous belief.
update_belief <- function(belief, Y_t, cfg) {
  newb <- ifelse(is.na(Y_t), belief, log(Y_t + 1e-4))
  newb <- (newb - mean(newb)) / (sd(newb) + 1e-8)
  (1 - cfg$belief_lr) * belief + cfg$belief_lr * newb
}

#' Fishery independent survey. Cells are sampled without replacement every
#' `every` months, either by simple random sampling (design "srs") or with
#' equal allocation across latitude bands (design "strat"). Both give equal
#' inclusion probability n_t / S, so the Horvitz Thompson total is the same
#' estimator; stratification changes its variance. Survey catchability and
#' encounter scale are fixed by protocol and do not respond to the event.
simulate_survey <- function(cfg, grid, N, pi_s, coverage = NULL, every = NULL,
                            design = NULL) {
  S <- nrow(grid); Tn <- ncol(pi_s)
  cov  <- if (is.null(coverage)) cfg$surv_cover else coverage
  ev   <- if (is.null(every)) cfg$surv_every else every
  des  <- if (is.null(design)) cfg$surv_design else design
  n_t  <- max(1L, round(cov * S))
  p_inc <- n_t / S

  V     <- matrix(NA_real_, S, Tn)
  samp  <- matrix(FALSE, S, Tn)
  sigma <- sqrt(log(1 + cfg$cv_surv^2))
  bands <- split(seq_len(S), grid$lat)
  per_band <- n_t %/% length(bands)

  for (t in seq_len(Tn)) {
    if ((t - 1) %% ev != 0) next
    if (des == "strat" && per_band >= 1) {
      idx <- unlist(lapply(bands, function(b) b[sample.int(length(b), per_band)]))
      if (length(idx) < n_t)                      # top up the remainder at random
        idx <- c(idx, sample(setdiff(seq_len(S), idx), n_t - length(idx)))
    } else idx <- sample.int(S, n_t)
    samp[idx, t] <- TRUE
    dens <- cfg$q_surv * N[t] * pi_s[idx, t] * S
    p <- 1 - exp(-cfg$a_surv * dens)
    p <- pmin(pmax(p, 1e-6), 1 - 1e-9)
    hit <- runif(length(idx)) < p
    vv <- rep(0, length(idx))
    if (any(hit))
      vv[hit] <- exp(rnorm(sum(hit), log(dens[hit] / p[hit]) - 0.5 * sigma^2, sigma))
    V[idx, t] <- vv
  }
  list(V = V, sampled = samp, p_inc = p_inc, n_per_month = n_t, every = ev,
       design = des)
}
