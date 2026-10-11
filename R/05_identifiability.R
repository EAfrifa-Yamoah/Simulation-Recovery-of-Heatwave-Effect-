# =============================================================================
# 05_identifiability.R
#
# Analytic statement (Section 2.1 of the manuscript). Fishery observations
# depend on the state only through the latent catch rate
#   mu_{s,t} = q_{s,t} N_t pi_{s,t}.
# For any positive sequence a_t the reparameterisation N_t* = a_t N_t,
# q_{s,t}* = q_{s,t} / a_t leaves every mu_{s,t} unchanged and therefore the
# whole fishery likelihood, including the effort allocation, which responds to
# observed catch rates. Total abundance and spatially UNIFORM catchability are
# structurally non identifiable from fishery data alone. A survey whose
# catchability does not carry a_t breaks the invariance.
#
# What the code below does, and does not, show. The reduced likelihood used
# here is written in log mu = base + (betaN - betaQ) I + gamma g_s I, so the
# singularity of the (betaN, betaQ) block is a direct consequence of the
# functional form and is an ILLUSTRATION of the analytic result, not an
# independent test of it. The informative part is the third parameter: a
# spatially structured catchability slope gamma multiplies a cell level weight
# g_s and is NOT annihilated by a scalar a_t. The fishery information matrix
# should therefore have exactly one zero eigenvalue (the uniform direction) and
# a strictly positive eigenvalue in the spatial direction. That is the
# practical content of the manuscript's claim that warm edge catchability is a
# separate mechanism from uniform catchability, and it is what is checked.
# =============================================================================

#' Reduced delta lognormal fishery log likelihood under
#'   log N_t = log N0 + betaN I_t,  log q_{s,t} = -betaQ I_t - gamma g_s I_t.
loglik_fishery <- function(par, dat, cfg) {
  betaN <- par[1]; betaQ <- par[2]; gamma <- if (length(par) > 2) par[3] else 0
  mu <- exp(dat$base + (betaN - betaQ) * dat$I - gamma * dat$g * dat$I)
  p  <- 1 - exp(-cfg$a0 * mu)
  p  <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  ll <- sum(log(1 - p[!dat$pos])) + sum(log(p[dat$pos]))
  m  <- mu[dat$pos] / p[dat$pos]
  ll + sum(dnorm(dat$ly[dat$pos],
                 log(m) - 0.5 * cfg$sigma_log^2, cfg$sigma_log, log = TRUE))
}

#' The survey observes q_surv N_t pi_{s,t}: betaN enters, betaQ and gamma do not.
loglik_survey <- function(par, sdat, cfg) {
  betaN <- par[1]
  sig <- sqrt(log(1 + cfg$cv_surv^2))
  mu  <- exp(sdat$base + betaN * sdat$I)
  p   <- 1 - exp(-cfg$a_surv * mu)
  p   <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  ll <- sum(log(1 - p[!sdat$pos])) + sum(log(p[sdat$pos]))
  m  <- mu[sdat$pos] / p[sdat$pos]
  ll + sum(dnorm(sdat$ly[sdat$pos], log(m) - 0.5 * sig^2, sig, log = TRUE))
}

#' Numerical observed information (central differences), symmetrised.
observed_information <- function(f, par, h = 1e-3) {
  k <- length(par); H <- matrix(0, k, k)
  for (i in seq_len(k)) for (j in seq_len(k)) {
    e_i <- e_j <- numeric(k); e_i[i] <- h; e_j[j] <- h
    H[i, j] <- (f(par + e_i + e_j) - f(par + e_i - e_j) -
                  f(par - e_i + e_j) + f(par - e_i - e_j)) / (4 * h^2)
  }
  -(H + t(H)) / 2
}

#' Build one data set under the operating model's spatial structure and
#' preferential effort, and return information matrices and profiles for the
#' two parameter (uniform) and three parameter (uniform + spatial) problems.
identifiability_experiment <- function(cfg, seed = 11, betaN = -0.30,
                                       betaQ = 0.30, gamma = 0.25, grid = NULL) {
  set.seed(seed)
  grid <- if (is.null(grid)) make_grid(cfg) else grid
  S <- nrow(grid); Tn <- cfg$n_months
  ev   <- draw_events(cfg)
  I    <- (ev$J - mean(ev$J)) / sd(ev$J)
  parts <- temperature_parts(cfg, grid, ev, sim_field(cfg, grid))
  phi  <- suitability(parts$base, cfg)
  pi_s <- sweep(phi, 2, colSums(phi), "/")
  h    <- pmax(1 - cfg$gamma_sp * grid$y, 0); g <- h / mean(h)   # warm edge weight

  base <- log(pi_s * S) + log(cfg$q0)
  Imat <- matrix(I, S, Tn, byrow = TRUE)
  Gmat <- matrix(g, S, Tn)
  mu   <- exp(base + (betaN - betaQ) * Imat - gamma * Gmat * Imat)

  ## fishery draw with preferential effort: cells fished with probability
  ## increasing in climatological density (not uniform random fishing)
  pref <- plogis(2 * as.numeric(scale(log(rowMeans(pi_s) * S))) - 0.6)
  fished <- matrix(runif(S * Tn) < pref, S, Tn)
  p  <- 1 - exp(-cfg$a0 * mu)
  pos <- (runif(S * Tn) < p) & fished
  ly <- matrix(NA_real_, S, Tn)
  ly[pos] <- rnorm(sum(pos), log(mu[pos] / p[pos]) - 0.5 * cfg$sigma_log^2,
                   cfg$sigma_log)
  dat <- list(base = base[fished], I = Imat[fished], g = Gmat[fished],
              pos = pos[fished], ly = ly[fished])

  ## survey draw at the default coverage
  sig <- sqrt(log(1 + cfg$cv_surv^2))
  ssel <- matrix(runif(S * Tn) < cfg$surv_cover, S, Tn)
  musv <- exp(base + betaN * Imat) * cfg$q_surv
  psv  <- 1 - exp(-cfg$a_surv * musv)
  spos <- (runif(S * Tn) < psv) & ssel
  sly  <- matrix(NA_real_, S, Tn)
  sly[spos] <- rnorm(sum(spos), log(musv[spos] / psv[spos]) - 0.5 * sig^2, sig)
  sdat <- list(base = base[ssel], I = Imat[ssel], pos = spos[ssel], ly = sly[ssel])

  ## two parameter problem (uniform catchability only)
  par2 <- c(betaN, betaQ)
  f2 <- function(p) loglik_fishery(c(p, gamma), dat, cfg)
  I2_fish <- observed_information(f2, par2)
  I2_both <- observed_information(function(p) f2(p) + loglik_survey(p, sdat, cfg), par2)

  ## three parameter problem (uniform + spatial catchability)
  par3 <- c(betaN, betaQ, gamma)
  I3_fish <- observed_information(function(p) loglik_fishery(p, dat, cfg), par3)
  I3_both <- observed_information(function(p) loglik_fishery(p, dat, cfg) +
                                    loglik_survey(p, sdat, cfg), par3)

  ## profile likelihood in betaN (gamma profiled out as well in the 3 par case)
  gridN <- seq(betaN - 0.45, betaN + 0.45, length.out = 41)
  prof <- function(withSurvey) {
    sapply(gridN, function(bn) {
      o <- optim(c(betaQ, gamma), function(z) {
        v <- loglik_fishery(c(bn, z[1], z[2]), dat, cfg)
        if (withSurvey) v <- v + loglik_survey(c(bn, z[1], z[2]), sdat, cfg)
        -v }, method = "L-BFGS-B", lower = c(-2, -2), upper = c(2, 2))
      -o$value
    })
  }
  pf <- prof(FALSE); pb <- prof(TRUE)

  ## profile in gamma from fishery data alone (should be curved: identified)
  gridG <- seq(gamma - 0.3, gamma + 0.3, length.out = 31)
  pg <- sapply(gridG, function(gm) {
    o <- optimize(function(d) -loglik_fishery(c(betaN, betaN - d, gm), dat, cfg),
                  interval = c(-2, 2))     # profile over the identified difference
    -o$objective
  })

  list(betaN = betaN, betaQ = betaQ, gamma = gamma,
       grid_betaN = gridN, grid_gamma = gridG,
       info2_fishery = I2_fish, info2_both = I2_both,
       info3_fishery = I3_fish, info3_both = I3_both,
       eig2_fishery = eigen(I2_fish, only.values = TRUE)$values,
       eig2_both    = eigen(I2_both, only.values = TRUE)$values,
       eig3_fishery = eigen(I3_fish, only.values = TRUE)$values,
       eig3_both    = eigen(I3_both, only.values = TRUE)$values,
       profile_fishery = pf - max(pf), profile_both = pb - max(pb),
       profile_gamma_fishery = pg - max(pg),
       # legacy names used by figures and tests
       eig_fishery = eigen(I2_fish, only.values = TRUE)$values,
       eig_both = eigen(I2_both, only.values = TRUE)$values)
}
