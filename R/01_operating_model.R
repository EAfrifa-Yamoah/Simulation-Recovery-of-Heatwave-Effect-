# =============================================================================
# 01_operating_model.R
# State process for the marine heatwave CPUE attribution framework.
#
# Revision notes (responses to review items):
#   P6  / 3.3  event severity J_t >= 0 drives every biological and fleet
#              response; the centred covariate I_t is used only in estimation.
#   P3.4       explicit event generating algorithm with recorded count,
#              duration, and category weighted cumulative intensity.
#   P3.5       temperature decomposed into climatology, chronic trend and
#              event anomaly, so counterfactual pairs are well defined.
#   P11 / 3.6  abundance is a state transition with density dependence,
#              fishing mortality and a lagged thermal dose.
#   P3.7       thermal dose is occupancy weighted, not domain mean.
#   P3.8       partial adjustment movement with an explicit outside frame state.
#   P3.9       the cross shelf axis carries a real depth gradient.
#   P3.10      spatial catchability weight normalised to domain mean one.
#   P3.16      grid dimensions, cell size and all units stated in config.
# =============================================================================

# ---- default configuration --------------------------------------------------

default_config <- function() {
  list(
    ## grid ---------------------------------------------------------------
    n_lat        = 12,      # latitude bands (rows), north = higher index
    n_shelf      = 12,      # cross shelf bands (cols), 1 = inshore/shallow
    cell_km      = 50,      # cell edge length, km
    n_months     = 96,      # months simulated
    ## temperature (degrees Celsius) --------------------------------------
    T_mean       = 18.0,    # domain mean climatological temperature
    g_lat        = 0.50,    # cooling per latitude band, degC per band
    g_depth      = 0.22,    # cooling per cross shelf band, degC per band
    seas_amp     = 3.0,     # seasonal amplitude
    seas_phase   = 0.6,     # seasonal phase, radians
    warm_trend   = 0.05,    # chronic warming, degC per year
    sd_field     = 0.40,    # marginal SD of the spatiotemporal noise field
    field_rho_s  = 0.55,    # spatial AR(1) correlation between adjacent cells
    field_rho_t  = 0.50,    # temporal AR(1) correlation between months
    ## event schedule ------------------------------------------------------
    n_events     = 6,       # heatwave episodes per replicate
    ev_dur_min   = 1,       # minimum duration, months
    ev_dur_max   = 4,       # maximum duration, months
    ev_peak_mean = 2.2,     # mean peak anomaly, degC
    ev_peak_sd   = 0.7,
    ev_thresh    = 1.0,     # category 1 threshold, degC above climatology
    ## thermal niche -------------------------------------------------------
    T_opt        = 20.0,
    sigma_th     = 1.9,
    T_crit       = 21.0,    # stress threshold for the thermal dose
    ## movement ------------------------------------------------------------
    rho_move     = 0.55,    # partial adjustment rate toward the ideal free target
    out_open     = TRUE,    # allow biomass to leave the sampling frame
    refuge_static = FALSE,  # TRUE: refuge suitability fixed in time (legacy); FALSE: refuge
                            #       is one band beyond the poleward edge and experiences the
                            #       same season, trend and event anomaly as the frame
    ## abundance state -----------------------------------------------------
    N0           = 1.0,
    r0           = 0.09,   # intrinsic monthly recovery rate
    K_cap        = 1.0,     # carrying capacity on the N scale
    M_nat        = 0.0,     # natural mortality beyond density dependence
    kappa        = 0.12,   # abundance sensitivity to thermal dose
    dose_lag     = 0,       # months between dose and abundance response
    F_scale      = 1.2e-4,  # fishing mortality per unit expected catch rate (q.pi.E.S), month^-1
    sd_proc      = 0.020,   # process noise SD on log N
    ## catchability --------------------------------------------------------
    q0           = 1.0,
    delta_dep    = 0.07,    # uniform catchability loss per unit event severity
    delta_sp     = 0.16,    # warm edge catchability loss
    gamma_sp     = 1.40,    # warm edge concentration
    ## observation ---------------------------------------------------------
    a0           = 0.9,     # encounter scale
    sigma_log    = 0.6,     # lognormal SD of positive catch
    ## fleet ---------------------------------------------------------------
    E_total      = 60,      # total effort units per month
    n_vessel     = 8,
    sd_vessel    = 0.25,    # vessel level log efficiency SD
    beta_mu      = 2.0,     # preference for high lagged catch rate
    beta_J       = 0.8,     # event driven contraction: proportional increase in the
                            #   fleet's targeting strength per unit severity, which
                            #   concentrates effort on fewer, better known cells
    beta_dist    = 0.15,    # travel cost against distance from port (inshore, south)
    belief_lr    = 0.5,     # weight on last month's observed catch rate in the fleet's belief
    ## survey --------------------------------------------------------------
    q_surv       = 1.0,
    a_surv       = 3.0,     # survey encounter scale
    cv_surv      = 0.45,    # survey lognormal CV on positive tows
    surv_cover   = 0.10,    # fraction of cells sampled per survey month
    surv_every   = 1,       # survey every k months (1 monthly, 3 quarterly, 6 semiannual)
    surv_design  = "srs"    # "srs" or "strat" (equal allocation across latitude bands)
  )
}

# ---- switch sets ------------------------------------------------------------
# The core design is a 2^3 factorial over abundance decline, niche tracking
# (redistribution) and uniform catchability loss, plus a sampling only cell and
# two spatially structured catchability cells (P9 / 4.1).

scenario_switches <- function(label) {
  # Core design: 2^4 factorial over abundance decline (A), redistribution (B),
  # uniform catchability loss (C) and the fleet's event driven contraction (D),
  # labelled S<A><B><C><D>. Warm edge catchability cells SPQ* carry the spatial
  # term; their last digit is the sampling switch.
  z <- list(decline = FALSE, tracking = FALSE, catch_uniform = FALSE,
            catch_spatial = FALSE, sampling = FALSE)
  if (grepl("^S[01]{4}$", label)) {
    b <- as.integer(strsplit(substr(label, 2, 5), "")[[1]]) == 1L
    return(modifyList(z, list(decline = b[1], tracking = b[2],
                              catch_uniform = b[3], sampling = b[4])))
  }
  switch(label,
    "SPQ0"  = modifyList(z, list(catch_spatial = TRUE)),
    "SPQ1"  = modifyList(z, list(catch_spatial = TRUE, sampling = TRUE)),
    "SPQT0" = modifyList(z, list(catch_spatial = TRUE, tracking = TRUE)),
    "SPQT1" = modifyList(z, list(catch_spatial = TRUE, tracking = TRUE,
                                 sampling = TRUE)),
    stop("unknown scenario label: ", label))
}

scenario_levels <- as.vector(outer(outer(outer(c("S0", "S1"), c("0", "1"), paste0),
                                   c("0", "1"), paste0), c("0", "1"), paste0))
scenario_levels <- sort(scenario_levels)      # S0000 ... S1111
scenario_warm   <- c("SPQ0", "SPQ1", "SPQT0", "SPQT1")

# ---- grid -------------------------------------------------------------------

make_grid <- function(cfg) {
  g <- expand.grid(shelf = seq_len(cfg$n_shelf), lat = seq_len(cfg$n_lat))
  g$cell <- seq_len(nrow(g))
  g$y    <- as.numeric(scale(g$lat))     # standardised northing
  g$x    <- as.numeric(scale(g$shelf))   # standardised offshore position
  # distance from port: ports sit inshore at the equatorward end
  g$dist <- sqrt((g$lat - 1)^2 + (g$shelf - 1)^2) / cfg$n_lat
  g
}

# ---- event schedule ---------------------------------------------------------
# Returns a non negative monthly severity series J_t (degC above the category 1
# threshold, category weighted) plus the event bookkeeping the review asks to
# be reported (number, duration, cumulative intensity).

draw_events <- function(cfg) {
  Tn <- cfg$n_months
  H  <- numeric(Tn)
  starts <- sort(sample.int(Tn - cfg$ev_dur_max, cfg$n_events))
  durs   <- sample(cfg$ev_dur_min:cfg$ev_dur_max, cfg$n_events, replace = TRUE)
  peaks  <- pmax(cfg$ev_thresh + 0.2,
                 rnorm(cfg$n_events, cfg$ev_peak_mean, cfg$ev_peak_sd))
  for (k in seq_len(cfg$n_events)) {
    idx <- starts[k]:min(Tn, starts[k] + durs[k] - 1L)
    # triangular within event profile peaking mid event
    w <- if (length(idx) == 1) 1 else {
      u <- seq(0, 1, length.out = length(idx))
      0.55 + 0.45 * sin(pi * u)
    }
    H[idx] <- pmax(H[idx], peaks[k] * w)
  }
  # category weighted severity: multiples of the category 1 threshold
  J <- pmax(H - cfg$ev_thresh, 0) / cfg$ev_thresh
  on <- J > 0
  n_real <- sum(on & !c(FALSE, on[-Tn]))     # separate runs after any merging
  list(H = H, J = J,
       n_events   = n_real,
       n_drawn    = cfg$n_events,
       mean_dur   = mean(durs),
       cum_int    = sum(J),
       event_month = J > 0)
}

# ---- correlated temperature noise field -------------------------------------

sim_field <- function(cfg, grid) {
  S <- nrow(grid); Tn <- cfg$n_months
  # separable AR(1) in space (nearest neighbour smoothing) and time
  raw <- matrix(rnorm(S * Tn), S, Tn)
  # spatial smoothing by a small kernel on the lattice
  lat <- grid$lat; shf <- grid$shelf
  W <- matrix(0, S, S)
  for (i in seq_len(S)) {
    d <- abs(lat - lat[i]) + abs(shf - shf[i])
    W[i, ] <- cfg$field_rho_s^d
  }
  W <- W / sqrt(rowSums(W^2))
  sm <- W %*% raw
  # temporal AR(1)
  out <- sm
  if (Tn > 1) for (t in 2:Tn)
    out[, t] <- cfg$field_rho_t * out[, t - 1] +
      sqrt(1 - cfg$field_rho_t^2) * sm[, t]
  out <- out / sd(as.numeric(out)) * cfg$sd_field
  out
}

# ---- temperature ------------------------------------------------------------
# T = climatology + chronic trend + event anomaly + noise.  Keeping the three
# deterministic pieces separate is what makes the paired counterfactual
# (event on / event off) well defined (P7, 3.5).

temperature_parts <- function(cfg, grid, events, field) {
  S <- nrow(grid); Tn <- cfg$n_months
  tt <- seq_len(Tn)
  clim_s <- cfg$T_mean - cfg$g_lat * (grid$lat - mean(grid$lat)) -
    cfg$g_depth * (grid$shelf - mean(grid$shelf))
  seas   <- cfg$seas_amp * sin(2 * pi * tt / 12 - cfg$seas_phase)
  trend  <- cfg$warm_trend * tt / 12
  base   <- outer(clim_s, seas + trend, "+") + field
  # climatology of the refuge: one latitude band beyond the poleward edge, mid shelf
  clim_ref <- cfg$T_mean - cfg$g_lat * ((cfg$n_lat + 1) - mean(grid$lat))
  list(base = base,                                    # no event
       event = matrix(events$H, S, Tn, byrow = TRUE),  # event anomaly
       refuge_base = clim_ref + seas + trend,          # refuge, no event
       refuge_event = events$H)                        # refuge anomaly
}

#' Refuge suitability series. Dynamic (default): the refuge is one band beyond
#' the poleward edge and experiences the same season, trend and event anomaly as
#' the frame. Static (legacy): fixed near optimal suitability, used only as a
#' sensitivity comparison, because a refuge that never warms acts as a
#' thermostat that draws biomass out of the frame regardless of conditions.
refuge_suitability <- function(cfg, parts, with_event = TRUE) {
  if (!cfg$out_open) return(rep(0, length(parts$refuge_base)))
  if (cfg$refuge_static) {
    v <- exp(-((cfg$T_opt - cfg$g_lat * 1.5) - cfg$T_opt)^2 / (2 * cfg$sigma_th^2))
    return(rep(v, length(parts$refuge_base)))
  }
  Tref <- parts$refuge_base + if (with_event) parts$refuge_event else 0
  exp(-(Tref - cfg$T_opt)^2 / (2 * cfg$sigma_th^2))
}

# ---- niche, movement, outside frame ----------------------------------------

suitability <- function(Tmat, cfg) {
  exp(-(Tmat - cfg$T_opt)^2 / (2 * cfg$sigma_th^2))
}

#' Distribution shares with partial adjustment and an outside frame state.
#' Returns a list with pi (S x T, summing to 1 - p_out each month) and p_out.
distribution_shares <- function(phi, cfg, tracking, phi_out_t,
                                phi_base = phi, phi_out_base = phi_out_t) {
  S <- nrow(phi); Tn <- ncol(phi)
  if (length(phi_out_t) == 1) phi_out_t <- rep(phi_out_t, Tn)
  if (length(phi_out_base) == 1) phi_out_base <- rep(phi_out_base, Tn)

  # the outside frame state represents one latitude band beyond the poleward
  # edge, so it carries the area weight of a band rather than of a single cell
  w_out <- cfg$n_shelf

  if (!tracking) {
    # fixed shares: the mean suitability over the series WITHOUT heatwaves
    # (phi_base, phi_out_base), so that the factual and counterfactual runs use
    # identical shares and no future event can enter the paired contrast
    base <- rowMeans(phi_base)
    tot  <- sum(base) + w_out * mean(phi_out_base)
    pi_s <- matrix(base / tot, S, Tn)
    p_out <- rep(1 - sum(base) / tot, Tn)
    return(list(pi = pi_s, p_out = p_out))
  }

  target <- rbind(phi, matrix(w_out * phi_out_t, 1, Tn, byrow = TRUE))
  target <- sweep(target, 2, colSums(target), "/")
  cur <- target[, 1]
  out <- matrix(0, S + 1, Tn)
  out[, 1] <- cur
  rho <- cfg$rho_move
  for (t in 2:Tn) {
    cur <- (1 - rho) * cur + rho * target[, t]
    cur <- cur / sum(cur)
    out[, t] <- cur
  }
  list(pi = out[seq_len(S), , drop = FALSE], p_out = out[S + 1, ])
}

# ---- abundance state process ------------------------------------------------
#' One month of the state transition on the log scale,
#'   log N_t = log N_{t-1} + r0 (1 - N_{t-1}/K) - M - F_t - kappa d_{t-L} + eta_t,
#' where F_t >= 0 is fishing mortality (see fishing_mortality) and d is the
#' occupancy weighted thermal dose. Called month by month from the sequential
#' simulation so that effort, removals and observations are mutually consistent.
abundance_step <- function(logN_prev, F_t, dose_lagged, eta_t, cfg, decline) {
  Nprev <- exp(logN_prev)
  logN_prev + cfg$r0 * (1 - Nprev / cfg$K_cap) - cfg$M_nat - F_t -
    (if (decline) cfg$kappa * dose_lagged else 0) + eta_t
}

#' Occupancy weighted thermal dose, d_t = sum_s pi_st max(T_st - T_crit, 0).
thermal_dose <- function(Tmat, pi_s, cfg) colSums(pi_s * pmax(Tmat - cfg$T_crit, 0))

#' Fishing mortality in month t: proportional to the expected catch per unit
#' abundance, S * sum_s q_st pi_st E_st, so it is non negative, rises with
#' fleet stock overlap and falls when catchability falls. Removal based, so
#' it is the same functional in the factual and counterfactual worlds.
fishing_mortality <- function(q_t, pi_t, E_t, cfg, S) {
  cfg$F_scale * S * sum(q_t * pi_t * E_t)
}

#' Severity scaled to unit standard deviation over event months; identically
#' zero when there are no events (the counterfactual world).
scale_severity <- function(Jt) {
  if (!any(Jt > 0)) return(rep(0, length(Jt)))
  Jt / max(sd(Jt[Jt > 0]), 1e-8)
}

# ---- catchability -----------------------------------------------------------

catchability_field <- function(cfg, grid, Jt, switches) {
  S <- nrow(grid); Tn <- length(Jt)
  Js <- scale_severity(Jt)                # non negative, event scaled severity
  q <- matrix(cfg$q0, S, Tn)
  if (switches$catch_uniform)
    q <- q * matrix(exp(-cfg$delta_dep * Js), S, Tn, byrow = TRUE)
  if (switches$catch_spatial) {
    h  <- pmax(1 - cfg$gamma_sp * grid$y, 0)
    gs <- h / mean(h)                    # domain mean one, as stated (3.10)
    q  <- q * exp(-cfg$delta_sp * outer(gs, Js, "*"))
  }
  q
}
