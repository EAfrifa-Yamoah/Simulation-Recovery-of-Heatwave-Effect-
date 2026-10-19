# =============================================================================
# scripts/make_figures.R
# All manuscript figures from results/, in base R. Rendered to PDF (vector)
# and rasterised to 300 dpi PNG with pdftoppm. Captions live in the manuscript,
# not in the figures; panels carry letters only. One meaning per colour across
# the whole suite (Okabe Ito palette); one shared legend per figure.
#   Rscript scripts/make_figures.R
# =============================================================================
Sys.setlocale("LC_CTYPE", "C.UTF-8")   # so the degree sign in labels is read as one character
for (f in list.files("R", full.names = TRUE)) source(f)
dir.create("figures", showWarnings = FALSE)

COL <- c(truth = "#000000", truth_in = "#7F7F7F", nominal = "#D55E00", naive = "#E69F00",
         delta = "#0072B2", delta_fished = "#56B4E9", deltaint = "#CC79A7",
         survey = "#009E73", event = "#F0E442")
LAB <- c(truth = "Truth (total abundance)", truth_in = "Truth (in-frame abundance)",
         nominal = "Nominal CPUE", naive = "Positive-only", delta = "Delta index",
         delta_fished = "Delta index, fished cells", deltaint = "Spatial delta index",
         survey = "Survey (design-based)")
## colours used only for elements that are not estimators, so that no estimator
## colour is reused with another meaning
COL_NEUTRAL <- "#5F5F5F"; COL_LAG <- "#CC79A7"
alpha <- function(col, a) { v <- col2rgb(col) / 255; rgb(v[1], v[2], v[3], a) }

fig_open <- function(nm, w, h, oma_top = 0) {
  cairo_pdf(file.path("figures", paste0(nm, ".pdf")), width = w, height = h, pointsize = 10)
  par(family = "sans", mgp = c(2.1, 0.6, 0), tcl = -0.3, las = 1, cex.axis = 0.85,
      oma = c(0, 0, oma_top, 0))
}
fig_close <- function(nm) {
  dev.off()
  system2("pdftoppm", c("-png", "-r", "300", "-singlefile",
                        shQuote(file.path("figures", paste0(nm, ".pdf"))),
                        shQuote(file.path("figures", nm))))
}
panel <- function(s, x = -0.02, y = 1.04) {
  u <- par("usr")
  px <- u[1] + x * diff(u[1:2]); py <- u[4] + (y - 1) * diff(u[3:4])
  if (par("xlog")) px <- 10^px; if (par("ylog")) py <- 10^py   # usr is log10 on log axes
  text(px, py, s, font = 2, cex = 1.15, xpd = NA, adj = c(1, 0))
}
shade_events <- function(month, event) {
  on <- as.logical(event); r <- rle(on); e <- cumsum(r$lengths); s <- e - r$lengths + 1
  for (k in which(r$values)) rect(month[s[k]] - 0.5, par("usr")[3], month[e[k]] + 0.5,
                                  par("usr")[4], col = alpha(COL["event"], 0.35), border = NA)
}
legend_top <- function(keys, lty = NULL, pch = NULL, lwd = 2, ncol = 4, inset = 0) {
  legend("top", legend = LAB[keys], col = COL[keys], lty = if (is.null(lty)) 1 else lty,
         pch = pch, lwd = lwd, bty = "n", horiz = FALSE, ncol = ncol, cex = 0.85,
         inset = inset, xpd = NA)
}
violin <- function(x, at, col, width = 0.38) {
  x <- x[is.finite(x)]; if (length(x) < 5) return(invisible())
  d <- density(x, n = 128); yy <- d$x; xx <- d$y / max(d$y) * width
  polygon(c(at + xx, rev(at - xx)), c(yy, rev(yy)), col = alpha(col, 0.55), border = col)
  segments(at - width * 0.7, mean(x), at + width * 0.7, mean(x), col = col, lwd = 2.5)
}

cfg  <- default_config(); grid <- make_grid(cfg)
## centre-of-gravity shifts are in units of the standardised northing; convert to km
KM <- cfg$cell_km * sd(grid$lat)
S <- nrow(grid); Tn <- cfg$n_months
fac <- readRDS("results/factorial.rds")
fs  <- read.csv("results/factorial_summary.csv")

## ---- example worlds for Figures 1 to 3 ---------------------------------------
example_world <- function(sc, seed) {
  set.seed(seed); sw <- scenario_switches(sc)
  ev <- draw_events(cfg); field <- sim_field(cfg, grid)
  eta <- rnorm(Tn, 0, cfg$sd_proc); fleet <- draw_fleet(cfg, grid)
  parts <- temperature_parts(cfg, grid, ev, field)
  w <- run_world(cfg, grid, sw, parts$base + parts$event, refuge_suitability(cfg, parts, TRUE),
                 ev$J, eta, fleet)
  wc <- run_world(cfg, grid, sw, parts$base, refuge_suitability(cfg, parts, FALSE),
                  rep(0, Tn), eta, fleet, E_fixed = w$E, observe = FALSE)
  I <- (ev$J - mean(ev$J)) / sd(ev$J)
  list(w = w, wc = wc, ev = ev, I = I, T = parts$base + parts$event, sw = sw)
}
to_grid <- function(v) { M <- matrix(NA, cfg$n_lat, cfg$n_shelf)
  M[cbind(grid$lat, grid$shelf)] <- v; t(M) }    # image() wants x = shelf, y = lat

# ---- Figure 1: spatial fields at a heatwave peak ------------------------------
fig1 <- function() {
  ex <- example_world("SPQT1", 20260724 + 7)      # all channels incl. warm edge catchability
  pk <- which.max(ex$ev$J)
  fig_open("fig1_spatial_fields", 9.2, 3.1)
  layout(matrix(1:6, 1, 6), widths = c(1, 0.18, 1, 0.18, 1, 0.18))
  par(mar = c(3.2, 3.2, 1.8, 0.4))
  draw <- function(v, pal, lab, letter, diverge = FALSE) {
    z <- to_grid(v)
    rng <- if (diverge) c(-1, 1) * max(abs(v)) else range(v)
    image(1:cfg$n_shelf, 1:cfg$n_lat, z, col = pal, zlim = rng, xlab = "cross shelf band",
          ylab = if (letter == "a") "latitude band (north up)" else "", axes = FALSE)
    box(); panel(letter)
    par(mar = c(3.2, 0.3, 1.8, 2.6))
    image(1, seq(rng[1], rng[2], length.out = 64), matrix(seq(rng[1], rng[2], length.out = 64), 1),
          col = pal, axes = FALSE, xlab = "", ylab = ""); axis(4, cex.axis = 0.75)
    mtext(lab, 4, line = 1.7, cex = 0.65, las = 0)
    par(mar = c(3.2, 3.2, 1.8, 0.4))
  }
  draw(ex$T[, pk], hcl.colors(64, "Inferno"), "temperature (°C)", "a")
  draw(ex$w$pi[, pk] * S, hcl.colors(64, "Viridis"), "stock density", "b")
  draw(ex$w$q[, pk] - cfg$q0, hcl.colors(64, "Blue-Red 3"), "catchability anomaly", "c", TRUE)
  fig_close("fig1_spatial_fields")
}

# ---- Figure 2: simulated observations ----------------------------------------
fig2 <- function() {
  ex <- example_world("S1111", 20260724 + 7)
  fig_open("fig2_data_distributions", 8.6, 5.8)
  par(mfrow = c(2, 2), mar = c(3.4, 3.6, 1.8, 0.8))
  tm <- colMeans(ex$T)
  plot(1:Tn, tm, type = "n", xlab = "month", ylab = "mean temperature (°C)")
  shade_events(1:Tn, ex$ev$J > 0); lines(1:Tn, tm, lwd = 1.4); panel("a")
  tt <- seq(min(ex$T) - 1, max(ex$T) + 1, length.out = 200)
  h <- hist(as.numeric(ex$T), breaks = 40, plot = FALSE)
  plot(tt, exp(-(tt - cfg$T_opt)^2 / (2 * cfg$sigma_th^2)), type = "n", xlab = "temperature (°C)",
       ylab = "thermal suitability")
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$density / max(h$density),
       col = alpha(COL["truth_in"], 0.35), border = NA)
  lines(tt, exp(-(tt - cfg$T_opt)^2 / (2 * cfg$sigma_th^2)), col = "black", lwd = 2.2)
  abline(v = cfg$T_opt, lty = 2, col = "grey40"); panel("b")
  y <- ex$w$Y[!is.na(ex$w$Y)]; fz <- mean(y == 0)
  hist(log(y[y > 0]), breaks = 40, col = alpha(COL_NEUTRAL, 0.75), border = "white",
       main = "", xlab = "log positive CPUE", ylab = "count")
  legend("topright", bty = "n", cex = 0.85, legend = sprintf("zeros on fished cells: %.0f%%", 100 * fz))
  panel("c")
  ff <- colMeans(ex$w$E > 0)
  plot(1:Tn, ff, type = "n", xlab = "month", ylab = "fraction of cells fished")
  shade_events(1:Tn, ex$ev$J > 0); lines(1:Tn, ff, lwd = 1.2); panel("d")
  fig_close("fig2_data_distributions")
}

# ---- Figure 3: index trajectories against both truths ------------------------
fig3 <- function() {
  scen <- c("S1001", "S0101", "S0011", "S1111")
  nm <- c("Abundance decline", "Redistribution", "Uniform catchability", "All three")
  fig_open("fig3_trajectories", 9, 6.8, oma_top = 2.4)
  par(mfrow = c(2, 2), mar = c(3.2, 3.6, 1.5, 0.6))
  for (k in seq_along(scen)) {
    ex <- example_world(scen[k], 20260724 + 50 + k)
    sim <- list(fish = list(Y = ex$w$Y, E = ex$w$E), grid = grid, J = ex$ev$J, I = ex$I)
    cen <- function(v) { v[is.finite(v)] <- v[is.finite(v)] - mean(v[is.finite(v)]); v }
    tr  <- cen(ex$w$logN - ex$wc$logN); tri <- cen(ex$w$logN_in - ex$wc$logN_in)
    nom <- cen(est_nominal(sim, ex$I)$index); del <- cen(est_delta_fe(sim, ex$I, grid)$index)
    sv  <- simulate_survey(cfg, grid, exp(ex$w$logN), ex$w$pi, coverage = 0.20)
    su  <- cen(est_survey(sv, ex$I, grid)$index)
    yl <- range(c(tr, tri, nom, del, su), na.rm = TRUE)
    plot(1:Tn, tr, type = "n", ylim = yl, xlab = if (k > 2) "month" else "",
         ylab = if (k %% 2 == 1) "log index (centred)" else "")
    shade_events(1:Tn, ex$ev$J > 0)
    lines(1:Tn, nom, col = COL["nominal"], lwd = 1.1); lines(1:Tn, del, col = COL["delta"], lwd = 1.1)
    lines(1:Tn, su, col = COL["survey"], lwd = 1.1)
    lines(1:Tn, tri, col = COL["truth_in"], lwd = 2, lty = 2); lines(1:Tn, tr, col = COL["truth"], lwd = 2.2)
    abline(h = 0, col = "grey70", lwd = 0.6)
    panel(letters[k])
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0.2, 0), new = TRUE)
  plot.new()
  legend("top", legend = LAB[c("truth", "truth_in", "nominal", "delta", "survey")],
         col = COL[c("truth", "truth_in", "nominal", "delta", "survey")],
         lty = c(1, 2, 1, 1, 1), lwd = c(2.2, 2, 1.1, 1.1, 1.1), bty = "n", ncol = 5, cex = 0.82, xpd = NA)
  fig_close("fig3_trajectories")
}

# ---- Figure 4: sampling distributions -----------------------------------------
fig4 <- function() {
  scen <- c("S1001", "S0101", "S0011", "S1111"); lab <- c("Decline", "Redistribute", "Catchability", "Mixed")
  fig_open("fig4_sampling_distributions", 9.4, 4.1, oma_top = 2.2)
  par(mfrow = c(1, 3), mar = c(4.6, 3.8, 2.6, 0.6))
  get <- function(sc, k) { v <- fac[[k]][fac$scenario == sc]; v[is.finite(v)] }
  pl <- function(key, ylab, letter, truth_key) {
    vals <- lapply(scen, get, paste0("beta_", key)); rng <- range(unlist(vals), fac$beta_cf[fac$scenario %in% scen], na.rm = TRUE)
    plot(NA, xlim = c(0.4, 4.6), ylim = rng, xaxt = "n", xlab = "", ylab = ylab)
    axis(1, at = 1:4, labels = lab, las = 2, cex.axis = 0.85)
    abline(h = 0, lty = 3, col = "grey50")
    for (i in 1:4) violin(vals[[i]], i, COL[key])
    for (i in 1:4) { points(i, mean(get(scen[i], "beta_cf")), pch = 4, lwd = 2.5, cex = 1.3)
      if (truth_key) points(i, mean(get(scen[i], "beta_cf_in")), pch = 1, lwd = 2, cex = 1.3, col = COL["truth_in"]) }
    panel(letter)
  }
  pl("nominal", "nominal CPUE coefficient", "a", TRUE)
  pl("survey", "survey index coefficient", "b", TRUE)
  vals <- lapply(lapply(scen, get, "dG_fishery"), function(v) v * KM)
  plot(NA, xlim = c(0.4, 4.6), ylim = range(unlist(vals), fac$dG_true * KM, na.rm = TRUE), xaxt = "n",
       xlab = "", ylab = "fishery centre of gravity shift (km)")
  axis(1, at = 1:4, labels = lab, las = 2, cex.axis = 0.85); abline(h = 0, lty = 3, col = "grey50")
  for (i in 1:4) violin(vals[[i]], i, COL["nominal"])
  for (i in 1:4) segments(i - 0.3, KM * mean(get(scen[i], "dG_true")), i + 0.3, KM * mean(get(scen[i], "dG_true")), lwd = 3)
  panel("c")
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0.2, 0), new = TRUE); plot.new()
  legend("top", legend = c("truth, total abundance", "truth, in-frame abundance", "true centre of gravity"),
         pch = c(4, 1, NA), lty = c(NA, NA, 1), lwd = c(2.5, 2, 3), col = c("black", COL["truth_in"], "black"),
         bty = "n", ncol = 3, cex = 0.85, xpd = NA)
  fig_close("fig4_sampling_distributions")
}

# ---- Figure 5: effect by case, 2^3 biological cells with the fleet contraction on
fig5 <- function() {
  scen <- c("S0001", "S1001", "S0101", "S0011", "S1101", "S1011", "S0111", "S1111")
  lab <- c("None", "Abund.", "Redist.", "Catch.", "Ab+Re", "Ab+Ca", "Re+Ca", "All 3")
  keys <- c("nominal", "delta", "delta_fished", "survey")
  cols <- c("beta_nominal", "beta_delta", "beta_delta_fished", "beta_survey")
  fig_open("fig5_effect_by_case", 9, 4.6)
  par(mar = c(3.4, 4.4, 1.2, 0.6))
  st <- lapply(cols, function(cc) t(sapply(scen, function(sc) { v <- fac[[cc]][fac$scenario == sc]; v <- v[is.finite(v)]
    c(mean(v), quantile(v, c(.025, .975))) })))
  rng <- range(unlist(st), fac$beta_cf[fac$scenario %in% scen], na.rm = TRUE)
  x <- seq_along(scen); w <- 0.17
  plot(NA, xlim = c(0.5, 8.5), ylim = rng, xaxt = "n", xlab = "", ylab = "heatwave coefficient")
  axis(1, at = x, labels = lab); abline(h = 0, lwd = 0.8)
  for (j in seq_along(keys)) {
    off <- (j - (length(keys) + 1) / 2) * w
    rect(x + off - w / 2, 0, x + off + w / 2, st[[j]][, 1], col = COL[keys[j]], border = NA)
    arrows(x + off, st[[j]][, 2], x + off, st[[j]][, 3], angle = 90, code = 3, length = 0.02, col = "grey30", lwd = 0.7)
  }
  tt <- sapply(scen, function(sc) mean(fac$beta_cf[fac$scenario == sc]))
  ti <- sapply(scen, function(sc) mean(fac$beta_cf_in[fac$scenario == sc]))
  segments(x - 0.42, tt, x + 0.42, tt, lwd = 2.6, col = COL["truth"])
  segments(x - 0.42, ti, x + 0.42, ti, lwd = 2, lty = 2, col = COL["truth_in"])
  legend("bottomleft", legend = c(LAB[keys], LAB["truth"], LAB["truth_in"]),
         fill = c(COL[keys], NA, NA), border = NA, lty = c(NA, NA, NA, NA, 1, 2),
         lwd = c(NA, NA, NA, NA, 2.6, 2), col = c(rep(NA, 4), COL["truth"], COL["truth_in"]),
         bty = "n", ncol = 2, cex = 0.8)
  fig_close("fig5_effect_by_case")
}

# ---- Figure 6: survey design curve (coverage x frequency) + power -------------
fig6 <- function() {
  cs <- read.csv("results/coverage_summary.csv")
  d0 <- cs[cs$cv_surv == cfg$cv_surv & cs$design == "srs", ]
  fig_open("fig6_survey_coverage", 9.2, 4.8, oma_top = 3)
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 2.6, 0.6))
  ltys <- c(`1` = 1, `3` = 2, `6` = 3)
  plot(NA, xlim = range(d0$tows_per_year), ylim = range(c(d0$rmse, d0$lo, d0$hi, d0$rmse_nominal, d0$rmse_delta), na.rm = TRUE),
       log = "xy", xlab = "survey tows per year", ylab = "error of the survey index (RMSE)")
  for (ev in c(1, 3, 6)) { d <- d0[d0$every == ev, ]; d <- d[order(d$tows_per_year), ]
    polygon(c(d$tows_per_year, rev(d$tows_per_year)), c(d$lo, rev(d$hi)), col = alpha(COL["survey"], 0.12), border = NA)
    lines(d$tows_per_year, d$rmse, col = COL["survey"], lwd = 2, lty = ltys[as.character(ev)], type = "b", pch = 19, cex = 0.8) }
  abline(h = mean(d0$rmse_nominal), col = COL["nominal"], lwd = 1.8, lty = 5)
  abline(h = mean(d0$rmse_delta), col = COL["delta"], lwd = 1.8, lty = 5)
  panel("a")
  plot(NA, xlim = range(d0$tows_per_year), ylim = c(0, 1), log = "x", xlab = "survey tows per year",
       ylab = "power: 95% CI excludes zero")
  for (ev in c(1, 3, 6)) { d <- d0[d0$every == ev, ]; d <- d[order(d$tows_per_year), ]
    lines(d$tows_per_year, d$power_survey, col = COL["survey"], lwd = 2, lty = ltys[as.character(ev)])
    points(d$tows_per_year, d$power_survey, col = COL["survey"], pch = ifelse(d$cover_survey >= 0.8, 19, 1), cex = 0.9, lwd = 1.5) }
  abline(h = mean(d0$power_nominal, na.rm = TRUE), col = COL["nominal"], lwd = 1.8, lty = 5)
  abline(h = mean(d0$power_delta, na.rm = TRUE), col = COL["delta"], lwd = 1.8, lty = 5)
  panel("b")
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0.2, 0), new = TRUE); plot.new()
  legend("top", legend = c("survey, monthly", "survey, quarterly", "survey, twice-yearly",
                           "nominal CPUE (no survey)", "delta index (no survey)", "open point: CI coverage below 0.8"),
         col = c(rep(COL["survey"], 3), COL["nominal"], COL["delta"], COL["survey"]), lty = c(1, 2, 3, 5, 5, NA),
         pch = c(19, 19, 19, NA, NA, 1), lwd = c(2, 2, 2, 1.8, 1.8, NA), bty = "n", ncol = 3, cex = 0.8, xpd = NA)
  fig_close("fig6_survey_coverage")
}

# ---- Figure 7: temporal recovery from the observed survey index ---------------
fig7 <- function() {
  tp <- readRDS("results/temporal.rds"); ex <- tp$example
  fig_open("fig7_temporal", 9.6, 4.3, oma_top = 2.2)
  layout(matrix(1:2, 1, 2), widths = c(2, 1)); par(mar = c(3.6, 3.8, 2.4, 0.6))
  yl <- range(c(ex$truth, ex$fit_acute_obs, ex$fit_dl_obs), na.rm = TRUE)
  plot(ex$month, ex$truth, type = "n", ylim = yl, xlab = "month", ylab = "heatwave effect on log abundance (centred)")
  shade_events(ex$month, ex$event); abline(h = 0, col = "grey70", lwd = 0.6)
  lines(ex$month, ex$fit_acute_obs, col = COL_NEUTRAL, lwd = 1.5, lty = 2)
  lines(ex$month, ex$fit_dl_obs, col = COL_LAG, lwd = 1.8)
  lines(ex$month, ex$truth, col = COL["truth"], lwd = 2.3)
  ca <- cor(ex$fit_acute_obs, ex$truth); cd <- cor(ex$fit_dl_obs, ex$truth)
  legend("bottomleft", bty = "n", cex = 0.8, legend = sprintf("correlation with truth:  same-month %.2f,  distributed lag %.2f", ca, cd))
  panel("a")
  w <- colMeans(tp$lag_weights_obs); q <- apply(tp$lag_weights_obs, 2, quantile, c(.25, .75))
  plot(0:tp$L, w, type = "n", ylim = range(q, 0), xlab = "lag (months)", ylab = "estimated weight on severity")
  abline(h = 0, col = "grey60", lwd = 0.6)
  polygon(c(0:tp$L, rev(0:tp$L)), c(q[1, ], rev(q[2, ])), col = alpha(COL_LAG, 0.25), border = NA)
  lines(0:tp$L, w, col = COL_LAG, lwd = 2); panel("b")
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0.2, 0), new = TRUE); plot.new()
  legend("top", legend = c("true heatwave effect", "same-month model", "distributed lag model"),
         col = c(COL["truth"], COL_NEUTRAL, COL_LAG), lty = c(1, 2, 1), lwd = c(2.3, 1.5, 1.8), bty = "n", ncol = 3, cex = 0.82, xpd = NA)
  fig_close("fig7_temporal")
}

# ---- Figure 8: warm edge catchability stress test ------------------------------
fig8 <- function() {
  scen <- c("S0101", "SPQ1", "SPQT1"); lab <- c("Redistribute", "Warm edge\ncatch.", "Warm edge\n+ redist.")
  fig_open("fig8_spatial_catchability", 9.4, 4.5, oma_top = 2.4)
  par(mfrow = c(1, 2), mar = c(4.4, 3.8, 2.6, 0.6))
  mci <- function(sc, k) { v <- fac[[k]][fac$scenario == sc]; v <- v[is.finite(v)]
    if (grepl("^dG", k)) v <- v * KM; c(mean(v), quantile(v, c(.025, .975))) }
  tmean <- function(sc, k) { m <- mean(fac[[k]][fac$scenario == sc]); if (grepl("^dG", k)) m * KM else m }
  bars <- function(keys, cols, ylab, letter, truth_cols = NULL) {
    st <- lapply(cols, function(cc) t(sapply(scen, mci, cc)))
    rng <- range(unlist(st), if (!is.null(truth_cols)) sapply(scen, tmean, truth_cols), na.rm = TRUE)
    plot(NA, xlim = c(0.5, 3.5), ylim = rng, xaxt = "n", xlab = "", ylab = ylab)
    axis(1, at = 1:3, labels = lab, padj = 0.6, cex.axis = 0.85); abline(h = 0, lwd = 0.8)
    for (j in seq_along(keys)) { off <- (j - 1.5) * 0.36
      rect(1:3 + off - 0.18, 0, 1:3 + off + 0.18, st[[j]][, 1], col = COL[keys[j]], border = NA)
      arrows(1:3 + off, st[[j]][, 2], 1:3 + off, st[[j]][, 3], angle = 90, code = 3, length = 0.02, col = "grey30", lwd = 0.7) }
    if (!is.null(truth_cols)) { tt <- sapply(scen, tmean, truth_cols)
      segments(1:3 - 0.42, tt, 1:3 + 0.42, tt, lwd = 2.6) }
    panel(letter)
  }
  bars(c("nominal", "survey"), c("dG_fishery", "dG_survey"), "poleward shift in centre of gravity (km)", "a", "dG_true")
  bars(c("nominal", "survey"), c("beta_nominal", "beta_survey"), "heatwave coefficient", "b", "beta_cf_in")
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0.2, 0), new = TRUE); plot.new()
  legend("top", legend = c("fishery (centre of gravity / nominal CPUE)", "survey (centre of gravity / index)", "truth (centre of gravity / in-frame effect)"),
         fill = c(COL["nominal"], COL["survey"], NA), border = NA, lty = c(NA, NA, 1), lwd = c(NA, NA, 2.6),
         col = c(NA, NA, "black"), bty = "n", ncol = 3, cex = 0.8, xpd = NA)
  fig_close("fig8_spatial_catchability")
}

# ---- Figure 9: matched calibration, identifiability, sensitivity --------------
fig9 <- function() {
  mm <- read.csv("results/matched.csv"); id <- readRDS("results/identifiability.rds")
  sens <- read.csv("results/sensitivity_design.csv")
  fig_open("fig9_sensitivity_identifiability", 9.6, 8.2)
  layout(matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE), widths = c(1, 1, 1.15))
  par(mar = c(3.6, 3.8, 2.4, 3.6))
  ks <- sort(unique(mm$kappa)); ds <- sort(unique(mm$delta_dep))
  surf <- function(col) { M <- matrix(NA, length(ks), length(ds))
    for (r in seq_len(nrow(mm))) M[match(mm$kappa[r], ks), match(mm$delta_dep[r], ds)] <- mm[[col]][r]; M }
  heat <- function(M, letter, label) {
    rng <- c(-1, 1) * max(abs(M), na.rm = TRUE)
    image(ks, ds, M, col = hcl.colors(64, "Blue-Red 3"), zlim = rng, xlab = expression("decline rate " * kappa),
          ylab = expression("catchability reduction " * delta[dep])); box(); panel(letter)
  }
  heat(surf("beta_nominal") - surf("beta_truth"), "a", "attribution gap (nominal minus truth)")
  heat(surf("beta_nominal"), "b", "observed nominal coefficient")
  contour(ks, ds, surf("beta_nominal"), levels = seq(-0.18, -0.03, by = 0.03), add = TRUE, lwd = 1, labcex = 0.6)
  ## c: profile in betaN, fishery vs both
  plot(id$grid_betaN, id$profile_both, type = "l", col = COL["survey"], lwd = 2.2, xlab = expression("abundance response " * beta[N]),
       ylab = "profile log likelihood (max = 0)")
  lines(id$grid_betaN, id$profile_fishery, col = COL["nominal"], lwd = 2.6); abline(v = id$betaN, lty = 2, col = "grey50")
  legend("bottomleft", legend = c("fishery only (flat)", "fishery + survey"), col = COL[c("nominal", "survey")], lwd = 2, bty = "n", cex = 0.8)
  panel("c")
  ## d: profile in gamma from the fishery alone
  plot(id$grid_gamma, id$profile_gamma_fishery, type = "l", col = COL["nominal"], lwd = 2.6,
       xlab = expression("warm-edge catchability slope " * gamma), ylab = "profile log likelihood (max = 0)")
  abline(v = id$gamma, lty = 2, col = "grey50"); panel("d")
  ## e, f: sensitivity SRC
  pretty <- c(beta_mu = "Targeting strength", beta_J = "Fleet contraction strength", a0 = "Fishery encounter scale", sigma_log = "Catch variability",
              cv_surv = "Survey CV", surv_cover = "Survey coverage", a_surv = "Survey encounter scale", belief_lr = "Fleet learning rate",
              rho_move = "Movement rate", n_events = "Number of events", warm_trend = "Warming trend", kappa = "Abundance sensitivity",
              delta_dep = "Uniform catchability loss", sigma_th = "Niche breadth", F_scale = "Fishing mortality scale")
  srcbar <- function(resp, letter, label) {
    z <- src(sens, resp); z <- z[order(abs(z$src)), ]
    par(mar = c(3.6, 8.6, 2.4, 0.8))
    barplot(z$src, horiz = TRUE, names.arg = pretty[z$parameter], las = 1, cex.names = 0.72,
            col = COL_NEUTRAL, border = NA, xlab = "standardised regression coefficient")
    abline(v = 0, col = "grey40"); panel(letter, x = -0.45)
  }
  srcbar("gap_nominal", "e", "nominal attribution gap")
  srcbar("bias_survey_total", "f", "survey bias for total abundance")
  fig_close("fig9_sensitivity_identifiability")
}

# ---- Figure 10: frame geometry --------------------------------------------------
fig10 <- function() {
  fr <- readRDS("results/frame.rds")
  scen <- c("S0101", "S1101"); lab <- c("Redistribution", "Abundance + redistribution")
  frames <- c("closed", "open_static", "open_dynamic"); fl <- c("closed frame", "open, static refuge", "open, warming refuge")
  fig_open("fig10_frame_geometry", 9.2, 4.5, oma_top = 2.4)
  par(mfrow = c(1, 2), mar = c(4.2, 3.8, 2.6, 0.6))
  for (k in 1:2) {
    st <- sapply(frames, function(fm) { d <- fr[fr$frame == fm & fr$scenario == scen[k], ]
      c(nom = mean(d$beta_nominal, na.rm = TRUE), del = mean(d$beta_delta, na.rm = TRUE), sur = mean(d$beta_survey, na.rm = TRUE),
        tt = mean(d$beta_cf), ti = mean(d$beta_cf_in)) })
    rng <- range(st, 0); x <- 1:3; w <- 0.24
    plot(NA, xlim = c(0.5, 3.5), ylim = rng, xaxt = "n", xlab = "", ylab = "heatwave coefficient")
    axis(1, at = x, labels = fl, cex.axis = 0.8); abline(h = 0, lwd = 0.8)
    rect(x - 1.5 * w, 0, x - 0.5 * w, st["nom", ], col = COL["nominal"], border = NA)
    rect(x - 0.5 * w, 0, x + 0.5 * w, st["del", ], col = COL["delta"], border = NA)
    rect(x + 0.5 * w, 0, x + 1.5 * w, st["sur", ], col = COL["survey"], border = NA)
    segments(x - 0.45, st["tt", ], x + 0.45, st["tt", ], lwd = 2.6, col = COL["truth"])
    segments(x - 0.45, st["ti", ], x + 0.45, st["ti", ], lwd = 2, lty = 2, col = COL["truth_in"])
    panel(letters[k])
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0.2, 0), new = TRUE); plot.new()
  legend("top", legend = c(LAB[c("nominal", "delta", "survey")], LAB["truth"], LAB["truth_in"]),
         fill = c(COL[c("nominal", "delta", "survey")], NA, NA), border = NA, lty = c(NA, NA, NA, 1, 2),
         lwd = c(NA, NA, NA, 2.6, 2), col = c(NA, NA, NA, COL["truth"], COL["truth_in"]), bty = "n", ncol = 5, cex = 0.78, xpd = NA)
  fig_close("fig10_frame_geometry")
}

# ---- Figure 11: identifiability (profiles and the matched decline x catchability grid)
PRETTY <- c(beta_mu = "Targeting strength", beta_J = "Fleet contraction strength",
            a0 = "Fishery encounter scale", sigma_log = "Catch variability",
            cv_surv = "Survey CV", surv_cover = "Survey coverage", a_surv = "Survey encounter scale",
            belief_lr = "Fleet learning rate", rho_move = "Movement rate", n_events = "Number of events",
            warm_trend = "Warming trend", kappa = "Abundance sensitivity",
            delta_dep = "Uniform catchability loss", sigma_th = "Niche breadth",
            F_scale = "Fishing mortality scale")
colkey <- function(rng, pal, lab) {
  par(mar = c(3.8, 0.3, 2.4, 3.4))
  image(1, seq(rng[1], rng[2], length.out = 64), matrix(seq(rng[1], rng[2], length.out = 64), 1),
        col = pal, axes = FALSE, xlab = "", ylab = "")
  axis(4, cex.axis = 0.75); box(); mtext(lab, 4, line = 2.2, cex = 0.65, las = 0)
}
fig11 <- function() {
  mm <- read.csv("results/matched.csv"); id <- readRDS("results/identifiability.rds")
  ks <- sort(unique(mm$kappa)); ds <- sort(unique(mm$delta_dep))
  surf <- function(v) { M <- matrix(NA, length(ks), length(ds))
    for (r in seq_len(nrow(mm))) M[match(mm$kappa[r], ks), match(mm$delta_dep[r], ds)] <- v[r]; M }
  fig_open("fig11_identifiability", 9.6, 7.8)
  layout(matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 2, 4, byrow = TRUE), widths = c(1, 0.24, 1, 0.24))
  par(mar = c(3.8, 4.0, 2.4, 1.0))
  plot(id$grid_betaN, id$profile_both, type = "l", col = COL["survey"], lwd = 2.2,
       ylim = range(c(id$profile_both, id$profile_fishery)),
       xlab = expression("abundance response " * beta[N]), ylab = "profile log likelihood (maximum = 0)")
  lines(id$grid_betaN, id$profile_fishery, col = COL["nominal"], lwd = 2.6)
  abline(v = id$betaN, lty = 2, col = "grey50")
  legend("bottomleft", legend = c("fishery data only", "fishery and survey data"),
         col = COL[c("nominal", "survey")], lwd = 2, bty = "n", cex = 0.85)
  panel("a")
  plot(id$grid_gamma, id$profile_gamma_fishery, type = "l", col = COL["nominal"], lwd = 2.6,
       xlab = expression("warm-edge catchability slope " * gamma), ylab = "profile log likelihood (maximum = 0)")
  abline(v = id$gamma, lty = 2, col = "grey50"); panel("b")
  pal <- rev(gray.colors(64, start = 0.30, end = 0.97))   # light to dark; reversed below so larger magnitudes are darker
  heat <- function(M, letter, lab, iso = FALSE) {
    rng <- range(M, na.rm = TRUE)
    par(mar = c(3.8, 4.0, 2.4, 0.4))
    image(ks, ds, M, col = rev(pal), zlim = rng, xlab = expression("decline rate " * kappa),
          ylab = expression("uniform catchability loss " * delta[dep]))
    box(); panel(letter)
    if (iso) {
      contour(ks, ds, M, levels = seq(-0.18, -0.03, by = 0.03), add = TRUE, lwd = 0.9,
              labcex = 0.65, col = "black")
      pts <- data.frame(k = c(0, 0.06, 0.12), d = c(0.07, 0.035, 0))
      points(pts$k, pts$d, pch = 21, bg = "white", cex = 1.6, lwd = 1.6)
      text(pts$k, pts$d, labels = 1:3, cex = 0.62, font = 2)
    }
    colkey(rng, rev(pal), lab)
  }
  heat(surf(mm$beta_nominal - mm$beta_truth), "c", "nominal coefficient minus truth")
  heat(surf(mm$beta_nominal), "d", "nominal coefficient", iso = TRUE)
  fig_close("fig11_identifiability")
}

# ---- Figure 12: global sensitivity, standardised regression coefficients -------
fig12 <- function() {
  sens <- read.csv("results/sensitivity_design.csv")
  fig_open("fig12_sensitivity", 9.6, 4.8)
  par(mfrow = c(1, 2))
  srcbar <- function(resp, letter) {
    z <- src(sens, resp); z <- z[order(abs(z$src)), ]
    par(mar = c(3.8, 10.2, 2.4, 0.8))
    barplot(z$src, horiz = TRUE, names.arg = PRETTY[z$parameter], las = 1, cex.names = 0.75,
            col = COL_NEUTRAL, border = NA, xlab = "standardised regression coefficient",
            xlim = c(-1, 1) * max(0.45, max(abs(z$src)) * 1.05))
    abline(v = 0, col = "grey30"); panel(letter, x = -0.50)
  }
  srcbar("gap_nominal", "a"); srcbar("bias_survey_total", "b")
  fig_close("fig12_sensitivity")
}

for (f in c("fig1", "fig2", "fig3", "fig4", "fig5", "fig6", "fig7", "fig8", "fig9", "fig10", "fig11", "fig12")) {
  message(f); get(f)()
}
message("figures written to figures/ (PDF + 300 dpi PNG)")
