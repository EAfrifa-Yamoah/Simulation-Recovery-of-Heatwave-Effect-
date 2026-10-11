# =============================================================================
# 03_estimators.R
# Estimator suite.
#
# Revision notes:
#   P3 / 3.2    zeros enter the primary index through a two part delta model.
#   P4.2        five estimators are benchmarked, not one.
#   P4.3        the demeaned index is retained only as a naive comparator and is
#               compared against a jointly fitted two way model.
#   P3.14       survey index and centroid are inclusion probability weighted.
#   P3.15       two centroids are reported: the raw fleet weighted catch
#               centroid and a fixed reference grid centroid from the delta
#               model predictions.
#
# The two way models are solved by alternating conditional maximisation, which
# is exact for the Gaussian part and converges quickly for the logistic part.
# This keeps a full delta model affordable across thousands of replicates.
# =============================================================================

# ---- grouped sums that keep every level -------------------------------------
gsum <- function(x, g, n) {
  out <- numeric(n)
  tmp <- rowsum(x, g, reorder = TRUE)
  out[as.integer(rownames(tmp))] <- tmp[, 1]
  out
}
gcount <- function(g, n) tabulate(g, nbins = n)

# ---- fast two way Gaussian: y = mu + a_i + b_j + e ---------------------------
twoway_gauss <- function(y, i, j, ni, nj, tol = 1e-9, maxit = 300) {
  a <- numeric(ni); b <- numeric(nj)
  mu <- mean(y)
  ci <- pmax(gcount(i, ni), 1); cj <- pmax(gcount(j, nj), 1)
  for (it in seq_len(maxit)) {
    a_new <- gsum(y - b[j] - mu, i, ni) / ci
    a_new <- a_new - mean(a_new)
    b_new <- gsum(y - a_new[i] - mu, j, nj) / cj
    if (max(abs(a_new - a), abs(b_new - b)) < tol) { a <- a_new; b <- b_new; break }
    a <- a_new; b <- b_new
  }
  res <- y - mu - a[i] - b[j]
  list(mu = mu, a = a, b = b, converged = it < maxit,
       sigma2 = sum(res^2) / max(1, length(y) - ni - nj))
}

# ---- fast two way logistic: logit(p) = a_i + b_j -----------------------------
# Diagonal Newton with a light ridge on the cell and month effects (penalised
# likelihood). The ridge is needed because cells with no encounters, or months
# with none, have no finite maximum likelihood estimate (complete separation);
# without it the solver cannot converge and the effect runs to minus infinity.
twoway_logit <- function(y, i, j, ni, nj, tol = 1e-6, maxit = 300, ridge = 0.05) {
  a <- rep(0, ni); b <- rep(0, nj)
  p0 <- min(max(mean(y), 1e-4), 1 - 1e-4)
  b[] <- log(p0 / (1 - p0))
  upd <- function(par, idx, other, n, pen) {
    eta  <- par[idx] + other
    p    <- 1 / (1 + exp(-eta))
    w    <- p * (1 - p) + 1e-8
    g    <- gsum(y - p, idx, n) - pen * par
    h    <- pmax(gsum(w, idx, n), 1e-6) + pen
    step <- g / h
    step[!is.finite(step)] <- 0
    par + pmax(pmin(step, 2), -2)
  }
  for (it in seq_len(maxit)) {
    a_new <- upd(a, i, b[j], ni, ridge); a_new <- a_new - mean(a_new)
    b_new <- upd(b, j, a_new[i], nj, 0)          # month effects unpenalised
    if (max(abs(a_new - a), abs(b_new - b)) < tol) { a <- a_new; b <- b_new; break }
    a <- a_new; b <- b_new
  }
  list(a = a, b = b, converged = it < maxit, iterations = it)
}

# ---- E_raw: nominal CPUE, total catch over total effort ----------------------
est_nominal <- function(sim, Ivec) {
  Y <- sim$fish$Y; E <- sim$fish$E
  Tn <- ncol(Y)
  num <- colSums(Y * E, na.rm = TRUE)
  den <- colSums(E * !is.na(Y), na.rm = TRUE)
  idx <- log(pmax(num / pmax(den, 1e-8), 1e-8))
  f <- slope_se(idx, Ivec)
  list(index = idx, beta = f$beta, se = f$se)
}

# ---- E0: naive positive only two way demeaned index -------------------------
est_naive <- function(sim, Ivec) {
  Y <- sim$fish$Y
  S <- nrow(Y); Tn <- ncol(Y)
  ok <- !is.na(Y) & Y > 0
  if (sum(ok) < 50) return(list(index = rep(NA, Tn), beta = NA_real_))
  lc  <- matrix(NA_real_, S, Tn); lc[ok] <- log(Y[ok])
  cs  <- rowMeans(lc, na.rm = TRUE)
  cs[!is.finite(cs)] <- 0
  moy <- ((seq_len(Tn) - 1) %% 12) + 1
  ms  <- tapply(as.numeric(lc), factor(rep(moy, each = S), levels = 1:12),
                mean, na.rm = TRUE)
  ms[!is.finite(ms)] <- 0
  msmat <- matrix(ms[moy], S, Tn, byrow = TRUE)
  dev <- lc - cs - msmat
  idx <- colMeans(dev, na.rm = TRUE)
  list(index = idx, beta = slope(idx, Ivec))
}

# ---- E1: two way fixed effect delta model, fixed reference grid --------------
est_delta_fe <- function(sim, Ivec, grid) {
  Y <- sim$fish$Y
  S <- nrow(Y); Tn <- ncol(Y)
  fished <- !is.na(Y)
  if (sum(fished) < 100) return(list(index = rep(NA, Tn), beta = NA_real_,
                                     pgrid = NULL))
  cell <- row(Y)[fished]; tim <- col(Y)[fished]
  yy   <- Y[fished]
  enc  <- as.numeric(yy > 0)

  fe_p <- twoway_logit(enc, cell, tim, S, Tn)
  posi <- yy > 0
  fe_m <- twoway_gauss(log(yy[posi]), cell[posi], tim[posi], S, Tn)

  # prediction on the complete reference grid at fixed reference effort.
  # Cells never fished in the series receive the mean cell effect: this is an
  # imputation, and its share is reported so that its influence can be judged.
  eta_p <- outer(fe_p$a, fe_p$b, "+")
  phat  <- 1 / (1 + exp(-eta_p))
  eta_m <- fe_m$mu + outer(fe_m$a, fe_m$b, "+")
  mhat  <- exp(eta_m + 0.5 * fe_m$sigma2)
  dens  <- phat * mhat
  idx   <- log(colMeans(dens))
  f     <- slope_se(idx, Ivec)

  # variant: average only over cells fished in that month (no imputation)
  dens_f <- dens; dens_f[!fished] <- NA
  idx_f  <- log(colMeans(dens_f, na.rm = TRUE))

  # imputation share during event and calm months
  J <- sim$J; ev_m <- J > 0
  unf <- 1 - colMeans(fished)
  list(index = idx, beta = f$beta, se = f$se, pgrid = dens,
       index_fished = idx_f, beta_fished = slope(idx_f, Ivec),
       converged = fe_p$converged && fe_m$converged,
       imputed_event = mean(unf[ev_m]), imputed_calm = mean(unf[!ev_m]))
}


# ---- E1b: delta model with a spatial severity interaction --------------------
# The additive two way model cannot move the centre of gravity: its spatial
# shape is fixed by construction. Adding a low dimensional space by severity
# interaction (a latitudinal and an offshore slope) gives a standardised
# density surface on fixed support that is free to move with the event, which
# is what the centroid diagnostic requires (review item 3.15).
est_delta_int <- function(sim, Ivec, grid, maxit = 12) {
  Y <- sim$fish$Y
  S <- nrow(Y); Tn <- ncol(Y)
  fished <- !is.na(Y)
  if (sum(fished) < 100)
    return(list(index = rep(NA, Tn), beta = NA_real_, dens = NULL))
  cell <- row(Y)[fished]; tim <- col(Y)[fished]; yy <- Y[fished]
  enc  <- as.numeric(yy > 0)
  Zi   <- cbind(grid$y[cell], grid$x[cell]) * Ivec[tim]

  ## encounter component
  gam_p <- c(0, 0)
  for (it in seq_len(maxit)) {
    off <- as.numeric(Zi %*% gam_p)
    fe  <- twoway_logit_off(enc, cell, tim, S, Tn, off)
    eta <- fe$a[cell] + fe$b[tim]
    pr  <- 1 / (1 + exp(-(eta + off)))
    w   <- pr * (1 - pr) + 1e-8
    z   <- (enc - pr) / w + off
    new <- tryCatch(as.numeric(coef(lm.wfit(Zi, z, w))), error = function(e) gam_p)
    new[!is.finite(new)] <- 0
    if (max(abs(new - gam_p)) < 1e-6) { gam_p <- new; break }
    gam_p <- new
  }
  off_p <- as.numeric(Zi %*% gam_p)
  fe_p  <- twoway_logit_off(enc, cell, tim, S, Tn, off_p)

  ## positive component
  posi <- yy > 0
  Zp   <- Zi[posi, , drop = FALSE]
  lv   <- log(yy[posi]); cp <- cell[posi]; tp <- tim[posi]
  th <- c(0, 0)
  for (it in seq_len(maxit)) {
    fe <- twoway_gauss_off(lv, cp, tp, S, Tn, as.numeric(Zp %*% th))
    r  <- lv - fe$mu - fe$a[cp] - fe$b[tp]
    new <- tryCatch(as.numeric(coef(.lm.fit(Zp, r))), error = function(e) th)
    new[!is.finite(new)] <- 0
    if (max(abs(new - th)) < 1e-8) { th <- new; break }
    th <- new
  }
  fe_m <- twoway_gauss_off(lv, cp, tp, S, Tn, as.numeric(Zp %*% th))

  ## fixed support prediction
  Zg <- cbind(grid$y, grid$x)
  int_p <- outer(as.numeric(Zg %*% gam_p), Ivec, "*")
  int_m <- outer(as.numeric(Zg %*% th),   Ivec, "*")
  phat <- 1 / (1 + exp(-(outer(fe_p$a, fe_p$b, "+") + int_p)))
  mhat <- exp(fe_m$mu + outer(fe_m$a, fe_m$b, "+") + int_m + 0.5 * fe_m$sigma2)
  dens <- phat * mhat
  idx  <- log(colMeans(dens))
  list(index = idx, beta = slope(idx, Ivec), dens = dens,
       gam_p = gam_p, theta = th)
}

twoway_gauss_off <- function(y, i, j, ni, nj, off, tol = 1e-9, maxit = 300) {
  yy <- y - off
  twoway_gauss(yy, i, j, ni, nj, tol, maxit)
}

twoway_logit_off <- function(y, i, j, ni, nj, off, tol = 1e-6, maxit = 300,
                             ridge = 0.05) {
  a <- rep(0, ni); b <- rep(0, nj)
  p0 <- min(max(mean(y), 1e-4), 1 - 1e-4); b[] <- log(p0 / (1 - p0))
  upd <- function(par, idx, other, n, pen) {
    eta <- par[idx] + other
    p <- 1 / (1 + exp(-eta)); w <- p * (1 - p) + 1e-8
    g <- gsum(y - p, idx, n) - pen * par; h <- pmax(gsum(w, idx, n), 1e-6) + pen
    st <- g / h; st[!is.finite(st)] <- 0
    par + pmax(pmin(st, 2), -2)
  }
  for (it in seq_len(maxit)) {
    a_new <- upd(a, i, b[j] + off, ni, ridge); a_new <- a_new - mean(a_new)
    b_new <- upd(b, j, a_new[i] + off, nj, 0)
    if (max(abs(a_new - a), abs(b_new - b)) < tol) { a <- a_new; b <- b_new; break }
    a <- a_new; b <- b_new
  }
  list(a = a, b = b, converged = it < maxit)
}

# ---- E3b: model based survey index -------------------------------------------
# A delta model fitted to the survey observations, predicted on the complete
# frame. Reported alongside the design based estimator so the two modes of
# inference are not conflated (review item 3.14).
est_survey_model <- function(surv, Ivec, grid) {
  V <- surv$V; S <- nrow(V); Tn <- ncol(V)
  ok <- !is.na(V)
  if (sum(ok) < 100) return(list(index = rep(NA, Tn), beta = NA_real_))
  cell <- row(V)[ok]; tim <- col(V)[ok]; vv <- V[ok]
  enc <- as.numeric(vv > 0)
  fe_p <- twoway_logit(enc, cell, tim, S, Tn)
  posi <- vv > 0
  if (sum(posi) < 50) return(list(index = rep(NA, Tn), beta = NA_real_))
  fe_m <- twoway_gauss(log(vv[posi]), cell[posi], tim[posi], S, Tn)
  phat <- 1 / (1 + exp(-outer(fe_p$a, fe_p$b, "+")))
  mhat <- exp(fe_m$mu + outer(fe_m$a, fe_m$b, "+") + 0.5 * fe_m$sigma2)
  idx  <- log(colMeans(phat * mhat))
  list(index = idx, beta = slope(idx, Ivec))
}

# ---- E2: spatiotemporal delta GAM -------------------------------------------
# Encounter and positive components with a two dimensional spatial smooth and a
# smooth event interaction, standing in for the geostatistical model class.
est_delta_gam <- function(sim, Ivec, grid, k_sp = 24) {
  if (!requireNamespace("mgcv", quietly = TRUE))
    return(list(index = NA, beta = NA_real_))
  Y <- sim$fish$Y
  S <- nrow(Y); Tn <- ncol(Y)
  fished <- !is.na(Y)
  cell <- row(Y)[fished]; tim <- col(Y)[fished]; yy <- Y[fished]
  d <- data.frame(enc = as.numeric(yy > 0), y = grid$y[cell], x = grid$x[cell],
                  tt = factor(tim, levels = seq_len(Tn)),
                  Iv = Ivec[tim], val = yy)
  ok <- tryCatch({
    g1 <- mgcv::bam(enc ~ s(x, y, k = k_sp) + s(x, y, by = Iv, k = k_sp) + tt,
                    family = binomial(), data = d, discrete = TRUE)
    dp <- subset(d, val > 0); dp$lv <- log(dp$val)
    g2 <- mgcv::bam(lv ~ s(x, y, k = k_sp) + s(x, y, by = Iv, k = k_sp) + tt,
                    data = dp, discrete = TRUE)
    nd <- expand.grid(cell = seq_len(S), tt = factor(seq_len(Tn), levels = seq_len(Tn)))
    nd$x <- grid$x[nd$cell]; nd$y <- grid$y[nd$cell]
    nd$Iv <- Ivec[as.integer(as.character(nd$tt))]
    p <- predict(g1, nd, type = "response")
    m <- exp(predict(g2, nd) + 0.5 * g2$sig2)
    idx <- log(tapply(p * m, nd$tt, mean))
    list(index = as.numeric(idx), beta = slope(as.numeric(idx), Ivec))
  }, error = function(e) list(index = NA, beta = NA_real_))
  ok
}

# ---- E3: design based survey index and centroid ------------------------------
est_survey <- function(surv, Ivec, grid) {
  V <- surv$V; S <- nrow(V); Tn <- ncol(V)
  tot <- rep(NA_real_, Tn); cg <- rep(NA_real_, Tn); cgx <- rep(NA_real_, Tn)
  for (t in seq_len(Tn)) {
    sel <- which(surv$sampled[, t])
    if (!length(sel)) next
    v <- V[sel, t]
    tot[t] <- sum(v / surv$p_inc)                   # Horvitz Thompson total
    if (sum(v) > 0) {                               # Hajek ratio centroids
      cg[t]  <- sum(grid$y[sel] * v) / sum(v)
      cgx[t] <- sum(grid$x[sel] * v) / sum(v)
    }
  }
  # a survey month with no catch at all has a total of zero; its log index is
  # floored at half the smallest positive total of the series (the usual
  # convention), so that one empty month cannot dominate the heatwave slope
  pos <- tot[is.finite(tot) & tot > 0]
  floor_v <- if (length(pos)) 0.5 * min(pos) else 1e-8
  idx <- log(pmax(tot, floor_v)); idx[!is.finite(tot)] <- NA
  f <- slope_se(idx, Ivec)
  list(index = idx, beta = f$beta, se = f$se, cg = cg, cgx = cgx,
       n_empty = sum(is.finite(tot) & tot == 0))
}

# ---- centroids ---------------------------------------------------------------
#' Raw fleet weighted catch centroid (a fleet use indicator) and the fixed
#' reference grid centroid from the delta model predictions.
centroids <- function(sim, delta_fit, grid, Ivec, tau = 0) {
  Y <- sim$fish$Y; Tn <- ncol(Y)
  cg_raw <- rep(NA_real_, Tn); cgx_raw <- rep(NA_real_, Tn)
  for (t in seq_len(Tn)) {
    v <- Y[, t]; ok <- !is.na(v) & v > 0
    if (sum(ok) > 1) {
      cg_raw[t]  <- sum(grid$y[ok] * v[ok]) / sum(v[ok])
      cgx_raw[t] <- sum(grid$x[ok] * v[ok]) / sum(v[ok])
    }
  }
  cg_grid <- rep(NA_real_, Tn); cgx_grid <- rep(NA_real_, Tn)
  if (!is.null(delta_fit$pgrid)) {
    d <- delta_fit$pgrid
    cg_grid  <- as.numeric(colSums(grid$y * d) / colSums(d))
    cgx_grid <- as.numeric(colSums(grid$x * d) / colSums(d))
  }
  list(raw = cg_raw, raw_x = cgx_raw, grid = cg_grid, grid_x = cgx_grid)
}

#' Event contrast of a monthly series: mean during event months minus mean
#' outside them.
event_contrast <- function(v, Jvec, tau = 0) {
  hi <- Jvec > tau; lo <- !hi
  if (!any(hi & !is.na(v)) || !any(lo & !is.na(v))) return(NA_real_)
  mean(v[hi], na.rm = TRUE) - mean(v[lo], na.rm = TRUE)
}

# ---- heatwave coefficient ----------------------------------------------------
#' Every index coefficient in the study is the partial slope on the centred
#' event severity covariate after controlling for a linear time trend and a
#' first order seasonal harmonic. Controlling for the chronic trend and
#' seasonality is what separates the event anomaly from background warming
#' (review item 3.5) and follows CPUE good practice (Hoyle et al., 2024).
#' The identical specification is applied to the estimand and to every
#' estimator, so all coefficients target the same quantity.
slope_se <- function(y, x, tvec = NULL, lag = 3) {
  n  <- length(y)
  tt <- if (is.null(tvec)) seq_len(n) else tvec
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) < 8) return(list(beta = NA_real_, se = NA_real_, rank = NA_integer_))
  X <- cbind(1, x[ok], tt[ok] / n,
             sin(2 * pi * tt[ok] / 12), cos(2 * pi * tt[ok] / 12))
  # drop collinear controls (a semiannual survey cannot separate the seasonal
  # sine and cosine); the severity column is kept whenever it is estimable
  q <- qr(X); keep <- sort(q$pivot[seq_len(q$rank)])
  if (!(2 %in% keep)) return(list(beta = NA_real_, se = NA_real_, rank = q$rank))
  X <- X[, keep, drop = FALSE]; j <- which(keep == 2)
  if (nrow(X) <= ncol(X) + 2) return(list(beta = NA_real_, se = NA_real_, rank = q$rank))
  fit <- try(.lm.fit(X, y[ok]), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(beta = NA_real_, se = NA_real_, rank = q$rank))
  e <- fit$residuals; m <- nrow(X)
  XtXi <- tryCatch(solve(crossprod(X)), error = function(z) NULL)
  if (is.null(XtXi)) return(list(beta = unname(coef(fit)[j]), se = NA_real_, rank = q$rank))
  # Newey West (Bartlett) long run covariance of the scores
  Xe <- X * e
  Om <- crossprod(Xe)
  for (l in seq_len(min(lag, m - 2))) {
    w <- 1 - l / (lag + 1)
    G <- crossprod(Xe[(l + 1):m, , drop = FALSE], Xe[1:(m - l), , drop = FALSE])
    Om <- Om + w * (G + t(G))
  }
  V <- XtXi %*% Om %*% XtXi
  list(beta = unname(coef(fit)[j]), se = sqrt(max(V[j, j], 0)), rank = q$rank)
}

slope <- function(y, x, tvec = NULL) slope_se(y, x, tvec)$beta
