# FedScore Module 3: federated logistic regression ----------------------------
#
# Two algorithms are available:
#   * ODAL2 (Duan et al., JAMIA 2020) for homogeneous sites. A lead site
#     combines its own likelihood with gradients/Hessians from every site.
#   * dCLR (Tong et al., npj Digit Med 2022) for heterogeneous sites. Pairwise
#     conditional likelihood removes site-specific intercepts (baseline risk).
# This file holds the computations; fed_protocol.R decides which site runs
# which piece and what is exchanged.

# Design matrices -------------------------------------------------------------

# Treatment-coded indicators for categorised data (first level = reference).
# The `map` attribute links each column back to its variable and category.
design_matrix <- function(cat_df, intercept = TRUE) {
  cols <- list()
  map <- data.frame(column = character(0), variable = character(0), level = character(0))
  for (v in names(cat_df)) {
    x <- cat_df[[v]]
    for (l in levels(x)[-1]) {
      nm <- paste0(v, " ", l)
      cols[[nm]] <- as.numeric(!is.na(x) & x == l)
      map[nrow(map) + 1, ] <- c(nm, v, l)
    }
  }
  X <- do.call(cbind, cols)
  if (is.null(X)) X <- matrix(numeric(0), nrow = nrow(cat_df), ncol = 0)
  if (intercept) {
    X <- cbind(`(Intercept)` = 1, X)
    map <- rbind(data.frame(column = "(Intercept)", variable = "(Intercept)", level = ""), map)
  }
  attr(X, "map") <- map
  X
}

# The columns design_matrix() produces for `vars`, worked out from the cut-offs
# alone (no data needed). The `levels` attribute lists every category,
# reference first, for building score tables.
design_terms <- function(cutoffs, vars, intercept = TRUE) {
  lv <- lapply(setNames(vars, vars), function(v) var_levels(cutoffs[[v]]))
  map <- do.call(rbind, c(
    if (intercept) list(data.frame(column = "(Intercept)", variable = "(Intercept)", level = "")),
    lapply(vars, function(v) data.frame(column = paste0(v, " ", lv[[v]][-1]), variable = rep(v, length(lv[[v]]) - 1),
                                        level = lv[[v]][-1]))))
  attr(map, "levels") <- lv
  map
}

# Categorise one data split and build its design matrix and outcome.
site_xy <- function(df, cutoffs, vars, intercept) {
  list(X = design_matrix(categorize(df, cutoffs, vars), intercept), y = df$label)
}

# Optimiser -------------------------------------------------------------------

safe_solve <- function(A, b) {
  out <- tryCatch(solve(A, b), error = function(e) NULL)
  if (is.null(out)) out <- MASS::ginv(A) %*% b
  drop(out)
}

safe_inv <- function(A) {
  out <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(out)) out <- MASS::ginv(A)
  out
}

# Maximise a concave objective. `f(beta, deriv)` returns list(value) and, when
# `deriv = TRUE`, also `grad` and `hess`.
newton_max <- function(f, beta0, maxit = 50, tol = 1e-7) {
  beta <- beta0
  cur <- f(beta, TRUE)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    step <- safe_solve(-cur$hess, cur$grad)
    t <- 1
    repeat {
      cand <- beta + t * step
      val <- f(cand, FALSE)$value
      if (is.finite(val) && val >= cur$value - 1e-12) break
      t <- t / 2
      if (t < 1e-6) break
    }
    if (t < 1e-6) { converged <- max(abs(cur$grad)) < 1e-6; break }
    beta <- cand
    prev <- cur$value
    cur <- f(beta, TRUE)
    if (max(abs(t * step)) < tol || abs(cur$value - prev) < 1e-12) { converged <- TRUE; break }
  }
  list(beta = beta, value = cur$value, grad = cur$grad, hess = cur$hess,
       iterations = it, converged = converged)
}

# Ordinary logistic likelihood (ODAL2) ----------------------------------------

# Mean log-likelihood with derivatives; `ridge` keeps the Hessian invertible
# when a category is empty at a site (the intercept is not penalised).
logistic_lik <- function(beta, X, y, deriv = TRUE, ridge = 1e-6) {
  eta <- drop(X %*% beta)
  pen <- ridge * c(0, rep(1, length(beta) - 1))
  out <- list(value = mean(y * eta - log1p(exp(-abs(eta))) - pmax(eta, 0)) - sum(pen * beta^2) / 2)
  if (deriv) {
    p <- stats::plogis(eta)
    out$grad <- drop(crossprod(X, y - p)) / length(y) - pen * beta
    out$hess <- -crossprod(X * (p * (1 - p)), X) / length(y) - diag(pen, length(beta))
  }
  out
}

logistic_mle <- function(X, y, ridge = 1e-6) {
  beta0 <- c(stats::qlogis(mean(y)), rep(0, ncol(X) - 1))
  newton_max(function(b, d) logistic_lik(b, X, y, d, ridge), beta0)
}

# Second-order surrogate likelihood shared by ODAL2 and dCLR: a site's own
# likelihood `lik`, corrected towards the all-site likelihood by the gap
# between the aggregated derivatives `total` and its own derivatives at
# beta_bar. Maximised by Newton steps from beta_bar.
surrogate_fit <- function(lik, beta_bar, total) {
  own <- lik(beta_bar, TRUE)
  dg <- total$grad - own$grad
  dH <- total$hess - own$hess
  f <- function(b, d) {
    l <- lik(b, d)
    diff <- b - beta_bar
    out <- list(value = l$value + sum(dg * diff) + 0.5 * drop(t(diff) %*% dH %*% diff))
    if (d) {
      out$grad <- l$grad + dg + drop(dH %*% diff)
      out$hess <- l$hess + dH
    }
    out
  }
  newton_max(f, beta_bar)
}

# Weighted sum of the sites' gradients and Hessians.
aggregate_derivs <- function(derivs, weights) {
  list(grad = Reduce(`+`, Map(function(d, w) w * d$grad, derivs, weights)),
       hess = Reduce(`+`, Map(function(d, w) w * d$hess, derivs, weights)))
}

# ODAL2 pieces. Design matrices include the intercept.
#   Lead site: odal2_init() gives beta_bar, shared with every site.
#   Every site: site_derivs("ODAL2", ...) at beta_bar.
#   Lead site: odal2_surrogate() with the aggregated derivatives; N is the
#   total training size across sites (for the standard errors).
odal2_init <- function(X, y, ridge = 1e-6) logistic_mle(X, y, ridge)$beta

odal2_surrogate <- function(X, y, beta_bar, total, N, ridge = 1e-6) {
  fit <- surrogate_fit(function(b, d) logistic_lik(b, X, y, d, ridge), beta_bar, total)
  list(beta = fit$beta, se = sqrt(pmax(diag(safe_inv(-N * fit$hess)), 0)), converged = fit$converged)
}

# Gradient and Hessian of a site's (mean) log-likelihood at `beta`.
site_derivs <- function(method, beta, X, y, ridge = 1e-6) {
  d <- if (method == "ODAL2") logistic_lik(beta, X, y, TRUE, ridge) else
    pairwise_lik(beta, collapse_patterns(X, y), TRUE, FALSE, ridge)
  list(grad = d$grad, hess = d$hess)
}

# Pairwise conditional likelihood (dCLR) ----------------------------------------
#
# For a site with n patients, the pairwise likelihood of Tong et al. is
#   l(b) = (1 / C(n,2)) * sum over case i, control k of log expit((x_i - x_k)' b),
# which matches the C implementation (mylik_odal.c and friends) up to an
# additive constant from same-outcome pairs. Rows are collapsed to unique
# covariate patterns, so the cost scales with patterns, not patients.

collapse_patterns <- function(X, y) {
  collapse_one <- function(Xg) {
    if (nrow(Xg) == 0) return(list(X = Xg, w = numeric(0)))
    key <- do.call(paste, c(as.data.frame(Xg), sep = "\r"))
    first <- !duplicated(key)
    list(X = Xg[first, , drop = FALSE], w = as.vector(table(factor(key, levels = key[first]))))
  }
  cases <- collapse_one(X[y == 1, , drop = FALSE])
  ctrls <- collapse_one(X[y == 0, , drop = FALSE])
  list(Xc = cases$X, wc = cases$w, Xk = ctrls$X, wk = ctrls$w, n = length(y), p = ncol(X))
}

# Value, derivatives and (optionally) the U-statistic sandwich variance.
pairwise_lik <- function(beta, pat, deriv = TRUE, variance = FALSE, ridge = 1e-6, block = 4e6) {
  p <- pat$p
  P <- pat$n * (pat$n - 1) / 2
  eta_c <- drop(pat$Xc %*% beta)
  eta_k <- drop(pat$Xk %*% beta)
  Xk_w <- pat$Xk * pat$wk
  value <- 0
  g_case <- numeric(p); col_W <- numeric(length(eta_k))
  H_cc <- matrix(0, p, p); H_ck <- matrix(0, p, p); col_V <- numeric(length(eta_k))
  u_c <- if (variance) matrix(0, length(eta_c), p) else NULL
  WtX <- if (variance) matrix(0, length(eta_k), p) else NULL
  rows <- max(1, floor(block / max(1, length(eta_k))))
  for (start in seq(1, length(eta_c), by = rows)) {
    i <- start:min(length(eta_c), start + rows - 1)
    D <- outer(eta_c[i], eta_k, "-")
    wc <- pat$wc[i]
    value <- value + sum(wc * (stats::plogis(D, log.p = TRUE) %*% pat$wk))
    if (!deriv) next
    W <- stats::plogis(-D)               # d/dD log expit(D)
    Ww <- drop(W %*% pat$wk)
    g_case <- g_case + drop(crossprod(pat$Xc[i, , drop = FALSE], wc * Ww))
    col_W <- col_W + drop(crossprod(W, wc))
    V <- W * (1 - W)
    Vw <- drop(V %*% pat$wk)
    Xci <- pat$Xc[i, , drop = FALSE]
    H_cc <- H_cc + crossprod(Xci * (wc * Vw), Xci)
    H_ck <- H_ck + crossprod(Xci * wc, V %*% Xk_w)
    col_V <- col_V + drop(crossprod(V, wc))
    if (variance) {
      u_c[i, ] <- Xci * Ww - W %*% Xk_w
      WtX <- WtX + crossprod(W * wc, Xci)
    }
  }
  pen <- ridge * sum(beta^2) / 2
  out <- list(value = value / P - pen)
  if (!deriv) return(out)
  out$grad <- (g_case - drop(crossprod(pat$Xk, pat$wk * col_W))) / P - ridge * beta
  H_kk <- crossprod(pat$Xk * (pat$wk * col_V), pat$Xk)
  out$hess <- -(H_cc + H_kk - H_ck - t(H_ck)) / P - diag(ridge, p)
  if (variance) {
    # Hajek projection of the U-statistic score: u_i = sum_k s_ik / (n - 1).
    u_k <- WtX - pat$Xk * col_W
    U <- rbind(u_c, u_k) / (pat$n - 1)
    w <- c(pat$wc, pat$wk)
    mu <- colSums(U * w) / sum(w)
    S1 <- crossprod(sweep(U, 2, mu) * sqrt(w)) / sum(w)
    A_inv <- safe_inv(out$hess)
    out$vcov <- A_inv %*% (4 * S1 / pat$n) %*% A_inv
  }
  out
}

pairwise_mle <- function(pat, beta0, ridge = 1e-6) {
  newton_max(function(b, d) pairwise_lik(b, pat, d, FALSE, ridge), beta0)
}

# Inverse-variance weighted average of site estimates.
iv_average <- function(betas, vcovs) {
  precs <- lapply(vcovs, safe_inv)
  P <- Reduce(`+`, precs)
  V <- safe_inv(P)
  list(beta = drop(V %*% Reduce(`+`, Map(function(Pr, b) Pr %*% b, precs, betas))), vcov = V)
}

# dCLR pieces. Design matrices have no intercept.
#   Every site: dclr_local() gives its estimate and sandwich variance; the
#   inverse-variance weighted average is beta_bar.
#   Every site: site_derivs("dCLR", ...) at beta_bar.
#   Every site: dclr_surrogate() with the aggregated derivatives; the
#   inverse-variance weighted average of these is the dCLR estimate.
# Repeating the last two steps around the new estimate ("surrogate passes")
# moves it towards the pooled-data fit, which helps when sites are very
# heterogeneous (e.g. a category without events at one site).
dclr_local <- function(X, y, ridge = 1e-6) {
  pat <- collapse_patterns(X, y)
  # Start from the slopes of an ordinary logistic fit, which are usually close.
  start <- logistic_mle(cbind(1, X), y, ridge)$beta[-1]
  fit <- pairwise_mle(pat, start, ridge)
  list(beta = fit$beta, vcov = pairwise_lik(fit$beta, pat, TRUE, TRUE, ridge)$vcov)
}

dclr_surrogate <- function(X, y, beta_bar, total, ridge = 1e-6) {
  pat <- collapse_patterns(X, y)
  fit <- surrogate_fit(function(b, d) pairwise_lik(b, pat, d, FALSE, ridge), beta_bar, total)
  list(beta = fit$beta, vcov = pairwise_lik(fit$beta, pat, TRUE, TRUE, ridge)$vcov, converged = fit$converged)
}
