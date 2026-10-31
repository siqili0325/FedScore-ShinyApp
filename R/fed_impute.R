# Missing values: FedIMPUTE (optional step before ranking) -----------------------
#
# Li et al., "FedIMPUTE: Privacy-preserving missing value imputation for
# multi-site heterogeneous electronic health records", JBI 2025. This is the
# DAC-MICE variant:
#   1. Each site completes its predictors locally with MICE (one imputation).
#      These values are only used as inputs to the imputation models.
#   2. For every variable to impute, a regression on the other predictors
#      (linear for numeric variables, logistic for binary ones) is fitted on
#      the rows where that variable is observed and federated across sites:
#      round 1 shares local coefficients, whose average is broadcast; each
#      further round shares gradients and Hessians at the broadcast
#      coefficients, and the lead takes one Newton step (divide-and-conquer
#      one-step estimation).
#   3. Each site fills its own missing values with the federated models,
#      optionally after recalibrating them on its own data (transfer learning).
# The outcome columns are never used, and imputed values never leave the site.

# Predictors that can be imputed, with the regression family used.
impute_family <- function(type, binary, levels) {
  if (type == "categorical") return(if (length(levels) == 2) "binomial" else NA_character_)
  if (isTRUE(binary)) "binomial" else "gaussian"
}

# Local completion of `df` (predictors only) with MICE; cached, because the
# same split is completed several times during a study.
.fill_cache <- new.env()
local_fill <- function(df, maxit = 5, seed = 2024) {
  if (!anyNA(df)) return(df)
  key <- digest::digest(list(df, maxit, seed))
  if (!is.null(.fill_cache[[key]])) return(.fill_cache[[key]])
  imp <- mice::mice(df, m = 1, maxit = maxit, seed = seed, printFlag = FALSE)
  out <- mice::complete(imp)
  # MICE leaves a column empty when it cannot model it (e.g. constant);
  # fall back to the median or most frequent value.
  for (v in names(out)) {
    miss <- is.na(out[[v]])
    if (!any(miss)) next
    out[[v]][miss] <- if (is.factor(out[[v]])) names(which.max(table(out[[v]]))) else stats::median(out[[v]], na.rm = TRUE)
  }
  assign(key, out, envir = .fill_cache)
  out
}

# Design matrix of the imputation model for `target`: intercept plus the other
# predictors (categorical ones as indicators of the study-wide levels).
impute_design <- function(df, target, predictors, levels) {
  cols <- list(`(Intercept)` = rep(1, nrow(df)))
  for (p in setdiff(predictors, target)) {
    lv <- levels[[p]]
    if (is.null(lv)) cols[[p]] <- as.numeric(df[[p]])
    else for (l in lv[-1]) cols[[paste0(p, " ", l)]] <- as.numeric(as.character(df[[p]]) == l)
  }
  do.call(cbind, cols)
}

# 0/1 or numeric response of an imputation model.
impute_response <- function(x, family, levels) {
  if (family == "gaussian") return(as.numeric(x))
  if (is.factor(x) || is.character(x)) return(as.numeric(as.character(x) == levels[2]))
  as.numeric(x)
}

impute_lik <- function(beta, X, y, family, deriv = TRUE, ridge = 1e-6) {
  if (family == "binomial") return(logistic_lik(beta, X, y, deriv, ridge))
  r <- y - drop(X %*% beta)
  out <- list(value = -sum(r^2) / (2 * length(y)))
  if (deriv) {
    out$grad <- drop(crossprod(X, r)) / length(y)
    out$hess <- -crossprod(X) / length(y) - diag(ridge, ncol(X))
  }
  out
}

# The rows a site can use for one imputation model, or NULL if there are too
# few (fewer than `min_count` observed values, or of either class).
impute_rows <- function(train, filled, model, req) {
  obs <- !is.na(train[[model$var]])
  y <- impute_response(train[[model$var]][obs], model$family, req$levels[[model$var]])
  X <- impute_design(filled[obs, , drop = FALSE], model$var, req$predictors, req$levels)
  k <- max(if (is.null(req$min_count)) 0 else req$min_count, ncol(X) + 1)
  small <- if (model$family == "binomial") min(sum(y), sum(1 - y)) else length(y)
  if (small < k) return(NULL)
  list(X = X, y = y)
}

# Site reply to round 1 (local coefficients) or later rounds (derivatives).
impute_reply <- function(req, data) {
  train <- data$train
  filled <- local_fill(train[, req$predictors, drop = FALSE], req$maxit, req$seed)
  list(models = lapply(req$models, function(m) {
    xy <- impute_rows(train, filled, m, req)
    if (is.null(xy)) return(list(n = 0))
    if (req$kind == "impute_local") {
      fit <- newton_max(function(b, d) impute_lik(b, xy$X, xy$y, m$family, d), rep(0, ncol(xy$X)))
      return(list(n = length(xy$y), beta = unname(fit$beta)))
    }
    d <- impute_lik(m$coef, xy$X, xy$y, m$family)
    list(n = length(xy$y), grad = unname(d$grad), hess = d$hess)
  }))
}

# Lead site: average the local coefficients (round 1) or take a Newton step
# from the aggregated derivatives (later rounds). Sites are weighted by their
# number of observed values, so the step targets the pooled-data fit.
impute_combine <- function(req, replies, m) {
  a <- lapply(replies, function(r) r$models[[m]])
  n <- vapply(a, function(x) as.numeric(x$n), numeric(1))
  if (!any(n > 0)) stop(sprintf("No site has enough observed values of %s to fit its imputation model.", req$models[[m]]$var))
  a <- a[n > 0]
  w <- n[n > 0] / sum(n)
  if (req$kind == "impute_local") return(drop(Reduce(`+`, Map(function(x, wk) wk * x$beta, a, w))))
  total <- aggregate_derivs(lapply(a, function(x) list(grad = x$grad, hess = as.matrix(x$hess))), w)
  req$models[[m]]$coef + safe_solve(-total$hess, total$grad)
}

# Fill the missing values of every split with the federated models. Returns
# the data with the imputed variables complete (other missing values are
# kept, and later become the "Unknown" category); attribute "imputed" counts
# the values filled per split and variable.
apply_imputation <- function(data, res) {
  preds <- unlist(res$predictors)
  out <- data
  counts <- list()
  for (s in names(SPLIT_ALIASES)) {
    raw <- data[[s]]
    filled <- local_fill(raw[, preds, drop = FALSE], res$maxit, res$seed)
    for (j in seq_along(res$models)) {
      m <- res$models[[j]]
      miss <- is.na(raw[[m$var]])
      counts[[s]][[m$var]] <- sum(miss)
      if (!any(miss)) next
      coef <- m$coef + calibration(data$train, m, res, s)
      eta <- drop(impute_design(filled[miss, , drop = FALSE], m$var, preds, res$levels) %*% coef)
      lv <- res$levels[[m$var]]
      val <- if (m$family == "gaussian") eta else {
        set.seed(res$seed + j)
        stats::rbinom(length(eta), 1, stats::plogis(eta))
      }
      if (!is.null(lv)) val <- factor(lv[val + 1], levels = levels(raw[[m$var]]))
      filled[[m$var]][miss] <- val
      out[[s]][[m$var]][miss] <- val
    }
  }
  attr(out, "outcome") <- attr(data, "outcome")
  attr(out, "imputed") <- counts
  out
}

# Optional transfer-learning recalibration: a lasso correction of the
# federated model on the site's own observed values (offset = federated
# prediction), computed once from the training set. Zero when switched off.
.calib_cache <- new.env()
calibration <- function(train, m, res, split) {
  zero <- numeric(length(m$coef))
  if (!isTRUE(res$calibrate)) return(zero)
  key <- digest::digest(list(train, m, res$seed, res$maxit))
  if (!is.null(.calib_cache[[key]])) return(.calib_cache[[key]])
  preds <- unlist(res$predictors)
  filled <- local_fill(train[, preds, drop = FALSE], res$maxit, res$seed)
  xy <- impute_rows(train, filled, m, list(predictors = preds, levels = res$levels, min_count = 0))
  delta <- zero
  if (!is.null(xy) && length(xy$y) >= 50) {
    set.seed(res$seed)
    cv <- tryCatch(glmnet::cv.glmnet(xy$X[, -1, drop = FALSE], xy$y, family = m$family,
                                     offset = drop(xy$X %*% m$coef), alpha = 1, nfolds = 5),
                   error = function(e) NULL)
    if (!is.null(cv)) delta <- as.numeric(stats::coef(cv, s = "lambda.min"))
  }
  assign(key, delta, envir = .calib_cache)
  delta
}
