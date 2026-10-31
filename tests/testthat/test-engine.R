# Run from the repository root with: Rscript -e 'testthat::test_dir("tests/testthat")'
# The app's functions are loaded by helper-study.R.

# Line-by-line R port of the original C code (examples/binary/hete_demo/C in the
# FedScore repository), used as the reference for the vectorised pairwise likelihood.
c_mylik_all <- function(n, y, X, par) {
  out <- 0
  for (j in 1:(n - 1)) for (l in (j + 1):n) {
    temp <- sum((X[j, ] - X[l, ]) * par)
    out <- out - log(1 + exp(-(y[j] - y[l]) * temp))
  }
  out
}
c_mylik_gradient <- function(n, y, X, par) {
  out <- numeric(length(par))
  for (j in 1:(n - 1)) for (l in (j + 1):n) {
    t1 <- sum((X[j, ] - X[l, ]) * par)
    temp <- exp(-(y[j] - y[l]) * t1)
    out <- out + (y[j] - y[l]) * (X[j, ] - X[l, ]) * (temp / (1 + temp))
  }
  out
}
c_mylik_gradient_second <- function(n, y, X, par) {
  p <- length(par)
  out <- matrix(0, p, p)
  for (j in 1:(n - 1)) for (l in (j + 1):n) {
    d <- X[j, ] - X[l, ]
    t1 <- sum(d * par)
    temp <- exp(-(y[j] - y[l]) * t1)
    out <- out - (y[j] - y[l])^2 * tcrossprod(d) * temp / (1 + temp)^2
  }
  out
}

sim_binary <- function(n, p, beta, intercept = -1, seed = 1) {
  set.seed(seed)
  X <- matrix(rbinom(n * p, 1, 0.4), n, p)
  y <- rbinom(n, 1, plogis(intercept + X %*% beta))
  list(X = X, y = y)
}

test_that("pairwise likelihood matches the original C implementation", {
  d <- sim_binary(60, 3, c(0.8, -0.5, 0.3))
  b <- c(0.4, -0.2, 0.1)
  n <- length(d$y)
  P <- n * (n - 1) / 2
  pat <- collapse_patterns(d$X, d$y)
  expect_lt(nrow(pat$Xc) + nrow(pat$Xk), n)   # rows really were collapsed
  got <- pairwise_lik(b, pat, ridge = 0)
  same_outcome_pairs <- choose(sum(d$y), 2) + choose(sum(1 - d$y), 2)
  expect_equal(got$value, (c_mylik_all(n, d$y, d$X, b) + same_outcome_pairs * log(2)) / P, tolerance = 1e-10)
  expect_equal(got$grad, c_mylik_gradient(n, d$y, d$X, b) / P, tolerance = 1e-10)
  expect_equal(got$hess, c_mylik_gradient_second(n, d$y, d$X, b) / P, tolerance = 1e-10)
})

test_that("chunked pairwise computation equals the unchunked one", {
  d <- sim_binary(400, 4, c(1, -1, 0.5, 0), seed = 2)
  pat <- collapse_patterns(d$X, d$y)
  b <- c(0.3, -0.3, 0.2, 0.1)
  a <- pairwise_lik(b, pat, variance = TRUE)
  c1 <- pairwise_lik(b, pat, variance = TRUE, block = 5)
  expect_equal(a, c1, tolerance = 1e-10)
})

test_that("logistic derivatives agree with finite differences", {
  d <- sim_binary(300, 3, c(1, -1, 0.5), seed = 3)
  X <- cbind(1, d$X)
  b <- c(-0.5, 0.2, 0.1, -0.3)
  an <- logistic_lik(b, X, d$y)
  num <- sapply(seq_along(b), function(k) {
    e <- replace(numeric(4), k, 1e-6)
    (logistic_lik(b + e, X, d$y, FALSE)$value - logistic_lik(b - e, X, d$y, FALSE)$value) / 2e-6
  })
  expect_equal(an$grad, num, tolerance = 1e-6)
  expect_equal(unname(logistic_mle(X, d$y, ridge = 0)$beta),
               unname(coef(glm(d$y ~ d$X, family = binomial))), tolerance = 1e-6)
})

test_that("surrogate fits with exact aggregated derivatives reach the pooled optimum", {
  beta <- c(1, -0.8, 0.5)
  a <- sim_binary(3000, 3, beta, intercept = -2, seed = 4)
  b <- sim_binary(1500, 3, beta, intercept = -0.5, seed = 5)
  Xs <- list(A = cbind(`(Intercept)` = 1, a$X), B = cbind(`(Intercept)` = 1, b$X))
  ys <- list(A = a$y, B = b$y)
  w <- site_weights(c(3000, 1500))
  # ODAL2 iterated to a fixed point is the size-weighted pooled MLE.
  bb <- odal2_init(Xs$A, ys$A)
  for (i in 1:5) {
    total <- aggregate_derivs(Map(function(X, y) site_derivs("ODAL2", bb, list(X = X, y = y)), Xs, ys), w)
    bb <- odal2_surrogate(Xs$A, ys$A, bb, total, N = 4500)$beta
  }
  expect_equal(unname(bb), unname(fit_pooled_logistic(Xs, ys)), tolerance = 1e-5)
})

test_that("categorisation gives identical levels at every site", {
  s1 <- data.frame(x = c(1:100), f = factor(rep(c("a", "b"), 50)), label = rep(0:1, 50))
  s2 <- data.frame(x = c(50:149), f = factor(rep(c("b", "c"), 50)), label = rep(0:1, 50))
  w <- c(0.5, 0.5)
  cuts <- federate_cutoffs(list(site_cut_stats(s1, c("x", "f")), site_cut_stats(s2, c("x", "f"))), w)
  expect_equal(cuts$f$levels, c("a", "b", "c"))
  c1 <- categorize(s1, cuts)
  c2 <- categorize(s2, cuts)
  expect_identical(levels(c1$x), levels(c2$x))
  expect_identical(levels(c1$f), levels(c2$f))
  expect_false(anyNA(c1$x))
  expect_identical(cut_labels(c(1, 2.5)), c("(,1)", "[1,2.5)", "[2.5,)"))
  expect_equal(parse_cuts("3, 1;2"), c(1, 2, 3))
})

test_that("score tables start at 0 per variable and total at most 100", {
  set.seed(8)
  s <- data.frame(x = rnorm(500), z = rnorm(500))
  s$label <- rbinom(500, 1, plogis(-1 + s$x - 0.5 * s$z))
  site <- list(train = s[1:300, ], validate = s[301:400, ], test = s[401:500, ])
  cuts <- federate_cutoffs(list(A = site_cut_stats(site$train, c("x", "z"))), 1)
  xy <- site_xy(site$train, cuts, c("x", "z"), TRUE)
  beta <- logistic_mle(xy$X, xy$y)$beta
  map <- design_terms(cuts, c("x", "z"))
  expect_identical(map$column, colnames(xy$X))
  tbl <- score_table(beta, map)
  expect_equal(as.vector(tapply(tbl$Points, tbl$Variable, min)), c(0, 0))
  expect_lte(sum(tapply(tbl$Points, tbl$Variable, max)), 101)
  sc <- compute_scores(s, tbl, cuts)
  expect_gt(auc_only(sc, s$label), 0.65)
  tl <- fit_tlscore(beta, "ODAL2", c("x", "z"), site, cuts)
  expect_equal(nrow(tl$table), nrow(tbl))
})

test_that("site data in all supported layouts is normalised", {
  df <- data.frame(a = 1:50, b = letters[rep(1:2, 25)], label = rep(c(TRUE, FALSE), 25))
  s <- normalize_site_data(df)
  expect_named(s, c("train", "validate", "test"))
  expect_true(is.factor(s$train$b))
  expect_true(all(s$train$label %in% 0:1))
  s2 <- normalize_site_data(list(non_cv = list(train = df, validate = df, test = df)))
  expect_equal(nrow(s2$test), 50)
  expect_error(normalize_site_data(list(train = df)), "validate")
  expect_equal(as_binary_label(factor(c("No", "Yes", "No"), levels = c("No", "Yes"))), c(0L, 1L, 0L))
  expect_equal(attr(s, "outcome"), "binary")
  surv <- data.frame(a = 1:40, label_time = c(0, 1:39), label_status = rep(c(TRUE, FALSE), 20))
  s3 <- normalize_site_data(list(train = surv, validate = surv, test = surv))
  expect_equal(attr(s3, "outcome"), "survival")
  expect_identical(s3$train$label_status, rep(c(1L, 0L), 20))
  expect_equal(predictor_names(s3$train), "a")
  bad <- surv
  bad$label_time[3] <- -1
  expect_error(normalize_site_data(list(train = bad, validate = surv, test = surv)), "negative")
  expect_error(normalize_site_data(list(train = surv, validate = df, test = surv)), "same outcome")
})

test_that("cut-offs at a site's minimum or maximum are dropped (no empty categories)", {
  set.seed(7)
  mk <- function(n) data.frame(k = c(rep(0, n * 0.3), round(runif(n * 0.7, 1, 60))), label = rbinom(n, 1, 0.2))
  s1 <- mk(400); s2 <- mk(200)
  cuts <- federate_cutoffs(list(site_cut_stats(s1, "k"), site_cut_stats(s2, "k")), c(2, 1) / 3)
  expect_gt(min(cuts$k$cuts), 0)
  tab <- table(categorize(rbind(s1, s2), cuts)$k)
  expect_true(all(tab > 0))
})

# Survival models ---------------------------------------------------------------

sim_surv <- function(n, beta, seed, rate = 0.1) {
  set.seed(seed)
  X <- cbind(a = rbinom(n, 1, 0.4), b = rbinom(n, 1, 0.5), c = rbinom(n, 1, 0.3))
  time <- round(rexp(n, rate * exp(drop(X %*% beta))), 1)       # rounding creates ties
  list(X = X, time = time, y = rbinom(n, 1, 0.75))
}

test_that("the Cox partial likelihood matches survival::coxph (Breslow ties)", {
  d <- sim_surv(400, c(0.7, -0.4, 0.3), 1)
  fit <- cox_mle(d$X, d$time, d$y, ridge = 0)
  ref <- survival::coxph(survival::Surv(d$time, d$y) ~ d$X, ties = "breslow")
  expect_equal(unname(fit$beta), unname(coef(ref)), tolerance = 1e-8)
  expect_equal(fit$value * 400, ref$loglik[2], tolerance = 1e-8)
  expect_equal(unname(safe_inv(-400 * fit$hess)), unname(vcov(ref)), tolerance = 1e-8)
  b <- c(0.2, -0.1, 0.4)
  an <- cox_lik(b, d$X, d$time, d$y, TRUE, 0)
  fd <- sapply(1:3, function(k) {
    e <- replace(numeric(3), k, 1e-6)
    (cox_lik(b + e, d$X, d$time, d$y, TRUE, 0)$grad - cox_lik(b - e, d$X, d$time, d$y, TRUE, 0)$grad) / 2e-6
  })
  expect_equal(unname(an$hess), unname(fd), tolerance = 1e-6)
})

# Line-by-line port of V_k() from the DiSAH code (fun.R), without tidyverse.
disah_V_k <- function(yobs, delta, X) {
  N <- length(yobs); p <- ncol(X)
  y_sort <- sort(yobs); y_rank <- rank(yobs, ties.method = "min")
  X_bar <- matrix(0, N, p)
  for (j in 1:N) { Y <- as.numeric(yobs >= y_sort[j]); X_bar[j, ] <- colSums(X * Y) / sum(Y) }
  dd <- diff(c(0, y_sort))
  A0 <- B0 <- matrix(0, p, p); d0 <- rep(0, p)
  for (i in 1:N) {
    Ki <- y_rank[i]
    Ri <- matrix(rep(X[i, ], Ki), nrow = Ki, byrow = TRUE) - X_bar[1:Ki, , drop = FALSE]
    Rdi <- Ri * sqrt(dd[1:Ki])
    A0 <- A0 + t(Rdi) %*% Rdi
    d0 <- d0 + (X[i, ] - X_bar[y_rank[i], ]) * delta[i]
    B0 <- B0 + delta[i] * (Ri[Ki, ] %*% t(Ri[Ki, ]))
  }
  list(A = A0, d = d0, B = B0)
}

test_that("the additive-hazards sums match the DiSAH reference code", {
  d <- sim_surv(250, c(0.7, -0.4, 0.3), 2)
  ref <- disah_V_k(d$time, d$y, d$X)
  got <- linying_stats(d$X, d$time, d$y)
  expect_equal(unname(got$A), unname(ref$A), tolerance = 1e-10)
  expect_equal(unname(got$d), unname(ref$d), tolerance = 1e-10)
  expect_equal(unname(got$B), unname(ref$B), tolerance = 1e-10)
  # Combining two sites equals solving the summed (stratified) equations.
  e <- sim_surv(150, c(0.7, -0.4, 0.3), 3)
  ref2 <- disah_V_k(e$time, e$y, e$X)
  fit <- linying_combine(list(got, linying_stats(e$X, e$time, e$y)))
  expect_equal(unname(fit$beta), unname(drop(solve(ref$A + ref2$A, ref$d + ref2$d))), tolerance = 1e-10)
  expect_true(all(fit$se > 0))
})

test_that("survival metrics match AutoScore-Survival", {
  set.seed(4)
  n <- 500
  status <- rbinom(n, 1, 0.7)
  time <- round(rexp(n, 0.05), 1) + 0.1
  score <- rnorm(n) - time / 20
  r <- evaluate_survival(score, time, status, c(5, 10, 20))
  val <- data.frame(label_time = time, label_status = status)
  expect_equal(r$metrics$iAUC, AutoScore:::eva_performance_iauc(score, val, print = FALSE), tolerance = 1e-12)
  ref_c <- survival::concordancefit(survival::Surv(time, status), -score)$concordance
  expect_equal(r$metrics$C_index, unname(ref_c))
  expect_gt(r$metrics$C_index, 0.6)
  expect_length(r$auc_t$auc, 3)
})

# Imputation ----------------------------------------------------------------------

test_that("one Newton round from the local average gives the pooled linear imputation model", {
  set.seed(5)
  make <- function(n, shift) {
    d <- data.frame(x1 = rnorm(n, shift), x2 = rnorm(n), g = factor(sample(c("a", "b", "c"), n, TRUE)))
    d$x3 <- 1 + 0.5 * d$x1 - d$x2 + (d$g == "b") + rnorm(n, sd = 0.5)
    d$x3[runif(n) < 0.3] <- NA
    d$x1[runif(n) < 0.1] <- NA
    d
  }
  sites <- list(A = make(400, 0), B = make(300, 1))
  levels <- list(g = c("a", "b", "c"))
  preds <- c("x1", "x2", "g", "x3")
  req <- list(kind = "impute_local", predictors = preds, levels = levels, maxit = 2, seed = 1, min_count = 10,
              models = list(list(var = "x3", family = "gaussian")))
  r1 <- lapply(sites, function(d) impute_reply(req, list(train = d)))
  b0 <- impute_combine(req, r1, 1)
  req2 <- req
  req2$kind <- "impute_derivs"
  req2$models <- list(list(var = "x3", family = "gaussian", coef = b0))
  r2 <- lapply(sites, function(d) impute_reply(req2, list(train = d)))
  b1 <- impute_combine(req2, r2, 1)
  # Pooled least squares on the same (locally completed) inputs.
  rows <- lapply(sites, function(d) impute_rows(d, local_fill(d[, preds], 2, 1), req$models[[1]], req))
  X <- do.call(rbind, lapply(rows, `[[`, "X")); y <- unlist(lapply(rows, `[[`, "y"))
  expect_equal(unname(b1), unname(drop(solve(crossprod(X), crossprod(X, y)))), tolerance = 1e-6)

  res <- list(predictors = preds, levels = levels, maxit = 2, seed = 1, calibrate = FALSE,
              models = list(list(var = "x3", family = "gaussian", coef = b1)))
  data <- list(train = sites$A, validate = sites$A, test = sites$B)
  out <- apply_imputation(data, res)
  expect_false(anyNA(out$train$x3))
  expect_true(anyNA(out$train$x1))                     # not chosen for imputation: left as is
  expect_equal(out$train$x3[!is.na(sites$A$x3)], sites$A$x3[!is.na(sites$A$x3)])
  expect_equal(attr(out, "imputed")$train$x3, sum(is.na(sites$A$x3)))
  expect_identical(apply_imputation(data, res)$test$x3, out$test$x3)   # deterministic
})
