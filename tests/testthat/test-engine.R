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
    total <- aggregate_derivs(Map(function(X, y) site_derivs("ODAL2", bb, X, y), Xs, ys), w)
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
