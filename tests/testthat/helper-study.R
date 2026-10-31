# Shared by the tests: load the app's functions and simulate sites taking part
# in a study through a real store.

app_dir <- normalizePath(file.path(testthat::test_path(), "..", ".."))
for (f in list.files(file.path(app_dir, "R"), full.names = TRUE)) source(f, local = TRUE)

study_state <- function(store, study) parse_study(store_snapshot(store, study), study)

# `site` answers its open request, if any. Returns TRUE if it wrote a reply.
respond_as <- function(store, study, site, data) {
  t <- site_task(study_state(store, study), site)
  if (is.null(t)) return(FALSE)
  store$write(round_file(study, t$step, t$round, site), site_reply(t$request, site, data))
  TRUE
}

# Run the open step to completion: every site replies, the lead combines.
finish_step <- function(store, study, sites, lead) {
  repeat {
    for (s in names(sites)) respond_as(store, study, s, sites[[s]])
    st <- study_state(store, study)
    lt <- lead_task(st)
    if (is.null(lt)) return(invisible(st))
    if (!lt$ready) stop("Step cannot continue: ", paste(c(lt$missing, lt$errors), collapse = "; "))
    out <- advance_step(st, lt$step, sites[[lead]])
    store$write(out$path, out$obj)
  }
}

run_step <- function(store, study, sites, lead, step, make_request) {
  st <- study_state(store, study)
  start_step(store, study, names(store_snapshot(store, study)), step, make_request(st), lead)
  finish_step(store, study, sites, lead)
}

# Benchmarks (not federated): the same models fitted on all sites' rows at once.
fit_pooled_logistic <- function(Xs, ys, ridge = 1e-6) {
  fit <- logistic_mle(do.call(rbind, Xs), unlist(ys), ridge)
  setNames(fit$beta, colnames(Xs[[1]]))
}

fit_pooled_pairwise <- function(Xs, ys, weights, ridge = 1e-6) {
  pats <- Map(collapse_patterns, Xs, ys)
  f <- function(b, d) {
    parts <- lapply(pats, function(pt) pairwise_lik(b, pt, d, FALSE, ridge))
    out <- list(value = sum(weights * vapply(parts, `[[`, numeric(1), "value")))
    if (d) {
      out$grad <- Reduce(`+`, Map(function(x, w) w * x$grad, parts, weights))
      out$hess <- Reduce(`+`, Map(function(x, w) w * x$hess, parts, weights))
    }
    out
  }
  setNames(newton_max(f, rep(0, ncol(Xs[[1]])))$beta, colnames(Xs[[1]]))
}

pooled_fit <- function(sites, cutoffs, vars, method, weights) {
  xy <- lapply(sites, function(s) site_xy(s$train, cutoffs, vars, method == "ODAL2"))
  Xs <- lapply(xy, `[[`, "X")
  ys <- lapply(xy, `[[`, "y")
  if (method == "ODAL2") fit_pooled_logistic(Xs, ys) else fit_pooled_pairwise(Xs, ys, weights[names(sites)])
}
