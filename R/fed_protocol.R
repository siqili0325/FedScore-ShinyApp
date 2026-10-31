# FedScore exchange protocol -----------------------------------------------------
#
# A study is one folder in the exchange store (see fed_store.R):
#   study.json                    settings, written by the lead site
#   <step>/round<k>.json          request k of a step, written by the lead site
#   <step>/round<k>_<site>.json   each site's reply to request k
#   <step>/result.json            outcome of the step, written by the lead site
# The lead site opens the steps in order. Every site, the lead included,
# answers each request from its own data only; the lead combines the replies
# into the next request or the step's result. Each request carries the id of
# its step run (`step_id`), so files left over from an earlier run are ignored.

PROTOCOL_VERSION <- 1L
STEPS <- c(summary = "Data summary", rank = "Ranking", cut = "Cut-offs",
           pars = "Model selection", final = "Final score")

now_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
new_id <- function(step) paste0(step, "-", format(Sys.time(), "%Y%m%dT%H%M%OS6", tz = "UTC"))

study_file <- function(study, ...) paste(c(study, ...), collapse = "/")
round_file <- function(study, step, k, site = NULL) {
  study_file(study, step, if (is.null(site)) sprintf("round%d.json", k) else sprintf("round%d_%s.json", k, site))
}

# What each kind of request asks sites to share (shown in the app).
shared_content <- function(req) {
  switch(req$kind,
         summary = "sample sizes, event counts and predictor names",
         rank = "rank position of each variable",
         cut = "quantiles of each continuous variable and the levels of categorical ones",
         local = "local estimate and its variance, for each model",
         derivs = "gradient and Hessian at the broadcast estimate, for each model",
         surrogate = "surrogate estimate and its variance, for each model",
         evaluate = if (isTRUE(req$full)) "test-set AUC with CI, operating point and ROC curve"
                    else "validation AUC of each model",
         req$kind)
}

# Reading a study ---------------------------------------------------------------

# Organise a store snapshot (path -> parsed JSON) into the study settings and,
# for each step, its current rounds (request plus replies so far) and result.
parse_study <- function(snap, study) {
  get <- function(...) snap[[study_file(study, ...)]]
  config <- get("study.json")
  if (is.null(config)) return(NULL)
  config$sites <- as.character(unlist(config$sites))
  steps <- lapply(names(STEPS), function(st) {
    first <- get(st, "round1.json")
    id <- first$step_id
    rounds <- list()
    k <- 1
    while (!is.null(id)) {
      req <- get(st, sprintf("round%d.json", k))
      if (is.null(req) || !identical(req$step_id, id)) break
      replies <- list()
      for (s in config$sites) {
        r <- get(st, sprintf("round%d_%s.json", k, s))
        if (!is.null(r) && identical(r$request, req$id)) replies[[s]] <- r
      }
      rounds[[k]] <- list(request = req, replies = replies)
      k <- k + 1
    }
    result <- get(st, "result.json")
    if (is.null(id) || !identical(result$step_id, id)) result <- NULL
    list(id = id, rounds = rounds, result = result)
  })
  names(steps) <- names(STEPS)
  list(config = config, steps = steps)
}

# The step that is waiting for replies or for the lead, if any.
open_step <- function(state) {
  for (st in names(STEPS)) {
    s <- state$steps[[st]]
    if (length(s$rounds) && is.null(s$result)) return(st)
  }
  NULL
}

# The request `site` still has to answer. A request it answered with an error
# is returned with `retry = TRUE`.
site_task <- function(state, site) {
  st <- open_step(state)
  if (is.null(st)) return(NULL)
  rounds <- state$steps[[st]]$rounds
  last <- rounds[[length(rounds)]]
  r <- last$replies[[site]]
  if (!is.null(r) && is.null(r$error)) return(NULL)
  list(step = st, round = length(rounds), request = last$request, retry = !is.null(r))
}

# What the lead is waiting for in the open step.
lead_task <- function(state) {
  st <- open_step(state)
  if (is.null(st)) return(NULL)
  rounds <- state$steps[[st]]$rounds
  last <- rounds[[length(rounds)]]
  errors <- Filter(function(r) !is.null(r$error), last$replies)
  missing <- setdiff(state$config$sites, names(last$replies))
  list(step = st, round = length(rounds), request = last$request, missing = missing,
       errors = vapply(errors, `[[`, "", "error"), ready = !length(missing) && !length(errors))
}

study_weights <- function(state) {
  w <- unlist(state$steps$summary$result$weights)
  w[state$config$sites]
}

# Lead site: creating a study and opening steps ---------------------------------

# `min_count`: the smallest number of patients a site may share statistics
# about (see check_min_count); 0 switches the check off.
create_study <- function(store, study, lead, sites, weighting = "size", min_count = 10) {
  if (!valid_name(study)) stop("Study names may contain only letters, digits, '-' and '_'.")
  bad <- sites[!vapply(sites, valid_name, logical(1))]
  if (length(bad)) stop("Invalid site name(s): ", paste(bad, collapse = ", "))
  if (anyDuplicated(sites)) stop("Site names must be unique.")
  if (length(sites) < 2) stop("A study needs at least two sites.")
  if (!lead %in% sites) stop("The lead site must be one of the sites.")
  if (!is.numeric(min_count) || length(min_count) != 1 || is.na(min_count) || min_count < 0)
    stop("The minimum count must be a number of at least 0.")
  store$write(study_file(study, "study.json"),
              list(protocol = PROTOCOL_VERSION, study = study, lead = lead, sites = sites,
                   min_count = as.integer(min_count), created = now_utc()))
  start_step(store, study, character(0), "summary", list(kind = "summary", weighting = weighting), lead)
}

# Write round 1 of `step`, first removing every file of this step and the
# later ones (their results depend on this step). `paths` are the study's
# current files.
start_step <- function(store, study, paths, step, request, lead) {
  later <- names(STEPS)[match(step, names(STEPS)):length(STEPS)]
  rel <- substring(paths, nchar(study) + 2)
  old <- paths[sub("/.*$", "", rel) %in% later & grepl("/", rel)]
  if (length(old)) store$delete(old)
  id <- new_id(step)
  req <- c(list(id = id, step_id = id, step = step, round = 1L, from = lead, created = now_utc()), request)
  store$write(round_file(study, step, 1), req)
  invisible(req)
}

# Round-1 request of each step, built from the earlier results.
rank_request <- function(state, method = "rf", ntree = 100, seed = 1234) {
  list(kind = "rank", predictors = state$steps$summary$result$predictors,
       method = method, ntree = ntree, seed = seed)
}

cut_request <- function(state, probs = DEFAULT_PROBS) {
  list(kind = "cut", predictors = state$steps$summary$result$predictors, probs = probs)
}

# Model fitting for the parsimony analysis (`models` = top 1..m variables,
# evaluated on the validation sets) or the final score (one model, test sets).
# ODAL2 starts from the lead site's own estimate, so it needs the lead's data.
fit_request <- function(state, models, method = c("ODAL2", "dCLR"), refinements = 1, split = "validate",
                        max_score = 100, lead_data = NULL) {
  method <- match.arg(method)
  cutoffs <- state$steps$cut$result$cutoffs
  intercept <- method == "ODAL2"
  specs <- lapply(models, function(vars) list(vars = vars, terms = design_terms(cutoffs, vars, intercept)$column))
  min_count <- if (is.null(state$config$min_count)) 0L else state$config$min_count
  base <- list(method = method, refinements = refinements, split = split, full = split == "test",
               max_score = max_score, cutoffs = cutoffs, weights = as.list(study_weights(state)),
               min_count = min_count)
  if (method == "dCLR") return(c(list(kind = "local"), base, list(models = specs)))
  if (is.null(lead_data)) stop("ODAL2 needs the lead site's data to start.")
  specs <- lapply(specs, function(m) {
    check_min_count(lead_data, list(kind = "derivs", cutoffs = cutoffs, min_count = min_count), m$vars)
    xy <- site_xy(lead_data$train, cutoffs, m$vars, TRUE)
    b <- unname(odal2_init(xy$X, xy$y))
    c(m, list(beta = b, init = b))
  })
  c(list(kind = "derivs", pass = 1L), base, list(models = specs))
}

# Every site: answering a request -------------------------------------------------

site_reply <- function(req, site, data) {
  body <- tryCatch(reply_body(req, data), error = function(e) list(error = conditionMessage(e)))
  c(list(request = req$id, step_id = req$step_id, site = site, created = now_utc()), body)
}

reply_body <- function(req, data) {
  switch(req$kind,
         summary = summary_reply(data),
         rank = list(ranking = names(site_rank(data, req$predictors, req$method, req$ntree, req$seed))),
         cut = list(stats = site_cut_stats(data$train, req$predictors, req$probs)),
         local = , derivs = , surrogate = , evaluate =
           list(models = lapply(req$models, function(m) model_reply(req, m, data))),
         stop("Unknown request: ", req$kind))
}

summary_reply <- function(data) {
  splits <- names(SPLIT_ALIASES)
  for (s in splits)
    if (length(unique(data[[s]]$label)) < 2) stop(sprintf("The %s set needs both outcome classes.", s))
  preds <- setdiff(names(data$train), "label")
  list(n = lapply(data[splits], nrow), events = lapply(data[splits], function(d) sum(d$label)),
       predictors = lapply(setNames(preds, preds), function(v)
         if (is.factor(data$train[[v]])) "categorical" else "numeric"))
}

# Small-cell protection. Model statistics computed from very few patients
# (a category with 3 patients, a test set with 2 events) come close to
# describing individuals, so a site refuses to share them when a non-empty
# category of the model's variables in its training data, or the events or
# non-events of the set being evaluated, number fewer than `req$min_count`.
check_min_count <- function(data, req, vars) {
  k <- if (is.null(req$min_count)) 0 else req$min_count
  if (k <= 0) return(invisible())
  if (req$kind == "evaluate") {
    y <- data[[req$split]]$label
    if (min(sum(y), sum(1 - y)) < k)
      stop(sprintf("The %s set has %d events and %d non-events; sharing its performance needs at least %d of each.",
                   req$split, sum(y), sum(1 - y), k))
    return(invisible())
  }
  counts <- lapply(categorize(data$train, req$cutoffs, vars), table)
  small <- unlist(lapply(names(counts), function(v) {
    tab <- counts[[v]]
    bad <- tab[tab > 0 & tab < k]
    if (length(bad)) sprintf("%s %s (%d)", v, names(bad), as.integer(bad))
  }))
  if (length(small))
    stop(sprintf("Too few patients to share statistics safely (minimum %d per category): %s. Merge these categories, or ask the lead site to change the cut-offs on step 3.",
                 k, paste(small, collapse = ", ")))
  invisible()
}

model_reply <- function(req, m, data) {
  check_min_count(data, req, m$vars)
  if (req$kind == "evaluate") {
    df <- data[[req$split]]
    score <- compute_scores(df, as.data.frame(m$table, stringsAsFactors = FALSE), req$cutoffs)
    if (isTRUE(req$full)) return(list(metrics = evaluate_scores(score, df$label), roc = roc_points(score, df$label)))
    return(list(n = nrow(df), auc = auc_only(score, df$label)))
  }
  xy <- site_xy(data$train, req$cutoffs, m$vars, req$method == "ODAL2")
  out <- switch(req$kind,
                local = dclr_local(xy$X, xy$y),
                derivs = c(site_derivs(req$method, m$beta, xy$X, xy$y), list(n = length(xy$y))),
                surrogate = dclr_surrogate(xy$X, xy$y, m$beta, list(grad = m$grad, hess = as.matrix(m$hess))))
  lapply(out, function(x) if (is.numeric(x) && !is.matrix(x)) unname(x) else x)
}

# Lead site: combining replies ----------------------------------------------------

# Combine the replies to the last round of `step` into the next request or
# the step's result. Returns list(path, obj) for the lead to write.
advance_step <- function(state, step, lead_data = NULL) {
  study <- state$config$study
  rounds <- state$steps[[step]]$rounds
  k <- length(rounds)
  req <- rounds[[k]]$request
  replies <- rounds[[k]]$replies[state$config$sites]
  out <- switch(step,
                summary = summary_result(req, replies),
                rank = rank_result(replies, study_weights(state)),
                cut = cut_result(replies, study_weights(state)),
                pars = , final = advance_fit(req, replies, lead_data))
  if (!is.null(out$kind)) {
    out <- c(list(id = new_id(step), step_id = req$step_id, step = step, round = k + 1L,
                  from = req$from, created = now_utc()), out)
    return(list(path = round_file(study, step, k + 1), obj = out))
  }
  list(path = study_file(study, step, "result.json"),
       obj = c(list(step_id = req$step_id, step = step, created = now_utc()), out))
}

summary_result <- function(req, replies) {
  sites <- names(replies)
  types <- lapply(replies, function(r) unlist(r$predictors))
  common <- Reduce(intersect, lapply(types, names))
  same_type <- vapply(common, function(v) length(unique(vapply(types, `[[`, "", v))) == 1, logical(1))
  warnings <- character(0)
  if (any(!same_type))
    warnings <- c(warnings, sprintf("Ignoring predictors stored as numeric at some sites and categorical at others: %s.",
                                    paste(common[!same_type], collapse = ", ")))
  common <- common[same_type]
  for (s in sites) {
    extra <- setdiff(names(types[[s]]), common)
    extra <- setdiff(extra, names(same_type)[!same_type])
    if (length(extra)) warnings <- c(warnings, sprintf("%s: ignoring columns not present at every site (%s).",
                                                       s, paste(extra, collapse = ", ")))
  }
  if (!length(common)) stop("The sites share no predictor columns.")
  n <- vapply(replies, function(r) as.numeric(r$n$train), numeric(1))
  ev <- vapply(replies, function(r) as.numeric(r$events$train), numeric(1))
  w <- setNames(site_weights(n, req$weighting), sites)
  table <- data.frame(
    Site = sites,
    Train = n,
    Validation = vapply(replies, function(r) as.numeric(r$n$validate), numeric(1)),
    Test = vapply(replies, function(r) as.numeric(r$n$test), numeric(1)),
    `Event rate (train)` = sprintf("%.1f%%", 100 * ev / n),
    Weight = sprintf("%.2f", w),
    check.names = FALSE)
  list(predictors = common, weighting = req$weighting, weights = as.list(w), table = table, warnings = warnings)
}

rank_result <- function(replies, weights) {
  ranks <- lapply(replies, function(r) setNames(seq_along(r$ranking), r$ranking))
  list(table = federate_rank(ranks, weights[names(replies)]))
}

cut_result <- function(replies, weights) {
  cutoffs <- federate_cutoffs(lapply(replies, `[[`, "stats"), weights[names(replies)])
  list(cutoffs = cutoffs, auto = cutoffs)
}

# The fitting rounds. ODAL2: derivs -> (lead fits) -> evaluate. dCLR: local ->
# derivs -> surrogate (-> derivs -> surrogate ...) -> evaluate.
advance_fit <- function(req, replies, lead_data) {
  w <- unlist(req$weights)[names(replies)]
  M <- length(req$models)
  answers <- function(m) lapply(replies, function(r) r$models[[m]])
  keep <- function(m, ...) {
    spec <- req$models[[m]]
    c(spec[intersect(c("vars", "terms", "init"), names(spec))], list(...))
  }
  # Score table of each model, ready for the sites to evaluate.
  evaluate <- function(models) {
    models <- lapply(models, function(m) {
      map <- design_terms(req$cutoffs, m$vars, req$method == "ODAL2")
      c(m, list(table = score_table(setNames(m$beta, m$terms), map, req$max_score)))
    })
    list(kind = "evaluate", models = models)
  }
  carry <- function(x) c(x, req[intersect(c("method", "refinements", "split", "full", "max_score", "cutoffs",
                                            "weights", "min_count"), names(req))])

  switch(req$kind,
    local = carry(list(kind = "derivs", pass = 1L, models = lapply(seq_len(M), function(m) {
      a <- answers(m)
      b <- iv_average(lapply(a, function(x) x$beta), lapply(a, function(x) as.matrix(x$vcov)))$beta
      keep(m, beta = b, init = b)
    }))),
    derivs = {
      totals <- lapply(seq_len(M), function(m)
        aggregate_derivs(lapply(answers(m), function(x) list(grad = x$grad, hess = as.matrix(x$hess))), w))
      if (req$method == "ODAL2") {
        if (is.null(lead_data)) stop("The lead site's data is needed to finish ODAL2.")
        N <- sum(vapply(replies, function(r) as.numeric(r$models[[1]]$n), numeric(1)))
        carry(evaluate(lapply(seq_len(M), function(m) {
          spec <- req$models[[m]]
          xy <- site_xy(lead_data$train, req$cutoffs, spec$vars, TRUE)
          fit <- odal2_surrogate(xy$X, xy$y, spec$beta, totals[[m]], N)
          keep(m, beta = unname(fit$beta), se = unname(fit$se), converged = fit$converged)
        })))
      } else {
        carry(list(kind = "surrogate", pass = req$pass, models = lapply(seq_len(M), function(m)
          keep(m, beta = req$models[[m]]$beta, grad = totals[[m]]$grad, hess = totals[[m]]$hess))))
      }
    },
    surrogate = {
      ests <- lapply(seq_len(M), function(m) {
        a <- answers(m)
        est <- iv_average(lapply(a, function(x) x$beta), lapply(a, function(x) as.matrix(x$vcov)))
        list(beta = est$beta, se = sqrt(pmax(diag(est$vcov), 0)),
             converged = all(vapply(a, function(x) isTRUE(x$converged), logical(1))))
      })
      if (req$pass < req$refinements) {
        carry(list(kind = "derivs", pass = req$pass + 1L,
                   models = lapply(seq_len(M), function(m) keep(m, beta = ests[[m]]$beta))))
      } else {
        carry(evaluate(lapply(seq_len(M), function(m)
          keep(m, beta = ests[[m]]$beta, se = ests[[m]]$se, converged = ests[[m]]$converged))))
      }
    },
    evaluate = fit_result(req, replies, w),
    stop("Unknown request: ", req$kind))
}

fit_result <- function(req, replies, w) {
  sites <- names(replies)
  base <- req[c("method", "refinements", "max_score", "cutoffs")]
  if (isTRUE(req$full)) {
    m <- req$models[[1]]
    perf <- do.call(rbind, lapply(sites, function(s)
      cbind(Site = s, as.data.frame(replies[[s]]$models[[1]]$metrics), stringsAsFactors = FALSE)))
    return(c(base, m[c("vars", "terms", "beta", "se", "init", "converged", "table")],
             list(perf = perf, roc = lapply(replies, function(r) r$models[[1]]$roc))))
  }
  auc <- do.call(rbind, lapply(seq_along(req$models), function(m) {
    vars <- req$models[[m]]$vars
    do.call(rbind, lapply(sites, function(s) data.frame(
      Model = m, Variable = vars[length(vars)], Site = s,
      N = replies[[s]]$models[[m]]$n, AUC = replies[[s]]$models[[m]]$auc, stringsAsFactors = FALSE)))
  }))
  c(base, list(models = lapply(req$models, function(m) m[setdiff(names(m), "table")]),
               auc = auc, summary = summarise_parsimony(auc, w)))
}

# Lead site: re-publish the cut-offs after fine-tuning. Later steps used the
# old cut-offs, so their files are removed.
publish_cutoffs <- function(store, state, paths, cutoffs) {
  study <- state$config$study
  res <- state$steps$cut$result
  rel <- substring(paths, nchar(study) + 2)
  old <- paths[sub("/.*$", "", rel) %in% c("pars", "final") & grepl("/", rel)]
  if (length(old)) store$delete(old)
  res$cutoffs <- cutoffs
  res$created <- now_utc()
  store$write(study_file(study, "cut", "result.json"), res)
}

# Tables arrive from JSON as lists of columns.
as_df <- function(x) as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)

count_numbers <- function(x) {
  if (is.list(x)) sum(vapply(x, count_numbers, numeric(1))) else if (is.numeric(x)) length(x) else 0
}

# Exchange log: one row per file in the study.
exchange_log <- function(snap, study) {
  if (!length(snap)) return(NULL)
  rel <- substring(names(snap), nchar(study) + 2)
  step <- ifelse(grepl("/", rel), sub("/.*$", "", rel), "")
  file <- sub("^.*/", "", rel)
  round <- suppressWarnings(as.integer(sub("^round([0-9]+).*$", "\\1", file)))
  from <- vapply(seq_along(snap), function(i) {
    obj <- snap[[i]]
    who <- c(obj[["site"]], obj[["from"]], obj[["lead"]])
    if (length(who)) who[1] else ""
  }, "")
  what <- vapply(seq_along(snap), function(i) {
    obj <- snap[[i]]
    if (file[i] == "study.json") return("study settings")
    if (file[i] == "result.json") return("result")
    if (!is.null(obj[["request"]])) return(if (!is.null(obj[["error"]])) paste("error:", obj[["error"]]) else "reply")
    paste("request:", shared_content(obj))
  }, "")
  numbers <- vapply(snap, count_numbers, numeric(1))
  out <- data.frame(Step = unname(ifelse(step == "", "", STEPS[step])), Round = round, From = from, Content = what,
                    Numbers = as.integer(numbers),
                    Written = vapply(snap, function(o) if (is.null(o[["created"]])) "" else o[["created"]], ""),
                    File = rel, stringsAsFactors = FALSE)
  out[order(match(step, c("", names(STEPS))), file == "result.json", out$Round, out$Written), ]
}
