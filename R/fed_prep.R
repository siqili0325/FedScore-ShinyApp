# FedScore Modules 1-2: federated variable ranking and unified cut-offs --------
#
# Every function named `site_*` runs on one site's own data and returns only
# summary statistics; functions named `federate_*` combine those summaries.

site_weights <- function(n, scheme = c("size", "equal")) {
  scheme <- match.arg(scheme)
  w <- if (scheme == "size") n else rep(1, length(n))
  w / sum(w)
}

# Module 1 --------------------------------------------------------------------

# Variable importance at one site (AutoScore Module 1). Returns a named vector,
# most important first.
site_rank <- function(site, predictors, method = c("rf", "auc"), ntree = 100, seed = 1234) {
  method <- match.arg(method)
  set.seed(seed)
  train <- site$train[, c(predictors, "label")]
  validate <- if (method == "auc") site$validate[, c(predictors, "label")] else NULL
  # AutoScore_rank prints and draws a plot; keep both out of the Shiny session.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  utils::capture.output(
    imp <- AutoScore::AutoScore_rank(train, validation_set = validate, method = method, ntree = ntree)
  )
  imp
}

# Combine site rankings into the federated ranking: the weighted mean of each
# variable's rank position across sites (lower = more important).
federate_rank <- function(site_ranks, weights) {
  vars <- names(site_ranks[[1]])
  pos <- sapply(site_ranks, function(r) match(vars, names(r)))
  rownames(pos) <- vars
  fed <- as.vector(pos %*% weights)
  out <- data.frame(Variable = vars, pos, `Federated rank` = round(fed, 2),
                    check.names = FALSE, stringsAsFactors = FALSE)
  names(out)[seq_along(site_ranks) + 1] <- paste("Rank at", names(site_ranks))
  out <- out[order(fed, pos[, 1]), ]
  rownames(out) <- NULL
  out
}

# Module 2 --------------------------------------------------------------------

DEFAULT_PROBS <- c(0.05, 0.2, 0.8, 0.95)

# Site-level cut-off statistics: quantiles for continuous variables, the
# observed levels for categorical ones.
site_cut_stats <- function(train, predictors, probs = DEFAULT_PROBS) {
  stats <- lapply(predictors, function(v) {
    x <- train[[v]]
    has_na <- anyNA(x)
    if (is.factor(x)) {
      return(list(type = "categorical", levels = levels(droplevels(x)), has_na = has_na))
    }
    ux <- sort(unique(na.omit(x)))
    if (length(ux) <= 2) {
      return(list(type = "categorical", levels = format_num(ux), has_na = has_na))
    }
    q <- unname(stats::quantile(x, probs, na.rm = TRUE))
    # Flags (not values) for quantiles equal to the site's min/max: such a
    # cut-off would create an empty category, so AutoScore drops it.
    list(type = "numeric", quantiles = q, at_min = q <= ux[1], at_max = q >= ux[length(ux)], has_na = has_na)
  })
  names(stats) <- predictors
  stats
}

# Unified cut-offs: weighted mean of site quantiles (rounded to 3 significant
# digits, as AutoScore does); union of levels for categorical variables.
federate_cutoffs <- function(site_stats, weights) {
  vars <- names(site_stats[[1]])
  out <- lapply(vars, function(v) {
    per_site <- lapply(site_stats, `[[`, v)
    has_na <- any(vapply(per_site, `[[`, logical(1), "has_na"))
    is_num <- vapply(per_site, function(s) s$type == "numeric", logical(1))
    if (!any(is_num)) {
      lv <- unique(unlist(lapply(per_site, `[[`, "levels")))
      if (all(!is.na(suppressWarnings(as.numeric(lv))))) lv <- lv[order(as.numeric(lv))]
      return(list(type = "categorical", levels = lv, has_na = has_na))
    }
    q <- sapply(per_site[is_num], `[[`, "quantiles")
    at_edge <- sapply(per_site[is_num], function(s) s$at_min | s$at_max)
    keep <- !apply(matrix(at_edge, nrow = nrow(as.matrix(q))), 1, all)
    if (!any(keep)) keep[] <- TRUE
    w <- weights[is_num] / sum(weights[is_num])
    cuts <- unique(signif(as.vector(as.matrix(q)[keep, , drop = FALSE] %*% w), 3))
    list(type = "numeric", cuts = sort(cuts), has_na = has_na)
  })
  names(out) <- vars
  out
}

format_num <- function(x) format(x, trim = TRUE, drop0trailing = TRUE, scientific = FALSE)

# Category labels in AutoScore style: "(,c1)", "[c1,c2)", ..., "[ck,)".
cut_labels <- function(cuts) {
  s <- format_num(cuts)
  k <- length(s)
  c(paste0("(,", s[1], ")"),
    if (k > 1) paste0("[", s[-k], ",", s[-1], ")"),
    paste0("[", s[k], ",)"))
}

var_levels <- function(cut) {
  lv <- if (cut$type == "numeric") cut_labels(cut$cuts) else cut$levels
  if (cut$has_na) c(lv, "Unknown") else lv
}

# Categorise `df[vars]` with the unified cut-offs. Every site gets the same
# factor levels, so design matrices line up across sites. Values that cannot
# be placed (a categorical level never seen in training) become NA.
categorize <- function(df, cutoffs, vars = names(cutoffs)) {
  out <- lapply(vars, function(v) {
    cut <- cutoffs[[v]]
    x <- df[[v]]
    lv <- var_levels(cut)
    if (cut$type == "numeric") {
      idx <- findInterval(x, cut$cuts) + 1L
      y <- lv[idx]
    } else {
      y <- as.character(if (is.factor(x)) x else format_num(x))
      y[!is.na(y) & !y %in% lv] <- NA
    }
    if (cut$has_na) y[is.na(x)] <- "Unknown"
    factor(y, levels = lv)
  })
  names(out) <- vars
  as.data.frame(out, check.names = FALSE)
}

cutoff_table <- function(cutoffs, order = names(cutoffs)) {
  data.frame(
    Variable = order,
    Type = vapply(cutoffs[order], `[[`, character(1), "type"),
    `Cut-offs` = vapply(cutoffs[order], function(c) {
      if (c$type == "numeric") paste(format_num(c$cuts), collapse = ", ") else paste(c$levels, collapse = " | ")
    }, character(1)),
    Categories = vapply(cutoffs[order], function(c) paste(var_levels(c), collapse = "  "), character(1)),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# Parse a user-edited cut-off string such as "90, 100, 120".
parse_cuts <- function(text) {
  vals <- suppressWarnings(as.numeric(trimws(strsplit(text, "[,;[:space:]]+")[[1]])))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) stop("Enter at least one numeric cut-off.")
  sort(unique(vals))
}
