# Loading and validating site data -------------------------------------------
#
# Each site is a list with `train`, `validate` and `test` data frames that share
# a binary `label` column. A single data frame is also accepted and is split
# 70/10/20 at load time.

SPLIT_ALIASES <- list(
  train    = c("train", "train_set", "training"),
  validate = c("validate", "validation", "validation_set", "valid"),
  test     = c("test", "test_set", "testing")
)

# Coerce a label vector (logical, factor, character or numeric) to integer 0/1.
as_binary_label <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.factor(x) || is.character(x)) {
    lv <- if (is.factor(x)) levels(droplevels(as.factor(x))) else sort(unique(na.omit(x)))
    if (length(lv) != 2) stop("`label` must have exactly two classes, found ", length(lv), ".")
    positive <- if (all(lv %in% c("0", "1"))) "1" else if (all(toupper(lv) %in% c("FALSE", "TRUE"))) lv[toupper(lv) == "TRUE"] else lv[2]
    return(as.integer(as.character(x) == positive))
  }
  vals <- sort(unique(na.omit(x)))
  if (!all(vals %in% c(0, 1))) stop("Numeric `label` must be coded 0/1.")
  as.integer(x)
}

clean_split <- function(df, split_name) {
  if (!is.data.frame(df)) stop("`", split_name, "` is not a data frame.")
  df <- as.data.frame(df)
  if (!"label" %in% names(df)) stop("`", split_name, "` has no `label` column.")
  df$label <- as_binary_label(df$label)
  if (anyNA(df$label)) stop("`label` in `", split_name, "` contains missing values.")
  for (v in setdiff(names(df), "label")) {
    if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  }
  df
}

split_site_data <- function(df, ratio = c(0.7, 0.1, 0.2), seed = 1234) {
  set.seed(seed)
  idx <- sample(c("train", "validate", "test"), nrow(df), replace = TRUE, prob = ratio)
  lapply(c(train = "train", validate = "validate", test = "test"), function(s) df[idx == s, , drop = FALSE])
}

# Normalise any supported object into list(train, validate, test).
normalize_site_data <- function(obj, seed = 1234) {
  if (is.list(obj) && !is.data.frame(obj) && !is.null(obj$non_cv)) obj <- obj$non_cv
  if (is.data.frame(obj)) obj <- split_site_data(obj, seed = seed)
  if (!is.list(obj)) stop("Expected a list of data frames or a single data frame.")
  out <- lapply(names(SPLIT_ALIASES), function(s) {
    hit <- intersect(SPLIT_ALIASES[[s]], names(obj))
    if (length(hit) == 0) stop("Missing the `", s, "` split (accepted names: ",
                               paste(SPLIT_ALIASES[[s]], collapse = ", "), ").")
    clean_split(obj[[hit[1]]], s)
  })
  names(out) <- names(SPLIT_ALIASES)
  out
}

read_site_file <- function(path, seed = 1234) {
  normalize_site_data(readRDS(path), seed = seed)
}

# Check that sites can be federated: shared predictors, both classes present.
# Returns a list of character vectors: `errors` block the analysis, `warnings` do not.
check_sites <- function(sites) {
  errors <- character(0)
  warnings <- character(0)
  if (length(sites) < 2) errors <- c(errors, "Select at least two sites.")
  preds <- lapply(sites, function(s) setdiff(names(s$train), "label"))
  common <- Reduce(intersect, preds)
  if (length(common) == 0) errors <- c(errors, "The sites share no predictor columns.")
  for (nm in names(sites)) {
    extra <- setdiff(preds[[nm]], common)
    if (length(extra)) warnings <- c(warnings, sprintf("%s: ignoring columns not present at every site (%s).",
                                                       nm, paste(extra, collapse = ", ")))
    for (s in names(SPLIT_ALIASES)) {
      if (length(unique(sites[[nm]][[s]]$label)) < 2)
        errors <- c(errors, sprintf("%s: the %s set needs both outcome classes.", nm, s))
    }
  }
  list(errors = errors, warnings = warnings, predictors = common)
}

site_summary <- function(sites) {
  do.call(rbind, lapply(names(sites), function(nm) {
    s <- sites[[nm]]
    data.frame(
      Site = nm,
      Train = nrow(s$train), Validation = nrow(s$validate), Test = nrow(s$test),
      Predictors = ncol(s$train) - 1,
      `Event rate (train)` = sprintf("%.1f%%", 100 * mean(s$train$label)),
      check.names = FALSE
    )
  }))
}
