# Loading and validating site data -------------------------------------------
#
# Each site is a list with `train`, `validate` and `test` data frames. The
# outcome is either binary (a `label` column) or survival (`label_time` and
# `label_status` columns, as in AutoScore-Survival). A single data frame is
# also accepted and is split 70/10/20 at load time.

SPLIT_ALIASES <- list(
  train    = c("train", "train_set", "training"),
  validate = c("validate", "validation", "validation_set", "valid"),
  test     = c("test", "test_set", "testing")
)

OUTCOME_COLUMNS <- c("label", "label_time", "label_status")

# "survival" when the data has label_time and label_status, "binary" when it
# has label.
data_outcome <- function(df) {
  if (all(c("label_time", "label_status") %in% names(df))) return("survival")
  if ("label" %in% names(df)) return("binary")
  stop("Expected a binary `label` column, or `label_time` and `label_status` columns for a survival outcome.")
}

predictor_names <- function(df) setdiff(names(df), OUTCOME_COLUMNS)

# 0/1 event indicator of a split: the binary label or the survival status.
event_indicator <- function(df) if ("label_status" %in% names(df)) df$label_status else df$label

# Coerce a label vector (logical, factor, character or numeric) to integer 0/1.
as_binary_label <- function(x, name = "label") {
  if (is.logical(x)) return(as.integer(x))
  if (is.factor(x) || is.character(x)) {
    lv <- if (is.factor(x)) levels(droplevels(as.factor(x))) else sort(unique(na.omit(x)))
    if (length(lv) != 2) stop("`", name, "` must have exactly two classes, found ", length(lv), ".")
    positive <- if (all(lv %in% c("0", "1"))) "1" else if (all(toupper(lv) %in% c("FALSE", "TRUE"))) lv[toupper(lv) == "TRUE"] else lv[2]
    return(as.integer(as.character(x) == positive))
  }
  vals <- sort(unique(na.omit(x)))
  if (!all(vals %in% c(0, 1))) stop("Numeric `", name, "` must be coded 0/1.")
  as.integer(x)
}

clean_split <- function(df, split_name, outcome) {
  if (!is.data.frame(df)) stop("`", split_name, "` is not a data frame.")
  df <- as.data.frame(df)
  if (data_outcome(df) != outcome) stop("`", split_name, "` does not have the same outcome columns as the training set.")
  if (outcome == "binary") {
    df$label <- as_binary_label(df$label)
    if (anyNA(df$label)) stop("`label` in `", split_name, "` contains missing values.")
  } else {
    df$label_status <- as_binary_label(df$label_status, "label_status")
    df$label_time <- suppressWarnings(as.numeric(df$label_time))
    if (anyNA(df$label_status) || anyNA(df$label_time))
      stop("`label_time` or `label_status` in `", split_name, "` contains missing or non-numeric values.")
    if (any(df$label_time < 0)) stop("`label_time` in `", split_name, "` must not be negative.")
    df$label <- NULL
  }
  for (v in predictor_names(df)) {
    if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  }
  df
}

split_site_data <- function(df, ratio = c(0.7, 0.1, 0.2), seed = 1234) {
  set.seed(seed)
  idx <- sample(c("train", "validate", "test"), nrow(df), replace = TRUE, prob = ratio)
  lapply(c(train = "train", validate = "validate", test = "test"), function(s) df[idx == s, , drop = FALSE])
}

# Normalise any supported object into list(train, validate, test). The
# outcome type is stored as attribute "outcome".
normalize_site_data <- function(obj, seed = 1234) {
  if (is.list(obj) && !is.data.frame(obj) && !is.null(obj$non_cv)) obj <- obj$non_cv
  if (is.data.frame(obj)) obj <- split_site_data(obj, seed = seed)
  if (!is.list(obj)) stop("Expected a list of data frames or a single data frame.")
  hits <- lapply(SPLIT_ALIASES, function(a) intersect(a, names(obj)))
  missing <- names(hits)[lengths(hits) == 0]
  if (length(missing)) stop("Missing the `", missing[1], "` split (accepted names: ",
                            paste(SPLIT_ALIASES[[missing[1]]], collapse = ", "), ").")
  if (!is.data.frame(obj[[hits$train[1]]])) stop("`train` is not a data frame.")
  outcome <- data_outcome(obj[[hits$train[1]]])
  out <- lapply(names(SPLIT_ALIASES), function(s) clean_split(obj[[hits[[s]][1]]], s, outcome))
  names(out) <- names(SPLIT_ALIASES)
  attr(out, "outcome") <- outcome
  out
}

read_site_file <- function(path, seed = 1234) {
  normalize_site_data(readRDS(path), seed = seed)
}
