# FedScore Modules 4-6: score tables, evaluation, and TLScore ------------------

# Turn coefficients into an AutoScore-style points table. Within each variable
# the lowest-risk category gets 0 points; points are then scaled so the
# maximum attainable total is `max_score`. `map` comes from design_terms().
score_table <- function(beta, map, max_score = 100) {
  levels <- attr(map, "levels")
  map <- map[map$variable != "(Intercept)", ]
  coef <- beta[map$column]
  vars <- names(levels)
  rows <- lapply(vars, function(v) {
    lv <- levels[[v]]
    data.frame(Variable = v, Category = lv,
               coef = c(0, unname(coef[map$variable == v]))[seq_along(lv)],
               stringsAsFactors = FALSE)
  })
  tbl <- do.call(rbind, rows)
  tbl$shift <- ave(tbl$coef, tbl$Variable, FUN = function(x) x - min(x))
  total <- sum(tapply(tbl$shift, tbl$Variable, max))
  tbl$Points <- if (total > 0) round(tbl$shift * max_score / total) else 0
  tbl[, c("Variable", "Category", "Points")]
}

# Total score for each row of `df` (raw data; categorised here).
compute_scores <- function(df, tbl, cutoffs) {
  vars <- unique(tbl$Variable)
  cat_df <- categorize(df, cutoffs, vars)
  total <- numeric(nrow(df))
  for (v in vars) {
    pts <- setNames(tbl$Points[tbl$Variable == v], tbl$Category[tbl$Variable == v])
    p <- pts[as.character(cat_df[[v]])]
    p[is.na(p)] <- 0
    total <- total + p
  }
  unname(total)
}

score_roc <- function(score, label) {
  pROC::roc(response = label, predictor = score, levels = c(0, 1), direction = "<", quiet = TRUE)
}

# AUC with DeLong CI plus sensitivity/specificity at the Youden threshold.
evaluate_scores <- function(score, label) {
  r <- score_roc(score, label)
  ci <- suppressWarnings(as.numeric(pROC::ci.auc(r, method = "delong")))
  best <- pROC::coords(r, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"),
                       transpose = FALSE)[1, ]
  data.frame(N = length(label), Events = sum(label), AUC = ci[2], Lower = ci[1], Upper = ci[3],
             Threshold = best$threshold, Sensitivity = best$sensitivity, Specificity = best$specificity)
}

auc_only <- function(score, label) as.numeric(pROC::auc(score_roc(score, label)))

# ROC curve as (1 - specificity, sensitivity) pairs at each distinct score.
roc_points <- function(score, label) {
  r <- score_roc(score, label)
  o <- order(1 - r$specificities, r$sensitivities)
  list(fpr = (1 - r$specificities)[o], tpr = r$sensitivities[o])
}

# Survival outcomes (as AutoScore-Survival) -------------------------------------

# Integrated AUC: Uno's time-dependent AUC at each event time, weighted by the
# drop in the Kaplan-Meier curve (AutoScore's eva_performance_iauc).
surv_iauc <- function(score, time, status) {
  y <- survival::Surv(time, status)
  km <- summary(survival::survfit(y ~ 1))
  keep <- seq_len(max(0, length(km$time) - 1))
  if (length(keep) < 1) return(NA_real_)
  auc <- survAUC::AUC.uno(y, y, lpnew = score, times = km$time[keep])$auc
  drop <- c(1, km$surv[keep][-length(keep)]) - km$surv[keep]
  sum(drop / sum(drop) * auc)
}

# iAUC, Harrell's C-index with a 95% CI, and AUC(t) at `times`.
evaluate_survival <- function(score, time, status, times) {
  y <- survival::Surv(time, status)
  cf <- survival::concordancefit(y, -score)
  se <- sqrt(cf$var[1])
  times <- times[times < max(time[status == 1])]
  auc_t <- if (length(times)) survAUC::AUC.uno(y, y, lpnew = score, times = times)$auc else numeric(0)
  list(metrics = data.frame(N = length(time), Events = sum(status), iAUC = surv_iauc(score, time, status),
                            C_index = cf$concordance, Lower = cf$concordance - 1.96 * se,
                            Upper = cf$concordance + 1.96 * se),
       auc_t = list(time = times, auc = auc_t))
}

# Mean validation AUC across sites, weighted by the site weights.
summarise_parsimony <- function(res, weights) {
  w <- weights[res$Site]
  agg <- aggregate(cbind(AUC = res$AUC * w, W = w) ~ Model + Variable, data = transform(res, w = w), FUN = sum)
  agg$AUC <- agg$AUC / agg$W
  agg[order(agg$Model), c("Model", "Variable", "AUC")]
}

# TLScore: adapt the federated model to one target site with a lasso-penalised
# correction fitted on that site's own training data (offset = federated
# linear predictor), i.e. beta_TL = beta_fed + delta. `beta` is named by the
# federated model's terms; only the target site's data is used. Available for
# the logistic models and for Cox (ODACH); not for additive hazards (DiSAH).
TL_METHODS <- c("ODAL2", "dCLR", "ODACH")

fit_tlscore <- function(beta, method, vars, site, cutoffs, max_score = 100, seed = 2024, nfolds = 5) {
  if (!method %in% TL_METHODS) stop("TLScore is not available for ", method, " models.")
  cat_df <- categorize(site$train, cutoffs, vars)
  X <- design_matrix(cat_df, intercept = FALSE)
  beta_fed <- beta[colnames(X)]
  eta <- drop(X %*% beta_fed)
  set.seed(seed)
  cv <- if (method == "ODACH") {
    # glmnet's Cox model rejects zero times; move them just above zero
    # (their order relative to the other times is unchanged).
    time <- site$train$label_time
    time[time <= 0] <- min(time[time > 0]) / 2
    y <- survival::Surv(time, site$train$label_status)
    glmnet::cv.glmnet(X, y, family = "cox", offset = eta, alpha = 1, nfolds = nfolds, maxit = 1e5)
  } else {
    # Recentre the offset on the site's event rate (dCLR has no intercept, and
    # ODAL2's comes from other sites); glmnet's own intercept absorbs the rest.
    offset <- eta - mean(eta) + stats::qlogis(mean(site$train$label))
    glmnet::cv.glmnet(X, site$train$label, family = "binomial", offset = offset,
                      alpha = 1, nfolds = nfolds, maxit = 1e5)
  }
  delta <- as.numeric(stats::coef(cv, s = "lambda.min"))
  if (method != "ODACH") delta <- delta[-1]           # drop glmnet's intercept
  names(delta) <- colnames(X)
  beta_tl <- beta
  beta_tl[names(delta)] <- beta_fed + delta
  list(beta = beta_tl, delta = delta, lambda = cv$lambda.min,
       table = score_table(beta_tl, design_terms(cutoffs, vars, method == "ODAL2"), max_score))
}
