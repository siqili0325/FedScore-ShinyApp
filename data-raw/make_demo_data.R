# How the demo data in data/ was made. Run from this repository's root with the
# path to a checkout of the FedScore repository:
#   Rscript data-raw/make_demo_data.R /path/to/FedScore
# site1.rds and site2.rds (binary outcome) are data/binary/site{1,2}.rds$non_cv
# from that repository, unchanged.

args <- commandArgs(TRUE)
fedscore <- if (length(args)) args[1] else "../FedScore_update"

# Survival outcome: the two heterogeneous sites of data/survival/demo_hete.rds.
# Their training sets are split 7:1 into train and validate (the test sets are
# kept), matching the 70/10/20 layout of the binary demo.
hete <- readRDS(file.path(fedscore, "data", "survival", "demo_hete.rds"))
set.seed(2024)
for (k in seq_along(hete)) {
  tr <- hete[[k]]$all_train
  val <- sample(nrow(tr), round(nrow(tr) / 8))
  site <- list(train = tr[-val, ], validate = tr[val, ], test = hete[[k]]$test)
  site <- lapply(site, function(d) { rownames(d) <- NULL; d })
  saveRDS(site, sprintf("data/surv_site%d.rds", k))
}

# Missing values: the binary demo with values removed at random, so that the
# sites differ in how much is missing (Lab_H, the top-ranked variable, is
# missing more often for older patients at site 1).
set.seed(2025)
for (k in 1:2) {
  site <- readRDS(file.path(fedscore, "data", "binary", sprintf("site%d.rds", k)))$non_cv
  site <- lapply(site, function(d) {
    n <- nrow(d)
    if (k == 1) {
      p <- plogis(-1.5 + 0.03 * (d$Age - mean(d$Age)))   # about 20%, more at older ages
      d$Lab_H[runif(n) < p] <- NA
      d$Lab_B[runif(n) < 0.10] <- NA
    } else {
      d$Lab_H[runif(n) < 0.05] <- NA
      d$Lab_B[runif(n) < 0.30] <- NA
      d$Vital_E[runif(n) < 0.15] <- NA
    }
    d
  })
  saveRDS(site, sprintf("data/missing_site%d.rds", k))
}
