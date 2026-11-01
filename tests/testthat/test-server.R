# End-to-end: the app runs as the lead site (site1, demo data) while site2 is
# simulated with the protocol functions, both using one folder store.

render_ok <- function(output, ids) {
  for (id in ids) testthat::expect_no_error(output[[id]], message = id)
}

for (alg in "dCLR") {
  test_that(sprintf("the lead site's app runs a whole %s study", alg), {
    root <- tempfile("store")
    site2 <- read_site_file(file.path(app_dir, "data", "site2.rds"))
    store2 <- folder_store(root)

    shiny::testServer(shiny::shinyAppDir(app_dir), {
      clicks <- 0
      # site2 answers whatever is open, then the lead's app re-reads the store.
      sync <- function(step) {
        for (i in 1:12) {
          respond_as(store2, "demo", "site2", site2)
          clicks <<- clicks + 1
          session$setInputs(refresh = clicks)
          if (!is.null(results[[step]])) return(invisible())
        }
        stop("step ", step, " did not finish")
      }

      session$setInputs(store_type = "folder", folder = root, study = "demo", me = "site1", auto = TRUE,
                        repo = "", branch = "main", gh_dir = "fedscore", token = "", demo = "site1")
      expect_equal(nrow(rv$data$train), 9800)
      session$setInputs(connect = 1)
      expect_null(config())
      expect_match(paste(output$study_panel$html), "Create study")

      session$setInputs(sites = "site1, site2", weighting = "size", min_count = 10, create = 1)
      expect_equal(config()$lead, "site1")
      expect_equal(config()$min_count, 10)
      expect_true(is_lead())
      sync("summary")
      expect_equal(results$summary$predictors, setdiff(names(site2$train), "label"))
      render_ok(output, c("site_table", "site_messages", "status_summary"))

      # Steps cannot start before their prerequisites.
      session$setInputs(var_mode = "top", n_vars = 4, max_score = 100, method = alg, refine = 2, start_final = 1)
      expect_null(state()$steps$final$id)

      session$setInputs(rank_method = "rf", ntree = 50, start_rank = 1)
      sync("rank")
      expect_equal(nrow(ranking()), 20)
      session$setInputs(quantiles = "0.05, 0.2, 0.8, 0.95", start_cut = 1)
      sync("cut")
      render_ok(output, c("rank_table", "rank_plot", "cut_table", "status_rank", "status_cut"))

      # Fine-tune one variable: the new cut-offs are published to the store.
      session$setInputs(tune_var = "Age", tune_cuts = "40, 60, 80", tune_apply = 1)
      expect_equal(results$cut$cutoffs$Age$cuts, c(40, 60, 80))
      expect_equal(store2$read("demo/cut/result.json")$cutoffs$Age$cuts, c(40, 60, 80))

      session$setInputs(max_vars = 4, start_pars = 1)
      sync("pars")
      expect_equal(nrow(as_df(results$pars$summary)), 4)
      render_ok(output, c("pars_plot", "pars_table", "status_pars"))

      session$setInputs(start_final = 2)
      sync("final")
      f <- results$final
      expect_equal(f$vars, head(ranking()$Variable, 4))
      expect_true(f$converged)
      render_ok(output, c("score_tbl", "perf_tbl", "perf_plot", "coef_tbl", "status_final"))
      expect_equal(as_df(f$perf)$Site, c("site1", "site2"))
      sites <- list(site1 = rv$data, site2 = site2)
      pooled <- pooled_fit(sites, f$cutoffs, f$vars, alg, study_weights(state()))
      expect_lt(max(abs(f$beta - pooled)), 0.01)

      session$setInputs(run_tl = 1)
      expect_equal(nrow(rv$tl$perf), 2)
      render_ok(output, c("tl_perf", "tl_tbl", "log_tbl"))
      lg <- log_df()
      expect_true(all(c("site1", "site2") %in% lg$From))
      session$setInputs(log_file = "final/result.json")
      expect_match(output$log_json, "\"perf\"")

      # Re-running the ranking clears every later step, for every site.
      session$setInputs(start_rank = 2)
      expect_null(results$cut)
      expect_null(results$final)
    })
  })
}

test_that("a site that is not the lead only answers requests", {
  root <- tempfile("store")
  store <- folder_store(root)
  site1 <- read_site_file(file.path(app_dir, "data", "site1.rds"))
  create_study(store, "demo", "site1", c("site1", "site2"))
  respond_as(store, "demo", "site1", site1)

  shiny::testServer(shiny::shinyAppDir(app_dir), {
    session$setInputs(store_type = "folder", folder = root, study = "demo", me = "site2", auto = FALSE,
                      repo = "", branch = "main", gh_dir = "fedscore", token = "")
    session$setInputs(connect = 1)
    expect_false(is_lead())
    expect_match(paste(output$status_summary$html), "Load your data")
    session$setInputs(demo = "site2")
    expect_match(paste(output$status_summary$html), "Send my reply")
    session$setInputs(do_reply = 1)
    expect_false(is.null(store$read("demo/summary/round1_site2.json")$n))
    # Only the lead may start steps.
    session$setInputs(rank_method = "rf", ntree = 10, start_rank = 1)
    expect_null(state()$steps$rank$id)
  })
})

# Drive the lead's app (site1) through a study while site2 is simulated.
# (testServer evaluates its expression among the server's own objects, so these
# argument names must not clash with them; the app has e.g. `outcome`.)
lead_session <- function(store_root, demo_id, outcome_type, body) {
  site2 <- read_site_file(file.path(app_dir, "data", sub("1", "2", DEMO_PATHS[[demo_id]])))
  store2 <- folder_store(store_root)
  shiny::testServer(shiny::shinyAppDir(app_dir), {
    clicks <- 0
    sync <- function(step) {
      for (i in 1:12) {
        respond_as(store2, "demo", "site2", site2)
        clicks <<- clicks + 1
        session$setInputs(refresh = clicks)
        if (!is.null(results[[step]])) return(invisible())
      }
      stop("step ", step, " did not finish")
    }
    session$setInputs(store_type = "folder", folder = store_root, study = "demo", me = "site1", auto = TRUE,
                      repo = "", branch = "main", gh_dir = "fedscore", token = "", demo = demo_id)
    session$setInputs(connect = 1)
    session$setInputs(sites = "site1, site2", weighting = "size", min_count = 10, outcome = outcome_type, create = 1)
    sync("summary")
    body(environment(), sync, site2)
  })
}
DEMO_PATHS <- c(surv_site1 = "surv_site1.rds", missing_site1 = "missing_site1.rds")

test_that("the lead site's app runs a survival study with ODACH", {
  lead_session(tempfile("store"), "surv_site1", "survival", function(env, sync, site2) with(env, {
    expect_equal(outcome(), "survival")
    expect_equal(config()$outcome, "survival")
    session$setInputs(rank_method = "rf", ntree = 10, start_rank = 1)
    sync("rank")
    session$setInputs(quantiles = "0.05, 0.2, 0.8, 0.95", start_cut = 1)
    sync("cut")
    session$setInputs(method = "ODACH", refine = 1, max_vars = 3, start_pars = 1)
    sync("pars")
    expect_true(all(as_df(results$pars$summary)$AUC > 0.6))
    render_ok(output, c("pars_plot", "pars_table"))
    session$setInputs(var_mode = "top", n_vars = 3, max_score = 100, times = "", start_final = 1)
    sync("final")
    f <- results$final
    expect_equal(f$method, "ODACH")
    expect_length(unlist(f$times), 3)                      # chosen from the lead's event times
    expect_true(all(as_df(f$perf)$C_index > 0.65))
    render_ok(output, c("score_tbl", "perf_tbl", "perf_plot", "coef_tbl"))
    expect_true(any(grepl("^AUC\\(", names(perf_table(f)))))
    session$setInputs(run_tl = 1)
    expect_equal(nrow(rv$tl$perf), 2)
    render_ok(output, "tl_perf")
  }))
})

test_that("the lead site's app imputes missing values with FedIMPUTE first", {
  lead_session(tempfile("store"), "missing_site1", "binary", function(env, sync, site2) with(env, {
    render_ok(output, c("missing_tbl", "impute_vars"))
    expect_setequal(candidates()$Variable, c("Lab_H", "Lab_B", "Vital_E"))
    session$setInputs(impute_sel = c("Lab_H", "Lab_B", "Vital_E"), impute_maxit = 5, impute_steps = 1,
                      impute_calibrate = FALSE, start_impute = 1)
    sync("impute")
    expect_length(results$impute$models, 3)
    expect_false(anyNA(my_data()$train[, c("Lab_H", "Lab_B")]))
    render_ok(output, c("impute_tbl", "impute_mine", "status_impute"))
    session$setInputs(rank_method = "rf", ntree = 20, start_rank = 1)
    sync("rank")
    session$setInputs(quantiles = "0.05, 0.2, 0.8, 0.95", start_cut = 1)
    sync("cut")
    expect_false(results$cut$cutoffs$Lab_H$has_na)
    session$setInputs(method = "dCLR", refine = 1, var_mode = "top", n_vars = 4, max_score = 100, start_final = 1)
    sync("final")
    expect_true(all(as_df(results$final$perf)$AUC > 0.75))
  }))
})
