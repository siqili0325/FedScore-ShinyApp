# The exchange protocol: stores, and complete studies run through a store.

sim_site <- function(n, intercept, seed, beta = c(1, -0.8, 0.5)) {
  set.seed(seed)
  make <- function(n) {
    X <- matrix(rbinom(n * 3, 1, 0.4), n, 3)
    y <- rbinom(n, 1, plogis(intercept + X %*% beta))
    data.frame(x1 = factor(X[, 1]), x2 = factor(X[, 2]), x3 = factor(X[, 3]), age = round(runif(n, 20, 90)), label = y)
  }
  list(train = make(n), validate = make(n / 2), test = make(n / 2))
}

test_that("a folder store round-trips JSON and lists, overwrites and deletes files", {
  store <- folder_store(tempfile("store"))
  obj <- list(m = matrix(c(1, 2, 3, 4) / 3, 2), v = c(0.1, NA), s = "x", l = c(TRUE, FALSE),
              t = data.frame(a = c("p", "q"), b = 1:2), one = matrix(0.5, 1, 1))
  store$write("st/a/round1.json", obj)
  store$write("st/a/round1_s1.json", list(x = 1))
  back <- store$read("st/a/round1.json")
  expect_equal(back$m, obj$m, tolerance = 1e-14)
  expect_equal(back$v, obj$v)
  expect_equal(as_df(back$t), obj$t)
  expect_equal(dim(back$one), c(1L, 1L))
  expect_setequal(store$list("st")$path, c("st/a/round1.json", "st/a/round1_s1.json"))
  expect_equal(nrow(store$list("other")), 0)
  v1 <- store$list("st")$version
  Sys.sleep(1.1)
  store$write("st/a/round1_s1.json", list(x = 2))
  expect_false(identical(store$list("st")$version, v1))
  expect_equal(store_snapshot(store, "st")[["st/a/round1_s1.json"]]$x, 2)
  store$delete("st/a/round1_s1.json")
  expect_equal(store$list("st")$path, "st/a/round1.json")
  expect_error(store$write("../escape.json", list()), "Invalid store path")
})

# A minimal in-memory stand-in for the GitHub REST API.
fake_github <- function() {
  files <- new.env()
  n <- 0
  json <- function(x, status = 200) {
    httr2::response(status, headers = list(`Content-Type` = "application/json"),
                    body = charToRaw(as.character(jsonlite::toJSON(x, auto_unbox = TRUE))))
  }
  handler <- function(req) {
    url <- sub("^https://api.github.com/repos/o/r", "", req$url)
    if (startsWith(url, "/git/trees/")) {
      return(json(list(truncated = FALSE, tree = lapply(ls(files), function(p) list(path = p, type = "blob", sha = files[[p]]$sha)))))
    }
    if (startsWith(url, "/git/blobs/")) {
      sha <- sub("/git/blobs/", "", url)
      hit <- Filter(function(p) files[[p]]$sha == sha, ls(files))
      return(httr2::response(200, body = charToRaw(files[[hit]]$content)))
    }
    if (startsWith(url, "/contents/")) {
      path <- utils::URLdecode(sub("/contents/", "", url))
      body <- req$body$data
      old <- if (exists(path, envir = files, inherits = FALSE)) files[[path]]$sha
      if (!identical(old, body$sha)) return(json(list(message = "sha does not match"), 409))
      if (req$method == "DELETE") {
        rm(list = path, envir = files)
        return(json(list(commit = list(sha = "c"))))
      }
      n <<- n + 1
      sha <- sprintf("sha%03d", n)
      assign(path, list(sha = sha, content = rawToChar(jsonlite::base64_dec(body$content))), envir = files)
      return(json(list(content = list(sha = sha)), if (is.null(old)) 201 else 200))
    }
    json(list(message = "Not Found"), 404)
  }
  list(handler = handler, files = files)
}

test_that("the GitHub store lists, reads, writes, retries on conflicts and deletes", {
  gh <- fake_github()
  httr2::local_mocked_responses(gh$handler)
  store <- github_store("o/r", "main", "fedscore", token = "t")
  expect_equal(nrow(store$list("st")), 0)
  store$write("st/study.json", list(lead = "s1", h = diag(2)))
  expect_true(exists("fedscore/st/study.json", envir = gh$files))
  files <- store$list("st")
  expect_equal(files$path, "st/study.json")
  expect_equal(store$read("st/study.json", files$version)$h, diag(2))
  # Another writer changes the file: the first upload conflicts, then succeeds.
  assign("fedscore/st/study.json", list(sha = "other", content = "{}"), envir = gh$files)
  store$write("st/study.json", list(lead = "s2"))
  expect_equal(store_snapshot(store, "st")[["st/study.json"]]$lead, "s2")
  store$delete("st/study.json")
  expect_equal(length(ls(gh$files)), 0)
  expect_error(github_store("not a repo"), "owner/name")
})

for (method in c("ODAL2", "dCLR")) {
  test_that(sprintf("a full %s study through the store matches the pooled fit", method), {
    sites <- list(A = sim_site(3000, -2, 4), B = sim_site(1500, -0.5, 5), C = sim_site(1000, -1, 6))
    store <- folder_store(tempfile("store"))
    create_study(store, "st", "A", names(sites), "size")
    st <- finish_step(store, "st", sites, "A")
    expect_equal(st$steps$summary$result$predictors, c("x1", "x2", "x3", "age"))
    w <- study_weights(st)
    expect_equal(unname(w), c(3000, 1500, 1000) / 5500)

    st <- run_step(store, "st", sites, "A", "rank", function(st) rank_request(st, "rf", 20))
    ranked <- as_df(st$steps$rank$result$table)$Variable
    expect_setequal(ranked, c("x1", "x2", "x3", "age"))
    st <- run_step(store, "st", sites, "A", "cut", cut_request)
    cutoffs <- st$steps$cut$result$cutoffs
    expect_equal(cutoffs$x1$levels, c("0", "1"))

    vars <- c("x1", "x2", "x3")
    passes <- if (method == "dCLR") 3 else 1
    st <- run_step(store, "st", sites, "A", "pars", function(st)
      fit_request(st, list("x1", c("x1", "x2"), vars), method, passes, lead_data = sites$A))
    pars <- st$steps$pars$result
    expect_equal(nrow(as_df(pars$auc)), 9)
    expect_equal(as_df(pars$summary)$Model, 1:3)

    st <- run_step(store, "st", sites, "A", "final", function(st)
      fit_request(st, list(vars), method, passes, "test", 100, sites$A))
    fin <- st$steps$final$result
    expect_true(fin$converged)
    # Rounds: ODAL2 = derivs + evaluate; dCLR = local + passes * (derivs + surrogate) + evaluate.
    expect_length(st$steps$final$rounds, if (method == "ODAL2") 2 else 2 + 2 * passes)
    pooled <- pooled_fit(sites, cutoffs, vars, method, w)
    expect_equal(fin$terms, names(pooled))
    expect_lt(max(abs(fin$beta - pooled)), if (method == "ODAL2") 0.02 else 1e-3)
    expect_lt(max(abs(fin$beta[fin$terms != "(Intercept)"] - c(1, -0.8, 0.5))), 0.3)
    expect_true(all(fin$se > 0))
    perf <- as_df(fin$perf)
    expect_equal(perf$Site, names(sites))
    expect_true(all(perf$AUC > 0.6))
    tbl <- as_df(fin$table)
    expect_equal(as.vector(tapply(tbl$Points, tbl$Variable, min)), c(0, 0, 0))

    lg <- exchange_log(store_snapshot(store, "st"), "st")
    expect_true(all(c("A", "B", "C") %in% lg$From))
    expect_equal(sum(lg$File == "final/result.json"), 1)
  })
}

test_that("errors, stale replies and restarts are handled", {
  sites <- list(A = sim_site(600, -1, 1), B = sim_site(400, -1, 2))
  store <- folder_store(tempfile("store"))
  create_study(store, "st", "A", names(sites))
  # Site B's data has one outcome class only: it replies with an error.
  bad <- sites$B
  bad$test$label <- 0L
  respond_as(store, "st", "A", sites$A)
  respond_as(store, "st", "B", bad)
  lt <- lead_task(study_state(store, "st"))
  expect_false(lt$ready)
  expect_match(lt$errors[["B"]], "both outcome classes")
  expect_true(site_task(study_state(store, "st"), "B")$retry)
  st <- finish_step(store, "st", sites, "A")   # B answers again with good data
  expect_false(is.null(st$steps$summary$result))

  st <- run_step(store, "st", sites, "A", "rank", function(st) rank_request(st, "rf", 10))
  st <- run_step(store, "st", sites, "A", "cut", cut_request)
  expect_false(is.null(st$steps$cut$result))
  # Re-running the ranking removes the later cut-off files.
  start_step(store, "st", names(store_snapshot(store, "st")), "rank", rank_request(st, "rf", 10), "A")
  st <- study_state(store, "st")
  expect_null(st$steps$cut$result)
  expect_false(any(grepl("^st/cut/", store$list("st")$path)))
  # A reply to an old request is ignored.
  store$write(round_file("st", "rank", 1, "B"), list(request = "rank-old", step_id = "rank-old", site = "B", ranking = "x1"))
  expect_null(study_state(store, "st")$steps$rank$rounds[[1]]$replies$B)

  expect_error(create_study(store, "bad name", "A", c("A", "B")), "Study names")
  expect_error(create_study(store, "st2", "A", c("A", "A")), "unique")
})

test_that("sites refuse to share statistics about very small groups", {
  sites <- list(A = sim_site(600, -1, 11), B = sim_site(400, -1, 12))
  # Site B has a category with only 3 patients.
  for (s in names(sites)) for (sp in names(sites[[s]])) {
    n <- nrow(sites[[s]][[sp]])
    sites[[s]][[sp]]$grp <- factor(c(rep("rare", if (s == "B") 3 else 60), rep("common", n - if (s == "B") 3 else 60)))
  }
  store <- folder_store(tempfile("store"))
  create_study(store, "st", "A", names(sites), min_count = 10)
  finish_step(store, "st", sites, "A")
  run_step(store, "st", sites, "A", "rank", function(st) rank_request(st, "rf", 10))
  st <- run_step(store, "st", sites, "A", "cut", cut_request)
  start_step(store, "st", names(store_snapshot(store, "st")), "final",
             fit_request(st, list(c("x1", "grp")), "dCLR", split = "test"), "A")
  respond_as(store, "st", "A", sites$A)
  respond_as(store, "st", "B", sites$B)
  lt <- lead_task(study_state(store, "st"))
  expect_false(lt$ready)
  expect_match(lt$errors[["B"]], "grp rare \\(3\\)")
  # Without the rare variable the model can be shared.
  start_step(store, "st", names(store_snapshot(store, "st")), "final",
             fit_request(st, list(c("x1", "x2")), "dCLR", split = "test"), "A")
  expect_false(is.null(finish_step(store, "st", sites, "A")$steps$final$result))

  # Evaluation needs enough events and non-events.
  few <- sites$A
  few$test$label <- c(1L, 1L, rep(0L, nrow(few$test) - 2))
  expect_error(check_min_count(few, list(kind = "evaluate", split = "test", min_count = 10), "x1"), "2 events")
  expect_silent(check_min_count(few, list(kind = "evaluate", split = "test", min_count = 0), "x1"))
  expect_error(create_study(store, "st3", "A", c("A", "B"), min_count = -1), "at least 0")
})
