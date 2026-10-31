# Exchange stores --------------------------------------------------------------
#
# Sites never see each other's data. They exchange small JSON files through a
# store that every site can read and write:
#   * folder_store(): a folder on a shared drive (network share, Dropbox,
#     OneDrive, Box, ...) or, for testing, a local folder;
#   * github_store(): a folder in a GitHub repository, e.g. a public repo so
#     anyone can audit what was shared.
# A store is a list of functions working on "/"-separated paths:
#   list(prefix)        data.frame(path, version) of the JSON files under prefix
#   read(path, version) parsed JSON
#   write(path, obj)    serialise obj to JSON and save it
#   delete(paths)       remove files
#   label               where the files live, for display

to_json <- function(obj) {
  jsonlite::toJSON(obj, auto_unbox = TRUE, digits = NA, null = "null", na = "null",
                   dataframe = "columns", rownames = FALSE)
}

from_json <- function(txt) {
  jsonlite::fromJSON(txt, simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = TRUE)
}

valid_name <- function(x) is.character(x) && length(x) == 1 && grepl("^[A-Za-z0-9_-]{1,64}$", x)

check_path <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  if (!length(parts) || any(!grepl("^[A-Za-z0-9_.-]+$", parts)) || any(parts %in% c(".", "..")))
    stop("Invalid store path: ", path)
  path
}

no_files <- function() data.frame(path = character(0), version = character(0), stringsAsFactors = FALSE)

folder_store <- function(root) {
  root <- normalizePath(path.expand(root), mustWork = FALSE)
  if (!dir.exists(root) && !dir.create(root, recursive = TRUE))
    stop("Cannot create the folder ", root)
  full <- function(path) file.path(root, check_path(path))
  list(
    label = root,
    list = function(prefix) {
      dir <- full(prefix)
      files <- list.files(dir, pattern = "\\.json$", recursive = TRUE)
      if (!length(files)) return(no_files())
      info <- file.info(file.path(dir, files))
      data.frame(path = paste(prefix, files, sep = "/"),
                 version = paste(as.numeric(info$mtime), info$size), stringsAsFactors = FALSE)
    },
    read = function(path, version = NULL) {
      from_json(paste(readLines(full(path), warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
    },
    write = function(path, obj) {
      f <- full(path)
      dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
      # Write then rename, so other sites never read a half-written file.
      tmp <- tempfile(".fedscore", tmpdir = dirname(f))
      writeLines(to_json(obj), tmp, useBytes = TRUE)
      if (!file.rename(tmp, f)) {
        file.copy(tmp, f, overwrite = TRUE)
        unlink(tmp)
      }
      invisible(path)
    },
    delete = function(paths) {
      for (p in paths) unlink(full(p))
      invisible(paths)
    }
  )
}

# GitHub: files live under `dir` on `branch` of `repo` ("owner/name"). Writing
# needs a token with read/write access to the repository contents (a
# fine-grained personal access token is enough). Reading a public repository
# works without a token, but GitHub then allows only 60 requests per hour.
github_store <- function(repo, branch = "main", dir = "fedscore", token = Sys.getenv("GITHUB_PAT"),
                         api = "https://api.github.com") {
  if (!requireNamespace("httr2", quietly = TRUE)) stop("The GitHub store needs the httr2 package.")
  if (!grepl("^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$", repo)) stop("Repository must look like owner/name.")
  dir <- gsub("^/+|/+$", "", dir)
  if (nzchar(dir)) check_path(dir)
  full <- function(path) paste0(if (nzchar(dir)) paste0(dir, "/"), check_path(path))
  shas <- new.env()                    # last known blob sha of each file

  call <- function(method, url, body = NULL, raw = FALSE) {
    req <- httr2::request(paste0(api, "/repos/", repo, url))
    req <- httr2::req_method(req, method)
    req <- httr2::req_headers(req, Accept = if (raw) "application/vnd.github.raw+json" else "application/vnd.github+json",
                              `X-GitHub-Api-Version` = "2022-11-28")
    if (nzchar(token)) req <- httr2::req_auth_bearer_token(req, token)
    req <- httr2::req_user_agent(req, "FedScore Shiny app")
    req <- httr2::req_error(req, is_error = function(resp) FALSE)
    if (!is.null(body)) req <- httr2::req_body_json(req, body, auto_unbox = TRUE)
    httr2::req_perform(req)
  }
  fail <- function(resp, what) {
    status <- httr2::resp_status(resp)
    detail <- tryCatch(httr2::resp_body_json(resp)$message, error = function(e) "")
    hint <- switch(as.character(status),
                   "401" = " Check the access token.",
                   "403" = " The token may lack write access, or the rate limit was reached.",
                   "404" = " Check the repository and branch names, and that the token can see the repository.",
                   "")
    stop(sprintf("GitHub %s failed (HTTP %d: %s).%s", what, status, detail, hint), call. = FALSE)
  }
  enc <- function(path) paste(vapply(strsplit(path, "/")[[1]], utils::URLencode, "", reserved = TRUE), collapse = "/")
  path_url <- function(path) paste0("/contents/", enc(full(path)))

  # One request lists the whole branch; remember each file's blob sha (needed
  # to update or delete it) and return the files under `prefix`.
  base <- if (nzchar(dir)) paste0(dir, "/") else ""
  refresh <- function() {
    resp <- call("GET", paste0("/git/trees/", utils::URLencode(branch, reserved = TRUE), "?recursive=1"))
    if (httr2::resp_status(resp) == 409) tree <- list()              # empty repository
    else if (httr2::resp_status(resp) != 200) fail(resp, "listing")
    else tree <- httr2::resp_body_json(resp)$tree
    keep <- Filter(function(e) e$type == "blob" && startsWith(e$path, base) && endsWith(e$path, ".json"), tree)
    rm(list = ls(shas, all.names = TRUE), envir = shas)
    for (e in keep) assign(substring(e$path, nchar(base) + 1), e$sha, envir = shas)
  }
  list_files <- function(prefix) {
    refresh()
    paths <- Filter(function(p) startsWith(p, paste0(prefix, "/")), ls(shas, all.names = TRUE))
    if (!length(paths)) return(no_files())
    data.frame(path = paths, version = unname(unlist(mget(paths, envir = shas))), stringsAsFactors = FALSE)
  }

  list(
    label = sprintf("github.com/%s (branch %s, folder %s)", repo, branch, if (nzchar(dir)) dir else "/"),
    list = list_files,
    read = function(path, version = NULL) {
      if (is.null(version)) version <- shas[[path]]
      if (is.null(version)) stop("Unknown file ", path)
      resp <- call("GET", paste0("/git/blobs/", version), raw = TRUE)
      if (httr2::resp_status(resp) != 200) fail(resp, "download")
      from_json(httr2::resp_body_string(resp))
    },
    write = function(path, obj) {
      content <- jsonlite::base64_enc(charToRaw(enc2utf8(as.character(to_json(obj)))))
      for (attempt in 1:5) {
        body <- list(message = paste("FedScore:", path), content = gsub("\n", "", content), branch = branch)
        if (!is.null(shas[[path]])) body$sha <- shas[[path]]
        resp <- call("PUT", path_url(path), body)
        status <- httr2::resp_status(resp)
        if (status %in% c(200, 201)) {
          assign(path, httr2::resp_body_json(resp)$content$sha, envir = shas)
          return(invisible(path))
        }
        # 409/422: another site committed at the same moment, or our sha is
        # stale. Refresh and try again.
        if (!status %in% c(409, 422) || attempt == 5) fail(resp, "upload")
        Sys.sleep(stats::runif(1, 0.5, 2))
        try(refresh(), silent = TRUE)
      }
    },
    delete = function(paths) {
      if (!all(vapply(paths, exists, NA, envir = shas, inherits = FALSE))) refresh()
      for (p in paths) {
        sha <- shas[[p]]
        if (is.null(sha)) next
        resp <- call("DELETE", path_url(p), list(message = paste("FedScore: remove", p), sha = sha, branch = branch))
        if (!httr2::resp_status(resp) %in% c(200, 404)) fail(resp, "delete")
        rm(list = p, envir = shas)
      }
      invisible(paths)
    }
  )
}

# Read every JSON file under `prefix`, re-downloading only files whose version
# changed since the last call (`cache` is an environment kept between calls).
# Returns a named list: path -> parsed JSON.
store_snapshot <- function(store, prefix, cache = new.env()) {
  files <- store$list(prefix)
  out <- list()
  for (i in seq_len(nrow(files))) {
    key <- files$path[i]
    hit <- cache[[key]]
    if (is.null(hit) || !identical(hit$version, files$version[i])) {
      obj <- tryCatch(store$read(key, files$version[i]), error = function(e) NULL)
      if (is.null(obj)) next
      hit <- list(version = files$version[i], obj = obj)
      assign(key, hit, envir = cache)
    }
    out[[key]] <- hit$obj
  }
  out
}
