# FedScore Shiny app: build a federated scoring system for a binary or a
# survival outcome, optionally imputing missing values first (FedIMPUTE).
#
# Every site runs this app on its own computer with its own data. Sites only
# exchange summary statistics, as JSON files in a shared store: a GitHub
# repository or a folder every site can reach. One site, the lead, creates the
# study and starts each step; the app answers the lead's requests from local
# data and shows what was shared.
#
# Run with shiny::runApp() from this folder. Files in R/ are sourced automatically.

library(shiny)
library(bslib)
library(ggplot2)

options(shiny.maxRequestSize = 200 * 1024^2)

DEMO_FILES <- c(site1 = "data/site1.rds", site2 = "data/site2.rds",
                surv_site1 = "data/surv_site1.rds", surv_site2 = "data/surv_site2.rds",
                missing_site1 = "data/missing_site1.rds", missing_site2 = "data/missing_site2.rds")
DEMO_CHOICES <- c("(none)" = "", "Binary outcome, site 1" = "site1", "Binary outcome, site 2" = "site2",
                  "Survival outcome, site 1" = "surv_site1", "Survival outcome, site 2" = "surv_site2",
                  "Binary, with missing values, site 1" = "missing_site1",
                  "Binary, with missing values, site 2" = "missing_site2")

ALGORITHMS <- list(
  binary = c("ODAL2: homogeneous sites" = "ODAL2",
             "dCLR: heterogeneous sites (site-specific baseline risk)" = "dCLR"),
  survival = c("ODACH: Cox model, site-specific baseline hazards" = "ODACH",
               "DiSAH: additive hazards, site-specific baseline hazards (one round)" = "DiSAH"))
RANK_METHODS <- list(
  binary = c("Random forest importance" = "rf", "Univariable AUC (uses validation set)" = "auc"),
  survival = c("Random survival forest importance" = "rf"))

# Connection defaults; set these environment variables (e.g. in .Renviron) to
# prefill the Study tab.
DEFAULTS <- list(
  store  = Sys.getenv("FEDSCORE_STORE", "folder"),
  folder = Sys.getenv("FEDSCORE_FOLDER", file.path("~", "FedScore-exchange")),
  repo   = Sys.getenv("FEDSCORE_REPO"),
  branch = Sys.getenv("FEDSCORE_BRANCH", "main"),
  dir    = Sys.getenv("FEDSCORE_DIR", "fedscore"),
  study  = Sys.getenv("FEDSCORE_STUDY", "my-study"),
  site   = Sys.getenv("FEDSCORE_SITE")
)

# Controls on the left, results on the right (stacks on narrow screens).
# A plain card is used instead of bslib::sidebar() so the app also works with
# older bslib releases.
step_layout <- function(controls, ...) {
  layout_columns(
    col_widths = breakpoints(sm = 12, lg = c(3, 9)),
    card(class = "step-controls", card_body(controls)),
    div(...)
  )
}

# Settings that only the lead site can change. They are shown or hidden via
# the hidden `lead_view` checkbox, which the server keeps in sync with the
# site's role (the server also checks the role before acting). A condition on
# an output (`output.x`) would be simpler but stops the Shiny client from
# starting with some shiny/bslib versions. The hidden `outcome_view` input
# does the same for settings that depend on the study's outcome.
lead_only <- function(...) {
  tagList(
    conditionalPanel("input.lead_view", ...),
    conditionalPanel("!input.lead_view", p(class = "step-hint", "The lead site chooses the settings and starts this step."))
  )
}

status_card <- function(step) card(card_header("Progress"), uiOutput(paste0("status_", step)))

ui <- page_navbar(
  title = "FedScore",
  theme = bs_theme(version = 5, bootswatch = "flatly", base_font = font_collection("system-ui", "sans-serif")),
  header = tagList(div(style = "display: none", checkboxInput("lead_view", NULL, FALSE),
                       textInput("outcome_view", NULL, "binary")), tags$head(tags$style(HTML("
    .shared-note { font-size: .85rem; color: var(--bs-secondary-color); }
    .table { font-size: .9rem; }
    .step-hint { color: var(--bs-secondary-color); font-style: italic; }
    .json-view { max-height: 480px; overflow: auto; font-size: .8rem; }
  ")))),

  nav_panel(
    "1 · Study",
    step_layout(
      tagList(
        h6("Exchange store"),
        radioButtons("store_type", NULL, c("Shared folder" = "folder", "GitHub repository" = "github"),
                     selected = DEFAULTS$store),
        conditionalPanel("input.store_type == 'folder'",
          textInput("folder", "Folder", DEFAULTS$folder),
          helpText("A folder every site can open: a network drive, or a Dropbox, OneDrive or Box folder.")),
        conditionalPanel("input.store_type == 'github'",
          textInput("repo", "Repository (owner/name)", DEFAULTS$repo),
          textInput("branch", "Branch", DEFAULTS$branch),
          textInput("gh_dir", "Folder in the repository", DEFAULTS$dir),
          passwordInput("token", "Access token"),
          helpText("Needs read and write access to the repository contents. Leave blank to use the",
                   "GITHUB_PAT environment variable. The token is used only on this computer.")),
        textInput("study", "Study name", DEFAULTS$study),
        textInput("me", "Your site name", DEFAULTS$site),
        actionButton("connect", "Connect", class = "btn-primary w-100"),
        hr(),
        h6("Your data"),
        fileInput("upload", NULL, accept = ".rds", placeholder = "Site data (.rds)"),
        selectInput("demo", "Or try the demo data", DEMO_CHOICES, selectize = FALSE),
        uiOutput("data_status"),
        hr(),
        checkboxInput("auto", "Reply and continue automatically", TRUE),
        helpText("When ticked, this app answers the lead's requests as soon as they appear (and, at the lead",
                 "site, moves on once every site has replied). Untick to review each step and click instead.")
      ),
      card(card_header("Study"), uiOutput("study_panel")),
      card(card_header("Sites"), uiOutput("status_summary"), tableOutput("site_table"), uiOutput("site_messages")),
      card(card_header("How it works"), markdown("
FedScore builds an [AutoScore](https://github.com/nliulab/AutoScore)-style point score across sites **without pooling patient-level data**, for a binary outcome or a survival (time-to-event) outcome.

* **Every site runs this app on its own computer**, with its own data. Data never leaves the site.
* Sites exchange only summary statistics, as JSON files in a **shared store**: a GitHub repository (a public one lets anyone audit the exchange) or a shared folder.
* **One site is the lead.** It creates the study, starts each step and combines the replies. The other sites answer its requests.

To try it on one computer, connect two browser tabs to the same folder: one as `site1` with a site-1 demo, which creates the study, and one as `site2` with the matching site-2 demo.

Steps: **missing values** (optional: FedIMPUTE federates imputation models; each site fills its own gaps) → **ranking** (each site ranks variables) → **cut-offs** (each site shares quantiles) → **model selection** (federated regression for the top 1…m variables: ODAL2 or dCLR for a binary outcome, ODACH or DiSAH for survival) → **final score** (points table, evaluated at every site) → **TLScore** (optional, local only). Every file exchanged is listed under **Exchange log**.
"))
    )
  ),

  nav_panel(
    "2 · Missing values",
    step_layout(
      tagList(
        lead_only(
          uiOutput("impute_vars"),
          numericInput("impute_maxit", "MICE iterations (local completion)", 5, min = 1, max = 20, step = 1),
          numericInput("impute_steps", "Newton rounds", 1, min = 0, max = 5, step = 1),
          helpText("0 = average of the sites' own models. Each round shares gradients and Hessians once more;",
                   "one round gives the pooled-data fit for numeric variables, binary ones may need 2-3."),
          checkboxInput("impute_calibrate", "Recalibrate the models at each site (transfer learning)", FALSE),
          actionButton("start_impute", "Start imputation", class = "btn-primary w-100")
        ),
        p(class = "shared-note mt-2", "Optional. Without it, missing values form an 'Unknown' category.",
          "Shared by each site: coefficients, then gradients and Hessians, of each imputation model.",
          "Imputed values never leave the site.")
      ),
      status_card("impute"),
      card(card_header("Missing values by site"), tableOutput("missing_tbl")),
      card(card_header("Federated imputation models"), tableOutput("impute_tbl"), uiOutput("impute_mine"))
    )
  ),

  nav_panel(
    "3 · Ranking",
    step_layout(
      tagList(
        lead_only(
          radioButtons("rank_method", "Ranking method", RANK_METHODS$binary),
          conditionalPanel("input.rank_method == 'rf'",
                           numericInput("ntree", "Number of trees", 100, min = 10, max = 2000, step = 10)),
          actionButton("start_rank", "Start ranking", class = "btn-primary w-100")
        ),
        p(class = "shared-note mt-2", "Shared by each site: the rank position of every variable.")
      ),
      status_card("rank"),
      layout_columns(
        col_widths = c(6, 6),
        card(full_screen = TRUE, card_header("Federated ranking"), tableOutput("rank_table")),
        card(full_screen = TRUE, card_header("Rank by site"), plotOutput("rank_plot", height = "520px"))
      )
    )
  ),

  nav_panel(
    "4 · Cut-offs",
    step_layout(
      tagList(
        lead_only(
          textInput("quantiles", "Quantiles for categorisation", "0.05, 0.2, 0.8, 0.95"),
          actionButton("start_cut", "Start cut-offs", class = "btn-primary w-100"),
          hr(),
          h6("Fine-tune a variable"),
          selectInput("tune_var", NULL, choices = NULL, selectize = FALSE),
          textInput("tune_cuts", "Cut-offs (comma separated)"),
          layout_columns(
            actionButton("tune_apply", "Publish", class = "btn-outline-primary"),
            actionButton("tune_reset", "Reset all", class = "btn-outline-secondary")
          ),
          helpText("Publishing new cut-offs clears model selection and the final score.")
        ),
        p(class = "shared-note mt-2", "Shared by each site: these quantiles for every continuous variable, and the observed levels of categorical ones.")
      ),
      status_card("cut"),
      card(full_screen = TRUE, card_header("Unified cut-offs (in federated rank order)"), tableOutput("cut_table"))
    )
  ),

  nav_panel(
    "5 · Model selection",
    step_layout(
      tagList(
        lead_only(
          radioButtons("method", "Federated algorithm", ALGORITHMS$binary),
          conditionalPanel("input.method == 'dCLR' || input.method == 'ODACH'",
            numericInput("refine", "Surrogate passes", 1, min = 1, max = 5, step = 1),
            helpText("1 = one-shot dCLR/ODACH as published. Each extra pass costs two more rounds of exchange",
                     "and moves the estimate closer to the pooled-data fit, which helps when a category is",
                     "sparse or event-free at some site.")),
          sliderInput("max_vars", "Variables to try", min = 1, max = 20, value = 10, step = 1),
          actionButton("start_pars", "Start model selection", class = "btn-primary w-100")
        ),
        p(class = "shared-note mt-2", "For each model size the federated model is fitted, and each site shares its",
          "validation AUC (binary) or integrated AUC (survival).")
      ),
      status_card("pars"),
      card(full_screen = TRUE, card_header("Parsimony plot"), plotOutput("pars_plot", height = "420px")),
      card(card_header("Validation AUC"), tableOutput("pars_table"))
    )
  ),

  nav_panel(
    "6 · Final score",
    step_layout(
      tagList(
        lead_only(
          radioButtons("var_mode", "Variables", c("Top-ranked" = "top", "Custom selection" = "custom")),
          conditionalPanel("input.var_mode == 'top'", numericInput("n_vars", "Number of variables", 6, min = 1, step = 1)),
          conditionalPanel("input.var_mode == 'custom'", uiOutput("custom_vars")),
          numericInput("max_score", "Maximum total score", 100, min = 10, step = 10),
          conditionalPanel("input.outcome_view == 'survival'",
            textInput("times", "Times for AUC(t)", ""),
            helpText("Comma separated, in the unit of label_time. Blank = the 25th, 50th and 75th",
                     "percentiles of the lead site's event times.")),
          p(class = "small", textOutput("final_method", inline = TRUE)),
          actionButton("start_final", "Build score", class = "btn-primary w-100")
        ),
        p(class = "shared-note mt-2", "Shared by each site: test-set performance. Binary: AUC with CI, best",
          "operating point, ROC curve. Survival: integrated AUC, C-index with CI, AUC at the chosen times.")
      ),
      status_card("final"),
      navset_card_tab(
        full_screen = TRUE,
        nav_panel("Score table", tableOutput("score_tbl"), downloadButton("dl_score", "Download score table", class = "btn-sm")),
        nav_panel("Test performance", tableOutput("perf_tbl"), plotOutput("perf_plot", height = "420px"),
                  downloadButton("dl_pred", "Download your test-set scores", class = "btn-sm")),
        nav_panel("Coefficients", tableOutput("coef_tbl"))
      )
    )
  ),

  nav_panel(
    "7 · TLScore",
    step_layout(
      tagList(
        actionButton("run_tl", "Fit TLScore for your site", class = "btn-primary w-100"),
        p(class = "shared-note mt-2", "Runs only on your computer and shares nothing. A lasso correction",
          "δ is fitted on your training data with the federated linear predictor as offset, giving",
          "β_TL = β_fed + δ. Available for ODAL2, dCLR and ODACH models.")
      ),
      card(card_header("Test performance at your site"), tableOutput("tl_perf")),
      card(full_screen = TRUE, card_header("FedScore vs. TLScore points"), tableOutput("tl_tbl"),
           downloadButton("dl_tl", "Download TLScore table", class = "btn-sm"))
    )
  ),

  nav_panel(
    "Exchange log",
    card(card_header("Every file in the study"),
         p(class = "shared-note", "These files are all that sites share. Numbers = count of numeric values in the file."),
         tableOutput("log_tbl"),
         downloadButton("dl_log", "Download all files (one JSON)", class = "btn-sm")),
    card(card_header("File contents"),
         selectInput("log_file", NULL, choices = NULL, selectize = FALSE, width = "100%"),
         div(class = "json-view", verbatimTextOutput("log_json")))
  ),

  nav_spacer(),
  nav_item(textOutput("whoami", inline = TRUE)),
  nav_item(actionButton("refresh", "Refresh", class = "btn-sm btn-outline-light")),
  nav_item(tags$a("Paper", href = "https://doi.org/10.1016/j.jbi.2023.104485", target = "_blank"))
)

server <- function(input, output, session) {
  rv <- reactiveValues(store = NULL, study = NULL, me = NULL, poll = 3000, data = NULL, data_name = NULL, tl = NULL)
  cache <- new.env()
  snap <- reactiveVal(NULL)            # path -> parsed JSON for every file in the study
  results <- reactiveValues()          # result of each step (changes only when its file does)
  tick <- reactiveVal(0)
  refresh_now <- function() tick(isolate(tick()) + 1)
  handled <- character(0)              # requests this session already acted on automatically

  run_safely <- function(expr) {
    tryCatch(expr, error = function(e) {
      message("FedScore error: ", conditionMessage(e))
      showNotification(conditionMessage(e), type = "error", duration = 10)
      NULL
    })
  }
  # Prerequisite checks for button handlers: tell the user, then stop quietly.
  guard <- function(ok, message) {
    if (!isTRUE(ok)) {
      showNotification(message, type = "warning")
      req(FALSE)
    }
  }

  # Connecting to the store ------------------------------------------------------
  observeEvent(input$connect, {
    me <- trimws(input$me)
    study <- trimws(input$study)
    guard(valid_name(me), "Site name: use letters, digits, '-' or '_' (up to 64 characters).")
    guard(valid_name(study), "Study name: use letters, digits, '-' or '_' (up to 64 characters).")
    store <- run_safely({
      s <- if (input$store_type == "github") {
        token <- if (nzchar(input$token)) input$token else Sys.getenv("GITHUB_PAT")
        github_store(trimws(input$repo), trimws(input$branch), trimws(input$gh_dir), token)
      } else folder_store(trimws(input$folder))
      s$list(study)
      s
    })
    req(store)
    rm(list = ls(cache, all.names = TRUE), envir = cache)
    snap(NULL)
    rv$store <- store
    rv$study <- study
    rv$me <- me
    rv$poll <- if (input$store_type == "github") 10000 else 3000
    refresh_now()
  })

  # Poll the store. snap() only changes when a file does.
  observe({
    tick()
    input$refresh
    req(rv$store)
    invalidateLater(rv$poll)
    s <- tryCatch(store_snapshot(rv$store, rv$study, cache), error = function(e) {
      message("FedScore: could not read the store: ", conditionMessage(e))
      NULL
    })
    if (!is.null(s)) snap(s)
  })

  state <- reactive({
    req(rv$store, !is.null(snap()))
    parse_study(snap(), rv$study)
  })
  config <- reactive(if (is.null(rv$store) || is.null(snap())) NULL else state()$config)
  is_lead <- reactive(!is.null(config()) && identical(config()$lead, rv$me))
  observe(updateCheckboxInput(session, "lead_view", value = is_lead()))
  # The study's outcome; a reactiveVal only notifies when the value changes, so
  # the lead's choices below are not reset each time a file arrives.
  outcome <- reactiveVal("binary")
  observe(outcome(if (is.null(config()) || is.null(config()$outcome)) "binary" else config()$outcome))
  observe({
    o <- outcome()
    updateTextInput(session, "outcome_view", value = o)
    updateRadioButtons(session, "method", choices = ALGORITHMS[[o]], selected = unname(ALGORITHMS[[o]][1]))
    updateRadioButtons(session, "rank_method", choices = RANK_METHODS[[o]], selected = "rf")
  })
  # Your data as used from the ranking onwards (after imputation, if any).
  my_data <- reactive({
    req(rv$data)
    if (is.null(results$impute)) rv$data else apply_imputation(rv$data, results$impute)
  })

  observe({
    st <- if (is.null(config())) NULL else state()
    for (s in names(STEPS)) results[[s]] <- st$steps[[s]]$result
  })

  output$whoami <- renderText({
    if (is.null(rv$store)) return("Not connected")
    sprintf("%s%s - %s", rv$me, if (is_lead()) " (lead)" else "", rv$study)
  })

  # Your data --------------------------------------------------------------------
  load_data <- function(path, name) {
    d <- run_safely(read_site_file(path))
    if (is.null(d)) return()
    rv$data <- d
    rv$data_name <- name
    rv$tl <- NULL
  }
  observeEvent(input$upload, load_data(input$upload$datapath, input$upload$name))
  observeEvent(input$demo, {
    req(input$demo)
    load_data(DEMO_FILES[[input$demo]], sprintf("demo data (%s)", input$demo))
    if (!nzchar(input$me)) updateTextInput(session, "me", value = sub("^.*_", "", input$demo))
  })
  output$data_status <- renderUI({
    if (is.null(rv$data)) return(p(class = "step-hint", "No data loaded."))
    d <- rv$data
    o <- data_outcome(d$train)
    n_na <- sum(vapply(d, function(x) sum(is.na(x[, predictor_names(x), drop = FALSE])), numeric(1)))
    div(class = "alert alert-success py-2 small mb-0",
        sprintf("%s: %s outcome; %d train / %d validation / %d test rows, %d predictors%s.", rv$data_name,
                o, nrow(d$train), nrow(d$validate), nrow(d$test), length(predictor_names(d$train)),
                if (n_na) sprintf("; %d missing values", n_na) else ""))
  })

  # Acting on requests -----------------------------------------------------------
  send_reply <- function(st, auto) {
    t <- site_task(st, rv$me)
    if (is.null(t) || (auto && (t$retry || t$request$id %in% handled))) return(FALSE)
    if (is.null(rv$data)) {
      if (!auto) showNotification("Load your data first (Study tab).", type = "warning")
      return(FALSE)
    }
    handled <<- c(handled, t$request$id)
    reply <- withProgress(message = sprintf("%s: preparing your reply", STEPS[[t$step]]), value = 0.3,
                          tryCatch(site_reply(t$request, rv$me, step_data(st, t$step, rv$data)),
                                   error = function(e) c(list(request = t$request$id, step_id = t$request$step_id,
                                                              site = rv$me, created = now_utc()),
                                                         list(error = conditionMessage(e)))))
    ok <- run_safely({
      rv$store$write(round_file(rv$study, t$step, t$round, rv$me), reply)
      TRUE
    })
    if (!is.null(reply$error))
      showNotification(paste("Your site could not answer:", reply$error), type = "error", duration = NULL)
    isTRUE(ok)
  }

  combine_replies <- function(st, auto) {
    if (!identical(st$config$lead, rv$me)) return(FALSE)
    lt <- lead_task(st)
    key <- paste0("lead:", lt$request$id)
    if (is.null(lt) || !lt$ready || (auto && key %in% handled)) return(FALSE)
    handled <<- c(handled, key)
    ok <- run_safely(withProgress(message = sprintf("%s: combining replies", STEPS[[lt$step]]), value = 0.5, {
      out <- advance_step(st, lt$step, step_data(st, lt$step, rv$data))
      rv$store$write(out$path, out$obj)
      TRUE
    }))
    isTRUE(ok)
  }

  observe({
    req(isTRUE(input$auto), config(), rv$me %in% config()$sites)
    st <- state()
    rv$data
    isolate(if (send_reply(st, TRUE) | combine_replies(st, TRUE)) refresh_now())
  })
  observeEvent(input$do_reply, if (send_reply(state(), FALSE)) refresh_now())
  observeEvent(input$do_advance, if (combine_replies(state(), FALSE)) refresh_now())

  # Lead site: starting steps --------------------------------------------------------
  start <- function(step, make_request, needs = character(0)) {
    guard(is_lead(), "Only the lead site can start steps.")
    st <- state()
    missing <- needs[vapply(needs, function(s) is.null(st$steps[[s]]$result), logical(1))]
    guard(!length(missing), sprintf("Finish %s first.", paste(STEPS[missing], collapse = " and ")))
    earlier <- names(STEPS)[seq_len(match(step, names(STEPS)) - 1)]
    open <- open_step(st)
    guard(is.null(open) || !open %in% earlier, sprintf("Wait until %s has finished.", STEPS[[open]]))
    ok <- run_safely({
      start_step(rv$store, rv$study, names(snap()), step, make_request(st), rv$me)
      TRUE
    })
    if (isTRUE(ok)) refresh_now()
  }

  output$study_panel <- renderUI({
    if (is.null(rv$store)) return(p(class = "step-hint", "Choose the exchange store, the study and your site name, then click Connect."))
    cfg <- config()
    if (is.null(cfg)) {
      return(tagList(
        p(sprintf("There is no study called '%s' in %s yet. Create it, and your site (%s) becomes its lead site.",
                  rv$study, rv$store$label, rv$me)),
        textInput("sites", "Sites taking part (comma separated, including yours)", value = rv$me, width = "100%"),
        radioButtons("outcome", "Outcome", c("Binary" = "binary", "Survival (time to event)" = "survival"),
                     selected = if (is.null(rv$data)) "binary" else data_outcome(rv$data$train), inline = TRUE),
        radioButtons("weighting", "Site weights", c("Proportional to training size" = "size", "Equal" = "equal"), inline = TRUE),
        numericInput("min_count", "Minimum patients per shared category", 10, min = 0, step = 1),
        helpText("Sites refuse to share model statistics about a category (or a validation/test set's events)",
                 "with fewer patients than this. 0 switches the check off."),
        actionButton("create", "Create study", class = "btn-primary")
      ))
    }
    lead <- identical(cfg$lead, rv$me)
    tagList(
      if (!rv$me %in% cfg$sites)
        div(class = "alert alert-danger py-2", sprintf("Your site name '%s' is not in this study. Its sites are: %s.",
                                                      rv$me, paste(cfg$sites, collapse = ", "))),
      tags$dl(class = "row mb-2",
        tags$dt(class = "col-sm-3", "Store"), tags$dd(class = "col-sm-9", rv$store$label),
        tags$dt(class = "col-sm-3", "Lead site"), tags$dd(class = "col-sm-9", cfg$lead),
        tags$dt(class = "col-sm-3", "Sites"), tags$dd(class = "col-sm-9", paste(cfg$sites, collapse = ", ")),
        tags$dt(class = "col-sm-3", "Your role"), tags$dd(class = "col-sm-9", if (lead) "lead site" else "site"),
        tags$dt(class = "col-sm-3", "Outcome"), tags$dd(class = "col-sm-9", outcome()),
        tags$dt(class = "col-sm-3", "Minimum count"),
        tags$dd(class = "col-sm-9", if (isTRUE(cfg$min_count > 0)) sprintf("%d patients per shared category", as.integer(cfg$min_count)) else "off"),
        tags$dt(class = "col-sm-3", "Progress"), tags$dd(class = "col-sm-9", progress_text(state()))),
      if (lead) tagList(
        radioButtons("weighting", "Site weights", c("Proportional to training size" = "size", "Equal" = "equal"),
                     selected = state()$steps$summary$rounds[[1]]$request$weighting, inline = TRUE),
        actionButton("restart_summary", "Re-run data summary", class = "btn-sm btn-outline-primary"),
        actionButton("delete_study", "Delete study", class = "btn-sm btn-outline-danger")
      )
    )
  })

  progress_text <- function(st) {
    done <- names(STEPS)[vapply(names(STEPS), function(s) !is.null(st$steps[[s]]$result), logical(1))]
    open <- open_step(st)
    paste0(if (length(done)) paste("finished:", paste(STEPS[done], collapse = ", ")) else "nothing finished yet",
           if (!is.null(open)) sprintf("; in progress: %s", STEPS[[open]]) else "")
  }

  observeEvent(input$create, {
    sites <- unique(trimws(strsplit(input$sites, ",")[[1]]))
    sites <- sites[nzchar(sites)]
    guard(rv$me %in% sites, "Include your own site in the list.")
    ok <- run_safely({
      create_study(rv$store, rv$study, rv$me, sites, input$weighting, input$min_count, input$outcome)
      TRUE
    })
    if (isTRUE(ok)) refresh_now()
  })
  observeEvent(input$restart_summary, start("summary", function(st) summary_request(outcome(), input$weighting)))
  observeEvent(input$delete_study, {
    showModal(modalDialog(
      title = "Delete study?",
      sprintf("This removes all %d files of '%s' from the store, for every site.", length(snap()), rv$study),
      footer = tagList(modalButton("Cancel"), actionButton("delete_confirm", "Delete", class = "btn-danger"))))
  })
  observeEvent(input$delete_confirm, {
    removeModal()
    guard(is_lead(), "Only the lead site can delete the study.")
    paths <- names(snap())
    run_safely(rv$store$delete(c(setdiff(paths, study_file(rv$study, "study.json")), study_file(rv$study, "study.json"))))
    refresh_now()
  })

  # Progress cards ---------------------------------------------------------------
  settings_text <- function(req) {
    switch(req$kind,
      summary = sprintf("%s outcome; site weights: %s.", if (identical(req$outcome, "survival")) "Survival" else "Binary",
                        if (identical(req$weighting, "equal")) "equal" else "proportional to training size"),
      impute_local = sprintf("Imputing %s; %d MICE iterations, %d Newton round(s)%s.",
                             paste(vapply(req$models, `[[`, "", "var"), collapse = ", "), as.integer(req$maxit),
                             as.integer(req$steps), if (isTRUE(req$calibrate)) ", recalibrated at each site" else ""),
      rank = if (req$method == "rf") sprintf("Random %sforest importance, %d trees.", if (outcome() == "survival") "survival " else "",
                                             as.integer(req$ntree)) else "Univariable AUC on the validation sets.",
      cut = sprintf("Quantiles %s.", paste(req$probs, collapse = ", ")),
      {
        alg <- if (req$method %in% c("dCLR", "ODACH")) sprintf("%s with %d surrogate pass(es)", req$method, as.integer(req$refinements)) else req$method
        if (isTRUE(req$full)) sprintf("%s, %d variables: %s.", alg, length(req$models[[1]]$vars), paste(req$models[[1]]$vars, collapse = ", "))
        else sprintf("%s, models with the top 1 to %d variables.", alg, length(req$models))
      })
  }

  status_ui <- function(step) {
    if (is.null(rv$store)) return(p(class = "step-hint", "Connect to a study on the Study tab."))
    st <- state()
    if (is.null(st)) return(p(class = "step-hint", "The study has not been created yet."))
    s <- st$steps[[step]]
    if (!length(s$rounds))
      return(p(class = "step-hint", if (is_lead()) "Not started. Choose the settings and start this step." else
        "Not started. Waiting for the lead site."))
    rounds <- lapply(seq_along(s$rounds), function(k) {
      r <- s$rounds[[k]]
      badges <- lapply(st$config$sites, function(x) {
        rep <- r$replies[[x]]
        if (is.null(rep)) span(class = "badge text-bg-light border me-1", paste(x, "waiting"))
        else if (!is.null(rep$error)) span(class = "badge text-bg-danger me-1", title = rep$error, paste(x, "error"))
        else span(class = "badge text-bg-success me-1", paste(x, "replied"))
      })
      errors <- Filter(Negate(is.null), lapply(r$replies, `[[`, "error"))
      div(class = "mb-1",
          strong(sprintf("Round %d. ", k)), sprintf("Sites share the %s. ", shared_content(r$request)), badges,
          lapply(names(errors), function(x) div(class = "text-danger small", sprintf("%s: %s", x, errors[[x]]))))
    })
    footer <- if (!is.null(s$result)) {
      div(class = "text-success mt-2", sprintf("Finished (%s).", s$result$created))
    } else {
      mine <- if (rv$me %in% st$config$sites) site_task(st, rv$me)
      lt <- lead_task(st)
      tagList(
        if (!is.null(mine) && is.null(rv$data))
          div(class = "alert alert-warning py-2 mt-2 mb-0", "Load your data on the Study tab to reply."),
        if (!is.null(mine) && !is.null(rv$data) && (mine$retry || !isTRUE(input$auto)))
          actionButton("do_reply", if (mine$retry) "Send my reply again" else "Send my reply", class = "btn-sm btn-primary mt-2"),
        if (is_lead() && isTRUE(lt$ready))
          (if (isTRUE(input$auto)) div(class = "step-hint mt-2", "All sites have replied; combining.") else
            actionButton("do_advance", "Combine replies and continue", class = "btn-sm btn-primary mt-2")),
        if (is.null(mine) && !isTRUE(lt$ready)) div(class = "step-hint mt-2", "Waiting for the other sites.")
      )
    }
    tagList(p(class = "small text-body-secondary mb-2", settings_text(s$rounds[[1]]$request)), rounds, footer)
  }
  for (step in names(STEPS)) local({
    st <- step
    output[[paste0("status_", st)]] <- renderUI(status_ui(st))
  })

  # Step 1: data summary -----------------------------------------------------------
  output$site_table <- renderTable({
    req(results$summary)
    as_df(results$summary$table)
  })
  output$site_messages <- renderUI({
    r <- results$summary
    req(r)
    tagList(
      lapply(unlist(r$warnings), function(m) div(class = "alert alert-warning py-2 mb-2", m)),
      div(class = "alert alert-success py-2 mb-0",
          sprintf("%d sites ready, %d shared predictors.", nrow(as_df(r$table)), length(r$predictors)))
    )
  })

  # Step 2: missing values (FedIMPUTE) ------------------------------------------------
  output$missing_tbl <- renderTable({
    validate(need(results$summary, "Waiting for the data summary."))
    miss <- as_df(results$summary$missing)
    validate(need(nrow(miss) > 0, "No missing values at any site. You can skip this step."))
    miss
  }, digits = 1)

  candidates <- reactive({
    req(results$summary)
    impute_candidates(results$summary)
  })
  output$impute_vars <- renderUI({
    validate(need(results$summary, "Finish the data summary first."))
    cand <- candidates()
    if (!nrow(cand)) return(p(class = "step-hint", "Nothing to impute."))
    checkboxGroupInput("impute_sel", "Variables to impute", selected = cand$Variable,
                       choiceNames = sprintf("%s (%s)", cand$Variable, ifelse(cand$Family == "gaussian", "numeric", "binary")),
                       choiceValues = cand$Variable)
  })
  observeEvent(input$start_impute, {
    guard(length(input$impute_sel) > 0, "Choose at least one variable to impute.")
    start("impute", function(st) impute_request(st, input$impute_sel, input$impute_maxit, input$impute_steps,
                                                input$impute_calibrate), "summary")
  })
  output$impute_tbl <- renderTable({
    validate(need(results$impute, "No imputation models yet (optional step)."))
    do.call(rbind, lapply(results$impute$models, function(m) {
      n <- unlist(m$n)
      cbind(data.frame(Variable = m$var, Model = if (m$family == "gaussian") "linear" else "logistic",
                       Coefficients = length(m$coef)),
            as.data.frame(as.list(setNames(as.integer(n), paste("Observed at", names(n)))), check.names = FALSE))
    }))
  })
  output$impute_mine <- renderUI({
    req(results$impute, rv$data)
    counts <- attr(my_data(), "imputed")
    total <- sum(unlist(counts))
    p(class = "small mt-2", sprintf("Your data: %d values imputed (%s).", total,
                                    paste(sprintf("%s %d", names(counts), vapply(counts, function(x) sum(unlist(x)), numeric(1))),
                                          collapse = ", ")))
  })

  # Step 3: ranking ----------------------------------------------------------------
  observeEvent(input$start_rank, start("rank", function(st) rank_request(st, input$rank_method, input$ntree), "summary"))

  ranking <- reactive({
    req(results$rank)
    as_df(results$rank$table)
  })
  observe({
    p <- nrow(ranking())
    updateSliderInput(session, "max_vars", max = p, value = min(10, p))
    updateNumericInput(session, "n_vars", max = p, value = min(6, p))
  })

  output$rank_table <- renderTable({
    validate(need(results$rank, "No federated ranking yet."))
    ranking()
  }, digits = 2)

  output$rank_plot <- renderPlot({
    d <- ranking()
    long <- do.call(rbind, lapply(grep("^Rank at ", names(d), value = TRUE), function(col)
      data.frame(Variable = d$Variable, Site = sub("^Rank at ", "", col), Rank = d[[col]])))
    long$Variable <- factor(long$Variable, levels = rev(d$Variable))
    fed <- data.frame(Variable = factor(d$Variable, levels = rev(d$Variable)), Rank = d$`Federated rank`)
    ggplot(long, aes(Rank, Variable)) +
      geom_col(data = fed, fill = "#2c3e50", alpha = .15, width = .7) +
      geom_point(aes(colour = Site), size = 3) +
      scale_x_continuous("Rank (1 = most important)", breaks = scales::breaks_pretty()) +
      labs(y = NULL, caption = "Bars: federated (weighted mean) rank") +
      theme_minimal(base_size = 14) + theme(legend.position = "top")
  })

  # Step 4: cut-offs ---------------------------------------------------------------
  observeEvent(input$start_cut, {
    probs <- run_safely({
      q <- parse_cuts(input$quantiles)
      if (any(q <= 0 | q >= 1)) stop("Quantiles must be strictly between 0 and 1.")
      q
    })
    req(probs)
    start("cut", function(st) cut_request(st, probs), "summary")
  })

  cut_order <- reactive({
    req(results$cut)
    vars <- names(results$cut$cutoffs)
    if (!is.null(results$rank)) intersect(as_df(results$rank$table)$Variable, vars) else vars
  })

  observe({
    req(results$cut)
    num <- Filter(function(v) results$cut$cutoffs[[v]]$type == "numeric", cut_order())
    updateSelectInput(session, "tune_var", choices = num, selected = isolate(input$tune_var))
  })
  observe({
    req(results$cut, input$tune_var %in% names(results$cut$cutoffs))
    updateTextInput(session, "tune_cuts", value = paste(format_num(results$cut$cutoffs[[input$tune_var]]$cuts), collapse = ", "))
  })
  publish <- function(cutoffs, note) {
    guard(is_lead(), "Only the lead site can change the cut-offs.")
    ok <- run_safely({
      publish_cutoffs(rv$store, state(), names(snap()), cutoffs)
      TRUE
    })
    if (isTRUE(ok)) {
      showNotification(note, type = "message")
      refresh_now()
    }
  }
  observeEvent(input$tune_apply, {
    req(results$cut, input$tune_var)
    cuts <- run_safely(parse_cuts(input$tune_cuts))
    req(cuts)
    new <- results$cut$cutoffs
    new[[input$tune_var]]$cuts <- cuts
    publish(new, sprintf("Published new cut-offs for %s.", input$tune_var))
  })
  observeEvent(input$tune_reset, {
    req(results$cut)
    publish(results$cut$auto, "Published the original cut-offs.")
  })

  output$cut_table <- renderTable({
    validate(need(results$cut, "No unified cut-offs yet."))
    cuts <- results$cut$cutoffs
    tbl <- cutoff_table(cuts, cut_order())
    edited <- vapply(tbl$Variable, function(v) !identical(cuts[[v]], results$cut$auto[[v]]), logical(1))
    tbl$Variable[edited] <- paste(tbl$Variable[edited], "(edited)")
    tbl
  })

  # Step 5: model selection ----------------------------------------------------------
  observeEvent(input$start_pars, {
    guard(!is.null(results$rank), "Finish Ranking first.")
    vars <- head(ranking()$Variable, input$max_vars)
    start("pars", function(st) fit_request(st, lapply(seq_along(vars), function(m) vars[seq_len(m)]),
                                           input$method, input$refine, "validate", 100, step_data(st, "pars", rv$data)),
          c("rank", "cut"))
  })

  auc_label <- reactive(if (outcome() == "survival") "iAUC" else "AUC")

  output$pars_plot <- renderPlot({
    validate(need(results$pars, "No model selection results yet."))
    sm <- as_df(results$pars$summary)
    raw <- as_df(results$pars$auc)
    lab <- factor(sm$Variable, levels = sm$Variable)
    raw$Variable <- factor(raw$Variable, levels = sm$Variable)
    lo <- max(0.5, floor(min(raw$AUC) * 20) / 20 - 0.05)
    ggplot(data.frame(Variable = lab, AUC = sm$AUC, Model = sm$Model), aes(Variable, AUC)) +
      geom_col(fill = "#2c3e50", width = .7) +
      geom_text(aes(label = Model), vjust = 1.6, colour = "white", size = 4) +
      geom_point(data = raw, aes(colour = Site), size = 2.5, position = position_dodge(width = .5)) +
      coord_cartesian(ylim = c(lo, 1)) +
      labs(x = "Variable added (in federated rank order)", y = paste("Validation", auc_label()),
           title = sprintf("%s: bars = weighted mean across sites", results$pars$method)) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
  })

  output$pars_table <- renderTable({
    req(results$pars)
    raw <- as_df(results$pars$auc)
    wide <- reshape(raw[, c("Model", "Site", "AUC")], idvar = "Model", timevar = "Site", direction = "wide")
    names(wide) <- sub("^AUC\\.", paste(auc_label(), "at "), names(wide))
    out <- merge(as_df(results$pars$summary), wide, by = "Model")
    names(out)[names(out) == "AUC"] <- auc_label()
    out
  }, digits = 3)

  # Step 6: final score ---------------------------------------------------------------
  output$custom_vars <- renderUI({
    validate(need(results$rank, "Finish the ranking first."))
    r <- ranking()
    checkboxGroupInput("custom_sel", NULL, selected = isolate(input$custom_sel),
                       choiceNames = sprintf("%s (rank %.1f)", r$Variable, r$`Federated rank`), choiceValues = r$Variable)
  })

  output$final_method <- renderText({
    if (input$method %in% c("dCLR", "ODACH")) sprintf("Algorithm: %s with %d surrogate pass(es) (set on step 5).", input$method, input$refine)
    else sprintf("Algorithm: %s (set on step 5).", input$method)
  })

  observeEvent(input$start_final, {
    guard(!is.null(results$rank), "Finish Ranking first.")
    all_vars <- ranking()$Variable
    vars <- if (input$var_mode == "top") head(all_vars, input$n_vars) else all_vars[all_vars %in% input$custom_sel]
    guard(length(vars) > 0, "Select at least one variable.")
    times <- NULL
    if (outcome() == "survival") {
      times <- if (nzchar(trimws(input$times))) run_safely(parse_cuts(input$times)) else if (!is.null(rv$data)) {
        tr <- rv$data$train
        unname(signif(stats::quantile(tr$label_time[tr$label_status == 1], c(.25, .5, .75)), 3))
      }
      guard(length(times) > 0, "Enter the times for AUC(t), or load your data so the app can choose them.")
    }
    start("final", function(st) fit_request(st, list(vars), input$method, input$refine, "test", input$max_score,
                                            step_data(st, "final", rv$data), times), c("rank", "cut"))
  })

  final <- reactive({
    validate(need(results$final, "No final score yet."))
    results$final
  })
  final_table <- reactive(as_df(final()$table))

  output$score_tbl <- renderTable({
    tbl <- final_table()
    tbl$Variable[duplicated(tbl$Variable)] <- ""
    tbl
  }, digits = 0)

  output$perf_tbl <- renderTable(perf_table(final()), digits = 3)

  # Test performance shared by every site, formatted for display.
  perf_table <- function(f) {
    out <- as_df(f$perf)
    if (!identical(f$outcome, "survival")) {
      out$`AUC (95% CI)` <- sprintf("%.3f (%.3f to %.3f)", out$AUC, out$Lower, out$Upper)
      return(out[, c("Site", "N", "Events", "AUC (95% CI)", "Threshold", "Sensitivity", "Specificity")])
    }
    out$`C-index (95% CI)` <- sprintf("%.3f (%.3f to %.3f)", out$C_index, out$Lower, out$Upper)
    out <- out[, c("Site", "N", "Events", "iAUC", "C-index (95% CI)")]
    for (s in seq_len(nrow(out))) {
      a <- f$auc_t[[out$Site[s]]]
      for (k in seq_along(unlist(a$time))) out[s, sprintf("AUC(%s)", format_num(unlist(a$time)[k]))] <- unlist(a$auc)[k]
    }
    out
  }

  output$perf_plot <- renderPlot({
    f <- final()
    perf <- as_df(f$perf)
    if (identical(f$outcome, "survival")) {
      d <- do.call(rbind, lapply(names(f$auc_t), function(s)
        data.frame(Site = sprintf("%s (iAUC %.3f)", s, perf$iAUC[perf$Site == s]),
                   Time = unlist(f$auc_t[[s]]$time), AUC = unlist(f$auc_t[[s]]$auc))))
      validate(need(nrow(d) > 0, "No AUC(t) values: the chosen times are beyond the follow-up of the test sets."))
      return(ggplot(d, aes(Time, AUC, colour = Site)) + geom_line(linewidth = 1) + geom_point(size = 3) +
               coord_cartesian(ylim = c(0.5, 1)) +
               labs(x = "Time", y = "AUC(t)", title = "Time-dependent AUC on each site's test set") +
               theme_minimal(base_size = 14) + theme(legend.position = "bottom", legend.direction = "vertical"))
    }
    d <- do.call(rbind, lapply(names(f$roc), function(s) {
      data.frame(Site = sprintf("%s (AUC %.3f)", s, perf$AUC[perf$Site == s]),
                 FPR = unlist(f$roc[[s]]$fpr), TPR = unlist(f$roc[[s]]$tpr))
    }))
    ggplot(d, aes(FPR, TPR, colour = Site)) +
      geom_abline(linetype = 2, colour = "grey60") + geom_path(linewidth = 1) +
      coord_equal() + labs(x = "1 - Specificity", y = "Sensitivity", title = "ROC on each site's test set") +
      theme_minimal(base_size = 14) + theme(legend.position = "bottom", legend.direction = "vertical")
  })

  output$coef_tbl <- renderTable({
    f <- final()
    scale <- switch(f$method, ODACH = "Federated (log hazard ratio)", DiSAH = "Federated (hazard difference)",
                    "Federated (log odds ratio)")
    out <- data.frame(Term = f$terms, Federated = f$beta, SE = f$se, check.names = FALSE)
    names(out)[2] <- scale
    if (length(f$init)) out[[if (f$method == "ODAL2") "Initial (lead site)" else "Initial (round 1 average)"]] <- f$init
    out
  }, digits = 4)

  my_test_scores <- reactive({
    f <- final()
    test <- my_data()$test
    data.frame(score = compute_scores(test, final_table(), f$cutoffs), test[, intersect(OUTCOME_COLUMNS, names(test)), drop = FALSE])
  })

  output$dl_score <- downloadHandler(
    filename = function() sprintf("fedscore_table_%s.csv", final()$method),
    content = function(file) write.csv(final_table(), file, row.names = FALSE)
  )
  output$dl_pred <- downloadHandler(
    filename = function() sprintf("fedscore_test_scores_%s.csv", rv$me),
    content = function(file) write.csv(my_test_scores(), file, row.names = FALSE)
  )

  # Step 7: TLScore (local only) -------------------------------------------------------
  observeEvent(results$final, rv$tl <- NULL, ignoreNULL = FALSE)
  observeEvent(input$run_tl, {
    guard(!is.null(results$final), "The final score (step 6) is not ready yet.")
    guard(!is.null(rv$data), "Load your data on the Study tab first.")
    f <- results$final
    guard(f$method %in% TL_METHODS, sprintf("TLScore needs a Cox (ODACH) or logistic model; this score uses %s.", f$method))
    d <- my_data()
    rv$tl <- run_safely(withProgress(message = "Fitting TLScore", value = .5, {
      tl <- fit_tlscore(setNames(f$beta, f$terms), f$method, f$vars, d, f$cutoffs, f$max_score)
      evaluate <- function(score) {
        if (identical(f$outcome, "survival"))
          evaluate_survival(score, d$test$label_time, d$test$label_status, numeric(0))$metrics
        else evaluate_scores(score, d$test$label)
      }
      tl$perf <- rbind(cbind(Score = "FedScore", evaluate(my_test_scores()$score)),
                       cbind(Score = "TLScore", evaluate(compute_scores(d$test, tl$table, f$cutoffs))))
      tl
    }))
  })

  output$tl_perf <- renderTable({
    validate(need(is.null(results$final) || results$final$method %in% TL_METHODS,
                  "TLScore is not available for additive-hazards (DiSAH) scores; it needs a logistic or Cox model."))
    validate(need(rv$tl, "When the final score is ready, fit TLScore for your site."))
    p <- rv$tl$perf
    adjusted <- c("", sprintf("%d of %d", sum(rv$tl$delta != 0), length(rv$tl$delta)))
    if ("C_index" %in% names(p)) {
      p$`C-index (95% CI)` <- sprintf("%.3f (%.3f to %.3f)", p$C_index, p$Lower, p$Upper)
      return(cbind(Site = rv$me, p[, c("Score", "N", "Events", "iAUC", "C-index (95% CI)")], `Coefficients adjusted` = adjusted))
    }
    p$`AUC (95% CI)` <- sprintf("%.3f (%.3f to %.3f)", p$AUC, p$Lower, p$Upper)
    cbind(Site = rv$me, p[, c("Score", "N", "Events", "AUC (95% CI)", "Sensitivity", "Specificity")],
          `Coefficients adjusted` = adjusted)
  }, digits = 3)

  tl_compare <- reactive({
    req(rv$tl)
    tbl <- final_table()
    data.frame(Variable = tbl$Variable, Category = tbl$Category, FedScore = tbl$Points, TLScore = rv$tl$table$Points)
  })
  output$tl_tbl <- renderTable({
    tbl <- tl_compare()
    tbl$Variable[duplicated(tbl$Variable)] <- ""
    tbl
  }, digits = 0)
  output$dl_tl <- downloadHandler(
    filename = function() sprintf("tlscore_table_%s.csv", rv$me),
    content = function(file) write.csv(tl_compare(), file, row.names = FALSE)
  )

  # Exchange log -------------------------------------------------------------------------
  log_df <- reactive({
    req(rv$store)
    exchange_log(snap(), rv$study)
  })
  output$log_tbl <- renderTable({
    validate(need(length(snap()) > 0, "Nothing has been exchanged yet."))
    log_df()
  }, na = "")
  observe({
    files <- if (length(snap())) log_df()$File else character(0)
    updateSelectInput(session, "log_file", choices = files, selected = isolate(input$log_file))
  })
  output$log_json <- renderText({
    req(input$log_file, rv$store)
    obj <- snap()[[study_file(rv$study, input$log_file)]]
    req(obj)
    txt <- jsonlite::prettify(to_json(obj), indent = 2)
    if (nchar(txt) > 50000) paste0(substr(txt, 1, 50000), "\n... (truncated; download the full log)") else txt
  })
  output$dl_log <- downloadHandler(
    filename = function() sprintf("fedscore_exchange_%s.json", rv$study),
    content = function(file) writeLines(to_json(snap()), file)
  )
}

shinyApp(ui, server)
