# FedScore Shiny app

An app for building a federated scoring system with FedScore, for a **binary** outcome or
a **survival** (time-to-event) outcome. It runs the whole workflow, from an optional
federated imputation of missing values to a final points table:

| Step | Binary outcome | Survival outcome |
|---|---|---|
| Missing values (optional) | FedIMPUTE | FedIMPUTE |
| Ranking | random forest or univariable AUC | random survival forest |
| Cut-offs | federated quantiles | federated quantiles |
| Model selection and final score | **ODAL2** (homogeneous sites) or **dCLR** (heterogeneous sites) | **ODACH** (Cox) or **DiSAH** (additive hazards), both with site-specific baseline hazards |
| Evaluation | AUC with CI, ROC curve | integrated AUC, C-index with CI, AUC(t) |
| TLScore (optional, local) | yes | ODACH scores |

The survival workflow follows AutoScore-Survival.

Every site runs the app on its own computer with its own data. Patient-level data never
leaves the site. Sites exchange only summary statistics, as small JSON files in a
**shared store**, which can be:

* a **GitHub repository**: a public repository lets anyone audit exactly what was
  shared; or
* a **shared folder**: a network drive, or a Dropbox, OneDrive or Box folder.

One site, the **lead**, creates the study, starts each step and combines the replies.
The other sites answer its requests.

The method is described in the FedScore paper (Li et al., *Journal of Biomedical
Informatics* 2023, <https://doi.org/10.1016/j.jbi.2023.104485>). The other methods are
cited under [Methods](#methods).

## Running the app

```r
# From the folder containing this repository
shiny::runApp("FedScore-ShinyApp")
```

Required packages: `shiny`, `bslib`, `ggplot2`, `AutoScore`, `pROC`, `glmnet`, `jsonlite`,
`mice`, `digest`, and `httr2` for the GitHub store. AutoScore brings in `randomForest`,
`randomForestSRC`, `survival` and `survAUC`. `MASS` and `scales` ship with standard R/ggplot2
installs.

```r
install.packages(c("shiny", "bslib", "ggplot2", "AutoScore", "pROC", "glmnet", "jsonlite",
                   "mice", "digest", "httr2"))
```

There is no C code to compile and no Python dependency. Tested with R 4.4, with both
shiny 1.10 / bslib 0.8 and shiny 1.13 / bslib 0.10. Each site should run the app
**locally**: deploying it to a hosting service such as shinyapps.io would upload the
site's data to that service.

### Try it on one computer

1. Start the app and open it in two browser tabs.
2. In tab 1, pick a site-1 demo (see below), keep the default folder and study name, and
   click **Connect**. Then enter `site1, site2` as the sites, check the outcome, and click
   **Create study**.
3. In tab 2, pick the matching site-2 demo and click **Connect**.
4. In tab 1 (the lead), start each step in turn. Tab 2 replies automatically.

There are three demo pairs:

* **Binary outcome:** the FedScore binary demo.
* **Survival outcome:** the FedScore heterogeneous survival demo.
* **Binary, with missing values:** the binary demo with values removed. Start with the
  Missing values step.

`data-raw/make_demo_data.R` shows how each demo file was made.

## Running a real study

### 1. Set up the store

**GitHub.** Create a repository with at least one commit (for example, one with a
README). Each site needs a token that can write to it:

* A **fine-grained personal access token**, limited to that repository, with
  *Contents: Read and write* permission.
* Collaborators need write access to the repository.

Enter the token in the app, or set `GITHUB_PAT` before starting R. The token is used only
by the app on that computer and is never written to the store. All study files go in one
folder of the repository (`fedscore` by default). Every exchange becomes a commit, so the
history is a complete audit trail.

**Shared folder.** Every site enters the path at which the folder appears on its own
computer.

### 2. Agree on names

Agree on a study name and on one name per site (letters, digits, `-` and `_`).

To prefill the connection fields, set these environment variables, for example in
`.Renviron`:

* `FEDSCORE_STORE` (`folder` or `github`)
* `FEDSCORE_FOLDER`
* `FEDSCORE_REPO`, `FEDSCORE_BRANCH` and `FEDSCORE_DIR`
* `FEDSCORE_STUDY` and `FEDSCORE_SITE`

### 3. Each site connects

Each site loads its data, enters the store details, the study name and its site name, and clicks **Connect**.

The lead connects first, then clicks **Create study**. It lists every site, chooses the
outcome (binary or survival), and sets the minimum count (see [Privacy](#privacy)).

### 4. The lead runs the steps

The lead starts each step and chooses its settings (variables to impute, ranking method,
quantiles, algorithm, number of variables, and so on). The Missing values step is
optional. Skip it, and missing values become an "Unknown" category of their variable. With *Reply and continue automatically*
ticked (the default), each site's app answers new requests as they appear, and the lead's
app moves on once every site has replied. Untick it to review each request and click
**Send my reply**, or at the lead, **Combine replies and continue**.

Each step's **Progress** card shows every round and which sites have replied. A site
whose data cannot answer a request, for example because its test set has no events,
replies with an error message. The lead then waits until that site fixes its data and
sends its reply again.

The app checks the store every 3 seconds for a shared folder, or every 10 seconds for
GitHub, and **Refresh** checks immediately. Re-running a step clears all later steps for
every site.

### 5. Results

Every site sees the federated ranking, cut-offs, parsimony plot, score table, coefficients
and every site's test performance. Each site can also run **TLScore** on its own data;
this step shares nothing.

## What is shared

| Step | What each site computes locally | What it shares |
|---|---|---|
| 1 · Study | checks its data | training/validation/test sizes, event counts, predictor names, types and categories, number of missing values per predictor |
| 2 · Missing values | completes its predictors with MICE; fits and later applies the imputation models | coefficients, then gradients and Hessians, of each imputation model (never imputed values) |
| 3 · Ranking | variable importance (random forest, random survival forest or univariable AUC) | rank position of each variable |
| 4 · Cut-offs | quantiles of each continuous variable | those quantiles; the levels of categorical variables |
| 5 · Model selection | federated fit for the top 1…m variables, then validation AUC (binary) or iAUC (survival) | estimates, variances, gradients and Hessians, or DiSAH's summed estimating equations (see below); validation AUC/iAUC |
| 6 · Final score | federated fit for the chosen variables, then test performance | the same for one model; binary: test AUC with CI, best operating point, ROC curve; survival: iAUC, C-index with CI, AUC at the chosen times |
| 7 · TLScore | lasso correction on its own data | nothing |

The **Exchange log** tab lists every file in the study, with the number of numeric values
in each, and shows any file's contents.

### Privacy

Summary statistics are much less sensitive than patient records, but they are not
risk-free. Statistics computed from very few patients can come close to describing
individuals. With categorical predictors, for example, a site's Hessian roughly encodes how
many of its patients fall in each category.

* **Minimum count.** The lead sets this when creating the study (default 10; 0 turns it
  off). A site refuses to share model statistics if any non-empty category of the model's
  variables has fewer patients than this in its training data. It also refuses to share
  validation or test performance if that set has fewer events or non-events than this. The
  refusal names the categories involved. The lead can merge sparse numeric categories by
  editing the cut-offs on step 4; sparse categorical levels must be merged in the site's
  data. For an imputation model, a site with fewer observed values than this simply does
  not take part in fitting that model.
* **Where the files live.** In a **public** GitHub repository, everyone can read every
  shared statistic. That suits demonstrations and open, audited studies. For real patient
  data, use a **private** repository that only the study sites can access, or a shared
  folder. Check this against your institution's data-governance and ethics approvals.
* **What is never shared:** patient-level rows, imputed values, the access token, and
  anything from TLScore.

### Files in the store

```
<study>/study.json                   lead site, sites (written by the lead)
<study>/<step>/round<k>.json         the lead's request k of a step
<study>/<step>/round<k>_<site>.json  each site's reply to request k
<study>/<step>/result.json           the outcome of the step (written by the lead)
```

Each site writes only its own reply files, so sites never overwrite each other's files.
Every request carries the id of its step run. Replies to an older run are ignored.

### Methods

#### Binary outcomes

* **ODAL2** (Duan et al., *JAMIA* 2020). The lead site fits a local logistic regression,
  and its estimate β̄ is broadcast in round 1. Every site returns its gradient and Hessian
  at β̄. The lead site maximises the second-order surrogate likelihood, and round 2
  evaluates the result.
* **dCLR** (Tong et al., *npj Digit Med* 2022). This uses each site's pairwise conditional
  likelihood, so site-specific baseline risks cancel:
  * Round 1: each site shares its local estimate and sandwich variance. These are combined
    into β̄ by inverse-variance weighting.
  * Round 2: each site shares its derivatives at β̄.
  * Round 3: each site maximises its surrogate likelihood. The surrogate estimates are
    combined by inverse-variance weighting.

  *Surrogate passes* (step 5) repeats rounds 2–3 around the updated estimate. One pass is
  the published one-shot dCLR. Extra passes move the result toward the pooled-data fit,
  which helps when sites are very heterogeneous. In the demo data, site 2 has no events in
  the youngest age group. One pass leaves Age coefficients about 0.24 from the pooled fit;
  two passes bring the gap to about 0.001.

#### Survival outcomes

Both survival models stratify by site: each site keeps its own baseline hazard, and the
covariate effects are shared.

* **ODACH** (Luo et al., *Sci Rep* 2022). This is the Cox model, using each site's Breslow
  partial likelihood. Its rounds are the same as dCLR's:
  1. Local Cox fits, combined into β̄ by inverse-variance weighting.
  2. Gradients and Hessians at β̄.
  3. Each site's surrogate estimate, combined by inverse-variance weighting.

  Surrogate passes work as for dCLR. Coefficients are log hazard ratios.
* **DiSAH** (Wang et al., arXiv:2601.14609, 2026). This is the additive hazards model
  λ_k(t) + xᵀβ, using Lin and Ying's estimating equations.
  * In one round, each site shares three sums: `A`, `d` and `B`.
  * The lead computes β = (ΣA)⁻¹Σd and the sandwich variance. This is exactly the
    pooled-data stratified estimate.
  * Coefficients are hazard differences per unit of time.

  The app implements DiSAH's stratified estimator. The unstratified estimator needs every
  patient's event time from every site, which this app does not share.

For both, the parsimony analysis uses the integrated AUC on the validation sets, as
AutoScore-Survival does. The final score reports:
* iAUC;
* Harrell's C-index with a 95% CI;
* Uno's AUC(t) at times chosen by the lead. By default these are the quartiles of the lead
  site's event times.

#### Missing values: FedIMPUTE

This is FedIMPUTE (Li et al., *J Biomed Inform* 2025), the DAC-MICE variant:

1. Each site completes its predictors with one MICE imputation. This uses only its own
   data, and never the outcome.
2. Each variable to impute gets an imputation model on the other predictors: linear for a
   numeric variable, logistic for a binary one. Each model is federated:
   * Round 1: each site shares local coefficients, and the lead averages them.
   * Each further *Newton round*: each site shares gradients and Hessians at the current
     coefficients, and the lead takes a Newton step.
3. Each site fills its own missing values with the federated models:
   * the predicted value for numeric variables;
   * a random draw for binary ones.

   Optionally, each site first recalibrates the models on its own data (*transfer
   learning*: a lasso correction with the federated prediction as offset).

Later steps use the completed data at every site. Variables with more than two categories
are not imputed; their missing values remain an "Unknown" category.

Differences from the FedIMPUTE code:
* Sites are weighted by their number of observed values, so one Newton round reproduces the
  pooled-data linear model exactly. FedIMPUTE averages sites equally.
* The adaptive-lasso step after several DAC iterations is not included.
* The transfer-learning recalibration uses a cross-validated lasso instead of FedIMPUTE's
  thresholding grid.

#### Site weights

Derivatives are aggregated with the site weights chosen when the study is created (by
training size or equal).

### Score table

Within each variable the lowest-risk category gets 0 points. Points are scaled so the
maximum attainable total equals *Maximum total score* (default 100), as in AutoScore.

## Site data format

A site's data is one `.rds` file, which can be:

* a list with `train`, `validate` and `test` data frames (the names `train_set`,
  `validation_set`, `validation` and `test_set` are also accepted);
* a list with a `non_cv` element in that format (as in `data/binary/` of the FedScore repository); or
* a single data frame, which is split 70/10/20.

Every split needs the outcome:

* **Binary:** a `label` column, coded as logical, 0/1 or a two-level factor.
* **Survival:** `label_time`, the follow-up time (≥ 0), and `label_status`, the event
  indicator (logical or 0/1). These are AutoScore-Survival's column names.

The other columns are predictors and may contain missing values. Character columns are
treated as categorical. Only columns present at every site, with the same type, are used.
Demo data is in `data/`.

## Code layout

| File | Contents |
|---|---|
| `app.R` | user interface and server |
| `R/fed_store.R` | folder and GitHub stores |
| `R/fed_protocol.R` | requests, replies and how the lead combines them |
| `R/fed_data.R` | loading and checking site data |
| `R/fed_prep.R` | ranking and cut-offs |
| `R/fed_impute.R` | FedIMPUTE: local completion, imputation models, filling missing values |
| `R/fed_model.R` | likelihoods and the ODAL2, dCLR, ODACH and DiSAH computations |
| `R/fed_score.R` | score tables, evaluation (binary and survival) and TLScore |
| `data-raw/make_demo_data.R` | how the demo data was made |

## Size and speed

A study with 5 sites and all 20 demo variables in the parsimony analysis took the
following time, added up over all sites. Real sites work in parallel, so each spends
about a fifth of this:

* ODAL2: about 20 seconds in total.
* dCLR (one pass): about 2.5 minutes in total.

The largest file in that study was about 1 MB. GitHub's blob API, which the app uses to
download files, supports files up to 100 MB.

## Tests

```r
testthat::test_dir("FedScore-ShinyApp/tests/testthat")
```

The tests check that:

* the vectorised pairwise likelihood, gradient and Hessian match a line-by-line R port of
  the original FedScore C code (`examples/binary/hete_demo/C` in the FedScore repository)
  to within 1e-10;
* both stores round-trip JSON. The GitHub store is tested against a mock of the GitHub API,
  including a write conflict between two sites;
* the Cox partial likelihood matches `survival::coxph` (Breslow ties), and the additive
  hazards sums match a port of DiSAH's own `V_k()` code, to within 1e-8;
* the survival metrics match AutoScore-Survival's;
* complete studies run through a store, and match the pooled-data fits:
  * ODAL2 and dCLR estimates match pooled logistic and pairwise fits;
  * ODACH matches pooled stratified `coxph` (within 0.001);
  * DiSAH matches the pooled stratified estimate exactly;
* FedIMPUTE's linear imputation models equal pooled least squares after one Newton round,
  and the imputed data is used by every later step;
* error replies, stale replies and restarted steps are handled;
* sites refuse to share statistics about groups smaller than the minimum count; and
* the app runs whole studies as the lead site, with a second site replying through the same
  folder (`shiny::testServer`): binary with ODAL2 and dCLR, survival with ODACH, and a study
  that imputes missing values first. A non-lead site can only reply.

The tests run on GitHub Actions for every push and pull request.

## Changes from the earlier prototype

This app replaces the prototype scripts (`shinyapp*.R`) on the `ShinyApp` branch of the
FedScore repository:

* **Sites now really are separate.** The prototype loaded every site's data into one R
  session and wrote intermediate files to its own working directory. Now each site runs
  its own copy and exchanges JSON through a shared store.
* **dCLR did not federate.** `get_meanG()` read the lead site's gradient file K times instead
  of averaging across sites. Each site also evaluated its gradient at its own estimate
  rather than at a shared β̄.
* **ODAL2 omitted the lead site from the global gradient and Hessian.**
* **The score table subtracted one global minimum across all coefficients, including the
  intercept.** It now shifts and scales points within each variable, AutoScore-style.
* **Unified cut-offs were floored to integers** (e.g. 36.8 °C became 36), which collapsed
  categories. They are now rounded to 3 significant digits. Cut-offs that equal a site's
  minimum or maximum are dropped, so no category is empty.
* **Sites could end up with different factor levels.** Every site now uses identical
  category levels, so design matrices always align.
* **TLScore read `site_beta.json`, which the dCLR path never writes, and added the full local
  β to the federated β.** It now fits a lasso correction δ with the federated linear
  predictor as an offset (β_TL = β_fed + δ).
* **The parsimony plot never scored the validation sets, and the "Run Regression" button did
  nothing.**
* **The C likelihoods (O(n²) loops, compiled by hand) and Nelder–Mead are gone.** Exact
  vectorised derivatives with Newton steps replace them. Rows are collapsed to unique
  covariate patterns, so a 10-variable dCLR fit on the demo data takes seconds.
* **The prototype required RStudio (`getSourceEditorContext`) and Python.**

The survival and imputation code was written from the papers and checked against the
DiSAH and FedIMPUTE code, not copied from it. Compared with FedScore's own survival scripts
(`code/R/fedscore_survival.R`), the surrogate step averages all sites' derivatives. Those
scripts read one site's file once per site, the same bug the old dCLR code had.

## License

GPL-3, the same as AutoScore and shiny (see `LICENSE`).

## Citation

Li, S., Ning, Y., Ong, M.E., Chakraborty, B., Hong, C., Xie, F., ... & Liu, N. (2023). FedScore: A privacy-preserving framework for federated scoring system development. *Journal of Biomedical Informatics*, 104485. https://doi.org/10.1016/j.jbi.2023.104485

Methods used in the app:

* Duan, R., et al. (2020). Learning from electronic health records across multiple sites: a communication-efficient and privacy-preserving distributed algorithm (ODAL2). *JAMIA*, 27(3), 376–385. https://doi.org/10.1093/jamia/ocz199
* Tong, J., et al. (2022). Distributed learning for heterogeneous clinical data with application to integrating COVID-19 data across 230 sites (dCLR). *npj Digital Medicine*, 5, 76. https://doi.org/10.1038/s41746-022-00615-8
* Luo, C., et al. (2022). ODACH: a one-shot distributed algorithm for Cox model with heterogeneous multi-center data. *Scientific Reports*, 12, 6627. https://doi.org/10.1038/s41598-022-09069-0
* Wang, Z., Li, S., Ong, M.E.H., & Liu, N. (2026). Communication-efficient distributed hazard difference estimation for heterogeneous multi-site survival data (DiSAH). arXiv:2601.14609.
* Li, S., et al. (2025). FedIMPUTE: Privacy-preserving missing value imputation for multi-site heterogeneous electronic health records. *Journal of Biomedical Informatics*. https://doi.org/10.1016/j.jbi.2025.104780
