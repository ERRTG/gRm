#' Exogenous score-effect diagnostics
#'
#' @param analysis A DIGRAM analysis object.
#' @param inference Inference mode. `"asymptotic"` maps to DIGRAM command 17,
#'   `"exact"` to command 2, and `"repeated"` to command 74.
#' @param nsim Number of simulations for exact branches.
#' @param seed Random seed for exact branches.
#' @param critlevel Optional DIGRAM exact command critical level on the
#'   per-1000 scale.
#' @param risk Optional DIGRAM repeated-Monte-Carlo risk on the per-1000 scale.
#' @param score_cap Highest score category retained before source collapsing.
#' @param jobs Reserved for future parallel implementations.
#' @param ... Reserved for S3 dispatch compatibility; ignored.
#' @return A `gRm_score_effects` result object.
#' @export
#' @examples
#' \donttest{
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0),
#'   site = c(0, 0, 1, 1, 0, 1, 0, 1)
#' )
#' analysis <- gRm(data, items = c("I1", "I2"), exogenous = "site", id = "ID")
#' effects <- score_effects(analysis)
#' summary(effects, which = "tests")
#' }
score_effects <- function(analysis,
                          inference = c("asymptotic", "exact", "repeated"),
                          nsim = 1000L,
                          seed = 9L,
                          critlevel = NULL,
                          risk = NULL,
                          score_cap = 56L,
                          jobs = 1L,
                          ...) {
  analysis <- as_public_gRm_analysis(analysis)
  inference <- match.arg(inference)
  exact_state <- gRm_exact_command_state_public(
    inference,
    nsim = nsim,
    seed = seed,
    critlevel = critlevel,
    risk = risk
  )
  jobs <- normalize_public_jobs(jobs)
  values <- exo_select_values(
    analysis$project,
    score_cap = score_cap,
    exact = exact_state$exact,
    repeated = exact_state$sequential,
    nsim = exact_state$nsim,
    seed = exact_state$seed,
    exact_state = exact_state
  )
  new_gRm_result(
    class = "gRm_score_effects",
    analysis = analysis,
    fit = NULL,
    values = values,
    result = "score_effects",
    metadata = list(
      inference = inference,
      nsim = exact_state$nsim,
      seed = if (exact_state$exact) exact_state$seed else NA_integer_,
      critlevel = if (is.null(critlevel)) NA_integer_ else as.integer(critlevel),
      risk = if (is.null(risk)) NA_integer_ else as.integer(risk),
      score_cap = as.integer(score_cap),
      jobs = jobs,
      exact_state = exact_state
    ),
    call = match.call()
  )
}

#' Item parameter estimates
#'
#' @param fit A fitted DIGRAM model.
#' @param ... Reserved for S3 dispatch compatibility; ignored.
#' @return A `gRm_item_parameters` result object.
#' @export
#' @examples
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0)
#' )
#' fit0 <- fit(gllrm(gRm(data, items = c("I1", "I2"), id = "ID")))
#' params <- item_parameters(fit0)
#' summary(params, which = "coefficients")
item_parameters <- function(fit, ...) {
  fit <- as_public_gRm_fit(fit)
  new_gRm_result(
    class = "gRm_item_parameters",
    analysis = fit$analysis %||% fit$spec$analysis,
    fit = fit,
    values = fit$values,
    result = "item_parameters",
    call = match.call()
  )
}

#' Item fit diagnostics
#'
#' Compute item-fit statistics for a fitted gRm model. Use
#' `summary(result, which = "tests")` for compact outfit, infit, and
#' item-restscore gamma tests. Use `summary(result, which = "items")` for the
#' per-item outfit/infit diagnostic summary produced by the extended item-fit
#' calculation.
#'
#' @param fit A fitted gRm model.
#' @param include_extended Whether to compute extended item-fit detail tables.
#'   When `FALSE`, `summary(result, which = "items")` returns an empty data
#'   frame because the per-item extended diagnostic summaries were not
#'   computed. The compact `summary(result, which = "tests")` table remains
#'   available.
#' @param ... Reserved for S3 dispatch compatibility; ignored.
#' @return A `gRm_item_fit` result object.
#' @details
#' The item-fit result has two different public summary layers. The
#' `summary(result, which = "tests")` layer is always available and returns one
#' row per item with the compact inferential statistics: outfit, infit,
#' item-restscore gamma, standard errors, p-values, FDR labels, and direction.
#' Its columns are `item_label`, `item_name`, `outfit`, `outfit_sd`,
#' `p_outfit`, `outfit_fdr`, `infit`, `infit_sd`, `p_infit`, `infit_fdr`,
#' `observed_gamma`, `expected_gamma`, `gamma_sd`, `p_gamma`, `gamma_fdr`,
#' and `direction`.
#'
#' The `summary(result, which = "items")` layer is also one row per item, but it
#' is not a second test table. It returns the extended outfit/infit
#' decomposition: aggregate observed and expected outfit components, the
#' resulting outfit value, and the observed, expected, variance, and ratio
#' components for infit. Its columns are `item_label`, `item_name`,
#' `outfit_total_n`, `outfit_total_observed`, `outfit_total_expected`,
#' `outfit_total_value`, `infit_observed`, `infit_expected`,
#' `infit_variance`, and `infit_value`.
#'
#' `summary(result, which = "bh")` returns the Benjamini-Hochberg critical
#' p-value limits used for the combined compact item-fit test family. Its
#' columns are `threshold` and `p_value`.
#'
#' The item-restscore gamma standard error follows the DIGRAM source
#' convention. For each item, gRm builds observed and fitted item-by-restscore
#' tables. The fitted table supplies both the expected gamma and the reference
#' spread used to place the observed gamma on a standard-error scale. The
#' reported `gamma_sd` and `p_gamma` therefore describe the DIGRAM
#' table-based comparison between observed and expected gamma; they are not
#' recalculated as a separate sample-only or resampling standard error.
#'
#' The `"items"` layer depends on the extended calculation. With the default
#' `include_extended = TRUE`, `item_fit()` computes that layer and
#' `summary(result, which = "items")` is populated. With
#' `include_extended = FALSE`, the extended calculation is skipped and
#' `summary(result, which = "items")` returns an empty data frame. The compact
#' `summary(result, which = "tests")` and `summary(result, which = "bh")`
#' sections remain available.
#' @export
#' @examples
#' \donttest{
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0)
#' )
#' fit0 <- fit(gllrm(gRm(data, items = c("I1", "I2"), id = "ID")))
#' item_tests <- item_fit(fit0)
#' summary(item_tests, which = "tests")
#' summary(item_tests, which = "items")
#' }
item_fit <- function(fit, include_extended = TRUE, ...) {
  fit <- as_public_gRm_fit(fit)
  values <- if (is_active_public_fit(fit)) {
    item_fits_values(fit, include_extended = include_extended)
  } else {
    item_fits_values(fit$project %||% fit$analysis$project, include_extended = include_extended)
  }
  new_gRm_result(
    class = "gRm_item_fit",
    analysis = fit$analysis %||% fit$spec$analysis,
    fit = fit,
    values = values,
    result = "item_fit",
    metadata = list(include_extended = isTRUE(include_extended)),
    call = match.call()
  )
}

#' Local-dependence diagnostics
#'
#' @param fit A fitted DIGRAM model.
#' @param max_step Maximum number of fitting iterations for candidate models.
#' @param max_delta Convergence threshold for candidate models.
#' @param jobs Number of parallel jobs. Defaults to source-stable serial work.
#' @param ... Reserved for S3 dispatch compatibility; ignored.
#' @return A `gRm_local_dependence` result object.
#' @export
#' @examples
#' \donttest{
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0)
#' )
#' fit0 <- fit(gllrm(gRm(data, items = c("I1", "I2"), id = "ID")))
#' ld <- local_dependence(fit0)
#' summary(ld, which = "tests")
#' }
local_dependence <- function(fit,
                             max_step = 5000L,
                             max_delta = 0.0001,
                             jobs = 1L,
                             ...) {
  fit <- as_public_gRm_fit(fit)
  jobs <- normalize_public_jobs(jobs)
  values <- if (is_active_public_fit(fit)) {
    active_gllrm_local_independence_values(
      fit,
      max_step = max_step,
      max_delta = max_delta,
      jobs = jobs
    )
  } else {
    local_independence_values(
      fit$project %||% fit$analysis$project,
      max_step = max_step,
      max_delta = max_delta,
      jobs = jobs
    )
  }
  new_gRm_result(
    class = "gRm_local_dependence",
    analysis = fit$analysis %||% fit$spec$analysis,
    fit = fit,
    values = values,
    result = "local_dependence",
    metadata = list(max_step = max_step, max_delta = max_delta, jobs = jobs),
    call = match.call()
  )
}

#' Differential item functioning diagnostics
#'
#' @param fit A fitted DIGRAM model.
#' @param max_step Maximum number of fitting iterations for candidate models.
#' @param max_delta Convergence threshold for candidate models.
#' @param jobs Number of parallel jobs. Defaults to source-stable serial work.
#' @param ... Reserved for S3 dispatch compatibility; ignored.
#' @return A `gRm_dif` result object.
#' @export
#' @examples
#' \donttest{
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0),
#'   site = c(0, 0, 1, 1, 0, 1, 0, 1)
#' )
#' analysis <- gRm(data, items = c("I1", "I2"), exogenous = "site", id = "ID")
#' fit0 <- fit(gllrm(analysis))
#' dif_tests <- dif(fit0)
#' summary(dif_tests, which = "tests")
#' }
dif <- function(fit,
                max_step = 5000L,
                max_delta = 0.0001,
                jobs = 1L,
                ...) {
  fit <- as_public_gRm_fit(fit)
  jobs <- normalize_public_jobs(jobs)
  values <- if (is_active_public_fit(fit)) {
    active_gllrm_dif_tests_values(
      fit,
      max_step = max_step,
      max_delta = max_delta,
      jobs = jobs
    )
  } else {
    dif_tests_values(
      fit$project %||% fit$analysis$project,
      max_step = max_step,
      max_delta = max_delta,
      jobs = jobs
    )
  }
  new_gRm_result(
    class = "gRm_dif",
    analysis = fit$analysis %||% fit$spec$analysis,
    fit = fit,
    values = values,
    result = "dif",
    metadata = list(max_step = max_step, max_delta = max_delta, jobs = jobs),
    call = match.call()
  )
}

#' Global homogeneity diagnostics
#'
#' @param fit A fitted DIGRAM model.
#' @param score_cuts Optional integer-like upper total-score cut values. When
#'   `NULL`, the analysis-level score cuts stored by [gRm()] or
#'   [read_digram_project()] are used. When supplied, these cuts override the
#'   analysis-level cuts for this global-homogeneity calculation only. The cuts
#'   define score groups as consecutive total-score intervals: the first group
#'   runs from the source-valid lowest score through the first cut, the next
#'   group starts at the following score and runs through the next cut, and so
#'   on.
#' @param max_step Maximum number of fitting iterations for group models.
#' @param max_delta Convergence threshold for group models.
#' @param jobs Reserved for future parallel implementations.
#' @param ... Reserved for S3 dispatch compatibility; ignored.
#' @return A `gRm_global_homogeneity` result object.
#' @export
#' @examples
#' \donttest{
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0)
#' )
#' fit0 <- fit(gllrm(gRm(data, items = c("I1", "I2"), id = "ID", score_cuts = c(1L, 2L))))
#' gh <- global_homogeneity(fit0)
#' summary(gh, which = "tests")
#' }
global_homogeneity <- function(fit,
                               score_cuts = NULL,
                               max_step = 5000L,
                               max_delta = 0.0001,
                               jobs = 1L,
                               ...) {
  fit <- as_public_gRm_fit(fit)
  analysis <- fit$analysis %||% fit$spec$analysis
  jobs <- normalize_public_jobs(jobs)
  score_cuts <- normalize_public_score_cuts(score_cuts %||% analysis$score_groups, analysis$project)
  values <- if (is_active_public_fit(fit)) {
    global_homogeneity_values(
      fit,
      score_cuts = score_cuts,
      max_step = max_step,
      max_delta = max_delta
    )
  } else {
    global_homogeneity_values(
      fit$project %||% analysis$project,
      score_cuts = score_cuts,
      max_step = max_step,
      max_delta = max_delta,
      bundle = fit$bundle %||% NULL,
      base_fit = fit$fit %||% NULL
    )
  }
  new_gRm_result(
    class = "gRm_global_homogeneity",
    analysis = analysis,
    fit = fit,
    values = values,
    result = "global_homogeneity",
    metadata = list(score_cuts = score_cuts, max_step = max_step, max_delta = max_delta, jobs = jobs),
    call = match.call()
  )
}

new_gRm_result <- function(class,
                              analysis,
                              fit,
                              values,
                              result,
                              metadata = list(),
                              call = NULL) {
  out <- list(
    analysis = analysis,
    project = analysis$project,
    fit = fit,
    values = values,
    result = result,
    metadata = metadata,
    source_trace = c(
      analysis$source_trace %||% character(),
      values$source_status %||% character(),
      api = result
    ),
    warnings = character(),
    unmodeled = character(),
    call = call
  )
  class(out) <- c(class, "list")
  out
}

as_public_gRm_analysis <- function(x) {
  if (inherits(x, "gRm_analysis") || inherits(x, "gRm_item_analysis")) {
    return(x)
  }
  if (inherits(x, "gRm_model") || inherits(x, "gRm_gllrm_spec")) {
    return(x$analysis)
  }
  if (inherits(x, "gRm_fit") || inherits(x, "gRm_gllrm_fit")) {
    return(x$analysis %||% x$spec$analysis)
  }
  if (inherits(x, "gRm_screen")) {
    return(x$analysis)
  }
  stop("Expected a DIGRAM analysis object.", call. = FALSE)
}

as_public_gRm_fit <- function(x) {
  if (inherits(x, "gRm_fit") || inherits(x, "gRm_gllrm_fit")) {
    return(x)
  }
  stop("Expected a fitted DIGRAM model.", call. = FALSE)
}

is_active_public_fit <- function(fit) {
  inherits(fit$values, "gRm_active_gllrm_values") ||
    nrow(fit$spec$ld %||% data.frame()) > 0L ||
    nrow(fit$spec$dif %||% data.frame()) > 0L
}

normalize_public_jobs <- function(jobs) {
  if (is.null(jobs)) {
    jobs <- 1L
  }
  if (length(jobs) != 1L || is.na(jobs) || jobs != as.integer(jobs) || jobs < 1L) {
    stop("`jobs` must be a positive integer.", call. = FALSE)
  }
  as.integer(jobs)
}

normalize_public_score_cuts <- function(score_cuts, project) {
  cuts <- as.integer(score_cuts %||% integer())
  if (length(cuts) < 2L) {
    cuts <- gRm_default_global_homogeneity_score_cuts(project)
  }
  if (length(cuts) < 2L || anyNA(cuts) || any(cuts != as.integer(cuts))) {
    stop("`score_cuts` must contain at least two non-missing integer-like score cuts.", call. = FALSE)
  }
  if (is.unsorted(cuts, strictly = TRUE)) {
    stop("`score_cuts` must be strictly increasing.", call. = FALSE)
  }
  cuts
}
