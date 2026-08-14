#' Exogenous score-effect diagnostics
#'
#' @param analysis A DIGRAM analysis object.
#' @param inference Inference mode. `"asymptotic"` maps to DIGRAM command 17,
#'   `"exact"` to command 2, and `"repeated"` to command 74.
#' @param nsim Single non-negative integer-like number of simulations for exact
#'   branches. The source convention `0` requests the DIGRAM default of `1000`.
#' @param seed Single integer-like random seed for exact branches.
#' @param critlevel Optional DIGRAM exact command critical level on the
#'   per-1000 scale; when supplied it must be a single non-negative
#'   integer-like value.
#' @param risk Optional DIGRAM repeated-Monte-Carlo risk on the per-1000 scale;
#'   when supplied it must be a single non-negative integer-like value.
#' @return A data-frame-like `gRm_score_effects` table. The unchanged
#'   `gRm_exo_select_values` backend object is available as
#'   `attr(result, "values")`; the selected rows and Benjamini-Hochberg
#'   threshold table are available as `attr(result, "selected")` and
#'   `attr(result, "bh")`.
#' @details
#' `score_effects()` returns its public test table directly. The source-shaped
#' backend values remain available as attributes rather than through separate
#' `summary()` views. Use `attr(result, "values")` for the full backend object,
#' `attr(result, "selected")` for source-selected rows, and `attr(result, "bh")`
#' for the Benjamini-Hochberg thresholds.
#'
#' The public table always includes `Exogenous`, `Hypothesis`, `Chisq`, `Df`,
#' `Pr(>Chisq)`, `Gamma`, `Pr(Gamma+)`, and `Pr(|Gamma|)`. For `inference =
#' "exact"` and `inference = "repeated"`, the table also includes exact or
#' simulation p-value columns and `Simulations` when those values are present.
#' Internal source-selection fields such as labels, markers, and selected flags
#' remain in `attr(result, "values")` or `attr(result, "selected")`.
#' @export
#' @examples
#' \donttest{
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0),
#'   site = c(0, 0, 1, 1, 0, 1, 0, 1)
#' )
#' analysis <- gRm(
#'   data,
#'   items = c("I1", "I2"),
#'   exogenous = "site",
#'   id = "ID",
#'   score_cuts = c(1L, 2L)
#' )
#' effects <- score_effects(analysis)
#' effects
#' }
score_effects <- function(analysis,
                          inference = c("asymptotic", "exact", "repeated"),
                          nsim = 1000L,
                          seed = 9L,
                          critlevel = NULL,
                          risk = NULL) {
  analysis <- as_public_gRm_analysis(analysis)
  inference <- match.arg(inference)
  exact_state <- gRm_exact_command_state_public(
    inference,
    nsim = nsim,
    seed = seed,
    critlevel = critlevel,
    risk = risk
  )
  values <- exo_select_values(
    analysis$project,
    score_cap = gRm_source_score_cap(),
    exact = exact_state$exact,
    repeated = exact_state$sequential,
    nsim = exact_state$nsim,
    seed = exact_state$seed,
    exact_state = exact_state
  )
  score_effects_public_table(
    values,
    metadata = list(
      inference = inference,
      nsim = exact_state$nsim,
      seed = if (exact_state$exact) exact_state$seed else NA_integer_,
      critlevel = if (is.null(critlevel)) NA_integer_ else as.integer(critlevel),
      risk = if (is.null(risk)) NA_integer_ else as.integer(risk),
      score_cap = gRm_source_score_cap(),
      exact_state = exact_state
    ),
    call = match.call()
  )
}

score_effects_public_table <- function(values, metadata, call) {
  bh <- score_effects_bh_table(values)
  table <- score_effects_tests_table(values, metadata$inference %||% "asymptotic")
  selected <- normalize_summary_table(values$selected %||% data.frame())
  out <- make_gRm_direct_table(
    table,
    class = "gRm_score_effects",
    values = values,
    bh = bh,
    which = "tests",
    title = "gRm: Score-effect tests",
    table_note = score_effects_bh_footer(bh),
    result = "score_effects",
    metadata = metadata,
    call = call
  )
  attr(out, "selected") <- selected
  out
}

score_effects_tests_table <- function(values, inference) {
  tests <- normalize_summary_table(values$screen %||% data.frame())
  if (!nrow(tests)) {
    return(data.frame(
      Exogenous = character(),
      Hypothesis = character(),
      Chisq = numeric(),
      Df = integer(),
      `Pr(>Chisq)` = numeric(),
      Gamma = numeric(),
      `Pr(Gamma+)` = numeric(),
      `Pr(|Gamma|)` = numeric(),
      check.names = FALSE
    ))
  }
  out <- data.frame(
    Exogenous = score_effects_column(tests, "exo_name", NA_character_),
    Hypothesis = score_effects_column(tests, "hypothesis", NA_character_),
    Chisq = score_effects_column(tests, "chi_square", NA_real_),
    Df = score_effects_column(tests, "df", NA_integer_),
    `Pr(>Chisq)` = score_effects_column(tests, "chi_p", NA_real_),
    Gamma = score_effects_column(tests, "gamma", NA_real_),
    `Pr(Gamma+)` = score_effects_column(tests, "gamma_p_one_sided", NA_real_),
    `Pr(|Gamma|)` = score_effects_column(tests, "gamma_p_two_sided", NA_real_),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  exact_columns <- c("exact_chi_p", "exact_gamma_p_one_sided", "exact_gamma_p_two_sided", "exact_nsim")
  has_exact_values <- any(vapply(
    exact_columns,
    function(column) column %in% names(tests) && any(!is.na(tests[[column]])),
    logical(1L)
  ))
  if (inference %in% c("exact", "repeated") && has_exact_values) {
    out[["Exact Pr(>Chisq)"]] <- score_effects_column(tests, "exact_chi_p", NA_real_)
    out[["Exact Pr(Gamma+)"]] <- score_effects_column(tests, "exact_gamma_p_one_sided", NA_real_)
    out[["Exact Pr(|Gamma|)"]] <- score_effects_column(tests, "exact_gamma_p_two_sided", NA_real_)
    out[["Simulations"]] <- score_effects_column(tests, "exact_nsim", NA_integer_)
  }
  out
}

score_effects_column <- function(tests, column, default) {
  if (column %in% names(tests)) {
    return(tests[[column]])
  }
  rep(default, nrow(tests))
}

#' Item fit diagnostics
#'
#' Run post-fit item diagnostics for a fitted gRm model. `item_fit()` is a
#' diagnostic procedure, not an item-parameter summary accessor: it checks the
#' fitted model against observed item behavior using outfit, infit, and
#' item-restscore gamma diagnostics. Use `which = "tests"` for compact outfit,
#' infit, and item-restscore gamma tests. Use `which = "items"` for the
#' per-item outfit/infit diagnostic summary produced by the extended item fit
#' calculation.
#'
#' @param fit A fitted gRm model.
#' @param which Item fit table to return. `"tests"` returns the compact
#'   inferential item fit diagnostics and is the default. `"items"` returns the
#'   extended per-item outfit/infit diagnostic decomposition.
#' @param include_extended Whether to compute extended item fit detail tables.
#'   When `FALSE`, `item_fit(fit, which = "items")` returns an empty data
#'   frame with the expected item-summary columns because the per-item extended
#'   diagnostic summaries were not computed. The compact `which = "tests"`
#'   table remains available.
#' @param ... Reserved for S3 dispatch compatibility; ignored.
#' @return A data-frame-like `gRm_item_fit` table. The unchanged
#'   `gRm_item_fits_values` backend object is available as
#'   `attr(result, "values")`; the Benjamini-Hochberg threshold table is
#'   available as `attr(result, "bh")`.
#' @details
#' The item fit result has two different public diagnostic tables. The
#' `which = "tests"` table is always available and returns one row per item
#' with the compact inferential statistics: outfit, infit, item-restscore
#' gamma, standard errors, p-values, diagnostic-specific FDR markers, and
#' gamma direction. Its public columns are `Item`, `Outfit`, `Outfit SE`,
#' `Pr(>Outfit)`, `Outfit FDR`, `Infit`, `Infit SE`, `Pr(>Infit)`,
#' `Infit FDR`, `Observed gamma`, `Expected gamma`, `Gamma SE`,
#' `Pr(>Gamma)`, `Gamma FDR`, and `Gamma direction`. The marker columns use
#' `""`, `"*"`, `"**"`, and `"***"` for no selection and selection at the
#' 5 percent, 1 percent, and 0.1 percent FDR levels, respectively. DIGRAM
#' computes these markers separately across items for outfit, infit, and
#' gamma; they are not obtained from the combined thresholds printed below
#' the table. When the result is printed, each marker is appended in a fixed
#' three-character field directly after its corresponding p-value, as in
#' DIGRAM. This preserves numeric alignment while the returned p-value columns
#' remain numeric. The source-facing backend table in
#' `attr(result, "values")$items` keeps the original DIGRAM-shaped columns,
#' including `item_label`, `item_name`, `outfit_fdr`, `infit_fdr`, and
#' `gamma_fdr`.
#'
#' The `which = "items"` table is also one row per item, but it is not a
#' second test table or an item-parameter table. It returns the extended
#' outfit/infit diagnostic decomposition: aggregate observed and expected
#' outfit components, the resulting outfit value, and the observed, expected,
#' variance, and ratio components for infit. Its public columns are `Item`,
#' `Outfit N`, `Outfit observed`, `Outfit expected`, `Outfit total`,
#' `Infit observed`, `Infit expected`, `Infit variance`, and `Infit ratio`.
#' The source-facing backend table in
#' `attr(result, "values")$extended$summaries` keeps the original DIGRAM-shaped
#' columns, including `item_label`, `item_name`, `outfit_total_value`, and
#' `infit_value`.
#'
#' The Benjamini-Hochberg critical p-value limits used for the combined compact
#' item fit test family are printed below the tests table and are also metadata
#' rather than a separate printed view. Use `attr(result, "bh")` to inspect all
#' available thresholds programmatically.
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
#' `item_fit(fit, which = "items")` is populated. With `include_extended =
#' FALSE`, the extended calculation is skipped and `item_fit(fit, which =
#' "items")` returns an empty data frame. The compact `which = "tests"` table
#' and BH metadata remain available.
#'
#' @section Incomplete item fit records:
#' The original DIGRAM source has additional item fit and item-restscore gamma
#' branches for incomplete response records when the runtime
#' `NincompleteRecs` state has been populated. The current gRm implementation
#' does not have a source-backed way to determine when that runtime list exists,
#' so public `item_fit()` results do not synthesize incomplete item fit records
#' from rows with missing item responses. Those rows are therefore outside the
#' source-faithful item fit scope of this package version unless a future
#' source/runtime trace establishes the corresponding DIGRAM state. Public
#' item fit values expose `incomplete_records_used = FALSE` and
#' `incomplete_records_status = "not_source_backed"` until the full source
#' include-incomplete path is implemented.
#' @export
#' @examples
#' \donttest{
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0)
#' )
#' fit0 <- fit(gllrm(gRm(
#'   data,
#'   items = c("I1", "I2"),
#'   id = "ID",
#'   score_cuts = c(1L, 2L)
#' )))
#' item_tests <- item_fit(fit0)
#' item_summaries <- item_fit(fit0, which = "items")
#' }
item_fit <- function(fit, which = c("tests", "items"), include_extended = TRUE, ...) {
  fit <- as_public_gRm_fit(fit)
  if (missing(which)) {
    which <- "tests"
  }
  if (!is.character(which) || length(which) != 1L || !which %in% c("tests", "items")) {
    stop("`which` must be one of \"tests\" or \"items\".", call. = FALSE)
  }
  values <- if (is_gllrm_public_fit(fit)) {
    item_fits_values(fit, include_extended = include_extended)
  } else {
    base_item_fits_values(fit$bundle, fit$fit, include_extended = include_extended)
  }
  item_fit_public_table(
    values,
    which = which,
    include_extended = isTRUE(include_extended),
    call = match.call()
  )
}

item_fit_public_table <- function(values,
                                  which,
                                  include_extended,
                                  call) {
  bh <- item_fit_bh_table(values)
  table <- switch(
    which,
    tests = item_fit_tests_table(values),
    items = item_fit_items_table(values)
  )
  note <- if (identical(which, "tests")) {
    paste(c(
      paste0(
        "Stars are based on separate Benjamini-Hochberg adjustments across ",
        "items for outfit, infit, and gamma: * = selected at 5% FDR, ",
        "** = selected at 1% FDR, *** = selected at 0.1% FDR. ",
        "The thresholds below are based on all item-fit p-values pooled ",
        "together and are reported separately."
      ),
      item_fit_bh_footer(bh)
    ), collapse = " ")
  } else {
    character()
  }
  title <- switch(
    which,
    tests = "gRm: Item fit tests",
    items = "gRm: Item fit item diagnostics"
  )
  make_gRm_direct_table(
    table,
    class = "gRm_item_fit",
    values = values,
    bh = bh,
    which = which,
    title = title,
    table_note = note,
    result = "item_fit",
    metadata = list(include_extended = include_extended),
    call = call
  )
}

item_fit_tests_table <- function(values) {
  table <- values$items %||% data.frame()
  if (!is.data.frame(table)) {
    table <- data.frame()
  }
  n <- nrow(table)
  data.frame(
    Item = item_fit_public_column(table, "item_name", character(), n),
    Outfit = item_fit_public_column(table, "outfit", numeric(), n),
    `Outfit SE` = item_fit_public_column(table, "outfit_sd", numeric(), n),
    `Pr(>Outfit)` = item_fit_public_column(table, "p_outfit", numeric(), n),
    `Outfit FDR` = item_fit_fdr_marker(
      item_fit_public_column(table, "outfit_fdr", integer(), n)
    ),
    Infit = item_fit_public_column(table, "infit", numeric(), n),
    `Infit SE` = item_fit_public_column(table, "infit_sd", numeric(), n),
    `Pr(>Infit)` = item_fit_public_column(table, "p_infit", numeric(), n),
    `Infit FDR` = item_fit_fdr_marker(
      item_fit_public_column(table, "infit_fdr", integer(), n)
    ),
    `Observed gamma` = item_fit_public_column(table, "observed_gamma", numeric(), n),
    `Expected gamma` = item_fit_public_column(table, "expected_gamma", numeric(), n),
    `Gamma SE` = item_fit_public_column(table, "gamma_sd", numeric(), n),
    `Pr(>Gamma)` = item_fit_public_column(table, "p_gamma", numeric(), n),
    `Gamma FDR` = item_fit_fdr_marker(
      item_fit_public_column(table, "gamma_fdr", integer(), n)
    ),
    `Gamma direction` = item_fit_public_column(table, "direction", character(), n),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Convert DIGRAM item-fit FDR risk grades to display markers
#'
#' @param risk Integer vector using the source grades 0, 1, 2, and 3.
#' @return A character vector containing `""`, `"*"`, `"**"`, or `"***"`.
#' @keywords internal
#' @noRd
item_fit_fdr_marker <- function(risk) {
  markers <- c("", "*", "**", "***")
  unname(markers[match(as.integer(risk), 0:3)])
}

item_fit_items_table <- function(values) {
  table <- values$extended$summaries %||% item_fit_empty_items_table()
  if (!is.data.frame(table) || !nrow(table)) {
    return(item_fit_empty_items_table())
  }
  n <- nrow(table)
  data.frame(
    Item = item_fit_public_column(table, "item_name", character(), n),
    `Outfit N` = item_fit_public_column(table, "outfit_total_n", integer(), n),
    `Outfit observed` = item_fit_public_column(table, "outfit_total_observed", numeric(), n),
    `Outfit expected` = item_fit_public_column(table, "outfit_total_expected", numeric(), n),
    `Outfit total` = item_fit_public_column(table, "outfit_total_value", numeric(), n),
    `Infit observed` = item_fit_public_column(table, "infit_observed", numeric(), n),
    `Infit expected` = item_fit_public_column(table, "infit_expected", numeric(), n),
    `Infit variance` = item_fit_public_column(table, "infit_variance", numeric(), n),
    `Infit ratio` = item_fit_public_column(table, "infit_value", numeric(), n),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

item_fit_empty_items_table <- function() {
  data.frame(
    Item = character(),
    `Outfit N` = integer(),
    `Outfit observed` = numeric(),
    `Outfit expected` = numeric(),
    `Outfit total` = numeric(),
    `Infit observed` = numeric(),
    `Infit expected` = numeric(),
    `Infit variance` = numeric(),
    `Infit ratio` = numeric(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

item_fit_public_column <- function(table, column, prototype, n) {
  if (column %in% names(table)) {
    return(table[[column]])
  }
  rep(item_fit_public_missing_value(prototype), n)
}

item_fit_public_missing_value <- function(prototype) {
  if (is.character(prototype)) {
    return(NA_character_)
  }
  if (is.integer(prototype)) {
    return(NA_integer_)
  }
  if (is.numeric(prototype)) {
    return(NA_real_)
  }
  if (is.logical(prototype)) {
    return(NA)
  }
  NA
}

item_fit_bh_table <- function(values) {
  table <- result_named_numeric_table(values$bh_limits %||% numeric(), "threshold", "p_value")
  if (!nrow(table)) {
    return(table)
  }
  table$fdr <- item_fit_bh_fdr_label(table$threshold)
  table[c("threshold", "fdr", "p_value")]
}

item_fit_bh_footer <- function(bh, digits = max(3L, getOption("digits") - 3L)) {
  if (!is.data.frame(bh) || !nrow(bh)) {
    return(character())
  }
  footer_rows <- bh[bh$fdr %in% c("0.05", "0.01"), , drop = FALSE]
  footer_rows <- footer_rows[match(c("0.05", "0.01"), footer_rows$fdr, nomatch = 0L), , drop = FALSE]
  if (!nrow(footer_rows)) {
    footer_rows <- bh
  }
  dig_tst <- max(1L, min(5L, digits - 1L))
  parts <- paste0(
    "FDR ",
    footer_rows$fdr,
    " = ",
    format.pval(footer_rows$p_value, digits = dig_tst, eps = .Machine$double.eps)
  )
  paste0("Benjamini-Hochberg thresholds: ", paste(parts, collapse = ", "))
}

score_effects_bh_table <- function(values) {
  bh <- values$bh %||% list()
  data.frame(
    threshold = c("fdr_05", "fdr_01"),
    fdr = c("0.05", "0.01"),
    p_value = c(bh$fdr_05 %||% NA_real_, bh$fdr_01 %||% NA_real_),
    stringsAsFactors = FALSE
  )
}

score_effects_bh_footer <- function(bh, digits = max(3L, getOption("digits") - 3L)) {
  item_fit_bh_footer(bh, digits = digits)
}

item_fit_bh_fdr_label <- function(threshold) {
  out <- as.character(threshold)
  out[out %in% c("fdr_5", "fdr_05")] <- "0.05"
  out[out %in% c("fdr_1", "fdr_01")] <- "0.01"
  out[out %in% c("fdr_001")] <- "0.001"
  out
}

make_gRm_direct_table <- function(table,
                                  class,
                                  values,
                                  bh,
                                  which,
                                  title = NULL,
                                  table_note = character(),
                                  analysis = NULL,
                                  fit = NULL,
                                  result = NULL,
                                  metadata = list(),
                                  call = NULL) {
  if (!is.data.frame(table)) {
    table <- data.frame()
  }
  rownames(table) <- NULL
  class(table) <- c(class, "gRm_direct_table", "data.frame")
  attr(table, "values") <- values
  attr(table, "bh") <- bh
  attr(table, "which") <- which
  attr(table, "title") <- title
  attr(table, "table_note") <- table_note
  attr(table, "analysis") <- analysis
  attr(table, "fit") <- fit
  attr(table, "result") <- result
  attr(table, "metadata") <- metadata
  attr(table, "call") <- call
  table
}

#' Local-dependence diagnostics
#'
#' @param fit A fitted DIGRAM model.
#' @param max_step Single positive integer-like maximum number of fitting
#'   iterations for candidate models.
#' @param max_delta Single positive finite convergence threshold for candidate
#'   models.
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
#' fit0 <- fit(gllrm(gRm(
#'   data,
#'   items = c("I1", "I2"),
#'   id = "ID",
#'   score_cuts = c(1L, 2L)
#' )))
#' ld <- local_dependence(fit0)
#' summary(ld)
#' }
local_dependence <- function(fit,
                             max_step = 5000L,
                             max_delta = 0.0001,
                             jobs = 1L,
                             ...) {
  fit <- as_public_gRm_fit(fit)
  controls <- normalize_public_fit_controls(max_step, max_delta)
  max_step <- controls$max_step
  max_delta <- controls$max_delta
  jobs <- normalize_public_jobs(jobs)
  values <- if (is_gllrm_public_fit(fit)) {
    gllrm_local_independence_values(
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
#' @param max_step Single positive integer-like maximum number of fitting
#'   iterations for candidate models.
#' @param max_delta Single positive finite convergence threshold for candidate
#'   models.
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
#' analysis <- gRm(
#'   data,
#'   items = c("I1", "I2"),
#'   exogenous = "site",
#'   id = "ID",
#'   score_cuts = c(1L, 2L)
#' )
#' fit0 <- fit(gllrm(analysis))
#' dif_tests <- dif(fit0)
#' summary(dif_tests)
#' }
dif <- function(fit,
                max_step = 5000L,
                max_delta = 0.0001,
                jobs = 1L,
                ...) {
  fit <- as_public_gRm_fit(fit)
  controls <- normalize_public_fit_controls(max_step, max_delta)
  max_step <- controls$max_step
  max_delta <- controls$max_delta
  jobs <- normalize_public_jobs(jobs)
  values <- if (is_gllrm_public_fit(fit)) {
    gllrm_dif_tests_values(
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
#'   on. Supplied cuts must contain at least two non-missing integer-like values,
#'   be strictly increasing, lie inside the possible score range, and leave at
#'   least two usable source score groups after boundary-score handling.
#' @param max_step Single positive integer-like maximum number of fitting
#'   iterations for group models.
#' @param max_delta Single positive finite convergence threshold for group
#'   models.
#' @param ... Reserved for S3 dispatch compatibility; ignored.
#' @return A `gRm_global_homogeneity` result object.
#' @details
#' `summary(result, which = "test")` reports the source-backed likelihood-ratio
#' test, degrees of freedom, and p-value. `summary(result, which =
#' "score_groups")` reports the fitted source score groups and convergence
#' status. `summary(result, which = "item_means")` reports the score-group item
#' means that gRm can derive source-faithfully from the available Pascal path. If
#' the fitted model contains LD and/or DIF terms, `summary(result)` also prints
#' the source-backed extended global-homogeneity sections for uniform local
#' dependence and uniform DIF. These sections are available directly with
#' `which = "uniform_ld"` and `which = "uniform_dif"`; their observed and
#' expected gamma columns are named by the score-group labels defined for the
#' particular analysis, followed by chi-square, degrees of freedom, and p-value
#' columns. Residual and marker cells that are not source-backed remain
#' available on the raw result values but are not printed in the public summary
#' tables.
#'
#' @section Global homogeneity residuals:
#' The item-level residual and marker cells in DIGRAM's global-homogeneity
#' report are a documented implementation restriction. The recovered Pascal
#' source identifies the residual/item fit report path, but not the hidden
#' runtime residual variance materialization needed to reproduce those cells
#' source-faithfully. gRm therefore leaves the item-level `residual` and
#' `marker` fields as `NA` and reports their status as `"not_source_backed"`
#' rather than inventing an unsupported formula.
#' @export
#' @examples
#' \donttest{
#' data <- expand.grid(I1 = 0:1, I2 = 0:1, I3 = 0:1, I4 = 0:1)
#' data$ID <- seq_len(nrow(data))
#' fit0 <- fit(gllrm(gRm(
#'   data,
#'   items = c("I1", "I2", "I3", "I4"),
#'   id = "ID",
#'   score_cuts = c(1L, 3L)
#' )))
#' gh <- global_homogeneity(fit0)
#' summary(gh, which = "test")
#' }
global_homogeneity <- function(fit,
                               score_cuts = NULL,
                               max_step = 5000L,
                               max_delta = 0.0001,
                               ...) {
  reject_public_dots(...)
  fit <- as_public_gRm_fit(fit)
  controls <- normalize_public_fit_controls(max_step, max_delta)
  max_step <- controls$max_step
  max_delta <- controls$max_delta
  analysis <- fit$analysis %||% fit$spec$analysis
  score_cuts <- normalize_public_score_cuts(
    score_cuts,
    analysis$project,
    default = analysis$score_groups,
    bundle = fit$bundle %||% NULL
  )
  values <- if (is_gllrm_public_fit(fit)) {
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
    metadata = list(score_cuts = score_cuts, max_step = max_step, max_delta = max_delta),
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
  if (inherits(x, "gRm_analysis")) {
    return(x)
  }
  if (inherits(x, "gRm_model")) {
    return(x$analysis)
  }
  if (inherits(x, "gRm_fit")) {
    return(x$analysis %||% x$spec$analysis)
  }
  if (inherits(x, "gRm_screen")) {
    return(x$analysis)
  }
  stop("Expected a DIGRAM analysis object.", call. = FALSE)
}

as_public_gRm_fit <- function(x) {
  if (inherits(x, "gRm_fit")) {
    return(x)
  }
  stop("Expected a fitted DIGRAM model.", call. = FALSE)
}

is_gllrm_public_fit <- function(fit) {
  inherits(fit$values, "gRm_gllrm_values") ||
    nrow(fit$spec$ld %||% data.frame()) > 0L ||
    nrow(fit$spec$dif %||% data.frame()) > 0L
}

normalize_public_jobs <- function(jobs) {
  if (is.null(jobs)) {
    jobs <- 1L
  }
  normalize_public_integer_like(
    jobs,
    "`jobs` must be a positive integer.",
    scalar = TRUE,
    lower = 1L
  )
}

reject_public_dots <- function(...) {
  dots <- list(...)
  if (length(dots)) {
    stop("The ... argument is reserved for future extensions and must be empty.", call. = FALSE)
  }
  invisible(NULL)
}

normalize_public_score_cuts <- function(score_cuts, project, default = NULL, bundle = NULL) {
  if (is.null(score_cuts)) {
    score_cuts <- default %||% integer()
    if (length(score_cuts) < 2L) {
      score_cuts <- gRm_default_global_homogeneity_score_cuts(project)
    }
  } else if (!is.numeric(score_cuts) && !is.integer(score_cuts)) {
    stop("`score_cuts` must be an integer-like vector of score cuts.", call. = FALSE)
  }

  cuts <- normalize_public_integer_like(
    score_cuts,
    "`score_cuts` must contain at least two non-missing integer-like score cuts.",
    min_length = 2L
  )
  if (is.unsorted(cuts, strictly = TRUE)) {
    stop("`score_cuts` must be strictly increasing.", call. = FALSE)
  }

  max_score <- sum(project$items$raw_max - 1L)
  if (any(cuts < 0L | cuts > max_score)) {
    stop("`score_cuts` must lie within the possible score range 0..", max_score, ".", call. = FALSE)
  }

  if (!is.null(bundle)) {
    groups <- global_homogeneity_score_groups(bundle, cuts)
    if (nrow(groups) < 2L) {
      stop("`score_cuts` must define at least two global-homogeneity score groups.", call. = FALSE)
    }
  }

  cuts
}
