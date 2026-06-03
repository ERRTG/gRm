#' Tidy a gRm object
#'
#' `tidy()` returns a compact primary table for inspection, filtering,
#' joining, or plotting. Use details() internally for complete structured numeric
#' result tables.
#'
#' @param x A clean gRm result object.
#' @param ... Reserved for S3 dispatch compatibility. Concrete method
#'   arguments are documented on the corresponding public result objects.
#' @return A data frame.
#' @noRd
tidy <- function(x, ...) {
  UseMethod("tidy")
}

#' Glance at a gRm object
#'
#' `glance()` returns exactly one row with object, run, model, or validation
#' status. It does not contain full detail tables.
#'
#' @param x A clean gRm result object.
#' @param ... Reserved for S3 dispatch compatibility; ignored.
#' @return A one-row data frame.
#' @noRd
glance <- function(x, ...) {
  UseMethod("glance")
}

#' Access detailed result tables
#'
#' `details()` returns source-faithful structured numeric tables already
#' calculated by the package. It does not render or parse DIGRAM text output.
#'
#' @param x A clean gRm object.
#' @param name Optional name of one detailed result table.
#' @param ... Reserved for S3 dispatch compatibility; ignored.
#' @return `details()` returns a `gRm_details` object when `name` is
#'   `NULL`, or one data frame when `name` is supplied. `detail_names()`
#'   returns a character vector of available detail table names.
#' @noRd
details <- function(x, name = NULL, ...) {
  UseMethod("details")
}

detail_names <- function(x, ...) {
  names(details(x, ...)$tables)
}

tidy.gRm_item_analysis <- function(x, ...) {
  item_analysis_item_summary(x)
}

glance.gRm_item_analysis <- function(x, ...) {
  item_analysis_glance_row(x)
}

tidy.gRm_gllrm_fit <- function(x,
                                  type = c("parameters", "items", "thresholds", "ld", "dif", "model_terms"),
                                  ...) {
  type <- match.arg(type)
  detail <- details(x)

  if (type == "items") {
    return(transform_for_tidy(x$values$item_statistics %||% data.frame()))
  }
  if (type == "thresholds") {
    return(detail$tables$thresholds %||% data.frame())
  }
  if (type == "ld") {
    return(detail$tables$ld_parameters %||% data.frame())
  }
  if (type == "dif") {
    return(detail$tables$dif_parameters %||% data.frame())
  }
  if (type == "model_terms") {
    return(detail$tables$model_terms %||% data.frame())
  }

  item_gamma <- detail$tables$item_gamma %||% data.frame()
  if (nrow(item_gamma)) {
    item_gamma$term_type <- "item"
    item_gamma$parameter_type <- "gamma"
    item_gamma$estimate <- item_gamma$gamma
    item_gamma$scale <- "multiplicative"
  }
  thresholds <- detail$tables$thresholds %||% data.frame()
  if (nrow(thresholds)) {
    thresholds$term_type <- "item"
    thresholds$parameter_type <- "threshold"
    thresholds$estimate <- thresholds$threshold
    thresholds$scale <- "threshold"
  }
  ld <- detail$tables$ld_parameters %||% data.frame()
  if (nrow(ld)) {
    ld$term_type <- "ld"
    ld$parameter_type <- "ld_gamma"
    ld$estimate <- ld$gamma
    ld$scale <- "multiplicative"
  }
  dif <- detail$tables$dif_parameters %||% data.frame()
  if (nrow(dif)) {
    dif$term_type <- "dif"
    dif$parameter_type <- "dif_gamma"
    dif$estimate <- dif$gamma
    dif$scale <- "multiplicative"
  }

  rbind_fill(item_gamma, thresholds, ld, dif)
}

glance.gRm_gllrm_fit <- function(x, ...) {
  values <- x$values %||% list()
  fit_summary <- item_parameter_fit_summary(values)
  data.frame(
    n_persons = nrow(x$analysis$project$raw_data %||% data.frame()),
    n_valid = values$likelihood_n %||% NA_integer_,
    n_items = length(x$analysis$items),
    n_exogenous = length(x$analysis$exogenous),
    n_ld = nrow(x$spec$ld),
    n_dif = nrow(x$spec$dif),
    n_parameters = values$n_parameters %||% NA_integer_,
    log_likelihood = values$log_likelihood %||% NA_real_,
    converged = x$convergence$converged %||% NA,
    iterations = x$convergence$iterations %||% values$n_step %||% NA_integer_,
    delta = x$convergence$delta %||% values$delta %||% NA_real_,
    score_min_estimated = fit_summary$score_min_estimated[[1L]] %||% NA_integer_,
    score_max_estimated = fit_summary$score_max_estimated[[1L]] %||% NA_integer_,
    warnings = length(x$warnings %||% character()),
    stringsAsFactors = FALSE
  )
}

tidy.gRm_screen <- function(x, candidates = c("selected", "all"), ...) {
  candidates <- match.arg(candidates)
  terms <- model_terms(x)
  rows <- rbind_fill(
    transform_for_tidy(terms$ld),
    transform_for_tidy(terms$dif),
    transform_for_tidy(terms$score_effects)
  )
  if (!nrow(rows)) {
    return(rows)
  }
  rows$inference <- x$inference
  if (candidates == "selected" && "status" %in% names(rows)) {
    rows <- rows[rows$status %in% "selected", , drop = FALSE]
  }
  rows
}

glance.gRm_screen <- function(x, ...) {
  detail <- details(x)
  bh <- detail$tables$bh_thresholds %||% data.frame()
  all_terms <- tidy.gRm_screen(x, candidates = "all")
  selected <- tidy.gRm_screen(x, candidates = "selected")
  data.frame(
    inference = x$inference,
    n_tests = first_or(bh$n_tests, NA_integer_),
    n_ld_candidates = count_term_type(all_terms, "ld"),
    n_dif_candidates = count_term_type(all_terms, "dif"),
    n_score_effect_candidates = count_term_type(all_terms, "score_effect"),
    n_ld_selected = count_term_type(selected, "ld"),
    n_dif_selected = count_term_type(selected, "dif"),
    n_score_effects_selected = count_term_type(selected, "score_effect"),
    bh_fdr_05 = first_or(bh$p_value[bh$fdr == "0.05"], NA_real_),
    bh_fdr_01 = first_or(bh$p_value[bh$fdr == "0.01"], NA_real_),
    bh_fdr_001 = first_or(bh$p_value[bh$fdr == "0.001"], NA_real_),
    warnings = length(x$warnings %||% character()),
    unmodeled = length(x$unmodeled %||% character()),
    stringsAsFactors = FALSE
  )
}

summary.gRm_item_analysis <- function(object, details = FALSE, detail = NULL, ...) {
  summary_detail_output(object, details, detail, glance(object))
}

summary.gRm_gllrm_fit <- function(object, details = FALSE, detail = NULL, ...) {
  summary_detail_output(object, details, detail, glance(object))
}

summary.gRm_screen <- function(object, details = FALSE, detail = NULL, ...) {
  summary_detail_output(object, details, detail, glance(object))
}

details.gRm_item_analysis <- function(x, name = NULL, ...) {
  project <- x$project
  score_data <- item_analysis_score_data(x)

  tables <- list(
    items = data.frame(
      name = project$items$name,
      label = project$items$label_code,
      categories = project$items$raw_max,
      stringsAsFactors = FALSE
    ),
    exogenous = data.frame(
      name = project$backgrounds$name,
      label = project$backgrounds$label_code,
      categories = project$backgrounds$raw_max,
      stringsAsFactors = FALSE
    ),
    score_distribution = score_distribution_table(
      score_data$total_score,
      score_data$max_score
    ),
    estimation_score_distribution = item_analysis_estimation_score_distribution(x)
  )
  select_detail(new_gRm_details(
    tables,
    object_class = class(x),
    metadata = result_metadata(x)
  ), name)
}

details.gRm_gllrm_fit <- function(x, name = NULL, ...) {
  values <- x$values %||% list()
  tables <- item_parameter_detail_tables(values)
  tables$model_terms <- rbind_fill(
    transform_for_tidy(model_terms(x)$ld),
    transform_for_tidy(model_terms(x)$dif)
  )
  if (inherits(values, "gRm_active_gllrm_values")) {
    tables <- c(tables, gllrm_active_detail_tables(values))
  }
  select_detail(new_gRm_details(
    tables,
    object_class = class(x),
    metadata = result_metadata(x)
  ), name)
}

details.gRm_screen <- function(x, name = NULL, ...) {
  terms <- model_terms(x)
  bh <- x$values$bh %||% list()
  value_tables <- if (inherits(x$values, "gRm_screen_j_values")) {
    details(x$values)$tables
  } else {
    list()
  }
  tables <- list(
    model_terms = tidy.gRm_screen(x, candidates = "all"),
    local_dependence_terms = transform_for_tidy(terms$ld),
    dif_terms = transform_for_tidy(terms$dif),
    score_effect_terms = transform_for_tidy(terms$score_effects),
    bh_thresholds = data.frame(
      fdr = c("0.05", "0.01", "0.001"),
      p_value = c(bh$fdr_05 %||% NA_real_, bh$fdr_01 %||% NA_real_, bh$fdr_001 %||% NA_real_),
      n_tests = bh$n_tests %||% NA_integer_,
      stringsAsFactors = FALSE
    )
  )
  value_tables <- value_tables[setdiff(names(value_tables), names(tables))]
  tables <- c(tables, value_tables)
  select_detail(new_gRm_details(
    tables,
    object_class = class(x),
    metadata = result_metadata(x)
  ), name)
}

details.gRm_dif_tests_values <- function(x, name = NULL, ...) {
  tables <- list(
    tests = transform_for_tidy(x$tests %||% data.frame()),
    selected = transform_for_tidy(detail_selected_by_bh(x$tests, x$bh_critical_p)),
    bh_thresholds = detail_bh_threshold_table(x, "missing_dif")
  )
  if (!is.null(x$active_tests)) {
    tables$active_tests <- transform_for_tidy(x$active_tests)
  }
  select_detail(new_gRm_details(tables, object_class = class(x)), name)
}

details.gRm_local_independence_values <- function(x, name = NULL, ...) {
  select_detail(new_gRm_details(list(
    tests = transform_for_tidy(x$tests %||% data.frame()),
    selected = transform_for_tidy(detail_selected_by_bh(x$tests, x$bh_critical_p)),
    bh_thresholds = detail_bh_threshold_table(x, "missing_ld")
  ), object_class = class(x)), name)
}

details.gRm_global_homogeneity_values <- function(x, name = NULL, ...) {
  tables <- list(
    summary = list_to_one_row(x$summary %||% list()),
    groups = transform_for_tidy(x$score_groups %||% data.frame()),
    items = transform_for_tidy(x$items %||% data.frame())
  )
  if (!is.null(x$item_parameters)) {
    tables$item_parameters <- details(x$item_parameters)$tables$item_statistics %||% data.frame()
  }
  select_detail(new_gRm_details(tables, object_class = class(x)), name)
}

details.gRm_exo_select_values <- function(x, name = NULL, ...) {
  select_detail(new_gRm_details(list(
    exogenous = transform_for_tidy(x$exogenous %||% data.frame()),
    missing = transform_for_tidy(x$missing %||% data.frame()),
    score_distribution = transform_for_tidy(x$score_distribution %||% data.frame()),
    score_summary = list_to_one_row(x$score_summary %||% list()),
    score_effect_tests = transform_for_tidy(x$screen %||% data.frame()),
    score_effect_selected = transform_for_tidy(x$selected %||% data.frame()),
    bh_thresholds = list_to_one_row(x$bh %||% list())
  ), object_class = class(x)), name)
}

details.gRm_item_parameters_values <- function(x, name = NULL, ...) {
  select_detail(new_gRm_details(item_parameter_detail_tables(x), object_class = class(x)), name)
}

details.gRm_active_gllrm_values <- function(x, name = NULL, ...) {
  tables <- c(item_parameter_detail_tables(x), gllrm_active_detail_tables(x))
  select_detail(new_gRm_details(tables, object_class = class(x)), name)
}

details.gRm_item_fits_values <- function(x, name = NULL, ...) {
  extended <- x$extended %||% list()
  tables <- list(
    statistics = transform_for_tidy(x$items %||% data.frame()),
    compact = transform_for_tidy(x$side_file %||% data.frame()),
    bh_thresholds = result_named_numeric_table(x$bh_limits %||% numeric(), "threshold", "p_value"),
    score_n = transform_for_tidy(extended$score_n %||% data.frame()),
    score_level_fit = transform_for_tidy(extended$scores %||% data.frame()),
    item_fit_summaries = transform_for_tidy(extended$summaries %||% data.frame()),
    restscore = transform_for_tidy(extended$restscore_tables %||% x$restscore_tables %||% data.frame()),
    local_restscore = transform_for_tidy(extended$local_restscore_tables %||% data.frame()),
    local_gamma = transform_for_tidy(extended$local_gamma$rows %||% data.frame()),
    local_gamma_bh_thresholds = transform_for_tidy(extended$local_gamma$limits %||% data.frame())
  )
  if (!is.null(x$component_gamma)) {
    tables$component_gamma <- transform_for_tidy(x$component_gamma)
  }
  if (!is.null(extended$component_restscore_tables)) {
    tables$component_restscore <- transform_for_tidy(extended$component_restscore_tables)
  }
  select_detail(new_gRm_details(tables, object_class = class(x)), name)
}

details.gRm_screen_j_values <- function(x, name = NULL, ...) {
  bh <- x$bh %||% list()
  select_detail(new_gRm_details(list(
    marginal_item = screen_j_marginal_item_detail(x),
    marginal_restscore = screen_j_marginal_restscore_detail(x),
    marginal_exogenous = screen_j_exogenous_detail(x, "marginal"),
    partial_item = screen_j_partial_item_detail(x),
    partial_exogenous = screen_j_exogenous_detail(x, "partial"),
    weighted_partial_gamma = screen_j_weighted_partial_gamma_detail(x),
    spurious_dif_iterations = screen_j_bias_iterations_detail(
      x$model$spurious_dif$rows %||% list(),
      "spurious_dif",
      x$items,
      x$backgrounds
    ),
    multiple_dif_iterations = screen_j_bias_iterations_detail(
      x$model$multiple_dif$rows %||% list(),
      "multiple_dif",
      x$items,
      x$backgrounds
    ),
    score_effect_tests = screen_j_score_effect_tests_detail(x),
    score_effect_selected = screen_j_score_effect_selected_detail(x),
    local_dependence_stepwise = screen_j_stepwise_ld_detail(x),
    dif_summary = screen_j_dif_summary_detail(x),
    bh_thresholds = data.frame(
      fdr = c("0.05", "0.01", "0.001"),
      p_value = c(bh$fdr_05 %||% NA_real_, bh$fdr_01 %||% NA_real_, bh$fdr_001 %||% NA_real_),
      n_tests = bh$n_tests %||% NA_integer_,
      stringsAsFactors = FALSE
    )
  ), object_class = class(x)), name)
}

details.gRm_items_select_values <- function(x, name = NULL, ...) {
  select_detail(new_gRm_details(list(
    items = transform_for_tidy(x$items %||% data.frame()),
    score_distribution = transform_for_tidy(x$score_distribution %||% data.frame()),
    score_summary = list_to_one_row(x$score_summary %||% list()),
    score_groups = transform_for_tidy(x$score_groups %||% data.frame())
  ), object_class = class(x)), name)
}

print.gRm_details <- function(x, ...) {
  cat("<gRm_details>\n")
  if (length(x$source_results)) {
    cat("  source results: ", paste(x$source_results, collapse = ", "), "\n", sep = "")
  }
  if (nrow(x$index)) {
    cat("  tables:\n")
    for (row in seq_len(nrow(x$index))) {
      cat(sprintf("    %-32s %5d x %-5d\n", x$index$name[[row]], x$index$n_rows[[row]], x$index$n_columns[[row]]))
    }
  } else {
    cat("  tables: none\n")
  }
  invisible(x)
}

select_detail <- function(collection, name = NULL) {
  if (is.null(name)) {
    return(collection)
  }
  if (length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-missing detail name.", call. = FALSE)
  }
  if (!name %in% names(collection$tables)) {
    stop(
      "Unknown detail `", name, "`. Available details: ",
      paste(names(collection$tables), collapse = ", "),
      call. = FALSE
    )
  }
  collection$tables[[name]]
}

summary_detail_output <- function(object, show_details, detail, compact) {
  if (!is.null(detail)) {
    return(details(object, detail))
  }
  if (isTRUE(show_details)) {
    return(details(object))
  }
  compact
}

new_gRm_details <- function(tables,
                                     object_class = character(),
                                     source_results = character(),
                                     metadata = list()) {
  tables <- normalize_detail_list(tables)
  out <- list(
    object_class = object_class,
    source_results = source_results,
    tables = tables,
    index = detail_index(tables),
    metadata = metadata
  )
  class(out) <- c("gRm_details", "list")
  out
}

normalize_detail_list <- function(tables) {
  if (is.null(tables)) {
    tables <- list()
  }
  if (!is.list(tables)) {
    stop("`tables` must be a named list of data frames.", call. = FALSE)
  }
  if (length(tables) && (is.null(names(tables)) || any(!nzchar(names(tables))))) {
    stop("Every result table must have a non-empty name.", call. = FALSE)
  }
  tables <- tables[vapply(tables, is.data.frame, logical(1L))]
  lapply(tables, function(x) {
    rownames(x) <- NULL
    x
  })
}

detail_index <- function(tables) {
  data.frame(
    name = names(tables),
    topic = detail_topic(names(tables)),
    n_rows = vapply(tables, nrow, integer(1L)),
    n_columns = vapply(tables, ncol, integer(1L)),
    status = vapply(tables, function(x) attr(x, "result_status", exact = TRUE) %||% "complete", character(1L)),
    stringsAsFactors = FALSE
  )
}

detail_topic <- function(names) {
  out <- rep("other", length(names))
  out[grepl("^item_|^items$|^exogenous$|^score_", names)] <- "model"
  out[grepl("^screen_|_terms$|^bh_thresholds$|^marginal_|^partial_|^weighted_partial|^spurious_dif|^multiple_dif|^local_dependence_stepwise$|^dif_summary$", names)] <- "screen"
  out[grepl("^missing_ld", names)] <- "missing_ld"
  out[grepl("^missing_dif", names)] <- "missing_dif"
  out[grepl("^global_homogeneity", names)] <- "global_homogeneity"
  out[grepl("^score_effect", names)] <- "score_effects"
  out[names %in% c("results", "warnings", "unmodeled")] <- "provenance"
  out
}

details_from_named_values <- function(values) {
  out <- list()
  for (id in names(values)) {
    value <- values[[id]]
    if (is.null(value)) {
      next
    }
    value_tables <- details(value)$tables
    names(value_tables) <- detail_public_names(id, names(value_tables))
    out <- c(out, value_tables)
  }
  out
}

detail_public_names <- function(result_name, table_names) {
  prefix <- switch(
    result_name,
    item_fit = "item_fit",
    missing_ld = "missing_ld",
    missing_dif = "missing_dif",
    global_homogeneity = "global_homogeneity",
    score_effects = "score_effect",
    screen = "screen",
    result_name
  )
  ifelse(grepl(paste0("^", prefix), table_names), table_names, paste(prefix, table_names, sep = "_"))
}

missing_detail_tables <- function(result_name, existing) {
  required <- switch(
    result_name,
    missing_ld = c("missing_ld_tests", "missing_ld_selected", "missing_ld_bh_thresholds"),
    missing_dif = c("missing_dif_tests", "missing_dif_selected", "missing_dif_bh_thresholds"),
    item_fit = c("item_fit_statistics"),
    global_homogeneity = c("global_homogeneity_summary", "global_homogeneity_groups", "global_homogeneity_items"),
    score_effects = c("score_effect_tests", "score_effect_selected", "score_effect_bh_thresholds"),
    character()
  )
  missing <- setdiff(required, existing)
  out <- rep(list(data.frame()), length(missing))
  out <- lapply(out, function(x) {
    attr(x, "result_status") <- "unavailable"
    x
  })
  stats::setNames(out, missing)
}

result_metadata <- function(x) {
  list(
    source_trace = x$source_trace %||% character(),
    warnings = x$warnings %||% character(),
    unmodeled = x$unmodeled %||% character()
  )
}

detail_bh_threshold_table <- function(value, result = NA_character_) {
  data.frame(
    result = result,
    fdr = "0.05",
    p_value = value$bh_critical_p %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

detail_selected_by_bh <- function(tests, threshold) {
  if (!is.data.frame(tests) || !"p_value" %in% names(tests) || is.na(threshold)) {
    return(data.frame())
  }
  tests[tests$p_value <= threshold, , drop = FALSE]
}

list_to_one_row <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (!is.list(x) || !length(x)) {
    return(data.frame())
  }
  data.frame(as.list(x), check.names = FALSE, stringsAsFactors = FALSE)
}

result_matrix_table <- function(x, key_name, value_name) {
  if (!is.matrix(x) || !length(x)) {
    return(data.frame())
  }
  idx <- which(!is.na(x), arr.ind = TRUE)
  data.frame(
    item = rownames(x)[idx[, 1L]],
    key = colnames(x)[idx[, 2L]],
    value = as.numeric(x[idx]),
    stringsAsFactors = FALSE
  ) |>
    stats::setNames(c("item", key_name, value_name))
}

result_named_numeric_table <- function(x, name_col, value_col) {
  if (!length(x)) {
    return(data.frame())
  }
  out <- data.frame(
    name = names(x),
    value = as.numeric(x),
    stringsAsFactors = FALSE
  )
  stats::setNames(out, c(name_col, value_col))
}

rbind_fill <- function(...) {
  parts <- list(...)
  parts <- parts[vapply(parts, nrow, integer(1L)) > 0L]
  if (!length(parts)) {
    return(data.frame())
  }
  cols <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(x) {
    missing <- setdiff(cols, names(x))
    for (col in missing) {
      x[[col]] <- NA_character_
    }
    x[cols]
  })
  do.call(rbind, parts)
}

transform_for_tidy <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) {
    return(data.frame())
  }
  rownames(x) <- NULL
  x
}
