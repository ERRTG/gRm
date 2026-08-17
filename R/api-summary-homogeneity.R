#' Internal public global homogeneity test helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param summary Internal `summary` value used by this helper.
#' @return The internal `public_global_homogeneity_test()` computation result.
#' @keywords internal
#' @noRd
public_global_homogeneity_test <- function(summary) {
  summary <- list_to_one_row(summary %||% list())
  data.frame(
    `Score groups` = diagnostic_column(summary, "n_groups", NA_integer_),
    Parameters = diagnostic_column(summary, "n_parameters", NA_integer_),
    `LogLik full` = diagnostic_column(summary, "full_log_likelihood", NA_real_),
    `LogLik groups` = diagnostic_column(summary, "subgroup_log_likelihood_sum", NA_real_),
    CLR = diagnostic_column(summary, "clr", NA_real_),
    Df = diagnostic_column(summary, "df", NA_integer_),
    `Pr(>CLR)` = diagnostic_column(summary, "p_value", NA_real_),
    check.names = FALSE
  )
}

#' Internal public global homogeneity score groups helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param groups Internal `groups` value used by this helper.
#' @return The internal `public_global_homogeneity_score_groups()` computation result.
#' @keywords internal
#' @noRd
public_global_homogeneity_score_groups <- function(groups) {
  if (!is.data.frame(groups)) {
    groups <- data.frame()
  }
  data.frame(
    `Score group` = global_homogeneity_score_group_labels(groups),
    Cases = diagnostic_column(groups, "n", NA_integer_),
    LogLik = diagnostic_column(groups, "log_likelihood", NA_real_),
    Converged = yes_no_display(diagnostic_column(groups, "converged", NA)),
    delta = diagnostic_column(groups, "delta", NA_real_),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Internal public global homogeneity item means helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param items Item selection or item metadata.
#' @param groups Internal `groups` value used by this helper.
#' @return The internal `public_global_homogeneity_item_means()` computation result.
#' @keywords internal
#' @noRd
public_global_homogeneity_item_means <- function(items, groups) {
  if (!is.data.frame(items)) {
    items <- data.frame()
  }
  data.frame(
    `Score group` = global_homogeneity_item_score_group_labels(items, groups),
    Item = diagnostic_column(items, "item_name", NA_character_),
    Cases = diagnostic_column(items, "n", NA_integer_),
    `Observed mean` = diagnostic_column(items, "observed_mean", NA_real_),
    `Expected mean` = diagnostic_column(items, "expected_mean", NA_real_),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Internal public global homogeneity uniform ld helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @param values Values to validate or transform.
#' @return The internal `public_global_homogeneity_uniform_ld()` computation result.
#' @keywords internal
#' @noRd
public_global_homogeneity_uniform_ld <- function(table, values) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(data.frame())
  }
  out <- data.frame(
    `Item 1` = global_homogeneity_label_names(
      diagnostic_column(table, "item1_label", NA_character_),
      global_homogeneity_context_variables(values, "items")
    ),
    `Item 2` = global_homogeneity_label_names(
      diagnostic_column(table, "item2_label", NA_character_),
      global_homogeneity_context_variables(values, "items")
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- cbind(
    out,
    global_homogeneity_uniform_gamma_columns(table, values),
    data.frame(
      Chisq = diagnostic_column(table, "chi_square", NA_real_),
      Df = diagnostic_column(table, "df", NA_integer_),
      `Pr(>Chisq)` = diagnostic_column(table, "p_value", NA_real_),
      check.names = FALSE
    )
  )
  rownames(out) <- NULL
  out
}

#' Internal public global homogeneity uniform dif helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @param values Values to validate or transform.
#' @return The internal `public_global_homogeneity_uniform_dif()` computation result.
#' @keywords internal
#' @noRd
public_global_homogeneity_uniform_dif <- function(table, values) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(data.frame())
  }
  out <- data.frame(
    Item = global_homogeneity_label_names(
      diagnostic_column(table, "item_label", NA_character_),
      global_homogeneity_context_variables(values, "items")
    ),
    Exogenous = global_homogeneity_label_names(
      diagnostic_column(table, "background_label", NA_character_),
      global_homogeneity_context_variables(values, "backgrounds")
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- cbind(
    out,
    global_homogeneity_uniform_gamma_columns(table, values),
    data.frame(
      Chisq = diagnostic_column(table, "chi_square", NA_real_),
      Df = diagnostic_column(table, "df", NA_integer_),
      `Pr(>Chisq)` = diagnostic_column(table, "p_value", NA_real_),
      check.names = FALSE
    )
  )
  rownames(out) <- NULL
  out
}

#' Internal global homogeneity default summary sections helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `global_homogeneity_default_summary_sections()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_default_summary_sections <- function(values) {
  sections <- c("test", "score_groups", "item_means")
  if (global_homogeneity_table_nrow(values$uniform_ld %||% data.frame()) > 0L) {
    sections <- c(sections, "uniform_ld")
  }
  if (global_homogeneity_table_nrow(values$uniform_dif %||% data.frame()) > 0L) {
    sections <- c(sections, "uniform_dif")
  }
  sections
}

#' Internal global homogeneity table nrow helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The internal `global_homogeneity_table_nrow()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_table_nrow <- function(table) {
  if (is.data.frame(table)) {
    nrow(table)
  } else {
    0L
  }
}

#' Internal global homogeneity uniform gamma columns helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @param values Values to validate or transform.
#' @return The internal `global_homogeneity_uniform_gamma_columns()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_uniform_gamma_columns <- function(table, values) {
  labels <- global_homogeneity_uniform_score_group_labels(table, values$score_groups %||% data.frame())
  n_groups <- length(labels)
  if (n_groups == 0L) {
    return(data.frame(row.names = seq_len(nrow(table))))
  }
  observed <- global_homogeneity_uniform_gamma_matrix(table, "observed_gamma", n_groups)
  expected <- global_homogeneity_uniform_gamma_matrix(table, "expected_gamma", n_groups)
  columns <- vector("list", 2L * n_groups)
  column_names <- character(2L * n_groups)
  for (group_index in seq_len(n_groups)) {
    observed_column <- (group_index - 1L) * 2L + 1L
    expected_column <- observed_column + 1L
    columns[[observed_column]] <- observed[, group_index]
    columns[[expected_column]] <- expected[, group_index]
    column_names[[observed_column]] <- paste("Obs gamma", labels[[group_index]])
    column_names[[expected_column]] <- paste("Exp gamma", labels[[group_index]])
  }
  out <- as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE)
  names(out) <- column_names
  out
}

#' Internal global homogeneity uniform score group labels helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @param groups Internal `groups` value used by this helper.
#' @return The internal `global_homogeneity_uniform_score_group_labels()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_uniform_score_group_labels <- function(table, groups) {
  labels <- global_homogeneity_score_group_labels(groups)
  n_groups <- max(
    length(labels),
    global_homogeneity_uniform_n_groups(table),
    0L
  )
  if (n_groups == 0L) {
    return(character())
  }
  if (length(labels) < n_groups) {
    labels <- c(labels, as.character(seq.int(length(labels) + 1L, n_groups)))
  }
  labels[seq_len(n_groups)]
}

#' Internal global homogeneity uniform n groups helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The internal `global_homogeneity_uniform_n_groups()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_uniform_n_groups <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(0L)
  }
  lengths <- integer()
  for (column in c("observed_gamma", "expected_gamma")) {
    if (column %in% names(table)) {
      lengths <- c(lengths, vapply(table[[column]], length, integer(1L)))
    }
  }
  if (length(lengths)) {
    max(lengths, na.rm = TRUE)
  } else {
    0L
  }
}

#' Internal global homogeneity uniform gamma matrix helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @param column Internal `column` value used by this helper.
#' @param n_groups Internal `n_groups` value used by this helper.
#' @return The internal `global_homogeneity_uniform_gamma_matrix()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_uniform_gamma_matrix <- function(table, column, n_groups) {
  out <- matrix(NA_real_, nrow = nrow(table), ncol = n_groups)
  if (!column %in% names(table) || n_groups == 0L) {
    return(out)
  }
  for (row_index in seq_len(nrow(table))) {
    values <- as.numeric(table[[column]][[row_index]])
    n <- min(length(values), n_groups)
    if (n > 0L) {
      out[row_index, seq_len(n)] <- values[seq_len(n)]
    }
  }
  out
}

#' Internal global homogeneity context variables helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @param kind Internal `kind` value used by this helper.
#' @return The internal `global_homogeneity_context_variables()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_context_variables <- function(values, kind = c("items", "backgrounds")) {
  kind <- match.arg(kind)
  context <- values$fit$context %||% list()
  variables <- context[[kind]] %||% NULL
  if (is.data.frame(variables)) {
    return(variables)
  }
  model <- values$bundle$model %||% list()
  variables <- model[[kind]] %||% NULL
  if (is.data.frame(variables)) {
    return(variables)
  }
  data.frame()
}

#' Internal global homogeneity label names helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param labels Source-facing labels.
#' @param variables Internal `variables` value used by this helper.
#' @return The internal `global_homogeneity_label_names()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_label_names <- function(labels, variables) {
  labels <- as.character(labels)
  if (!is.data.frame(variables) || !all(c("label_code", "name") %in% names(variables))) {
    return(labels)
  }
  index <- match(labels, as.character(variables$label_code))
  out <- as.character(variables$name[index])
  missing <- is.na(out)
  out[missing] <- labels[missing]
  out
}

#' Internal global homogeneity score group labels helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param groups Internal `groups` value used by this helper.
#' @return The internal `global_homogeneity_score_group_labels()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_score_group_labels <- function(groups) {
  if (!is.data.frame(groups) || !nrow(groups)) {
    return(character())
  }
  if (all(c("from_score", "to_score") %in% names(groups))) {
    return(mapply(score_range_label, groups$from_score, groups$to_score, USE.NAMES = FALSE))
  }
  if ("label" %in% names(groups)) {
    return(gsub("[[:space:]]+-[[:space:]]+", "-", groups$label))
  }
  as.character(diagnostic_column(groups, "group", NA_integer_))
}

#' Internal global homogeneity item score group labels helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param items Item selection or item metadata.
#' @param groups Internal `groups` value used by this helper.
#' @return The internal `global_homogeneity_item_score_group_labels()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_item_score_group_labels <- function(items, groups) {
  if (!is.data.frame(items) || !nrow(items)) {
    return(character())
  }
  item_groups <- diagnostic_column(items, "group", NA_integer_)
  if (!is.data.frame(groups) || !nrow(groups) || !"group" %in% names(groups)) {
    return(as.character(item_groups))
  }
  group_labels <- global_homogeneity_score_group_labels(groups)
  names(group_labels) <- as.character(groups$group)
  out <- unname(group_labels[as.character(item_groups)])
  missing <- is.na(out)
  out[missing] <- as.character(item_groups[missing])
  out
}

#' Internal diagnostic column helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param tests Diagnostic test rows.
#' @param column Internal `column` value used by this helper.
#' @param default Internal `default` value used by this helper.
#' @return The internal `diagnostic_column()` computation result.
#' @keywords internal
#' @noRd
diagnostic_column <- function(tests, column, default) {
  if (column %in% names(tests)) {
    return(tests[[column]])
  }
  rep(default, nrow(tests))
}
