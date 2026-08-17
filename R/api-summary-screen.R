#' Internal screen summary header helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @param tables Internal `tables` value used by this helper.
#' @return The internal `screen_summary_header()` computation result.
#' @keywords internal
#' @noRd
screen_summary_header <- function(object, tables) {
  exact_state <- object$exact_state %||% list()
  bh <- object$values$bh %||% list()
  header <- c(
    "",
    "Screening",
    paste0("  Inference: ", object$inference %||% "asymptotic")
  )
  if (isTRUE(exact_state$exact)) {
    header <- c(
      header,
      paste0("  Simulations: ", summary_scalar(object$nsim %||% exact_state$nsim %||% NA_integer_))
    )
    if (isTRUE(exact_state$sequential)) {
      header <- c(
        header,
        paste0("  Sequential limit: ", summary_scalar(exact_state$seq_limit %||% NA_integer_))
      )
    }
    header <- c(
      header,
      paste0("  Seed: ", summary_scalar(object$seed %||% exact_state$seed %||% NA_integer_))
    )
  }
  c(
    header,
    paste0("  Tested relations: ", summary_scalar(bh$n_tests %||% NA_integer_)),
    paste0("  Tested local-dependence pairs: ", nrow(tables$local_dependence)),
    paste0("  Tested DIF relations: ", nrow(tables$dif)),
    paste0("  Tested score effects: ", nrow(tables$score_effects))
  )
}

#' Internal public screen bh table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @return The internal `public_screen_bh_table()` computation result.
#' @keywords internal
#' @noRd
public_screen_bh_table <- function(object) {
  bh <- object$values$bh %||% list()
  data.frame(
    fdr = c("0.05", "0.01", "0.001"),
    p_value = c(bh$fdr_05 %||% NA_real_, bh$fdr_01 %||% NA_real_, bh$fdr_001 %||% NA_real_),
    n_tests = bh$n_tests %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}

#' Internal public screen summary tables helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @return The internal `public_screen_summary_tables()` computation result.
#' @keywords internal
#' @noRd
public_screen_summary_tables <- function(object) {
  values <- object$values %||% list()
  terms <- public_screen_terms(object)
  model_terms <- rbind_fill(terms$ld, terms$dif)
  selected_ld <- public_selected_rows(terms$ld)
  selected_dif <- public_selected_rows(terms$dif)
  list(
    local_dependence = public_screen_local_dependence_tests(values),
    dif = public_screen_dif_tests(values),
    score_effects = public_screen_score_effect_tests(values),
    selected = rbind_fill(selected_ld, selected_dif),
    selected_ld = selected_ld,
    selected_dif = selected_dif,
    model_terms = model_terms,
    bh = public_screen_bh_table(object)
  )
}

#' Internal public screen local dependence tests helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `public_screen_local_dependence_tests()` computation result.
#' @keywords internal
#' @noRd
public_screen_local_dependence_tests <- function(values) {
  p <- values$partial$item_p %||% matrix(numeric(), nrow = 0L, ncol = 0L)
  gamma <- values$partial$item_gamma %||% matrix(numeric(), nrow = 0L, ncol = 0L)
  wpg <- values$partial$weighted_gamma %||% matrix(numeric(), nrow = 0L, ncol = 0L)
  items <- values$items %||% data.frame()
  n_items <- nrow(items)
  if (!is.matrix(p) || n_items < 2L) {
    return(data.frame(
      `Item 1` = character(),
      `Item 2` = character(),
      `Gamma 1->2` = numeric(),
      `Pr(>|Gamma 1->2|)` = numeric(),
      `Gamma 2->1` = numeric(),
      `Pr(>|Gamma 2->1|)` = numeric(),
      WPG = numeric(),
      `Gamma sum` = numeric(),
      Decision = character(),
      ` ` = character(),
      check.names = FALSE
    ))
  }
  selected <- values$model$local_dependence$matrix %||%
    matrix(FALSE, nrow = n_items, ncol = n_items)
  stepwise <- values$model$local_dependence$stepwise_matrix %||% selected
  rows <- which(upper.tri(matrix(FALSE, nrow = n_items, ncol = n_items)), arr.ind = TRUE)
  reverse_rows <- cbind(rows[, "col"], rows[, "row"])
  item_names <- screen_item_names(items, n_items)
  gamma_forward <- screen_matrix_value(gamma, rows)
  gamma_reverse <- screen_matrix_value(gamma, reverse_rows)
  included <- screen_matrix_value(selected, rows, default = FALSE) %in% TRUE
  provisional <- screen_matrix_value(stepwise, rows, default = FALSE) %in% TRUE
  out <- data.frame(
    `Item 1` = item_names[rows[, "row"]],
    `Item 2` = item_names[rows[, "col"]],
    `Gamma 1->2` = gamma_forward,
    `Pr(>|Gamma 1->2|)` = screen_matrix_value(p, rows),
    `Gamma 2->1` = gamma_reverse,
    `Pr(>|Gamma 2->1|)` = screen_matrix_value(p, reverse_rows),
    WPG = screen_matrix_value(wpg, rows),
    `Gamma sum` = gamma_forward + gamma_reverse,
    Decision = ifelse(
      included,
      "included",
      ifelse(provisional, "negative LD; not included", "")
    ),
    ` ` = ifelse(included, "*", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  for (column in c("Pr(>|Gamma 1->2|)", "Pr(>|Gamma 2->1|)")) {
    out[!is.na(out[[column]]) & out[[column]] > 1, column] <- NA_real_
  }
  out
}

#' Internal public screen dif tests helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `public_screen_dif_tests()` computation result.
#' @keywords internal
#' @noRd
public_screen_dif_tests <- function(values) {
  p <- values$partial$exo_p %||% matrix(numeric(), nrow = 0L, ncol = 0L)
  stat <- values$partial$exo_stat %||% matrix(numeric(), nrow = 0L, ncol = 0L)
  kinds <- values$partial$exo_kind %||% character()
  items <- values$items %||% data.frame()
  backgrounds <- values$backgrounds %||% data.frame()
  n_items <- nrow(items)
  n_exo <- nrow(backgrounds)
  if (!is.matrix(p) || n_items == 0L || n_exo == 0L) {
    return(data.frame(
      Item = character(),
      Exogenous = character(),
      Chisq = numeric(),
      Df = integer(),
      `Pr(>Chisq)` = numeric(),
      Gamma = numeric(),
      `Pr(>|Gamma|)` = numeric(),
      ` ` = character(),
      check.names = FALSE
    ))
  }
  rows <- expand.grid(item = seq_len(n_items), exogenous = seq_len(n_exo))
  item_names <- screen_item_names(items, n_items)
  exo_names <- screen_item_names(backgrounds, n_exo)
  selected <- values$model$item_bias %||% matrix(FALSE, nrow = n_items, ncol = n_exo)
  kind <- rep_len(as.character(kinds), n_exo)
  row_kind <- kind[rows$exogenous]
  use_gamma <- row_kind == "Gamma"
  row_index <- cbind(rows$item, rows$exogenous)
  statistic <- stat[row_index]
  p_value <- p[row_index]
  data.frame(
    Item = item_names[rows$item],
    Exogenous = exo_names[rows$exogenous],
    Chisq = ifelse(use_gamma, NA_real_, statistic),
    Df = NA_integer_,
    `Pr(>Chisq)` = ifelse(use_gamma, NA_real_, p_value),
    Gamma = ifelse(use_gamma, statistic, NA_real_),
    `Pr(>|Gamma|)` = ifelse(use_gamma, p_value, NA_real_),
    ` ` = ifelse(screen_matrix_value(selected, row_index, default = FALSE) %in% TRUE, "*", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Internal public screen score effect tests helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param values Values to validate or transform.
#' @return The internal `public_screen_score_effect_tests()` computation result.
#' @keywords internal
#' @noRd
public_screen_score_effect_tests <- function(values) {
  rows <- values$model$score_effects$rows %||% data.frame()
  if (!is.data.frame(rows) || !nrow(rows)) {
    return(data.frame(
      Exogenous = character(),
      Chisq = numeric(),
      Df = integer(),
      `Pr(>Chisq)` = numeric(),
      Gamma = numeric(),
      `Pr(>|Gamma|)` = numeric(),
      ` ` = character(),
      check.names = FALSE
    ))
  }
  out <- data.frame(
    Exogenous = diagnostic_column(rows, "name", NA_character_),
    Chisq = diagnostic_column(rows, "chi_square", NA_real_),
    Df = diagnostic_column(rows, "df", NA_integer_),
    `Pr(>Chisq)` = diagnostic_column(rows, "p_chi", NA_real_),
    Gamma = diagnostic_column(rows, "gamma", NA_real_),
    `Pr(>|Gamma|)` = diagnostic_column(rows, "p_gamma", NA_real_),
    ` ` = ifelse(diagnostic_column(rows, "selected", FALSE) %in% TRUE, "*", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  exact_nsim <- diagnostic_column(rows, "exact_nsim", NA_integer_)
  if (any(!is.na(exact_nsim) & exact_nsim > 0L)) {
    out$Simulations <- exact_nsim
  }
  out
}

#' Internal screen item names helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @param n Internal `n` value used by this helper.
#' @return The internal `screen_item_names()` computation result.
#' @keywords internal
#' @noRd
screen_item_names <- function(table, n) {
  names <- if (is.data.frame(table) && "name" %in% names(table)) {
    as.character(table$name)
  } else {
    character()
  }
  if (length(names) < n) {
    names <- c(names, as.character(seq.int(length(names) + 1L, n)))
  }
  names
}

#' Internal screen matrix value helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param matrix Internal `matrix` value used by this helper.
#' @param index One-based internal index.
#' @param default Internal `default` value used by this helper.
#' @return The internal `screen_matrix_value()` computation result.
#' @keywords internal
#' @noRd
screen_matrix_value <- function(matrix, index, default = NA_real_) {
  if (!is.matrix(matrix) || !nrow(index)) {
    return(rep(default, nrow(index)))
  }
  valid <- index[, 1L] >= 1L & index[, 1L] <= nrow(matrix) &
    index[, 2L] >= 1L & index[, 2L] <= ncol(matrix)
  out <- rep(default, nrow(index))
  out[valid] <- matrix[index[valid, , drop = FALSE]]
  out
}

#' Internal public local dependence tests helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param tests Diagnostic test rows.
#' @param bh_critical_p Internal `bh_critical_p` value used by this helper.
#' @return The internal `public_local_dependence_tests()` computation result.
#' @keywords internal
#' @noRd
public_local_dependence_tests <- function(tests, bh_critical_p = NA_real_) {
  if (!is.data.frame(tests)) {
    tests <- data.frame()
  }
  converged <- diagnostic_column(tests, "converged", NA)
  p_value <- diagnostic_column(tests, "p_value", NA_real_)
  threshold <- bh_critical_p %||% NA_real_
  threshold <- threshold[[1L]]
  data.frame(
    `Item 1` = diagnostic_column(tests, "item1_name", NA_character_),
    `Item 2` = diagnostic_column(tests, "item2_name", NA_character_),
    Chisq = diagnostic_column(tests, "chi_square", NA_real_),
    Df = diagnostic_column(tests, "degrees_of_freedom", NA_integer_),
    `Pr(>Chisq)` = p_value,
    WPG = diagnostic_column(tests, "wpg_gamma", NA_real_),
    Converged = yes_no_display(converged),
    delta = diagnostic_column(tests, "delta", NA_real_),
    ` ` = ifelse(!is.na(p_value) & !is.na(threshold) & p_value <= threshold, "*", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Internal public dif tests helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param tests Diagnostic test rows.
#' @param bh_critical_p Internal `bh_critical_p` value used by this helper.
#' @return The internal `public_dif_tests()` computation result.
#' @keywords internal
#' @noRd
public_dif_tests <- function(tests, bh_critical_p = NA_real_) {
  if (!is.data.frame(tests)) {
    tests <- data.frame()
  }
  converged <- diagnostic_column(tests, "converged", NA)
  stable <- diagnostic_column(tests, "output_stable", FALSE)
  p_value <- diagnostic_column(tests, "p_value", NA_real_)
  threshold <- bh_critical_p %||% NA_real_
  threshold <- threshold[[1L]]
  data.frame(
    Item = diagnostic_column(tests, "item_name", NA_character_),
    Exogenous = diagnostic_column(tests, "background_name", NA_character_),
    Chisq = diagnostic_column(tests, "chi_square", NA_real_),
    Df = diagnostic_column(tests, "degrees_of_freedom", NA_integer_),
    `Pr(>Chisq)` = p_value,
    Gamma = diagnostic_column(tests, "gamma", NA_real_),
    Converged = yes_no_display(converged),
    Stable = yes_no_display(stable),
    delta = diagnostic_column(tests, "delta", NA_real_),
    ` ` = ifelse(!is.na(p_value) & !is.na(threshold) & p_value <= threshold, "*", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
