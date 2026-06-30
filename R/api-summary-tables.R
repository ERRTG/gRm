empty_result_table <- function(columns) {
  columns <- as.character(columns)
  out <- stats::setNames(rep(list(logical()), length(columns)), columns)
  as.data.frame(out, stringsAsFactors = FALSE)
}

result_matrix_long_table <- function(mat, row_name, col_name, value_name, keep_diagonal = TRUE) {
  if (!is.matrix(mat) || !length(mat)) {
    return(empty_result_table(c(row_name, col_name, value_name)))
  }
  idx <- which(!is.na(mat), arr.ind = TRUE)
  if (!keep_diagonal && nrow(idx)) {
    idx <- idx[idx[, 1L] != idx[, 2L], , drop = FALSE]
  }
  if (!nrow(idx)) {
    return(empty_result_table(c(row_name, col_name, value_name)))
  }
  row_keys <- rownames(mat)
  col_keys <- colnames(mat)
  if (is.null(row_keys)) row_keys <- as.character(seq_len(nrow(mat)))
  if (is.null(col_keys)) col_keys <- as.character(seq_len(ncol(mat)))
  out <- data.frame(
    row = row_keys[idx[, 1L]],
    col = col_keys[idx[, 2L]],
    value = as.numeric(mat[idx]),
    stringsAsFactors = FALSE
  )
  stats::setNames(out, c(row_name, col_name, value_name))
}

result_named_vector_table <- function(x, name_col, value_col) {
  if (!length(x)) {
    return(empty_result_table(c(name_col, value_col)))
  }
  names_x <- names(x)
  if (is.null(names_x)) names_x <- as.character(seq_along(x))
  out <- data.frame(
    name = names_x,
    value = as.numeric(x),
    stringsAsFactors = FALSE
  )
  stats::setNames(out, c(name_col, value_col))
}

item_parameter_fit_summary <- function(x) {
  range <- x$observed_score_range %||% c(NA_integer_, NA_integer_)
  data.frame(
    n_step = x$n_step %||% NA_integer_,
    delta = x$delta %||% NA_real_,
    log_likelihood = x$log_likelihood %||% NA_real_,
    likelihood_n = x$likelihood_n %||% NA_integer_,
    n_parameters = x$n_parameters %||% NA_integer_,
    score_min_estimated = range[[1L]] %||% NA_integer_,
    score_max_estimated = range[[2L]] %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}

item_gamma_result_table <- function(x) {
  rows <- result_matrix_long_table(x$item_gamma %||% matrix(numeric(), nrow = 0L), "item", "item_score", "gamma")
  if (nrow(rows)) rows$item_score <- as.integer(rows$item_score)
  rows
}

item_threshold_result_table <- function(x) {
  rows <- result_matrix_table(x$thresholds %||% matrix(numeric(), nrow = 0L), "score", "threshold")
  if (nrow(rows)) rows$score <- as.integer(rows$score)
  rows
}

item_score_matrix_result_table <- function(mat, value_name) {
  rows <- result_matrix_long_table(mat %||% matrix(numeric(), nrow = 0L), "item", "item_score", value_name)
  if (nrow(rows)) rows$item_score <- as.integer(rows$item_score)
  rows
}

item_effect_result_tables <- function(x) {
  ice <- result_named_vector_table(x$ice_item_effect %||% numeric(), "item", "ice_item_effect")
  mice <- result_named_vector_table(x$mice_item_effect %||% numeric(), "item", "mice_item_effect")
  merge(ice, mice, by = "item", all = TRUE, sort = FALSE)
}

item_parameter_result_tables <- function(x) {
  location_rows <- data.frame(
    item = names(x$locations %||% numeric()),
    location = as.numeric(x$locations %||% numeric()),
    stringsAsFactors = FALSE
  )
  tables <- list(
    fit_summary = item_parameter_fit_summary(x),
    input_stats = list_to_one_row(x$input_stats %||% list()),
    item_gamma = item_gamma_result_table(x),
    thresholds = item_threshold_result_table(x),
    locations = location_rows,
    ice = item_score_matrix_result_table(x$ice, "ice"),
    mice = item_score_matrix_result_table(x$mice, "mice"),
    item_effects = item_effect_result_tables(x),
    item_statistics = normalize_summary_table(x$item_statistics %||% data.frame())
  )
  if (is.data.frame(x$ice_fields)) {
    tables$ice_fields <- normalize_summary_table(x$ice_fields)
  }
  tables
}
