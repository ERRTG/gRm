first_or <- function(x, default) {
  if (is.null(x) || length(x) == 0L) {
    return(default)
  }
  x[[1L]]
}

count_term_type <- function(rows, type) {
  if (!is.data.frame(rows) || !nrow(rows)) {
    return(0L)
  }
  term_col <- if ("term_type" %in% names(rows)) rows$term_type else rows$type
  sum(term_col %in% type, na.rm = TRUE)
}

empty_detail <- function(columns) {
  columns <- as.character(columns)
  out <- stats::setNames(rep(list(logical()), length(columns)), columns)
  as.data.frame(out, stringsAsFactors = FALSE)
}

bind_detail_rows <- function(rows, columns = character()) {
  rows <- rows[vapply(rows, is.data.frame, logical(1L))]
  rows <- rows[vapply(rows, nrow, integer(1L)) > 0L]
  if (!length(rows)) {
    return(empty_detail(columns))
  }
  do.call(rbind_fill, rows)
}

detail_matrix_long <- function(mat, row_name, col_name, value_name, keep_diagonal = TRUE) {
  if (!is.matrix(mat) || !length(mat)) {
    return(empty_detail(c(row_name, col_name, value_name)))
  }
  idx <- which(!is.na(mat), arr.ind = TRUE)
  if (!keep_diagonal && nrow(idx)) {
    idx <- idx[idx[, 1L] != idx[, 2L], , drop = FALSE]
  }
  if (!nrow(idx)) {
    return(empty_detail(c(row_name, col_name, value_name)))
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

detail_named_vector <- function(x, name_col, value_col) {
  if (!length(x)) {
    return(empty_detail(c(name_col, value_col)))
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

item_analysis_score_data <- function(x) {
  project <- x$project
  items <- project$items
  raw <- project$raw_data
  n <- nrow(raw)
  scores <- matrix(NA_integer_, nrow = n, ncol = nrow(items))
  colnames(scores) <- items$name
  valid <- matrix(FALSE, nrow = n, ncol = nrow(items))
  colnames(valid) <- items$name

  for (item_index in seq_len(nrow(items))) {
    values <- as.integer(raw[, items$position[[item_index]]])
    valid[, item_index] <- !is.na(values) &
      values >= 1L &
      values <= items$raw_max[[item_index]]
    scores[valid[, item_index], item_index] <- values[valid[, item_index]] - 1L
  }

  complete_items <- if (ncol(valid)) apply(valid, 1L, all) else rep(TRUE, n)
  total_score <- rep(NA_integer_, n)
  total_score[complete_items] <- rowSums(scores[complete_items, , drop = FALSE])

  list(
    scores = scores,
    item_valid = valid,
    complete_items = complete_items,
    total_score = total_score,
    max_score = sum(items$raw_max - 1L)
  )
}

item_analysis_exogenous_complete <- function(x) {
  backgrounds <- x$project$backgrounds
  raw <- x$project$raw_data
  if (!nrow(backgrounds)) {
    return(rep(TRUE, nrow(raw)))
  }
  complete <- rep(TRUE, nrow(raw))
  for (exo_index in seq_len(nrow(backgrounds))) {
    values <- as.integer(raw[, backgrounds$position[[exo_index]]])
    complete <- complete &
      !is.na(values) &
      values >= 1L &
      values <= backgrounds$raw_max[[exo_index]]
  }
  complete
}

score_distribution_table <- function(scores, max_score) {
  keep <- !is.na(scores) & scores >= 0L & scores <= max_score
  data.frame(
    score = seq.int(0L, max_score),
    n = as.integer(tabulate(scores[keep] + 1L, nbins = max_score + 1L)),
    stringsAsFactors = FALSE
  )
}

item_analysis_estimation_score_distribution <- function(x) {
  bundle <- build_item_parameters_bundle(x$project)
  max_score <- sum(x$project$items$raw_max - 1L)
  score_distribution_table(bundle$data$score[bundle$data$status == 1L], max_score)
}

item_analysis_item_summary <- function(x) {
  project <- x$project
  score_data <- item_analysis_score_data(x)
  items <- project$items
  complete <- score_data$complete_items

  rows <- vector("list", nrow(items))
  for (item_index in seq_len(nrow(items))) {
    item_scores <- score_data$scores[, item_index]
    valid <- score_data$item_valid[, item_index]
    rows[[item_index]] <- data.frame(
      item = items$name[[item_index]],
      item_label = items$label_code[[item_index]],
      item_name = items$name[[item_index]],
      categories = items$raw_max[[item_index]],
      score_min = 0L,
      score_max = items$raw_max[[item_index]] - 1L,
      n_valid_item = sum(valid),
      mean_item_score = if (any(valid)) mean(item_scores[valid]) else NA_real_,
      n_complete_items = sum(complete),
      mean_item_score_complete = if (any(complete)) mean(item_scores[complete]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

item_analysis_glance_row <- function(x) {
  score_data <- item_analysis_score_data(x)
  exo_complete <- item_analysis_exogenous_complete(x)
  complete_items <- score_data$complete_items
  observed <- score_data$total_score[complete_items]

  data.frame(
    n_rows = nrow(x$project$raw_data),
    n_items = length(x$items),
    n_exogenous = length(x$exogenous),
    n_complete_items = sum(complete_items),
    n_complete_item_exogenous = sum(complete_items & exo_complete),
    n_missing_items = sum(!complete_items),
    n_missing_exogenous = sum(!exo_complete),
    score_type = x$score$type %||% NA_character_,
    score_min_possible = 0L,
    score_max_possible = score_data$max_score,
    score_min_observed = if (length(observed)) min(observed) else NA_integer_,
    score_max_observed = if (length(observed)) max(observed) else NA_integer_,
    score_groups = paste(x$score_groups %||% integer(), collapse = ","),
    stringsAsFactors = FALSE
  )
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

item_gamma_detail <- function(x) {
  rows <- detail_matrix_long(x$item_gamma %||% matrix(numeric(), nrow = 0L), "item", "item_score", "gamma")
  if (nrow(rows)) rows$item_score <- as.integer(rows$item_score)
  rows
}

item_threshold_detail <- function(x) {
  rows <- result_matrix_table(x$thresholds %||% matrix(numeric(), nrow = 0L), "score", "threshold")
  if (nrow(rows)) rows$score <- as.integer(rows$score)
  rows
}

item_score_matrix_detail <- function(mat, value_name) {
  rows <- detail_matrix_long(mat %||% matrix(numeric(), nrow = 0L), "item", "item_score", value_name)
  if (nrow(rows)) rows$item_score <- as.integer(rows$item_score)
  rows
}

item_effect_detail <- function(x) {
  ice <- detail_named_vector(x$ice_item_effect %||% numeric(), "item", "ice_item_effect")
  mice <- detail_named_vector(x$mice_item_effect %||% numeric(), "item", "mice_item_effect")
  merge(ice, mice, by = "item", all = TRUE, sort = FALSE)
}

item_parameter_detail_tables <- function(x) {
  location_rows <- data.frame(
    item = names(x$locations %||% numeric()),
    location = as.numeric(x$locations %||% numeric()),
    stringsAsFactors = FALSE
  )
  tables <- list(
    fit_summary = item_parameter_fit_summary(x),
    input_stats = list_to_one_row(x$input_stats %||% list()),
    item_gamma = item_gamma_detail(x),
    thresholds = item_threshold_detail(x),
    locations = location_rows,
    ice = item_score_matrix_detail(x$ice, "ice"),
    mice = item_score_matrix_detail(x$mice, "mice"),
    item_effects = item_effect_detail(x),
    item_statistics = transform_for_tidy(x$item_statistics %||% data.frame())
  )
  if (is.data.frame(x$ice_fields)) {
    tables$ice_fields <- transform_for_tidy(x$ice_fields)
  }
  tables
}

screen_inference_label <- function(x) {
  if (isTRUE(x$exact)) "exact" else "asymptotic"
}

screen_j_marginal_item_detail <- function(x) {
  items <- x$items
  rows <- list()
  for (row_item in seq_len(nrow(items))) {
    for (col_item in seq_len(nrow(items))) {
      if (row_item == col_item) next
      rows[[length(rows) + 1L]] <- data.frame(
        item = items$name[[row_item]],
        item_label = items$label_code[[row_item]],
        item_name = items$name[[row_item]],
        item2 = items$name[[col_item]],
        item2_label = items$label_code[[col_item]],
        item2_name = items$name[[col_item]],
        gamma = x$marginal$item_gamma[row_item, col_item],
        p_value = x$marginal$item_p[row_item, col_item],
        inference = screen_inference_label(x),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_detail_rows(rows, c("item", "item2", "gamma", "p_value", "inference"))
}

screen_j_marginal_restscore_detail <- function(x) {
  items <- x$items
  data.frame(
    item = items$name,
    item_label = items$label_code,
    item_name = items$name,
    gamma = as.numeric(x$marginal$rest_gamma),
    p_value = as.numeric(x$marginal$rest_p),
    inference = screen_inference_label(x),
    stringsAsFactors = FALSE
  )
}

screen_j_exogenous_detail <- function(x, section = c("marginal", "partial")) {
  section <- match.arg(section)
  values <- x[[section]]
  items <- x$items
  backgrounds <- x$backgrounds
  rows <- list()
  if (!nrow(backgrounds)) {
    return(empty_detail(c("section", "target", "item", "exo", "statistic", "p_value")))
  }

  for (exo_index in seq_len(nrow(backgrounds))) {
    for (item_index in seq_len(nrow(items))) {
      status <- if (!is.null(x$model$item_bias_status) && length(x$model$item_bias_status)) {
        x$model$item_bias_status[item_index, exo_index]
      } else {
        0L
      }
      kind <- values$exo_kind[[exo_index]] %||% x$marginal$exo_kind[[exo_index]]
      rows[[length(rows) + 1L]] <- data.frame(
        section = section,
        target = "item",
        item = items$name[[item_index]],
        item_label = items$label_code[[item_index]],
        item_name = items$name[[item_index]],
        exo = backgrounds$name[[exo_index]],
        exo_label = backgrounds$label_code[[exo_index]],
        exo_name = backgrounds$name[[exo_index]],
        statistic_name = kind,
        statistic = values$exo_stat[item_index, exo_index],
        p_value = values$exo_p[item_index, exo_index],
        ppq = if (!is.null(values$exo_ppq)) values$exo_ppq[item_index, exo_index] else NA_real_,
        pmq = if (!is.null(values$exo_pmq)) values$exo_pmq[item_index, exo_index] else NA_real_,
        selected = status > 0L,
        status = status,
        inference = screen_inference_label(x),
        stringsAsFactors = FALSE
      )
    }
    if (section == "marginal") {
      rows[[length(rows) + 1L]] <- data.frame(
        section = section,
        target = "score",
        item = NA_character_,
        item_label = NA_character_,
        item_name = NA_character_,
        exo = backgrounds$name[[exo_index]],
        exo_label = backgrounds$label_code[[exo_index]],
        exo_name = backgrounds$name[[exo_index]],
        statistic_name = x$marginal$exo_kind[[exo_index]],
        statistic = x$marginal$score_stat[[exo_index]],
        p_value = x$marginal$score_p[[exo_index]],
        ppq = NA_real_,
        pmq = NA_real_,
        selected = FALSE,
        status = 0L,
        inference = screen_inference_label(x),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_detail_rows(rows)
}

screen_j_partial_item_detail <- function(x) {
  items <- x$items
  rows <- list()
  for (row_item in seq_len(nrow(items))) {
    for (col_item in seq_len(nrow(items))) {
      if (row_item == col_item) next
      selected <- if (!is.null(x$model$local_dependence$matrix)) {
        x$model$local_dependence$matrix[row_item, col_item]
      } else {
        FALSE
      }
      rows[[length(rows) + 1L]] <- data.frame(
        item = items$name[[row_item]],
        item_label = items$label_code[[row_item]],
        item_name = items$name[[row_item]],
        item2 = items$name[[col_item]],
        item2_label = items$label_code[[col_item]],
        item2_name = items$name[[col_item]],
        gamma = x$partial$item_gamma[row_item, col_item],
        p_value = x$partial$item_p[row_item, col_item],
        ppq = x$partial$item_ppq[row_item, col_item],
        pmq = x$partial$item_pmq[row_item, col_item],
        selected = selected,
        inference = screen_inference_label(x),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_detail_rows(rows)
}

screen_j_weighted_partial_gamma_detail <- function(x) {
  items <- x$items
  rows <- list()
  if (nrow(items) < 2L) {
    return(empty_detail(c("item", "item2", "wpg_gamma")))
  }
  for (row_item in seq_len(nrow(items) - 1L)) {
    for (col_item in seq.int(row_item + 1L, nrow(items))) {
      selected <- if (!is.null(x$model$local_dependence$matrix)) {
        x$model$local_dependence$matrix[row_item, col_item]
      } else {
        FALSE
      }
      rows[[length(rows) + 1L]] <- data.frame(
        item = items$name[[row_item]],
        item_label = items$label_code[[row_item]],
        item_name = items$name[[row_item]],
        item2 = items$name[[col_item]],
        item2_label = items$label_code[[col_item]],
        item2_name = items$name[[col_item]],
        gamma_item_to_item2 = x$partial$item_gamma[row_item, col_item],
        p_value_item_to_item2 = x$partial$item_p[row_item, col_item],
        gamma_item2_to_item = x$partial$item_gamma[col_item, row_item],
        p_value_item2_to_item = x$partial$item_p[col_item, row_item],
        wpg_gamma = x$partial$weighted_gamma[row_item, col_item],
        selected = selected,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_detail_rows(rows)
}

screen_j_bias_iterations_detail <- function(entries, kind, items, backgrounds) {
  rows <- list()
  for (entry_index in seq_along(entries)) {
    entry <- entries[[entry_index]]
    iterations <- entry$iterations %||% list()
    for (iteration_index in seq_along(iterations)) {
      iteration <- iterations[[iteration_index]]
      if (!is.data.frame(iteration$rows) || !nrow(iteration$rows)) next
      rows_df <- iteration$rows
      rows_df$analysis <- kind
      rows_df$entry_index <- entry_index
      rows_df$iteration <- iteration_index
      rows_df$excluded <- if (length(iteration$excluded)) iteration$excluded else NA_integer_
      rows_df$excluded_label <- iteration$excluded_label %||% NA_character_
      rows_df$excluded_name <- iteration$excluded_name %||% NA_character_
      if (identical(kind, "spurious_dif") && !is.null(entry$background)) {
        rows_df$exo <- backgrounds$name[[entry$background]]
        rows_df$exo_label <- backgrounds$label_code[[entry$background]]
        rows_df$exo_name <- backgrounds$name[[entry$background]]
      }
      if (identical(kind, "multiple_dif") && !is.null(entry$item)) {
        rows_df$item <- items$name[[entry$item]]
        rows_df$item_label <- items$label_code[[entry$item]]
        rows_df$item_name <- items$name[[entry$item]]
      }
      rows[[length(rows) + 1L]] <- rows_df
    }
  }
  bind_detail_rows(rows)
}

screen_j_dif_summary_detail <- function(x) {
  items <- x$items
  backgrounds <- x$backgrounds
  status <- x$model$item_bias_status
  rows <- list()
  if (!nrow(backgrounds) || is.null(status) || !length(status)) {
    return(empty_detail(c("item", "exo", "statistic", "p_value", "status")))
  }
  for (item_index in seq_len(nrow(items))) {
    for (exo_index in seq_len(nrow(backgrounds))) {
      if (status[item_index, exo_index] <= 0L) next
      statistic_name <- x$partial$exo_kind[[exo_index]]
      rows[[length(rows) + 1L]] <- data.frame(
        item = items$name[[item_index]],
        item_label = items$label_code[[item_index]],
        item_name = items$name[[item_index]],
        exo = backgrounds$name[[exo_index]],
        exo_label = backgrounds$label_code[[exo_index]],
        exo_name = backgrounds$name[[exo_index]],
        statistic_name = statistic_name,
        statistic = x$partial$exo_stat[item_index, exo_index],
        p_value = x$partial$exo_p[item_index, exo_index],
        gamma = if (identical(statistic_name, "Gamma")) x$partial$exo_stat[item_index, exo_index] else NA_real_,
        status = status[item_index, exo_index],
        possible_dif_source = status[item_index, exo_index] > 1L,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_detail_rows(rows)
}

screen_j_score_effect_tests_detail <- function(x) {
  rows <- transform_for_tidy(x$model$score_effects$rows %||% data.frame())
  if (!nrow(rows)) {
    return(rows)
  }
  rows$term_type <- "score_effect"
  rows$inference <- screen_inference_label(x)
  rows
}

screen_j_score_effect_selected_detail <- function(x) {
  rows <- screen_j_score_effect_tests_detail(x)
  if (!nrow(rows) || !"selected" %in% names(rows)) {
    return(rows[0L, , drop = FALSE])
  }
  rows[rows$selected %in% TRUE, , drop = FALSE]
}

screen_j_stepwise_ld_detail <- function(x) {
  rows <- transform_for_tidy(x$model$local_dependence$rows %||% data.frame())
  if (!nrow(rows)) {
    return(rows)
  }
  rows$term_type <- "ld"
  rows$selection_rule <- "screen_j_stepwise_weighted_partial_gamma"
  rows
}

details_to_test_rows <- function(tables) {
  candidates <- tables[grepl("_(tests|selected|active_tests)$", names(tables))]
  rows <- lapply(names(candidates), function(table_name) {
    tab <- transform_for_tidy(candidates[[table_name]])
    if (!nrow(tab)) {
      return(tab)
    }
    tab$result <- sub("_(tests|selected|active_tests)$", "", table_name)
    tab
  })
  bind_detail_rows(rows)
}

detail_public_prefix <- function(result_name) {
  switch(
    result_name,
    item_fit = "item_fit",
    missing_ld = "missing_ld",
    missing_dif = "missing_dif",
    global_homogeneity = "global_homogeneity",
    score_effects = "score_effect",
    screen = "screen",
    result_name
  )
}
