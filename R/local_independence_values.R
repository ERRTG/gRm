local_independence_values_cache <- new.env(parent = emptyenv())

local_independence_cache_key <- function(project, max_step, max_delta, jobs) {
  raw_sum <- sum(as.matrix(project$raw_data), na.rm = TRUE)
  input_dir <- if (!is.null(project$paths$input_dir)) project$paths$input_dir else ""
  item_signature <- paste(
    project$items$name,
    project$items$position,
    project$items$raw_max,
    collapse = "|"
  )
  paste(
    normalizePath(input_dir, mustWork = FALSE),
    nrow(project$raw_data),
    ncol(project$raw_data),
    item_signature,
    raw_sum,
    max_step,
    max_delta,
    jobs,
    sep = "\r"
  )
}

#' Count observed margins for one active local-dependence candidate
#'
#' Counts the observed item-by-item margin used when DIGRAM tests a missing
#' local-dependence term between two items.
#'
#' @param bundle Source-shaped bundle.
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @return Observed pair margin matrix.
#' @keywords internal
active_ld_counts <- function(bundle, item1, item2) {
  items <- bundle$model$items
  data <- bundle$data
  item1_name <- items$name[[item1]]
  item2_name <- items$name[[item2]]
  item1_max <- items$raw_max[[item1]]
  item2_max <- items$raw_max[[item2]]
  observed_ld <- matrix(
    0L,
    nrow = item1_max,
    ncol = item2_max,
    dimnames = list(
      as.character(seq.int(0L, item1_max - 1L)),
      as.character(seq.int(0L, item2_max - 1L))
    )
  )

  for (row_index in which(data$status == 1L)) {
    score1 <- data[[item1_name]][[row_index]]
    score2 <- data[[item2_name]][[row_index]]
    observed_ld[score1 + 1L, score2 + 1L] <-
      observed_ld[score1 + 1L, score2 + 1L] + 1L
  }

  list(observed_ld = observed_ld)
}

#' Source Goodman-Kruskal gamma count totals
#'
#' @param tab Two-way integer table.
#' @return Gamma, PPQ, and PMQ totals.
#' @keywords internal
local_independence_source_gamma_counts <- function(tab) {
  n_rows <- nrow(tab)
  n_cols <- ncol(tab)
  cumulative <- t(apply(apply(tab, 2L, cumsum), 1L, cumsum))
  if (n_rows == 1L) {
    cumulative <- matrix(cumulative, nrow = 1L)
  }
  total <- sum(tab)
  p <- 0
  q <- 0

  cell_sum <- function(row_to, col_to) {
    if (row_to <= 0L || col_to <= 0L) {
      0
    } else {
      cumulative[row_to, col_to]
    }
  }

  for (row in seq_len(n_rows)) {
    for (col in seq_len(n_cols)) {
      less_less <- cell_sum(row - 1L, col - 1L)
      greater_greater <- total - cell_sum(row, n_cols) - cell_sum(n_rows, col) + cell_sum(row, col)
      less_greater <- cell_sum(row - 1L, n_cols) - cell_sum(row - 1L, col)
      greater_less <- cell_sum(n_rows, col - 1L) - cell_sum(row, col - 1L)
      p <- p + tab[row, col] * (less_less + greater_greater)
      q <- q + tab[row, col] * (less_greater + greater_less)
    }
  }

  ppq <- p + q
  pmq <- p - q
  list(gamma = if (ppq > 0) pmq / ppq else 0, ppq = ppq, pmq = pmq)
}

#' Calculate source-shaped local-independence WPG gamma
#'
#' @param project DIGRAM project.
#' @return Raw item matrix and source score vectors.
#' @keywords internal
build_local_independence_gamma_context <- function(project) {
  items <- project$items
  item_matrix <- matrix(0L, nrow = nrow(project$raw_data), ncol = nrow(items))
  complete_items <- rep(TRUE, nrow(item_matrix))
  for (item_index in seq_len(nrow(items))) {
    item_matrix[, item_index] <- as.integer(project$raw_data[, items$position[[item_index]]])
    complete_items <- complete_items &
      item_matrix[, item_index] >= 1L &
      item_matrix[, item_index] <= items$raw_max[[item_index]]
  }
  item_score <- rowSums(sweep(item_matrix, 2L, 1L, "-"))

  list(
    item_matrix = item_matrix,
    complete_items = complete_items,
    item_score = item_score
  )
}

#' Calculate directed source-shaped local-independence gamma counts
#'
#' @param context Context from [build_local_independence_gamma_context()].
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @param items Item metadata.
#' @return PPQ and PMQ totals.
#' @keywords internal
local_independence_directed_gamma_counts <- function(context, item1, item2, items) {
  rest_score <- context$item_score - (context$item_matrix[, item1] - 1L)
  x <- context$item_matrix[, item1]
  y <- context$item_matrix[, item2]
  valid <- context$complete_items
  ppq <- 0
  pmq <- 0

  for (score in sort(unique(rest_score[valid]))) {
    in_stratum <- valid & rest_score == score
    tab <- matrix(
      0L,
      nrow = items$raw_max[[item1]],
      ncol = items$raw_max[[item2]]
    )
    index <- x[in_stratum] + (y[in_stratum] - 1L) * nrow(tab)
    tab[] <- tabulate(index, nbins = length(tab))
    stats <- local_independence_source_gamma_counts(tab)
    ppq <- ppq + stats$ppq
    pmq <- pmq + stats$pmq
  }

  list(ppq = ppq, pmq = pmq)
}

#' Calculate source-shaped local-independence WPG gamma
#'
#' Ports the weighted partial gamma used by `DGRirtD.MissingLD`: the two
#' directed item-screening partial gamma count totals are pooled before PMQ/PPQ.
#'
#' @param project DIGRAM project.
#' @param context Context from [build_local_independence_gamma_context()].
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @return Source weighted partial gamma.
#' @keywords internal
local_independence_wpg_gamma <- function(project, context, item1, item2) {
  items <- project$items
  forward <- local_independence_directed_gamma_counts(context, item1, item2, items)
  backward <- local_independence_directed_gamma_counts(context, item2, item1, items)
  ppq <- forward$ppq + backward$ppq
  pmq <- forward$pmq + backward$pmq
  if (ppq > 0) pmq / ppq else 0
}

#' Build a fast source-shaped local-dependence estimation context
#'
#' @param bundle Source-shaped bundle.
#' @return Cached integer matrices and source bounds for active-LD fitting.
#' @keywords internal
build_active_ld_context <- function(bundle) {
  items <- bundle$model$items
  data <- bundle$data

  item_matrix <- matrix(0L, nrow = nrow(data), ncol = nrow(items))
  for (item_index in seq_len(nrow(items))) {
    item_matrix[, item_index] <- as.integer(data[[items$name[[item_index]]]])
  }

  item_raw_max <- as.integer(items$raw_max)
  valid_rows <- which(data$status == 1L)
  list(
    n_items = nrow(items),
    item_raw_max = item_raw_max,
    max_total_score = as.integer(bundle$model$max_total_score),
    item_score_columns = lapply(item_raw_max, seq_len),
    item_score_values = lapply(item_raw_max, function(raw_max) seq.int(0L, raw_max - 1L)),
    item_score_reference = gllrm_item_score_reference(item_matrix, valid_rows, items),
    item_matrix = item_matrix,
    score = as.integer(data$score),
    valid_rows = valid_rows
  )
}

#' Count observed margins for one active LD candidate using cached data
#'
#' @param context Active LD context from [build_active_ld_context()].
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @return Observed pair margin matrix.
#' @keywords internal
active_ld_counts_context <- function(context, item1, item2) {
  item1_max <- context$item_raw_max[[item1]]
  item2_max <- context$item_raw_max[[item2]]
  observed_ld <- matrix(
    0L,
    nrow = item1_max,
    ncol = item2_max,
    dimnames = list(
      as.character(seq.int(0L, item1_max - 1L)),
      as.character(seq.int(0L, item2_max - 1L))
    )
  )

  score1 <- context$item_matrix[context$valid_rows, item1]
  score2 <- context$item_matrix[context$valid_rows, item2]
  linear <- score1 + score2 * item1_max + 1L
  observed_ld[] <- tabulate(linear, nbins = length(observed_ld))
  list(observed_ld = observed_ld)
}

#' Convolve component weights excluding one component
#'
#' @param components List of score-weight vectors.
#' @param excluded One-based component index to exclude.
#' @param max_total_score Maximum score.
#' @return Score polynomial for all non-excluded components.
#' @keywords internal
convolve_components_except <- function(components, excluded, max_total_score) {
  result <- c(1, numeric(max_total_score))
  for (component_index in seq_along(components)) {
    if (component_index == excluded) {
      next
    }
    result <- convolve_score_vectors(result, components[[component_index]], max_total_score)
  }
  result
}

#' Build active local-dependence component weights
#'
#' @param bundle Source-shaped bundle.
#' @param item_gamma Item gamma matrix.
#' @param ld_gamma Candidate LD gamma matrix.
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @return Component metadata and score polynomials.
#' @keywords internal
build_active_ld_components <- function(bundle, item_gamma, ld_gamma, item1, item2) {
  items <- bundle$model$items
  max_total_score <- bundle$model$max_total_score
  components <- list()
  component_item <- integer(0)
  pair_component <- NA_integer_

  for (item_index in seq_len(nrow(items))) {
    if (item_index %in% c(item1, item2)) {
      next
    }
    scores <- seq.int(0L, items$raw_max[[item_index]] - 1L)
    components[[length(components) + 1L]] <-
      as.numeric(item_gamma[item_index, as.character(scores)])
    component_item <- c(component_item, item_index)
  }

  pair_weights <- numeric(max_total_score + 1L)
  for (score1 in seq.int(0L, items$raw_max[[item1]] - 1L)) {
    for (score2 in seq.int(0L, items$raw_max[[item2]] - 1L)) {
      # Source trace: an active LD candidate multiplies the two item score
      # parameters by the item-by-item LD cell parameter before conditioning on
      # the total score.
      pair_weights[[score1 + score2 + 1L]] <- pair_weights[[score1 + score2 + 1L]] +
        item_gamma[item1, as.character(score1)] *
        item_gamma[item2, as.character(score2)] *
        ld_gamma[score1 + 1L, score2 + 1L]
    }
  }
  components[[length(components) + 1L]] <- pair_weights
  component_item <- c(component_item, NA_integer_)
  pair_component <- length(components)

  list(
    components = components,
    component_item = component_item,
    pair_component = pair_component,
    full_gamma = Reduce(
      function(left, right) convolve_score_vectors(left, right, max_total_score),
      components,
      init = c(1, numeric(max_total_score))
    )
  )
}

#' Calculate expected margins for one active local-dependence candidate
#'
#' @param bundle Source-shaped bundle.
#' @param base_counts Base Rasch counts.
#' @param candidate_counts Candidate LD counts.
#' @param item_gamma Item gamma matrix.
#' @param ld_gamma Candidate LD gamma matrix.
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @return Expected item and LD margins.
#' @keywords internal
calculate_active_ld_expected <- function(bundle, base_counts, candidate_counts, item_gamma, ld_gamma, item1, item2) {
  items <- bundle$model$items
  max_total_score <- bundle$model$max_total_score
  expected_items <- item_gamma
  expected_items[,] <- 0
  expected_ld <- ld_gamma
  expected_ld[,] <- 0
  score_counts <- base_counts$score_counts
  component_info <- build_active_ld_components(bundle, item_gamma, ld_gamma, item1, item2)
  full_gamma <- component_info$full_gamma

  without_component <- vector("list", length(component_info$components))
  for (component_index in seq_along(component_info$components)) {
    without_component[[component_index]] <- convolve_components_except(
      component_info$components, component_index, max_total_score
    )
  }

  for (score in seq.int(0L, max_total_score)) {
    n_score <- score_counts[[score + 1L]]
    denominator <- full_gamma[[score + 1L]]
    if (n_score == 0L || denominator <= 0) {
      next
    }

    for (component_index in seq_along(component_info$components)) {
      item_index <- component_info$component_item[[component_index]]
      if (is.na(item_index)) {
        next
      }
      rest_gamma <- without_component[[component_index]]
      for (item_score in seq.int(0L, items$raw_max[[item_index]] - 1L)) {
        if (score < item_score) {
          next
        }
        numerator <- item_gamma[item_index, as.character(item_score)] *
          rest_gamma[[score - item_score + 1L]]
        expected_items[item_index, as.character(item_score)] <-
          expected_items[item_index, as.character(item_score)] +
          n_score * numerator / denominator
      }
    }

    rest_pair_gamma <- without_component[[component_info$pair_component]]
    for (score1 in seq.int(0L, items$raw_max[[item1]] - 1L)) {
      for (score2 in seq.int(0L, items$raw_max[[item2]] - 1L)) {
        pair_score <- score1 + score2
        if (score < pair_score) {
          next
        }
        numerator <- item_gamma[item1, as.character(score1)] *
          item_gamma[item2, as.character(score2)] *
          ld_gamma[score1 + 1L, score2 + 1L] *
          rest_pair_gamma[[score - pair_score + 1L]]
        expected <- n_score * numerator / denominator
        expected_items[item1, as.character(score1)] <-
          expected_items[item1, as.character(score1)] + expected
        expected_items[item2, as.character(score2)] <-
          expected_items[item2, as.character(score2)] + expected
        expected_ld[score1 + 1L, score2 + 1L] <-
          expected_ld[score1 + 1L, score2 + 1L] + expected
      }
    }
  }

  list(expected_items = expected_items, expected_ld = expected_ld)
}

#' Calculate expected margins for one active LD candidate using cached bounds
#'
#' @param context Active LD context from [build_active_ld_context()].
#' @param base_counts Base Rasch counts.
#' @param item_gamma Item gamma matrix.
#' @param ld_gamma Candidate LD gamma matrix.
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @return Expected item and LD margins.
#' @keywords internal
calculate_active_ld_expected_context <- function(context, base_counts, item_gamma, ld_gamma, item1, item2) {
  max_total_score <- context$max_total_score
  expected_items <- item_gamma
  expected_items[,] <- 0
  expected_ld <- ld_gamma
  expected_ld[,] <- 0

  components <- list()
  component_item <- integer(0)
  for (item_index in seq_len(context$n_items)) {
    if (item_index %in% c(item1, item2)) {
      next
    }
    components[[length(components) + 1L]] <-
      as.numeric(item_gamma[item_index, context$item_score_columns[[item_index]]])
    component_item <- c(component_item, item_index)
  }

  pair_weights <- numeric(max_total_score + 1L)
  for (score1 in context$item_score_values[[item1]]) {
    for (score2 in context$item_score_values[[item2]]) {
      pair_weights[[score1 + score2 + 1L]] <- pair_weights[[score1 + score2 + 1L]] +
        item_gamma[item1, score1 + 1L] *
        item_gamma[item2, score2 + 1L] *
        ld_gamma[score1 + 1L, score2 + 1L]
    }
  }
  components[[length(components) + 1L]] <- pair_weights
  component_item <- c(component_item, NA_integer_)
  pair_component <- length(components)

  prefix <- vector("list", length(components) + 1L)
  suffix <- vector("list", length(components) + 1L)
  prefix[[1L]] <- c(1, numeric(max_total_score))
  for (component_index in seq_along(components)) {
    prefix[[component_index + 1L]] <- convolve_score_vectors(
      prefix[[component_index]],
      components[[component_index]],
      max_total_score
    )
  }
  suffix[[length(components) + 1L]] <- c(1, numeric(max_total_score))
  for (component_index in rev(seq_along(components))) {
    suffix[[component_index]] <- convolve_score_vectors(
      components[[component_index]],
      suffix[[component_index + 1L]],
      max_total_score
    )
  }
  full_gamma <- prefix[[length(components) + 1L]]

  score_counts <- base_counts$score_counts
  active_score_indices <- which(score_counts > 0)
  if (length(active_score_indices) == 0L) {
    return(list(expected_items = expected_items, expected_ld = expected_ld))
  }
  denominator <- full_gamma[active_score_indices]
  usable <- denominator > 0
  if (!any(usable)) {
    return(list(expected_items = expected_items, expected_ld = expected_ld))
  }
  active_score_indices <- active_score_indices[usable]
  scores <- active_score_indices - 1L
  score_weight <- score_counts[active_score_indices] / denominator[usable]

  for (component_index in seq_along(components)) {
    item_index <- component_item[[component_index]]
    rest_gamma <- convolve_score_vectors(
      prefix[[component_index]],
      suffix[[component_index + 1L]],
      max_total_score
    )
    if (!is.na(item_index)) {
      for (item_score in context$item_score_values[[item_index]]) {
        score_usable <- scores >= item_score
        if (!any(score_usable)) {
          next
        }
        item_col <- item_score + 1L
        expected_items[item_index, item_col] <- expected_items[item_index, item_col] +
          item_gamma[item_index, item_col] *
          sum(score_weight[score_usable] * rest_gamma[scores[score_usable] - item_score + 1L])
      }
      next
    }

    for (score1 in context$item_score_values[[item1]]) {
      for (score2 in context$item_score_values[[item2]]) {
        pair_score <- score1 + score2
        score_usable <- scores >= pair_score
        if (!any(score_usable)) {
          next
        }
        weight <- item_gamma[item1, score1 + 1L] *
          item_gamma[item2, score2 + 1L] *
          ld_gamma[score1 + 1L, score2 + 1L]
        expected <- weight * sum(
          score_weight[score_usable] * rest_gamma[scores[score_usable] - pair_score + 1L]
        )
        expected_items[item1, score1 + 1L] <- expected_items[item1, score1 + 1L] + expected
        expected_items[item2, score2 + 1L] <- expected_items[item2, score2 + 1L] + expected
        expected_ld[score1 + 1L, score2 + 1L] <- expected_ld[score1 + 1L, score2 + 1L] + expected
      }
    }
  }

  list(expected_items = expected_items, expected_ld = expected_ld)
}

#' Fit one active local-dependence candidate model
#'
#' @param bundle Source-shaped bundle.
#' @param base_counts Base Rasch counts.
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @param max_step Maximum IPF iterations.
#' @param max_delta Convergence threshold.
#' @param initial_item_gamma_matrix Optional starting item gamma matrix.
#' @return Candidate fit state.
#' @keywords internal
fit_active_ld_candidate <- function(bundle, base_counts, item1, item2, max_step = 5000L, max_delta = 0.0001, initial_item_gamma_matrix = NULL, context = NULL, candidate_counts = NULL) {
  items <- bundle$model$items
  if (is.null(context)) {
    context <- build_active_ld_context(bundle)
  }
  item_gamma <- if (is.null(initial_item_gamma_matrix)) initial_item_gamma(bundle) else initial_item_gamma_matrix
  ld_gamma <- matrix(
    1,
    nrow = items$raw_max[[item1]],
    ncol = items$raw_max[[item2]],
    dimnames = list(
      as.character(seq.int(0L, items$raw_max[[item1]] - 1L)),
      as.character(seq.int(0L, items$raw_max[[item2]] - 1L))
    )
  )
  if (is.null(candidate_counts)) {
    candidate_counts <- active_ld_counts_context(context, item1, item2)
  }
  delta <- 0
  converged <- FALSE

  for (step in seq_len(max_step)) {
    expected <- calculate_active_ld_expected_context(
      context, base_counts, item_gamma, ld_gamma, item1, item2
    )
    delta <- 0

    for (item_index in seq_len(nrow(items))) {
      for (item_score in seq.int(0L, items$raw_max[[item_index]] - 1L)) {
        observed <- base_counts$item_counts[item_index, as.character(item_score)]
        fitted <- expected$expected_items[item_index, as.character(item_score)]
        if (observed > 0 && fitted > 0) {
          # Source trace: IPF updates multiply each item score parameter by the
          # observed/fitted sufficient-margin ratio.
          ratio <- observed / fitted
          item_gamma[item_index, as.character(item_score)] <-
            item_gamma[item_index, as.character(item_score)] * ratio
          delta <- max(delta, abs(fitted - observed))
        } else if (observed == 0) {
          item_gamma[item_index, as.character(item_score)] <- 0
        }
      }
    }

    for (score1 in seq.int(0L, items$raw_max[[item1]] - 1L)) {
      for (score2 in seq.int(0L, items$raw_max[[item2]] - 1L)) {
        observed <- candidate_counts$observed_ld[score1 + 1L, score2 + 1L]
        fitted <- expected$expected_ld[score1 + 1L, score2 + 1L]
        if (observed > 0 && fitted > 0) {
          # Source trace: the active IJ/LD parameter is updated by the same
          # observed/fitted margin ratio as the Pascal IPF step.
          ratio <- observed / fitted
          ld_gamma[score1 + 1L, score2 + 1L] <-
            ld_gamma[score1 + 1L, score2 + 1L] * ratio
          delta <- max(delta, abs(fitted - observed))
        } else if (observed == 0) {
          ld_gamma[score1 + 1L, score2 + 1L] <- 0
        }
      }
    }

    ld_reference <- as.integer(context$item_score_reference %||% 0L) + 1L
    ld_gamma <- adjust_ld_gamma_source_reference_details(
      candidate_counts$observed_ld,
      ld_gamma,
      i_ref = ld_reference,
      j_ref = ld_reference,
      preserve_current_ties = TRUE
    )$adjusted
    item_gamma <- adjust_item_gammas_source_scale(bundle, item_gamma)
    if (delta <= max_delta) {
      converged <- TRUE
      break
    }
  }

  list(
    item_gamma = item_gamma,
    ld_gamma = ld_gamma,
    candidate_counts = candidate_counts,
    delta = delta,
    converged = converged
  )
}

#' Apply source reference adjustment to one LD gamma table
#'
#' Ports the LD part of `SourceRaschCore.pas::AdjustDependencyParameters`,
#' choosing the densest observed row/column as reference and converting raw IJ
#' cell parameters to source reference-relative parameters.
#'
#' @param observed_ld Observed item-pair margin.
#' @param ld_gamma Current LD gamma matrix.
#' @return Adjusted LD gamma matrix.
#' @keywords internal
adjust_ld_gamma_source_reference <- function(observed_ld, ld_gamma) {
  adjust_ld_gamma_source_reference_details(observed_ld, ld_gamma)$adjusted
}

adjust_ld_gamma_source_reference_details <- function(observed_ld,
                                                     ld_gamma,
                                                     i_ref = 1L,
                                                     j_ref = 1L,
                                                     preserve_current_ties = TRUE) {
  rows <- nrow(ld_gamma)
  cols <- ncol(ld_gamma)
  i_ref <- min(max(as.integer(i_ref), 1L), rows)
  j_ref <- min(max(as.integer(j_ref), 1L), cols)
  i_cells <- rowSums(observed_ld > 0)
  j_cells <- colSums(observed_ld > 0)

  if (i_cells[[i_ref]] < cols) {
    if (i_cells[[rows]] == cols) {
      i_ref <- rows
    } else if (isTRUE(preserve_current_ties)) {
      for (candidate in seq_len(rows)) {
        if (i_cells[[candidate]] > i_cells[[i_ref]]) {
          i_ref <- candidate
        }
      }
    } else {
      i_ref <- which.max(i_cells)
    }
  }
  if (j_cells[[j_ref]] < rows) {
    if (j_cells[[cols]] == rows) {
      j_ref <- cols
    } else if (isTRUE(preserve_current_ties)) {
      for (candidate in seq_len(cols)) {
        if (j_cells[[candidate]] > j_cells[[j_ref]]) {
          j_ref <- candidate
        }
      }
    } else {
      j_ref <- which.max(j_cells)
    }
  }

  ix_ref <- ld_gamma[i_ref, j_ref]
  positive_ref <- ix_ref > 0
  complete_ref <- i_cells[[i_ref]] == cols && j_cells[[j_ref]] == rows
  i_pos <- sum(i_cells > 0)
  j_pos <- sum(j_cells > 0)

  i_first <- ld_gamma[, j_ref]
  for (score1 in seq_len(rows)) {
    if (i_first[[score1]] == 0) {
      observed_cols <- which(observed_ld[score1, ] > 0)
      if (length(observed_cols) > 0L) {
        i_first[[score1]] <- ld_gamma[score1, observed_cols[[1L]]]
      }
    }
  }
  j_first <- ld_gamma[i_ref, ]
  for (score2 in seq_len(cols)) {
    if (j_first[[score2]] == 0) {
      observed_rows <- which(observed_ld[, score2] > 0)
      if (length(observed_rows) > 0L) {
        j_first[[score2]] <- ld_gamma[observed_rows[[1L]], score2]
      }
    }
  }

  old_ld <- ld_gamma
  adjusted <- ld_gamma
  adjusted_with_item_factors <- FALSE
  if (complete_ref || (i_pos > 1L && j_pos > 1L && positive_ref)) {
    adjusted_with_item_factors <- TRUE
    for (score1 in seq_len(rows)) {
      for (score2 in seq_len(cols)) {
        if (observed_ld[score1, score2] > 0) {
          # Source trace: adjusted cell = old_cell * reference /
          # (first_in_column * first_in_row).
          denominator <- j_first[[score2]] * i_first[[score1]]
          if (denominator > 0) {
            adjusted[score1, score2] <- old_ld[score1, score2] * ix_ref / denominator
          } else {
            adjusted[score1, score2] <- 0
          }
        } else {
          adjusted[score1, score2] <- 0
        }
      }
    }
  } else {
    adjusted[,] <- 0
    adjusted[observed_ld > 0] <- 1
  }
  list(
    adjusted = adjusted,
    i_first = i_first,
    j_first = j_first,
    i_ref = i_ref,
    j_ref = j_ref,
    reference = ix_ref,
    adjusted_with_item_factors = adjusted_with_item_factors
  )
}

#' Calculate active local-dependence negative log likelihood
#'
#' @param bundle Source-shaped bundle.
#' @param fit Candidate fit from [fit_active_ld_candidate()].
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @return Negative conditional log likelihood.
#' @keywords internal
active_ld_loglike <- function(bundle, fit, item1, item2) {
  data <- bundle$data
  items <- bundle$model$items
  item1_name <- items$name[[item1]]
  item2_name <- items$name[[item2]]
  full_gamma <- build_active_ld_components(
    bundle, fit$item_gamma, fit$ld_gamma, item1, item2
  )$full_gamma

  loglike <- 0
  for (row_index in which(data$status == 1L)) {
    score <- data$score[[row_index]]
    product_gamma <- 1
    for (item_index in seq_len(nrow(items))) {
      item_score <- data[[items$name[[item_index]]]][[row_index]]
      product_gamma <- product_gamma * fit$item_gamma[item_index, as.character(item_score)]
    }
    product_gamma <- product_gamma *
      fit$ld_gamma[data[[item1_name]][[row_index]] + 1L, data[[item2_name]][[row_index]] + 1L]
    probability <- product_gamma / full_gamma[[score + 1L]]
    if (probability > 0) {
      loglike <- loglike - log(probability)
    }
  }
  loglike
}

#' Calculate active local-dependence negative log likelihood using cached data
#'
#' @param bundle Source-shaped bundle.
#' @param context Active LD context from [build_active_ld_context()].
#' @param fit Candidate fit from [fit_active_ld_candidate()].
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @return Negative conditional log likelihood.
#' @keywords internal
active_ld_loglike_context <- function(bundle, context, fit, item1, item2) {
  full_gamma <- build_active_ld_components(
    bundle, fit$item_gamma, fit$ld_gamma, item1, item2
  )$full_gamma

  valid_rows <- context$valid_rows
  product_gamma <- rep(1, length(valid_rows))
  for (item_index in seq_len(context$n_items)) {
    item_scores <- context$item_matrix[valid_rows, item_index]
    product_gamma <- product_gamma * fit$item_gamma[cbind(item_index, item_scores + 1L)]
  }
  product_gamma <- product_gamma *
    fit$ld_gamma[cbind(
      context$item_matrix[valid_rows, item1] + 1L,
      context$item_matrix[valid_rows, item2] + 1L
    )]
  probability <- product_gamma / full_gamma[context$score[valid_rows] + 1L]
  probability <- probability[probability > 0]
  -sum(log(probability))
}

#' Compute DIGRAM local-independence candidate test values
#'
#' Computes the item-pair likelihood-ratio tests used in DIGRAM's
#' `check-local-independence.txt` and
#' `check-local-independence-extended.txt` reports. The implementation follows
#' `docs/example_LOCAL_INDEPENDENCE_SOURCE_TRACE.md`, especially the
#' `DGRirtD.pas` `MissingLD` branch.
#'
#' @param project A source-shaped DIGRAM project list, such as the `project`
#'   component returned by [gRm()] or [read_digram_project()].
#' @param max_step Maximum source Rasch/GLLRM estimation iterations.
#' @param max_delta Convergence threshold.
#' @param jobs Number of parallel candidate LD fits.
#' @return A `gRm_local_independence_values` object.
#' @keywords internal
local_independence_values <- function(project, max_step = 5000L, max_delta = 0.0001, jobs = min(32L, parallel::detectCores(logical = TRUE), 128L)) {
  if ((inherits(project, "gRm_fit") || inherits(project, "gRm_gllrm_fit")) && inherits(project$values, "gRm_active_gllrm_values")) {
    return(active_gllrm_local_independence_values(project, max_step = max_step, max_delta = max_delta, jobs = jobs))
  }
  if (!is.list(project)) {
    stop("`project` must be a DIGRAM project list.", call. = FALSE)
  }
  if (is.null(project$items) || nrow(project$items) < 2L) {
    stop("Local-independence tests require at least two item variables.", call. = FALSE)
  }
  jobs <- max(1L, min(as.integer(jobs), 128L))
  cache_key <- local_independence_cache_key(project, max_step, max_delta, jobs)
  if (exists(cache_key, envir = local_independence_values_cache, inherits = FALSE)) {
    return(get(cache_key, envir = local_independence_values_cache, inherits = FALSE))
  }

  bundle <- build_item_parameters_bundle(project)
  context <- build_active_ld_context(bundle)
  gamma_context <- build_local_independence_gamma_context(project)
  base_fit <- fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)
  base_counts <- base_fit$counts
  base_loglike <- base_rasch_loglike(bundle, base_fit$item_gamma)
  items <- project$items
  candidates <- expand.grid(
    item1_index = seq_len(nrow(items)),
    item2_index = seq_len(nrow(items))
  )
  candidates <- candidates[candidates$item1_index < candidates$item2_index, , drop = FALSE]
  candidates <- candidates[order(candidates$item1_index, candidates$item2_index), , drop = FALSE]
  candidate_counts <- vector("list", nrow(candidates))
  for (candidate_row in seq_len(nrow(candidates))) {
    candidate_counts[[candidate_row]] <- active_ld_counts_context(
      context,
      candidates$item1_index[[candidate_row]],
      candidates$item2_index[[candidate_row]]
    )
  }

  fit_one <- function(candidate_row) {
    item1 <- candidates$item1_index[[candidate_row]]
    item2 <- candidates$item2_index[[candidate_row]]
    candidate_fit <- fit_active_ld_candidate(
      bundle, base_counts, item1, item2,
      max_step = max_step, max_delta = max_delta,
      initial_item_gamma_matrix = base_fit$item_gamma,
      context = context,
      candidate_counts = candidate_counts[[candidate_row]]
    )
    candidate_loglike <- active_ld_loglike_context(bundle, context, candidate_fit, item1, item2)

    # Source trace: DGRirtD MissingLD reports
    # lr := 2*abs(Raschloglike-Raschloglike1).
    clr <- 2 * abs(base_loglike - candidate_loglike)
    # Source trace: for example four-category items, IJ df is the candidate
    # parameter increment (4 - 1) * (4 - 1) = 9.
    df <- (items$raw_max[[item1]] - 1L) * (items$raw_max[[item2]] - 1L)

    data.frame(
      pair_label = paste0(items$label_code[[item1]], items$label_code[[item2]]),
      item1_label = items$label_code[[item1]],
      item2_label = items$label_code[[item2]],
      item1_name = items$name[[item1]],
      item2_name = items$name[[item2]],
      chi_square = clr,
      degrees_of_freedom = df,
      # Source trace: p := pfchi(df, lr).
      p_value = source_pfchi(df, clr),
      wpg_gamma = local_independence_wpg_gamma(project, gamma_context, item1, item2),
      converged = candidate_fit$converged,
      delta = candidate_fit$delta,
      stringsAsFactors = FALSE
    )
  }

  if (.Platform$OS.type == "unix" && jobs > 1L) {
    fit_rows <- parallel::mclapply(seq_len(nrow(candidates)), fit_one, mc.cores = jobs, mc.preschedule = FALSE)
  } else {
    fit_rows <- lapply(seq_len(nrow(candidates)), fit_one)
  }
  tests <- do.call(rbind, fit_rows)
  bh_critical <- source_bh_critical(tests$p_value, 0.05)

  result <- structure(
    list(
      tests = tests,
      bh_critical_p = bh_critical,
      suggested_ld = tests$pair_label[tests$p_value <= bh_critical],
      max_step = as.integer(max_step),
      max_delta = max_delta
    ),
    class = "gRm_local_independence_values"
  )
  assign(cache_key, result, envir = local_independence_values_cache)
  result
}

active_gllrm_local_independence_values <- function(fit, max_step = 5000L, max_delta = 0.0001, jobs = min(32L, parallel::detectCores(logical = TRUE), 128L)) {
  context <- fit$fit$context
  items <- context$items
  candidates <- active_gllrm_li_candidates(context)

  if (nrow(candidates) == 0L) {
    empty <- data.frame(
      pair_label = character(),
      item1_label = character(),
      item2_label = character(),
      item1_name = character(),
      item2_name = character(),
      chi_square = numeric(),
      degrees_of_freedom = integer(),
      p_value = numeric(),
      wpg_gamma = numeric(),
      converged = logical(),
      delta = numeric(),
      stringsAsFactors = FALSE
    )
    return(structure(
      list(
        tests = empty,
        bh_critical_p = 0,
        suggested_ld = character(),
        max_step = as.integer(max_step),
        max_delta = max_delta
      ),
      class = "gRm_local_independence_values"
    ))
  }

  jobs <- max(1L, min(as.integer(jobs), 128L, nrow(candidates)))
  base_loglike <- fit$fit$log_likelihood
  fit_one <- function(candidate_row) {
    item1 <- candidates$item1_index[[candidate_row]]
    item2 <- candidates$item2_index[[candidate_row]]
    candidate_fit <- fit_gllrm_with_added_ld(
      fit,
      item1 = item1,
      item2 = item2,
      max_step = max_step,
      max_delta = max_delta
    )
    candidate_loglike <- candidate_fit$fit$log_likelihood
    clr <- 2 * abs(base_loglike - candidate_loglike)
    ld_index <- active_gllrm_context_ld_index(candidate_fit$fit$context, item1, item2)
    df <- source_ij_observed_df(
      candidate_fit$fit$context$observed_ld[[ld_index]]
    )
    data.frame(
      pair_label = paste0(items$label_code[[item1]], items$label_code[[item2]]),
      item1_label = items$label_code[[item1]],
      item2_label = items$label_code[[item2]],
      item1_name = items$name[[item1]],
      item2_name = items$name[[item2]],
      chi_square = clr,
      degrees_of_freedom = df,
      p_value = source_pfchi(df, clr),
      converged = isTRUE(candidate_fit$convergence$converged),
      delta = candidate_fit$fit$report_delta %||% NA_real_,
      stringsAsFactors = FALSE
    )
  }

  if (.Platform$OS.type == "unix" && jobs > 1L) {
    fit_rows <- parallel::mclapply(seq_len(nrow(candidates)), fit_one, mc.cores = jobs, mc.preschedule = FALSE)
  } else {
    fit_rows <- lapply(seq_len(nrow(candidates)), fit_one)
  }
  tests <- do.call(rbind, fit_rows)
  bh_critical <- source_bh_critical(tests$p_value, 0.05)
  structure(
    list(
      tests = tests,
      bh_critical_p = bh_critical,
      suggested_ld = tests$pair_label[tests$p_value <= bh_critical],
      max_step = as.integer(max_step),
      max_delta = max_delta
    ),
    class = "gRm_local_independence_values"
  )
}

active_gllrm_context_ld_index <- function(context, item1, item2) {
  key <- active_gllrm_ld_key(item1, item2)
  hit <- which(vapply(context$ld_specs, function(spec) {
    identical(active_gllrm_ld_key(spec$item1, spec$item2), key)
  }, logical(1L)))
  if (length(hit) != 1L) {
    stop("Could not identify the candidate LD term in the fitted GLLRM context.", call. = FALSE)
  }
  hit[[1L]]
}

active_gllrm_li_candidates <- function(context, components = NULL) {
  rows <- expand.grid(
    item1_index = seq_len(context$n_items),
    item2_index = seq_len(context$n_items)
  )
  rows <- rows[rows$item1_index < rows$item2_index, , drop = FALSE]
  active <- active_gllrm_ld_lookup(context)
  rows <- rows[
    !vapply(
      seq_len(nrow(rows)),
      function(i) active[[active_gllrm_ld_key(rows$item1_index[[i]], rows$item2_index[[i]])]] %||% FALSE,
      logical(1L)
    ),
    ,
    drop = FALSE
  ]
  rows[order(rows$item1_index, rows$item2_index), , drop = FALSE]
}

active_gllrm_ld_lookup <- function(context) {
  out <- list()
  for (spec in context$ld_specs) {
    out[[active_gllrm_ld_key(spec$item1, spec$item2)]] <- TRUE
  }
  out
}

active_gllrm_ld_key <- function(item1, item2) {
  paste(min(item1, item2), max(item1, item2), sep = ":")
}

source_ij_observed_df <- function(observed_ij) {
  nonzero_item1_scores <- sum(rowSums(observed_ij) > 0)
  nonzero_item2_scores <- sum(colSums(observed_ij) > 0)
  as.integer(max(0L, nonzero_item1_scores - 1L) * max(0L, nonzero_item2_scores - 1L))
}

active_gllrm_li_pair_test <- function(context, components, item1, item2) {
  component_of <- components$component_of
  component_to_remove <- component_of[[item1]]
  score_items <- which(component_of != component_to_remove)
  rest_score <- if (length(score_items) > 0L) {
    rowSums(context$item_matrix[, score_items, drop = FALSE])
  } else {
    rep(0L, nrow(context$item_matrix))
  }
  max_score <- sum(context$item_raw_max[score_items] - 1L)
  bias_backgrounds <- active_gllrm_li_bias_backgrounds(context, score_items)

  condition_values <- matrix(nrow = nrow(context$item_matrix), ncol = 0L)
  condition_dims <- integer()
  if (length(bias_backgrounds) > 0L) {
    condition_values <- cbind(condition_values, context$background_matrix[, bias_backgrounds, drop = FALSE])
    condition_dims <- c(condition_dims, context$background_raw_max[bias_backgrounds])
  }
  condition_values <- cbind(condition_values, rest_score)
  condition_dims <- c(condition_dims, max_score - 1L)

  valid <- seq_len(nrow(context$item_matrix)) %in% context$valid_rows
  valid <- valid & rest_score > 0L & rest_score < max_score
  screen_j_conditional_bias_test(
    x = context$item_matrix[, item1] + 1L,
    y = context$item_matrix[, item2] + 1L,
    x_dim = context$item_raw_max[[item1]],
    y_dim = context$item_raw_max[[item2]],
    condition_values = condition_values,
    condition_dims = condition_dims,
    valid = valid,
    exact = FALSE,
    native = FALSE
  )
}

active_gllrm_li_bias_backgrounds <- function(context, score_items) {
  backgrounds <- integer()
  for (spec in context$dif_specs) {
    if (spec$item %in% score_items) {
      backgrounds <- c(backgrounds, spec$background)
    }
  }
  sort(unique(backgrounds))
}
