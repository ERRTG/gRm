#' Internal screen j exo values helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param project Encoded gRm project.
#' @param item_matrix Internal `item_matrix` value used by this helper.
#' @param item_score Internal `item_score` value used by this helper.
#' @param complete_items Internal `complete_items` value used by this helper.
#' @param items Item selection or item metadata.
#' @param backgrounds Internal `backgrounds` value used by this helper.
#' @param max_score Internal `max_score` value used by this helper.
#' @param partial Internal `partial` value used by this helper.
#' @return The internal `screen_j_exo_values()` computation result.
#' @keywords internal
#' @noRd
screen_j_exo_values <- function(project,
                                item_matrix,
                                item_score,
                                complete_items,
                                items,
                                backgrounds,
                                max_score,
                                partial) {
  n_items <- nrow(items)
  n_exo <- nrow(backgrounds)
  stat <- matrix(0, nrow = n_items, ncol = n_exo)
  p <- matrix(2, nrow = n_items, ncol = n_exo)
  ppq <- matrix(0, nrow = n_items, ncol = n_exo)
  pmq <- matrix(0, nrow = n_items, ncol = n_exo)
  score_stat <- rep(0, n_exo)
  score_p <- rep(2, n_exo)
  kind <- rep("Gamma", n_exo)
  exo_matrix <- matrix(0L, nrow = nrow(project$raw_data), ncol = n_exo)
  if (n_exo == 0L) {
    return(list(stat = stat, p = p, ppq = ppq, pmq = pmq, score_stat = score_stat, score_p = score_p, kind = kind, exo_values = exo_matrix))
  }
  dimnames(stat) <- dimnames(p) <- dimnames(ppq) <- dimnames(pmq) <- list(items$label_code, backgrounds$label_code)
  names(score_stat) <- names(score_p) <- names(kind) <- backgrounds$label_code

  for (exo_index in seq_len(n_exo)) {
    exo_values <- as.integer(project$raw_data[, backgrounds$position[[exo_index]]])
    exo_matrix[, exo_index] <- exo_values
    exo_dim <- backgrounds$raw_max[[exo_index]]
    exo_valid <- exo_values >= 1L & exo_values <= exo_dim
    use_gamma <- backgrounds$vtype[[exo_index]] > 2L || exo_dim == 2L
    kind[[exo_index]] <- if (use_gamma) "Gamma" else "CHI*2"

    score_tab <- screen_j_pair_table(item_score + 1L, exo_values, max_score + 1L, exo_dim, complete_items & exo_valid)
    if (use_gamma) {
      score_stats <- screen_rc_gamma(score_tab)
      score_stat[[exo_index]] <- score_stats$gamma
      score_p[[exo_index]] <- score_stats$p_value
    } else {
      score_stats <- screen_rc_chi(score_tab)
      score_stat[[exo_index]] <- score_stats$chi_square
      score_p[[exo_index]] <- score_stats$p_value
    }

    for (item_index in seq_len(n_items)) {
      item_values <- item_matrix[, item_index]
      rest_score <- item_score - item_values
      if (!partial) {
        valid_item <- item_values >= 1L & item_values <= items$raw_max[[item_index]]
        tab <- screen_j_pair_table(
          item_values, exo_values, items$raw_max[[item_index]], exo_dim,
          valid_item & exo_valid
        )
        if (use_gamma) {
          stats <- screen_rc_gamma(tab)
          stat[item_index, exo_index] <- stats$gamma
          p[item_index, exo_index] <- stats$p_value
        } else {
          stats <- screen_rc_chi(tab)
          stat[item_index, exo_index] <- stats$chi_square
          p[item_index, exo_index] <- stats$p_value
        }
      } else {
        strata <- screen_j_strata_table(
          item_values, exo_values, item_score + 1L,
          items$raw_max[[item_index]], exo_dim, max_score + 1L,
          complete_items & exo_valid & item_score > 0L & item_score < max_score
        )
        if (use_gamma) {
          stats <- screen_j_partial_gamma(strata)
          stat[item_index, exo_index] <- stats$gamma
          ppq[item_index, exo_index] <- stats$ppq
          pmq[item_index, exo_index] <- stats$pmq
          p[item_index, exo_index] <- if (isTRUE(attr(project, "screen_j_exact"))) {
            chi_stats <- screen_j_partial_chi(strata)
            screen_j_exact_partial_gamma(
              strata,
              stats$gamma,
              attr(project, "screen_j_nsim"),
              attr(project, "screen_j_seed"),
              chi_stats$stat,
              sequential = isTRUE(attr(project, "screen_j_repeated")),
              seq_limit = attr(project, "screen_j_seq_limit"),
              seq_p0 = attr(project, "screen_j_seq_p0"),
              seq_boundary = attr(project, "screen_j_seq_boundary")
            )
          } else {
            stats$p_value
          }
        } else {
          stats <- screen_j_partial_chi(strata)
          stat[item_index, exo_index] <- stats$stat
          p[item_index, exo_index] <- if (isTRUE(attr(project, "screen_j_exact"))) {
            screen_j_exact_partial_chi(
              strata,
              stats$stat,
              attr(project, "screen_j_nsim"),
              attr(project, "screen_j_seed"),
              sequential = isTRUE(attr(project, "screen_j_repeated")),
              seq_limit = attr(project, "screen_j_seq_limit"),
              seq_p0 = attr(project, "screen_j_seq_p0"),
              seq_boundary = attr(project, "screen_j_seq_boundary")
            )
          } else {
            stats$p_value
          }
        }
      }
    }
  }
  list(stat = stat, p = p, ppq = ppq, pmq = pmq, score_stat = score_stat, score_p = score_p, kind = kind, exo_values = exo_matrix)
}

#' Internal screen j weighted partial gamma helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param ppq Internal `ppq` value used by this helper.
#' @param pmq Internal `pmq` value used by this helper.
#' @return The internal `screen_j_weighted_partial_gamma()` computation result.
#' @keywords internal
#' @noRd
screen_j_weighted_partial_gamma <- function(ppq, pmq) {
  n_items <- nrow(ppq)
  weighted <- matrix(0, nrow = n_items, ncol = n_items, dimnames = dimnames(ppq))
  diag(weighted) <- 1
  for (i in seq_len(n_items - 1L)) {
    for (j in seq.int(i + 1L, n_items)) {
      denom <- ppq[i, j] + ppq[j, i]
      value <- if (denom > 0) (pmq[i, j] + pmq[j, i]) / denom else 0
      weighted[i, j] <- value
      weighted[j, i] <- value
    }
  }
  weighted
}

#' Internal screen j average abs partial gamma helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param gamma Internal `gamma` value used by this helper.
#' @return The internal `screen_j_average_abs_partial_gamma()` computation result.
#' @keywords internal
#' @noRd
screen_j_average_abs_partial_gamma <- function(gamma) {
  n_items <- nrow(gamma)
  if (n_items < 2L) {
    return(0)
  }
  total <- 0
  for (i in seq_len(n_items - 1L)) {
    for (j in seq.int(i + 1L, n_items)) {
      total <- total + abs(gamma[i, j] + gamma[j, i])
    }
  }
  total / (n_items * (n_items - 1L))
}

#' Internal screen j problem counts helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param gamma Internal `gamma` value used by this helper.
#' @param p Internal `p` value used by this helper.
#' @param fdr05 Internal `fdr05` value used by this helper.
#' @param i Internal `i` value used by this helper.
#' @param j Internal `j` value used by this helper.
#' @param temp_p Internal `temp_p` value used by this helper.
#' @return The internal `screen_j_problem_counts()` computation result.
#' @keywords internal
#' @noRd
screen_j_problem_counts <- function(gamma, p, fdr05, i, j, temp_p = p) {
  positive <- 0L
  negative <- 0L
  if (temp_p[i, j] <= fdr05) {
    if (gamma[i, j] > 0) positive <- positive + 1L else negative <- negative + 1L
  }
  if (temp_p[j, i] <= fdr05) {
    if (gamma[j, i] > 0) positive <- positive + 1L else negative <- negative + 1L
  }
  c(positive = positive, negative = negative)
}

#' Select and finalize SCREEN J local-dependence evidence
#'
#' Reproduce the four-stage greedy selection in
#' `SKbias7.stepwise_elimination`, then apply the separate command-level
#' positive-local-dependence filter in `DIGRAM1f.pas`. The provisional result
#' is retained for diagnostics, whereas `matrix` and `rows` contain only the
#' terms that DIGRAM permits in the screen model.
#'
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param gamma Directed partial-gamma matrix.
#' @param p Directed partial-gamma p-value matrix.
#' @param weighted_gamma Symmetric weighted-partial-gamma matrix.
#' @param fdr05 Global SCREEN J Benjamini-Hochberg boundary at FDR 0.05.
#' @return A list containing final model `matrix` and `rows`, provisional
#'   `stepwise_matrix` and `stepwise_rows`, and excluded `negative_matrix` and
#'   `negative_rows`.
#' @keywords internal
#' @noRd
screen_j_stepwise_local_dependence <- function(gamma, p, weighted_gamma, fdr05) {
  n_items <- nrow(gamma)
  temp_p <- p
  selected <- matrix(FALSE, nrow = n_items, ncol = n_items, dimnames = dimnames(gamma))
  rows <- list()

  add_selected <- function(i, j, value, stage) {
    selected[i, j] <<- TRUE
    selected[j, i] <<- TRUE
    temp_p[i, ] <<- 2
    temp_p[j, ] <<- 2
    rows[[length(rows) + 1L]] <<- data.frame(
      row = i,
      col = j,
      pair = paste0(rownames(gamma)[[i]], colnames(gamma)[[j]]),
      value = value,
      stage = stage,
      stringsAsFactors = FALSE
    )
  }

  choose_pair <- function(target, sign) {
    best_i <- 0L
    best_j <- 0L
    best_value <- 0
    for (i in seq_len(n_items - 1L)) {
      for (j in seq.int(i + 1L, n_items)) {
        counts <- screen_j_problem_counts(gamma, p, fdr05, i, j, p)
        adjusted <- screen_j_problem_counts(gamma, p, fdr05, i, j, temp_p)
        if (counts[[target]] == sign$count && adjusted[[target]] > 0L) {
          value <- weighted_gamma[i, j]
          if ((sign$direction == "positive" && value > best_value) ||
              (sign$direction == "negative" && value < best_value)) {
            best_i <- i
            best_j <- j
            best_value <- value
          }
        }
      }
    }
    if (best_i > 0L) {
      add_selected(best_i, best_j, best_value, paste0(sign$direction, sign$count))
      TRUE
    } else {
      FALSE
    }
  }

  while (choose_pair("positive", list(count = 2L, direction = "positive"))) {}
  while (choose_pair("positive", list(count = 1L, direction = "positive"))) {}
  while (choose_pair("negative", list(count = 2L, direction = "negative"))) {}
  while (choose_pair("negative", list(count = 1L, direction = "negative"))) {}

  if (length(rows) == 0L) {
    rows_df <- data.frame(
      row = integer(),
      col = integer(),
      pair = character(),
      value = numeric(),
      stage = character(),
      stringsAsFactors = FALSE
    )
  } else {
    rows_df <- do.call(rbind, rows)
  }

  screen_j_finalize_local_dependence(gamma, selected, rows_df)
}

#' Apply DIGRAM's final SCREEN J local-dependence model rule
#'
#' The greedy SCREEN J procedure reports both positive and negative evidence.
#' Since DIGRAM version 3.37, a provisional pair becomes a screen-model term
#' only when the sum of its two directed partial gammas is strictly positive.
#'
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param gamma Directed partial-gamma matrix.
#' @param stepwise_matrix Symmetric logical matrix of provisional selections.
#' @param stepwise_rows Ordered provisional-selection table.
#' @return A local-dependence selection list with separate provisional, final,
#'   and excluded-negative representations.
#' @keywords internal
#' @noRd
screen_j_finalize_local_dependence <- function(gamma, stepwise_matrix, stepwise_rows) {
  rows <- stepwise_rows
  if (!is.data.frame(rows)) {
    rows <- data.frame()
  }

  # Source trace: DIGRAM1f command 3, parameters I/J, removes a provisional
  # LocalDependence[i,j] unless IJXgamma[i,j] + IJXgamma[j,i] > 0. This is an
  # unweighted sum of the two directed partial gammas; it is not the WPG sign
  # and it is not implied by the provisional stage label.
  if (nrow(rows)) {
    index_forward <- cbind(rows$row, rows$col)
    index_reverse <- cbind(rows$col, rows$row)
    rows$directed_gamma_sum <- gamma[index_forward] + gamma[index_reverse]
    rows$included <- !is.na(rows$directed_gamma_sum) & rows$directed_gamma_sum > 0
  } else {
    rows$directed_gamma_sum <- numeric()
    rows$included <- logical()
  }

  final_matrix <- matrix(
    FALSE,
    nrow = nrow(stepwise_matrix),
    ncol = ncol(stepwise_matrix),
    dimnames = dimnames(stepwise_matrix)
  )
  if (nrow(rows) && any(rows$included)) {
    included_rows <- rows[rows$included, , drop = FALSE]
    final_matrix[cbind(included_rows$row, included_rows$col)] <- TRUE
    final_matrix[cbind(included_rows$col, included_rows$row)] <- TRUE
  }

  negative_matrix <- stepwise_matrix & !final_matrix
  included_rows <- rows[rows$included, , drop = FALSE]
  negative_rows <- rows[!rows$included, , drop = FALSE]
  rownames(rows) <- NULL
  rownames(included_rows) <- NULL
  rownames(negative_rows) <- NULL

  list(
    matrix = final_matrix,
    rows = included_rows,
    stepwise_matrix = stepwise_matrix,
    stepwise_rows = rows,
    negative_matrix = negative_matrix,
    negative_rows = negative_rows
  )
}

#' Internal screen j post screen dif helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param item_matrix Internal `item_matrix` value used by this helper.
#' @param exo_values Internal `exo_values` value used by this helper.
#' @param item_score Internal `item_score` value used by this helper.
#' @param complete_items Internal `complete_items` value used by this helper.
#' @param items Item selection or item metadata.
#' @param backgrounds Internal `backgrounds` value used by this helper.
#' @param max_score Internal `max_score` value used by this helper.
#' @param partial_exo Internal `partial_exo` value used by this helper.
#' @param initial_item_bias Internal `initial_item_bias` value used by this helper.
#' @param fdr_01 Internal `fdr_01` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param repeated Whether to use repeated sequential simulation.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_post_screen_dif()` computation result.
#' @keywords internal
#' @noRd
screen_j_post_screen_dif <- function(item_matrix,
                                     exo_values,
                                     item_score,
                                     complete_items,
                                     items,
                                     backgrounds,
                                     max_score,
                                     partial_exo,
                                     initial_item_bias,
                                     fdr_01,
                                     exact = FALSE,
                                     nsim = 1000L,
                                     seed = NULL,
                                     repeated = FALSE,
                                     seq_limit = nsim,
                                     seq_p0 = 0.05,
                                     seq_boundary = 1.058) {
  n_items <- nrow(items)
  n_exo <- nrow(backgrounds)
  item_bias <- initial_item_bias
  item_bias_status <- matrix(
    0L,
    nrow = n_items,
    ncol = n_exo,
    dimnames = dimnames(initial_item_bias)
  )
  item_bias_status[initial_item_bias] <- 1L
  spurious_rows <- list()
  multiple_rows <- list()

  if (n_exo == 0L) {
    return(list(
      item_bias = item_bias,
      item_bias_status = item_bias_status,
      spurious_dif = list(rows = spurious_rows),
      multiple_dif = list(rows = multiple_rows)
    ))
  }

  for (exo_index in seq_len(n_exo)) {
    biased_items <- which(item_bias[, exo_index])
    entry <- NULL

    if (length(biased_items) == n_items) {
      retained <- which(partial_exo$p[, exo_index] <= fdr_01)
      item_bias[, exo_index] <- FALSE
      item_bias[retained, exo_index] <- TRUE
      biased_items <- retained
      entry <- list(
        background = exo_index,
        all_items_warning = TRUE,
        fdr01_reduction = TRUE,
        still_all_warning = length(biased_items) == n_items,
        iterations = list(),
        remaining = biased_items
      )
    }

    if (length(biased_items) > 1L && length(biased_items) < n_items) {
      analysis <- screen_j_spurious_item_bias_analysis(
        item_matrix = item_matrix,
        exo_values = exo_values[, exo_index],
        item_score = item_score,
        complete_items = complete_items,
        items = items,
        background = backgrounds[exo_index, , drop = FALSE],
        exo_index = exo_index,
        biased_items = biased_items,
        max_score = max_score,
        exact = exact,
        nsim = nsim,
        seed = seed,
        repeated = repeated,
        seq_limit = seq_limit,
        seq_p0 = seq_p0,
        seq_boundary = seq_boundary
      )
      if (!is.null(entry)) {
        analysis$all_items_warning <- entry$all_items_warning
        analysis$fdr01_reduction <- entry$fdr01_reduction
        analysis$still_all_warning <- entry$still_all_warning
      }
      spurious_rows[[length(spurious_rows) + 1L]] <- analysis
      item_bias[, exo_index] <- FALSE
      item_bias[analysis$remaining, exo_index] <- TRUE
      for (item_index in analysis$remaining) {
        item_bias_status[item_index, exo_index] <- 2L
      }
    } else {
      if (!is.null(entry)) {
        spurious_rows[[length(spurious_rows) + 1L]] <- entry
      }
      if (length(biased_items) == 1L) {
        item_bias_status[biased_items, exo_index] <- 2L
      }
    }
  }

  for (item_index in seq_len(n_items)) {
    biased_exo <- which(item_bias[item_index, ])
    if (length(biased_exo) > 1L) {
      analysis <- screen_j_multiple_item_bias_analysis(
        item_matrix = item_matrix,
        exo_values = exo_values,
        item_score = item_score,
        complete_items = complete_items,
        items = items,
        backgrounds = backgrounds,
        item_index = item_index,
        biased_exo = biased_exo,
        max_score = max_score,
        exact = exact,
        nsim = nsim,
        seed = seed,
        repeated = repeated,
        seq_limit = seq_limit,
        seq_p0 = seq_p0,
        seq_boundary = seq_boundary
      )
      multiple_rows[[length(multiple_rows) + 1L]] <- analysis
      for (exo_index in biased_exo) {
        if (exo_index %in% analysis$remaining) {
          item_bias_status[item_index, exo_index] <- 2L
        } else if (item_bias[item_index, exo_index]) {
          item_bias_status[item_index, exo_index] <- 1L
        }
      }
      item_bias[item_index, ] <- FALSE
      item_bias[item_index, analysis$remaining] <- TRUE
    } else if (length(biased_exo) == 1L) {
      item_bias_status[item_index, biased_exo] <- 2L
    }
  }

  list(
    item_bias = item_bias,
    item_bias_status = item_bias_status,
    spurious_dif = list(rows = spurious_rows),
    multiple_dif = list(rows = multiple_rows)
  )
}

#' Internal screen j spurious item bias analysis helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param item_matrix Internal `item_matrix` value used by this helper.
#' @param exo_values Internal `exo_values` value used by this helper.
#' @param item_score Internal `item_score` value used by this helper.
#' @param complete_items Internal `complete_items` value used by this helper.
#' @param items Item selection or item metadata.
#' @param background One-based exogenous-variable index.
#' @param exo_index Internal `exo_index` value used by this helper.
#' @param biased_items Internal `biased_items` value used by this helper.
#' @param max_score Internal `max_score` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param repeated Whether to use repeated sequential simulation.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_spurious_item_bias_analysis()` computation result.
#' @keywords internal
#' @noRd
screen_j_spurious_item_bias_analysis <- function(item_matrix,
                                                exo_values,
                                                item_score,
                                                complete_items,
                                                items,
                                                background,
                                                exo_index,
                                                biased_items,
                                                max_score,
                                                exact = FALSE,
                                                nsim = 1000L,
                                                seed = NULL,
                                                repeated = FALSE,
                                                seq_limit = nsim,
                                                seq_p0 = 0.05,
                                                seq_boundary = 1.058) {
  current <- biased_items
  iterations <- list()
  use_gamma <- background$vtype[[1L]] > 2L || background$raw_max[[1L]] == 2L

  repeat {
    rows <- vector("list", length(current))
    p_min <- rep(0, length(current))
    for (candidate_index in seq_along(current)) {
      item_index <- current[[candidate_index]]
      other_items <- setdiff(current, item_index)
      condition_values <- cbind(
        if (length(other_items) > 0L) item_matrix[, other_items, drop = FALSE] else NULL,
        item_score
      )
      condition_dims <- c(items$raw_max[other_items], max_score - 1L)
      test <- screen_j_conditional_bias_test(
        x = exo_values,
        y = item_matrix[, item_index],
        x_dim = background$raw_max[[1L]],
        y_dim = items$raw_max[[item_index]],
        condition_values = condition_values,
        condition_dims = condition_dims,
        valid = complete_items &
          exo_values >= 1L & exo_values <= background$raw_max[[1L]] &
          item_score > 0L & item_score < max_score,
        exact = exact,
        nsim = nsim,
        seed = seed,
        repeated = repeated,
        native = screen_j_conditional_native_allowed(repeated, seq_p0, seq_boundary),
        seq_limit = seq_limit,
        seq_p0 = seq_p0,
        seq_boundary = seq_boundary
      )
      p_min[[candidate_index]] <- screen_j_source_stepwise_p_min(test, use_gamma, exact)
      rows[[candidate_index]] <- data.frame(
        index = candidate_index,
        hypothesis = screen_j_hypothesis_label(
          background$label_code[[1L]],
          items$label_code[[item_index]],
          c(items$label_code[other_items], "#")
        ),
        candidate = item_index,
        candidate_label = items$label_code[[item_index]],
        candidate_name = items$name[[item_index]],
        chi_square = test$chi_square,
        df = test$df,
        p_chi = test$p_chi,
        p_chi_asymp = test$p_chi_asymp,
        p_chi_exact = test$p_chi_exact,
        gamma = test$gamma,
        p_gamma = test$p_gamma,
        p_gamma_asymp = test$p_gamma_asymp,
        p_gamma_exact = test$p_gamma_exact,
        p_min = p_min[[candidate_index]],
        use_gamma = use_gamma,
        n = test$exact_nsim,
        stringsAsFactors = FALSE
      )
    }

    rows_df <- do.call(rbind, rows)
    rows_df$bh_05 <- source_bh_critical(screen_j_exa_p_values(rows_df), 0.05)
    rows_df$bh_01 <- source_bh_critical(screen_j_exa_p_values(rows_df), 0.01)
    excluded <- integer()
    excluded_label <- NA_character_
    excluded_name <- NA_character_
    if (length(p_min) > 0L && max(p_min) > 0.05) {
      excluded_index <- which.max(p_min)
      excluded <- current[[excluded_index]]
      excluded_label <- items$label_code[[excluded]]
      excluded_name <- items$name[[excluded]]
    }
    iterations[[length(iterations) + 1L]] <- list(
      rows = rows_df,
      excluded = excluded,
      excluded_label = excluded_label,
      excluded_name = excluded_name
    )
    if (length(excluded) == 0L) {
      break
    }
    current <- current[current != excluded]
    if (length(current) == 0L) {
      break
    }
  }

  list(
    background = exo_index,
    all_items_warning = FALSE,
    fdr01_reduction = FALSE,
    still_all_warning = FALSE,
    iterations = iterations,
    remaining = current
  )
}

#' Internal screen j multiple item bias analysis helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param item_matrix Internal `item_matrix` value used by this helper.
#' @param exo_values Internal `exo_values` value used by this helper.
#' @param item_score Internal `item_score` value used by this helper.
#' @param complete_items Internal `complete_items` value used by this helper.
#' @param items Item selection or item metadata.
#' @param backgrounds Internal `backgrounds` value used by this helper.
#' @param item_index One-based item index.
#' @param biased_exo Internal `biased_exo` value used by this helper.
#' @param max_score Internal `max_score` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param repeated Whether to use repeated sequential simulation.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_multiple_item_bias_analysis()` computation result.
#' @keywords internal
#' @noRd
screen_j_multiple_item_bias_analysis <- function(item_matrix,
                                                exo_values,
                                                item_score,
                                                complete_items,
                                                items,
                                                backgrounds,
                                                item_index,
                                                biased_exo,
                                                max_score,
                                                exact = FALSE,
                                                nsim = 1000L,
                                                seed = NULL,
                                                repeated = FALSE,
                                                seq_limit = nsim,
                                                seq_p0 = 0.05,
                                                seq_boundary = 1.058) {
  current <- biased_exo
  iterations <- list()

  repeat {
    rows <- vector("list", length(current))
    p_min <- rep(0, length(current))
    for (candidate_index in seq_along(current)) {
      exo_index <- current[[candidate_index]]
      other_exo <- setdiff(current, exo_index)
      condition_values <- cbind(
        if (length(other_exo) > 0L) exo_values[, other_exo, drop = FALSE] else NULL,
        item_score
      )
      condition_dims <- c(backgrounds$raw_max[other_exo], max_score - 1L)
      valid <- complete_items & item_score > 0L & item_score < max_score
      for (valid_exo in current) {
        valid <- valid &
          exo_values[, valid_exo] >= 1L &
          exo_values[, valid_exo] <= backgrounds$raw_max[[valid_exo]]
      }
      test <- screen_j_conditional_bias_test(
        x = item_matrix[, item_index],
        y = exo_values[, exo_index],
        x_dim = items$raw_max[[item_index]],
        y_dim = backgrounds$raw_max[[exo_index]],
        condition_values = condition_values,
        condition_dims = condition_dims,
        valid = valid,
        exact = exact,
        nsim = nsim,
        seed = seed,
        repeated = repeated,
        native = screen_j_conditional_native_allowed(repeated, seq_p0, seq_boundary),
        seq_limit = seq_limit,
        seq_p0 = seq_p0,
        seq_boundary = seq_boundary
      )
      use_gamma <- backgrounds$vtype[[exo_index]] > 2L || backgrounds$raw_max[[exo_index]] == 2L
      p_min[[candidate_index]] <- screen_j_source_stepwise_p_min(test, use_gamma, exact)
      rows[[candidate_index]] <- data.frame(
        index = candidate_index,
        hypothesis = screen_j_hypothesis_label(
          items$label_code[[item_index]],
          backgrounds$label_code[[exo_index]],
          c(backgrounds$label_code[other_exo], "#")
        ),
        candidate = exo_index,
        candidate_label = backgrounds$label_code[[exo_index]],
        candidate_name = backgrounds$name[[exo_index]],
        chi_square = test$chi_square,
        df = test$df,
        p_chi = test$p_chi,
        p_chi_asymp = test$p_chi_asymp,
        p_chi_exact = test$p_chi_exact,
        gamma = test$gamma,
        p_gamma = test$p_gamma,
        p_gamma_asymp = test$p_gamma_asymp,
        p_gamma_exact = test$p_gamma_exact,
        p_min = p_min[[candidate_index]],
        use_gamma = use_gamma,
        n = test$exact_nsim,
        stringsAsFactors = FALSE
      )
    }

    rows_df <- do.call(rbind, rows)
    rows_df$bh_05 <- source_bh_critical(screen_j_exa_p_values(rows_df), 0.05)
    rows_df$bh_01 <- source_bh_critical(screen_j_exa_p_values(rows_df), 0.01)
    excluded <- integer()
    excluded_label <- NA_character_
    excluded_name <- NA_character_
    if (length(p_min) > 0L && max(p_min) > 0.05) {
      excluded_index <- which.max(p_min)
      excluded <- current[[excluded_index]]
      excluded_label <- backgrounds$label_code[[excluded]]
      excluded_name <- backgrounds$name[[excluded]]
    }
    iterations[[length(iterations) + 1L]] <- list(
      rows = rows_df,
      excluded = excluded,
      excluded_label = excluded_label,
      excluded_name = excluded_name
    )
    if (length(excluded) == 0L) {
      break
    }
    current <- current[current != excluded]
    if (length(current) == 0L) {
      break
    }
  }

  list(item = item_index, iterations = iterations, remaining = current)
}
