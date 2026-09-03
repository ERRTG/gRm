#' Source trace: source/digram_source_20260817/skunits/SKbias13.pas::StepwiseScoreScreening builds
#' score-effect screening rows from the same conditional bias-test machinery
#' used by screen J. The R helper returns those rows as numeric values rather
#' than writing DIGRAM text.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param project Encoded gRm project.
#' @param item_matrix Internal `item_matrix` value used by this helper.
#' @param item_score Internal `item_score` value used by this helper.
#' @param complete_items Internal `complete_items` value used by this helper.
#' @param backgrounds Internal `backgrounds` value used by this helper.
#' @param max_score Internal `max_score` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param repeated Whether to use repeated sequential simulation.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_score_effect_values()` computation result.
#' @keywords internal
#' @noRd
screen_j_score_effect_values <- function(project,
                                         item_matrix,
                                         item_score,
                                         complete_items,
                                         backgrounds,
                                         max_score,
                                         exact = FALSE,
                                         nsim = 1000L,
                                         seed = NULL,
                                         repeated = FALSE,
                                         seq_limit = nsim,
                                         seq_p0 = 0.05,
                                         seq_boundary = 1.058) {
  n_exo <- nrow(backgrounds)
  maxdim <- 57L
  score_top <- maxdim - 1L
  score_category <- pmin(item_score, score_top) + 1L
  score_valid <- complete_items
  rows <- vector("list", n_exo)
  if (n_exo == 0L) {
    return(list(
      rows = data.frame(),
      fdr_05 = 0,
      fdr_01 = 0,
      selected = character(),
      score_top = score_top,
      exact = exact,
      screening = list(rows = data.frame(), iterations = list())
    ))
  }
  p_values <- numeric(0)
  nsim <- as.integer(nsim)
  marginal_p_min <- rep(1, n_exo)
  exo_matrix <- matrix(0L, nrow = nrow(project$raw_data), ncol = n_exo)
  for (exo_index in seq_len(n_exo)) {
    exo_values <- as.integer(project$raw_data[, backgrounds$position[[exo_index]]])
    exo_matrix[, exo_index] <- exo_values
    exo_dim <- backgrounds$raw_max[[exo_index]]
    exo_valid <- exo_values >= 1L & exo_values <= exo_dim
    test <- screen_j_conditional_bias_test(
      x = score_category,
      y = exo_values,
      x_dim = score_top + 1L,
      y_dim = exo_dim,
      condition_values = matrix(integer(), nrow = nrow(item_matrix), ncol = 0L),
      condition_dims = integer(),
      valid = score_valid & exo_valid,
      exact = exact,
      nsim = nsim,
      seed = seed,
      repeated = repeated,
      native = screen_j_conditional_native_allowed(repeated, seq_p0, seq_boundary),
      seq_limit = seq_limit,
      seq_p0 = seq_p0,
      seq_boundary = seq_boundary
    )
    p_chi <- test$p_chi
    p_gamma <- test$p_gamma
    p_chi_exact <- test$p_chi_exact
    p_gamma_exact <- test$p_gamma_exact
    exact_nsim <- test$exact_nsim
    marginal_p_min[[exo_index]] <- min(p_chi, p_gamma)
    p_values <- c(p_values, p_chi, p_gamma)
    rows[[exo_index]] <- data.frame(
      label = backgrounds$label_code[[exo_index]],
      name = backgrounds$name[[exo_index]],
      chi_square = test$chi_square,
      df = test$df,
      p_chi = p_chi,
      p_chi_asymp = test$p_chi_asymp,
      p_chi_exact = p_chi_exact,
      p_chi_exact_low = if (isTRUE(exact)) screen_j_conflimit99_field(exact_nsim, p_chi_exact, "low") else NA_real_,
      p_chi_exact_high = if (isTRUE(exact)) screen_j_conflimit99_field(exact_nsim, p_chi_exact, "high") else NA_real_,
      gamma = test$gamma,
      p_gamma = p_gamma,
      p_gamma_asymp = test$p_gamma_asymp,
      p_gamma_exact = p_gamma_exact,
      p_gamma_exact_low = if (isTRUE(exact)) screen_j_conflimit99_field(exact_nsim, p_gamma_exact, "low") else NA_real_,
      p_gamma_exact_high = if (isTRUE(exact)) screen_j_conflimit99_field(exact_nsim, p_gamma_exact, "high") else NA_real_,
      exact_nsim = exact_nsim,
      selected = FALSE,
      stringsAsFactors = FALSE
    )
  }
  rows_df <- do.call(rbind, rows)
  screening <- screen_j_stepwise_score_screening(
    score_category = score_category,
    score_valid = score_valid,
    score_in_range = item_score <= score_top,
    exo_values = exo_matrix,
    backgrounds = backgrounds,
    marginal_p_min = marginal_p_min,
    score_top = score_top,
    exact = exact,
    nsim = nsim,
    seed = seed,
    repeated = repeated,
    seq_limit = seq_limit,
    seq_p0 = seq_p0,
    seq_boundary = seq_boundary
  )
  rows_df$selected <- seq_len(n_exo) %in% screening$remaining
  report_rows <- screen_j_score_report_rows(rows_df, screening$rows, exact)
  p_values <- c(p_values, screen_j_score_screening_p_values(screening$rows))
  list(
    rows = report_rows,
    fdr_05 = source_bh_critical(p_values, 0.05),
    fdr_01 = source_bh_critical(p_values, 0.01),
    marginal_selected = rows_df$label[marginal_p_min <= 0.05],
    selected = rows_df$label[rows_df$selected],
    score_top = score_top,
    exact = exact,
    screening = screening
  )
}

#' Source trace: source/digram_source_20260817/skunits/SKbias13.pas::StepwiseScoreScreening runs
#' the marginal and conditional score-screening loop. The R function keeps the
#' same candidate order and source random-table convention when exact or
#' repeated inference is requested.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param score_category Internal `score_category` value used by this helper.
#' @param score_valid Internal `score_valid` value used by this helper.
#' @param score_in_range Internal `score_in_range` value used by this helper.
#' @param exo_values Internal `exo_values` value used by this helper.
#' @param backgrounds Internal `backgrounds` value used by this helper.
#' @param marginal_p_min Internal `marginal_p_min` value used by this helper.
#' @param score_top Internal `score_top` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param repeated Whether to use repeated sequential simulation.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_stepwise_score_screening()` computation result.
#' @keywords internal
#' @noRd
screen_j_stepwise_score_screening <- function(score_category,
                                              score_valid,
                                              score_in_range,
                                              exo_values,
                                              backgrounds,
                                              marginal_p_min,
                                              score_top,
                                              exact = FALSE,
                                              nsim = 1000L,
                                              seed = NULL,
                                              repeated = FALSE,
                                              seq_limit = nsim,
                                              seq_p0 = 0.05,
                                              seq_boundary = 1.058) {
  n_exo <- nrow(backgrounds)
  current <- which(marginal_p_min <= 0.05)
  iterations <- list()
  if (length(current) <= 1L) {
    return(list(rows = data.frame(), iterations = iterations, remaining = current))
  }

  while (length(current) > 1L) {
    rows <- vector("list", length(current))
    p_min <- rep(0, length(current))
    for (candidate_index in seq_along(current)) {
      exo_index <- current[[candidate_index]]
      other_exo <- setdiff(current, exo_index)
      condition_values <- if (length(other_exo) > 0L) {
        exo_values[, other_exo, drop = FALSE]
      } else {
        matrix(integer(), nrow = nrow(exo_values), ncol = 0L)
      }
      valid <- score_valid & score_in_range
      for (valid_exo in current) {
        valid <- valid &
          exo_values[, valid_exo] >= 1L &
          exo_values[, valid_exo] <= backgrounds$raw_max[[valid_exo]]
      }
      test <- screen_j_conditional_bias_test(
        x = score_category,
        y = exo_values[, exo_index],
        x_dim = score_top + 1L,
        y_dim = backgrounds$raw_max[[exo_index]],
        condition_values = condition_values,
        condition_dims = backgrounds$raw_max[other_exo],
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
      p_min[[candidate_index]] <- screen_j_source_stepwise_p_min(test, use_gamma = TRUE, exact = exact)
      rows[[candidate_index]] <- data.frame(
        index = candidate_index,
        hypothesis = screen_j_hypothesis_label(
          "#",
          backgrounds$label_code[[exo_index]],
          backgrounds$label_code[other_exo]
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
        n = test$exact_nsim,
        stringsAsFactors = FALSE
      )
    }

    rows_df <- do.call(rbind, rows)
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
    if (length(current) <= 1L) {
      break
    }
  }

  rows_df <- if (length(iterations) == 0L) data.frame() else do.call(rbind, lapply(iterations, `[[`, "rows"))
  list(rows = rows_df, iterations = iterations, remaining = current)
}

#' Internal screen j score screening p values helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param rows Rows used by the computation.
#' @return The internal `screen_j_score_screening_p_values()` computation result.
#' @keywords internal
#' @noRd
screen_j_score_screening_p_values <- function(rows) {
  if (is.null(rows) || nrow(rows) == 0L) {
    return(numeric(0))
  }
  c(rows$p_chi, rows$p_gamma)
}

#' Internal screen j score report rows helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param marginal_rows Internal `marginal_rows` value used by this helper.
#' @param screening_rows Internal `screening_rows` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @return The internal `screen_j_score_report_rows()` computation result.
#' @keywords internal
#' @noRd
screen_j_score_report_rows <- function(marginal_rows, screening_rows, exact) {
  marginal_rows$hypothesis <- paste0("#&", marginal_rows$label)
  if (is.null(screening_rows) || nrow(screening_rows) == 0L) {
    return(marginal_rows)
  }
  screening <- data.frame(
    label = screening_rows$candidate_label,
    name = screening_rows$candidate_name,
    chi_square = screening_rows$chi_square,
    df = screening_rows$df,
    p_chi = screening_rows$p_chi,
    p_chi_asymp = screening_rows$p_chi_asymp,
    p_chi_exact = screening_rows$p_chi_exact,
    p_chi_exact_low = if (isTRUE(exact)) screen_j_conflimit99_field(screening_rows$n, screening_rows$p_chi_exact, "low") else NA_real_,
    p_chi_exact_high = if (isTRUE(exact)) screen_j_conflimit99_field(screening_rows$n, screening_rows$p_chi_exact, "high") else NA_real_,
    gamma = screening_rows$gamma,
    p_gamma = screening_rows$p_gamma,
    p_gamma_asymp = screening_rows$p_gamma_asymp,
    p_gamma_exact = screening_rows$p_gamma_exact,
    p_gamma_exact_low = if (isTRUE(exact)) screen_j_conflimit99_field(screening_rows$n, screening_rows$p_gamma_exact, "low") else NA_real_,
    p_gamma_exact_high = if (isTRUE(exact)) screen_j_conflimit99_field(screening_rows$n, screening_rows$p_gamma_exact, "high") else NA_real_,
    exact_nsim = screening_rows$n,
    selected = FALSE,
    hypothesis = screening_rows$hypothesis,
    stringsAsFactors = FALSE
  )
  rbind(marginal_rows, screening[, names(marginal_rows)])
}

#' Internal screen j conflimit99 field helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param n Internal `n` value used by this helper.
#' @param p Internal `p` value used by this helper.
#' @param field Internal `field` value used by this helper.
#' @return The internal `screen_j_conflimit99_field()` computation result.
#' @keywords internal
#' @noRd
screen_j_conflimit99_field <- function(n, p, field) {
  n <- rep_len(as.integer(n), length(p))
  p <- as.numeric(p)
  vapply(seq_along(p), function(i) {
    source_conflimit99(n[[i]], p[[i]])[[field]]
  }, numeric(1L))
}
