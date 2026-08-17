#' Internal gRm source score cap helper
#'
#' Supports the exo select values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @return The internal `gRm_source_score_cap()` computation result.
#' @keywords internal
#' @noRd
gRm_source_score_cap <- function() {
  # Source trace: source/PAS_skunits/SKTypes.pas defines MAXDIM = 57, and
  # source/PAS_skunits/SKbias13.pas collapses scores above maxdim - 1.
  56L
}

#' Derive DIGRAM exogenous-selection numeric values
#'
#' Computes source-shaped exogenous overview, missing-exogenous diagnostics,
#' score distribution, and score/exogenous screening values. The source path is
#' `SKbias2.pas` for exogenous overview, missing-exogenous diagnostics and
#' score distribution, then `SKbias13.pas::StepwiseScoreScreening` for the
#' score/exogenous screening.
#'
#' Production R computes from the encoded project object and its `raw_data`
#' matrix; Pascal and supplied DIGRAM output files are external validation
#' inputs only.
#'
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param project A parsed DIGRAM project from [read_digram_project()].
#' @param score_cap Highest printed score category for the score/exogenous
#'   screening table. Source `maxdim - 1` gives `56` for this source tree.
#' @param exact Logical; when `TRUE`, run the source-shaped Monte Carlo exact
#'   table branch used by `SKexa2`/`SKrandom.GENTAB1`.
#' @param nsim Number of random tables for exact p-values.
#' @param seed Seed for the R implementation of the source Monte Carlo branch.
#' @return A `gRm_exo_select_values` object.
#' @examples
#' \dontrun{
#' project <- read_digram_project("path/to/DIGRAM")
#' values <- exo_select_values(project)
#' values$screen
#' }
#' @param repeated Whether to use repeated sequential simulation.
#' @param exact_state Internal `exact_state` value used by this helper.
#' @keywords internal
#' @noRd
exo_select_values <- function(project,
                              score_cap = gRm_source_score_cap(),
                              exact = FALSE,
                              repeated = FALSE,
                              nsim = 1000L,
                              seed = 9L,
                              exact_state = NULL) {
  exact <- isTRUE(exact)
  repeated <- isTRUE(repeated)
  nsim <- as.integer(nsim)
  seed <- as.integer(seed)
  if (is.null(exact_state)) {
    exact_state <- gRm_exact_state_from_flags(exact, repeated, nsim = nsim, seed = seed)
  }
  exact <- isTRUE(exact_state$exact)
  repeated <- isTRUE(exact_state$sequential)
  nsim <- as.integer(exact_state$nsim)
  seed <- as.integer(exact_state$seed)
  items <- project$items
  backgrounds <- project$backgrounds
  raw <- project$raw_data
  item_matrix <- raw[, items$position, drop = FALSE]
  max_values <- items$raw_max
  valid_item <- sweep(item_matrix, 2L, max_values, `<=`) & item_matrix >= 1L
  complete_item <- rowSums(valid_item) == ncol(item_matrix)

  # Source trace: DIGRAM item responses are stored as 1..item_dim, while the
  # raw score used by SKbias2/SKbias13 subtracts one from each item response.
  scores <- rowSums(item_matrix - 1L)

  exo_matrix <- raw[, backgrounds$position, drop = FALSE]
  valid_exo <- sweep(exo_matrix, 2L, backgrounds$raw_max, `<=`) & exo_matrix >= 1L
  complete_item_exo <- complete_item & rowSums(valid_exo) == ncol(exo_matrix)

  score_distribution <- exo_select_score_distribution(
    scores,
    known = complete_item_exo,
    missing_known = rep(TRUE, length(scores)),
    min_score = 0L,
    max_score = sum(items$raw_max - 1L),
    largest_possible_score = sum(items$raw_max - 1L)
  )
  missing <- exo_select_missing_diagnostics(project, scores, complete_item, valid_exo)
  screen <- exo_select_score_screen(
    project,
    scores,
    complete_item,
    score_cap,
    exact = exact,
    repeated = repeated,
    nsim = nsim,
    seed = seed,
    exact_state = exact_state
  )
  marginal_selected <- screen$exo_label[screen$marginal_selected]
  selected <- screen$exo_label[screen$selected]

  structure(
    list(
      project = project,
      exogenous = data.frame(
        exo_label = backgrounds$label_code,
        exo_name = backgrounds$name,
        categories = backgrounds$raw_max,
        recursive_level = 1L,
        stringsAsFactors = FALSE
      ),
      complete_item_cases = sum(complete_item),
      complete_item_exo_cases = sum(complete_item_exo),
      missing = missing,
      score_distribution = score_distribution$distribution,
      score_summary = score_distribution$summary,
      recursive_line = paste0(paste(items$label_code, collapse = ""), paste(backgrounds$label_code, collapse = ""), "# <- \u00a4"),
      score_cap = as.integer(score_cap),
      exact = isTRUE(exact),
      repeated = isTRUE(repeated),
      nsim = if (isTRUE(exact)) as.integer(nsim) else 0L,
      seed = if (isTRUE(exact)) as.integer(seed) else NA_integer_,
      exact_state = exact_state,
      screen = screen,
      # Source trace: SKbias13.StepwiseScoreScreening calls
      # SKexa1.EXA_SUMMARY1_2 after Fillchar(Hypotese, ..., 3), so each
      # marginal score/exogenous hypothesis contributes both the chi-square
      # p-value and the two-sided gamma p-value. SKexa1 then calls the source
      # BenjaminiHochberg routine over the actual p-value vector, not a fixed
      # Bonferroni divisor; the resulting critical level depends on which
      # ordered p-values cross the BH line.
      bh = exo_select_bh_thresholds(screen),
      marginal_selected = marginal_selected,
      selected = selected
    ),
    class = "gRm_exo_select_values"
  )
}

#' Missing-exogenous diagnostics for complete item cases
#'
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param project Parsed project.
#' @param scores Raw scores for all records.
#' @param complete_item Complete-item mask.
#' @param valid_exo Valid exogenous-value matrix.
#' @return Data frame of source-shaped missing diagnostics.
#' @keywords internal
#' @noRd
exo_select_missing_diagnostics <- function(project, scores, complete_item, valid_exo) {
  backgrounds <- project$backgrounds
  rows <- vector("list", nrow(backgrounds))
  for (exo_index in seq_len(nrow(backgrounds))) {
    missing_mask <- complete_item & !valid_exo[, exo_index]
    known_mask <- complete_item & valid_exo[, exo_index]
    missing_scores <- scores[missing_mask]
    known_scores <- scores[known_mask]
    if (length(missing_scores) == 0L || length(known_scores) == 0L) {
      rows[[exo_index]] <- data.frame(
        exo_label = backgrounds$label_code[[exo_index]],
        exo_name = backgrounds$name[[exo_index]],
        count = length(missing_scores),
        missing_mean = NA_real_,
        known_mean = NA_real_,
        t_stat = NA_real_,
        p_value = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    missing_mean <- mean(missing_scores)
    known_mean <- mean(known_scores)
    # Source trace: SKbias2.pas divides the within-group population variance by
    # group n before adding the two standard errors in the denominator.
    sd_missing <- sqrt((mean(missing_scores^2) - missing_mean^2) / length(missing_scores))
    sd_known <- sqrt((mean(known_scores^2) - known_mean^2) / length(known_scores))
    t_stat <- (-known_mean + missing_mean) / (sd_missing + sd_known)
    rows[[exo_index]] <- data.frame(
      exo_label = backgrounds$label_code[[exo_index]],
      exo_name = backgrounds$name[[exo_index]],
      count = length(missing_scores),
      missing_mean = missing_mean,
      known_mean = known_mean,
      t_stat = t_stat,
      p_value = 2 * source_tail_norm(abs(t_stat), TRUE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

#' Score distribution for exo-selection item-score cases
#'
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param scores Raw scores for all records.
#' @param known Logical mask for records with source-known item scores.
#' @param min_score Lower score limit.
#' @param max_score Upper score limit.
#' @param largest_possible_score Largest source-possible score.
#' @return Distribution and summary list.
#' @param missing_known Internal `missing_known` value used by this helper.
#' @keywords internal
#' @noRd
exo_select_score_distribution <- function(scores,
                                          known = rep(TRUE, length(scores)),
                                          missing_known = known,
                                          min_score = 0L,
                                          max_score = NULL,
                                          largest_possible_score = NULL) {
  known <- as.logical(known)
  known[is.na(known)] <- FALSE
  missing_known <- as.logical(missing_known)
  missing_known[is.na(missing_known)] <- FALSE
  if (length(known) != length(scores)) {
    stop("`known` must have the same length as `scores`.", call. = FALSE)
  }
  if (length(missing_known) != length(scores)) {
    stop("`missing_known` must have the same length as `scores`.", call. = FALSE)
  }

  min_score <- as.integer(min_score)
  if (is.null(max_score)) {
    max_score <- if (any(known)) max(scores[known]) else min_score
  }
  if (is.null(largest_possible_score)) {
    largest_possible_score <- max_score
  }
  max_score <- as.integer(max_score)
  largest_possible_score <- as.integer(largest_possible_score)
  known_scores <- as.integer(scores[known])
  eligible_scores <- known_scores[known_scores >= min_score & known_scores <= max_score]
  n <- length(eligible_scores)
  missing <- sum(!missing_known)
  below <- sum(known_scores < min_score)
  above <- sum(known_scores > max_score & known_scores <= largest_possible_score)

  # Source trace: SKbias2.pas::SHOW_SCOREDISTRIBUTION prints a dedicated
  # "No cases with known scores" branch when the eligible score count is zero.
  if (n == 0L) {
    return(list(
      distribution = data.frame(
        score = integer(),
        count = integer(),
        percent = numeric(),
        cumulative = numeric(),
        stringsAsFactors = FALSE
      ),
      summary = list(
        n = 0L,
        mean = NA_real_,
        variance = NA_real_,
        sd = NA_real_,
        skewness = NA_real_,
        below = as.integer(below),
        above = as.integer(above),
        missing = as.integer(missing)
      )
    ))
  }

  highest_score <- max(eligible_scores)
  score_values <- seq.int(min_score, highest_score)
  score_factor <- factor(eligible_scores, levels = score_values)
  counts <- as.integer(tabulate(score_factor, nbins = length(score_values)))
  mean_score <- sum(score_values * counts) / n
  # Source trace: SKbias2.pas::SHOW_SCOREDISTRIBUTION prints the sample
  # variance after accumulating score sums over the displayed distribution.
  variance <- if (n > 1L) {
    (sum(score_values^2 * counts) / n - mean_score^2) * n / (n - 1L)
  } else {
    NA_real_
  }
  sd_score <- sqrt(variance)
  skewness <- if (n > 2L && is.finite(sd_score) && sd_score > 0) {
    third_moment <- sum((score_values - mean_score)^3 * counts / n)
    third_moment * n^2 / ((n - 1L) * (n - 2L)) / sd_score^3
  } else {
    NA_real_
  }
  percent <- 100 * counts / n
  distribution <- data.frame(
    score = score_values,
    count = counts,
    percent = percent,
    cumulative = cumsum(percent),
    stringsAsFactors = FALSE
  )
  list(
    distribution = distribution,
    summary = list(
      n = n,
      mean = mean_score,
      variance = variance,
      sd = sd_score,
      skewness = skewness,
      below = as.integer(below),
      above = as.integer(above),
      missing = as.integer(missing)
    )
  )
}

#' Source score/exogenous screening values
#'
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param project Parsed project.
#' @param scores Raw scores.
#' @param complete_item Complete-item mask.
#' @param score_cap Maximum score category after source collapsing.
#' @return Data frame of screening rows.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param repeated Whether to use repeated sequential simulation.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param exact_state Internal `exact_state` value used by this helper.
#' @keywords internal
#' @noRd
exo_select_score_screen <- function(project,
                                    scores,
                                    complete_item,
                                    score_cap = gRm_source_score_cap(),
                                    exact = FALSE,
                                    repeated = FALSE,
                                    nsim = 1000L,
                                    seed = 9L,
                                    exact_state = NULL) {
  if (is.null(exact_state)) {
    exact_state <- gRm_exact_state_from_flags(exact, repeated, nsim = nsim, seed = seed)
  }
  exact <- isTRUE(exact_state$exact)
  repeated <- isTRUE(exact_state$sequential)
  nsim <- as.integer(exact_state$nsim)
  seed <- as.integer(exact_state$seed)
  backgrounds <- project$backgrounds
  raw <- project$raw_data
  rows <- vector("list", nrow(backgrounds))
  score_category <- pmin(scores, score_cap) + 1L
  exo_matrix <- matrix(0L, nrow = nrow(raw), ncol = nrow(backgrounds))
  marginal_p_min <- rep(1, nrow(backgrounds))
  for (exo_index in seq_len(nrow(backgrounds))) {
    exo_values <- raw[, backgrounds$position[[exo_index]]]
    exo_matrix[, exo_index] <- exo_values
    valid <- complete_item & exo_values >= 1L & exo_values <= backgrounds$raw_max[[exo_index]]
    collapsed <- pmin(scores[valid], score_cap)
    tab <- table(
      factor(collapsed, levels = seq.int(0L, score_cap)),
      factor(exo_values[valid], levels = seq_len(backgrounds$raw_max[[exo_index]]))
    )
    source_tab <- tab
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    test <- screen_rc_chi(tab)
    gamma <- exo_select_gamma_test(tab)
    chi_p <- test$p_value
    gamma_p_two_sided <- 2 * gamma$p_one_sided
    exact_result <- if (isTRUE(exact)) {
      # Source trace: SKbias13.StepwiseScoreScreening marginal screening calls
      # Inexpensive_bt_tests once per exogenous variable with NHYP=1; SKexa1
      # INIT_EXACT_TEST runs inside each call and SKrandom.GENTAB1 then draws
      # from Pascal Random. Reset the source stream for each marginal table.
      exo_select_exact_test(
        source_tab,
        test$chi_square,
        gamma$gamma,
        nsim = nsim,
        random_draw = screen_j_source_random_stream(seed),
        seed = seed,
        sequential = repeated,
        sequential_count_cutoff = exo_select_repeated_count_cutoff(project, exact_state = exact_state),
        seq_p0 = exact_state$seq_p0,
        seq_boundary = exact_state$seq_b
      )
    } else {
        list(chi_p = NA_real_, gamma_p_one_sided = NA_real_, gamma_p_two_sided = NA_real_, nsim = 0L)
    }
    selected_chi_p <- if (isTRUE(exact)) exact_result$chi_p else chi_p
    selected_gamma_p <- if (isTRUE(exact)) exact_result$gamma_p_two_sided else gamma_p_two_sided
    marginal_p_min[[exo_index]] <- min(selected_chi_p, selected_gamma_p)
    rows[[exo_index]] <- data.frame(
      exo_label = backgrounds$label_code[[exo_index]],
      exo_name = backgrounds$name[[exo_index]],
      hypothesis = paste0("#&", backgrounds$label_code[[exo_index]]),
      chi_square = test$chi_square,
      df = test$df,
      chi_p = chi_p,
      gamma = gamma$gamma,
      gamma_p_one_sided = gamma$p_one_sided,
      gamma_p_two_sided = gamma_p_two_sided,
      exact_chi_p = exact_result$chi_p,
      exact_gamma_p_one_sided = exact_result$gamma_p_one_sided,
      exact_gamma_p_two_sided = exact_result$gamma_p_two_sided,
      exact_nsim = exact_result$nsim,
      selected = FALSE,
      chi_marker = "",
      gamma_marker = "",
      stringsAsFactors = FALSE
    )
  }
  screen <- do.call(rbind, rows)
  screening <- exo_select_stepwise_score_screening(
    score_category = score_category,
    score_valid = complete_item,
    score_in_range = scores <= score_cap,
    exo_values = exo_matrix,
    backgrounds = backgrounds,
    marginal_p_min = marginal_p_min,
    score_cap = score_cap,
    exact = exact,
    repeated = repeated,
    nsim = nsim,
    seed = seed,
    exact_state = exact_state
  )
  if (nrow(screening$rows) > 0L) {
    screen <- rbind(screen, screening$rows[, names(screen)])
  }
  bh <- exo_select_bh_thresholds(screen)
  chi_marker_p <- if (isTRUE(exact)) screen$exact_chi_p else screen$chi_p
  gamma_marker_p <- if (isTRUE(exact)) screen$exact_gamma_p_two_sided else screen$gamma_p_one_sided
  selected_chi_p <- if (isTRUE(exact)) screen$exact_chi_p else screen$chi_p
  selected_gamma_p <- if (isTRUE(exact)) screen$exact_gamma_p_two_sided else screen$gamma_p_two_sided
  screen$selected <- FALSE
  screen$marginal_selected <- FALSE
  screen$marginal_selected[seq_len(nrow(backgrounds))] <- marginal_p_min <= 0.05
  screen$selected[seq_len(nrow(backgrounds))] <- seq_len(nrow(backgrounds)) %in% screening$remaining
  screen$chi_marker <- mapply(
    exo_select_marker,
    p_value = chi_marker_p,
    statistic = screen$chi_p,
    MoreArgs = list(bh05 = bh$fdr_05, bh01 = bh$fdr_01, positive = TRUE, chi = TRUE),
    USE.NAMES = FALSE
  )
  screen$gamma_marker <- mapply(
    exo_select_marker,
    # Source/runtime trace: the supplied DIGRAM score-screening reports print
    # the EXA_SUMMARY1_2 two-sided gamma p-value, but the significance marker
    # matches the undoubled Results[*,6] value used elsewhere in the SKexa1
    # EXA summary family. example's first score/exo row is the boundary case:
    # printed two-sided p = 0.002 and printed FDR 0.01 threshold = 0.002, while
    # the runtime marker is "++".
    p_value = gamma_marker_p,
    statistic = screen$gamma,
    positive = screen$gamma >= 0,
    MoreArgs = list(bh05 = bh$fdr_05, bh01 = bh$fdr_01, chi = FALSE),
    USE.NAMES = FALSE
  )
  screen
}

#' Internal exo select stepwise score screening helper
#'
#' Supports the exo select values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param score_category Internal `score_category` value used by this helper.
#' @param score_valid Internal `score_valid` value used by this helper.
#' @param score_in_range Internal `score_in_range` value used by this helper.
#' @param exo_values Internal `exo_values` value used by this helper.
#' @param backgrounds Internal `backgrounds` value used by this helper.
#' @param marginal_p_min Internal `marginal_p_min` value used by this helper.
#' @param score_cap Internal `score_cap` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param repeated Whether to use repeated sequential simulation.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param exact_state Internal `exact_state` value used by this helper.
#' @return The internal `exo_select_stepwise_score_screening()` computation result.
#' @keywords internal
#' @noRd
exo_select_stepwise_score_screening <- function(score_category,
                                                score_valid,
                                                score_in_range,
                                                exo_values,
                                                backgrounds,
                                                marginal_p_min,
                                                score_cap,
                                                exact = FALSE,
                                                repeated = FALSE,
                                                nsim = 1000L,
                                                seed = NULL,
                                                exact_state = NULL) {
  if (is.null(exact_state)) {
    exact_state <- gRm_exact_state_from_flags(exact, repeated, nsim = nsim, seed = seed %||% 9L)
  }
  exact <- isTRUE(exact_state$exact)
  repeated <- isTRUE(exact_state$sequential)
  nsim <- as.integer(exact_state$nsim)
  seed <- as.integer(exact_state$seed)
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
        x_dim = score_cap + 1L,
        y_dim = backgrounds$raw_max[[exo_index]],
        condition_values = condition_values,
        condition_dims = backgrounds$raw_max[other_exo],
        valid = valid,
        exact = exact,
        nsim = nsim,
        seed = seed,
        repeated = repeated,
        seq_limit = as.integer(exact_state$seq_limit[[1L]]),
        seq_p0 = exact_state$seq_p0,
        seq_boundary = exact_state$seq_b
      )
      gamma_p_two_sided <- test$p_gamma
      p_min[[candidate_index]] <- min(test$p_chi, gamma_p_two_sided)
      rows[[candidate_index]] <- data.frame(
        exo_label = backgrounds$label_code[[exo_index]],
        exo_name = backgrounds$name[[exo_index]],
        hypothesis = screen_j_hypothesis_label("#", backgrounds$label_code[[exo_index]], backgrounds$label_code[other_exo]),
        chi_square = test$chi_square,
        df = test$df,
        chi_p = test$p_chi_asymp,
        gamma = test$gamma,
        gamma_p_one_sided = test$p_gamma_asymp / 2,
        gamma_p_two_sided = test$p_gamma_asymp,
        exact_chi_p = test$p_chi_exact,
        exact_gamma_p_one_sided = test$p_gamma_exact / 2,
        exact_gamma_p_two_sided = test$p_gamma_exact,
        exact_nsim = test$exact_nsim,
        selected = FALSE,
        chi_marker = "",
        gamma_marker = "",
        stringsAsFactors = FALSE
      )
    }

    rows_df <- do.call(rbind, rows)
    excluded <- integer()
    if (length(p_min) > 0L && max(p_min) > 0.05) {
      excluded <- current[[which.max(p_min)]]
    }
    iterations[[length(iterations) + 1L]] <- list(rows = rows_df, excluded = excluded)
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

#' Internal exo select repeated count cutoff helper
#'
#' Supports the exo select values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param project Encoded gRm project.
#' @param exact_state Internal `exact_state` value used by this helper.
#' @return The internal `exo_select_repeated_count_cutoff()` computation result.
#' @keywords internal
#' @noRd
exo_select_repeated_count_cutoff <- function(project, exact_state = NULL) {
  if (is.null(exact_state)) {
    exact_state <- gRm_exact_command_state("repeated")
  }
  as.integer(exact_state$seq_limit[[1L]])
}

#' Source Benjamini-Hochberg thresholds for score/exogenous screening
#'
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param screen Screening rows from `exo_select_score_screen()`.
#' @return List with FDR 0.05 and 0.01 source critical p-values.
#' @keywords internal
#' @noRd
exo_select_bh_thresholds <- function(screen) {
  exact_available <- "exact_chi_p" %in% names(screen) &&
    all(!is.na(screen$exact_chi_p)) &&
    all(!is.na(screen$exact_gamma_p_two_sided))
  p_values <- if (exact_available) {
    c(screen$exact_chi_p, screen$exact_gamma_p_two_sided)
  } else {
    c(screen$chi_p, screen$gamma_p_two_sided)
  }
  list(
    fdr_05 = source_bh_critical(p_values, 0.05),
    fdr_01 = source_bh_critical(p_values, 0.01)
  )
}

#' Source-shaped Monte Carlo exact score/exogenous test
#'
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param tab Observed score by exogenous table after source score collapsing.
#' @param observed_chi Observed chi-square statistic.
#' @param observed_gamma Observed RC gamma statistic.
#' @param nsim Number of random tables.
#' @return List of simulated exact p-values and effective simulation count.
#' @param random_draw Internal `random_draw` value used by this helper.
#' @param seed Random-stream seed.
#' @param sequential Internal `sequential` value used by this helper.
#' @param sequential_count_cutoff Internal `sequential_count_cutoff` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @keywords internal
#' @noRd
exo_select_exact_test <- function(tab,
                                  observed_chi,
                                  observed_gamma,
                                  nsim = 1000L,
                                  random_draw = NULL,
                                  seed = NULL,
                                  sequential = FALSE,
                                  sequential_count_cutoff = TRUE,
                                  seq_p0 = 0.05,
                                  seq_boundary = 1.058) {
  nsim <- as.integer(nsim)
  seq_limit <- if (isTRUE(sequential_count_cutoff)) {
    20L
  } else if (isFALSE(sequential_count_cutoff) || is.null(sequential_count_cutoff)) {
    0L
  } else {
    as.integer(sequential_count_cutoff[[1L]])
  }
  nchi <- 0L
  ngamma <- 0L
  nabsgamma <- 0L
  chi_status <- FALSE
  gamma_status <- FALSE
  seq_p0 <- as.numeric(seq_p0[[1L]])
  seq_boundary <- as.numeric(seq_boundary[[1L]])
  prepared <- exo_select_prepare_gentab1(tab)
  for (sim in seq_len(nsim)) {
    generated <- exo_select_gentab1_prepared(prepared, random_draw = random_draw)
    generated_chi <- screen_rc_chi_square(generated)
    generated_gamma <- screen_rc_gamma_counts(generated)$gamma
    # Source trace: SKexa2 counts random tables whose chi-square statistic is
    # at least as large as the observed statistic.
    if (generated_chi >= observed_chi) {
      nchi <- nchi + 1L
    }
    # Source trace: Results[*,7] is directional. The random gamma must be at
    # least as extreme in the observed sign direction, with gamma = 0 counting
    # every simulated table exactly as in the Pascal branch.
    if ((observed_gamma > 0 && generated_gamma >= observed_gamma) ||
        (observed_gamma < 0 && generated_gamma <= observed_gamma) ||
        observed_gamma == 0) {
      ngamma <- ngamma + 1L
    }
    # Source trace: Results[*,9], the value printed by EXA_SUMMARY1_2 for this
    # two-sided exact gamma column, counts absolute gamma exceedances.
    if (abs(generated_gamma) >= abs(observed_gamma)) {
      nabsgamma <- nabsgamma + 1L
    }
    if (isTRUE(sequential)) {
      chi_stop <- screen_j_seq_t(nchi, sim, seq_p0) >= seq_boundary
      gamma_stop <- screen_j_seq_t(nabsgamma, sim, seq_p0) >= seq_boundary
      if (seq_limit > 0L) {
        chi_stop <- chi_stop || nchi >= seq_limit
        gamma_stop <- gamma_stop || nabsgamma >= seq_limit
      }
      if (chi_stop) {
        chi_status <- TRUE
      }
      if (gamma_stop) {
        gamma_status <- TRUE
      }
      if (chi_status && gamma_status) {
        nsim <- sim
        break
      }
    }
  }
  list(
    chi_p = nchi / nsim,
    gamma_p_one_sided = ngamma / nsim,
    gamma_p_two_sided = nabsgamma / nsim,
    nsim = nsim
  )
}

#' Generate one conditional random table with fixed margins
#'
#' This is the R port of `source/PAS_scd/SKrandom.pas::GENTAB1`. It fills the
#' free cells sequentially using the source hypergeometric probability ordering
#' around the rounded expected cell and preserves the observed margins.
#'
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param tab Observed two-way table.
#' @return Generated two-way table with the same margins.
#' @param random_draw Internal `random_draw` value used by this helper.
#' @keywords internal
#' @noRd
exo_select_gentab1 <- function(tab, random_draw = NULL) {
  exo_select_gentab1_prepared(exo_select_prepare_gentab1(tab), random_draw = random_draw)
}

#' Internal exo select prepare gentab1 helper
#'
#' Supports the exo select values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param tab Internal `tab` value used by this helper.
#' @return The internal `exo_select_prepare_gentab1()` computation result.
#' @keywords internal
#' @noRd
exo_select_prepare_gentab1 <- function(tab) {
  tab <- as.matrix(tab)
  row_total <- as.integer(rowSums(tab))
  col_total <- as.integer(colSums(tab))
  grand_total <- sum(row_total)
  log_factorial <- numeric(grand_total + 1L)
  if (grand_total >= 2L) {
    source_limit <- min(grand_total, 1000L)
    log_factorial[seq.int(2L, source_limit) + 1L] <- cumsum(log(seq.int(2L, source_limit)))
    if (grand_total > 1000L) {
      n <- seq.int(1001L, grand_total)
      log_factorial[n + 1L] <- log((2 * pi)^2) + (n + 0.5) * log(n) - n + 1 / (24 * n) - 2.756
    }
  }
  list(
    cdim = nrow(tab),
    rdim = ncol(tab),
    row_total = row_total,
    col_total = col_total,
    grand_total = grand_total,
    log_factorial = log_factorial
  )
}

#' Internal exo select gentab1 prepared helper
#'
#' Supports the exo select values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param prepared Internal `prepared` value used by this helper.
#' @param random_draw Internal `random_draw` value used by this helper.
#' @return The internal `exo_select_gentab1_prepared()` computation result.
#' @keywords internal
#' @noRd
exo_select_gentab1_prepared <- function(prepared, random_draw = NULL) {
  cdim <- prepared$cdim
  rdim <- prepared$rdim
  row_total <- prepared$row_total
  col_total <- prepared$col_total
  grand_total <- prepared$grand_total
  log_factorial <- prepared$log_factorial
  newtab <- matrix(0L, nrow = cdim, ncol = rdim)
  generated_rows <- integer(rdim)
  generated_total <- 0L
  for (ic in seq_len(cdim - 1L)) {
    if (row_total[[ic]] == 0L) {
      next
    }
    generated_column <- 0L
    free_n <- grand_total - generated_total
    for (ir in seq_len(rdim - 1L)) {
      if ((col_total[[ir]] - generated_rows[[ir]]) == 0L || free_n == 0L) {
        t11 <- 0L
      } else {
        column1 <- row_total[[ic]] - generated_column
        row1 <- col_total[[ir]] - generated_rows[[ir]]
        column2 <- free_n - column1
        row2 <- free_n - row1
        tmax <- min(column1, row1)
        tmin <- max(column1 + row1 - free_n, 0L)
        expected <- round(column1 * (row1 / free_n))
        draw <- if (is.null(random_draw)) stats::runif(1L) else random_draw()
        cumulative <- 0
        cell_probability <- function(value) {
          t12 <- column1 - value
          t21 <- row1 - value
          t22 <- free_n - row1 - t12
          # Source trace: SKrandom.GENTAB1 uses the hypergeometric probability
          # for this free cell conditional on the already generated margins.
          exp(
            log_factorial[[column1 + 1L]] + log_factorial[[column2 + 1L]] +
              log_factorial[[row1 + 1L]] + log_factorial[[row2 + 1L]] -
              log_factorial[[value + 1L]] - log_factorial[[t12 + 1L]] -
              log_factorial[[t21 + 1L]] - log_factorial[[t22 + 1L]] -
              log_factorial[[free_n + 1L]]
          )
        }
        step_down <- expected - tmin
        step_up <- tmax - expected
        step_min <- min(step_down, step_up)
        t11 <- expected
        cumulative <- cumulative + cell_probability(t11)
        if (cumulative < draw) {
          found <- FALSE
          if (step_min >= 1L) {
            for (step in seq_len(step_min)) {
              t11 <- expected + step
              cumulative <- cumulative + cell_probability(t11)
              if (cumulative >= draw) {
                found <- TRUE
                break
              }
              t11 <- expected - step
              cumulative <- cumulative + cell_probability(t11)
              if (cumulative >= draw) {
                found <- TRUE
                break
              }
            }
          }
          if (!found) {
            if (step_down > step_up) {
              for (step in seq.int(step_min + 1L, step_down)) {
                t11 <- expected - step
                cumulative <- cumulative + cell_probability(t11)
                if (cumulative >= draw) break
              }
            } else if (step_up >= step_min + 1L) {
              for (step in seq.int(step_min + 1L, step_up)) {
                t11 <- expected + step
                cumulative <- cumulative + cell_probability(t11)
                if (cumulative >= draw) break
              }
            }
          }
        }
      }
      newtab[ic, ir] <- t11
      generated_total <- generated_total + t11
      generated_column <- generated_column + t11
      free_n <- free_n - col_total[[ir]] + generated_rows[[ir]]
      generated_rows[[ir]] <- generated_rows[[ir]] + t11
    }
    t11 <- row_total[[ic]] - generated_column
    newtab[ic, rdim] <- t11
    generated_total <- generated_total + t11
    generated_rows[[rdim]] <- generated_rows[[rdim]] + t11
  }
  newtab[cdim, ] <- col_total - generated_rows
  newtab
}

#' Wilson-style 99 percent source confidence limits
#'
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param n Number of simulations.
#' @param p Simulated p-value.
#' @return Numeric vector with lower and upper limits.
#' @keywords internal
#' @noRd
source_conflimit99 <- function(n, p) {
  z <- 2.5758
  zsquare <- z * z
  r <- p * n
  q <- 1 - p
  a <- 2 * r + zsquare
  b <- if (zsquare + 4 * r * q > 0) z * sqrt(zsquare + 4 * r * q) else 0
  c_value <- 2 * (n + zsquare)
  c(low = (a - b) / c_value, high = (a + b) / c_value)
}

#' Source gamma test for a score/exogenous table
#'
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param tab Contingency table.
#' @return List with gamma and one-sided p-value.
#' @keywords internal
#' @noRd
exo_select_gamma_test <- function(tab) {
  tab <- as.matrix(tab)
  cells <- gamma_cell_tables(tab)
  ppq <- cells$p + cells$q
  pmq <- cells$p - cells$q
  if (ppq <= 0) {
    return(list(gamma = 0, p_one_sided = 1))
  }
  n <- sum(tab)
  # Source trace: SKbias13/Inexpensive_bt_tests stores gamma in Results[1,5]
  # and the asymptotic one-sided gamma p-value in Results[1,6]. This variance
  # is the SKbias13 RC gamma statistic, not the fitted-gamma variance used by
  # item fit residual reports.
  s <- -pmq * (pmq / n)
  for (row in seq_len(nrow(tab))) {
    for (col in seq_len(ncol(tab))) {
      m <- cells$aij[row, col] - cells$dij[row, col]
      s <- s + tab[row, col] * m * m
    }
  }
  s <- 4 * s
  gamma <- pmq / ppq
  p_one <- if (s > 0) source_tail_norm(abs(gamma / (sqrt(s) / ppq)), TRUE) else 1
  list(gamma = gamma, p_one_sided = p_one)
}

#' Source significance marker
#'
#' Source trace: `source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening`.
#' @param p_value P-value used for the marker.
#' @param statistic Statistic sign source.
#' @param bh05 FDR 0.05 threshold.
#' @param bh01 FDR 0.01 threshold.
#' @param positive Direction for gamma markers.
#' @param chi Whether to use chi-square marker style.
#' @return Marker text.
#' @keywords internal
#' @noRd
exo_select_marker <- function(p_value, statistic, bh05, bh01, positive = TRUE, chi = FALSE) {
  if (p_value <= bh01) {
    if (chi) return("xx")
    return(if (positive) "++" else "--")
  }
  if (p_value <= bh05) {
    if (chi) return("x")
    return(if (positive) "+" else "-")
  }
  ""
}
