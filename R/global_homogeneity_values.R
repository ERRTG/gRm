#' Build DIGRAM score groups from source score cuts
#'
#' @param bundle Source-shaped bundle from [build_item_parameters_bundle()].
#' @param score_cuts Integer vector of upper score cuts.
#' @return A data frame with active score-group bounds.
#' @keywords internal
global_homogeneity_score_groups <- function(bundle, score_cuts) {
  # Source trace: DGRirtD.pas global homogeneity loop after scorecuts[0] :=
  # minscore - 1. The active group bounds are clipped to LeastScore/LargestScore
  # and skipped when outside the source-valid interior score range.
  if (length(score_cuts) < 2L) {
    stop("At least two score cuts are required for global homogeneity.", call. = FALSE)
  }

  previous_cut <- bundle$model$least_score - 1L
  rows <- list()
  for (cut_index in seq_along(score_cuts)) {
    from_score <- previous_cut + 1L
    if (bundle$model$least_score > from_score) {
      from_score <- bundle$model$least_score
    }
    to_score <- score_cuts[[cut_index]]
    if (bundle$model$largest_score < to_score) {
      to_score <- bundle$model$largest_score
    }
    if (
      from_score <= to_score &&
        from_score < bundle$model$max_total_score &&
        to_score > 0L
    ) {
      rows[[length(rows) + 1L]] <- data.frame(
        group = length(rows) + 1L,
        from_score = as.integer(from_score),
        to_score = as.integer(to_score),
        label = if (to_score > from_score) {
          sprintf("%d - %d", from_score, to_score)
        } else {
          sprintf("%d", from_score)
        },
        stringsAsFactors = FALSE
      )
    }
    previous_cut <- score_cuts[[cut_index]]
  }

  if (length(rows) == 0L) {
    stop("No active global homogeneity score groups remain after source clipping.", call. = FALSE)
  }
  do.call(rbind, rows)
}

global_homogeneity_score_group_lookup <- function(groups, max_score) {
  lookup <- rep(NA_integer_, max_score + 1L)
  if (nrow(groups) == 0L) {
    return(lookup)
  }
  for (group_index in seq_len(nrow(groups))) {
    from_score <- max(0L, groups$from_score[[group_index]])
    to_score <- min(max_score, groups$to_score[[group_index]])
    if (from_score > to_score) {
      next
    }
    positions <- seq.int(from_score, to_score) + 1L
    missing <- is.na(lookup[positions])
    lookup[positions[missing]] <- group_index
  }
  lookup
}

global_homogeneity_lookup_score <- function(lookup, score) {
  if (score < 0L || score >= length(lookup)) {
    return(NA_integer_)
  }
  lookup[[score + 1L]]
}

gRm_default_global_homogeneity_score_cuts <- function(project) {
  cuts <- tryCatch(
    as.integer(items_select_values(project)$score_groups$to_score),
    error = function(e) integer()
  )
  if (length(cuts) >= 2L) {
    return(cuts)
  }

  item_matrix <- as.matrix(project$raw_data[, project$items$position, drop = FALSE])
  complete <- rep(TRUE, nrow(item_matrix))
  for (item_index in seq_len(nrow(project$items))) {
    complete <- complete &
      item_matrix[, item_index] >= 1L &
      item_matrix[, item_index] <= project$items$raw_max[[item_index]]
  }
  scores <- rowSums(sweep(item_matrix, 2L, 1L, "-"))
  max_score <- if (any(complete)) max(scores[complete]) else sum(project$items$raw_max - 1L)
  unique_scores <- sort(unique(scores[complete]))
  first_cut <- if (length(unique_scores) >= 2L) {
    unique_scores[[max(1L, floor(length(unique_scores) / 2L))]]
  } else {
    floor(max_score / 2L)
  }
  as.integer(c(first_cut, max_score))
}

#' Subset a DIGRAM bundle to one score interval
#'
#' @param bundle Source-shaped bundle.
#' @param from_score Lower inclusive score bound.
#' @param to_score Upper inclusive score bound.
#' @return Bundle with rows outside the interval marked invalid.
#' @keywords internal
subset_bundle_to_score_group <- function(bundle, from_score, to_score) {
  # Source trace: DGRirtD.pas passes fra/til into Estimate_GLLRM. The data
  # records remain the same source rows, but only valid complete records whose
  # total score is inside the group contribute to fitting and margins.
  group_bundle <- bundle
  keep <- bundle$data$status == 1L &
    bundle$data$score >= from_score &
    bundle$data$score <= to_score
  group_bundle$data$status[!keep] <- 0L
  group_bundle$manifest$nvalid <- sum(keep)
  group_bundle$model$least_score <- as.integer(from_score)
  group_bundle$model$largest_score <- as.integer(to_score)
  group_bundle
}

#' Summarize DIGRAM item means and unmodeled residual cells
#'
#' @param bundle Score-group bundle.
#' @param fit Score-group Rasch fit.
#' @param group One-based score-group index.
#' @return Data frame of observed/expected item mean rows. The `residual`,
#'   hidden expected-variance, and marker cells are `NA`: the available Pascal
#'   source identifies the `skbias15.pas` item-mean residual path, but the
#'   historical runtime residual boundary has not been reproduced without
#'   empirical report-specific correction factors.
#' @keywords internal
global_homogeneity_item_mean_rows <- function(bundle, fit, group, expected_item_gamma) {
  # Source trace: skbias15.pas::Calculate_residuals_and_item_fits output(22).
  # CalculateMeans stores n, mean, and expected score variance at item_max+1,
  # item_max+2, and item_max+3; Meanres then standardizes the mean difference.
  # We can source-back n, observed mean, and expected mean from the available
  # Pascal path. The runtime's hidden variance materialization still differs in
  # some rounded example cells, so the printed residual and marker are left NA
  # instead of using report-derived correction factors.
  #
  items <- bundle$model$items
  counts <- fit$counts$item_counts
  score_item_n <- global_homogeneity_score_item_n(bundle)
  expected_tables <- global_homogeneity_expected_item_margin_tables(
    bundle,
    fit$counts,
    expected_item_gamma
  )
  rows <- vector("list", nrow(items))

  for (item_index in seq_len(nrow(items))) {
    item_max <- items$raw_max[[item_index]] - 1L
    item_scores <- seq.int(0L, item_max)
    score_names <- as.character(item_scores)
    observed_counts <- as.numeric(counts[item_index, score_names])
    n <- sum(observed_counts)

    observed_mean <- if (n > 0) {
      sum(item_scores * observed_counts) / n
    } else {
      0
    }
    expected_counts <- as.numeric(colSums(expected_tables[[item_index]])[score_names])
    expected_n <- sum(expected_counts)
    expected_mean <- if (expected_n > 0) {
      sum(item_scores * expected_counts) / expected_n
    } else {
      0
    }
    expected_variance <- NA_real_
    residual <- NA_real_
    marker <- NA_character_
    rows[[item_index]] <- data.frame(
      group = as.integer(group),
      item_label = items$label_code[[item_index]],
      item_name = items$name[[item_index]],
      n = as.integer(n),
      observed_mean = observed_mean,
      expected_mean = expected_mean,
      expected_variance = expected_variance,
      visible_source_residual = residual,
      residual = residual,
      marker = marker,
      residual_runtime_source_backed = FALSE,
      marker_runtime_source_backed = FALSE,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

#' Calculate DIGRAM's score-weighted expected item mean variances
#'
#' @param bundle Score-group bundle.
#' @param counts Count list from [rasch_counts()].
#' @param item_gamma Item gamma matrix used for expected margins.
#' @param expected_tables Optional precomputed expected item margin tables.
#' @param score_item_n Optional observed score-by-item row totals.
#' @return Numeric vector with one variance per item.
#' @keywords internal
global_homogeneity_expected_score_variances <- function(bundle,
                                                        counts,
                                                        item_gamma,
                                                        expected_tables = NULL,
                                                        score_item_n = NULL) {
  items <- bundle$model$items
  if (is.null(score_item_n)) {
    score_item_n <- global_homogeneity_score_item_n(bundle)
  }

  result <- numeric(nrow(items))
  max_total_score <- bundle$model$max_total_score
  if (max_total_score < 2L) {
    return(result)
  }

  if (is.null(expected_tables)) {
    expected_tables <- global_homogeneity_expected_item_margin_tables(
      bundle,
      counts,
      item_gamma
    )
  }
  # Source trace: skbias15.pas::CalculateExpectedValues builds
  # ExpectedItemMargTables from full and item-excluded gamma arrays:
  # e = g_item * g_without_item * n_s / g_full_s. CalculateMeans then applies
  # SummarizeTal to those accumulated expected tables, and SummarizeVar combines
  # the per-score variances with observed score-count weights.

  for (item_index in seq_len(nrow(items))) {
    item_max <- items$raw_max[[item_index]] - 1L
    expected_table <- expected_tables[[item_index]]
    total_summary <- global_homogeneity_summarize_tal(
      colSums(expected_table),
      item_max
    )
    total_expected_n <- total_summary$n
    if (total_expected_n <= 0) {
      next
    }

    for (score in seq.int(1L, max_total_score - 1L)) {
      observed_score_n <- score_item_n[score + 1L, item_index]
      if (observed_score_n <= 0) {
        next
      }
      score_summary <- global_homogeneity_summarize_tal(
        expected_table[score + 1L, ],
        item_max
      )
      result[[item_index]] <- result[[item_index]] +
        (observed_score_n / total_expected_n) * score_summary$variance
    }
  }

  result
}

#' Calculate global homogeneity expected item summaries in the native C++ kernel
#'
#' @param bundle Score-group bundle.
#' @param counts Count list from [rasch_counts()].
#' @param item_gamma Item gamma matrix used for expected margins.
#' @param score_item_n Score-by-item observed row totals.
#' @return List with expected variance, n, and mean vectors, or `NULL`.
#' @keywords internal
global_homogeneity_expected_summary_cpp <- function(bundle,
                                                    counts,
                                                    item_gamma,
                                                    score_item_n = NULL) {
  if (is.null(score_item_n)) {
    score_item_n <- global_homogeneity_score_item_n(bundle)
  }
  tryCatch(
    .Call(
      "gRm_global_homogeneity_expected_summary",
      item_gamma,
      counts$score_counts,
      score_item_n,
      bundle$model$items$raw_max,
      bundle$model$least_score,
      bundle$model$largest_score,
      PACKAGE = "gRm"
    ),
    error = function(e) NULL
  )
}

#' Count observed item-margin row totals by score
#'
#' @param bundle Score-group bundle.
#' @return Numeric score-by-item matrix of observed item row totals.
#' @keywords internal
global_homogeneity_score_item_n <- function(bundle) {
  # Source trace: skbias15.pas::SummarizeVar uses
  # ItemMargTables[score,item,a] as the per-score numerator. This equals the
  # score count for complete data, but remains item-specific when missing item
  # responses are present.
  items <- bundle$model$items
  max_total_score <- bundle$model$max_total_score
  result <- matrix(
    0,
    nrow = max_total_score + 1L,
    ncol = nrow(items),
    dimnames = list(
      as.character(seq.int(0L, max_total_score)),
      items$name
    )
  )
  valid <- bundle$data$status == 1L
  for (item_index in seq_len(nrow(items))) {
    item_values <- bundle$data[[items$name[[item_index]]]]
    item_valid <- valid & !is.na(item_values)
    result[, item_index] <- tabulate(
      bundle$data$score[item_valid] + 1L,
      nbins = max_total_score + 1L
    )
  }
  result
}

#' Build Pascal-shaped expected item margin tables for global homogeneity
#'
#' @param bundle Score-group bundle.
#' @param counts Count list from [rasch_counts()].
#' @param item_gamma Item gamma matrix used for expected margins.
#' @return List of score-by-item-score expected margin matrices.
#' @keywords internal
global_homogeneity_expected_item_margin_tables <- function(bundle, counts, item_gamma) {
  items <- bundle$model$items
  max_total_score <- bundle$model$max_total_score
  use_items <- rep(TRUE, nrow(items))
  full_gamma <- global_homogeneity_score_gamma(bundle, item_gamma, use_items)
  tables <- vector("list", nrow(items))
  source_item_max <- max(items$raw_max) - 1L
  scfra <- if (bundle$model$least_score > 1L) bundle$model$least_score else 1L
  sctil <- if (bundle$model$largest_score < max_total_score - 1L) {
    bundle$model$largest_score
  } else {
    max_total_score - 1L
  }

  for (item_index in seq_len(nrow(items))) {
    item_max <- items$raw_max[[item_index]] - 1L
    table <- matrix(
      0,
      nrow = max_total_score + 1L,
      ncol = item_max + 1L,
      dimnames = list(
        as.character(seq.int(0L, max_total_score)),
        as.character(seq.int(0L, item_max))
      )
    )
    use_items[[item_index]] <- FALSE
    without_item <- global_homogeneity_score_gamma(bundle, item_gamma, use_items)
    use_items[[item_index]] <- TRUE

    for (item_score in seq.int(0L, item_max)) {
      first_score <- item_score
      if (first_score < scfra) {
        first_score <- scfra
      }
      last_score <- item_score + (nrow(items) - 1L) * source_item_max
      if (last_score > sctil) {
        last_score <- sctil
      }
      if (first_score > last_score) {
        next
      }
      for (score in seq.int(first_score, last_score)) {
        score_n <- counts$score_counts[[score + 1L]]
        g_full <- full_gamma[[score + 1L]]
        rest_score <- score - item_score
        if (score_n <= 0 || g_full <= 0 || rest_score < 0L) {
          next
        }
        # skbias15.pas computes e := g1 * g2 * (n / g). The order is
        # observable for sparse score rows that sit on the tal[a] > 1 branch.
        ratio <- score_n / g_full
        numerator <- item_gamma[item_index, as.character(item_score)] *
          without_item[[rest_score + 1L]]
        table[score + 1L, item_score + 1L] <-
          table[score + 1L, item_score + 1L] + numerator * ratio
      }
    }

    tables[[item_index]] <- table
  }

  tables
}

#' Build the component gamma array used by skbias15 expected item margins
#'
#' @param bundle Source-shaped bundle.
#' @param item_gamma Item gamma matrix.
#' @param excluded_item Optional one-based item index removed from the
#'   convolution.
#' @return Numeric vector indexed by total score plus one.
#' @keywords internal
global_homogeneity_component_score_gamma <- function(bundle, item_gamma, excluded_item = NA_integer_) {
  # Source trace: skbias15.pas::CalculateExpectedValues removes one item
  # component with CalculateBiasedGammaValues2 before RunThroughItemComponent.
  # For the examples' no-LD/no-DIF models, each component is one item, so this
  # is the direct score-polynomial convolution with that item optionally
  # omitted.
  items <- bundle$model$items
  max_total_score <- bundle$model$max_total_score
  gamma_values <- numeric(max_total_score + 1L)
  gamma_values[[1L]] <- 1
  current_max <- 0L
  excluded_item <- as.integer(excluded_item)[[1L]]

  for (item_index in seq_len(nrow(items))) {
    if (!is.na(excluded_item) && item_index == excluded_item) {
      next
    }
    next_values <- numeric(max_total_score + 1L)
    next_max <- current_max + items$raw_max[[item_index]] - 1L
    for (score in seq.int(0L, current_max)) {
      current_weight <- gamma_values[[score + 1L]]
      if (current_weight == 0) {
        next
      }
      for (item_score in seq.int(0L, items$raw_max[[item_index]] - 1L)) {
        target_score <- score + item_score
        next_values[[target_score + 1L]] <-
          next_values[[target_score + 1L]] +
          current_weight * item_gamma[item_index, as.character(item_score)]
      }
    }
    gamma_values <- next_values
    current_max <- next_max
  }

  gamma_values
}

#' Build the full score gamma array used by ExpectedItemMargTables
#'
#' @param bundle Score-group bundle.
#' @param item_gamma Item gamma matrix used for expected margins.
#' @return Numeric vector indexed by total score plus one.
#' @keywords internal
global_homogeneity_full_score_gamma <- function(bundle, item_gamma) {
  global_homogeneity_score_gamma(
    bundle,
    item_gamma,
    use_items = rep(TRUE, nrow(bundle$model$items))
  )
}

#' Build a score gamma array using the skbias12 Inexpensive_Gamma_Calculation order
#'
#' @param bundle Score-group bundle.
#' @param item_gamma Item gamma matrix used for expected margins.
#' @param use_items Logical vector selecting items included in the score gamma.
#' @return Numeric vector indexed by total score plus one.
#' @keywords internal
global_homogeneity_score_gamma <- function(bundle, item_gamma, use_items) {
  build_source_score_gamma(bundle, item_gamma, use_items)
}

#' Run skbias15.pas SummarizeTal on item-score cells
#'
#' @param cells Numeric expected item-score cells.
#' @param item_max Maximum item score for the item.
#' @return List containing Pascal `tal[a]`, `tal[b]`, and `tal[c]`.
#' @keywords internal
global_homogeneity_summarize_tal <- function(cells, item_max) {
  n <- 0
  mean_value <- 0
  second_moment <- 0

  for (item_score in seq.int(0L, item_max)) {
    cell <- cells[[item_score + 1L]]
    n <- n + cell
    mean_value <- mean_value + item_score * cell
    second_moment <- second_moment + item_score * item_score * cell
  }

  if (n > 0) {
    mean_value <- mean_value / n
    if (n > 1) {
      second_moment <- second_moment / n
    }
  }

  variance <- if (n > 1) second_moment - mean_value * mean_value else 0
  list(n = n, mean = mean_value, variance = variance)
}

#' Derive DIGRAM global homogeneity numeric values
#'
#' Computes non-GUI DIGRAM global homogeneity values from raw DIGRAM
#' input using native R code.
#'
#' The item mean rows are computed from the shared `skbias15.pas`-shaped helper
#' used by both global homogeneity and global invariance. Row counts, observed
#' means, expected means, and summary CLR/df/p values are protected by tests.
#' The printed item-mean residual and marker cells are deliberately `NA`
#' because the available source does not fully explain the historical runtime's
#' hidden residual variance materialization.
#'
#' @param project A parsed DIGRAM project from [read_digram_project()].
#' @param score_cuts Integer upper score cuts. For the supplied validation runtime
#'   example this is `c(30, 87)`, corresponding to `CUT 30 87`.
#' @param max_step Maximum number of Rasch IPF iterations.
#' @param max_delta Convergence threshold for Rasch IPF.
#' @return A `gRm_global_homogeneity_values` object with the full-model fit,
#'   score-group fits, item mean rows, and CLR summary.
#' @examples
#' \dontrun{
#' project <- read_digram_project("path/to/DIGRAM")
#' values <- global_homogeneity_values(project, score_cuts = c(30, 87))
#' values$summary
#' }
#' @keywords internal
global_homogeneity_values <- function(project,
                                      score_cuts,
                                      max_step = 5000L,
                                      max_delta = 0.0001,
                                      bundle = NULL,
                                      base_fit = NULL,
                                      fit = NULL) {
  if ((inherits(project, "gRm_fit") || inherits(project, "gRm_gllrm_fit")) && inherits(project$values, "gRm_active_gllrm_values")) {
    return(active_gllrm_global_homogeneity_values(
      project,
      score_cuts = score_cuts,
      max_step = max_step,
      max_delta = max_delta
    ))
  }
  bundle <- bundle %||% build_item_parameters_bundle(project)
  base_fit <- base_fit %||% fit
  full_fit <- base_fit %||% fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)
  item_parameters <- item_parameters_values(full_fit, bundle)
  expected_item_gamma <- full_fit$item_gamma
  groups <- global_homogeneity_score_groups(bundle, as.integer(score_cuts))
  groups$n <- integer(nrow(groups))
  groups$log_likelihood <- numeric(nrow(groups))
  groups$converged <- logical(nrow(groups))
  groups$delta <- numeric(nrow(groups))

  group_values <- vector("list", nrow(groups))
  item_rows <- vector("list", nrow(groups))
  subgroup_loglike_sum <- 0

  for (group_index in seq_len(nrow(groups))) {
    group_bundle <- subset_bundle_to_score_group(
      bundle,
      groups$from_score[[group_index]],
      groups$to_score[[group_index]]
    )
    group_fit <- fit_rasch_base(group_bundle, max_step = max_step, max_delta = max_delta)
    group_loglike <- base_rasch_loglike(group_bundle, group_fit$item_gamma)
    subgroup_loglike_sum <- subgroup_loglike_sum + group_loglike
    group_values[[group_index]] <- list(
      bundle = group_bundle,
      fit = group_fit,
      log_likelihood = group_loglike
    )
    item_rows[[group_index]] <- global_homogeneity_item_mean_rows(
      group_bundle,
      group_fit,
      groups$group[[group_index]],
      expected_item_gamma
    )
    groups$n[[group_index]] <- group_fit$counts$n_valid
    groups$log_likelihood[[group_index]] <- group_loglike
    groups$converged[[group_index]] <- group_fit$converged
    groups$delta[[group_index]] <- group_fit$delta
  }

  full_loglike <- base_rasch_loglike(bundle, full_fit$item_gamma)
  n_parameters <- calculate_source_n_parameters(full_fit$counts$item_counts)
  clr <- 2 * abs(full_loglike - subgroup_loglike_sum)
  df <- (nrow(groups) - 1L) * n_parameters
  p_value <- source_pfchi(df, clr)

  result <- list(
    item_parameters = item_parameters,
    bundle = bundle,
    fit = full_fit,
    score_groups = groups,
    group_values = group_values,
    items = do.call(rbind, item_rows),
    summary = list(
      n_groups = nrow(groups),
      n_parameters = n_parameters,
      full_log_likelihood = full_loglike,
      subgroup_log_likelihood_sum = subgroup_loglike_sum,
      clr = clr,
      df = df,
      p_value = p_value
    )
  )
  class(result) <- c("gRm_global_homogeneity_values", class(result))
  result
}

active_gllrm_global_homogeneity_values <- function(fit,
                                                   score_cuts,
                                                   max_step = 5000L,
                                                   max_delta = 0.0001) {
  bundle <- fit$bundle %||% fit$fit$context$bundle
  groups <- global_homogeneity_score_groups(bundle, as.integer(score_cuts))
  groups$n <- integer(nrow(groups))
  groups$log_likelihood <- numeric(nrow(groups))
  groups$converged <- logical(nrow(groups))
  groups$delta <- numeric(nrow(groups))

  group_values <- vector("list", nrow(groups))
  item_rows <- vector("list", nrow(groups))
  subgroup_loglike_sum <- 0
  full_state <- fit$fit
  full_state$context <- NULL

  for (group_index in seq_len(nrow(groups))) {
    group_bundle <- subset_bundle_to_score_group(
      bundle,
      groups$from_score[[group_index]],
      groups$to_score[[group_index]]
    )
    group_active <- fit_gllrm_active(
      fit$spec,
      max_step = max_step,
      max_delta = max_delta,
      bundle = group_bundle
    )
    group_values_object <- gllrm_active_values(group_active, fit$spec)
    group_loglike <- group_active$state$log_likelihood
    subgroup_loglike_sum <- subgroup_loglike_sum + group_loglike
    group_values[[group_index]] <- list(
      bundle = group_bundle,
      fit = list(
        context = group_active$context,
        state = group_active$state,
        values = group_values_object
      ),
      log_likelihood = group_loglike
    )
    item_rows[[group_index]] <- active_global_homogeneity_item_mean_rows(
      group_active$context,
      full_state,
      groups$group[[group_index]],
      probability_cache = new_active_gllrm_probability_cache(group_active$context, full_state)
    )
    groups$n[[group_index]] <- group_active$context$counts$n_valid
    groups$log_likelihood[[group_index]] <- group_loglike
    groups$converged[[group_index]] <- group_active$state$converged
    groups$delta[[group_index]] <- group_active$state$report_delta
  }

  full_loglike <- fit$values$log_likelihood
  n_parameters <- fit$values$n_parameters
  clr <- 2 * abs(full_loglike - subgroup_loglike_sum)
  df <- (nrow(groups) - 1L) * n_parameters
  p_value <- source_pfchi(df, clr)
  full_probability_cache <- new_active_gllrm_probability_cache(fit$fit$context, fit$fit)

  result <- list(
    item_parameters = fit$values,
    bundle = bundle,
    fit = fit$fit,
    score_groups = groups,
    group_values = group_values,
    items = do.call(rbind, item_rows),
    uniform_ld = active_global_homogeneity_uniform_ld(
      fit,
      groups,
      probability_cache = full_probability_cache
    ),
    uniform_dif = active_global_homogeneity_uniform_dif(
      fit,
      groups,
      probability_cache = full_probability_cache
    ),
    summary = list(
      n_groups = nrow(groups),
      n_parameters = n_parameters,
      full_log_likelihood = full_loglike,
      subgroup_log_likelihood_sum = subgroup_loglike_sum,
      clr = clr,
      df = df,
      p_value = p_value
    )
  )
  class(result) <- c("gRm_global_homogeneity_values", class(result))
  result
}

active_global_homogeneity_uniform_ld <- function(fit, groups, probability_cache = NULL) {
  context <- fit$fit$context
  state <- fit$fit
  if (length(context$ld_specs) == 0L) {
    return(data.frame())
  }

  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  probability_cache <- probability_cache %||%
    new_active_gllrm_probability_cache(context, state, components = components)
  tables_by_ld <- active_uniform_ld_scoregroup_tables_all(context, groups, probability_cache)
  rows <- vector("list", length(context$ld_specs))
  for (ld_index in seq_along(context$ld_specs)) {
    spec <- context$ld_specs[[ld_index]]
    rows[[ld_index]] <- active_uniform_ld_summary_row(context, groups, spec, tables_by_ld[[ld_index]])
  }
  do.call(rbind, rows)
}

active_uniform_ld_scoregroup_tables <- function(context, groups, spec, ld_index, probability_cache) {
  active_uniform_ld_scoregroup_tables_all(context, groups, probability_cache)[[ld_index]]
}

active_uniform_ld_scoregroup_tables_all <- function(context, groups, probability_cache) {
  tables_by_ld <- lapply(context$ld_specs, function(spec) {
    dimensions <- c(context$item_raw_max[[spec$item1]], context$item_raw_max[[spec$item2]], nrow(groups))
    list(
      observed = array(0, dim = dimensions),
      expected = array(0, dim = dimensions)
    )
  })
  if (length(tables_by_ld) == 0L) {
    return(tables_by_ld)
  }

  score_group_lookup <- global_homogeneity_score_group_lookup(groups, context$max_total_score)
  for (row in context$valid_rows) {
    group_index <- global_homogeneity_lookup_score(score_group_lookup, context$score[[row]])
    if (is.na(group_index)) {
      next
    }
    for (ld_index in seq_along(context$ld_specs)) {
      spec <- context$ld_specs[[ld_index]]
      score1 <- context$item_matrix[row, spec$item1] + 1L
      score2 <- context$item_matrix[row, spec$item2] + 1L
      tables_by_ld[[ld_index]]$observed[score1, score2, group_index] <-
        tables_by_ld[[ld_index]]$observed[score1, score2, group_index] + 1
    }
  }

  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    score <- group$score[[1L]]
    homogeneity_group <- global_homogeneity_lookup_score(score_group_lookup, score)
    if (is.na(homogeneity_group)) {
      next
    }
    background_values <- gllrm_group_background_values(context, group)
    probabilities_by_ld <- active_gllrm_cached_ld_probabilities(
      probability_cache,
      total_score = score,
      background_values = background_values
    )
    for (ld_index in seq_along(context$ld_specs)) {
      tables_by_ld[[ld_index]]$expected[, , homogeneity_group] <-
        tables_by_ld[[ld_index]]$expected[, , homogeneity_group] +
        group$count[[1L]] * probabilities_by_ld[[ld_index]]
    }
  }

  tables_by_ld
}

active_uniform_ld_summary_row <- function(context, groups, spec, tables) {
  n_groups <- nrow(groups)
  observed_gamma <- numeric(n_groups)
  expected_gamma <- numeric(n_groups)
  chi_square <- 0
  df <- 0L
  x_dim <- dim(tables$observed)[[1L]]
  y_dim <- dim(tables$observed)[[2L]]

  for (group_index in seq_len(n_groups)) {
    observed <- tables$observed[, , group_index, drop = FALSE][, , 1L]
    expected <- tables$expected[, , group_index, drop = FALSE][, , 1L]
    positive <- expected > 0
    if (any(positive)) {
      chi_square <- chi_square + sum((observed[positive] - expected[positive])^2 / expected[positive])
    }
    standardized <- active_uniform_dif_standardize_expected(observed, expected)
    observed_gamma[[group_index]] <- source_gamma_from_table(observed)
    expected_gamma[[group_index]] <- source_gamma_from_table(standardized$table)
    df <- df + standardized$df
  }
  df <- df - (x_dim - 1L) * (y_dim - 1L)
  if (df < 1L) {
    df <- 1L
  }

  data.frame(
    item1 = spec$item1,
    item2 = spec$item2,
    item1_label = context$items$label_code[[spec$item1]],
    item2_label = context$items$label_code[[spec$item2]],
    observed_gamma = I(list(observed_gamma)),
    expected_gamma = I(list(expected_gamma)),
    chi_square = chi_square,
    df = as.integer(df),
    p_value = source_pfchi(df, chi_square),
    stringsAsFactors = FALSE
  )
}

active_gllrm_group_ld_probabilities <- function(context,
                                                state,
                                                total_score,
                                                background_values,
                                                components = NULL) {
  components <- components %||% context$ld_components_items %||% gllrm_ld_components(context)$items
  out <- lapply(context$ld_specs, function(spec) {
    matrix(0, nrow = context$item_raw_max[[spec$item1]], ncol = context$item_raw_max[[spec$item2]])
  })
  if (total_score < 0L || total_score > context$max_total_score) {
    return(out)
  }

  component_gamma <- lapply(components, function(component_items) {
    one <- gllrm_component_gamma(context, state, component_items, background_values)
    one$gamma * one$scale
  })
  convolutions <- gllrm_component_convolutions(component_gamma, context$max_total_score)
  denominator <- convolutions$full[[total_score + 1L]]
  if (denominator <= 0) {
    return(out)
  }

  for (component_index in seq_along(components)) {
    component_items <- components[[component_index]]
    key <- gllrm_component_key(component_items)
    ld_local <- context$component_ld_local_matrices[[key]] %||%
      context$component_ld_local_indices[[key]]
    if (is.null(ld_local) || nrow(ld_local) == 0L) {
      next
    }
    rest_gamma <- convolutions$rest[[component_index]]
    configs <- context$component_config_matrices[[key]]
    config_scores <- context$component_config_scores[[key]]
    weights <- gllrm_component_config_weights_fast(
      context,
      state,
      component_items,
      configs,
      background_values,
      key = key
    )
    valid <- config_scores <= total_score & weights > 0
    if (!any(valid)) {
      next
    }
    valid_index <- which(valid)
    rest_values <- rest_gamma[total_score - config_scores[valid_index] + 1L]
    probabilities <- weights[valid_index] * rest_values / denominator
    positive <- probabilities > 0
    if (!any(positive)) {
      next
    }
    valid_index <- valid_index[positive]
    probabilities <- probabilities[positive]

    for (ld_row in seq_len(nrow(ld_local))) {
      ld_index <- ld_local[ld_row, 1L]
      score1 <- configs[valid_index, ld_local[ld_row, 2L]] + 1L
      score2 <- configs[valid_index, ld_local[ld_row, 3L]] + 1L
      n_score1 <- nrow(out[[ld_index]])
      pair_index <- score1 + (score2 - 1L) * n_score1
      by_pair <- rowsum(probabilities, pair_index, reorder = FALSE)
      out[[ld_index]][as.integer(rownames(by_pair))] <-
        out[[ld_index]][as.integer(rownames(by_pair))] + as.numeric(by_pair[, 1L])
    }
  }

  out
}

active_global_homogeneity_uniform_dif <- function(fit, groups, probability_cache = NULL) {
  context <- fit$fit$context
  state <- fit$fit
  if (length(context$dif_specs) == 0L) {
    return(data.frame())
  }

  rows <- vector("list", length(context$dif_specs))
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  probability_cache <- probability_cache %||%
    new_active_gllrm_probability_cache(context, state, components = components)
  for (dif_index in seq_along(context$dif_specs)) {
    spec <- context$dif_specs[[dif_index]]
    tables <- active_uniform_dif_scoregroup_tables(context, groups, spec, probability_cache)
    rows[[dif_index]] <- active_uniform_dif_summary_row(context, groups, spec, tables)
  }
  do.call(rbind, rows)
}

active_uniform_dif_scoregroup_tables <- function(context, groups, spec, probability_cache) {
  item_max <- context$item_raw_max[[spec$item]] - 1L
  background_max <- context$background_raw_max[[spec$background]]
  dimensions <- c(item_max + 1L, background_max, nrow(groups))
  observed <- array(0, dim = dimensions)
  expected <- array(0, dim = dimensions)

  score_group_lookup <- global_homogeneity_score_group_lookup(groups, context$max_total_score)

  for (row in context$valid_rows) {
    group_index <- global_homogeneity_lookup_score(score_group_lookup, context$score[[row]])
    if (is.na(group_index)) {
      next
    }
    item_score <- context$item_matrix[row, spec$item] + 1L
    background_value <- context$background_matrix[row, spec$background]
    observed[item_score, background_value, group_index] <-
      observed[item_score, background_value, group_index] + 1
  }

  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    score <- group$score[[1L]]
    homogeneity_group <- global_homogeneity_lookup_score(score_group_lookup, score)
    if (is.na(homogeneity_group)) {
      next
    }
    background_values <- gllrm_group_background_values(context, group)
    background_value <- background_values[[spec$background]]
    probabilities <- active_gllrm_cached_item_probabilities(
      probability_cache,
      total_score = score,
      background_values = background_values
    )[[spec$item]]
    expected[seq_along(probabilities), background_value, homogeneity_group] <-
      expected[seq_along(probabilities), background_value, homogeneity_group] +
        group$count[[1L]] * probabilities
  }

  list(observed = observed, expected = expected)
}

active_uniform_dif_summary_row <- function(context, groups, spec, tables) {
  n_groups <- nrow(groups)
  observed_gamma <- numeric(n_groups)
  expected_gamma <- numeric(n_groups)
  chi_square <- 0
  df <- 0L
  x_dim <- dim(tables$observed)[[1L]]
  y_dim <- dim(tables$observed)[[2L]]

  for (group_index in seq_len(n_groups)) {
    observed <- tables$observed[, , group_index, drop = FALSE][, , 1L]
    expected <- tables$expected[, , group_index, drop = FALSE][, , 1L]
    positive <- expected > 0
    if (any(positive)) {
      chi_square <- chi_square + sum((observed[positive] - expected[positive])^2 / expected[positive])
    }
    standardized <- active_uniform_dif_standardize_expected(observed, expected)
    observed_gamma[[group_index]] <- source_gamma_from_table(observed)
    expected_gamma[[group_index]] <- source_gamma_from_table(standardized$table)
    df <- df + standardized$df
  }
  df <- df - (x_dim - 1L) * (y_dim - 1L)
  if (df < 1L) {
    df <- 1L
  }

  data.frame(
    item = spec$item,
    background = spec$background,
    item_label = context$items$label_code[[spec$item]],
    background_label = context$backgrounds$label_code[[spec$background]],
    observed_gamma = I(list(observed_gamma)),
    expected_gamma = I(list(expected_gamma)),
    chi_square = chi_square,
    df = as.integer(df),
    p_value = source_pfchi(df, chi_square),
    stringsAsFactors = FALSE
  )
}

active_uniform_dif_standardize_expected <- function(observed, expected) {
  out <- expected
  row_margins <- rowSums(observed)
  col_margins <- colSums(observed)
  for (step in seq_len(30L)) {
    row_sums <- rowSums(out)
    for (row in seq_along(row_sums)) {
      if (row_sums[[row]] > 0) {
        out[row, ] <- out[row, ] * row_margins[[row]] / row_sums[[row]]
      }
    }
    col_sums <- colSums(out)
    for (col in seq_along(col_sums)) {
      if (col_sums[[col]] > 0) {
        out[, col] <- out[, col] * col_margins[[col]] / col_sums[[col]]
      }
    }
  }
  list(
    table = out,
    df = max(0L, (sum(row_margins > 0) - 1L) * (sum(col_margins > 0) - 1L))
  )
}

active_global_homogeneity_item_mean_rows <- function(context,
                                                     full_state,
                                                     group,
                                                     probability_cache = NULL) {
  items <- context$items
  counts <- context$counts$item_counts
  expected_tables <- active_global_homogeneity_expected_item_margin_tables(
    context,
    full_state,
    probability_cache = probability_cache
  )
  expected_variance <- active_global_homogeneity_item_variances(
    context,
    expected_tables = expected_tables
  )
  rows <- vector("list", nrow(items))

  for (item_index in seq_len(nrow(items))) {
    item_max <- items$raw_max[[item_index]] - 1L
    item_scores <- seq.int(0L, item_max)
    score_cols <- item_scores + 1L
    observed_counts <- as.numeric(counts[item_index, as.character(item_scores)])
    n <- sum(observed_counts)
    observed_mean <- if (n > 0) {
      sum(item_scores * observed_counts) / n
    } else {
      0
    }
    expected_counts <- colSums(expected_tables[[item_index]])[score_cols]
    expected_n <- sum(expected_counts)
    expected_mean <- if (expected_n > 0) {
      sum(item_scores * expected_counts) / expected_n
    } else {
      0
    }
    item_expected_variance <- expected_variance[[item_index]]
    # The available Pascal source backs the active-model counts, expected
    # means, score-group likelihoods, and CLR, but it does not reproduce the
    # runtime denominator used for this compact residual column. Historical
    # harness traces in this repository end at the same boundary, so the
    # residual and marker are deliberately exposed as unmodeled values.
    residual <- NA_real_
    marker <- NA_character_
    rows[[item_index]] <- data.frame(
      group = as.integer(group),
      item_label = items$label_code[[item_index]],
      item_name = items$name[[item_index]],
      n = as.integer(n),
      observed_mean = observed_mean,
      expected_mean = expected_mean,
      expected_variance = item_expected_variance,
      visible_source_residual = residual,
      residual = residual,
      marker = marker,
      residual_runtime_source_backed = FALSE,
      marker_runtime_source_backed = FALSE,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

active_global_homogeneity_expected_item_margin_tables <- function(context, full_state, probability_cache = NULL) {
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  probability_cache <- probability_cache %||%
    new_active_gllrm_probability_cache(context, full_state, components = components)
  tables <- lapply(seq_len(context$n_items), function(item_index) {
    item_max <- context$item_raw_max[[item_index]] - 1L
    matrix(
      0,
      nrow = context$max_total_score + 1L,
      ncol = item_max + 1L,
      dimnames = list(
        as.character(seq.int(0L, context$max_total_score)),
        as.character(seq.int(0L, item_max))
      )
    )
  })

  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    total_score <- group$score[[1L]]
    background_values <- gllrm_group_background_values(context, group)
    probabilities <- active_gllrm_cached_item_probabilities(
      probability_cache,
      total_score = total_score,
      background_values = background_values
    )
    for (item_index in seq_len(context$n_items)) {
      p <- probabilities[[item_index]]
      tables[[item_index]][total_score + 1L, seq_along(p)] <-
        tables[[item_index]][total_score + 1L, seq_along(p)] +
        group$count[[1L]] * p
    }
  }

  tables
}

active_global_homogeneity_item_variances <- function(context, expected_tables) {
  out <- numeric(context$n_items)
  score_item_n <- active_global_homogeneity_score_item_n(context)
  scfra <- if (min(context$score[context$valid_rows]) == 0L) {
    1L
  } else {
    min(context$score[context$valid_rows])
  }
  sctil <- if (max(context$score[context$valid_rows]) == context$max_total_score) {
    context$max_total_score - 1L
  } else {
    max(context$score[context$valid_rows])
  }

  for (item_index in seq_len(context$n_items)) {
    item_max <- context$item_raw_max[[item_index]] - 1L
    expected_table <- expected_tables[[item_index]]
    total_summary <- global_homogeneity_summarize_tal(
      colSums(expected_table),
      item_max
    )
    if (total_summary$n <= 0 || scfra > sctil) {
      next
    }
    for (score in seq.int(scfra, sctil)) {
      observed_score_n <- score_item_n[score + 1L, item_index]
      if (observed_score_n <= 0) {
        next
      }
      score_summary <- global_homogeneity_summarize_tal(
        expected_table[score + 1L, ],
        item_max
      )
      out[[item_index]] <- out[[item_index]] +
        (observed_score_n / total_summary$n) * score_summary$variance
    }
  }

  out
}

active_global_homogeneity_score_item_n <- function(context) {
  out <- matrix(
    0,
    nrow = context$max_total_score + 1L,
    ncol = context$n_items,
    dimnames = list(
      as.character(seq.int(0L, context$max_total_score)),
      context$items$name
    )
  )
  for (item_index in seq_len(context$n_items)) {
    out[, item_index] <- tabulate(
      context$score[context$valid_rows] + 1L,
      nbins = context$max_total_score + 1L
    )
  }
  out
}
