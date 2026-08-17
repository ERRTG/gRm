#' Build DIGRAM score groups from source score cuts
#'
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param bundle Source-shaped bundle from `build_item_parameters_bundle()`.
#' @param score_cuts Integer vector of upper score cuts.
#' @return A data frame with retained score-group bounds.
#' @keywords internal
#' @noRd
global_homogeneity_score_groups <- function(bundle, score_cuts) {
  # Source trace: DGRirtD.pas global homogeneity loop after scorecuts[0] :=
  # minscore - 1. The retained group bounds are clipped to LeastScore/LargestScore
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
    stop("No global homogeneity score groups remain after source clipping.", call. = FALSE)
  }
  do.call(rbind, rows)
}

#' Internal global homogeneity score group lookup helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param groups Internal `groups` value used by this helper.
#' @param max_score Internal `max_score` value used by this helper.
#' @return The internal `global_homogeneity_score_group_lookup()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal global homogeneity uniform score group lookup helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param groups Internal `groups` value used by this helper.
#' @param max_score Internal `max_score` value used by this helper.
#' @return The internal `global_homogeneity_uniform_score_group_lookup()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_uniform_score_group_lookup <- function(groups, max_score) {
  # Source trace: the extended uniform LD/DIF block calls
  # Count_IJX_tabel/Count_IXZ_tabel before the subgroup refits. Those table
  # counters use ScoreGruppe(score, scoredim, minscore, maxscore, scorecuts),
  # where cutpoint-defined score groups keep minscore = 0. The subgroup
  # homogeneity refits still clip displayed groups to LeastScore = 1.
  source_groups <- groups
  if (nrow(source_groups) > 0L) {
    source_groups$from_score[[1L]] <- 0L
  }
  global_homogeneity_score_group_lookup(source_groups, max_score)
}

#' Internal global homogeneity lookup score helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param lookup Internal `lookup` value used by this helper.
#' @param score Zero-based total or item score.
#' @return The internal `global_homogeneity_lookup_score()` computation result.
#' @keywords internal
#' @noRd
global_homogeneity_lookup_score <- function(lookup, score) {
  if (score < 0L || score >= length(lookup)) {
    return(NA_integer_)
  }
  lookup[[score + 1L]]
}

#' Internal gllrm uniform complete rows helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param context Prepared GLLRM computation context.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @return The internal `gllrm_uniform_complete_rows()` computation result.
#' @keywords internal
#' @noRd
gllrm_uniform_complete_rows <- function(context, score_group_lookup) {
  item_complete <- rowSums(context$item_matrix < 0L) == 0L
  background_complete <- if (context$n_backgrounds > 0L) {
    rowSums(context$background_matrix < 1L) == 0L
  } else {
    rep(TRUE, nrow(context$item_matrix))
  }
  rows <- which(item_complete & background_complete)
  if (!length(rows)) {
    return(integer())
  }
  in_group <- vapply(
    context$score[rows],
    function(score) !is.na(global_homogeneity_lookup_score(score_group_lookup, score)),
    logical(1L)
  )
  rows[in_group]
}

#' Internal gRm default global homogeneity score cuts helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param project Encoded gRm project.
#' @return The internal `gRm_default_global_homogeneity_score_cuts()` computation result.
#' @keywords internal
#' @noRd
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
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param bundle Source-shaped bundle.
#' @param from_score Lower inclusive score bound.
#' @param to_score Upper inclusive score bound.
#' @return Bundle with rows outside the interval marked invalid.
#' @keywords internal
#' @noRd
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
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param bundle Score-group bundle.
#' @param fit Score-group Rasch fit.
#' @param group One-based score-group index.
#' @return Data frame of observed/expected item mean rows. The `residual`,
#'   hidden expected-variance, and marker cells are `NA`: the available Pascal
#'   source identifies the `skbias15.pas` item-mean residual path, but the
#'   historical runtime residual boundary has not been reproduced without
#'   empirical report-specific correction factors.
#' @param expected_item_gamma Internal `expected_item_gamma` value used by this helper.
#' @keywords internal
#' @noRd
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

#' Build Pascal-shaped expected item margin tables for global homogeneity
#'
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param bundle Score-group bundle.
#' @param counts Count list from `rasch_counts()`.
#' @param item_gamma Item gamma matrix used for expected margins.
#' @return List of score-by-item-score expected margin matrices.
#' @keywords internal
#' @noRd
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

#' Build a score gamma array using the skbias12 Inexpensive_Gamma_Calculation order
#'
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param bundle Score-group bundle.
#' @param item_gamma Item gamma matrix used for expected margins.
#' @param use_items Logical vector selecting items included in the score gamma.
#' @return Numeric vector indexed by total score plus one.
#' @keywords internal
#' @noRd
global_homogeneity_score_gamma <- function(bundle, item_gamma, use_items) {
  build_source_score_gamma(bundle, item_gamma, use_items)
}

#' Run skbias15.pas SummarizeTal on item-score cells
#'
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param cells Numeric expected item-score cells.
#' @param item_max Maximum item score for the item.
#' @return List containing Pascal `tal[a]`, `tal[b]`, and `tal[c]`.
#' @keywords internal
#' @noRd
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
