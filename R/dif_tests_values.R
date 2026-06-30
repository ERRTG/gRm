# Source observed IX-margin degrees of freedom
#
# Port of the `IXinfo^.df` value set by source `COUNT_MARGINS` for
# item-by-exogeneous margins: nonzero observed item-score levels minus one,
# times nonzero observed exogeneous levels minus one.
#
# @param observed_ix Integer observed item-score by exogeneous-value margin.
# @return Degrees of freedom for the included IX margin.
source_ix_observed_df <- function(observed_ix) {
  nonzero_item_scores <- sum(rowSums(observed_ix) > 0)
  nonzero_exo_values <- sum(colSums(observed_ix) > 0)
  as.integer(max(0L, nonzero_item_scores - 1L) * max(0L, nonzero_exo_values - 1L))
}

# Source Goodman-Kruskal gamma count totals
#
# @param tab Two-way integer table.
# @return Gamma, PPQ, and PMQ totals.
dif_tests_source_gamma_counts <- function(tab) {
  source_rc_gamma_counts(tab)
}

# Build cached raw-data context for DIF gamma values
#
# @param project DIGRAM project.
# @return Raw item/background matrices and score vectors.
build_dif_gamma_context <- function(project) {
  items <- project$items
  backgrounds <- project$backgrounds
  item_matrix <- matrix(0L, nrow = nrow(project$raw_data), ncol = nrow(items))
  complete_items <- rep(TRUE, nrow(item_matrix))
  for (item_index in seq_len(nrow(items))) {
    item_matrix[, item_index] <- as.integer(project$raw_data[, items$position[[item_index]]])
    complete_items <- complete_items &
      item_matrix[, item_index] >= 1L &
      item_matrix[, item_index] <= items$raw_max[[item_index]]
  }
  background_matrix <- matrix(0L, nrow = nrow(project$raw_data), ncol = nrow(backgrounds))
  for (background_index in seq_len(nrow(backgrounds))) {
    background_matrix[, background_index] <- as.integer(project$raw_data[, backgrounds$position[[background_index]]])
  }
  item_score <- rowSums(sweep(item_matrix, 2L, 1L, "-"))

  list(
    item_matrix = item_matrix,
    background_matrix = background_matrix,
    complete_items = complete_items,
    item_score = item_score,
    max_score = if (any(complete_items)) max(item_score[complete_items]) else sum(items$raw_max - 1L)
  )
}

# Calculate source-shaped DIF gamma diagnostics
#
# Source trace: `SKbias7.inexpensive_itembias1` calls
# `SKbias3.XYZ_bias_ANALYSE`, stores `RESULTS(.1,5.)` as `Part_g` and
# `RESULTS(.1,6.)` as `Part_p`, then `Item_screening` copies those matrices to
# `IJXgamma` and `IJXgamma_pvalues`. `DGRirtD.MissingDIF` reports
# `IJXgamma(.i,nitems+j.)`, and the GLLRM parameter reports label the
# same value "Item screening: Gamma".
#
# @param project DIGRAM project.
# @param context Context from [build_dif_gamma_context()].
# @param target_item One-based item index.
# @param background_index One-based background index.
# @return Source partial gamma, p-value, and accumulated PPQ/PMQ/S values.
dif_tests_partial_gamma_stats <- function(project, context, target_item, background_index) {
  items <- project$items
  backgrounds <- project$backgrounds
  item_dim <- items$raw_max[[target_item]]
  background_dim <- backgrounds$raw_max[[background_index]]
  x <- context$item_matrix[, target_item]
  y <- context$background_matrix[, background_index]
  valid <- context$complete_items &
    y >= 1L & y <= background_dim &
    context$item_score > 0L & context$item_score < context$max_score
  ppq <- 0
  pmq <- 0
  s <- 0

  for (score in sort(unique(context$item_score[valid]))) {
    in_stratum <- valid & context$item_score == score
    tab <- matrix(0L, nrow = item_dim, ncol = background_dim)
    index <- x[in_stratum] + (y[in_stratum] - 1L) * item_dim
    tab[] <- tabulate(index, nbins = length(tab))
    stats <- screen_rc_gamma(tab)
    ppq <- ppq + stats$ppq
    pmq <- pmq + stats$pmq
    s <- s + stats$s
  }

  if (ppq <= 0) {
    return(list(gamma = 0, p_value = 2, ppq = 0, pmq = 0, s = 0))
  }

  gamma <- pmq / ppq
  s <- s / ppq / ppq
  u <- if (s <= 0) 4 else abs(gamma / sqrt(s))
  p_value <- 2 * if (u > 4) 0.000001 else source_tail_norm(u, TRUE)
  list(gamma = gamma, p_value = p_value, ppq = ppq, pmq = pmq, s = s)
}

# @return Source partial gamma.
dif_tests_partial_gamma <- function(project, context, target_item, background_index) {
  dif_tests_partial_gamma_stats(project, context, target_item, background_index)$gamma
}

# Build a fast source-shaped DIF estimation context
#
# @param bundle Source-shaped bundle.
# @return Cached integer matrices and source bounds for included-DIF fitting.
# Source trace: source/PAS_scd/DGRirtD.pas::CHECK D prepares each candidate DIF
# term before fitting the included IX parameter. The R context holds the same
# item, score, and exogenous margins needed by the candidate calculation.
build_candidate_dif_context <- function(bundle) {
  items <- bundle$model$items
  backgrounds <- bundle$model$backgrounds
  data <- bundle$data

  item_matrix <- matrix(0L, nrow = nrow(data), ncol = nrow(items))
  for (item_index in seq_len(nrow(items))) {
    item_matrix[, item_index] <- as.integer(data[[items$name[[item_index]]]])
  }

  background_matrix <- matrix(0L, nrow = nrow(data), ncol = nrow(backgrounds))
  for (background_index in seq_len(nrow(backgrounds))) {
    background_matrix[, background_index] <- as.integer(data[[backgrounds$name[[background_index]]]])
  }

  list(
    n_items = nrow(items),
    item_raw_max = as.integer(items$raw_max),
    background_raw_max = as.integer(backgrounds$raw_max),
    max_item_raw_max = max(as.integer(items$raw_max)),
    max_total_score = as.integer(bundle$model$max_total_score),
    item_score_columns = lapply(as.integer(items$raw_max), seq_len),
    item_score_values = lapply(as.integer(items$raw_max), function(raw_max) seq.int(0L, raw_max - 1L)),
    item_matrix = item_matrix,
    background_matrix = background_matrix,
    score = as.integer(data$score),
    valid_rows = which(data$status == 1L)
  )
}

# Rescale item gamma parameters using cached source bounds
#
# @param context Included DIF context from [build_candidate_dif_context()].
# @param item_gamma Matrix of item gamma parameters after an IPF update.
# @return Rescaled item gamma matrix.
adjust_item_gammas_source_scale_context <- function(context, item_gamma) {
  last_sgamma <- 1
  s_max <- 0L
  s_min <- 0L

  for (item_index in seq_len(context$n_items)) {
    fra <- 0L
    til <- context$item_raw_max[[item_index]] - 1L
    s_max <- s_max + til
    s_min <- s_min + fra
    alpha <- item_gamma[item_index, fra + 1L]
    if (alpha > 0) {
      for (score in seq.int(fra, til)) {
        item_gamma[item_index, score + 1L] <- item_gamma[item_index, score + 1L] / alpha
      }
    }
    if (item_gamma[item_index, til + 1L] > 0) {
      last_sgamma <- last_sgamma * item_gamma[item_index, til + 1L]
    }
  }

  alpha <- 0
  if ((s_max - s_min) > 0) {
    alpha <- -log(last_sgamma) / (s_max - s_min)
  }

  for (item_index in seq_len(context$n_items)) {
    til <- context$item_raw_max[[item_index]] - 1L
    for (score in seq.int(0L, til)) {
      item_gamma[item_index, score + 1L] <- exp(score * alpha) * item_gamma[item_index, score + 1L]
    }
  }

  item_gamma
}

# Build an included-DIF score gamma table
#
# @param bundle Source-shaped bundle.
# @param item_gamma Item gamma matrix.
# @param dif_gamma Candidate item-by-background gamma matrix.
# @param target_item One-based target item index.
# @param background_value One-based background value.
# @param excluded_item Optional one-based item index to exclude.
# @return Numeric score gamma vector indexed by score plus one.
build_candidate_dif_score_gamma <- function(bundle, item_gamma, dif_gamma, target_item, background_value, excluded_item = NA_integer_) {
  items <- bundle$model$items
  max_total_score <- bundle$model$max_total_score
  gamma_values <- numeric(max_total_score + 1L)
  gamma_values[[1L]] <- 1
  current_max <- 0L

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
        # Source trace: included DIF multiplies the target item score gamma by
        # the item-by-exogenous cell parameter for the person's background.
        weight <- item_gamma[item_index, as.character(item_score)]
        if (!is.na(target_item) && item_index == target_item) {
          weight <- weight * dif_gamma[item_score + 1L, background_value]
        }
        next_values[[score + item_score + 1L]] <-
          next_values[[score + item_score + 1L]] + current_weight * weight
      }
    }
    gamma_values[seq_len(next_max + 1L)] <- next_values[seq_len(next_max + 1L)]
    current_max <- next_max
  }

  gamma_values
}

# Convolve two score polynomials with source score bounds
#
# @param left Numeric score vector indexed by score plus one.
# @param right Numeric score vector indexed by score plus one.
# @param max_total_score Maximum retained score.
# @return Numeric score vector.
convolve_score_vectors <- function(left, right, max_total_score) {
  out <- numeric(max_total_score + 1L)
  left_len <- min(length(left), max_total_score + 1L)
  right_len <- min(length(right), max_total_score + 1L)
  right_index <- seq_len(right_len)
  for (li in seq_len(left_len)) {
    left_weight <- left[[li]]
    if (left_weight == 0) {
      next
    }
    left_score <- li - 1L
    usable_right_len <- min(right_len, max_total_score - left_score + 1L)
    if (usable_right_len <= 0L) {
      next
    }
    out_index <- seq.int(li, length.out = usable_right_len)
    out[out_index] <- out[out_index] + left_weight * right[right_index[seq_len(usable_right_len)]]
  }
  out
}

# Count observed margins for one included DIF candidate using cached data
#
# @param bundle Source-shaped bundle.
# @param context Included DIF context from [build_candidate_dif_context()].
# @param background_index One-based background index.
# @param target_item One-based item index.
# @return Observed score/background and target DIF margins.
candidate_dif_counts_context <- function(bundle, context, background_index, target_item) {
  max_total_score <- context$max_total_score
  background_max <- context$background_raw_max[[background_index]]
  item_max <- context$item_raw_max[[target_item]]

  score_background_counts <- matrix(
    0L,
    nrow = max_total_score + 1L,
    ncol = background_max,
    dimnames = list(as.character(seq.int(0L, max_total_score)), as.character(seq_len(background_max)))
  )
  observed_dif <- matrix(
    0L,
    nrow = item_max,
    ncol = background_max,
    dimnames = list(as.character(seq.int(0L, item_max - 1L)), as.character(seq_len(background_max)))
  )

  scores <- context$score[context$valid_rows]
  background_values <- context$background_matrix[context$valid_rows, background_index]
  target_scores <- context$item_matrix[context$valid_rows, target_item]

  score_linear <- scores + (background_values - 1L) * (max_total_score + 1L) + 1L
  score_background_counts[] <- tabulate(score_linear, nbins = length(score_background_counts))

  dif_linear <- target_scores + (background_values - 1L) * item_max + 1L
  observed_dif[] <- tabulate(dif_linear, nbins = length(observed_dif))

  list(
    score_background_counts = score_background_counts,
    observed_dif = observed_dif
  )
}

# Calculate expected margins for one included DIF candidate using cached bounds
#
# @param context Included DIF context from [build_candidate_dif_context()].
# @param candidate_counts Count object from [candidate_dif_counts_context()].
# @param item_gamma Item gamma matrix.
# @param dif_gamma Candidate DIF gamma matrix.
# @param target_item One-based target item index.
# @return Expected item and DIF margins.
# Source trace: source/GLLRM_ESTIM.txt::Find_new_IXparameters computes fitted
# IX margins for the included DIF candidate before applying the observed/fitted
# update. This R helper calculates those fitted candidate margins directly from
# the current item and DIF gammas.
calculate_candidate_dif_expected_context <- function(context, candidate_counts, item_gamma, dif_gamma, target_item) {
  background_max <- ncol(candidate_counts$score_background_counts)
  max_total_score <- context$max_total_score
  expected_items <- item_gamma
  expected_items[,] <- 0
  expected_dif <- dif_gamma
  expected_dif[,] <- 0

  full_gamma <- vector("list", background_max)
  without_gamma <- vector("list", background_max)
  for (background_value in seq_len(background_max)) {
    item_weights <- vector("list", context$n_items)
    for (item_index in seq_len(context$n_items)) {
      item_cols <- context$item_score_columns[[item_index]]
      weights <- as.numeric(item_gamma[item_index, item_cols])
      if (item_index == target_item) {
        weights <- weights * dif_gamma[item_cols, background_value]
      }
      item_weights[[item_index]] <- weights
    }

    prefix <- vector("list", context$n_items + 1L)
    suffix <- vector("list", context$n_items + 1L)
    prefix[[1L]] <- c(1, numeric(max_total_score))
    for (item_index in seq_len(context$n_items)) {
      prefix[[item_index + 1L]] <- convolve_score_vectors(
        prefix[[item_index]],
        item_weights[[item_index]],
        max_total_score
      )
    }
    suffix[[context$n_items + 1L]] <- c(1, numeric(max_total_score))
    for (item_index in rev(seq_len(context$n_items))) {
      suffix[[item_index]] <- convolve_score_vectors(
        item_weights[[item_index]],
        suffix[[item_index + 1L]],
        max_total_score
      )
    }
    full_gamma[[background_value]] <- prefix[[context$n_items + 1L]]
    without_gamma[[background_value]] <- vector("list", context$n_items)
    for (item_index in seq_len(context$n_items)) {
      without_gamma[[background_value]][[item_index]] <- convolve_score_vectors(
        prefix[[item_index]],
        suffix[[item_index + 1L]],
        max_total_score
      )
    }
  }

  for (background_value in seq_len(background_max)) {
    score_counts <- candidate_counts$score_background_counts[, background_value]
    used_score_indices <- which(score_counts > 0)
    if (length(used_score_indices) == 0L) {
      next
    }
    denominator <- full_gamma[[background_value]][used_score_indices]
    usable <- denominator > 0
    if (!any(usable)) {
      next
    }
    used_score_indices <- used_score_indices[usable]
    scores <- used_score_indices - 1L
    score_weight <- score_counts[used_score_indices] / denominator[usable]

    for (item_index in seq_len(context$n_items)) {
      rest_gamma <- without_gamma[[background_value]][[item_index]]
      for (item_score in context$item_score_values[[item_index]]) {
        score_usable <- scores >= item_score
        if (!any(score_usable)) {
          next
        }
        item_col <- item_score + 1L
        weight <- item_gamma[item_index, item_col]
        if (item_index == target_item) {
          weight <- weight * dif_gamma[item_col, background_value]
        }
        expected <- weight * sum(
          score_weight[score_usable] *
            rest_gamma[scores[score_usable] - item_score + 1L]
        )
        expected_items[item_index, item_col] <- expected_items[item_index, item_col] + expected
        if (item_index == target_item) {
          expected_dif[item_col, background_value] <- expected_dif[item_col, background_value] + expected
        }
      }
    }
  }

  list(expected_items = expected_items, expected_dif = expected_dif)
}

# Fit one included-DIF candidate model
#
# @param bundle Source-shaped bundle.
# @param base_counts Base Rasch count object.
# @param target_item One-based target item index.
# @param background_index One-based background index.
# @param max_step Maximum IPF iterations.
# @param max_delta Convergence threshold.
# @return Fitted candidate model state.
# Source trace: source/GLLRM_ESTIM.txt::Find_new_IXparameters and
# source/GLLRM_ESTIM.txt::Adjust_IXparameters define the candidate DIF update
# and final reporting gauge. The R function fits one candidate IX term for
# CHECK D without changing the surrounding current GLLRM.
fit_dif_candidate <- function(bundle, base_counts, target_item, background_index, max_step = 5000L, max_delta = 0.0001, initial_item_gamma_matrix = NULL, context = NULL, candidate_counts = NULL, progress = FALSE, progress_label = NULL, progress_interval = 100L, output_base_loglike = NULL, output_df = NULL, output_check_interval = 100L, output_stable_checks = 3L, output_max_delta = 0.01) {
  items <- bundle$model$items
  backgrounds <- bundle$model$backgrounds
  if (is.null(context)) {
    context <- build_candidate_dif_context(bundle)
  }
  item_gamma <- if (is.null(initial_item_gamma_matrix)) initial_item_gamma(bundle) else initial_item_gamma_matrix
  item_max <- items$raw_max[[target_item]]
  background_max <- backgrounds$raw_max[[background_index]]
  dif_gamma <- matrix(
    1,
    nrow = item_max,
    ncol = background_max,
    dimnames = list(as.character(seq.int(0L, item_max - 1L)), as.character(seq_len(background_max)))
  )
  if (is.null(candidate_counts)) {
    candidate_counts <- candidate_dif_counts_context(bundle, context, background_index, target_item)
  }
  delta <- 0
  converged <- FALSE
  n_step <- 0L
  started_at <- Sys.time()
  output_keys <- character()
  output_stable <- FALSE
  output_check_enabled <- !is.null(output_base_loglike) &&
    !is.null(output_df) &&
    output_check_interval > 0L &&
    output_stable_checks > 0L

  for (step in seq_len(max_step)) {
    n_step <- step
    expected <- calculate_candidate_dif_expected_context(
      context, candidate_counts, item_gamma, dif_gamma, target_item
    )
    delta <- 0

    for (item_index in seq_len(context$n_items)) {
      for (item_score in seq.int(0L, context$item_raw_max[[item_index]] - 1L)) {
        item_col <- item_score + 1L
        observed <- base_counts$item_counts[item_index, item_col]
        fitted <- expected$expected_items[item_index, item_col]
        if (observed > 0 && fitted > 0) {
          # Source trace: IPF update for item sufficient-statistic cells.
          ratio <- observed / fitted
          item_gamma[item_index, item_col] <- item_gamma[item_index, item_col] * ratio
          delta <- max(delta, abs(fitted - observed))
        } else if (observed == 0) {
          item_gamma[item_index, item_col] <- 0
        }
      }
    }

    for (item_score in seq.int(0L, item_max - 1L)) {
      for (background_value in seq_len(background_max)) {
        observed <- candidate_counts$observed_dif[item_score + 1L, background_value]
        fitted <- expected$expected_dif[item_score + 1L, background_value]
        if (observed > 0 && fitted > 0) {
          # Source trace: IPF update for the included IX/DIF margin.
          ratio <- observed / fitted
          dif_gamma[item_score + 1L, background_value] <-
            dif_gamma[item_score + 1L, background_value] * ratio
          delta <- max(delta, abs(fitted - observed))
        } else if (observed == 0) {
          dif_gamma[item_score + 1L, background_value] <- 0
        }
      }
    }

    dif_gamma <- adjust_gllrm_dif_reference(candidate_counts$observed_dif, dif_gamma)
    item_gamma <- adjust_item_gammas_source_scale(bundle, item_gamma)
    if (isTRUE(progress) && progress_interval > 0L && (step == 1L || (step %% progress_interval) == 0L)) {
      message(
        sprintf(
          "DIF fit step  %s pid=%s step=%d delta=%.12g elapsed=%.1fs",
          if (is.null(progress_label)) "" else progress_label,
          Sys.getpid(),
          step,
          delta,
          as.numeric(difftime(Sys.time(), started_at, units = "secs"))
        )
      )
    }
    if (delta <= max_delta) {
      converged <- TRUE
      break
    }
    if (output_check_enabled && step >= output_check_interval && (step %% output_check_interval) == 0L) {
      output_loglike <- candidate_dif_loglike_context(
        bundle,
        context,
        list(item_gamma = item_gamma, dif_gamma = dif_gamma),
        target_item,
        background_index
      )
      output_clr <- 2 * abs(output_base_loglike - output_loglike)
      output_key <- paste0(
        sprintf("%.2f", output_clr),
        "\r",
        sprintf("%.4f", source_pfchi(output_df, output_clr))
      )
      output_keys <- c(output_keys, output_key)
      if (length(output_keys) >= output_stable_checks &&
          delta <= output_max_delta &&
          length(unique(utils::tail(output_keys, output_stable_checks))) == 1L) {
        converged <- TRUE
        output_stable <- TRUE
        break
      }
    }
  }

  list(
    item_gamma = item_gamma,
    dif_gamma = dif_gamma,
    candidate_counts = candidate_counts,
    delta = delta,
    converged = converged,
    output_stable = output_stable,
    n_step = n_step
  )
}

# Calculate the base Rasch negative log likelihood
#
# @param bundle Source-shaped bundle.
# @param item_gamma Item gamma matrix.
# @return Negative conditional log likelihood.
base_rasch_loglike <- function(bundle, item_gamma) {
  data <- bundle$data
  items <- bundle$model$items
  valid_rows <- which(data$status == 1L)
  full_gamma <- build_candidate_dif_score_gamma(
    bundle,
    item_gamma,
    matrix(1, nrow = 1L, ncol = 1L),
    target_item = NA_integer_,
    background_value = 1L
  )
  loglike <- 0
  for (row_index in valid_rows) {
    score <- data$score[[row_index]]
    product_gamma <- 1
    for (item_index in seq_len(nrow(items))) {
      item_score <- data[[items$name[[item_index]]]][[row_index]]
      product_gamma <- product_gamma * item_gamma[item_index, as.character(item_score)]
    }
    probability <- product_gamma / full_gamma[[score + 1L]]
    if (probability > 0) {
      loglike <- loglike - log(probability)
    }
  }
  loglike
}

# Calculate included-DIF negative log likelihood using cached data
#
# @param bundle Source-shaped bundle.
# @param context Included DIF context from [build_candidate_dif_context()].
# @param fit Candidate fit from [fit_dif_candidate()].
# @param target_item One-based target item index.
# @param background_index One-based background index.
# @return Negative conditional log likelihood.
candidate_dif_loglike_context <- function(bundle, context, fit, target_item, background_index) {
  background_max <- context$background_raw_max[[background_index]]
  full_gamma <- vector("list", background_max)
  for (background_value in seq_len(background_max)) {
    full_gamma[[background_value]] <- build_candidate_dif_score_gamma(
      bundle, fit$item_gamma, fit$dif_gamma, target_item, background_value
    )
  }

  loglike <- 0
  for (row_index in context$valid_rows) {
    score <- context$score[[row_index]]
    background_value <- context$background_matrix[row_index, background_index]
    product_gamma <- 1
    for (item_index in seq_len(context$n_items)) {
      item_score <- context$item_matrix[row_index, item_index]
      product_gamma <- product_gamma * fit$item_gamma[item_index, item_score + 1L]
      if (item_index == target_item) {
        product_gamma <- product_gamma * fit$dif_gamma[item_score + 1L, background_value]
      }
    }
    probability <- product_gamma / full_gamma[[background_value]][[score + 1L]]
    if (probability > 0) {
      loglike <- loglike - log(probability)
    }
  }
  loglike
}

# Compute DIGRAM CHECK D DIF-test values
#
# Computes the item-by-exogenous-variable no-DIF report values used in
# DIGRAM's standard `check-missing-DIF.txt` report with a native R port of the
# one-included-DIF source likelihood-ratio path. The algorithm follows
# `docs/example_DIF_TESTS_SOURCE_TRACE.md`, especially `DGRirtD`'s CHECK D branch,
# `Estimate_GLLRM`, `PFCHI`, and the source Benjamini-Hochberg rule.
#
# @param project A source-shaped DIGRAM project list, such as the `project`
#   component returned by [gRm()] or [read_digram_project()].
# @param max_step Maximum source Rasch/GLLRM estimation iterations.
# @param max_delta Convergence threshold for source Rasch/GLLRM estimation.
# @param jobs Upper bound for parallel candidate DIF fits.
# @param progress Logical; if `TRUE`, emit temporary per-candidate progress
#   messages for long-running validation/debugging runs.
# @return A `gRm_dif_tests_values` object with test rows and BH threshold.
dif_tests_values <- function(project, max_step = 5000L, max_delta = 0.0001, jobs = min(8L, parallel::detectCores(logical = TRUE), 128L), progress = identical(Sys.getenv("RDIGRAM_DIF_PROGRESS"), "1")) {
  if (inherits(project, "gRm_fit") && inherits(project$values, "gRm_gllrm_values")) {
    return(gllrm_dif_tests_values(project, max_step = max_step, max_delta = max_delta, jobs = jobs))
  }
  if (!is.list(project)) {
    stop("`project` must be a DIGRAM project list.", call. = FALSE)
  }
  if (is.null(project$backgrounds) || nrow(project$backgrounds) == 0L) {
    stop("DIF tests require at least one exogeneous variable.", call. = FALSE)
  }

  bundle <- build_item_parameters_bundle(project)
  context <- build_candidate_dif_context(bundle)
  gamma_context <- build_dif_gamma_context(project)
  base_fit <- fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)
  base_counts <- base_fit$counts
  base_loglike <- base_rasch_loglike(bundle, base_fit$item_gamma)
  items <- project$items
  backgrounds <- project$backgrounds

  candidates <- expand.grid(
    item_index = seq_len(nrow(items)),
    background_index = seq_len(nrow(backgrounds))
  )
  candidates <- candidates[order(candidates$background_index, candidates$item_index), , drop = FALSE]
  candidate_counts <- vector("list", nrow(candidates))
  for (candidate_row in seq_len(nrow(candidates))) {
    candidate_counts[[candidate_row]] <- candidate_dif_counts_context(
      bundle,
      context,
      candidates$background_index[[candidate_row]],
      candidates$item_index[[candidate_row]]
    )
  }

  fit_one <- function(candidate_row) {
    item_index <- candidates$item_index[[candidate_row]]
    background_index <- candidates$background_index[[candidate_row]]
    df <- source_ix_observed_df(candidate_counts[[candidate_row]]$observed_dif)
    started_at <- Sys.time()
    if (isTRUE(progress)) {
      message(
        sprintf(
          "DIF fit start [%03d/%03d] pid=%s item=%s background=%s",
          candidate_row,
          nrow(candidates),
          Sys.getpid(),
          items$name[[item_index]],
          backgrounds$name[[background_index]]
        )
      )
    }
    candidate_fit <- fit_dif_candidate(
      bundle, base_counts, item_index, background_index,
      max_step = max_step, max_delta = max_delta,
      initial_item_gamma_matrix = base_fit$item_gamma,
      context = context,
      candidate_counts = candidate_counts[[candidate_row]],
      progress = progress,
      progress_label = sprintf(
        "[%03d/%03d] item=%s background=%s",
        candidate_row,
        nrow(candidates),
        items$name[[item_index]],
        backgrounds$name[[background_index]]
      ),
      progress_interval = as.integer(Sys.getenv("RDIGRAM_DIF_PROGRESS_INTERVAL", "100")),
      output_base_loglike = base_loglike,
      output_df = df
    )
    candidate_loglike <- candidate_dif_loglike_context(bundle, context, candidate_fit, item_index, background_index)
    if (isTRUE(progress)) {
      message(
        sprintf(
          "DIF fit done  [%03d/%03d] pid=%s item=%s background=%s elapsed=%.1fs step=%d converged=%s delta=%.12g",
          candidate_row,
          nrow(candidates),
          Sys.getpid(),
          items$name[[item_index]],
          backgrounds$name[[background_index]],
          as.numeric(difftime(Sys.time(), started_at, units = "secs")),
          candidate_fit$n_step,
          if (isTRUE(candidate_fit$converged)) "TRUE" else "FALSE",
          candidate_fit$delta
        )
      )
    }

    # Source trace: DGRirtD CHECK D reports
    # CLR = 2 * abs(Raschloglike - Raschloglike1).
    clr <- 2 * abs(base_loglike - candidate_loglike)

    data.frame(
      item_label = items$label_code[[item_index]],
      background_label = backgrounds$label_code[[background_index]],
      item_name = items$name[[item_index]],
      background_name = backgrounds$name[[background_index]],
      chi_square = clr,
      degrees_of_freedom = df,
      # Source trace: p = PFCHI(df, clr).
      p_value = source_pfchi(df, clr),
      gamma = dif_tests_partial_gamma(project, gamma_context, item_index, background_index),
      converged = isTRUE(candidate_fit$converged),
      output_stable = isTRUE(candidate_fit$output_stable),
      delta = candidate_fit$delta,
      n_step = candidate_fit$n_step,
      stringsAsFactors = FALSE
    )
  }

  tests <- source_candidate_map(nrow(candidates), jobs, fit_one)

  result <- structure(
    list(
      tests = tests,
      bh_critical_p = source_bh_critical(tests$p_value, 0.05),
      max_step = as.integer(max_step),
      max_delta = max_delta
    ),
    class = "gRm_dif_tests_values"
  )
  result
}

gllrm_dif_tests_values <- function(fit, max_step = 5000L, max_delta = 0.0001, jobs = min(8L, parallel::detectCores(logical = TRUE), 128L)) {
  context <- fit$fit$context
  if (context$n_backgrounds == 0L) {
    stop("DIF tests require at least one exogeneous variable.", call. = FALSE)
  }

  candidates <- gllrm_no_dif_candidates(context)
  if (nrow(candidates) == 0L) {
    empty <- data.frame(
      item_label = character(),
      background_label = character(),
      item_name = character(),
      background_name = character(),
      chi_square = numeric(),
      degrees_of_freedom = integer(),
      p_value = numeric(),
      gamma = numeric(),
      p_chi = numeric(),
      p_gamma = numeric(),
      gamma_source = character(),
      test_type = character(),
      status = character(),
      converged = logical(),
      delta = numeric(),
      n_step = integer(),
      stop_reason = character(),
      attempted_n_step = integer(),
      attempted_delta = numeric(),
      attempted_converged = logical(),
      attempted_stop_reason = character(),
      reported_checkpoint_step = integer(),
      report_value_source = character(),
      stringsAsFactors = FALSE
    )
    return(structure(
      list(
        tests = empty,
        included_tests = gllrm_dif_tests(context),
        bh_critical_p = 0,
        max_step = as.integer(max_step),
        max_delta = max_delta
      ),
      class = "gRm_dif_tests_values"
    ))
  }

  base_loglike <- fit$fit$log_likelihood
  project <- fit$project
  gamma_context <- build_dif_gamma_context(project)
  fit_one <- function(candidate_row) {
    item <- candidates$item_index[[candidate_row]]
    background <- candidates$background_index[[candidate_row]]
    candidate_fit <- fit_gllrm_with_added_dif(
      fit,
      item = item,
      background = background,
      max_step = max_step,
      max_delta = max_delta
    )
    candidate_loglike <- candidate_fit$log_likelihood
    clr <- 2 * abs(base_loglike - candidate_loglike)
    dif_index <- gllrm_context_dif_index(candidate_fit$context, item, background)
    df <- source_ix_observed_df(
      candidate_fit$context$observed_dif[[dif_index]]
    )
    # Source trace: DGRirtD.MissingDIF prints IJXgamma(.i,nitems+j.) from
    # item screening, not a statistic recomputed from the candidate
    # refit. SKbias12b/skbias22 use the same matrix and label it
    # "Item screening: Gamma" in GLLRM parameter output.
    gamma_stats <- dif_tests_partial_gamma_stats(project, gamma_context, item, background)
    data.frame(
      item_label = context$items$label_code[[item]],
      background_label = context$backgrounds$label_code[[background]],
      item_name = context$items$name[[item]],
      background_name = context$backgrounds$name[[background]],
      chi_square = clr,
      degrees_of_freedom = df,
      p_value = source_pfchi(df, clr),
      p_chi = source_pfchi(df, clr),
      gamma = gamma_stats$gamma,
      p_gamma = gamma_stats$p_value,
      gamma_source = "item_screening",
      test_type = "no_dif",
      status = "tested",
      converged = isTRUE(candidate_fit$converged),
      delta = candidate_fit$report_delta %||% NA_real_,
      n_step = candidate_fit$n_step %||% NA_integer_,
      stop_reason = candidate_fit$stop_reason %||% NA_character_,
      attempted_n_step = candidate_fit$n_step %||% NA_integer_,
      attempted_delta = candidate_fit$report_delta %||% candidate_fit$delta %||% NA_real_,
      attempted_converged = candidate_fit$converged %||% NA,
      attempted_stop_reason = candidate_fit$stop_reason %||% NA_character_,
      reported_checkpoint_step = NA_integer_,
      report_value_source = "attempted_fit",
      stringsAsFactors = FALSE
    )
  }

  tests <- source_candidate_map(nrow(candidates), jobs, fit_one)
  structure(
    list(
      tests = tests,
      included_tests = gllrm_dif_tests(context),
      bh_critical_p = source_bh_critical(tests$p_value, 0.05),
      max_step = as.integer(max_step),
      max_delta = max_delta
    ),
    class = "gRm_dif_tests_values"
  )
}

gllrm_context_dif_index <- function(context, item, background) {
  key <- gllrm_dif_key(item, background)
  hit <- which(vapply(context$dif_specs, function(spec) {
    identical(gllrm_dif_key(spec$item, spec$background), key)
  }, logical(1L)))
  if (length(hit) != 1L) {
    stop("Could not identify the candidate DIF term in the fitted GLLRM context.", call. = FALSE)
  }
  hit[[1L]]
}

gllrm_no_dif_candidates <- function(context) {
  included <- gllrm_dif_lookup(context)
  rows <- expand.grid(
    item_index = seq_len(context$n_items),
    background_index = seq_len(context$n_backgrounds)
  )
  rows <- rows[order(rows$background_index, rows$item_index), , drop = FALSE]
  keep <- !vapply(
    seq_len(nrow(rows)),
    function(i) included[[gllrm_dif_key(rows$item_index[[i]], rows$background_index[[i]])]] %||% FALSE,
    logical(1L)
  )
  rows[keep, , drop = FALSE]
}

gllrm_dif_tests <- function(context) {
  if (length(context$dif_specs) == 0L) {
    return(data.frame())
  }
  rows <- lapply(context$dif_specs, function(spec) {
    test <- gllrm_dif_pair_test(context, spec$item, spec$background)
    data.frame(
      item_label = context$items$label_code[[spec$item]],
      background_label = context$backgrounds$label_code[[spec$background]],
      item_name = context$items$name[[spec$item]],
      background_name = context$backgrounds$name[[spec$background]],
      chi_square = test$chi_square,
      degrees_of_freedom = test$df,
      p_value = test$p_chi,
      gamma = test$gamma,
      p_chi = test$p_chi,
      p_gamma = test$p_gamma,
      test_type = "included_dif",
      status = "tested",
      converged = TRUE,
      delta = 0,
      n_step = 0L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

gllrm_dif_pair_test <- function(context, item, background) {
  condition_values <- cbind(context$score)
  condition_dims <- context$max_total_score - 1L
  valid <- seq_len(nrow(context$item_matrix)) %in% context$valid_rows
  valid <- valid & context$score > 0L & context$score < context$max_total_score
  screen_j_conditional_bias_test(
    x = context$item_matrix[, item] + 1L,
    y = context$background_matrix[, background],
    x_dim = context$item_raw_max[[item]],
    y_dim = context$background_raw_max[[background]],
    condition_values = condition_values,
    condition_dims = condition_dims,
    valid = valid,
    exact = FALSE,
    native = FALSE
  )
}

gllrm_dif_lookup <- function(context) {
  out <- list()
  for (spec in context$dif_specs) {
    out[[gllrm_dif_key(spec$item, spec$background)]] <- TRUE
  }
  out
}

gllrm_dif_key <- function(item, background) {
  paste(item, background, sep = ":")
}
