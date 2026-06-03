#' Source normal tail probability
#'
#' Port of `SkStat.TAILNORM`, used by the source chi-square tail approximation.
#'
#' @param value Normal deviate.
#' @param upper Logical; `TRUE` returns the upper tail.
#' @return Tail probability.
#' @keywords internal
source_tail_norm <- function(value, upper = TRUE) {
  normal_c <- 0.3989422804014327
  big_x <- 170
  if (value == 0) {
    return(0.5)
  }
  if (value < 0) {
    upper <- !upper
  }
  value <- abs(value)
  x2 <- value * value
  if ((x2 / 2) < big_x) {
    y <- normal_c * exp(-0.5 * x2)
  } else {
    y <- 0
  }
  n <- y / value
  if (!upper && n == 0) {
    return(1)
  }
  if (upper && n == 0) {
    return(0)
  }
  if ((upper && value > 2.32) || (!upper && value > 3.5)) {
    q1 <- value
    p2 <- y * value
    i <- 1
    p1 <- y
    q2 <- x2 + i
    if (upper) {
      s <- p1 / q1
      m <- s
      t <- p2 / q2
    } else {
      s <- i - p1 / q1
      m <- s
      t <- i - p2 / q2
    }
    while (m != t && s != t) {
      i <- i + 1
      s <- value * p2 + i * p1
      p1 <- p2
      p2 <- s
      s <- value * q2 + i * q1
      q1 <- q2
      q2 <- s
      s <- m
      m <- t
      if (upper) {
        t <- p2 / q2
      } else {
        t <- 1 - p2 / q2
      }
    }
    return(t)
  }

  s <- y * value
  term <- y * value
  i <- 1
  t <- 0
  while (s != t) {
    i <- i + 2
    t <- s
    term <- term * x2 / i
    s <- s + term
  }
  if (upper) {
    0.5 - s
  } else {
    0.5 + s
  }
}

#' Source chi-square upper-tail probability
#'
#' Port of `SkStat.PFCHI`.
#'
#' @param df Degrees of freedom.
#' @param x Chi-square statistic.
#' @return Upper-tail probability.
#' @keywords internal
source_pfchi <- function(df, x) {
  source_co <- 0.3989422804014327
  big_x <- 170
  df_div_2 <- df %/% 2L
  if (df == 0L || x <= 0) {
    return(1)
  }
  if (df > 100L) {
    transformed <- sqrt(4.5 * df) * (exp(log(x / df) / 3) + 1 / (4.5 * df) - 1)
    return(source_tail_norm(transformed, TRUE))
  }
  if ((df %% 2L) == 0L) {
    if (x < big_x) {
      p <- exp(-0.5 * x)
    } else {
      p <- 0
    }
    c_value <- p
    if (df_div_2 > 1L) {
      for (i in seq_len(df_div_2 - 1L)) {
        c_value <- c_value * 0.5 * x / i
        p <- p + c_value
      }
    }
    return(p)
  }

  p <- 2 * source_tail_norm(sqrt(x), TRUE)
  if (x < big_x) {
    c_value <- exp(-0.5 * x) * 2 * source_co / sqrt(x)
  } else {
    c_value <- 0
  }
  if (df_div_2 > 0L) {
    for (i in seq_len(df_div_2)) {
      c_value <- c_value * x / (2 * i - 1)
      p <- p + c_value
    }
  }
  p
}

#' Source Benjamini-Hochberg critical p-value
#'
#' Port of `SourceBenjaminiHochbergCritical`, matching the current Pascal
#' harness behavior.
#'
#' @param p_values Numeric p-values.
#' @param alpha False discovery rate target.
#' @return Critical p-value.
#' @keywords internal
source_bh_critical <- function(p_values, alpha = 0.05) {
  n <- length(p_values)
  if (n <= 0L) {
    return(0)
  }
  sorted <- sort(p_values)
  result <- alpha / n
  for (index in rev(seq_len(n))) {
    candidate <- (index / n) * alpha
    result <- candidate
    if (sorted[[index]] <= candidate) {
      break
    }
  }
  result
}

dif_tests_values_cache <- new.env(parent = emptyenv())

dif_tests_cache_key <- function(project, max_step, max_delta, jobs) {
  raw_sum <- sum(as.matrix(project$raw_data), na.rm = TRUE)
  input_dir <- if (!is.null(project$paths$input_dir)) project$paths$input_dir else ""
  paste(
    normalizePath(input_dir, mustWork = FALSE),
    nrow(project$raw_data),
    nrow(project$items),
    nrow(project$backgrounds),
    raw_sum,
    max_step,
    max_delta,
    jobs,
    sep = "\r"
  )
}

#' Source observed IX-margin degrees of freedom
#'
#' Port of the `IXinfo^.df` value set by source `COUNT_MARGINS` for
#' item-by-exogeneous margins: nonzero observed item-score levels minus one,
#' times nonzero observed exogeneous levels minus one.
#'
#' @param observed_ix Integer observed item-score by exogeneous-value margin.
#' @return Degrees of freedom for the active IX margin.
#' @keywords internal
source_ix_observed_df <- function(observed_ix) {
  nonzero_item_scores <- sum(rowSums(observed_ix) > 0)
  nonzero_exo_values <- sum(colSums(observed_ix) > 0)
  as.integer(max(0L, nonzero_item_scores - 1L) * max(0L, nonzero_exo_values - 1L))
}

#' Source Goodman-Kruskal gamma count totals
#'
#' @param tab Two-way integer table.
#' @return Gamma, PPQ, and PMQ totals.
#' @keywords internal
dif_tests_source_gamma_counts <- function(tab) {
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

#' Build cached raw-data context for DIF gamma values
#'
#' @param project DIGRAM project.
#' @return Raw item/background matrices and score vectors.
#' @keywords internal
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

#' Calculate source-shaped DIF gamma
#'
#' Ports the `SKbias7.inexpensive_itembias1` partial gamma used by
#' `DGRirtD.MissingDIF`: item by exogenous variable, stratified by total score.
#'
#' @param project DIGRAM project.
#' @param context Context from [build_dif_gamma_context()].
#' @param target_item One-based item index.
#' @param background_index One-based background index.
#' @return Source partial gamma.
#' @keywords internal
dif_tests_partial_gamma <- function(project, context, target_item, background_index) {
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

  for (score in sort(unique(context$item_score[valid]))) {
    in_stratum <- valid & context$item_score == score
    tab <- matrix(0L, nrow = item_dim, ncol = background_dim)
    index <- x[in_stratum] + (y[in_stratum] - 1L) * item_dim
    tab[] <- tabulate(index, nbins = length(tab))
    stats <- dif_tests_source_gamma_counts(tab)
    ppq <- ppq + stats$ppq
    pmq <- pmq + stats$pmq
  }

  if (ppq > 0) pmq / ppq else 0
}

#' Build a fast source-shaped DIF estimation context
#'
#' @param bundle Source-shaped bundle.
#' @return Cached integer matrices and source bounds for active-DIF fitting.
#' @keywords internal
build_active_dif_context <- function(bundle) {
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

#' Rescale item gamma parameters using cached source bounds
#'
#' @param context Active DIF context from [build_active_dif_context()].
#' @param item_gamma Matrix of item gamma parameters after an IPF update.
#' @return Rescaled item gamma matrix.
#' @keywords internal
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

#' Build an active-DIF score gamma table
#'
#' @param bundle Source-shaped bundle.
#' @param item_gamma Item gamma matrix.
#' @param dif_gamma Candidate item-by-background gamma matrix.
#' @param target_item One-based target item index.
#' @param background_value One-based background value.
#' @param excluded_item Optional one-based item index to exclude.
#' @return Numeric score gamma vector indexed by score plus one.
#' @keywords internal
build_active_dif_score_gamma <- function(bundle, item_gamma, dif_gamma, target_item, background_value, excluded_item = NA_integer_) {
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
        # Source trace: active DIF multiplies the target item score gamma by
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

#' Convolve two score polynomials with source score bounds
#'
#' @param left Numeric score vector indexed by score plus one.
#' @param right Numeric score vector indexed by score plus one.
#' @param max_total_score Maximum retained score.
#' @return Numeric score vector.
#' @keywords internal
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

#' Count observed margins for one active DIF candidate
#'
#' @param bundle Source-shaped bundle.
#' @param background_index One-based background index.
#' @param target_item One-based item index.
#' @return Observed score/background and target DIF margins.
#' @keywords internal
active_dif_counts <- function(bundle, background_index, target_item) {
  items <- bundle$model$items
  backgrounds <- bundle$model$backgrounds
  data <- bundle$data
  valid <- data$status == 1L
  max_total_score <- bundle$model$max_total_score
  background_name <- backgrounds$name[[background_index]]
  target_name <- items$name[[target_item]]
  background_max <- backgrounds$raw_max[[background_index]]
  item_max <- items$raw_max[[target_item]]

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

  valid_rows <- which(valid)
  for (row_index in valid_rows) {
    score <- data$score[[row_index]]
    background_value <- data[[background_name]][[row_index]]
    target_score <- data[[target_name]][[row_index]]
    score_background_counts[score + 1L, background_value] <-
      score_background_counts[score + 1L, background_value] + 1L
    observed_dif[target_score + 1L, background_value] <-
      observed_dif[target_score + 1L, background_value] + 1L
  }

  list(
    score_background_counts = score_background_counts,
    observed_dif = observed_dif
  )
}

#' Count observed margins for one active DIF candidate using cached data
#'
#' @param bundle Source-shaped bundle.
#' @param context Active DIF context from [build_active_dif_context()].
#' @param background_index One-based background index.
#' @param target_item One-based item index.
#' @return Observed score/background and target DIF margins.
#' @keywords internal
active_dif_counts_context <- function(bundle, context, background_index, target_item) {
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

#' Calculate expected margins for one active DIF candidate
#'
#' @param bundle Source-shaped bundle.
#' @param counts Base Rasch count object.
#' @param candidate_counts Count object from [active_dif_counts()].
#' @param item_gamma Item gamma matrix.
#' @param dif_gamma Candidate DIF gamma matrix.
#' @param target_item One-based target item index.
#' @return Expected item and DIF margins.
#' @keywords internal
calculate_active_dif_expected <- function(bundle, counts, candidate_counts, item_gamma, dif_gamma, target_item) {
  items <- bundle$model$items
  background_max <- ncol(candidate_counts$score_background_counts)
  max_total_score <- bundle$model$max_total_score
  expected_items <- item_gamma
  expected_items[,] <- 0
  expected_dif <- dif_gamma
  expected_dif[,] <- 0

  full_gamma <- vector("list", background_max)
  without_gamma <- vector("list", background_max)
  for (background_value in seq_len(background_max)) {
    item_weights <- vector("list", nrow(items))
    for (item_index in seq_len(nrow(items))) {
      scores <- seq.int(0L, items$raw_max[[item_index]] - 1L)
      weights <- as.numeric(item_gamma[item_index, as.character(scores)])
      if (item_index == target_item) {
        weights <- weights * dif_gamma[scores + 1L, background_value]
      }
      item_weights[[item_index]] <- weights
    }

    prefix <- vector("list", nrow(items) + 1L)
    suffix <- vector("list", nrow(items) + 1L)
    prefix[[1L]] <- c(1, numeric(max_total_score))
    for (item_index in seq_len(nrow(items))) {
      prefix[[item_index + 1L]] <- convolve_score_vectors(
        prefix[[item_index]],
        item_weights[[item_index]],
        max_total_score
      )
    }
    suffix[[nrow(items) + 1L]] <- c(1, numeric(max_total_score))
    for (item_index in rev(seq_len(nrow(items)))) {
      suffix[[item_index]] <- convolve_score_vectors(
        item_weights[[item_index]],
        suffix[[item_index + 1L]],
        max_total_score
      )
    }
    full_gamma[[background_value]] <- prefix[[nrow(items) + 1L]]
    without_gamma[[background_value]] <- vector("list", nrow(items))
    for (item_index in seq_len(nrow(items))) {
      without_gamma[[background_value]][[item_index]] <- convolve_score_vectors(
        prefix[[item_index]],
        suffix[[item_index + 1L]],
        max_total_score
      )
    }
  }

  for (background_value in seq_len(background_max)) {
    score_counts <- candidate_counts$score_background_counts[, background_value]
    for (score in seq.int(0L, bundle$model$max_total_score)) {
      n_score <- score_counts[[score + 1L]]
      denominator <- full_gamma[[background_value]][[score + 1L]]
      if (n_score == 0L || denominator <= 0) {
        next
      }
      for (item_index in seq_len(nrow(items))) {
        rest_gamma <- without_gamma[[background_value]][[item_index]]
        for (item_score in seq.int(0L, items$raw_max[[item_index]] - 1L)) {
          if (score < item_score) {
            next
          }
          weight <- item_gamma[item_index, as.character(item_score)]
          if (item_index == target_item) {
            weight <- weight * dif_gamma[item_score + 1L, background_value]
          }
          numerator <- weight * rest_gamma[[score - item_score + 1L]]
          expected <- n_score * numerator / denominator
          expected_items[item_index, as.character(item_score)] <-
            expected_items[item_index, as.character(item_score)] + expected
          if (item_index == target_item) {
            # Source trace: the active IX margin is the item score by the
            # selected exogenous value, summed over total-score strata.
            expected_dif[item_score + 1L, background_value] <-
              expected_dif[item_score + 1L, background_value] + expected
          }
        }
      }
    }
  }

  list(expected_items = expected_items, expected_dif = expected_dif)
}

#' Calculate expected margins for one active DIF candidate using cached bounds
#'
#' @param context Active DIF context from [build_active_dif_context()].
#' @param candidate_counts Count object from [active_dif_counts_context()].
#' @param item_gamma Item gamma matrix.
#' @param dif_gamma Candidate DIF gamma matrix.
#' @param target_item One-based target item index.
#' @return Expected item and DIF margins.
#' @keywords internal
calculate_active_dif_expected_context <- function(context, candidate_counts, item_gamma, dif_gamma, target_item) {
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
    active_score_indices <- which(score_counts > 0)
    if (length(active_score_indices) == 0L) {
      next
    }
    denominator <- full_gamma[[background_value]][active_score_indices]
    usable <- denominator > 0
    if (!any(usable)) {
      next
    }
    active_score_indices <- active_score_indices[usable]
    scores <- active_score_indices - 1L
    score_weight <- score_counts[active_score_indices] / denominator[usable]

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

#' Fit one active-DIF candidate model
#'
#' @param bundle Source-shaped bundle.
#' @param base_counts Base Rasch count object.
#' @param target_item One-based target item index.
#' @param background_index One-based background index.
#' @param max_step Maximum IPF iterations.
#' @param max_delta Convergence threshold.
#' @return Fitted candidate model state.
#' @keywords internal
fit_active_dif_candidate <- function(bundle, base_counts, target_item, background_index, max_step = 5000L, max_delta = 0.0001, initial_item_gamma_matrix = NULL, context = NULL, candidate_counts = NULL, progress = FALSE, progress_label = NULL, progress_interval = 100L, output_base_loglike = NULL, output_df = NULL, output_check_interval = 100L, output_stable_checks = 3L, output_max_delta = 0.01) {
  items <- bundle$model$items
  backgrounds <- bundle$model$backgrounds
  if (is.null(context)) {
    context <- build_active_dif_context(bundle)
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
    candidate_counts <- active_dif_counts_context(bundle, context, background_index, target_item)
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
    expected <- calculate_active_dif_expected_context(
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
          # Source trace: IPF update for the active IX/DIF margin.
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
      output_loglike <- active_dif_loglike_context(
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

#' Calculate the base Rasch negative log likelihood
#'
#' @param bundle Source-shaped bundle.
#' @param item_gamma Item gamma matrix.
#' @return Negative conditional log likelihood.
#' @keywords internal
base_rasch_loglike <- function(bundle, item_gamma) {
  data <- bundle$data
  items <- bundle$model$items
  valid_rows <- which(data$status == 1L)
  full_gamma <- build_active_dif_score_gamma(
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

#' Calculate the active-DIF negative log likelihood
#'
#' @param bundle Source-shaped bundle.
#' @param fit Candidate fit from [fit_active_dif_candidate()].
#' @param target_item One-based target item index.
#' @param background_index One-based background index.
#' @return Negative conditional log likelihood.
#' @keywords internal
active_dif_loglike <- function(bundle, fit, target_item, background_index) {
  data <- bundle$data
  items <- bundle$model$items
  backgrounds <- bundle$model$backgrounds
  background_name <- backgrounds$name[[background_index]]
  valid_rows <- which(data$status == 1L)
  full_gamma <- vector("list", backgrounds$raw_max[[background_index]])
  for (background_value in seq_len(backgrounds$raw_max[[background_index]])) {
    full_gamma[[background_value]] <- build_active_dif_score_gamma(
      bundle, fit$item_gamma, fit$dif_gamma, target_item, background_value
    )
  }

  loglike <- 0
  for (row_index in valid_rows) {
    score <- data$score[[row_index]]
    background_value <- data[[background_name]][[row_index]]
    product_gamma <- 1
    for (item_index in seq_len(nrow(items))) {
      item_score <- data[[items$name[[item_index]]]][[row_index]]
      product_gamma <- product_gamma * fit$item_gamma[item_index, as.character(item_score)]
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

#' Calculate active-DIF negative log likelihood using cached data
#'
#' @param bundle Source-shaped bundle.
#' @param context Active DIF context from [build_active_dif_context()].
#' @param fit Candidate fit from [fit_active_dif_candidate()].
#' @param target_item One-based target item index.
#' @param background_index One-based background index.
#' @return Negative conditional log likelihood.
#' @keywords internal
active_dif_loglike_context <- function(bundle, context, fit, target_item, background_index) {
  background_max <- context$background_raw_max[[background_index]]
  full_gamma <- vector("list", background_max)
  for (background_value in seq_len(background_max)) {
    full_gamma[[background_value]] <- build_active_dif_score_gamma(
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

#' Compute DIGRAM CHECK D DIF-test values
#'
#' Computes the item-by-exogenous-variable no-DIF report values used in
#' DIGRAM's standard `check-missing-DIF.txt` report with a native R port of the
#' one-active-DIF source likelihood-ratio path. The algorithm follows
#' `docs/example_DIF_TESTS_SOURCE_TRACE.md`, especially `DGRirtD`'s CHECK D branch,
#' `Estimate_GLLRM`, `PFCHI`, and the source Benjamini-Hochberg rule.
#'
#' @param project A source-shaped DIGRAM project list, such as the `project`
#'   component returned by [gRm()] or [read_digram_project()].
#' @param max_step Maximum source Rasch/GLLRM estimation iterations.
#' @param max_delta Convergence threshold for source Rasch/GLLRM estimation.
#' @param jobs Upper bound for parallel candidate DIF fits.
#' @param progress Logical; if `TRUE`, emit temporary per-candidate progress
#'   messages for long-running validation/debugging runs.
#' @return A `gRm_dif_tests_values` object with test rows and BH threshold.
#' @keywords internal
dif_tests_values <- function(project, max_step = 5000L, max_delta = 0.0001, jobs = min(8L, parallel::detectCores(logical = TRUE), 128L), progress = identical(Sys.getenv("RDIGRAM_DIF_PROGRESS"), "1")) {
  if ((inherits(project, "gRm_fit") || inherits(project, "gRm_gllrm_fit")) && inherits(project$values, "gRm_active_gllrm_values")) {
    return(active_gllrm_dif_tests_values(project, max_step = max_step, max_delta = max_delta, jobs = jobs))
  }
  if (!is.list(project)) {
    stop("`project` must be a DIGRAM project list.", call. = FALSE)
  }
  if (is.null(project$backgrounds) || nrow(project$backgrounds) == 0L) {
    stop("DIF tests require at least one exogeneous variable.", call. = FALSE)
  }
  cache_key <- dif_tests_cache_key(project, max_step, max_delta, jobs)
  if (exists(cache_key, envir = dif_tests_values_cache, inherits = FALSE)) {
    return(get(cache_key, envir = dif_tests_values_cache, inherits = FALSE))
  }

  bundle <- build_item_parameters_bundle(project)
  context <- build_active_dif_context(bundle)
  gamma_context <- build_dif_gamma_context(project)
  base_fit <- fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)
  base_counts <- base_fit$counts
  base_loglike <- base_rasch_loglike(bundle, base_fit$item_gamma)
  items <- project$items
  backgrounds <- project$backgrounds
  tests <- data.frame(
    item_label = character(nrow(items) * nrow(backgrounds)),
    background_label = character(nrow(items) * nrow(backgrounds)),
    item_name = character(nrow(items) * nrow(backgrounds)),
    background_name = character(nrow(items) * nrow(backgrounds)),
    chi_square = numeric(nrow(items) * nrow(backgrounds)),
    degrees_of_freedom = integer(nrow(items) * nrow(backgrounds)),
    p_value = numeric(nrow(items) * nrow(backgrounds)),
    gamma = numeric(nrow(items) * nrow(backgrounds)),
    converged = logical(nrow(items) * nrow(backgrounds)),
    delta = numeric(nrow(items) * nrow(backgrounds)),
    n_step = integer(nrow(items) * nrow(backgrounds)),
    stringsAsFactors = FALSE
  )

  candidates <- expand.grid(
    item_index = seq_len(nrow(items)),
    background_index = seq_len(nrow(backgrounds))
  )
  candidates <- candidates[order(candidates$background_index, candidates$item_index), , drop = FALSE]
  jobs <- max(1L, min(as.integer(jobs), 128L, nrow(candidates)))
  candidate_counts <- vector("list", nrow(candidates))
  for (candidate_row in seq_len(nrow(candidates))) {
    candidate_counts[[candidate_row]] <- active_dif_counts_context(
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
    candidate_fit <- fit_active_dif_candidate(
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
    candidate_loglike <- active_dif_loglike_context(bundle, context, candidate_fit, item_index, background_index)
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

  if (.Platform$OS.type == "unix" && jobs > 1L) {
    fit_rows <- parallel::mclapply(seq_len(nrow(candidates)), fit_one, mc.cores = jobs, mc.preschedule = FALSE)
  } else {
    fit_rows <- lapply(seq_len(nrow(candidates)), fit_one)
  }
  tests <- do.call(rbind, fit_rows)

  result <- structure(
    list(
      tests = tests,
      bh_critical_p = source_bh_critical(tests$p_value, 0.05),
      max_step = as.integer(max_step),
      max_delta = max_delta
    ),
    class = "gRm_dif_tests_values"
  )
  assign(cache_key, result, envir = dif_tests_values_cache)
  result
}

active_gllrm_dif_tests_values <- function(fit, max_step = 5000L, max_delta = 0.0001, jobs = min(8L, parallel::detectCores(logical = TRUE), 128L)) {
  context <- fit$fit$context
  if (context$n_backgrounds == 0L) {
    stop("DIF tests require at least one exogeneous variable.", call. = FALSE)
  }

  candidates <- active_gllrm_no_dif_candidates(context)
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
      test_type = character(),
      status = character(),
      converged = logical(),
      delta = numeric(),
      n_step = integer(),
      stringsAsFactors = FALSE
    )
    return(structure(
      list(
        tests = empty,
        active_tests = active_gllrm_active_dif_tests(context),
        bh_critical_p = 0,
        max_step = as.integer(max_step),
        max_delta = max_delta
      ),
      class = "gRm_dif_tests_values"
    ))
  }

  jobs <- max(1L, min(as.integer(jobs), 128L, nrow(candidates)))
  base_loglike <- fit$fit$log_likelihood
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
    candidate_loglike <- candidate_fit$fit$log_likelihood
    clr <- 2 * abs(base_loglike - candidate_loglike)
    dif_index <- active_gllrm_context_dif_index(candidate_fit$fit$context, item, background)
    df <- source_ix_observed_df(
      candidate_fit$fit$context$observed_dif[[dif_index]]
    )
    data.frame(
      item_label = context$items$label_code[[item]],
      background_label = context$backgrounds$label_code[[background]],
      item_name = context$items$name[[item]],
      background_name = context$backgrounds$name[[background]],
      chi_square = clr,
      degrees_of_freedom = df,
      p_value = source_pfchi(df, clr),
      p_chi = source_pfchi(df, clr),
      p_gamma = NA_real_,
      test_type = "no_dif",
      status = "tested",
      converged = isTRUE(candidate_fit$convergence$converged),
      delta = candidate_fit$fit$report_delta %||% NA_real_,
      n_step = candidate_fit$fit$n_step %||% NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  if (.Platform$OS.type == "unix" && jobs > 1L) {
    fit_rows <- parallel::mclapply(seq_len(nrow(candidates)), fit_one, mc.cores = jobs, mc.preschedule = FALSE)
  } else {
    fit_rows <- lapply(seq_len(nrow(candidates)), fit_one)
  }
  tests <- do.call(rbind, fit_rows)
  structure(
    list(
      tests = tests,
      active_tests = active_gllrm_active_dif_tests(context),
      bh_critical_p = source_bh_critical(tests$p_value, 0.05),
      max_step = as.integer(max_step),
      max_delta = max_delta
    ),
    class = "gRm_dif_tests_values"
  )
}

active_gllrm_context_dif_index <- function(context, item, background) {
  key <- active_gllrm_dif_key(item, background)
  hit <- which(vapply(context$dif_specs, function(spec) {
    identical(active_gllrm_dif_key(spec$item, spec$background), key)
  }, logical(1L)))
  if (length(hit) != 1L) {
    stop("Could not identify the candidate DIF term in the fitted GLLRM context.", call. = FALSE)
  }
  hit[[1L]]
}

active_gllrm_no_dif_candidates <- function(context) {
  active <- active_gllrm_dif_lookup(context)
  rows <- expand.grid(
    item_index = seq_len(context$n_items),
    background_index = seq_len(context$n_backgrounds)
  )
  rows <- rows[order(rows$background_index, rows$item_index), , drop = FALSE]
  keep <- !vapply(
    seq_len(nrow(rows)),
    function(i) active[[active_gllrm_dif_key(rows$item_index[[i]], rows$background_index[[i]])]] %||% FALSE,
    logical(1L)
  )
  rows[keep, , drop = FALSE]
}

active_gllrm_active_dif_tests <- function(context) {
  if (length(context$dif_specs) == 0L) {
    return(data.frame())
  }
  rows <- lapply(context$dif_specs, function(spec) {
    test <- active_gllrm_active_dif_pair_test(context, spec$item, spec$background)
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
      test_type = "active_dif",
      status = "tested",
      converged = TRUE,
      delta = 0,
      n_step = 0L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

active_gllrm_no_dif_pair_test <- function(context, item, background) {
  condition <- active_gllrm_no_dif_conditioning(context, item, background)
  valid <- seq_len(nrow(context$item_matrix)) %in% context$valid_rows
  valid <- valid & context$score > 0L & context$score < context$max_total_score
  screen_j_conditional_bias_test(
    x = context$item_matrix[, item] + 1L,
    y = context$background_matrix[, background],
    x_dim = context$item_raw_max[[item]],
    y_dim = context$background_raw_max[[background]],
    condition_values = condition$values,
    condition_dims = condition$dims,
    valid = valid,
    exact = FALSE,
    native = FALSE
  )
}

active_gllrm_active_dif_pair_test <- function(context, item, background) {
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

active_gllrm_no_dif_conditioning <- function(context, item, background) {
  source_backgrounds <- active_gllrm_item_dif_backgrounds(context, item)
  active_backgrounds <- active_gllrm_model_dif_backgrounds(context)
  xtra_backgrounds <- setdiff(active_backgrounds, c(background, source_backgrounds))
  biased_items <- active_gllrm_background_dif_items(context, background)
  biased_items <- setdiff(biased_items, item)

  values <- cbind(context$score)
  dims <- context$max_total_score - 1L
  for (background_index in c(source_backgrounds, xtra_backgrounds)) {
    values <- cbind(values, context$background_matrix[, background_index])
    dims <- c(dims, context$background_raw_max[[background_index]])
  }
  for (item_index in biased_items) {
    values <- cbind(values, context$item_matrix[, item_index] + 1L)
    dims <- c(dims, context$item_raw_max[[item_index]])
  }
  list(values = values, dims = dims)
}

active_gllrm_dif_lookup <- function(context) {
  out <- list()
  for (spec in context$dif_specs) {
    out[[active_gllrm_dif_key(spec$item, spec$background)]] <- TRUE
  }
  out
}

active_gllrm_dif_key <- function(item, background) {
  paste(item, background, sep = ":")
}

active_gllrm_item_dif_backgrounds <- function(context, item) {
  sort(unique(vapply(
    context$dif_specs[vapply(context$dif_specs, function(spec) identical(spec$item, item), logical(1L))],
    function(spec) spec$background,
    integer(1L)
  )))
}

active_gllrm_model_dif_backgrounds <- function(context) {
  sort(unique(vapply(context$dif_specs, function(spec) spec$background, integer(1L))))
}

active_gllrm_background_dif_items <- function(context, background) {
  sort(unique(vapply(
    context$dif_specs[vapply(context$dif_specs, function(spec) identical(spec$background, background), logical(1L))],
    function(spec) spec$item,
    integer(1L)
  )))
}
