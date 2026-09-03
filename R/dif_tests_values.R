#' Source observed IX-margin degrees of freedom
#'
#' Port of the `IXinfo^.df` value set by source `COUNT_MARGINS` for
#' item-by-exogeneous margins: nonzero observed item-score levels minus one,
#' times nonzero observed exogeneous levels minus one.
#'
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param observed_ix Integer observed item-score by exogeneous-value margin.
#' @return Degrees of freedom for the included IX margin.
#' @keywords internal
#' @noRd
source_ix_observed_df <- function(observed_ix) {
  nonzero_item_scores <- sum(rowSums(observed_ix) > 0)
  nonzero_exo_values <- sum(colSums(observed_ix) > 0)
  as.integer(max(0L, nonzero_item_scores - 1L) * max(0L, nonzero_exo_values - 1L))
}

#' Source Goodman-Kruskal gamma count totals
#'
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param tab Two-way integer table.
#' @return Gamma, PPQ, and PMQ totals.
#' @keywords internal
#' @noRd
dif_tests_source_gamma_counts <- function(tab) {
  source_rc_gamma_counts(tab)
}

#' Build cached raw-data context for DIF gamma values
#'
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param project DIGRAM project.
#' @return Raw item/background matrices and score vectors.
#' @keywords internal
#' @noRd
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

#' Calculate source-shaped DIF gamma diagnostics
#'
#' Source trace: `SKbias7.inexpensive_itembias1` calls
#' `SKbias3.XYZ_bias_ANALYSE`, stores `RESULTS(.1,5.)` as `Part_g` and
#' `RESULTS(.1,6.)` as `Part_p`, then `Item_screening` copies those matrices to
#' `IJXgamma` and `IJXgamma_pvalues`. `DGRirtD.MissingDIF` reports
#' `IJXgamma(.i,nitems+j.)`, and the GLLRM parameter reports label the
#' same value "Item screening: Gamma".
#'
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param project DIGRAM project.
#' @param context Context from `build_dif_gamma_context()`.
#' @param target_item One-based item index.
#' @param background_index One-based background index.
#' @return Source partial gamma, p-value, and accumulated PPQ/PMQ/S values.
#' @keywords internal
#' @noRd
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

#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @return Source partial gamma.
#' @param project Encoded gRm project.
#' @param context Prepared GLLRM computation context.
#' @param target_item Internal `target_item` value used by this helper.
#' @param background_index Internal `background_index` value used by this helper.
#' @keywords internal
#' @noRd
dif_tests_partial_gamma <- function(project, context, target_item, background_index) {
  dif_tests_partial_gamma_stats(project, context, target_item, background_index)$gamma
}

#' Convolve two score polynomials with source score bounds
#'
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param left Numeric score vector indexed by score plus one.
#' @param right Numeric score vector indexed by score plus one.
#' @param max_total_score Maximum retained score.
#' @return Numeric score vector.
#' @keywords internal
#' @noRd
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

#' Calculate the base Rasch negative log likelihood
#'
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param bundle Source-shaped bundle.
#' @param item_gamma Item gamma matrix.
#' @return Negative conditional log likelihood.
#' @keywords internal
#' @noRd
base_rasch_loglike <- function(bundle, item_gamma) {
  data <- bundle$data
  items <- bundle$model$items
  valid_rows <- which(data$status == 1L)
  # Source trace: source/digram_source_20260817/skunits/skbias12.pas::
  # Inexpensive_Gamma_Calculation builds the same denominator used by the base
  # Rasch estimator. Candidate diagnostics now share this canonical routine.
  full_gamma <- build_source_score_gamma(
    bundle,
    item_gamma,
    use_items = rep(TRUE, nrow(items))
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

#' Compute DIGRAM CHECK D DIF-test values
#'
#' Computes the item-by-exogenous-variable no-DIF report values used in
#' DIGRAM's standard `check-missing-DIF.txt` report with a native R port of the
#' one-included-DIF source likelihood-ratio path. The algorithm follows
#' `docs/source-traces/BIRT_DIF_TESTS_SOURCE_TRACE.md`, especially `DGRirtD`'s CHECK D branch,
#' `Estimate_GLLRM`, `PFCHI`, and the source Benjamini-Hochberg rule.
#'
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param project A source-shaped DIGRAM project list, such as the `project`
#'   component returned by [gRm()] or [read_digram_project()].
#' @param max_step Maximum source Rasch/GLLRM estimation iterations.
#' @param max_delta Convergence threshold for source Rasch/GLLRM estimation.
#' @param jobs Upper bound for parallel candidate DIF fits.
#' @param progress Logical; if `TRUE`, emit temporary per-candidate progress
#'   messages for long-running validation/debugging runs.
#' @return A `gRm_dif_tests_values` object with test rows and BH threshold.
#' @keywords internal
#' @noRd
dif_tests_values <- function(project, max_step = 5000L, max_delta = 0.0001, jobs = min(8L, parallel::detectCores(logical = TRUE), 128L), progress = identical(Sys.getenv("RDIGRAM_DIF_PROGRESS"), "1")) {
  if (inherits(project, "gRm_fit") && inherits(project$values, "gRm_gllrm_values")) {
    return(gllrm_dif_tests_values(project, max_step = max_step, max_delta = max_delta, jobs = jobs))
  }
  candidate_object <- if (inherits(project, "gRm_fit")) project else NULL
  if (inherits(project, "gRm_fit")) {
    project <- project$project %||% project$analysis$project
  }
  if (!is.list(project)) {
    stop("`project` must be a DIGRAM project list.", call. = FALSE)
  }
  if (is.null(project$backgrounds) || nrow(project$backgrounds) == 0L) {
    stop("DIF tests require at least one exogeneous variable.", call. = FALSE)
  }

  candidate_object <- candidate_object %||% gllrm_candidate_base_model(project)
  bundle <- candidate_object$bundle %||% build_item_parameters_bundle(project)
  gamma_context <- build_dif_gamma_context(project)
  base_fit <- if (inherits(candidate_object, "gRm_fit")) {
    candidate_object$fit
  } else {
    fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)
  }
  base_loglike <- base_rasch_loglike(bundle, base_fit$item_gamma)
  items <- project$items
  backgrounds <- project$backgrounds

  candidates <- expand.grid(
    item_index = seq_len(nrow(items)),
    background_index = seq_len(nrow(backgrounds))
  )
  candidates <- candidates[order(candidates$background_index, candidates$item_index), , drop = FALSE]

  fit_one <- function(candidate_row) {
    item_index <- candidates$item_index[[candidate_row]]
    background_index <- candidates$background_index[[candidate_row]]
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
    # Source trace: DGRirtD.MissingDIF reaches Find_new_IXparameters and
    # Adjust_IXparameters through the shared Estimate_GLLRM candidate engine.
    candidate_fit <- fit_gllrm_candidate_dif(
      candidate_object,
      item = item_index,
      background = background_index,
      max_step = max_step,
      max_delta = max_delta
    )
    candidate_loglike <- candidate_fit$log_likelihood
    dif_index <- gllrm_context_dif_index(
      candidate_fit$context,
      item_index,
      background_index
    )
    df <- source_ix_observed_df(candidate_fit$context$observed_dif[[dif_index]])
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
          candidate_fit$report_delta
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
      # The former formatted CLR/p-value repetition rule was not a Pascal stop
      # condition. Retain the compatibility column without using it to stop.
      output_stable = FALSE,
      delta = candidate_fit$report_delta,
      n_step = candidate_fit$n_step,
      stop_reason = candidate_fit$stop_reason,
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

#' Internal gllrm dif tests values helper
#'
#' Supports the dif tests values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param fit Fitted gRm model.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @param jobs Requested worker count.
#' @return The internal `gllrm_dif_tests_values()` computation result.
#' @keywords internal
#' @noRd
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
    # Source trace: DGRirtD.MissingDIF reaches Find_new_IXparameters and
    # Adjust_IXparameters through the shared Estimate_GLLRM candidate engine.
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

#' Internal gllrm context dif index helper
#'
#' Supports the dif tests values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param context Prepared GLLRM computation context.
#' @param item One-based item index.
#' @param background One-based exogenous-variable index.
#' @return The internal `gllrm_context_dif_index()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal gllrm no dif candidates helper
#'
#' Supports the dif tests values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_no_dif_candidates()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal gllrm dif tests helper
#'
#' Supports the dif tests values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_dif_tests()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal gllrm dif pair test helper
#'
#' Supports the dif tests values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param context Prepared GLLRM computation context.
#' @param item One-based item index.
#' @param background One-based exogenous-variable index.
#' @return The internal `gllrm_dif_pair_test()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal gllrm dif lookup helper
#'
#' Supports the dif tests values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_dif_lookup()` computation result.
#' @keywords internal
#' @noRd
gllrm_dif_lookup <- function(context) {
  out <- list()
  for (spec in context$dif_specs) {
    out[[gllrm_dif_key(spec$item, spec$background)]] <- TRUE
  }
  out
}

#' Internal gllrm dif key helper
#'
#' Supports the dif tests values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/scd/DGRirtD.pas::MissingDIF`.
#' @param item One-based item index.
#' @param background One-based exogenous-variable index.
#' @return The internal `gllrm_dif_key()` computation result.
#' @keywords internal
#' @noRd
gllrm_dif_key <- function(item, background) {
  paste(item, background, sep = ":")
}
