#' Source Goodman-Kruskal gamma count totals
#'
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param tab Two-way integer table.
#' @return Gamma, PPQ, and PMQ totals.
#' @keywords internal
#' @noRd
local_independence_source_gamma_counts <- function(tab) {
  source_rc_gamma_counts(tab)
}

#' Calculate source-shaped local-independence WPG gamma
#'
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param project DIGRAM project.
#' @return Raw item matrix and source score vectors.
#' @keywords internal
#' @noRd
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
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param context Context from `build_local_independence_gamma_context()`.
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @param items Item metadata.
#' @return PPQ and PMQ totals.
#' @keywords internal
#' @noRd
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
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param project DIGRAM project.
#' @param context Context from `build_local_independence_gamma_context()`.
#' @param item1 One-based first item index.
#' @param item2 One-based second item index.
#' @return Source weighted partial gamma.
#' @keywords internal
#' @noRd
local_independence_wpg_gamma <- function(project, context, item1, item2) {
  items <- project$items
  forward <- local_independence_directed_gamma_counts(context, item1, item2, items)
  backward <- local_independence_directed_gamma_counts(context, item2, item1, items)
  ppq <- forward$ppq + backward$ppq
  pmq <- forward$pmq + backward$pmq
  if (ppq > 0) pmq / ppq else 0
}

#' Internal adjust ld gamma source reference details helper
#'
#' Supports the local independence values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param observed_ld Internal `observed_ld` value used by this helper.
#' @param ld_gamma Internal `ld_gamma` value used by this helper.
#' @param i_ref Internal `i_ref` value used by this helper.
#' @param j_ref Internal `j_ref` value used by this helper.
#' @param preserve_current_ties Internal `preserve_current_ties` value used by this helper.
#' @return The internal `adjust_ld_gamma_source_reference_details()` computation result.
#' @keywords internal
#' @noRd
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

#' Compute DIGRAM local-independence candidate test values
#'
#' Computes the item-pair likelihood-ratio tests used in DIGRAM's
#' `check-local-independence.txt` and
#' `check-local-independence-extended.txt` reports. The implementation follows
#' `docs/source-traces/BIRT_LOCAL_INDEPENDENCE_SOURCE_TRACE.md`, especially the
#' `DGRirtD.pas` `MissingLD` branch.
#'
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param project A source-shaped DIGRAM project list, such as the `project`
#'   component returned by [gRm()] or [read_digram_project()].
#' @param max_step Maximum source Rasch/GLLRM estimation iterations.
#' @param max_delta Convergence threshold.
#' @param jobs Number of parallel candidate LD fits.
#' @return A `gRm_local_independence_values` object.
#' @keywords internal
#' @noRd
local_independence_values <- function(project, max_step = 5000L, max_delta = 0.0001, jobs = min(32L, parallel::detectCores(logical = TRUE), 128L)) {
  if (inherits(project, "gRm_fit") && inherits(project$values, "gRm_gllrm_values")) {
    return(gllrm_local_independence_values(project, max_step = max_step, max_delta = max_delta, jobs = jobs))
  }
  candidate_object <- if (inherits(project, "gRm_fit")) project else NULL
  if (inherits(project, "gRm_fit")) {
    project <- project$project %||% project$analysis$project
  }
  if (!is.list(project)) {
    stop("`project` must be a DIGRAM project list.", call. = FALSE)
  }
  if (is.null(project$items) || nrow(project$items) < 2L) {
    stop("Local-independence tests require at least two item variables.", call. = FALSE)
  }

  candidate_object <- candidate_object %||% gllrm_candidate_base_model(project)
  bundle <- candidate_object$bundle %||% build_item_parameters_bundle(project)
  gamma_context <- build_local_independence_gamma_context(project)
  base_fit <- if (inherits(candidate_object, "gRm_fit")) {
    candidate_object$fit
  } else {
    fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)
  }
  base_loglike <- base_rasch_loglike(bundle, base_fit$item_gamma)
  items <- project$items
  candidates <- expand.grid(
    item1_index = seq_len(nrow(items)),
    item2_index = seq_len(nrow(items))
  )
  candidates <- candidates[candidates$item1_index < candidates$item2_index, , drop = FALSE]
  candidates <- candidates[order(candidates$item1_index, candidates$item2_index), , drop = FALSE]

  fit_one <- function(candidate_row) {
    item1 <- candidates$item1_index[[candidate_row]]
    item2 <- candidates$item2_index[[candidate_row]]
    # Source trace: DGRirtD.MissingLD reaches Find_new_IJparameters and
    # Adjust_IJparameters through the shared Estimate_GLLRM candidate engine.
    candidate_fit <- fit_gllrm_candidate_ld(
      candidate_object,
      item1 = item1,
      item2 = item2,
      max_step = max_step,
      max_delta = max_delta
    )
    candidate_loglike <- candidate_fit$log_likelihood

    # Source trace: DGRirtD MissingLD reports
    # lr := 2*abs(Raschloglike-Raschloglike1).
    clr <- 2 * abs(base_loglike - candidate_loglike)
    # Source trace: IJ df follows the observed nonzero item margins stored with
    # the candidate counts, not merely the declared raw category dimensions.
    ld_index <- gllrm_context_ld_index(candidate_fit$context, item1, item2)
    df <- source_ij_observed_df(candidate_fit$context$observed_ld[[ld_index]])

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
      delta = candidate_fit$report_delta,
      n_step = candidate_fit$n_step,
      stop_reason = candidate_fit$stop_reason,
      stringsAsFactors = FALSE
    )
  }

  tests <- source_candidate_map(nrow(candidates), jobs, fit_one)
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
  result
}

#' Internal gllrm local independence values helper
#'
#' Supports the local independence values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param fit Fitted gRm model.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @param jobs Requested worker count.
#' @return The internal `gllrm_local_independence_values()` computation result.
#' @keywords internal
#' @noRd
gllrm_local_independence_values <- function(fit, max_step = 5000L, max_delta = 0.0001, jobs = min(32L, parallel::detectCores(logical = TRUE), 128L)) {
  context <- fit$fit$context
  items <- context$items
  candidates <- gllrm_li_candidates(context)

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
      wpg_gamma_source = character(),
      converged = logical(),
      delta = numeric(),
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
        bh_critical_p = 0,
        suggested_ld = character(),
        max_step = as.integer(max_step),
        max_delta = max_delta
      ),
      class = "gRm_local_independence_values"
    ))
  }

  base_loglike <- fit$fit$log_likelihood
  project <- fit$project
  gamma_context <- build_local_independence_gamma_context(project)
  fit_one <- function(candidate_row) {
    item1 <- candidates$item1_index[[candidate_row]]
    item2 <- candidates$item2_index[[candidate_row]]
    # Source trace: DGRirtD.MissingLD reaches Find_new_IJparameters and
    # Adjust_IJparameters through the shared Estimate_GLLRM candidate engine.
    candidate_fit <- fit_gllrm_with_added_ld(
      fit,
      item1 = item1,
      item2 = item2,
      max_step = max_step,
      max_delta = max_delta
    )
    attempted_fit <- candidate_fit
    reported_checkpoint_step <- NA_integer_
    report_value_source <- "attempted_fit"
    # Source trace: source/PAS_scd/DGRirtD.pas::MissingLD calls
    # Estimate_GLLRM once and reports that fit. Post-50 recurrence is handled
    # inside the common estimator; there is no universal step-51 replacement.
    candidate_loglike <- candidate_fit$log_likelihood
    clr <- 2 * abs(base_loglike - candidate_loglike)
    ld_index <- gllrm_context_ld_index(candidate_fit$context, item1, item2)
    df <- source_ij_observed_df(
      candidate_fit$context$observed_ld[[ld_index]]
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
      # Source trace: DGRirtD.MissingLD prints WPGgamma from the
      # item-screening WPG matrix, not from the candidate refit.
      wpg_gamma = local_independence_wpg_gamma(project, gamma_context, item1, item2),
      wpg_gamma_source = "item_screening",
      converged = isTRUE(candidate_fit$converged),
      delta = candidate_fit$report_delta %||% NA_real_,
      stop_reason = candidate_fit$stop_reason %||% NA_character_,
      attempted_n_step = attempted_fit$n_step %||% NA_integer_,
      attempted_delta = attempted_fit$report_delta %||% attempted_fit$delta %||% NA_real_,
      attempted_converged = attempted_fit$converged %||% NA,
      attempted_stop_reason = attempted_fit$stop_reason %||% NA_character_,
      reported_checkpoint_step = reported_checkpoint_step,
      report_value_source = report_value_source,
      stringsAsFactors = FALSE
    )
  }

  tests <- source_candidate_map(nrow(candidates), jobs, fit_one)
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

#' Internal gllrm context ld index helper
#'
#' Supports the local independence values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param context Prepared GLLRM computation context.
#' @param item1 Internal `item1` value used by this helper.
#' @param item2 Internal `item2` value used by this helper.
#' @return The internal `gllrm_context_ld_index()` computation result.
#' @keywords internal
#' @noRd
gllrm_context_ld_index <- function(context, item1, item2) {
  key <- gllrm_ld_key(item1, item2)
  hit <- which(vapply(context$ld_specs, function(spec) {
    identical(gllrm_ld_key(spec$item1, spec$item2), key)
  }, logical(1L)))
  if (length(hit) != 1L) {
    stop("Could not identify the candidate LD term in the fitted GLLRM context.", call. = FALSE)
  }
  hit[[1L]]
}

#' Internal gllrm li candidates helper
#'
#' Supports the local independence values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param context Prepared GLLRM computation context.
#' @param components Internal `components` value used by this helper.
#' @return The internal `gllrm_li_candidates()` computation result.
#' @keywords internal
#' @noRd
gllrm_li_candidates <- function(context, components = NULL) {
  rows <- expand.grid(
    item1_index = seq_len(context$n_items),
    item2_index = seq_len(context$n_items)
  )
  rows <- rows[rows$item1_index < rows$item2_index, , drop = FALSE]
  included <- gllrm_ld_lookup(context)
  rows <- rows[
    !vapply(
      seq_len(nrow(rows)),
      function(i) included[[gllrm_ld_key(rows$item1_index[[i]], rows$item2_index[[i]])]] %||% FALSE,
      logical(1L)
    ),
    ,
    drop = FALSE
  ]
  rows[order(rows$item1_index, rows$item2_index), , drop = FALSE]
}

#' Internal gllrm ld lookup helper
#'
#' Supports the local independence values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_ld_lookup()` computation result.
#' @keywords internal
#' @noRd
gllrm_ld_lookup <- function(context) {
  out <- list()
  for (spec in context$ld_specs) {
    out[[gllrm_ld_key(spec$item1, spec$item2)]] <- TRUE
  }
  out
}

#' Internal gllrm ld key helper
#'
#' Supports the local independence values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param item1 Internal `item1` value used by this helper.
#' @param item2 Internal `item2` value used by this helper.
#' @return The internal `gllrm_ld_key()` computation result.
#' @keywords internal
#' @noRd
gllrm_ld_key <- function(item1, item2) {
  paste(min(item1, item2), max(item1, item2), sep = ":")
}

#' Internal source ij observed df helper
#'
#' Supports the local independence values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::MissingLD`.
#' @param observed_ij Internal `observed_ij` value used by this helper.
#' @return The internal `source_ij_observed_df()` computation result.
#' @keywords internal
#' @noRd
source_ij_observed_df <- function(observed_ij) {
  nonzero_item1_scores <- sum(rowSums(observed_ij) > 0)
  nonzero_item2_scores <- sum(colSums(observed_ij) > 0)
  as.integer(max(0L, nonzero_item1_scores - 1L) * max(0L, nonzero_item2_scores - 1L))
}
