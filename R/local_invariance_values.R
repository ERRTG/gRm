#' Subset a DIGRAM bundle to a contiguous exogenous category range
#'
#' @param bundle Source-shaped bundle.
#' @param background_index One-based background variable index.
#' @param from_value First one-based category value.
#' @param to_value Last one-based category value.
#' @return Bundle with rows outside the category range marked invalid.
#' @keywords internal
subset_bundle_to_background_range <- function(bundle, background_index, from_value, to_value) {
  backgrounds <- bundle$model$backgrounds
  background_name <- backgrounds$name[[background_index]]
  group_bundle <- bundle
  keep <- bundle$data$status == 1L &
    bundle$data[[background_name]] >= from_value &
    bundle$data[[background_name]] <= to_value
  group_bundle$data$status[!keep] <- 0L
  group_bundle$manifest$nvalid <- sum(keep)
  group_bundle
}

#' DIGRAM source label for a collapsed ordinal group
#'
#' @param from_value First one-based category value.
#' @param to_value Last one-based category value.
#' @return Source-style group label.
#' @keywords internal
source_group_label <- function(from_value, to_value) {
  if (from_value == to_value) {
    as.character(from_value)
  } else {
    paste0(from_value, "+", to_value)
  }
}

#' Derive DIGRAM local-invariance MCA report values
#'
#' Computes the restricted MCA analysis printed by DIGRAM's
#' `test-of-local-invariance.txt` report. The source path is
#' `DGRirtD.pas`'s `LocalDIF` branch calling
#' `skbias14.pas::MCAanalysis_of_LL_ratingscale`: DIGRAM compares adjacent
#' exogenous categories, refits Rasch models for each group and their union,
#' computes `CLR = |2*LL(union) - 2*LL(left) - 2*LL(right)|`, converts the
#' statistic with `SkStat.PFCHI`, and applies the source Benjamini-Hochberg
#' cutoff at each step.
#'
#' @param project A parsed DIGRAM project from [read_digram_project()] or a
#'   source-shaped project list.
#' @param max_step Maximum number of Rasch IPF iterations.
#' @param max_delta Convergence threshold for Rasch IPF.
#' @return A `gRm_local_invariance_values` object.
#' @examples
#' \dontrun{
#' project <- read_digram_project("path/to/DIGRAM",
#'   items = c(paste0("A", 1:5), paste0("C", 1:5), paste0("E", 1:5),
#'             paste0("N", 1:5), paste0("O", 1:5)),
#'   exogenous = c("gender", "education"))
#' values <- local_invariance_values(project)
#' }
#' @keywords internal
local_invariance_values <- function(project, max_step = 5000L, max_delta = 0.0001) {
  bundle <- build_item_parameters_bundle(project)
  backgrounds <- bundle$model$backgrounds
  n_parameters <- calculate_source_n_parameters(
    fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)$counts$item_counts
  )

  analyses <- list()
  analysis_index <- 0L

  for (background_index in seq_len(nrow(backgrounds))) {
    background <- backgrounds[background_index, , drop = FALSE]
    n_categories <- background$raw_max[[1L]]
    if (n_categories <= 2L) {
      next
    }

    fit_cache <- new.env(parent = emptyenv())
    fit_range <- function(from_value, to_value) {
      key <- paste(from_value, to_value, sep = ":")
      if (!exists(key, envir = fit_cache, inherits = FALSE)) {
        group_bundle <- subset_bundle_to_background_range(
          bundle, background_index, from_value, to_value
        )
        # Source trace: skbias14.pas::EstimateModel applies the exofilter and
        # then calls Estimate_GLLRM. This R implementation refits the same
        # Rasch base model to the rows inside the current category interval.
        fit <- fit_rasch_base(group_bundle, max_step = max_step, max_delta = max_delta)
        loglike <- base_rasch_loglike(group_bundle, fit$item_gamma)
        assign(key, list(
          from = from_value,
          to = to_value,
          label = source_group_label(from_value, to_value),
          n = fit$counts$n_valid,
          log_likelihood = loglike,
          converged = fit$converged,
          delta = fit$delta
        ), envir = fit_cache)
      }
      get(key, envir = fit_cache, inherits = FALSE)
    }

    category_counts <- vapply(seq_len(n_categories), function(value) {
      fit_range(value, value)$n
    }, integer(1))

    groups <- lapply(seq_len(n_categories), function(value) {
      list(from = value, to = value, active = TRUE)
    })
    initial_pairs <- NULL
    steps <- list()
    step_index <- 0L

    repeat {
      active_indices <- which(vapply(groups, `[[`, logical(1), "active"))
      if (length(active_indices) < 2L) {
        break
      }

      pair_rows <- vector("list", length(active_indices) - 1L)
      for (pair_index in seq_len(length(active_indices) - 1L)) {
        left_index <- active_indices[[pair_index]]
        right_index <- active_indices[[pair_index + 1L]]
        left_group <- groups[[left_index]]
        right_group <- groups[[right_index]]
        left_fit <- fit_range(left_group$from, left_group$to)
        right_fit <- fit_range(right_group$from, right_group$to)
        union_fit <- fit_range(left_group$from, right_group$to)
        # skbias14 stores min2loglike as 2*RaschLoglike. The CLR below is
        # therefore the likelihood-ratio contrast between the joined group and
        # the two separate adjacent groups.
        clr <- abs(2 * union_fit$log_likelihood -
          2 * left_fit$log_likelihood -
          2 * right_fit$log_likelihood)
        pair_rows[[pair_index]] <- data.frame(
          left_index = left_index,
          right_index = right_index,
          left_label = source_group_label(left_group$from, left_group$to),
          right_label = source_group_label(right_group$from, right_group$to),
          left_from = left_group$from,
          left_to = left_group$to,
          right_from = right_group$from,
          right_to = right_group$to,
          chi_square = clr,
          degrees_of_freedom = n_parameters,
          p_value = source_pfchi(n_parameters, clr),
          stringsAsFactors = FALSE
        )
      }

      pairs <- do.call(rbind, pair_rows)
      fdr05 <- source_bh_critical(pairs$p_value, 0.05)
      max_row_index <- which.max(pairs$p_value)
      max_row <- pairs[max_row_index, , drop = FALSE]
      step_index <- step_index + 1L
      steps[[step_index]] <- cbind(
        data.frame(
          step = step_index,
          n_pvalues = nrow(pairs),
          fdr05 = fdr05,
          collapsed = max_row$p_value >= fdr05,
          stringsAsFactors = FALSE
        ),
        max_row
      )

      if (is.null(initial_pairs)) {
        initial_pairs <- pairs
        initial_fdr05 <- fdr05
      }

      # Source trace: skbias14 prints the current max-p row first. Only after
      # that does it collapse the groups when maxp >= FDR05.
      if (max_row$p_value >= fdr05) {
        groups[[max_row$left_index]]$to <- groups[[max_row$right_index]]$to
        groups[[max_row$right_index]]$active <- FALSE
      } else {
        break
      }
    }

    final_groups <- do.call(rbind, lapply(which(vapply(groups, `[[`, logical(1), "active")), function(index) {
      group <- groups[[index]]
      data.frame(
        from = group$from,
        to = group$to,
        label = source_group_label(group$from, group$to),
        stringsAsFactors = FALSE
      )
    }))

    final_pairs <- NULL
    if (nrow(final_groups) >= 2L) {
      final_pairs <- do.call(rbind, lapply(seq_len(nrow(final_groups) - 1L), function(pair_index) {
        left_group <- final_groups[pair_index, ]
        right_group <- final_groups[pair_index + 1L, ]
        left_fit <- fit_range(left_group$from, left_group$to)
        right_fit <- fit_range(right_group$from, right_group$to)
        union_fit <- fit_range(left_group$from, right_group$to)
        clr <- abs(2 * union_fit$log_likelihood -
          2 * left_fit$log_likelihood -
          2 * right_fit$log_likelihood)
        data.frame(
          left_label = left_group$label,
          right_label = right_group$label,
          chi_square = clr,
          degrees_of_freedom = n_parameters,
          p_value = source_pfchi(n_parameters, clr),
          stringsAsFactors = FALSE
        )
      }))
    }

    fits <- as.list(fit_cache)
    fits <- do.call(rbind, lapply(fits, as.data.frame, stringsAsFactors = FALSE))

    analysis_index <- analysis_index + 1L
    analyses[[analysis_index]] <- list(
      background_index = background_index,
      background_name = background$name[[1L]],
      background_display_name = substr(background$name[[1L]], 1L, 8L),
      background_label = background$label_code[[1L]],
      n_categories = n_categories,
      category_counts = category_counts,
      n_parameters = n_parameters,
      initial_pairs = initial_pairs,
      initial_fdr05 = initial_fdr05,
      steps = do.call(rbind, steps),
      final_groups = final_groups,
      final_pairs = final_pairs,
      fits = fits
    )
  }

  structure(
    list(bundle = bundle, analyses = analyses),
    class = "gRm_local_invariance_values"
  )
}
