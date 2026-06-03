#' Subset a DIGRAM bundle to one exogenous variable category
#'
#' @param bundle Source-shaped bundle.
#' @param background_index One-based background variable index.
#' @param background_value One-based background category value.
#' @return Bundle with rows outside the category marked invalid.
#' @keywords internal
subset_bundle_to_background_value <- function(bundle, background_index, background_value) {
  backgrounds <- bundle$model$backgrounds
  background_name <- backgrounds$name[[background_index]]
  group_bundle <- bundle
  keep <- bundle$data$status == 1L & bundle$data[[background_name]] == background_value
  group_bundle$data$status[!keep] <- 0L
  group_bundle$manifest$nvalid <- sum(keep)
  group_bundle
}

#' DIGRAM source name for an exogenous category
#'
#' @param value One-based category value.
#' @return Source-style category name.
#' @keywords internal
source_exo_value_name <- function(value) {
  names <- c(
    "one", "two", "three", "four", "five", "six", "seven", "eight",
    "nine", "ten"
  )
  if (value >= 1L && value <= length(names)) {
    names[[value]]
  } else {
    as.character(value)
  }
}

#' Derive DIGRAM global invariance report values
#'
#' Computes the values printed by DIGRAM's
#' `check-global-invariance.txt` report. The source path is
#' `DGRirtD.pas`'s `TestGlobalDIF` branch: for each selected exogenous
#' variable category, DIGRAM refits the model for that category to accumulate
#' the category likelihood, then prints `skbias15.pas` item-mean residual rows
#' using the current full-model item parameters filtered to the category.
#'
#' The item mean residual rows reuse the same source-shaped helper as global
#' homogeneity. The exact printed residual and marker cells are left as `NA`:
#' the Pascal source identifies the `skbias15` residual path, but the historical
#' runtime residual boundary has not been reproduced source-faithfully without
#' empirical report-specific correction factors.
#'
#' @param project A parsed DIGRAM project from [read_digram_project()].
#' @param max_step Maximum number of Rasch IPF iterations.
#' @param max_delta Convergence threshold for Rasch IPF.
#' @return A `gRm_global_invariance_values` object.
#' @examples
#' \dontrun{
#' project <- read_digram_project("path/to/DIGRAM")
#' values <- global_invariance_values(project)
#' values$tests
#' }
#' @keywords internal
global_invariance_values <- function(project, max_step = 5000L, max_delta = 0.0001) {
  bundle <- build_item_parameters_bundle(project)
  backgrounds <- bundle$model$backgrounds
  full_fit <- fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)
  full_loglike <- base_rasch_loglike(bundle, full_fit$item_gamma)
  n_parameters <- calculate_source_n_parameters(full_fit$counts$item_counts)

  sections <- list()
  tests <- list()
  section_index <- 0L

  for (background_index in seq_len(nrow(backgrounds))) {
    background <- backgrounds[background_index, , drop = FALSE]
    subgroup_loglike_sum <- 0
    max_delta_seen <- 0
    n_invalid_groups <- 0L

    for (background_value in seq_len(background$raw_max[[1L]])) {
      group_bundle <- subset_bundle_to_background_value(bundle, background_index, background_value)
      group_fit <- fit_rasch_base(group_bundle, max_step = max_step, max_delta = max_delta)
      if (group_fit$counts$n_valid > 0L) {
        group_loglike <- base_rasch_loglike(group_bundle, group_fit$item_gamma)
        subgroup_loglike_sum <- subgroup_loglike_sum + group_loglike
        if (!group_fit$converged) {
          max_delta_seen <- max(max_delta_seen, group_fit$delta)
        }
        item_rows <- global_homogeneity_item_mean_rows(
          group_bundle,
          group_fit,
          background_value,
          full_fit$item_gamma
        )
        # Source/parity note: global invariance calls the same item-mean
        # residual helper as global homogeneity because DGRirtD.pas dispatches
        # both report families through skbias15.pas::Calculate_residuals_and_item_fits.
        # The helper deliberately leaves residual and marker as NA because the
        # exact historical runtime residual boundary is not source-backed.
      } else {
        group_loglike <- 0
        n_invalid_groups <- n_invalid_groups + 1L
        item_rows <- data.frame()
      }

      section_index <- section_index + 1L
      sections[[section_index]] <- list(
        background_index = background_index,
        background_label = background$label_code[[1L]],
        background_name = background$name[[1L]],
        background_display_name = substr(background$name[[1L]], 1L, 8L),
        background_value = background_value,
        background_value_name = source_exo_value_name(background_value),
        n = group_fit$counts$n_valid,
        log_likelihood = group_loglike,
        converged = group_fit$converged,
        delta = group_fit$delta,
        items = item_rows
      )
    }

    n_groups <- background$raw_max[[1L]]
    # Source trace: DGRirtD.pas computes df as
    # (ngrp - 1 - nInvalidGroups) * Nparameters, then subtracts existing DIF
    # parameter contributions. The example model has no active DIF terms here.
    df <- (n_groups - 1L - n_invalid_groups) * n_parameters
    clr <- 2 * abs(full_loglike - subgroup_loglike_sum)
    p_value <- source_pfchi(df, clr)

    tests[[background_index]] <- data.frame(
      background_index = background_index,
      background_label = background$label_code[[1L]],
      background_name = background$name[[1L]],
      background_display_name = substr(background$name[[1L]], 1L, 8L),
      n_groups = n_groups,
      n_invalid_groups = n_invalid_groups,
      chi_square = clr,
      degrees_of_freedom = df,
      p_value = p_value,
      max_delta = max_delta_seen,
      stringsAsFactors = FALSE
    )
  }

  tests <- do.call(rbind, tests)

  # Source trace: the DGRirtD summary loop stores clrtable[0,3] instead of
  # clrtable[i,3] for each exogenous variable. With no simultaneous global
  # homogeneity test in the supplied example report, this behaves as zero p-values
  # for the BH summary, so the displayed critical levels remain the nominal
  # alpha levels.
  summary_bh_p_values <- rep(0, nrow(tests))

  result <- structure(
    list(
      bundle = bundle,
      fit = full_fit,
      sections = sections,
      tests = tests,
      bh = list(
        fdr_05 = source_bh_critical(summary_bh_p_values, 0.05),
        fdr_01 = source_bh_critical(summary_bh_p_values, 0.01),
        fdr_001 = source_bh_critical(summary_bh_p_values, 0.001)
      )
    ),
    class = "gRm_global_invariance_values"
  )
  result
}
