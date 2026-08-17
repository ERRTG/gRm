#' Derive DIGRAM global homogeneity numeric values
#'
#' Computes non-GUI DIGRAM global homogeneity values from raw DIGRAM
#' input using native R code.
#'
#' The item mean rows are computed from the `skbias15.pas`-shaped global
#' homogeneity helper. Row counts, observed means, expected means, and summary
#' CLR/df/p values are protected by tests.
#' The printed item-mean residual and marker cells are deliberately `NA`
#' because the available source does not fully explain the historical runtime's
#' hidden residual variance materialization.
#'
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
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
#' @param bundle Source-shaped analysis bundle.
#' @param base_fit Internal `base_fit` value used by this helper.
#' @param fit Fitted gRm model.
#' @keywords internal
#' @noRd
global_homogeneity_values <- function(project,
                                      score_cuts,
                                      max_step = 5000L,
                                      max_delta = 0.0001,
                                      bundle = NULL,
                                      base_fit = NULL,
                                      fit = NULL) {
  if (inherits(project, "gRm_fit") && inherits(project$values, "gRm_gllrm_values")) {
    return(gllrm_global_homogeneity_values(
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

#' Internal gllrm global homogeneity values helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param fit Fitted gRm model.
#' @param score_cuts Resolved total-score cut values.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @return The internal `gllrm_global_homogeneity_values()` computation result.
#' @keywords internal
#' @noRd
gllrm_global_homogeneity_values <- function(fit,
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
    group_gllrm_fit <- fit_gllrm(
      fit$spec,
      max_step = max_step,
      max_delta = max_delta,
      bundle = group_bundle
    )
    group_values_object <- gllrm_values(group_gllrm_fit, fit$spec)
    group_loglike <- group_gllrm_fit$state$log_likelihood
    subgroup_loglike_sum <- subgroup_loglike_sum + group_loglike
    group_values[[group_index]] <- list(
      bundle = group_bundle,
      fit = list(
        context = group_gllrm_fit$context,
        state = group_gllrm_fit$state,
        values = group_values_object
      ),
      log_likelihood = group_loglike
    )
    item_rows[[group_index]] <- gllrm_global_homogeneity_item_mean_rows(
      group_gllrm_fit$context,
      full_state,
      groups$group[[group_index]],
      probability_cache = new_gllrm_probability_cache(group_gllrm_fit$context, full_state)
    )
    groups$n[[group_index]] <- group_gllrm_fit$context$counts$n_valid
    groups$log_likelihood[[group_index]] <- group_loglike
    groups$converged[[group_index]] <- group_gllrm_fit$state$converged
    groups$delta[[group_index]] <- group_gllrm_fit$state$report_delta
  }

  full_loglike <- fit$values$log_likelihood
  n_parameters <- fit$values$n_parameters
  clr <- 2 * abs(full_loglike - subgroup_loglike_sum)
  df <- (nrow(groups) - 1L) * n_parameters
  p_value <- source_pfchi(df, clr)
  full_probability_cache <- new_gllrm_probability_cache(fit$fit$context, fit$fit)

  result <- list(
    item_parameters = fit$values,
    bundle = bundle,
    fit = fit$fit,
    score_groups = groups,
    group_values = group_values,
    items = do.call(rbind, item_rows),
    uniform_ld = gllrm_global_homogeneity_uniform_ld(
      fit,
      groups,
      probability_cache = full_probability_cache
    ),
    uniform_dif = gllrm_global_homogeneity_uniform_dif(
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

#' Internal gllrm global homogeneity uniform ld helper
#'
#' Supports the global homogeneity values implementation while preserving its internal contract.
#' Source trace: `source/PAS_scd/DGRirtD.pas::GlobalHomogeneity`.
#' @param fit Fitted gRm model.
#' @param groups Internal `groups` value used by this helper.
#' @param probability_cache Internal `probability_cache` value used by this helper.
#' @return The internal `gllrm_global_homogeneity_uniform_ld()` computation result.
#' @keywords internal
#' @noRd
gllrm_global_homogeneity_uniform_ld <- function(fit, groups, probability_cache = NULL) {
  context <- fit$fit$context
  state <- fit$fit
  if (length(context$ld_specs) == 0L) {
    return(data.frame())
  }

  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  probability_cache <- probability_cache %||%
    new_gllrm_probability_cache(context, state, components = components)
  tables_by_ld <- gllrm_uniform_ld_scoregroup_tables_all(context, groups, probability_cache)
  rows <- vector("list", length(context$ld_specs))
  for (ld_index in seq_along(context$ld_specs)) {
    spec <- context$ld_specs[[ld_index]]
    rows[[ld_index]] <- gllrm_uniform_ld_summary_row(context, groups, spec, tables_by_ld[[ld_index]])
  }
  do.call(rbind, rows)
}
