#' Internal item fits values helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param project Encoded gRm project.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @param include_extended Internal `include_extended` value used by this helper.
#' @return The internal `item_fits_values()` computation result.
#' @keywords internal
#' @noRd
item_fits_values <- function(project, max_step = 5000L, max_delta = 0.0001, include_extended = TRUE) {
  if (inherits(project, "gRm_fit") && inherits(project$values, "gRm_gllrm_values")) {
    return(gllrm_item_fits_values(project, include_extended = include_extended))
  }
  bundle <- build_item_parameters_bundle(project)
  fit <- fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)
  base_item_fits_values(bundle, fit, include_extended = include_extended)
}
#' Internal assemble item fits values helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @param fit_like Internal `fit_like` value used by this helper.
#' @param item_fit Internal `item_fit` value used by this helper.
#' @param gamma_fit Internal `gamma_fit` value used by this helper.
#' @param extended Internal `extended` value used by this helper.
#' @param component_gamma Internal `component_gamma` value used by this helper.
#' @param gllrm_copy_exists Internal `gllrm_copy_exists` value used by this helper.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
assemble_item_fits_values <- function(bundle,
                                      fit_like,
                                      item_fit,
                                      gamma_fit,
                                      extended,
                                      component_gamma = NULL,
                                      gllrm_copy_exists = FALSE) {
  rows <- data.frame(
    item_label = bundle$model$items$label_code,
    item_name = bundle$model$items$name,
    outfit = item_fit$outfit,
    outfit_sd = item_fit$outfit_sd,
    p_outfit = item_fit$p_outfit,
    outfit_fdr = item_fdr_risk(item_fit$p_outfit),
    infit = item_fit$infit,
    infit_sd = item_fit$infit_sd,
    p_infit = item_fit$p_infit,
    infit_fdr = item_fdr_risk(item_fit$p_infit),
    observed_gamma = gamma_fit$observed_gamma,
    expected_gamma = gamma_fit$expected_gamma,
    gamma_sd = gamma_fit$gamma_sd,
    p_gamma = gamma_fit$p_gamma,
    gamma_fdr = item_fdr_risk(gamma_fit$p_gamma),
    stringsAsFactors = FALSE
  )
  rows$direction <- item_fit_direction(rows)

  side_file <- data.frame(
    item = rows$item_label,
    outfit = rows$outfit,
    p_outfit = rows$p_outfit,
    infit = rows$infit,
    p_infit = rows$p_infit,
    ObsGamma = rows$observed_gamma,
    ExpGamma = rows$expected_gamma,
    p_gamma = rows$p_gamma,
    stringsAsFactors = FALSE
  )

  component_rows <- if (is.null(component_gamma)) {
    NULL
  } else {
    component_gamma$rows
  }
  all_p <- c(rows$p_infit, rows$p_outfit, rows$p_gamma)
  if (is.data.frame(component_rows) && nrow(component_rows)) {
    all_p <- c(all_p, component_rows$p_gamma)
  }

  result <- list(
    items = rows,
    side_file = side_file,
    bh_limits = c(
      fdr_5 = source_bh_critical(all_p, 0.05),
      fdr_1 = source_bh_critical(all_p, 0.01)
    ),
    fit = fit_like,
    extended = extended,
    incomplete_records_used = FALSE,
    incomplete_records_status = "not_source_backed"
  )
  if (!is.null(component_rows)) {
    result$component_gamma <- component_rows
  }
  if (isTRUE(gllrm_copy_exists)) {
    result$gllrm_copy_exists <- TRUE
  }
  class(result) <- c("gRm_item_fits_values", class(result))
  result
}

#' Internal base item fits values helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param bundle Source-shaped analysis bundle.
#' @param fit Fitted gRm model.
#' @param include_extended Internal `include_extended` value used by this helper.
#' @return The internal `base_item_fits_values()` computation result.
#' @keywords internal
#' @noRd
base_item_fits_values <- function(bundle, fit, include_extended = TRUE) {
  conditional <- item_conditional_moments(bundle, fit$item_gamma, include_probabilities = TRUE)
  # Source-faithfulness guard: Pascal item fits consume incomplete records only
  # after the include-incomplete analysis path has populated NincompleteRecs
  # and related runtime arrays via Execute_incomplete_GLLRM_estimates. Do not
  # wire collect_source_incomplete_records() into public item_fit() until that
  # gating and completion path is implemented end to end.
  incomplete <- empty_source_incomplete_records(bundle)

  item_fit <- calculate_conditional_item_fit_values(bundle, fit, conditional = conditional, incomplete = incomplete)
  gamma_fit <- calculate_item_restscore_gamma_values(bundle, fit, conditional = conditional, incomplete = incomplete)
  extended <- if (include_extended) {
    calculate_extended_item_fit_values(bundle, fit, conditional = conditional)
  } else {
    NULL
  }

  assemble_item_fits_values(
    bundle = bundle,
    fit_like = fit,
    item_fit = item_fit,
    gamma_fit = gamma_fit,
    extended = extended
  )
}

#' Internal gllrm item fits values helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param fit Fitted gRm model.
#' @param include_extended Internal `include_extended` value used by this helper.
#' @return The internal `gllrm_item_fits_values()` computation result.
#' @keywords internal
#' @noRd
gllrm_item_fits_values <- function(fit, include_extended = TRUE) {
  context <- fit$fit$context
  state <- fit$fit
  state$context <- NULL
  bundle <- context$bundle
  fit_like <- list(
    model = "gllrm",
    item_gamma = state$item_gamma,
    counts = context$counts,
    context = context,
    state = state
  )
  conditional <- gllrm_item_conditional_moments(
    context,
    state,
    include_probabilities = TRUE,
    probability_cache = new_gllrm_probability_cache(context, state)
  )

  # Source-faithfulness guard: Pascal item fits consume incomplete records only
  # after the include-incomplete analysis path has populated NincompleteRecs
  # and related runtime arrays via Execute_incomplete_GLLRM_estimates. Do not
  # wire collect_source_incomplete_records() into public item_fit() until that
  # gating and completion path is implemented end to end.
  incomplete <- empty_source_incomplete_records(bundle)
  item_fit <- calculate_conditional_item_fit_values(bundle, fit_like, conditional = conditional, incomplete = incomplete)
  gamma_fit <- calculate_item_restscore_gamma_values(bundle, fit_like, conditional = conditional, incomplete = incomplete)
  component_gamma <- gllrm_component_restscore_values(context, state)
  extended <- if (include_extended) {
    out <- calculate_extended_item_fit_values(bundle, fit_like, conditional = conditional)
    out$component_restscore_tables <- component_gamma$table_rows
    out
  } else {
    NULL
  }

  assemble_item_fits_values(
    bundle = bundle,
    fit_like = fit_like,
    item_fit = item_fit,
    gamma_fit = gamma_fit,
    extended = extended,
    component_gamma = component_gamma,
    gllrm_copy_exists = TRUE
  )
}

#' Internal gllrm component restscore values helper
#'
#' Supports the item fits values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`.
#' @param context Prepared GLLRM computation context.
#' @param state Current fitted or iterative parameter state.
#' @return The internal `gllrm_component_restscore_values()` computation result.
#' @keywords internal
#' @noRd
gllrm_component_restscore_values <- function(context, state) {
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  dependent_components <- components[lengths(components) > 1L]
  if (length(dependent_components) == 0L) {
    empty_rows <- data.frame(
      component = character(),
      representative_item = integer(),
      representative_item_label = character(),
      observed_gamma = numeric(),
      expected_gamma = numeric(),
      gamma_sd = numeric(),
      p_gamma = numeric(),
      stringsAsFactors = FALSE
    )
    empty_tables <- data.frame(
      component = character(),
      representative_item_label = character(),
      restscore_table = character(),
      restscore = integer(),
      component_score = integer(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
    return(list(rows = empty_rows, table_rows = empty_tables))
  }

  rows <- list()
  table_rows <- list()
  for (component_items in dependent_components) {
    tables <- gllrm_component_restscore_tables(context, state, component_items)
    observed_gamma <- goodman_kruskal_gamma(tables$observed)
    fitted <- fitted_gamma_stats(tables$expected)
    expected_gamma <- fitted$gamma
    gamma_sd <- sqrt(fitted$variance)
    p_gamma <- if (gamma_sd > 0) {
      two_sided_source_normal_p((observed_gamma - expected_gamma) / gamma_sd)
    } else {
      1
    }
    component <- paste(context$items$label_code[component_items], collapse = "")
    representative_item <- component_items[[length(component_items)]]
    representative_item_label <- context$items$label_code[[representative_item]]
    rows[[length(rows) + 1L]] <- data.frame(
      component = component,
      representative_item = representative_item,
      representative_item_label = representative_item_label,
      observed_gamma = observed_gamma,
      expected_gamma = expected_gamma,
      gamma_sd = gamma_sd,
      p_gamma = p_gamma,
      stringsAsFactors = FALSE
    )
    table_rows[[length(table_rows) + 1L]] <- component_restscore_table_rows(
      component = component,
      representative_item_label = representative_item_label,
      observed_table = tables$observed,
      expected_table = tables$expected
    )
  }

  list(
    rows = type.convert(do.call(rbind, rows), as.is = TRUE),
    table_rows = type.convert(do.call(rbind, table_rows), as.is = TRUE)
  )
}
