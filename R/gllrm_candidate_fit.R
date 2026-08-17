# GLLRM candidate fits for LD and DIF diagnostics.
#
# Source trace: CHECK LID and CHECK DIF candidate tests compare the current
# current GLLRM with a source-ordered model containing exactly one additional
# IJ or IX term. These helpers run the same GLLRM estimation core used by
# public fit(), but skip public item-parameter output construction that the
# candidate likelihood-ratio tests do not use.

#' Internal gllrm candidate base model helper
#'
#' Supports the gllrm candidate fit implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param object Object dispatched to this helper.
#' @return The internal `gllrm_candidate_base_model()` computation result.
#' @keywords internal
#' @noRd
gllrm_candidate_base_model <- function(object) {
  if (inherits(object, "gRm_model")) {
    return(object)
  }
  if (inherits(object, "gRm_fit")) {
    return(object$model %||% object$spec)
  }
  if (inherits(object, "gRm_project") ||
      (is.list(object) && is.data.frame(object$items))) {
    analysis <- list(
      project = object,
      items = as.character(object$items$name),
      exogenous = as.character(object$backgrounds$name %||% character()),
      source_trace = object$source_trace %||% character()
    )
    class(analysis) <- c("gRm_analysis", "list")
    return(new_gRm_model(
      analysis = analysis,
      ld = empty_ld_terms(),
      dif = empty_dif_terms(),
      call = match.call()
    ))
  }
  gllrm(object)
}

#' Internal gllrm candidate ld spec helper
#'
#' Supports the gllrm candidate fit implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param object Object dispatched to this helper.
#' @param item1 Internal `item1` value used by this helper.
#' @param item2 Internal `item2` value used by this helper.
#' @return The internal `gllrm_candidate_ld_spec()` computation result.
#' @keywords internal
#' @noRd
gllrm_candidate_ld_spec <- function(object, item1, item2) {
  model <- gllrm_candidate_base_model(object)
  items <- model$project$items
  candidate <- canonical_ld_term(
    vars = c(items$name[[item1]], items$name[[item2]]),
    items = items$name,
    label = paste(items$name[[item1]], items$name[[item2]], sep = ":"),
    source = "user"
  )
  ld_terms <- source_order_ld_table(
    items,
    rbind_fill(model$ld %||% empty_ld_terms(), candidate)
  )
  dif_terms <- source_order_dif_table(
    items,
    model$project$backgrounds,
    model$dif %||% empty_dif_terms()
  )
  new_gRm_model_from_canonical_terms(
    model,
    ld = ld_terms,
    dif = dif_terms,
    call = match.call()
  )
}

#' Internal gllrm candidate dif spec helper
#'
#' Supports the gllrm candidate fit implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param object Object dispatched to this helper.
#' @param item One-based item index.
#' @param background One-based exogenous-variable index.
#' @return The internal `gllrm_candidate_dif_spec()` computation result.
#' @keywords internal
#' @noRd
gllrm_candidate_dif_spec <- function(object, item, background) {
  model <- gllrm_candidate_base_model(object)
  items <- model$project$items
  backgrounds <- model$project$backgrounds
  candidate <- canonical_dif_term(
    vars = c(items$name[[item]], backgrounds$name[[background]]),
    items = items$name,
    exogenous = backgrounds$name,
    label = paste(items$name[[item]], backgrounds$name[[background]], sep = ":"),
    source = "user"
  )
  ld_terms <- source_order_ld_table(
    items,
    model$ld %||% empty_ld_terms()
  )
  dif_terms <- source_order_dif_table(
    items,
    backgrounds,
    rbind_fill(model$dif %||% empty_dif_terms(), candidate)
  )
  new_gRm_model_from_canonical_terms(
    model,
    ld = ld_terms,
    dif = dif_terms,
    call = match.call()
  )
}

#' Internal fit gllrm candidate core helper
#'
#' Supports the gllrm candidate fit implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param spec GLLRM model specification.
#' @param bundle Source-shaped analysis bundle.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @return The internal `fit_gllrm_candidate_core()` computation result.
#' @keywords internal
#' @noRd
fit_gllrm_candidate_core <- function(spec, bundle, max_step, max_delta) {
  fit_gllrm(
    spec,
    max_step = max_step,
    max_delta = max_delta,
    bundle = bundle
  )
}

#' Internal new gllrm candidate fit helper
#'
#' Supports the gllrm candidate fit implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param gllrm_fit Internal `gllrm_fit` value used by this helper.
#' @param spec GLLRM model specification.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_gllrm_candidate_fit <- function(gllrm_fit, spec, max_step, max_delta) {
  convergence <- gRm_gllrm_convergence_state(
    gllrm_fit$state,
    max_step = max_step,
    max_delta = max_delta
  )
  out <- list(
    spec = spec,
    context = gllrm_fit$context,
    state = gllrm_fit$state,
    log_likelihood = gllrm_fit$state$log_likelihood,
    convergence = convergence,
    n_step = convergence$iterations,
    report_delta = convergence$report_delta,
    final_delta = convergence$final_delta,
    converged = convergence$converged,
    stop_reason = convergence$stop_reason,
    max_step = max_step,
    max_delta = max_delta
  )
  class(out) <- c("gRm_gllrm_candidate_fit", "list")
  out
}

#' Internal fit gllrm candidate bundle helper
#'
#' Supports the gllrm candidate fit implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param object Object dispatched to this helper.
#' @param spec GLLRM model specification.
#' @return The internal `fit_gllrm_candidate_bundle()` computation result.
#' @keywords internal
#' @noRd
fit_gllrm_candidate_bundle <- function(object, spec) {
  object$bundle %||%
    object$fit$bundle %||%
    build_item_parameters_bundle(spec$project)
}

#' Internal fit gllrm candidate ld helper
#'
#' Supports the gllrm candidate fit implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param object Object dispatched to this helper.
#' @param item1 Internal `item1` value used by this helper.
#' @param item2 Internal `item2` value used by this helper.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @return The internal `fit_gllrm_candidate_ld()` computation result.
#' @keywords internal
#' @noRd
fit_gllrm_candidate_ld <- function(object, item1, item2, max_step, max_delta) {
  spec <- gllrm_candidate_ld_spec(object, item1 = item1, item2 = item2)
  bundle <- fit_gllrm_candidate_bundle(object, spec)
  gllrm_fit <- fit_gllrm_candidate_core(
    spec,
    bundle = bundle,
    max_step = max_step,
    max_delta = max_delta
  )
  new_gllrm_candidate_fit(gllrm_fit, spec, max_step = max_step, max_delta = max_delta)
}

#' Internal fit gllrm candidate dif helper
#'
#' Supports the gllrm candidate fit implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/skbias12b.pas::Estimate_GLLRM`.
#' @param object Object dispatched to this helper.
#' @param item One-based item index.
#' @param background One-based exogenous-variable index.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @return The internal `fit_gllrm_candidate_dif()` computation result.
#' @keywords internal
#' @noRd
fit_gllrm_candidate_dif <- function(object, item, background, max_step, max_delta) {
  spec <- gllrm_candidate_dif_spec(object, item = item, background = background)
  bundle <- fit_gllrm_candidate_bundle(object, spec)
  gllrm_fit <- fit_gllrm_candidate_core(
    spec,
    bundle = bundle,
    max_step = max_step,
    max_delta = max_delta
  )
  new_gllrm_candidate_fit(gllrm_fit, spec, max_step = max_step, max_delta = max_delta)
}
