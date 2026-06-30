# GLLRM candidate fits for LD and DIF diagnostics.
#
# Source trace: CHECK LID and CHECK DIF candidate tests compare the current
# current GLLRM with a source-ordered model containing exactly one additional
# IJ or IX term. These helpers run the same GLLRM estimation core used by
# public fit(), but skip public item-parameter output construction that the
# candidate likelihood-ratio tests do not use.

gllrm_candidate_ld_spec <- function(object, item1, item2) {
  model <- object$model %||% object$spec
  context <- object$fit$context
  candidate <- canonical_ld_term(
    vars = c(context$items$name[[item1]], context$items$name[[item2]]),
    items = context$items$name,
    label = paste(context$items$name[[item1]], context$items$name[[item2]], sep = ":"),
    source = "user"
  )
  ld_terms <- source_order_ld_table(
    context$items,
    rbind_fill(model$ld %||% empty_ld_terms(), candidate)
  )
  dif_terms <- source_order_dif_table(
    context$items,
    context$backgrounds,
    model$dif %||% empty_dif_terms()
  )
  new_gRm_model_from_canonical_terms(
    model,
    ld = ld_terms,
    dif = dif_terms,
    call = match.call()
  )
}

gllrm_candidate_dif_spec <- function(object, item, background) {
  model <- object$model %||% object$spec
  context <- object$fit$context
  candidate <- canonical_dif_term(
    vars = c(context$items$name[[item]], context$backgrounds$name[[background]]),
    items = context$items$name,
    exogenous = context$backgrounds$name,
    label = paste(context$items$name[[item]], context$backgrounds$name[[background]], sep = ":"),
    source = "user"
  )
  ld_terms <- source_order_ld_table(
    context$items,
    model$ld %||% empty_ld_terms()
  )
  dif_terms <- source_order_dif_table(
    context$items,
    context$backgrounds,
    rbind_fill(model$dif %||% empty_dif_terms(), candidate)
  )
  new_gRm_model_from_canonical_terms(
    model,
    ld = ld_terms,
    dif = dif_terms,
    call = match.call()
  )
}

fit_gllrm_candidate_core <- function(spec, bundle, max_step, max_delta) {
  fit_gllrm(
    spec,
    max_step = max_step,
    max_delta = max_delta,
    bundle = bundle
  )
}

new_gllrm_candidate_fit <- function(gllrm_fit, spec, max_step, max_delta) {
  out <- list(
    spec = spec,
    context = gllrm_fit$context,
    state = gllrm_fit$state,
    log_likelihood = gllrm_fit$state$log_likelihood,
    n_step = gllrm_fit$state$n_step,
    report_delta = gllrm_fit$state$report_delta,
    converged = gllrm_fit$state$converged,
    stop_reason = gllrm_fit$state$stop_reason,
    max_step = max_step,
    max_delta = max_delta
  )
  class(out) <- c("gRm_gllrm_candidate_fit", "list")
  out
}

fit_gllrm_candidate_bundle <- function(object, spec) {
  object$bundle %||%
    object$fit$bundle %||%
    build_item_parameters_bundle(spec$project)
}

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
