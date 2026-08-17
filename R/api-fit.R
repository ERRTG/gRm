#' Fit a DIGRAM model
#'
#' Estimate the parameters of a model specified by [gllrm()].
#'
#' @usage
#' fit(object, ...)
#' \method{fit}{gRm_model}(object, max_step = 5000L, max_delta = 0.0001, ...)
#' @aliases fit.gRm_model fit.gRm_analysis fit.gRm_screen fit.default
#' @param object A `gRm_model` object created by [gllrm()]. Calling
#'   `fit()` directly on an analysis object is an error; specify the model
#'   first.
#' @param max_step Single positive integer-like maximum number of fitting
#'   iterations.
#' @param max_delta Single positive finite strict convergence tolerance for the
#'   maximum absolute observed-versus-fitted sufficient-count discrepancy.
#' @param ... Reserved for S3 dispatch compatibility; no additional fitting
#'   arguments are currently supported.
#' @return A `gRm_fit` object containing the fitted model, convergence
#'   metadata, source-shaped sufficient statistics, and result values used by
#'   downstream accessors.
#' @details
#' For Rasch models, `fit()` uses the source-faithful base Rasch fitting path.
#' For models with LD or DIF terms, it uses the GLLRM fitting path. The
#' returned fit object is the input to [item_fit()], [local_dependence()],
#' [dif()], and [global_homogeneity()].
#'
#' Printing the fitted object shows compact fit status. `summary(fit)` shows
#' fitted item parameters, thresholds, and fitted LD/DIF terms. Use
#' `summary(fit, which = "parameters")` for item-parameter summaries only and
#' `summary(fit, which = "thresholds")` for threshold rows only.
#' @examples
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0)
#' )
#' analysis <- gRm(data, items = c("I1", "I2"), id = "ID", score_cuts = c(1L, 2L))
#' model <- gllrm(analysis)
#' fitted <- fit(model)
#' fitted
#' summary(fitted)
#' @seealso [gllrm()], [item_fit()], [local_dependence()], [dif()],
#'   [global_homogeneity()]
#' @export
fit <- function(object, ...) {
  UseMethod("fit")
}

#' @export
fit.gRm_analysis <- function(object,
                                ld = NULL,
                                dif = NULL,
                                max_step = 5000L,
                                max_delta = 0.0001,
                                ...) {
  stop("fit() requires a DIGRAM model specification created by gllrm().", call. = FALSE)
}

#' @export
fit.gRm_model <- function(object,
                             max_step = 5000L,
                             max_delta = 0.0001,
                             ...) {
  reject_public_dots(...)
  controls <- normalize_public_fit_controls(max_step, max_delta)
  max_step <- controls$max_step
  max_delta <- controls$max_delta
  bundle <- build_item_parameters_bundle(object$analysis$project)
  assert_gRm_analysis_identity(object$analysis, bundle)
  assert_public_estimable_fit_bundle(bundle, object)
  if (nrow(object$ld) > 0L || nrow(object$dif) > 0L) {
    gllrm_fit <- fit_gllrm(object, max_step = max_step, max_delta = max_delta, bundle = bundle)
    output_state <- gllrm_output_parameter_state(gllrm_fit$context, gllrm_fit$state)
    output_fit <- list(context = gllrm_fit$context, state = output_state, bundle = gllrm_fit$bundle)
    values <- gllrm_values(output_fit, object)
    assert_public_fit_values(values)
    return(new_gRm_fit(
      model = object,
      bundle = gllrm_fit$bundle,
      fit = c(gllrm_fit$state, list(context = gllrm_fit$context)),
      values = values,
      convergence = gRm_gllrm_convergence_state(
        gllrm_fit$state,
        max_step = max_step,
        max_delta = max_delta
      ),
      source_trace = c(object$source_trace %||% character(), fit = "fit_gllrm"),
      call = match.call()
    ))
  }

  base_fit <- fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)
  values <- item_parameters_values(base_fit, bundle)
  assert_public_fit_values(values)
  new_gRm_fit(
    model = object,
    bundle = bundle,
    fit = base_fit,
    values = values,
    convergence = gRm_rasch_convergence_state(
      base_fit,
      max_step = max_step,
      max_delta = max_delta
    ),
    source_trace = c(object$source_trace %||% character(), fit = "fit_rasch_base"),
    call = match.call()
  )
}

#' @export
fit.gRm_screen <- function(object, ...) {
  stop(
    "fit() requires a GLLRM model specification. Use gllrm(screen_obj) to create one from selected screening terms.",
    call. = FALSE
  )
}

#' @export
fit.default <- function(object, ...) {
  stop("fit() requires a DIGRAM model specification created by gllrm().", call. = FALSE)
}

#' Internal normalize public fit controls helper
#'
#' Supports the api fit implementation while preserving its internal contract.
#' @param max_step Maximum fitting iteration.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
normalize_public_fit_controls <- function(max_step, max_delta) {
  list(
    max_step = normalize_public_max_step(max_step),
    max_delta = normalize_public_max_delta(max_delta)
  )
}

#' Internal normalize public max step helper
#'
#' Supports the api fit implementation while preserving its internal contract.
#' @param max_step Maximum fitting iteration.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
normalize_public_max_step <- function(max_step) {
  if (
    length(max_step) != 1L ||
      !is.numeric(max_step) ||
      is.na(max_step) ||
      !is.finite(max_step) ||
      max_step <= 0 ||
      max_step != floor(max_step) ||
      max_step > .Machine$integer.max
  ) {
    stop("`max_step` must be a single positive integer-like value.", call. = FALSE)
  }
  as.integer(max_step)
}

#' Internal normalize public max delta helper
#'
#' Supports the api fit implementation while preserving its internal contract.
#' @param max_delta Sufficient-count discrepancy tolerance.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
normalize_public_max_delta <- function(max_delta) {
  if (
    length(max_delta) != 1L ||
      !is.numeric(max_delta) ||
      is.na(max_delta) ||
      !is.finite(max_delta) ||
      max_delta <= 0
  ) {
    stop("`max_delta` must be a single positive finite number.", call. = FALSE)
  }
  as.numeric(max_delta)
}

#' Internal new gRm fit helper
#'
#' Supports the api fit implementation while preserving its internal contract.
#' @param model gRm model specification.
#' @param bundle Source-shaped analysis bundle.
#' @param fit Fitted gRm model.
#' @param values Values to validate or transform.
#' @param convergence Internal `convergence` value used by this helper.
#' @param source_trace Internal `source_trace` value used by this helper.
#' @param call Captured R call.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_gRm_fit <- function(model, bundle, fit, values, convergence, source_trace, call) {
  out <- list(
    analysis = model$analysis,
    project = model$analysis$project,
    analysis_fingerprint = model$analysis_fingerprint %||% model$analysis$analysis_fingerprint,
    likelihood_sample = model$likelihood_sample %||% model$analysis$likelihood_sample,
    model = model,
    model_type = model$model_type,
    bundle = bundle,
    fit = fit,
    values = values,
    parameters = values,
    convergence = convergence,
    source_trace = source_trace,
    unmodeled = character(),
    warnings = character(),
    call = call
  )
  out$spec <- out$model
  class(out) <- c("gRm_fit", "list")
  out
}

#' Internal stop gRm not implemented helper
#'
#' Supports the api fit implementation while preserving its internal contract.
#' @param feature Internal `feature` value used by this helper.
#' @param validation_needed Internal `validation_needed` value used by this helper.
#' @return The internal `stop_gRm_not_implemented()` computation result.
#' @keywords internal
#' @noRd
stop_gRm_not_implemented <- function(feature, validation_needed = NULL) {
  message <- paste0(feature, " is not implemented in the source-faithful R API yet.")
  if (!is.null(validation_needed)) {
    message <- paste0(message, " Required before enabling: ", validation_needed, ".")
  }
  condition <- structure(
    list(
      message = message,
      feature = feature,
      validation_needed = validation_needed
    ),
    class = c("gRm_not_implemented", "error", "condition")
  )
  stop(condition)
}
