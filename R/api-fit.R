#' Fit a DIGRAM model
#'
#' Estimate the parameters of a model specified by [gllrm()].
#'
#' @usage
#' fit(object, ...)
#' \method{fit}{gRm_model}(object, max_step = 5000L, max_delta = 0.0001, ...)
#' @aliases fit.gRm_model
#' @param object A `gRm_model` object created by [gllrm()]. Calling
#'   `fit()` directly on an analysis object is an error; specify the model
#'   first.
#' @param max_step Single positive integer-like maximum number of fitting
#'   iterations.
#' @param max_delta Single positive finite convergence tolerance for the maximum
#'   absolute parameter update.
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
  controls <- normalize_public_fit_controls(max_step, max_delta)
  max_step <- controls$max_step
  max_delta <- controls$max_delta
  bundle <- build_item_parameters_bundle(object$analysis$project)
  assert_public_estimable_fit_bundle(bundle)
  if (nrow(object$ld) > 0L || nrow(object$dif) > 0L) {
    gllrm_fit <- fit_gllrm(object, max_step = max_step, max_delta = max_delta, bundle = bundle)
    output_state <- gllrm_output_parameter_state(gllrm_fit$context, gllrm_fit$state)
    output_fit <- list(context = gllrm_fit$context, state = output_state, bundle = gllrm_fit$bundle)
    values <- gllrm_values(output_fit, object)
    return(new_gRm_fit(
      model = object,
      bundle = gllrm_fit$bundle,
      fit = c(gllrm_fit$state, list(context = gllrm_fit$context)),
      values = values,
      convergence = list(
        converged = gllrm_fit$state$converged,
        iterations = gllrm_fit$state$n_step,
        max_step = max_step,
        max_delta = max_delta
      ),
      source_trace = c(object$source_trace %||% character(), fit = "fit_gllrm"),
      call = match.call()
    ))
  }

  base_fit <- fit_rasch_base(bundle, max_step = max_step, max_delta = max_delta)
  values <- item_parameters_values(base_fit, bundle)
  new_gRm_fit(
    model = object,
    bundle = bundle,
    fit = base_fit,
    values = values,
    convergence = list(
      converged = base_fit$converged %||% values$converged %||% NA,
      iterations = base_fit$n_step %||% base_fit$iteration %||% values$iteration %||% NA,
      max_step = max_step,
      max_delta = max_delta
    ),
    source_trace = c(object$source_trace %||% character(), fit = "fit_rasch_base"),
    call = match.call()
  )
}

assert_public_estimable_fit_bundle <- function(bundle) {
  # Source trace: source/PAS_skunits/skbias21.pas exits report sections when
  # Nvalid < 1. The R fitting API should therefore not manufacture a converged
  # model object when the source-shaped estimating margins are empty.
  if (isTRUE(bundle$manifest$nvalid < 1L)) {
    stop(
      paste(
        "Rasch fitting requires at least one source-valid complete response pattern",
        "inside the DIGRAM score window."
      ),
      call. = FALSE
    )
  }
  invisible(bundle)
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

normalize_public_fit_controls <- function(max_step, max_delta) {
  list(
    max_step = normalize_public_max_step(max_step),
    max_delta = normalize_public_max_delta(max_delta)
  )
}

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

new_gRm_fit <- function(model, bundle, fit, values, convergence, source_trace, call) {
  out <- list(
    analysis = model$analysis,
    project = model$analysis$project,
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
