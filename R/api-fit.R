#' Fit a DIGRAM model
#'
#' Estimate the parameters of a model specified by [gllrm()].
#'
#' @usage fit(object, max_step = 5000L, max_delta = 0.0001, ...)
#' @param object A `gRm_model` object created by [gllrm()]. Calling
#'   `fit()` directly on an analysis object is an error; specify the model
#'   first.
#' @param max_step Maximum number of fitting iterations.
#' @param max_delta Convergence tolerance for the maximum absolute parameter
#'   update.
#' @param ... Reserved for S3 dispatch compatibility; no additional fitting
#'   arguments are currently supported.
#' @return A `gRm_fit` object containing the fitted model, convergence
#'   metadata, source-shaped sufficient statistics, and result values used by
#'   downstream accessors.
#' @details
#' For Rasch models, `fit()` uses the source-faithful base Rasch fitting path.
#' For models with LD or DIF terms, it uses the active GLLRM fitting path. The
#' returned fit object is the input to [item_parameters()], [item_fit()],
#' [local_dependence()], [dif()], and [global_homogeneity()].
#'
#' Use `summary(fit, which = "fit")` for convergence and likelihood summaries,
#' `summary(fit, which = "parameters")` for item-parameter summaries, and
#' `summary(fit, which = "terms")` for fitted model terms.
#' @examples
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0)
#' )
#' analysis <- gRm(data, items = c("I1", "I2"), id = "ID")
#' model <- gllrm(analysis)
#' fitted <- fit(model, max_step = 100L)
#' summary(fitted, which = "fit")
#' @seealso [gllrm()], [item_parameters()], [item_fit()],
#'   [local_dependence()], [dif()], [global_homogeneity()]
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
fit.gRm_item_analysis <- function(object,
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
  if (nrow(object$ld) > 0L || nrow(object$dif) > 0L) {
    active_fit <- fit_gllrm_active(object, max_step = max_step, max_delta = max_delta)
    output_state <- gllrm_output_parameter_state(active_fit$context, active_fit$state)
    output_fit <- list(context = active_fit$context, state = output_state, bundle = active_fit$bundle)
    values <- gllrm_active_values(output_fit, object)
    return(new_gRm_fit(
      model = object,
      bundle = active_fit$bundle,
      fit = c(active_fit$state, list(context = active_fit$context)),
      values = values,
      convergence = list(
        converged = active_fit$state$converged,
        iterations = active_fit$state$n_step,
        max_step = max_step,
        max_delta = max_delta
      ),
      source_trace = c(object$source_trace %||% character(), fit = "fit_gllrm_active"),
      call = match.call()
    ))
  }

  bundle <- build_item_parameters_bundle(object$analysis$project)
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

#' @export
fit.gRm_gllrm_spec <- function(object,
                                  max_step = 5000L,
  max_delta = 0.0001,
                                  ...) {
  fit.gRm_model(object, max_step = max_step, max_delta = max_delta, ...)
}

fit.gRm_screen <- function(object, ...) {
  stop(
    "fit() requires a GLLRM model specification. Use gllrm(screen_obj) to create one from selected screening terms.",
    call. = FALSE
  )
}

fit.default <- function(object, ...) {
  stop("fit() requires a DIGRAM model specification created by gllrm().", call. = FALSE)
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

#' @export
print.gRm_fit <- function(x, ...) {
  label <- if (identical(x$model_type, "rasch")) "DIGRAM Rasch fit" else "DIGRAM GLLRM fit"
  model <- x$model %||% x$spec
  cat("<", label, ">\n", sep = "")
  cat("  items: ", length(x$analysis$items), "\n", sep = "")
  cat("  local dependence terms: ", nrow(model$ld), "\n", sep = "")
  cat("  DIF terms: ", nrow(model$dif), "\n", sep = "")
  cat("  converged: ", as.character(x$convergence$converged), "\n", sep = "")
  cat("  use summary() for results\n", sep = "")
  invisible(x)
}

#' @export
print.gRm_gllrm_fit <- function(x, ...) {
  print.gRm_fit(x, ...)
}
