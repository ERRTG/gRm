#' Likelihood methods for fitted gRm models
#'
#' Compare fitted gRm models with a likelihood-ratio table that follows the
#' shape of R model-comparison output while using likelihood-ratio terminology.
#'
#' @name gRm-likelihood
#' @aliases logLik.gRm_fit anova.gRm_fit
#' @usage
#' \method{logLik}{gRm_fit}(object, ...)
#' \method{anova}{gRm_fit}(object, ..., test = "Chisq")
#' @param object A fitted model returned by [fit()].
#' @param ... For `anova.gRm_fit()`, additional `gRm_fit` objects fitted to the
#'   same gRm analysis. For `logLik.gRm_fit()`, reserved for S3 compatibility
#'   and must be empty.
#' @param test Test column to include in the model-comparison table. `"Chisq"`
#'   and `"LRT"` request the standard likelihood-ratio chi-square p-value
#'   computed with [stats::pchisq()]. Use `NULL`, `"none"`, or `FALSE` to omit
#'   the p-value column.
#' @return `logLik.gRm_fit()` returns a standard `"logLik"` object using R's
#'   sign convention. gRm stores DIGRAM's source-style negative log likelihood
#'   internally: if the fitted object stores `L = -log(likelihood)`,
#'   `logLik(fit)` returns `-L`.
#'
#'   `anova.gRm_fit()` returns an `"anova"` data frame with columns:
#'   `"Model Df"`, `"-2 logLik"`, `"Df"`, `"Chisq"`, and, when requested,
#'   `"Pr(>Chisq)"`.
#' @details
#' gRm stores the source-faithful DIGRAM quantity `-log(likelihood)` on fitted
#'   objects. `logLik.gRm_fit()` converts that value to R's standard log
#'   likelihood sign, matching the sign convention used by [stats::logLik()].
#'   `anova.gRm_fit()` reports `-2 logLik` for each model and the
#'   likelihood-ratio change between consecutive rows. Its printed format
#'   follows the same style as `anova(glm_fit0, glm_fit1, test = "Chisq")`,
#'   but uses likelihood-ratio column names instead of residual-deviance names.
#'   The test statistic remains the source-faithful likelihood-ratio statistic,
#'   while the p-value is computed with R's standard
#'   `stats::pchisq(..., lower.tail = FALSE)` when the row has a positive
#'   degrees-of-freedom increase and non-negative likelihood-ratio statistic.
#'   Reversed or otherwise invalid comparison rows keep the signed diagnostic
#'   `Df` and `Chisq` values, but `Pr(>Chisq)` is set to `NA`.
#'
#' The supplied models must have nested LD/DIF term sets and must be fitted to
#'   exactly the same analysis and likelihood rows. The comparison checks the
#'   deterministic analysis fingerprint, its canonical payload, the exact
#'   likelihood row mask, and term-set inclusion before calculating a table.
#' @examples
#' data <- expand.grid(I1 = 1:3, I2 = 1:3, I3 = 1:3)
#' data$ID <- seq_len(nrow(data))
#' data$site <- rep(1:3, length.out = nrow(data))
#' analysis <- gRm(
#'   data,
#'   items = c("I1", "I2", "I3"),
#'   exogenous = "site",
#'   id = "ID",
#'   score_cuts = "auto"
#' )
#' fit0 <- fit(gllrm(analysis))
#' fit1 <- fit(gllrm(analysis, ld = ~ I1:I2))
#'
#' # gRm stores DIGRAM's negative log likelihood internally, but logLik()
#' # returns the usual R signed log likelihood.
#' logLik(fit0)
#'
#' # Compare nested fitted models with an anova-style likelihood-ratio table.
#' anova(fit0, fit1, test = "Chisq")
#' @seealso [fit()], [gllrm()]
NULL

#' @export
logLik.gRm_fit <- function(object, ...) {
  reject_public_dots(...)
  info <- gRm_fit_likelihood_info(object)
  out <- -info$negative_log_likelihood
  attr(out, "df") <- info$n_parameters
  attr(out, "nobs") <- info$likelihood_n
  class(out) <- "logLik"
  out
}

#' @export
anova.gRm_fit <- function(object, ..., test = "Chisq") {
  fits <- list(object, ...)
  if (length(fits) < 2L) {
    stop("anova.gRm_fit() requires at least two fitted gRm models.", call. = FALSE)
  }
  if (!all(vapply(fits, inherits, logical(1L), what = "gRm_fit"))) {
    stop("All objects supplied to anova.gRm_fit() must be fitted gRm models.", call. = FALSE)
  }
  validate_gRm_fit_likelihood_comparison(fits)

  log_likelihoods <- lapply(fits, stats::logLik)
  model_df <- vapply(log_likelihoods, attr, numeric(1L), which = "df")
  neg2_loglik <- -2 * vapply(log_likelihoods, as.numeric, numeric(1L))
  delta_df <- c(NA_real_, diff(model_df))
  lr_statistic <- c(NA_real_, -diff(neg2_loglik))

  out <- data.frame(
    check.names = FALSE,
    "Model Df" = model_df,
    "-2 logLik" = neg2_loglik,
    "Df" = delta_df,
    "Chisq" = lr_statistic
  )
  row.names(out) <- as.character(seq_len(nrow(out)))

  if (!is.null(normalize_gRm_likelihood_test(test))) {
    p_value <- rep(NA_real_, length(fits))
    test_rows <- seq_along(fits)[-1L]
    # Source-defined DIGRAM reports compute non-negative LR statistics for
    # specific current/previous or candidate-addition comparisons; see
    # source/PAS_scd/DGRirtD.pas. anova.gRm_fit() accepts arbitrary R row order,
    # so keep signed row diagnostics and suppress only invalid p-values.
    usable <- !is.na(delta_df[test_rows]) &
      !is.na(lr_statistic[test_rows]) &
      delta_df[test_rows] > 0 &
      lr_statistic[test_rows] >= 0
    p_value[test_rows[usable]] <- stats::pchisq(
      lr_statistic[test_rows[usable]],
      df = delta_df[test_rows[usable]],
      lower.tail = FALSE
    )
    if (any(!usable)) {
      warning(
        "Some model comparisons have non-positive df changes or negative likelihood-ratio statistics; p-values set to NA.",
        call. = FALSE
      )
    }
    out[["Pr(>Chisq)"]] <- p_value
  }

  class(out) <- c("anova", "data.frame")
  attr(out, "heading") <- c(
    "Likelihood ratio tests for gRm fits\n",
    paste(vapply(seq_along(fits), gRm_fit_anova_heading, character(1L), fits = fits), collapse = "\n")
  )
  out
}

#' Internal gRm fit likelihood info helper
#'
#' Supports the api likelihood comparison implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @return The internal `gRm_fit_likelihood_info()` computation result.
#' @keywords internal
#' @noRd
gRm_fit_likelihood_info <- function(object) {
  values <- object$values %||% list()
  negative_log_likelihood <- values$negative_log_likelihood %||%
    values$log_likelihood %||%
    object$fit$negative_log_likelihood %||%
    object$fit$log_likelihood %||%
    NA_real_
  n_parameters <- values$n_parameters %||% object$fit$n_parameters %||% NA_real_
  likelihood_n <- values$likelihood_n %||% object$fit$likelihood_n %||% NA_real_

  if (
    length(negative_log_likelihood) != 1L ||
      !is.numeric(negative_log_likelihood) ||
      is.na(negative_log_likelihood) ||
      !is.finite(negative_log_likelihood)
  ) {
    stop("Fitted gRm model does not contain a finite negative log-likelihood value.", call. = FALSE)
  }
  if (
    length(n_parameters) != 1L ||
      !is.numeric(n_parameters) ||
      is.na(n_parameters) ||
      !is.finite(n_parameters) ||
      n_parameters < 0 ||
      n_parameters != floor(n_parameters)
  ) {
    stop("Fitted gRm model does not contain a non-negative integer-like parameter count.", call. = FALSE)
  }
  if (
    length(likelihood_n) != 1L ||
      !is.numeric(likelihood_n) ||
      is.na(likelihood_n) ||
      !is.finite(likelihood_n) ||
      likelihood_n < 0 ||
      likelihood_n != floor(likelihood_n)
  ) {
    likelihood_n <- NA_real_
  }

  list(
    negative_log_likelihood = as.numeric(negative_log_likelihood),
    n_parameters = as.numeric(n_parameters),
    likelihood_n = as.numeric(likelihood_n)
  )
}

#' Internal validate gRm fit likelihood comparison helper
#'
#' Supports the api likelihood comparison implementation while preserving its internal contract.
#' @param fits Internal `fits` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
validate_gRm_fit_likelihood_comparison <- function(fits) {
  fingerprints <- lapply(fits, gRm_fit_analysis_fingerprint)
  if (any(vapply(fingerprints, is.null, logical(1L)))) {
    stop(
      "Fitted gRm models must contain an analysis fingerprint to be compared.",
      call. = FALSE
    )
  }
  reference_fingerprint <- fingerprints[[1L]]
  same_fingerprint <- vapply(
    fingerprints[-1L],
    identical,
    logical(1L),
    y = reference_fingerprint
  )
  if (!all(same_fingerprint)) {
    stop("Fitted gRm models must use the same analysis data and encoding to be compared.", call. = FALSE)
  }

  identities <- lapply(fits, gRm_fit_analysis_identity)
  known_identities <- !vapply(identities, is.null, logical(1L))
  if (any(known_identities) && !all(known_identities)) {
    stop("Fitted gRm models have incomplete analysis identity metadata.", call. = FALSE)
  }
  if (all(known_identities)) {
    reference_identity <- identities[[1L]]
    same_identity <- vapply(
      identities[-1L],
      identical,
      logical(1L),
      y = reference_identity
    )
    if (!all(same_identity)) {
      stop("Fitted gRm models must use identical canonical analysis inputs.", call. = FALSE)
    }
  }

  samples <- lapply(fits, gRm_fit_likelihood_sample)
  if (any(vapply(samples, is.null, logical(1L)))) {
    stop(
      "Fitted gRm models must contain an exact likelihood row mask to be compared.",
      call. = FALSE
    )
  }
  reference_sample <- samples[[1L]]
  same_sample <- vapply(samples[-1L], identical, logical(1L), y = reference_sample)
  if (!all(same_sample)) {
    stop("Fitted gRm models must use the same likelihood rows to be compared.", call. = FALSE)
  }

  sample_sizes <- vapply(fits, function(fit) {
    gRm_fit_likelihood_info(fit)$likelihood_n
  }, numeric(1L))
  known <- sample_sizes[!is.na(sample_sizes)]
  if (length(known) > 1L && any(known != known[[1L]])) {
    stop("Fitted gRm models must use the same likelihood sample size.", call. = FALSE)
  }

  sample_counts <- vapply(samples, function(sample) sample$n, integer(1L))
  if (any(!is.na(sample_sizes) & sample_sizes != sample_counts)) {
    stop("Fitted gRm model likelihood row metadata disagrees with its sample size.", call. = FALSE)
  }

  term_sets <- lapply(fits, gRm_fit_term_set)
  for (index in seq.int(2L, length(term_sets))) {
    previous <- term_sets[[index - 1L]]
    current <- term_sets[[index]]
    nested <- all(previous %in% current) || all(current %in% previous)
    if (!nested) {
      stop(
        "Consecutive fitted gRm models must have nested LD/DIF term sets; models ",
        index - 1L,
        " and ",
        index,
        " are non-nested.",
        call. = FALSE
      )
    }
  }

  invisible(fits)
}

#' Extract a fitted analysis fingerprint
#'
#' @param fit A `gRm_fit` object.
#' @return Scalar fingerprint or `NULL`.
#' @keywords internal
gRm_fit_analysis_fingerprint <- function(fit) {
  analysis <- fit$analysis %||% (fit$model %||% fit$spec)$analysis %||% NULL
  fingerprint <- fit$analysis_fingerprint %||% analysis$analysis_fingerprint %||% NULL
  if (
    is.null(fingerprint) ||
      !is.character(fingerprint) ||
      length(fingerprint) != 1L ||
      is.na(fingerprint) ||
      !nzchar(fingerprint)
  ) {
    return(NULL)
  }
  fingerprint
}

#' Extract a canonical fitted analysis identity
#'
#' @param fit A `gRm_fit` object.
#' @return Canonical identity list or `NULL`.
#' @keywords internal
gRm_fit_analysis_identity <- function(fit) {
  analysis <- fit$analysis %||% (fit$model %||% fit$spec)$analysis %||% NULL
  if (!is.list(analysis)) {
    return(NULL)
  }
  analysis$analysis_identity %||% NULL
}

#' Extract and validate fitted likelihood-sample identity
#'
#' @param fit A `gRm_fit` object.
#' @return Likelihood-sample identity list or `NULL`.
#' @keywords internal
gRm_fit_likelihood_sample <- function(fit) {
  analysis <- fit$analysis %||% (fit$model %||% fit$spec)$analysis %||% NULL
  sample <- fit$likelihood_sample %||% analysis$likelihood_sample %||% NULL
  required <- c("schema", "fingerprint", "row_mask", "row_indices", "n")
  if (!is.list(sample) || !all(required %in% names(sample))) {
    return(NULL)
  }
  if (
    !is.logical(sample$row_mask) ||
      anyNA(sample$row_mask) ||
      !identical(which(sample$row_mask), sample$row_indices) ||
      !identical(as.integer(sum(sample$row_mask)), sample$n)
  ) {
    return(NULL)
  }
  sample
}

#' Encode a model-term component without delimiter ambiguity
#'
#' @param x Scalar term component.
#' @return Length-prefixed UTF-8 string.
#' @keywords internal
gRm_fit_term_component <- function(x) {
  x <- enc2utf8(as.character(x))
  paste0(nchar(x, type = "bytes"), ":", x)
}

#' Extract the canonical LD/DIF term set from a fit
#'
#' @param fit A `gRm_fit` object.
#' @return Sorted unique character term keys.
#' @keywords internal
gRm_fit_term_set <- function(fit) {
  spec <- fit$model %||% fit$spec %||% list()
  ld <- spec$ld %||% data.frame()
  dif <- spec$dif %||% data.frame()
  keys <- character()
  if (is.data.frame(ld) && nrow(ld)) {
    keys <- c(keys, vapply(seq_len(nrow(ld)), function(index) {
      paste0(
        "LD/",
        gRm_fit_term_component(ld$item1[[index]]),
        "/",
        gRm_fit_term_component(ld$item2[[index]])
      )
    }, character(1L)))
  }
  if (is.data.frame(dif) && nrow(dif)) {
    keys <- c(keys, vapply(seq_len(nrow(dif)), function(index) {
      paste0(
        "DIF/",
        gRm_fit_term_component(dif$item[[index]]),
        "/",
        gRm_fit_term_component(dif$exogenous[[index]])
      )
    }, character(1L)))
  }
  sort(unique(keys))
}

#' Internal normalize gRm likelihood test helper
#'
#' Supports the api likelihood comparison implementation while preserving its internal contract.
#' @param test Internal `test` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
normalize_gRm_likelihood_test <- function(test) {
  if (is.null(test) || identical(test, FALSE)) {
    return(NULL)
  }
  if (!is.character(test) || length(test) != 1L || is.na(test)) {
    stop("`test` must be \"Chisq\", \"LRT\", \"none\", NULL, or FALSE.", call. = FALSE)
  }
  key <- tolower(test)
  if (key %in% c("none", "no", "false")) {
    return(NULL)
  }
  if (key %in% c("chisq", "chi-square", "chi square", "lrt")) {
    return("Chisq")
  }
  stop("`test` must be \"Chisq\", \"LRT\", \"none\", NULL, or FALSE.", call. = FALSE)
}

#' Internal gRm fit anova heading helper
#'
#' Supports the api likelihood comparison implementation while preserving its internal contract.
#' @param index One-based internal index.
#' @param fits Internal `fits` value used by this helper.
#' @return The internal `gRm_fit_anova_heading()` computation result.
#' @keywords internal
#' @noRd
gRm_fit_anova_heading <- function(index, fits) {
  label <- gRm_fit_model_label(fits[[index]])
  paste0("Model ", index, ": ", label)
}

#' Internal gRm fit model label helper
#'
#' Supports the api likelihood comparison implementation while preserving its internal contract.
#' @param fit Fitted gRm model.
#' @return The internal `gRm_fit_model_label()` computation result.
#' @keywords internal
#' @noRd
gRm_fit_model_label <- function(fit) {
  spec <- fit$model %||% fit$spec %||% list()
  type <- public_model_type(spec)
  terms <- model_terms_from_spec_for_label(spec)
  if (!nzchar(terms)) {
    return(type)
  }
  paste(type, terms)
}

#' Internal model terms from spec for label helper
#'
#' Supports the api likelihood comparison implementation while preserving its internal contract.
#' @param spec GLLRM model specification.
#' @return The internal `model_terms_from_spec_for_label()` computation result.
#' @keywords internal
#' @noRd
model_terms_from_spec_for_label <- function(spec) {
  ld <- spec$ld %||% data.frame()
  dif <- spec$dif %||% data.frame()
  pieces <- character()
  if (is.data.frame(ld) && nrow(ld)) {
    pieces <- c(pieces, paste0("LD(", paste(paste(ld$item1, ld$item2, sep = ":"), collapse = ", "), ")"))
  }
  if (is.data.frame(dif) && nrow(dif)) {
    pieces <- c(pieces, paste0("DIF(", paste(paste(dif$item, dif$exogenous, sep = ":"), collapse = ", "), ")"))
  }
  paste(pieces, collapse = " ")
}
