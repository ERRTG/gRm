#' gRm: source-faithful R implementation of non-GUI DIGRAM slices
#'
#' `gRm` is the active R implementation track for selected non-GUI DIGRAM
#' computations. The implementation grows against source-shaped Pascal harnesses
#' in this repository and exposes native numeric R results.
#'
#' The package is currently experimental. The Pascal harness remains the
#' executable source-faithful reference, while `gRm` implements matching
#' computations independently in R and is validated by parity tests.
#'
#' Public workflows read DIGRAM projects, build item-analysis objects, specify
#' and fit Rasch/GLLRM models, screen candidate terms, run explicit diagnostic
#' functions, and inspect results through [summary()].
#' The installed package does not generate DIGRAM runtime output files or parse
#' historical DIGRAM output. Repository validation against those files belongs
#' under `validation/digram_oracle/`.
#'
#' @useDynLib gRm, .registration = TRUE
#' @importFrom stats convolve
#' @importFrom utils type.convert
#' @examples
#' \dontrun{
#' analysis <- gRm(data, items = c("I1", "I2", "I3"), exogenous = "site")
#' fit0 <- fit(gllrm(analysis))
#' summary(fit0)
#' summary(item_parameters(fit0), which = "coefficients")
#' summary(item_fit(fit0), which = "tests")
#' }
#'
#' @keywords internal
"_PACKAGE"

#' gRm object shapes
#'
#' Main `gRm` functions exchange S3 lists with compact public classes such
#' as `gRm_analysis`, `gRm_model`, `gRm_fit`, and result classes for
#' each diagnostic family. This topic documents the expected internal
#' structures for package contributors and tests.
#'
#' @section Project object:
#' Returned by [read_digram_project()]. Contains `paths`, `variables`, `items`,
#' `backgrounds`, and `raw_data`. `variables`, `items`, and `backgrounds` are
#' data frames with source label code, raw data position, raw maximum category,
#' variable name, and item/background flag.
#'
#' @section Bundle object:
#' Returned by [build_item_parameters_bundle()]. Contains `model`, `manifest`,
#' and `data`. `model` stores item/background metadata and maximum total score.
#' `manifest` stores source-style read and validity counts. `data` stores
#' recoded zero-based item scores, background values, total score, status, and
#' missing/invalid flags.
#'
#' @section Fit object:
#' Returned by [fit_rasch_base()]. Contains iteration diagnostics, convergence
#' state, fitted item gamma parameters, expected item counts, update ratios, and
#' observed count summaries.
#'
#' @section Values object:
#' Returned by numeric value helpers. Contains item labels, gamma values,
#' thresholds, locations, ICE/MICE effects, item statistics, observed score
#' range, estimated parameter count, and convergence diagnostics for public
#' accessors.
#'
#' @name gRm-object-shapes
#' @keywords internal
NULL
