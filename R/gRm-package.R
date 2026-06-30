#' gRm: source-faithful R implementation of non-GUI DIGRAM slices
#'
#' `gRm` is the maintained R implementation track for selected non-GUI DIGRAM
#' computations. The implementation grows against source-shaped Pascal harnesses
#' in this repository and exposes native numeric R results.
#'
#' The package is currently experimental. The Pascal harness remains the
#' executable source-faithful reference, while `gRm` implements matching
#' computations independently in R and is validated by parity tests.
#'
#' @section Current ordinal-only restriction:
#' The current `gRm` package interface supports ordinal item and exogeneous
#' variables only. Internally, variables are represented with DIGRAM's ordinal
#' source type. Nominal and mixed variable-type behavior from the original
#' DIGRAM source is not implemented in this version. In particular,
#' multi-category nominal exogeneous variables are outside the source-faithful
#' scope of the current package and should not be interpreted as validated
#' nominal-variable analyses.
#'
#' @section Legacy DIGRAM import category-code restriction:
#' [read_digram_project()] supports the simple ordinal DIGRAM import subset
#' where category codes in `DIGRAM.imv` are contiguous one-based codes matching
#' the values in `DIGRAM.csv`. The reader follows `DIGRAM.imp` for the source
#' directory and project prefix, and treats `.imv` recursive-level marker rows
#' as structural metadata rather than variables. Historical DIGRAM projects with
#' zero-based, non-contiguous, or otherwise recoded `.imv` category values are
#' not source-faithfully implemented in this package version. The package does
#' not currently preserve separate source category codes and internal one-based
#' analysis categories.
#'
#' @section Incomplete item fit records:
#' Incomplete item fit records are a documented implementation restriction.
#' DIGRAM's Pascal source can add such records to item fit and item-restscore
#' gamma calculations only when the runtime `NincompleteRecs` state has been
#' populated. gRm does not currently have source-backed evidence for when that
#' runtime state exists, so incomplete item fit records are not synthesized from
#' missing item-response rows.
#'
#' @section Global homogeneity residuals:
#' Global homogeneity residuals and marker cells are a documented
#' implementation restriction. gRm computes the source-backed likelihood-ratio
#' test, p-value, score groups, and item means, but the item-level residual and
#' marker cells remain `NA` because the recovered source does not make the
#' runtime residual variance materialization source-backed. Public summaries
#' report those cells as not source-backed.
#'
#' @section Workflow scope:
#' Public workflows read DIGRAM projects, build item-analysis objects, specify
#' and fit Rasch/GLLRM models, screen candidate terms, run explicit diagnostic
#' functions, and inspect results through [summary()].
#' The installed package does not generate DIGRAM runtime output files or parse
#' historical DIGRAM output. Repository validation against those files belongs
#' under `validation/digram_oracle/`.
#'
#' @section Canonical workflow:
#' A `gRm_analysis` stores the encoded DIGRAM project and score-group setup.
#' From there, either specify a model manually with [gllrm()] or run [screen()]
#' to discover candidate local-dependence and DIF terms. [screen()] returns a
#' screening result, not a fitted model; pass that result to [gllrm()] to create
#' the model selected by screening. [fit()] is the boundary between model
#' specification and model diagnostics.
#'
#' ```
#' analysis <- gRm(data, items = items, exogenous = exogenous)
#'
#' manual_model <- gllrm(analysis, ld = ~ item1:item2, dif = ~ item3:group)
#' screened <- screen(analysis)
#' screened_model <- gllrm(screened)
#'
#' fitted <- fit(screened_model)
#'
#' summary(fitted)
#' item_fit(fitted)
#' local_dependence(fitted)
#' dif(fitted)
#' global_homogeneity(fitted)
#' ```
#'
#' [score_effects()] is analysis-level and can be called directly on a
#' `gRm_analysis` object. Model diagnostics require a fitted `gRm_fit` object.
#'
#' @useDynLib gRm, .registration = TRUE
#' @importFrom stats convolve
#' @importFrom utils type.convert
#' @examples
#' \dontrun{
#' analysis <- gRm(data, items = c("I1", "I2", "I3"), exogenous = "site")
#' fit0 <- fit(gllrm(analysis))
#' fit0
#' summary(fit0)
#' item_fit(fit0)
#' }
#'
"_PACKAGE"

# gRm object shapes
#
# Main `gRm` functions exchange S3 lists with compact public classes such
# as `gRm_analysis`, `gRm_model`, `gRm_fit`, and result classes for
# each diagnostic family. This topic documents the expected internal
# structures for package contributors and tests.
#
# @section Project object:
# Returned by [read_digram_project()]. Contains `paths`, `variables`, `items`,
# `backgrounds`, and `raw_data`. `variables`, `items`, and `backgrounds` are
# data frames with source label code, raw data position, raw maximum category,
# variable name, and item/background flag.
#
# @section Bundle object:
# Returned by [build_item_parameters_bundle()]. Contains `model`, `manifest`,
# and `data`. `model` stores item/background metadata and maximum total score.
# `manifest` stores source-style read and validity counts. `data` stores
# recoded zero-based item scores, background values, total score, status, and
# missing/invalid flags.
#
# @section Fit object:
# Returned by [fit_rasch_base()]. Contains iteration diagnostics, convergence
# state, fitted item gamma parameters, expected item counts, update ratios, and
# observed count summaries.
#
# @section Values object:
# Returned by numeric value helpers. Contains item labels, gamma values,
# thresholds, locations, ICE/MICE effects, item statistics, observed score
# range, estimated parameter count, and convergence diagnostics for public
# accessors.
#
# @name gRm-object-shapes
NULL
