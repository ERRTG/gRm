#' Summary tables for gRm objects
#'
#' Return R-native summary tables for objects produced by the public gRm
#' API.
#'
#' @name gRm-summary
#' @aliases summary.gRm_analysis summary.gRm_item_analysis
#'   summary.gRm_model summary.gRm_gllrm_spec summary.gRm_fit
#'   summary.gRm_gllrm_fit summary.gRm_screen
#'   summary.gRm_score_effects summary.gRm_item_parameters
#'   summary.gRm_item_fit summary.gRm_local_dependence
#'   summary.gRm_dif summary.gRm_global_homogeneity print.summary.gRm
#' @usage
#' \method{summary}{gRm_analysis}(object, which = "data", ...)
#' \method{summary}{gRm_model}(object, which = c("model", "ld", "dif"), ...)
#' \method{summary}{gRm_fit}(object, which = c("fit", "parameters", "terms"), ...)
#' \method{summary}{gRm_screen}(object, which = c("tests", "selected", "all", "score_effects", "bh"), ...)
#' \method{summary}{gRm_score_effects}(object, which = c("selected", "tests", "bh"), ...)
#' \method{summary}{gRm_item_parameters}(object, which = c("tests", "coefficients", "items", "thresholds", "fit"), ...)
#' \method{summary}{gRm_item_fit}(object, which = c("tests", "items", "bh"), ...)
#' \method{summary}{gRm_local_dependence}(object, which = c("selected", "tests", "bh"), ...)
#' \method{summary}{gRm_dif}(object, which = c("selected", "tests", "active", "bh"), ...)
#' \method{summary}{gRm_global_homogeneity}(object, which = c("tests", "summary", "groups", "items"), ...)
#' \method{print}{summary.gRm}(x, ...)
#' @param object An gRm analysis, model, fit, screen, or result object.
#' @param which Named summary section to return. Available values depend on the
#'   object class and are shown in the usage block.
#' @param x A `summary.gRm` object returned by `summary()`.
#' @param ... Reserved for S3 dispatch compatibility; ignored by all gRm
#'   summary and summary-print methods.
#' @return A `summary.gRm` list. The requested tables are available both by
#'   section name and through the `tables` component. `print.summary.gRm()`
#'   returns the input object invisibly.
#' @details
#' The public output layer is based on `summary()` rather than DIGRAM
#' fixed-width text reports. Users request semantic sections with `which`; the
#' returned object contains R data frames suitable for inspection, downstream
#' analysis, or formatting by user code.
#'
#' Available sections are class-specific:
#'
#' * `gRm_analysis`: `"data"`.
#' * `gRm_model`: `"model"`, `"ld"`, `"dif"`.
#' * `gRm_fit`: `"fit"`, `"parameters"`, `"terms"`.
#' * `gRm_screen`: `"tests"`, `"selected"`, `"all"`, `"score_effects"`, `"bh"`.
#' * `gRm_score_effects`: `"selected"`, `"tests"`, `"bh"`.
#' * `gRm_item_parameters`: `"tests"`, `"coefficients"`, `"items"`, `"thresholds"`, `"fit"`.
#' * `gRm_item_fit`: `"tests"` returns compact inferential item-fit
#'   statistics for outfit, infit, and item-restscore gamma; `"items"`
#'   returns the extended per-item outfit/infit diagnostic summary when the
#'   item-fit object was created with `include_extended = TRUE`; `"bh"`
#'   returns Benjamini-Hochberg thresholds for the compact test family.
#' * `gRm_local_dependence`: `"selected"`, `"tests"`, `"bh"`.
#' * `gRm_dif`: `"selected"`, `"tests"`, `"active"`, `"bh"`.
#' * `gRm_global_homogeneity`: `"tests"`, `"summary"`, `"groups"`, `"items"`.
#'
#' For `gRm_item_fit`, the distinction between `"tests"` and `"items"` is
#' intentional. The `"tests"` section is the compact inferential output and is
#' always available. The `"items"` section is the extended per-item diagnostic
#' decomposition and is populated only when the item-fit object was created by
#' `item_fit(..., include_extended = TRUE)`, which is the default. If
#' `include_extended = FALSE`, the `"items"` section is an empty data frame
#' rather than a fallback copy of the test table.
#' @examples
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0),
#'   site = c(0, 0, 1, 1, 0, 1, 0, 1)
#' )
#' analysis <- gRm(data, items = c("I1", "I2"), exogenous = "site", id = "ID")
#' model <- gllrm(analysis)
#' summary(analysis, which = "data")
#' summary(model, which = "model")
#' @seealso [gRm()], [gllrm()], [fit()], [screen()],
#'   [score_effects()], [item_parameters()], [item_fit()],
#'   [local_dependence()], [dif()], [global_homogeneity()]
NULL

#' Print gRm objects
#'
#' Print compact status information for public gRm objects.
#'
#' @name gRm-print
#' @aliases print.gRm_analysis print.gRm_model print.gRm_fit
#'   print.gRm_screen print.gRm_score_effects
#'   print.gRm_item_parameters print.gRm_item_fit
#'   print.gRm_local_dependence print.gRm_dif
#'   print.gRm_global_homogeneity
#' @usage
#' \method{print}{gRm_analysis}(x, ...)
#' \method{print}{gRm_model}(x, ...)
#' \method{print}{gRm_fit}(x, ...)
#' \method{print}{gRm_screen}(x, ...)
#' \method{print}{gRm_score_effects}(x, ...)
#' \method{print}{gRm_item_parameters}(x, ...)
#' \method{print}{gRm_item_fit}(x, ...)
#' \method{print}{gRm_local_dependence}(x, ...)
#' \method{print}{gRm_dif}(x, ...)
#' \method{print}{gRm_global_homogeneity}(x, ...)
#' @param x An gRm analysis, model, fit, screen, or result object.
#' @param ... Reserved for S3 dispatch compatibility; ignored by all gRm
#'   print methods.
#' @return The input object, invisibly.
#' @details
#' The print methods are intentionally compact. They identify the object type
#' and key counts or status fields, but they do not expose the full numerical
#' result surface. Use `summary()` for R-native result tables.
#' @examples
#' data <- data.frame(
#'   ID = 1:4,
#'   I1 = c(0, 1, 0, 1),
#'   I2 = c(1, 0, 1, 0)
#' )
#' analysis <- gRm(data, items = c("I1", "I2"), id = "ID")
#' print(analysis)
#' print(gllrm(analysis))
#' @seealso [summary()], [gRm()], [gllrm()], [fit()]
NULL
