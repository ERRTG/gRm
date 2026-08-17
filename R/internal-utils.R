#' Internal %||% helper
#'
#' Supports the internal utils implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @param y Internal `y` value used by this helper.
#' @return The internal `%||%()` computation result.
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
