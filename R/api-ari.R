#' DIGRAM ARI item score curves
#'
#' Compute the DIGRAM ARI item-by-total-score table for a fitted gRm model.
#'
#' @param fit A fitted gRm model.
#' @param ... Reserved for future extensions and must be empty.
#' @return A data frame with class `gRm_ari`. Rows follow DIGRAM's
#'   item-major, score-minor order. The columns are `ItemNo`, `Item`, `Score`,
#'   `n`, `Obs0..ObsK`, `ObsMean`, `ObsVar`, `Exp0..ExpK`, `ExpMean`,
#'   `ExpVar`, and `z`, where `K` is the global maximum item score.
#' @details
#' ARI is the item response distribution by raw total score produced by the
#' DIGRAM `ITA 18` path. The returned observed and expected category columns are
#' probabilities, not raw counts; the raw row count is stored in `n`. Total
#' score `0` and the maximum possible total score are excluded, following the
#' source call to `Calculate_Ari` with `LeastScore = 1` and
#' `LargestScore = highest_possible_score - 1`.
#'
#' Values are returned at native R precision. The function does not write
#' `Ari_dot.csv` or `Ari_comma.csv`.
#' @aliases print.gRm_ari
#' @export
#' @examples
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0)
#' )
#' fitted <- fit(gllrm(gRm(
#'   data,
#'   items = c("I1", "I2"),
#'   id = "ID",
#'   score_cuts = c(1L, 2L)
#' )))
#' ari_table <- ari(fitted)
#' print(ari_table)
#' head(ari_table)
#' plot(ari_table)
ari <- function(fit, ...) {
  reject_public_dots(...)
  fit <- as_public_gRm_fit(fit)
  new_gRm_ari(
    ari_values(fit),
    metadata = list(
      model_type = fit$model_type %||% fit$spec$model_type %||% NA_character_,
      nitems = nrow(fit$bundle$model$items),
      score_range = c(1L, fit$bundle$model$max_total_score - 1L),
      source_trace = c(
        fit$source_trace %||% character(),
        ari = "source/PAS_scd/DGRirtD.pas::TargetSlut -> source/PAS_skunits/skbias15.pas::Calculate_Ari"
      )
    )
  )
}

#' Internal new gRm ari helper
#'
#' Supports the api ari implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @param metadata Internal `metadata` value used by this helper.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
new_gRm_ari <- function(x, metadata = list()) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  attr(x, "metadata") <- metadata
  class(x) <- c("gRm_ari", "data.frame")
  x
}

#' @export
print.gRm_ari <- function(x, ...) {
  reject_public_dots(...)
  metadata <- attr(x, "metadata") %||% list()
  score_range <- metadata$score_range %||% range(x$Score, na.rm = TRUE)
  if (!length(score_range) || anyNA(score_range) || any(!is.finite(score_range))) {
    score_text <- "none"
  } else {
    score_text <- paste(score_range, collapse = "-")
  }
  cat("gRm: ARI item score curves\n")
  cat("Rows: ", nrow(x), "\n", sep = "")
  cat("Items: ", metadata$nitems %||% length(unique(x$Item)), "\n", sep = "")
  cat("Scores: ", score_text, "\n", sep = "")
  if (!is.null(metadata$model_type) && !is.na(metadata$model_type)) {
    cat("Model: ", metadata$model_type, "\n", sep = "")
  }
  if (nrow(x)) {
    print(utils::head(as.data.frame(x)))
  }
  invisible(x)
}
