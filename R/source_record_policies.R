#' Return rows eligible for CML/GLLRM estimation
#'
#' This policy retains the source score-window and complete-background status
#' used by `Estimate_GLLRM`. Diagnostic routines must select their own record
#' policy instead of reusing this set.
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Count_IJtable`.
#' @param context A fitted Rasch/GLLRM context.
#' @return An integer vector of source-row indices in input order.
#' @keywords internal
source_estimation_rows <- function(context) {
  as.integer(context$estimation_rows %||% context$valid_rows %||% integer())
}

#' Return rows with complete, in-range item responses
#'
#' Mirrors `Get_Items` in `source/digram_source_20260817/skunits/SKbias2.pas`: every fitted item
#' must be present and inside its source category range. CM2/CM3's focal
#' `UseItems` set does not narrow this record reader. No estimation score
#' window or exogenous-data condition is applied.
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Count_IJtable`.
#' @param context A fitted Rasch/GLLRM context.
#' @return An integer vector of source-row indices in input order.
#' @keywords internal
source_complete_item_rows <- function(context) {
  if (!is.null(context$complete_item_rows)) {
    return(as.integer(context$complete_item_rows))
  }
  item_matrix <- context$item_matrix
  if (is.null(item_matrix) || nrow(item_matrix) == 0L) {
    return(integer())
  }
  if (ncol(item_matrix) == 0L) {
    return(seq_len(nrow(item_matrix)))
  }

  complete <- rep(TRUE, nrow(item_matrix))
  for (item in seq_len(ncol(item_matrix))) {
    values <- item_matrix[, item]
    complete <- complete & !is.na(values) & values >= 0L &
      values < as.integer(context$item_raw_max[[item]])
  }
  which(complete)
}

#' Return rows with complete, in-range exogenous values
#'
#' Mirrors `get_exogene` in `source/digram_source_20260817/skunits/skbias14.pas`: every
#' exogenous variable must be present and inside its one-based source category
#' range. Projects without exogenous variables accept every source row.
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Count_IJtable`.
#' @param context A fitted Rasch/GLLRM context.
#' @return An integer vector of source-row indices in input order.
#' @keywords internal
source_complete_background_rows <- function(context) {
  if (!is.null(context$complete_background_rows)) {
    return(as.integer(context$complete_background_rows))
  }
  background_matrix <- context$background_matrix
  n_rows <- if (!is.null(context$item_matrix)) nrow(context$item_matrix) else nrow(background_matrix)
  if (is.null(background_matrix) || ncol(background_matrix) == 0L) {
    return(seq_len(n_rows))
  }

  complete <- rep(TRUE, nrow(background_matrix))
  for (background in seq_len(ncol(background_matrix))) {
    values <- background_matrix[, background]
    complete <- complete & !is.na(values) & values >= 1L &
      values <= as.integer(context$background_raw_max[[background]])
  }
  which(complete)
}

#' Return rows accepted by ordinary CM2/CM3 record readers
#'
#' The observed and expected IJ/IX/score and exogenous three-way routines call
#' `Get_Items` and then `get_exogene` without applying the CML score window.
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Count_IJtable`.
#' @param context A fitted Rasch/GLLRM context.
#' @return An integer vector of source-row indices in input order.
#' @keywords internal
source_complete_item_exogenous_rows <- function(context) {
  if (!is.null(context$complete_item_exogenous_rows)) {
    return(as.integer(context$complete_item_exogenous_rows))
  }
  item_rows <- source_complete_item_rows(context)
  background_rows <- source_complete_background_rows(context)
  item_rows[item_rows %in% background_rows]
}

#' Return rows accepted by the observed CM3 item-triple counter
#'
#' `Count_IJK` in `source/digram_source_20260817/skunits/skbias14.pas` calls `get_exogene`, but
#' its failure branch is commented out. Consequently item-complete rows remain
#' eligible even when exogenous data are missing. This compatibility policy is
#' intentionally different from the expected-IJK policy.
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Count_IJK`.
#' @param context A fitted Rasch/GLLRM context.
#' @return An integer vector of source-row indices in input order.
#' @keywords internal
source_cm3_observed_ijk_rows <- function(context) {
  as.integer(context$cm3_observed_ijk_rows %||% source_complete_item_rows(context))
}
