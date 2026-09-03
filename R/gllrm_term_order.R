#' GLLRM source-order canonicalization.
#'
#' Source trace: source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo
#' fills IJinfo by scanning `i := 1..nitems-1`, `j := i+1..nitems`, and
#' fills IXinfo by scanning `i := 1..nitems`, `j := 1..nexogene`.
#' source/digram_source_20260817/skunits/skbias22.pas then reports included LD and DIF terms by iterating
#' k := 1..Nij and k := 1..Nix. R keeps the same model-term order so formula
#' order does not leak into source-shaped model arrays or public term tables.
#' @param analysis Prepared gRm analysis.
#' @param ld Internal `ld` value used by this helper.
#' @param dif Internal `dif` value used by this helper.
#' @return The internal `source_order_model_terms()` computation result.
#' @keywords internal
#' @noRd
source_order_model_terms <- function(analysis, ld, dif) {
  list(
    ld = source_order_ld_table(analysis$items, ld),
    dif = source_order_dif_table(analysis$items, analysis$exogenous, dif)
  )
}

#' Internal source order ld table helper
#'
#' Supports the gllrm term order implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param items Item selection or item metadata.
#' @param terms Internal `terms` value used by this helper.
#' @return The internal `source_order_ld_table()` computation result.
#' @keywords internal
#' @noRd
source_order_ld_table <- function(items, terms) {
  if (is.null(terms) || nrow(terms) == 0L) {
    return(empty_ld_terms())
  }
  item_names <- source_order_names(items)
  out <- terms
  index1 <- match(out$item1, item_names)
  index2 <- match(out$item2, item_names)
  if (anyNA(index1) || anyNA(index2)) {
    stop("LD terms contain item names that are not in the GLLRM context.", call. = FALSE)
  }
  lo <- pmin(index1, index2)
  hi <- pmax(index1, index2)
  out$item1 <- item_names[lo]
  out$item2 <- item_names[hi]
  order_index <- order(lo, hi)
  key <- paste(lo, hi, sep = "\r")
  out <- out[order_index, names(empty_ld_terms()), drop = FALSE]
  key <- key[order_index]
  out <- out[!duplicated(key), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Internal source order dif table helper
#'
#' Supports the gllrm term order implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param items Item selection or item metadata.
#' @param backgrounds Internal `backgrounds` value used by this helper.
#' @param terms Internal `terms` value used by this helper.
#' @return The internal `source_order_dif_table()` computation result.
#' @keywords internal
#' @noRd
source_order_dif_table <- function(items, backgrounds, terms) {
  if (is.null(terms) || nrow(terms) == 0L) {
    return(empty_dif_terms())
  }
  item_names <- source_order_names(items)
  background_names <- source_order_names(backgrounds)
  out <- terms
  item_index <- match(out$item, item_names)
  background_index <- match(out$exogenous, background_names)
  if (anyNA(item_index) || anyNA(background_index)) {
    stop("DIF terms contain names that are not in the GLLRM context.", call. = FALSE)
  }
  order_index <- order(item_index, background_index)
  key <- paste(item_index, background_index, sep = "\r")
  out <- out[order_index, names(empty_dif_terms()), drop = FALSE]
  key <- key[order_index]
  out <- out[!duplicated(key), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Internal source order names helper
#'
#' Supports the gllrm term order implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param x Object or value to process.
#' @return The internal `source_order_names()` computation result.
#' @keywords internal
#' @noRd
source_order_names <- function(x) {
  if (is.null(x)) {
    return(character())
  }
  if (is.data.frame(x)) {
    return(as.character(x$name %||% character()))
  }
  as.character(x)
}
