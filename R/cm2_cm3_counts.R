# CM2/CM3 observed table counters.

#' Build a resolved diagnostic score lookup for CM2/CM3 count tables.
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Fitted GLLRM context.
#' @param score_cuts Integer score cut vector already resolved by the caller.
#' @return Integer lookup indexed by zero-based total score plus one, with
#'   `score_cuts` and `score_groups` attributes.
#' @keywords internal
#' @noRd
cm2_cm3_score_group_lookup <- function(context, score_cuts) {
  groups <- global_homogeneity_score_groups(context$bundle, as.integer(score_cuts))
  max_score <- as.integer(context$max_total_score %||% context$bundle$model$max_total_score)
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::Count_IStable,
  # Count_IJStable, and Count_IXStable call Scoregruppe(score, scoredim,
  # Minscore, maxscore, scorecuts). Reuse the existing source-shaped lookup
  # helper so CM2/CM3 score rows keep the same scorecuts boundary semantics.
  lookup <- global_homogeneity_uniform_score_group_lookup(groups, max_score)
  attr(lookup, "score_cuts") <- as.integer(score_cuts)
  attr(lookup, "score_groups") <- groups
  lookup
}

#' Count one observed CM2/CM3 table selected by a prepared margin spec.
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Fitted GLLRM context.
#' @param spec CM2/CM3 margin spec from `cm2_cm3_prepare_margins()`.
#' @param score_group_lookup Resolved lookup from `cm2_cm3_score_group_lookup()`.
#' @return Matrix or array of observed counts.
#' @keywords internal
#' @noRd
cm2_cm3_count_observed <- function(context, spec, score_group_lookup) {
  switch(
    spec$kind,
    item_item = cm2_cm3_count_item_item(context, spec$item1, spec$item2),
    item_exogenous = cm2_cm3_count_item_exogenous(context, spec$item, spec$exogenous),
    item_score_group = cm2_cm3_count_item_score_group(context, spec$item, score_group_lookup),
    item_item_item = cm2_cm3_count_item_item_item(context, spec$item1, spec$item2, spec$item3),
    item_item_exogenous = cm2_cm3_count_item_item_exogenous(
      context,
      spec$item1,
      spec$item2,
      spec$exogenous
    ),
    item_item_score_group = cm2_cm3_count_item_item_score_group(
      context,
      spec$item1,
      spec$item2,
      score_group_lookup
    ),
    item_exogenous_exogenous = cm2_cm3_count_item_exogenous_exogenous(
      context,
      spec$item,
      spec$exogenous1,
      spec$exogenous2
    ),
    item_exogenous_score_group = cm2_cm3_count_item_exogenous_score_group(
      context,
      spec$item,
      spec$exogenous,
      score_group_lookup
    ),
    stop("Unknown CM2/CM3 margin kind.", call. = FALSE)
  )
}

#' Source trace: source/digram_source_20260817/skunits/skbias14.pas::Count_IJtable calls Get_Items
#' and get_exogene over every source record without an estimation-score test.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param item1 Internal `item1` value used by this helper.
#' @param item2 Internal `item2` value used by this helper.
#' @return The internal `cm2_cm3_count_item_item()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_count_item_item <- function(context, item1, item2) {
  dims <- c(cm2_cm3_item_raw_max(context, item1), cm2_cm3_item_raw_max(context, item2))
  out <- matrix(
    0L,
    nrow = dims[[1L]],
    ncol = dims[[2L]],
    dimnames = list(cm2_cm3_item_labels(dims[[1L]]), cm2_cm3_item_labels(dims[[2L]]))
  )
  rows <- source_complete_item_exogenous_rows(context)
  if (!length(rows)) {
    return(out)
  }

  score1 <- context$item_matrix[rows, item1]
  score2 <- context$item_matrix[rows, item2]
  cm2_cm3_fill_observed(out, score1 + 1L, score2 + 1L)
}

#' Source trace: source/digram_source_20260817/skunits/skbias14.pas::Count_IXtable calls Get_Items
#' and get_exogene over every source record without an estimation-score test.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param item One-based item index.
#' @param background One-based exogenous-variable index.
#' @return The internal `cm2_cm3_count_item_exogenous()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_count_item_exogenous <- function(context, item, background) {
  dims <- c(cm2_cm3_item_raw_max(context, item), cm2_cm3_background_raw_max(context, background))
  out <- matrix(
    0L,
    nrow = dims[[1L]],
    ncol = dims[[2L]],
    dimnames = list(cm2_cm3_item_labels(dims[[1L]]), cm2_cm3_background_labels(dims[[2L]]))
  )
  rows <- source_complete_item_exogenous_rows(context)
  if (!length(rows)) {
    return(out)
  }

  item_score <- context$item_matrix[rows, item]
  background_value <- context$background_matrix[rows, background]
  cm2_cm3_fill_observed(out, item_score + 1L, background_value)
}

#' Source trace: source/digram_source_20260817/skunits/skbias14.pas::Count_IStable counts observed
#' item-by-score-group tables using resolved diagnostic score cuts.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param item One-based item index.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @return The internal `cm2_cm3_count_item_score_group()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_count_item_score_group <- function(context, item, score_group_lookup) {
  dims <- c(cm2_cm3_item_raw_max(context, item), cm2_cm3_score_group_count(score_group_lookup))
  out <- matrix(
    0L,
    nrow = dims[[1L]],
    ncol = dims[[2L]],
    dimnames = list(cm2_cm3_item_labels(dims[[1L]]), cm2_cm3_score_group_labels(score_group_lookup))
  )
  rows <- source_complete_item_exogenous_rows(context)
  if (!length(rows)) {
    return(out)
  }

  item_score <- context$item_matrix[rows, item]
  # Count_IStable recomputes the total from Get_Items after record acceptance.
  total_score <- rowSums(context$item_matrix[rows, , drop = FALSE])
  score_group <- cm2_cm3_lookup_scores(score_group_lookup, total_score)
  cm2_cm3_fill_observed(out, item_score + 1L, score_group)
}

#' Source trace: source/digram_source_20260817/skunits/skbias14.pas::Count_IJK retains every row
#' accepted by Get_Items. Its get_exogene failure branch is inside a commented
#' block, so missing exogenous values do not reject an observed item triple.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param item1 Internal `item1` value used by this helper.
#' @param item2 Internal `item2` value used by this helper.
#' @param item3 Internal `item3` value used by this helper.
#' @return The internal `cm2_cm3_count_item_item_item()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_count_item_item_item <- function(context, item1, item2, item3) {
  dims <- c(
    cm2_cm3_item_raw_max(context, item1),
    cm2_cm3_item_raw_max(context, item2),
    cm2_cm3_item_raw_max(context, item3)
  )
  out <- array(0L, dim = dims, dimnames = lapply(dims, cm2_cm3_item_labels))
  rows <- source_cm3_observed_ijk_rows(context)
  if (!length(rows)) {
    return(out)
  }

  score1 <- context$item_matrix[rows, item1]
  score2 <- context$item_matrix[rows, item2]
  score3 <- context$item_matrix[rows, item3]
  cm2_cm3_fill_observed(out, score1 + 1L, score2 + 1L, score3 + 1L)
}

#' Source trace: source/digram_source_20260817/skunits/skbias14.pas::Count_IJXtable counts observed
#' item-by-item-by-exogeneous tables over source-valid records.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param item1 Internal `item1` value used by this helper.
#' @param item2 Internal `item2` value used by this helper.
#' @param background One-based exogenous-variable index.
#' @return The internal `cm2_cm3_count_item_item_exogenous()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_count_item_item_exogenous <- function(context, item1, item2, background) {
  dims <- c(
    cm2_cm3_item_raw_max(context, item1),
    cm2_cm3_item_raw_max(context, item2),
    cm2_cm3_background_raw_max(context, background)
  )
  out <- array(
    0L,
    dim = dims,
    dimnames = list(
      cm2_cm3_item_labels(dims[[1L]]),
      cm2_cm3_item_labels(dims[[2L]]),
      cm2_cm3_background_labels(dims[[3L]])
    )
  )
  rows <- source_complete_item_exogenous_rows(context)
  if (!length(rows)) {
    return(out)
  }

  score1 <- context$item_matrix[rows, item1]
  score2 <- context$item_matrix[rows, item2]
  background_value <- context$background_matrix[rows, background]
  cm2_cm3_fill_observed(out, score1 + 1L, score2 + 1L, background_value)
}

#' Source trace: source/digram_source_20260817/skunits/skbias14.pas::Count_IJStable counts observed
#' item-by-item-by-score-group tables using resolved diagnostic score cuts.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param item1 Internal `item1` value used by this helper.
#' @param item2 Internal `item2` value used by this helper.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @return The internal `cm2_cm3_count_item_item_score_group()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_count_item_item_score_group <- function(context, item1, item2, score_group_lookup) {
  dims <- c(
    cm2_cm3_item_raw_max(context, item1),
    cm2_cm3_item_raw_max(context, item2),
    cm2_cm3_score_group_count(score_group_lookup)
  )
  out <- array(
    0L,
    dim = dims,
    dimnames = list(
      cm2_cm3_item_labels(dims[[1L]]),
      cm2_cm3_item_labels(dims[[2L]]),
      cm2_cm3_score_group_labels(score_group_lookup)
    )
  )
  rows <- source_complete_item_exogenous_rows(context)
  if (!length(rows)) {
    return(out)
  }

  score1 <- context$item_matrix[rows, item1]
  score2 <- context$item_matrix[rows, item2]
  # Count_IJStable applies Scoregruppe only after Get_Items/get_exogene succeed.
  total_score <- rowSums(context$item_matrix[rows, , drop = FALSE])
  score_group <- cm2_cm3_lookup_scores(score_group_lookup, total_score)
  cm2_cm3_fill_observed(out, score1 + 1L, score2 + 1L, score_group)
}

#' Source trace: source/digram_source_20260817/skunits/skbias14.pas::Count_IXZtable counts observed
#' item-by-exogeneous-by-exogeneous tables over source-valid records.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param item One-based item index.
#' @param background1 Internal `background1` value used by this helper.
#' @param background2 Internal `background2` value used by this helper.
#' @return The internal `cm2_cm3_count_item_exogenous_exogenous()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_count_item_exogenous_exogenous <- function(context, item, background1, background2) {
  dims <- c(
    cm2_cm3_item_raw_max(context, item),
    cm2_cm3_background_raw_max(context, background1),
    cm2_cm3_background_raw_max(context, background2)
  )
  out <- array(
    0L,
    dim = dims,
    dimnames = list(
      cm2_cm3_item_labels(dims[[1L]]),
      cm2_cm3_background_labels(dims[[2L]]),
      cm2_cm3_background_labels(dims[[3L]])
    )
  )
  rows <- source_complete_item_exogenous_rows(context)
  if (!length(rows)) {
    return(out)
  }

  item_score <- context$item_matrix[rows, item]
  background_value1 <- context$background_matrix[rows, background1]
  background_value2 <- context$background_matrix[rows, background2]
  cm2_cm3_fill_observed(out, item_score + 1L, background_value1, background_value2)
}

#' Source trace: source/digram_source_20260817/skunits/skbias14.pas::Count_IXStable counts observed
#' item-by-exogeneous-by-score-group tables using resolved diagnostic score cuts.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param item One-based item index.
#' @param background One-based exogenous-variable index.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @return The internal `cm2_cm3_count_item_exogenous_score_group()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_count_item_exogenous_score_group <- function(context, item, background, score_group_lookup) {
  dims <- c(
    cm2_cm3_item_raw_max(context, item),
    cm2_cm3_background_raw_max(context, background),
    cm2_cm3_score_group_count(score_group_lookup)
  )
  out <- array(
    0L,
    dim = dims,
    dimnames = list(
      cm2_cm3_item_labels(dims[[1L]]),
      cm2_cm3_background_labels(dims[[2L]]),
      cm2_cm3_score_group_labels(score_group_lookup)
    )
  )
  rows <- source_complete_item_exogenous_rows(context)
  if (!length(rows)) {
    return(out)
  }

  item_score <- context$item_matrix[rows, item]
  background_value <- context$background_matrix[rows, background]
  # Count_IXStable applies Scoregruppe only after Get_Items/get_exogene succeed.
  total_score <- rowSums(context$item_matrix[rows, , drop = FALSE])
  score_group <- cm2_cm3_lookup_scores(score_group_lookup, total_score)
  cm2_cm3_fill_observed(out, item_score + 1L, background_value, score_group)
}

#' Internal cm2 cm3 fill observed helper
#'
#' Supports the cm2 cm3 counts implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param out Internal `out` value used by this helper.
#' @param index1 Internal `index1` value used by this helper.
#' @param index2 Internal `index2` value used by this helper.
#' @param index3 Internal `index3` value used by this helper.
#' @return The internal `cm2_cm3_fill_observed()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_fill_observed <- function(out, index1, index2, index3 = NULL) {
  dims <- dim(out)
  if (is.null(index3)) {
    valid <- !is.na(index1) & !is.na(index2) &
      index1 >= 1L & index1 <= dims[[1L]] &
      index2 >= 1L & index2 <= dims[[2L]]
    if (any(valid)) {
      linear_index <- index1[valid] + (index2[valid] - 1L) * dims[[1L]]
      out[] <- tabulate(linear_index, nbins = length(out))
    }
    return(out)
  }

  valid <- !is.na(index1) & !is.na(index2) & !is.na(index3) &
    index1 >= 1L & index1 <= dims[[1L]] &
    index2 >= 1L & index2 <= dims[[2L]] &
    index3 >= 1L & index3 <= dims[[3L]]
  if (any(valid)) {
    linear_index <- index1[valid] +
      (index2[valid] - 1L) * dims[[1L]] +
      (index3[valid] - 1L) * dims[[1L]] * dims[[2L]]
    out[] <- tabulate(linear_index, nbins = length(out))
  }
  out
}

#' Internal cm2 cm3 lookup scores helper
#'
#' Supports the cm2 cm3 counts implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @param scores Zero-based score values.
#' @return The internal `cm2_cm3_lookup_scores()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_lookup_scores <- function(score_group_lookup, scores) {
  vapply(
    as.integer(scores),
    function(score) global_homogeneity_lookup_score(score_group_lookup, score),
    integer(1L)
  )
}

#' Internal cm2 cm3 item raw max helper
#'
#' Supports the cm2 cm3 counts implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param item One-based item index.
#' @return The internal `cm2_cm3_item_raw_max()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_item_raw_max <- function(context, item) {
  raw_max <- context$item_raw_max %||% context$items$raw_max %||% context$project$items$raw_max
  as.integer(raw_max[[item]])
}

#' Internal cm2 cm3 background raw max helper
#'
#' Supports the cm2 cm3 counts implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param background One-based exogenous-variable index.
#' @return The internal `cm2_cm3_background_raw_max()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_background_raw_max <- function(context, background) {
  raw_max <- context$background_raw_max %||%
    context$backgrounds$raw_max %||%
    context$project$backgrounds$raw_max
  as.integer(raw_max[[background]])
}

#' Internal cm2 cm3 score group count helper
#'
#' Supports the cm2 cm3 counts implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @return The internal `cm2_cm3_score_group_count()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_score_group_count <- function(score_group_lookup) {
  groups <- attr(score_group_lookup, "score_groups")
  if (!is.null(groups)) {
    return(nrow(groups))
  }
  max(score_group_lookup, na.rm = TRUE)
}

#' Internal cm2 cm3 item labels helper
#'
#' Supports the cm2 cm3 counts implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param raw_max Internal `raw_max` value used by this helper.
#' @return The internal `cm2_cm3_item_labels()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_item_labels <- function(raw_max) {
  as.character(seq.int(0L, raw_max - 1L))
}

#' Internal cm2 cm3 background labels helper
#'
#' Supports the cm2 cm3 counts implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param raw_max Internal `raw_max` value used by this helper.
#' @return The internal `cm2_cm3_background_labels()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_background_labels <- function(raw_max) {
  as.character(seq_len(raw_max))
}

#' Internal cm2 cm3 score group labels helper
#'
#' Supports the cm2 cm3 counts implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param score_group_lookup Internal `score_group_lookup` value used by this helper.
#' @return The internal `cm2_cm3_score_group_labels()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_score_group_labels <- function(score_group_lookup) {
  groups <- attr(score_group_lookup, "score_groups")
  if (!is.null(groups)) {
    return(as.character(groups$label))
  }
  as.character(seq_len(cm2_cm3_score_group_count(score_group_lookup)))
}
