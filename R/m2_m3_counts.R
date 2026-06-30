# M2/M3 observed table counters.

# Build a resolved diagnostic score lookup for M2/M3 count tables.
#
# @param context Fitted GLLRM context.
# @param score_cuts Integer score cut vector already resolved by the caller.
# @return Integer lookup indexed by zero-based total score plus one, with
#   `score_cuts` and `score_groups` attributes.
# @keywords internal
m2_m3_score_group_lookup <- function(context, score_cuts) {
  groups <- global_homogeneity_score_groups(context$bundle, as.integer(score_cuts))
  max_score <- as.integer(context$max_total_score %||% context$bundle$model$max_total_score)
  # Source trace: source/PAS_skunits/skbias14.pas::Count_IStable,
  # Count_IJStable, and Count_IXStable call Scoregruppe(score, scoredim,
  # Minscore, maxscore, scorecuts). Reuse the existing source-shaped lookup
  # helper so CM2/CM3 score rows keep the same scorecuts boundary semantics.
  lookup <- global_homogeneity_uniform_score_group_lookup(groups, max_score)
  attr(lookup, "score_cuts") <- as.integer(score_cuts)
  attr(lookup, "score_groups") <- groups
  lookup
}

# Count one observed M2/M3 table selected by a prepared margin spec.
#
# @param context Fitted GLLRM context.
# @param spec M2/M3 margin spec from [m2_m3_prepare_margins()].
# @param score_group_lookup Resolved lookup from [m2_m3_score_group_lookup()].
# @return Matrix or array of observed counts.
# @keywords internal
m2_m3_count_observed <- function(context, spec, score_group_lookup) {
  switch(
    spec$kind,
    item_item = m2_m3_count_item_item(context, spec$item1, spec$item2),
    item_exogenous = m2_m3_count_item_exogenous(context, spec$item, spec$exogenous),
    item_score_group = m2_m3_count_item_score_group(context, spec$item, score_group_lookup),
    item_item_item = m2_m3_count_item_item_item(context, spec$item1, spec$item2, spec$item3),
    item_item_exogenous = m2_m3_count_item_item_exogenous(
      context,
      spec$item1,
      spec$item2,
      spec$exogenous
    ),
    item_item_score_group = m2_m3_count_item_item_score_group(
      context,
      spec$item1,
      spec$item2,
      score_group_lookup
    ),
    item_exogenous_exogenous = m2_m3_count_item_exogenous_exogenous(
      context,
      spec$item,
      spec$exogenous1,
      spec$exogenous2
    ),
    item_exogenous_score_group = m2_m3_count_item_exogenous_score_group(
      context,
      spec$item,
      spec$exogenous,
      score_group_lookup
    ),
    stop("Unknown M2/M3 margin kind.", call. = FALSE)
  )
}

# Source trace: source/PAS_skunits/skbias14.pas::Count_IJtable counts observed
# item-by-item tables over source-valid records.
m2_m3_count_item_item <- function(context, item1, item2) {
  dims <- c(m2_m3_item_raw_max(context, item1), m2_m3_item_raw_max(context, item2))
  out <- matrix(
    0L,
    nrow = dims[[1L]],
    ncol = dims[[2L]],
    dimnames = list(m2_m3_item_labels(dims[[1L]]), m2_m3_item_labels(dims[[2L]]))
  )
  rows <- context$valid_rows
  if (!length(rows)) {
    return(out)
  }

  score1 <- context$item_matrix[rows, item1]
  score2 <- context$item_matrix[rows, item2]
  m2_m3_fill_observed(out, score1 + 1L, score2 + 1L)
}

# Source trace: source/PAS_skunits/skbias14.pas::Count_IXtable counts observed
# item-by-exogeneous tables over source-valid records.
m2_m3_count_item_exogenous <- function(context, item, background) {
  dims <- c(m2_m3_item_raw_max(context, item), m2_m3_background_raw_max(context, background))
  out <- matrix(
    0L,
    nrow = dims[[1L]],
    ncol = dims[[2L]],
    dimnames = list(m2_m3_item_labels(dims[[1L]]), m2_m3_background_labels(dims[[2L]]))
  )
  rows <- context$valid_rows
  if (!length(rows)) {
    return(out)
  }

  item_score <- context$item_matrix[rows, item]
  background_value <- context$background_matrix[rows, background]
  m2_m3_fill_observed(out, item_score + 1L, background_value)
}

# Source trace: source/PAS_skunits/skbias14.pas::Count_IStable counts observed
# item-by-score-group tables using resolved diagnostic score cuts.
m2_m3_count_item_score_group <- function(context, item, score_group_lookup) {
  dims <- c(m2_m3_item_raw_max(context, item), m2_m3_score_group_count(score_group_lookup))
  out <- matrix(
    0L,
    nrow = dims[[1L]],
    ncol = dims[[2L]],
    dimnames = list(m2_m3_item_labels(dims[[1L]]), m2_m3_score_group_labels(score_group_lookup))
  )
  rows <- context$valid_rows
  if (!length(rows)) {
    return(out)
  }

  item_score <- context$item_matrix[rows, item]
  score_group <- m2_m3_lookup_scores(score_group_lookup, context$score[rows])
  m2_m3_fill_observed(out, item_score + 1L, score_group)
}

# Source trace: source/PAS_skunits/skbias14.pas::Count_IJK counts observed
# item-by-item-by-item tables over source-valid records.
m2_m3_count_item_item_item <- function(context, item1, item2, item3) {
  dims <- c(
    m2_m3_item_raw_max(context, item1),
    m2_m3_item_raw_max(context, item2),
    m2_m3_item_raw_max(context, item3)
  )
  out <- array(0L, dim = dims, dimnames = lapply(dims, m2_m3_item_labels))
  rows <- context$valid_rows
  if (!length(rows)) {
    return(out)
  }

  score1 <- context$item_matrix[rows, item1]
  score2 <- context$item_matrix[rows, item2]
  score3 <- context$item_matrix[rows, item3]
  m2_m3_fill_observed(out, score1 + 1L, score2 + 1L, score3 + 1L)
}

# Source trace: source/PAS_skunits/skbias14.pas::Count_IJXtable counts observed
# item-by-item-by-exogeneous tables over source-valid records.
m2_m3_count_item_item_exogenous <- function(context, item1, item2, background) {
  dims <- c(
    m2_m3_item_raw_max(context, item1),
    m2_m3_item_raw_max(context, item2),
    m2_m3_background_raw_max(context, background)
  )
  out <- array(
    0L,
    dim = dims,
    dimnames = list(
      m2_m3_item_labels(dims[[1L]]),
      m2_m3_item_labels(dims[[2L]]),
      m2_m3_background_labels(dims[[3L]])
    )
  )
  rows <- context$valid_rows
  if (!length(rows)) {
    return(out)
  }

  score1 <- context$item_matrix[rows, item1]
  score2 <- context$item_matrix[rows, item2]
  background_value <- context$background_matrix[rows, background]
  m2_m3_fill_observed(out, score1 + 1L, score2 + 1L, background_value)
}

# Source trace: source/PAS_skunits/skbias14.pas::Count_IJStable counts observed
# item-by-item-by-score-group tables using resolved diagnostic score cuts.
m2_m3_count_item_item_score_group <- function(context, item1, item2, score_group_lookup) {
  dims <- c(
    m2_m3_item_raw_max(context, item1),
    m2_m3_item_raw_max(context, item2),
    m2_m3_score_group_count(score_group_lookup)
  )
  out <- array(
    0L,
    dim = dims,
    dimnames = list(
      m2_m3_item_labels(dims[[1L]]),
      m2_m3_item_labels(dims[[2L]]),
      m2_m3_score_group_labels(score_group_lookup)
    )
  )
  rows <- context$valid_rows
  if (!length(rows)) {
    return(out)
  }

  score1 <- context$item_matrix[rows, item1]
  score2 <- context$item_matrix[rows, item2]
  score_group <- m2_m3_lookup_scores(score_group_lookup, context$score[rows])
  m2_m3_fill_observed(out, score1 + 1L, score2 + 1L, score_group)
}

# Source trace: source/PAS_skunits/skbias14.pas::Count_IXZtable counts observed
# item-by-exogeneous-by-exogeneous tables over source-valid records.
m2_m3_count_item_exogenous_exogenous <- function(context, item, background1, background2) {
  dims <- c(
    m2_m3_item_raw_max(context, item),
    m2_m3_background_raw_max(context, background1),
    m2_m3_background_raw_max(context, background2)
  )
  out <- array(
    0L,
    dim = dims,
    dimnames = list(
      m2_m3_item_labels(dims[[1L]]),
      m2_m3_background_labels(dims[[2L]]),
      m2_m3_background_labels(dims[[3L]])
    )
  )
  rows <- context$valid_rows
  if (!length(rows)) {
    return(out)
  }

  item_score <- context$item_matrix[rows, item]
  background_value1 <- context$background_matrix[rows, background1]
  background_value2 <- context$background_matrix[rows, background2]
  m2_m3_fill_observed(out, item_score + 1L, background_value1, background_value2)
}

# Source trace: source/PAS_skunits/skbias14.pas::Count_IXStable counts observed
# item-by-exogeneous-by-score-group tables using resolved diagnostic score cuts.
m2_m3_count_item_exogenous_score_group <- function(context, item, background, score_group_lookup) {
  dims <- c(
    m2_m3_item_raw_max(context, item),
    m2_m3_background_raw_max(context, background),
    m2_m3_score_group_count(score_group_lookup)
  )
  out <- array(
    0L,
    dim = dims,
    dimnames = list(
      m2_m3_item_labels(dims[[1L]]),
      m2_m3_background_labels(dims[[2L]]),
      m2_m3_score_group_labels(score_group_lookup)
    )
  )
  rows <- context$valid_rows
  if (!length(rows)) {
    return(out)
  }

  item_score <- context$item_matrix[rows, item]
  background_value <- context$background_matrix[rows, background]
  score_group <- m2_m3_lookup_scores(score_group_lookup, context$score[rows])
  m2_m3_fill_observed(out, item_score + 1L, background_value, score_group)
}

m2_m3_fill_observed <- function(out, index1, index2, index3 = NULL) {
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

m2_m3_lookup_scores <- function(score_group_lookup, scores) {
  vapply(
    as.integer(scores),
    function(score) global_homogeneity_lookup_score(score_group_lookup, score),
    integer(1L)
  )
}

m2_m3_item_raw_max <- function(context, item) {
  raw_max <- context$item_raw_max %||% context$items$raw_max %||% context$project$items$raw_max
  as.integer(raw_max[[item]])
}

m2_m3_background_raw_max <- function(context, background) {
  raw_max <- context$background_raw_max %||%
    context$backgrounds$raw_max %||%
    context$project$backgrounds$raw_max
  as.integer(raw_max[[background]])
}

m2_m3_score_group_count <- function(score_group_lookup) {
  groups <- attr(score_group_lookup, "score_groups")
  if (!is.null(groups)) {
    return(nrow(groups))
  }
  max(score_group_lookup, na.rm = TRUE)
}

m2_m3_item_labels <- function(raw_max) {
  as.character(seq.int(0L, raw_max - 1L))
}

m2_m3_background_labels <- function(raw_max) {
  as.character(seq_len(raw_max))
}

m2_m3_score_group_labels <- function(score_group_lookup) {
  groups <- attr(score_group_lookup, "score_groups")
  if (!is.null(groups)) {
    return(as.character(groups$label))
  }
  as.character(seq_len(m2_m3_score_group_count(score_group_lookup)))
}
