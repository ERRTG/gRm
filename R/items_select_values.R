#' Derive DIGRAM item-selection report values
#'
#' Computes the deterministic data-summary values printed by DIGRAM's
#' `items-select.txt` report. The source path is `DGRexe.pas` item selection,
#' which calls `SKbias2.pas::SHOW_ITEMS`, `Sort_items`, `Calculate_scores`,
#' `Calculate_ChronbachsAlpha`, and `Cut_scores`.
#'
#' Production R computes from `DIGRAM.var` and `DIGRAM.dat`; Pascal and the
#' supplied DIGRAM report are test oracles only.
#'
#' @param project A parsed DIGRAM project from [read_digram_project()].
#' @return A `gRm_items_select_values` object.
#' @examples
#' \dontrun{
#' project <- read_digram_project("path/to/DIGRAM")
#' values <- items_select_values(project)
#' values$score_summary
#' }
#' @keywords internal
items_select_values <- function(project) {
  items <- project$items
  raw <- project$raw_data
  item_matrix <- raw[, items$position, drop = FALSE]
  max_values <- items$raw_max
  valid_item <- sweep(item_matrix, 2L, max_values, `<=`) & item_matrix >= 1L
  complete <- rowSums(valid_item) == ncol(item_matrix)

  scores_zero <- item_matrix - 1L
  complete_scores <- rowSums(scores_zero[complete, , drop = FALSE])
  obtainable_min <- 0L
  obtainable_max <- sum(max_values - 1L)
  score_counts <- tabulate(complete_scores + 1L, nbins = obtainable_max + 1L)
  score_values <- seq.int(0L, obtainable_max)
  n_complete <- sum(complete)

  item_rows <- data.frame(
    item_label = items$label_code,
    item_name = items$name,
    n = integer(nrow(items)),
    mean = numeric(nrow(items)),
    mean_complete = numeric(nrow(items)),
    min_score = integer(nrow(items)),
    max_score = integer(nrow(items)),
    stringsAsFactors = FALSE
  )

  for (item_index in seq_len(nrow(items))) {
    valid <- valid_item[, item_index]
    values <- item_matrix[valid, item_index] - 1L
    complete_values <- item_matrix[complete, item_index] - 1L
    item_rows$n[[item_index]] <- length(values)
    item_rows$mean[[item_index]] <- mean(values)
    item_rows$mean_complete[[item_index]] <- mean(complete_values)
    item_rows$min_score[[item_index]] <- min(complete_values)
    item_rows$max_score[[item_index]] <- max(complete_values)
  }

  score_mean <- sum(score_values * score_counts) / n_complete
  # Source trace: SKbias2.pas::SHOWSCOREDISTRIBUTION computes the displayed
  # variance as the sample variance from the complete-item score distribution.
  score_variance <- (sum(score_values^2 * score_counts) / n_complete - score_mean^2) *
    n_complete / (n_complete - 1L)
  score_sd <- sqrt(score_variance)
  third_moment <- sum((score_values - score_mean)^3 * score_counts / n_complete)
  score_skewness <- third_moment * n_complete^2 /
    ((n_complete - 1L) * (n_complete - 2L)) / score_sd^3

  complete_item_scores <- scores_zero[complete, , drop = FALSE]
  item_variances <- apply(complete_item_scores, 2L, function(x) {
    mean(x^2) - mean(x)^2
  })
  score_var_population <- sum(complete_scores^2) / n_complete - score_mean^2
  # Source trace: SKbias2.pas::Calculate_ChronbachsAlpha uses population item
  # variances and population total-score variance, then multiplies by k/(k-1).
  alpha <- nrow(items) / (nrow(items) - 1L) *
    (1 - sum(item_variances) / score_var_population)

  cut <- item_selection_default_cut(score_counts)
  group1 <- sum(score_counts[seq.int(obtainable_min, cut) + 1L])
  group2 <- n_complete - group1
  score_groups <- data.frame(
    from_score = c(obtainable_min, cut + 1L),
    to_score = c(cut, max(complete_scores)),
    count = c(group1, group2),
    percent = 100 * c(group1, group2) / n_complete,
    cumulative = c(100 * group1 / n_complete, 100),
    stringsAsFactors = FALSE
  )

  distribution <- data.frame(
    score = score_values,
    count = score_counts,
    percent = 100 * score_counts / n_complete,
    cumulative = cumsum(100 * score_counts / n_complete),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      project = project,
      items = item_rows,
      item_labels = paste(items$label_code, collapse = ""),
      score_distribution = distribution,
      score_summary = list(
        n = n_complete,
        missing = nrow(raw) - n_complete,
        obtainable_min = obtainable_min,
        obtainable_max = obtainable_max,
        observed_min = min(complete_scores),
        observed_max = max(complete_scores),
        mean = score_mean,
        variance = score_variance,
        sd = score_sd,
        skewness = score_skewness,
        alpha = alpha
      ),
      score_groups = score_groups
    ),
    class = "gRm_items_select_values"
  )
}

#' Source default score cut for item selection
#'
#' @param score_counts Integer counts indexed by score plus one.
#' @return Integer cutpoint.
#' @keywords internal
item_selection_default_cut <- function(score_counts) {
  # Source trace: DGRexe.pas chooses a median-like cut by accumulating
  # non-extreme score counts and then moving one score down if the previous
  # cumulative count is closer to half of the non-extreme N.
  n <- sum(score_counts)
  cut <- 0L
  previous <- 0L
  current <- 0L
  for (score in seq.int(1L, length(score_counts) - 2L)) {
    previous <- current
    current <- current + score_counts[[score + 1L]]
    if (cut == 0L && current >= 0.5 * n) {
      cut <- score
      sum_before <- previous
      sum_after <- current
    }
  }
  if (cut > 1L) {
    if ((0.5 * n - sum_before) < (sum_after - 0.5 * n)) {
      cut <- cut - 1L
    }
  }
  cut
}
