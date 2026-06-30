ari_validate_plot_input <- function(x, confidence = NULL) {
  if (!inherits(x, "gRm_ari")) {
    stop("ARI plotting requires a gRm_ari object returned by ari().", call. = FALSE)
  }
  required <- c("ItemNo", "Item", "Score", "n", "ObsMean", "ExpMean", "ExpVar")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(
      "ARI table is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!nrow(x)) {
    stop("ARI plotting requires a non-empty ARI table.", call. = FALSE)
  }

  numeric_columns <- c("ItemNo", "Score", "n", "ObsMean", "ExpMean", "ExpVar")
  for (column in numeric_columns) {
    value <- x[[column]]
    if (!is.numeric(value) && !is.integer(value)) {
      stop("ARI column ", column, " must be numeric.", call. = FALSE)
    }
    if (any(!is.finite(value))) {
      stop("ARI column ", column, " contains non-finite values.", call. = FALSE)
    }
  }
  if (any(x$n < 0)) {
    stop("ARI row counts cannot be negative.", call. = FALSE)
  }

  if (!is.null(confidence)) {
    ari_validate_probability_scalar(confidence, "confidence")
  }

  invisible(x)
}

ari_validate_whole_number <- function(x, name) {
  if (!is.numeric(x) && !is.integer(x)) {
    stop(name, " must be a positive whole number.", call. = FALSE)
  }
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0 || x != floor(x)) {
    stop(name, " must be a positive whole number.", call. = FALSE)
  }
  as.integer(x)
}

ari_validate_probability_scalar <- function(x, name) {
  if (!is.numeric(x) && !is.integer(x)) {
    stop(name, " must be a number between 0 and 1.", call. = FALSE)
  }
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0 || x >= 1) {
    stop(name, " must be a number between 0 and 1.", call. = FALSE)
  }
  as.numeric(x)
}

ari_score_intervals <- function(x, class_size = 40L) {
  ari_validate_plot_input(x)
  class_size <- ari_validate_whole_number(class_size, "class_size")

  first_item <- min(x$ItemNo)
  score_distribution <- x[x$ItemNo == first_item, c("Score", "n"), drop = FALSE]
  if (anyDuplicated(score_distribution$Score)) {
    stop("ARI table has duplicate Score rows for the first item.", call. = FALSE)
  }
  score_distribution <- score_distribution[order(score_distribution$Score), , drop = FALSE]
  rownames(score_distribution) <- NULL

  intervals <- integer(nrow(score_distribution))
  current_interval <- 1L
  running_frequency <- 0
  for (row in seq_len(nrow(score_distribution))) {
    if (running_frequency >= class_size) {
      current_interval <- current_interval + 1L
      running_frequency <- 0
    }
    running_frequency <- running_frequency + score_distribution$n[[row]]
    intervals[[row]] <- current_interval
  }

  final_interval <- max(intervals)
  if (final_interval > 1L && sum(score_distribution$n[intervals == final_interval]) < class_size) {
    intervals[intervals == final_interval] <- final_interval - 1L
  }

  data.frame(
    Score = as.integer(score_distribution$Score),
    interval = as.integer(intervals),
    stringsAsFactors = FALSE
  )
}

ari_confidence_multiplier <- function(confidence) {
  confidence <- ari_validate_probability_scalar(confidence, "confidence")
  stats::qnorm((1 + confidence) / 2)
}

ari_plot_data <- function(x, class_size = 40L, confidence = 0.95) {
  ari_validate_plot_input(x, confidence = confidence)
  class_size <- ari_validate_whole_number(class_size, "class_size")
  multiplier <- ari_confidence_multiplier(confidence)
  intervals <- ari_score_intervals(x, class_size = class_size)

  table <- as.data.frame(x, stringsAsFactors = FALSE)
  table$interval <- intervals$interval[match(table$Score, intervals$Score)]
  if (anyNA(table$interval)) {
    stop("ARI score intervals could not be matched back to every ARI row.", call. = FALSE)
  }

  keys <- unique(table[c("ItemNo", "Item", "interval")])
  keys <- keys[order(keys$ItemNo, keys$interval), , drop = FALSE]
  rows <- vector("list", nrow(keys))
  for (row in seq_len(nrow(keys))) {
    selected <- table$ItemNo == keys$ItemNo[[row]] &
      table$Item == keys$Item[[row]] &
      table$interval == keys$interval[[row]]
    group <- table[selected, , drop = FALSE]
    n <- sum(group$n)
    if (n <= 0) {
      stop("ARI plot groups must have positive total row counts.", call. = FALSE)
    }
    observed_mean <- sum(group$n * group$ObsMean) / n
    expected_mean <- sum(group$n * group$ExpMean) / n
    expected_variance <- sum(group$n * group$ExpVar) / n
    standard_error <- sqrt(expected_variance / n)

    rows[[row]] <- data.frame(
      ItemNo = as.integer(keys$ItemNo[[row]]),
      Item = as.character(keys$Item[[row]]),
      interval = as.integer(keys$interval[[row]]),
      O = observed_mean,
      N = n,
      E = expected_mean,
      V = expected_variance,
      lower = expected_mean - multiplier * standard_error,
      upper = expected_mean + multiplier * standard_error,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ari_filter_plot_items <- function(x, items = NULL) {
  if (is.null(items)) {
    return(x)
  }
  if (is.list(items) || is.matrix(items) || is.data.frame(items)) {
    stop("items must be item numbers or item names.", call. = FALSE)
  }
  if (is.factor(items)) {
    items <- as.character(items)
  }
  if (!length(items)) {
    stop("items must select at least one item.", call. = FALSE)
  }

  if (is.numeric(items) || is.integer(items)) {
    item_numbers <- vapply(items, ari_validate_whole_number, integer(1L), name = "items")
    unknown <- setdiff(item_numbers, unique(x$ItemNo))
    if (length(unknown)) {
      stop("Unknown items: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    out <- x[x$ItemNo %in% item_numbers, , drop = FALSE]
  } else if (is.character(items)) {
    numeric_like <- grepl("^[0-9]+$", items)
    if (any(numeric_like) && any(!numeric_like)) {
      stop("items must be all item numbers or all item names.", call. = FALSE)
    }
    unknown <- setdiff(items, unique(x$Item))
    if (length(unknown)) {
      stop("Unknown items: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    out <- x[x$Item %in% items, , drop = FALSE]
  } else {
    stop("items must be item numbers or item names.", call. = FALSE)
  }

  out <- out[order(out$ItemNo, out$interval), , drop = FALSE]
  rownames(out) <- NULL
  out
}

ari_validate_facet_layout <- function(plot_data, rows = NULL, columns = NULL) {
  nrow_value <- NULL
  ncol_value <- NULL
  if (!is.null(rows)) {
    nrow_value <- ari_validate_whole_number(rows, "rows")
  }
  if (!is.null(columns)) {
    ncol_value <- ari_validate_whole_number(columns, "columns")
  }
  item_count <- length(unique(plot_data$ItemNo))
  if (!is.null(nrow_value) && !is.null(ncol_value) && nrow_value * ncol_value < item_count) {
    stop(
      "The requested rows and columns layout has fewer panels than selected items.",
      call. = FALSE
    )
  }
  list(nrow = nrow_value, ncol = ncol_value)
}

#' Plot DIGRAM ARI item score curves
#'
#' Draw the DIGRAM ARI item mean curves for a `gRm_ari` table returned by
#' [ari()]. The plot collapses raw total scores into adjacent class intervals,
#' then shows the observed item mean score and a model-expected confidence band
#' for each item.
#'
#' @param x A `gRm_ari` table returned by [ari()]. The table must contain the
#'   item number, item name, total score, row count, observed mean, expected
#'   mean, and expected variance columns produced by [ari()].
#' @param ... Reserved for future extensions. Supplying any argument through
#'   `...` is an error so misspelled plotting arguments are reported instead of
#'   silently ignored.
#' @param class_size Positive whole number giving the minimum target number of
#'   persons in each displayed score interval. Adjacent raw total scores are
#'   collapsed until this target is reached. Defaults to `40`, matching the
#'   reference DIGRAM ARI plotting macro.
#' @param rows Optional positive whole number giving the number of rows in the
#'   faceted plot layout. If both `rows` and `columns` are supplied, their
#'   product must be large enough for the selected items.
#' @param columns Optional positive whole number giving the number of columns in
#'   the faceted plot layout. If both `rows` and `columns` are supplied, their
#'   product must be large enough for the selected items.
#' @param items Optional item selector used to draw only a subset of items. Supply
#'   either item numbers from the `ItemNo` column or item names from the `Item`
#'   column. Numeric and character selectors cannot be mixed, and the source item
#'   order is preserved after filtering.
#' @param confidence Number between 0 and 1 giving the confidence level for the
#'   expected band. The critical value is calculated with
#'   `stats::qnorm((1 + confidence) / 2)` for every confidence level.
#' @param show_expected Single logical value. If `TRUE`, draw the expected mean
#'   as a dashed line in addition to the observed mean and expected confidence
#'   band. The default `FALSE` mirrors the reference DIGRAM/SAS plot.
#' @return A `ggplot` object.
#' @details
#' The observed line is `sum(n * ObsMean) / sum(n)` within each item and score
#' interval. The expected band is centered at `sum(n * ExpMean) / sum(n)` and
#' uses the weighted expected variance `sum(n * ExpVar) / sum(n)`. Item names
#' are shown as ordinary ggplot2 facet-strip titles. No files are read or
#' written.
#' @importFrom rlang .data
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
#' plot(ari_table)
plot.gRm_ari <- function(x,
                         ...,
                         class_size = 40L,
                         rows = NULL,
                         columns = NULL,
                         items = NULL,
                         confidence = 0.95,
                         show_expected = FALSE) {
  dots <- list(...)
  if (length(dots)) {
    stop("The ... argument is reserved for future extensions and must be empty.", call. = FALSE)
  }
  if (!is.logical(show_expected) || length(show_expected) != 1L || is.na(show_expected)) {
    stop("show_expected must be TRUE or FALSE.", call. = FALSE)
  }

  plot_data <- ari_plot_data(x, class_size = class_size, confidence = confidence)
  plot_data <- ari_filter_plot_items(plot_data, items = items)
  layout <- ari_validate_facet_layout(plot_data, rows = rows, columns = columns)
  item_labels <- plot_data[!duplicated(plot_data$ItemNo), c("ItemNo", "Item"), drop = FALSE]
  item_labeller <- ggplot2::as_labeller(stats::setNames(
    as.character(item_labels$Item),
    as.character(item_labels$ItemNo)
  ))

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["interval"]])) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]], fill = "95% CI"),
      alpha = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data[["O"]]),
      linewidth = 0.7,
      colour = "grey40"
    )

  if (show_expected) {
    plot <- plot +
      ggplot2::geom_line(
        ggplot2::aes(y = .data[["E"]]),
        linewidth = 0.55,
        linetype = "dashed",
        colour = "grey25"
      )
  }

  plot +
    ggplot2::facet_wrap(
      ggplot2::vars(.data[["ItemNo"]]),
      nrow = layout$nrow,
      ncol = layout$ncol,
      labeller = item_labeller
    ) +
    ggplot2::scale_x_continuous(breaks = seq_len(max(plot_data$interval))) +
    ggplot2::scale_fill_manual(NULL, values = c("95% CI" = "grey70")) +
    ggplot2::labs(x = "class interval", y = "Mean item score") +
    ggplot2::theme_minimal()
}
