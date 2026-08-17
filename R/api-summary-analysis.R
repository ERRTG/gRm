#' Internal validate summary which helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param which Internal `which` value used by this helper.
#' @param allowed Internal `allowed` value used by this helper.
#' @return The normalized or validated internal value.
#' @keywords internal
#' @noRd
validate_summary_which <- function(which, allowed) {
  if (missing(which) || is.null(which)) {
    return(allowed)
  }
  if (!is.character(which) || anyNA(which) || !length(which)) {
    stop("`which` must name one or more summary sections.", call. = FALSE)
  }
  unknown <- setdiff(which, allowed)
  if (length(unknown)) {
    stop(
      "Unknown `which` summary section",
      if (length(unknown) > 1L) "s" else "",
      ": ",
      paste(unknown, collapse = ", "),
      ". Available sections: ",
      paste(allowed, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  unique(which)
}

#' Internal reject summary which helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param ... Additional internal arguments passed through this helper.
#' @return The internal `reject_summary_which()` computation result.
#' @keywords internal
#' @noRd
reject_summary_which <- function(...) {
  dots <- list(...)
  if ("which" %in% names(dots)) {
    stop("This summary has one public view and does not accept `which`.", call. = FALSE)
  }
  reject_public_dots(...)
  invisible(NULL)
}

#' Internal analysis summary table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @return The internal `analysis_summary_table()` computation result.
#' @keywords internal
#' @noRd
analysis_summary_table <- function(object) {
  data.frame(
    n_rows = nrow(object$project$raw_data),
    n_items = length(object$items),
    n_exogenous = length(object$exogenous),
    items = paste(object$items, collapse = ", "),
    exogenous = paste(object$exogenous, collapse = ", "),
    score_groups = paste(object$score_groups %||% integer(), collapse = ", "),
    stringsAsFactors = FALSE
  )
}

#' Internal analysis summary header helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @return The internal `analysis_summary_header()` computation result.
#' @keywords internal
#' @noRd
analysis_summary_header <- function(object) {
  data_name <- object$data_name %||% object$name %||% "<unnamed>"
  items <- object$items %||% character()
  exogenous <- object$exogenous %||% character()
  item_levels <- analysis_summary_level_counts(object$project$items)
  exogenous_levels <- analysis_summary_level_counts(object$project$backgrounds)
  id <- object$id %||% "none"
  score_groups <- analysis_summary_score_group_label(object)
  c(
    "",
    "Data",
    paste0("  Source: ", data_name),
    paste0("  Rows: ", nrow(object$project$raw_data)),
    paste0("  ID: ", id),
    "",
    "Variables",
    paste0("  Items: ", length(items), " (", summary_header_names(items), ")"),
    paste0("  Item levels: ", summary_header_names(item_levels, empty = "none")),
    paste0(
      "  Exogenous: ",
      length(exogenous),
      " (",
      summary_header_names(exogenous, empty = "none"),
      ")"
    ),
    paste0("  Exogenous levels: ", summary_header_names(exogenous_levels, empty = "none")),
    "",
    "Score groups",
    paste0("  Groups: ", score_groups)
  )
}

#' Internal analysis summary score group label helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @return The internal `analysis_summary_score_group_label()` computation result.
#' @keywords internal
#' @noRd
analysis_summary_score_group_label <- function(object) {
  distribution <- tryCatch(
    analysis_score_group_distribution(object),
    error = function(e) data.frame()
  )
  if (is.data.frame(distribution) && nrow(distribution)) {
    groups <- distribution[distribution$score != "Total", , drop = FALSE]
    if (nrow(groups)) {
      return(paste0(nrow(groups), " (", paste(groups$score, collapse = ", "), ")"))
    }
  }
  summary_header_names(object$score_groups %||% integer(), empty = "none")
}

#' Internal analysis score group distribution helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @return The internal `analysis_score_group_distribution()` computation result.
#' @keywords internal
#' @noRd
analysis_score_group_distribution <- function(object) {
  values <- items_select_values(object$project)
  score_summary <- values$score_summary
  score_distribution <- values$score_distribution
  cuts <- as.integer(object$score_groups %||% integer())
  if (!length(cuts)) {
    return(data.frame())
  }

  n_complete <- as.integer(score_summary$n)
  rows <- list()
  from_score <- 0L
  cumulative_cases <- 0L
  for (group_index in seq_along(cuts)) {
    to_score <- cuts[[group_index]]
    if (group_index == length(cuts)) {
      to_score <- min(to_score, score_summary$observed_max)
    }
    if (from_score <= to_score) {
      in_group <- score_distribution$score >= from_score & score_distribution$score <= to_score
      cases <- as.integer(sum(score_distribution$count[in_group]))
      cumulative_cases <- cumulative_cases + cases
      rows[[length(rows) + 1L]] <- data.frame(
        score = score_range_label(from_score, to_score),
        count = cases,
        percent = 100 * cases / n_complete,
        cumulative = 100 * cumulative_cases / n_complete,
        stringsAsFactors = FALSE
      )
    }
    from_score <- cuts[[group_index]] + 1L
  }

  rows[[length(rows) + 1L]] <- data.frame(
    score = "Total",
    count = n_complete,
    percent = 100,
    cumulative = 100,
    stringsAsFactors = FALSE
  )
  out <- do.call(rbind, rows)
  attr(out, "observed_min") <- as.integer(score_summary$observed_min)
  attr(out, "observed_max") <- as.integer(score_summary$observed_max)
  attr(out, "missing_item_score_rows") <- as.integer(score_summary$missing)
  out
}

#' Internal score range label helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param from_score Internal `from_score` value used by this helper.
#' @param to_score Internal `to_score` value used by this helper.
#' @return The internal `score_range_label()` computation result.
#' @keywords internal
#' @noRd
score_range_label <- function(from_score, to_score) {
  if (identical(as.integer(from_score), as.integer(to_score))) {
    return(as.character(as.integer(from_score)))
  }
  paste0(as.integer(from_score), "-", as.integer(to_score))
}

#' Internal print analysis score group distribution helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param table Numeric contingency or result table.
#' @return The input object, invisibly, or a graphics result.
#' @keywords internal
#' @noRd
print_analysis_score_group_distribution <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(invisible(NULL))
  }
  observed_min <- attr(table, "observed_min", exact = TRUE)
  observed_max <- attr(table, "observed_max", exact = TRUE)
  missing <- attr(table, "missing_item_score_rows", exact = TRUE)

  cat("\nScore group distribution\n")
  cat("Observed score range: ", observed_min, "-", observed_max, "\n", sep = "")
  cat("Score-group cases: ", sum(table$count[table$score != "Total"]), "\n", sep = "")
  cat("Missing item-score rows: ", missing, "\n\n", sep = "")

  display <- table
  display$percent <- sprintf("%.1f", display$percent)
  display$cumulative <- sprintf("%.1f", display$cumulative)
  names(display) <- c("Score", "Count", "Percent", "Cumulative")
  print(display, row.names = FALSE)
  invisible(NULL)
}

#' Internal analysis summary level counts helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param variables Internal `variables` value used by this helper.
#' @return The internal `analysis_summary_level_counts()` computation result.
#' @keywords internal
#' @noRd
analysis_summary_level_counts <- function(variables) {
  if (!is.data.frame(variables) || !nrow(variables)) {
    return(character())
  }
  paste0(variables$name, "=", variables$raw_max)
}

#' Internal summary header names helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param x Object or value to process.
#' @param empty Internal `empty` value used by this helper.
#' @return The internal `summary_header_names()` computation result.
#' @keywords internal
#' @noRd
summary_header_names <- function(x, empty = "none") {
  if (!length(x)) {
    return(empty)
  }
  paste(as.character(x), collapse = ", ")
}

#' Internal model summary table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @return The internal `model_summary_table()` computation result.
#' @keywords internal
#' @noRd
model_summary_table <- function(object) {
  data.frame(
    model_type = public_model_type(object),
    n_items = length(object$analysis$items),
    n_exogenous = length(object$analysis$exogenous),
    n_ld = nrow(object$ld %||% data.frame()),
    n_dif = nrow(object$dif %||% data.frame()),
    stringsAsFactors = FALSE
  )
}

#' Internal fit summary table helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @return The internal `fit_summary_table()` computation result.
#' @keywords internal
#' @noRd
fit_summary_table <- function(object) {
  values <- object$values %||% list()
  data.frame(
    model_type = public_model_type(object$spec),
    log_likelihood = values$log_likelihood %||% NA_real_,
    n_parameters = values$n_parameters %||% NA_integer_,
    likelihood_n = values$likelihood_n %||% NA_integer_,
    converged = object$convergence$converged %||% NA,
    iterations = object$convergence$iterations %||% values$n_step %||% NA_integer_,
    delta = object$convergence$report_delta %||% object$convergence$delta %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Internal public model type helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @return The internal `public_model_type()` computation result.
#' @keywords internal
#' @noRd
public_model_type <- function(object) {
  object$model_type %||% if (
    nrow(object$ld %||% data.frame()) > 0L ||
      nrow(object$dif %||% data.frame()) > 0L
  ) {
    "gllrm"
  } else {
    "rasch"
  }
}

#' Internal public model label helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @param noun Internal `noun` value used by this helper.
#' @return The internal `public_model_label()` computation result.
#' @keywords internal
#' @noRd
public_model_label <- function(object, noun) {
  paste("gRm:", public_model_type_label(public_model_type(object)), noun)
}

#' Internal public model type label helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param model_type Internal `model_type` value used by this helper.
#' @return The internal `public_model_type_label()` computation result.
#' @keywords internal
#' @noRd
public_model_type_label <- function(model_type) {
  if (identical(model_type, "rasch")) {
    "Rasch model"
  } else {
    "Graphical Log-Linear Rasch Model"
  }
}

#' Internal public screen terms helper
#'
#' Supports the api summary implementation while preserving its internal contract.
#' @param object Object dispatched to this helper.
#' @return The internal `public_screen_terms()` computation result.
#' @keywords internal
#' @noRd
public_screen_terms <- function(object) {
  if (!is.null(object$terms)) {
    return(object$terms)
  }
  model_terms(object)
}
