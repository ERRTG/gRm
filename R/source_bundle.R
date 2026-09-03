#' Build the source-shaped estimation bundle
#'
#' Converts an internal gRm project into the source-shaped bundle used by base
#' Rasch, item-parameter, item fit, local-dependence, DIF, and global
#' homogeneity calculations. Item responses are recoded from DIGRAM's one-based
#' raw categories to zero-based scores, row scores and validity flags are
#' computed, and source-style manifest counts are assembled.
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Count_Margins`.
#' @param project A project list returned by [read_digram_project()].
#' @return A list with components:
#'   \describe{
#'     \item{`model`}{Scenario metadata, item/background specifications, and the
#'       maximum total score.}
#'     \item{`manifest`}{Source-style read, completeness, validity, and missing
#'       count summaries.}
#'     \item{`data`}{A data frame containing recoded item scores, background
#'       values, total score, status, and missing/invalid flags.}
#'   }
#' @examples
#' \dontrun{
#' project <- read_digram_project("path/to/DIGRAM")
#' bundle <- build_item_parameters_bundle(project)
#' bundle$manifest
#' }
#' @keywords internal
#' @noRd
build_item_parameters_bundle <- function(project) {
  items <- project$items
  backgrounds <- project$backgrounds
  raw <- project$raw_data
  encoded <- source_encode_project_rows(items, backgrounds, raw)
  max_total_score <- sum(items$raw_max - 1L)
  classified <- source_classify_bundle_rows(
    encoded,
    estimation_largest_score = max_total_score - 1L
  )

  data <- cbind(
    encoded$item_data,
    encoded$background_data,
    data.frame(
      score = classified$score,
      status = classified$status,
      missing_items = encoded$missing_items,
      invalid_items = integer(nrow(raw)),
      missing_backgrounds = encoded$missing_backgrounds,
      invalid_backgrounds = integer(nrow(raw))
    )
  )

  list(
    model = list(
      scenario = "DIGRAM",
      items = items,
      backgrounds = backgrounds,
      max_total_score = max_total_score,
      least_score = classified$least_score,
      largest_score = classified$largest_score
    ),
    manifest = list(
      scenario = "DIGRAM",
      nitems = nrow(items),
      nbackgrounds = nrow(backgrounds),
      nread = nrow(raw),
      ncomplete_items = encoded$n_complete_items,
      ncomplete_backgrounds = encoded$n_complete_backgrounds,
      ncomplete_item_backgrounds = encoded$n_complete_item_backgrounds,
      nvalid = classified$n_valid,
      # DIGRAM report label: "missing items"; source Count_Margins counter:
      # Nincomplete, incremented for incomplete item rows.
      nmissing_items = classified$n_incomplete,
      ninvalid_items = 0L,
      # DIGRAM report label: "missing backgrounds"; source Count_Margins
      # counter: Nuseless, incremented for rows unusable for
      # background-conditioned margins, including missing exogeneous values
      # and rows with both missing item and background data.
      nmissing_backgrounds = classified$n_useless,
      ninvalid_backgrounds = 0L
    ),
    data = data
  )
}

#' Encode source project records into item/background state
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Count_Margins`.
#' Mathematical step: traverse records, items, and backgrounds in source order;
#' recode valid item categories to zero-based scores and preserve one-based
#' background values while recording completeness and raw row scores.
#' @param items Parsed item metadata.
#' @param backgrounds Parsed background metadata.
#' @param raw Raw one-based DIGRAM records.
#' @return Encoded data, record flags, score state, and completeness counters.
#' @keywords internal
#' @noRd
source_encode_project_rows <- function(items, backgrounds, raw) {
  item_data <- data.frame(matrix(nrow = nrow(raw), ncol = nrow(items)))
  names(item_data) <- items$name
  background_data <- data.frame(matrix(nrow = nrow(raw), ncol = nrow(backgrounds)))
  names(background_data) <- backgrounds$name
  score <- integer(nrow(raw))
  status <- integer(nrow(raw))
  missing_items <- integer(nrow(raw))
  missing_backgrounds <- integer(nrow(raw))
  complete_item_flags <- logical(nrow(raw))
  complete_background_flags <- logical(nrow(raw))
  row_scores <- integer(nrow(raw))
  complete_item_scores <- integer(nrow(raw))
  n_complete_items <- 0L
  n_complete_backgrounds <- 0L
  n_complete_item_backgrounds <- 0L

  for (row_index in seq_len(nrow(raw))) {
    row_score <- 0L
    complete_items <- TRUE
    complete_backgrounds <- TRUE
    for (item_index in seq_len(nrow(items))) {
      value <- raw[row_index, items$position[[item_index]]]
      if (value < 1L || value > items$raw_max[[item_index]]) {
        item_data[[item_index]][[row_index]] <- -1L
        complete_items <- FALSE
        missing_items[[row_index]] <- 1L
      } else {
        recoded <- value - 1L
        item_data[[item_index]][[row_index]] <- recoded
        row_score <- row_score + recoded
      }
    }
    for (background_index in seq_len(nrow(backgrounds))) {
      value <- raw[row_index, backgrounds$position[[background_index]]]
      if (value < 1L || value > backgrounds$raw_max[[background_index]]) {
        background_data[[background_index]][[row_index]] <- -1L
        complete_backgrounds <- FALSE
        missing_backgrounds[[row_index]] <- 1L
      } else {
        background_data[[background_index]][[row_index]] <- value
      }
    }
    if (complete_items) {
      n_complete_items <- n_complete_items + 1L
      complete_item_scores[[n_complete_items]] <- row_score
    }
    if (complete_backgrounds) {
      n_complete_backgrounds <- n_complete_backgrounds + 1L
    }
    complete_item_flags[[row_index]] <- complete_items
    complete_background_flags[[row_index]] <- complete_backgrounds
    row_scores[[row_index]] <- row_score
    if (complete_items && complete_backgrounds) {
      n_complete_item_backgrounds <- n_complete_item_backgrounds + 1L
      score[[row_index]] <- row_score
    } else {
      status[[row_index]] <- 0L
      score[[row_index]] <- -1L
    }
  }

  list(
    item_data = item_data,
    background_data = background_data,
    score = score,
    status = status,
    missing_items = missing_items,
    missing_backgrounds = missing_backgrounds,
    complete_item_flags = complete_item_flags,
    complete_background_flags = complete_background_flags,
    row_scores = row_scores,
    complete_item_scores = complete_item_scores,
    n_complete_items = n_complete_items,
    n_complete_backgrounds = n_complete_backgrounds,
    n_complete_item_backgrounds = n_complete_item_backgrounds
  )
}

#' Classify encoded records for source estimation and manifest counters
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Count_Margins`.
#' Mathematical step: apply the CML score window before exogenous usability,
#' retaining DIGRAM's distinct `Nincomplete` and `Nuseless` branches.
#' @param encoded Result from `source_encode_project_rows()`.
#' @param estimation_largest_score Upper score bound supplied to the source
#'   GLLRM estimator, normally `highest_possible_score - 1`.
#' @return Status/score vectors, score bounds, and source manifest counters.
#' @keywords internal
#' @noRd
source_classify_bundle_rows <- function(encoded, estimation_largest_score) {
  score <- encoded$score
  status <- encoded$status
  least_score <- 1L
  estimation_largest_score <- as.integer(estimation_largest_score)
  largest_score <- if (encoded$n_complete_items > 0L) {
    max(encoded$complete_item_scores[seq_len(encoded$n_complete_items)])
  } else {
    0L
  }
  n_valid <- 0L
  n_incomplete <- 0L
  n_useless <- 0L
  for (row_index in seq_along(score)) {
    # skbias22.GLLRM_estim supplies the source CML bounds as
    # 1..highest_possible_score-1. The separately retained largest_score is
    # observed-range metadata and must not admit a maximum-score record.
    if (
      score[[row_index]] >= least_score &&
        score[[row_index]] <= estimation_largest_score
    ) {
      status[[row_index]] <- 1L
      n_valid <- n_valid + 1L
    }
  }
  for (row_index in seq_along(score)) {
    if (encoded$complete_item_flags[[row_index]]) {
      # Count_Margins checks the complete item score range before reading the
      # exogenous values, so boundary rows do not increment Nuseless.
      if (
        encoded$row_scores[[row_index]] < least_score ||
          encoded$row_scores[[row_index]] > estimation_largest_score
      ) {
        next
      }
      if (!encoded$complete_background_flags[[row_index]]) {
        n_useless <- n_useless + 1L
      }
    } else if (encoded$complete_background_flags[[row_index]]) {
      n_incomplete <- n_incomplete + 1L
    } else {
      n_useless <- n_useless + 1L
    }
  }
  list(
    score = score,
    status = status,
    least_score = least_score,
    largest_score = largest_score,
    n_valid = n_valid,
    n_incomplete = n_incomplete,
    n_useless = n_useless
  )
}

#' Write a source-shaped DIGRAM bundle
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Count_Margins`.
#' @param bundle Bundle returned by `build_item_parameters_bundle()`.
#' @param output_dir Directory to receive `model.tsv`, `manifest.tsv`, and
#'   `GLLRMdata.txt`.
#' @param extra_manifest Named list of additional manifest keys.
#' @return `output_dir`, invisibly.
#' @keywords internal
#' @noRd
write_source_bundle <- function(bundle, output_dir, extra_manifest = list()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  items <- bundle$model$items
  backgrounds <- bundle$model$backgrounds
  model_row_list <- vector("list", 1L + nrow(items) + nrow(backgrounds))
  model_row_index <- 1L
  model_row_list[[model_row_index]] <- data.frame(
    record_type = "scenario",
    position = 0L,
    name = bundle$model$scenario,
    max_score = bundle$model$max_total_score,
    pair_id = "",
    item1_pos = "",
    item2_pos = "",
    item1 = "",
    item2 = "",
    stringsAsFactors = FALSE
  )
  for (item_index in seq_len(nrow(items))) {
    model_row_index <- model_row_index + 1L
    model_row_list[[model_row_index]] <- data.frame(
      record_type = "item",
      position = item_index,
      name = items$name[[item_index]],
      max_score = items$raw_max[[item_index]] - 1L,
      pair_id = "",
      item1_pos = "",
      item2_pos = "",
      item1 = "",
      item2 = "",
      stringsAsFactors = FALSE
    )
  }
  for (background_index in seq_len(nrow(backgrounds))) {
    model_row_index <- model_row_index + 1L
    model_row_list[[model_row_index]] <- data.frame(
      record_type = "background",
      position = background_index,
      name = backgrounds$name[[background_index]],
      max_score = backgrounds$raw_max[[background_index]],
      pair_id = "",
      item1_pos = "",
      item2_pos = "",
      item1 = "",
      item2 = "",
      stringsAsFactors = FALSE
    )
  }
  model_rows <- do.call(rbind, model_row_list)

  manifest <- bundle$manifest
  manifest$score_min_possible <- 0L
  manifest$score_max_possible <- bundle$model$max_total_score
  manifest$least_score <- bundle$model$least_score
  manifest$largest_score <- bundle$model$largest_score
  manifest$score_min_estimated <- bundle$model$least_score
  manifest$score_max_estimated <- bundle$model$largest_score
  manifest$nld <- 0L
  manifest <- c(manifest, extra_manifest)
  manifest_rows <- data.frame(
    key = names(manifest),
    value = unlist(manifest, use.names = FALSE),
    stringsAsFactors = FALSE
  )

  utils::write.table(
    model_rows,
    file.path(output_dir, "model.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
  utils::write.table(
    manifest_rows,
    file.path(output_dir, "manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
  utils::write.table(
    bundle$data,
    file.path(output_dir, "GLLRMdata.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  invisible(output_dir)
}
