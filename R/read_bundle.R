#' Build the item-parameters analysis bundle
#'
#' Converts a parsed DIGRAM project into the source-shaped bundle used by the
#' source item-parameters report implementation. Item responses are recoded from
#' DIGRAM's one-based raw categories to zero-based scores, row scores and
#' validity flags are computed, and source-style manifest counts are assembled.
#'
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
build_item_parameters_bundle <- function(project) {
  items <- project$items
  backgrounds <- project$backgrounds
  raw <- project$raw_data

  item_data <- data.frame(matrix(nrow = nrow(raw), ncol = nrow(items)))
  names(item_data) <- items$name
  background_data <- data.frame(matrix(nrow = nrow(raw), ncol = nrow(backgrounds)))
  names(background_data) <- backgrounds$name

  score <- integer(nrow(raw))
  status <- integer(nrow(raw))
  missing_items <- integer(nrow(raw))
  invalid_items <- integer(nrow(raw))
  missing_backgrounds <- integer(nrow(raw))
  invalid_backgrounds <- integer(nrow(raw))
  complete_item_flags <- logical(nrow(raw))
  complete_background_flags <- logical(nrow(raw))
  row_scores <- integer(nrow(raw))

  n_complete_items <- 0L
  n_complete_backgrounds <- 0L
  n_valid <- 0L
  n_incomplete <- 0L
  n_useless <- 0L
  complete_item_background_scores <- integer()

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
    }
    if (complete_backgrounds) {
      n_complete_backgrounds <- n_complete_backgrounds + 1L
    }
    complete_item_flags[[row_index]] <- complete_items
    complete_background_flags[[row_index]] <- complete_backgrounds
    row_scores[[row_index]] <- row_score

    if (complete_items && complete_backgrounds) {
      score[[row_index]] <- row_score
      complete_item_background_scores <- c(complete_item_background_scores, row_score)
    } else {
      status[[row_index]] <- 0L
      score[[row_index]] <- -1L
    }
  }

  least_score <- 1L
  # Source trace: CML estimation excludes the lower boundary total score.
  # The upper assumption is the largest complete item+background score observed
  # for the declared data set, matching the BFI 1..118 and example 1..87 reports.
  largest_score <- if (length(complete_item_background_scores) > 0L) {
    max(complete_item_background_scores)
  } else {
    0L
  }
  for (row_index in seq_len(nrow(raw))) {
    if (score[[row_index]] >= least_score && score[[row_index]] <= largest_score) {
      status[[row_index]] <- 1L
      n_valid <- n_valid + 1L
    }
  }
  for (row_index in seq_len(nrow(raw))) {
    if (complete_item_flags[[row_index]]) {
      # Source trace: Count_Margins in source/PAS_skunits/skbias12b.pas
      # checks the complete item score range before reading/checking exogenous
      # values. Boundary complete-item rows therefore do not count as Nuseless.
      if (row_scores[[row_index]] < least_score || row_scores[[row_index]] > largest_score) {
        next
      }
      if (!complete_background_flags[[row_index]]) {
        n_useless <- n_useless + 1L
      }
    } else if (complete_background_flags[[row_index]]) {
      n_incomplete <- n_incomplete + 1L
    } else {
      n_useless <- n_useless + 1L
    }
  }

  data <- cbind(
    item_data,
    background_data,
    data.frame(
      score = score,
      status = status,
      missing_items = missing_items,
      invalid_items = invalid_items,
      missing_backgrounds = missing_backgrounds,
      invalid_backgrounds = invalid_backgrounds
    )
  )

  list(
    model = list(
      scenario = "DIGRAM",
      items = items,
      backgrounds = backgrounds,
      max_total_score = sum(items$raw_max - 1L),
      least_score = least_score,
      largest_score = largest_score
    ),
    manifest = list(
      scenario = "DIGRAM",
      nitems = nrow(items),
      nbackgrounds = nrow(backgrounds),
      nread = nrow(raw),
      ncomplete_items = n_complete_items,
      ncomplete_backgrounds = n_complete_backgrounds,
      nvalid = n_valid,
      nmissing_items = n_incomplete,
      ninvalid_items = 0L,
      nmissing_backgrounds = n_useless,
      ninvalid_backgrounds = 0L
    ),
    data = data
  )
}

#' Write a source-shaped DIGRAM bundle
#'
#' @param bundle Bundle returned by [build_item_parameters_bundle()].
#' @param output_dir Directory to receive `model.tsv`, `manifest.tsv`, and
#'   `GLLRMdata.txt`.
#' @param extra_manifest Named list of additional manifest keys.
#' @return `output_dir`, invisibly.
#' @keywords internal
write_source_bundle <- function(bundle, output_dir, extra_manifest = list()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  items <- bundle$model$items
  backgrounds <- bundle$model$backgrounds
  model_rows <- data.frame(
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
    model_rows <- rbind(
      model_rows,
      data.frame(
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
    )
  }
  for (background_index in seq_len(nrow(backgrounds))) {
    model_rows <- rbind(
      model_rows,
      data.frame(
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
    )
  }

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
