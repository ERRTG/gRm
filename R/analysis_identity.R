#' Create a deterministic fingerprint for a canonical R object
#'
#' The payload is serialized with R's version-2 format and hashed with MD5.
#' Equality-sensitive public checks also compare the canonical payload itself,
#' so the compact hash is an identifier rather than the sole collision guard.
#'
#' @param x Canonical object without calls, environments, or transient paths.
#' @param kind Fingerprint namespace.
#' @return A scalar character fingerprint.
#' @keywords internal
gRm_stable_fingerprint <- function(x, kind) {
  path <- tempfile(pattern = "gRm-fingerprint-", fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  writeBin(serialize(x, NULL, version = 2L), path)
  paste0(kind, ":", unname(tools::md5sum(path)))
}

#' Canonical category metadata for analysis identity
#'
#' @param project Internal `gRm_project` object.
#' @return A list containing ordered variable declarations and category maps.
#' @keywords internal
gRm_analysis_category_identity <- function(project) {
  metadata_fields <- intersect(
    c("label_code", "position", "raw_max", "vtype", "name", "is_item"),
    names(project$variables)
  )
  variables <- project$variables[, metadata_fields, drop = FALSE]
  rownames(variables) <- NULL

  category_levels <- project$category_levels %||% list(
    items = stats::setNames(
      lapply(project$items$raw_max, seq_len),
      project$items$name
    ),
    backgrounds = stats::setNames(
      lapply(project$backgrounds$raw_max, seq_len),
      project$backgrounds$name
    )
  )
  list(variables = variables, category_levels = category_levels)
}

#' Build the canonical analysis identity payload
#'
#' The payload includes the selected encoded values in row order, explicit row
#' identity, source missingness/validity state, category mappings, variable
#' roles, and score cuts. Those are precisely the inputs that can change the
#' conditional likelihood sample or interpretation while leaving column names
#' unchanged.
#'
#' @param project Internal `gRm_project` object.
#' @param data User-facing source data retained by the analysis.
#' @param id Optional identifier column name.
#' @param score_groups Normalized analysis score cuts.
#' @param bundle Optional already-built source estimation bundle.
#' @return Canonical analysis-identity list.
#' @keywords internal
gRm_analysis_identity <- function(project,
                                  data,
                                  id,
                                  score_groups,
                                  bundle = NULL) {
  bundle <- bundle %||% build_item_parameters_bundle(project)
  encoded <- unclass(as.matrix(project$raw_data))
  storage.mode(encoded) <- "integer"
  dimnames(encoded) <- NULL

  state_fields <- intersect(
    c(
      "score", "status", "missing_items", "invalid_items",
      "missing_backgrounds", "invalid_backgrounds"
    ),
    names(bundle$data)
  )
  row_state <- bundle$data[, state_fields, drop = FALSE]
  rownames(row_state) <- NULL

  row_identity <- if (
    !is.null(id) && is.data.frame(data) && id %in% names(data)
  ) {
    data[[id]]
  } else {
    seq_len(nrow(encoded))
  }

  list(
    schema = "gRm-analysis-identity-v1",
    selected = list(
      items = as.character(project$items$name),
      exogenous = as.character(project$backgrounds$name),
      id = id
    ),
    score_cuts = as.integer(score_groups),
    categories = gRm_analysis_category_identity(project),
    row_position = seq_len(nrow(encoded)),
    row_names = if (is.data.frame(data)) rownames(data) else as.character(seq_len(nrow(encoded))),
    row_identity = row_identity,
    encoded_values = encoded,
    row_state = row_state
  )
}

#' Build exact likelihood-sample identity
#'
#' @param analysis_fingerprint Scalar analysis fingerprint.
#' @param bundle Source-shaped estimation bundle.
#' @return Structured row mask, row indices, count, and fingerprint.
#' @keywords internal
gRm_likelihood_sample_identity <- function(analysis_fingerprint, bundle) {
  row_mask <- unname(as.logical(bundle$data$status == 1L))
  payload <- list(
    schema = "gRm-likelihood-sample-v1",
    analysis_fingerprint = analysis_fingerprint,
    row_mask = row_mask
  )
  list(
    schema = payload$schema,
    fingerprint = gRm_stable_fingerprint(payload, "gRm-likelihood-sample-v1"),
    row_mask = row_mask,
    row_indices = which(row_mask),
    n = as.integer(sum(row_mask))
  )
}

#' Materialize analysis and likelihood identity fields
#'
#' @param project Internal `gRm_project` object.
#' @param data User-facing source data retained by the analysis.
#' @param id Optional identifier column name.
#' @param score_groups Normalized analysis score cuts.
#' @param bundle Optional already-built source estimation bundle.
#' @return List containing canonical identity, its fingerprint, and likelihood
#'   sample identity.
#' @keywords internal
gRm_analysis_identity_fields <- function(project,
                                         data,
                                         id,
                                         score_groups,
                                         bundle = NULL) {
  bundle <- bundle %||% build_item_parameters_bundle(project)
  identity <- gRm_analysis_identity(
    project,
    data,
    id,
    score_groups,
    bundle = bundle
  )
  fingerprint <- gRm_stable_fingerprint(identity, "gRm-analysis-identity-v1")
  list(
    identity = identity,
    fingerprint = fingerprint,
    likelihood_sample = gRm_likelihood_sample_identity(fingerprint, bundle)
  )
}

#' Verify that an analysis identity still matches its project
#'
#' @param analysis A `gRm_analysis` object.
#' @param bundle Source-shaped estimation bundle built for fitting.
#' @return Current likelihood-sample identity, invisibly.
#' @keywords internal
assert_gRm_analysis_identity <- function(analysis, bundle) {
  current <- gRm_analysis_identity_fields(
    analysis$project,
    analysis$data,
    analysis$id,
    analysis$score_groups,
    bundle = bundle
  )
  if (
    !identical(analysis$analysis_fingerprint, current$fingerprint) ||
      !identical(analysis$analysis_identity, current$identity)
  ) {
    stop(
      "The analysis data or encoding metadata changed after construction; rebuild the gRm analysis before fitting.",
      call. = FALSE
    )
  }
  if (!identical(analysis$likelihood_sample, current$likelihood_sample)) {
    stop(
      "The analysis likelihood row mask changed after construction; rebuild the gRm analysis before fitting.",
      call. = FALSE
    )
  }
  invisible(current$likelihood_sample)
}
