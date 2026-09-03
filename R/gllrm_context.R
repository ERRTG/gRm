#' GLLRM context construction.
#'
#' Source trace: source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo
#' prepares the model dimensions and source/digram_source_20260817/skunits/skbias12b.pas::
#' Estimate_GLLRM receives the item, LD, DIF, score, and exogenous structures
#' used by the GLLRM fitting loop. The R context keeps those same
#' sufficient statistics in list/data-frame form instead of Pascal global
#' arrays.
#' @param spec GLLRM model specification.
#' @param bundle Source-shaped analysis bundle.
#' @param max_joint_configs Internal `max_joint_configs` value used by this helper.
#' @return A newly assembled internal object or table.
#' @keywords internal
#' @noRd
build_gllrm_context <- function(spec, bundle, max_joint_configs = 200000L) {
  items <- bundle$model$items
  backgrounds <- bundle$model$backgrounds
  data <- bundle$data
  item_matrix <- gllrm_item_matrix(data, items)
  background_matrix <- gllrm_background_matrix(data, backgrounds)
  valid_rows <- which(data$status == 1L)
  item_score_reference <- gllrm_item_score_reference(item_matrix, valid_rows, items)
  ld_specs <- gllrm_ld_specs(spec$ld, items)
  dif_specs <- gllrm_dif_specs(spec$dif, items, backgrounds)

  context <- list(
    spec = spec,
    bundle = bundle,
    items = items,
    backgrounds = backgrounds,
    n_items = nrow(items),
    n_backgrounds = nrow(backgrounds),
    max_total_score = as.integer(bundle$model$max_total_score),
    item_raw_max = as.integer(items$raw_max),
    background_raw_max = as.integer(backgrounds$raw_max %||% integer()),
    item_score_values = lapply(as.integer(items$raw_max), function(raw_max) seq.int(0L, raw_max - 1L)),
    item_score_columns = lapply(as.integer(items$raw_max), seq_len),
    item_score_reference = item_score_reference,
    background_values = lapply(as.integer(backgrounds$raw_max %||% integer()), seq_len),
    item_matrix = item_matrix,
    background_matrix = background_matrix,
    score = as.integer(data$score),
    estimation_rows = valid_rows,
    valid_rows = valid_rows,
    counts = rasch_counts(bundle),
    ld_specs = ld_specs,
    dif_specs = dif_specs,
    dif_background_indices = sort(unique(vapply(dif_specs, `[[`, integer(1L), "background"))),
    max_joint_configs = as.integer(max_joint_configs)
  )
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::Count_IJtable through
  # Count_IXStable use Get_Items/get_exogene over all records, while Count_IJK
  # retains item-complete records after the exogenous failure check was
  # commented out. Keep these diagnostic policies separate from Estimate_GLLRM.
  context$complete_item_rows <- source_complete_item_rows(context)
  context$complete_background_rows <- source_complete_background_rows(context)
  context$complete_item_exogenous_rows <- source_complete_item_exogenous_rows(context)
  context$cm3_observed_ijk_rows <- source_cm3_observed_ijk_rows(context)
  components <- gllrm_ld_components(context)
  context$ld_components_items <- components$items
  context$ld_component_of <- components$component_of
  context$dif_map <- gllrm_dif_map(context)
  context$dif_by_item <- gllrm_dif_by_item(context)
  context$dif_by_item_matrices <- lapply(context$dif_by_item, gllrm_metadata_matrix)
  context$component_keys <- vapply(
    context$ld_components_items,
    gllrm_component_key,
    character(1L)
  )
  context$component_configurations <- stats::setNames(
    lapply(context$ld_components_items, function(component_items) {
      gllrm_build_component_configurations(context, component_items)
    }),
    context$component_keys
  )
  context$component_config_matrices <- stats::setNames(
    lapply(context$component_configurations, function(configs) {
      out <- as.matrix(configs)
      storage.mode(out) <- "integer"
      out
    }),
    context$component_keys
  )
  context$component_config_scores <- stats::setNames(
    lapply(context$component_config_matrices, rowSums),
    context$component_keys
  )
  context$component_lookup <- stats::setNames(
    lapply(context$ld_components_items, function(component_items) {
      lookup <- logical(context$n_items)
      lookup[component_items] <- TRUE
      lookup
    }),
    context$component_keys
  )
  context$component_ld_indices <- stats::setNames(
    lapply(context$ld_components_items, function(component_items) {
      lookup <- logical(context$n_items)
      lookup[component_items] <- TRUE
      which(vapply(context$ld_specs, function(spec) {
        isTRUE(lookup[[spec$item1]]) && isTRUE(lookup[[spec$item2]])
      }, logical(1L)))
    }),
    context$component_keys
  )
  context$component_ld_local_indices <- stats::setNames(
    Map(function(component_items, key) {
      ld_indices <- context$component_ld_indices[[key]]
      if (length(ld_indices) == 0L) {
        return(data.frame(
          ld_index = integer(),
          item1_pos = integer(),
          item2_pos = integer()
        ))
      }
      data.frame(
        ld_index = ld_indices,
        item1_pos = match(vapply(context$ld_specs[ld_indices], `[[`, integer(1L), "item1"), component_items),
        item2_pos = match(vapply(context$ld_specs[ld_indices], `[[`, integer(1L), "item2"), component_items),
        stringsAsFactors = FALSE
      )
    }, context$ld_components_items, context$component_keys),
    context$component_keys
  )
  context$component_ld_local_matrices <- stats::setNames(
    lapply(context$component_ld_local_indices, gllrm_metadata_matrix),
    context$component_keys
  )
  context$observed_ld <- gllrm_observed_ld(context)
  context$observed_dif <- gllrm_observed_dif(context)
  context$score_exo_groups <- gllrm_score_exo_groups(context)
  gllrm_check_component_complexity(context)
  context
}

#' Internal gllrm item score reference helper
#'
#' Supports the gllrm context implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param item_matrix Internal `item_matrix` value used by this helper.
#' @param valid_rows Internal `valid_rows` value used by this helper.
#' @param items Item selection or item metadata.
#' @return The internal `gllrm_item_score_reference()` computation result.
#' @keywords internal
#' @noRd
gllrm_item_score_reference <- function(item_matrix, valid_rows, items) {
  max_score <- max(as.integer(items$raw_max), 0L) - 1L
  if (max_score < 0L || !length(valid_rows)) {
    return(0L)
  }
  scores <- as.integer(item_matrix[valid_rows, , drop = FALSE])
  counts <- tabulate(scores + 1L, nbins = max_score + 1L)
  which.max(counts) - 1L
}

#' Internal gllrm metadata matrix helper
#'
#' Supports the gllrm context implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param x Object or value to process.
#' @return The internal `gllrm_metadata_matrix()` computation result.
#' @keywords internal
#' @noRd
gllrm_metadata_matrix <- function(x) {
  out <- as.matrix(x)
  storage.mode(out) <- "integer"
  dimnames(out) <- NULL
  out
}

#' Internal gllrm ld specs helper
#'
#' Supports the gllrm context implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param ld_terms Internal `ld_terms` value used by this helper.
#' @param items Item selection or item metadata.
#' @return The internal `gllrm_ld_specs()` computation result.
#' @keywords internal
#' @noRd
gllrm_ld_specs <- function(ld_terms, items) {
  if (is.null(ld_terms) || nrow(ld_terms) == 0L) {
    return(list())
  }
  lapply(seq_len(nrow(ld_terms)), function(i) {
    item1 <- match(ld_terms$item1[[i]], items$name)
    item2 <- match(ld_terms$item2[[i]], items$name)
    if (is.na(item1) || is.na(item2)) {
      stop("LD term contains an item that is absent from the source bundle.", call. = FALSE)
    }
    if (item2 < item1) {
      tmp <- item1
      item1 <- item2
      item2 <- tmp
    }
    list(
      term = paste(items$name[[item1]], items$name[[item2]], sep = ":"),
      item1 = item1,
      item2 = item2,
      item1_name = items$name[[item1]],
      item2_name = items$name[[item2]]
    )
  })
}

#' Internal gllrm dif specs helper
#'
#' Supports the gllrm context implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param dif_terms Internal `dif_terms` value used by this helper.
#' @param items Item selection or item metadata.
#' @param backgrounds Internal `backgrounds` value used by this helper.
#' @return The internal `gllrm_dif_specs()` computation result.
#' @keywords internal
#' @noRd
gllrm_dif_specs <- function(dif_terms, items, backgrounds) {
  if (is.null(dif_terms) || nrow(dif_terms) == 0L) {
    return(list())
  }
  lapply(seq_len(nrow(dif_terms)), function(i) {
    item <- match(dif_terms$item[[i]], items$name)
    background <- match(dif_terms$exogenous[[i]], backgrounds$name)
    if (is.na(item) || is.na(background)) {
      stop("DIF term contains a variable that is absent from the source bundle.", call. = FALSE)
    }
    list(
      term = paste(items$name[[item]], backgrounds$name[[background]], sep = ":"),
      item = item,
      background = background,
      item_name = items$name[[item]],
      background_name = backgrounds$name[[background]]
    )
  })
}

#' Internal gllrm item matrix helper
#'
#' Supports the gllrm context implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param data Input data for the computation.
#' @param items Item selection or item metadata.
#' @return The internal `gllrm_item_matrix()` computation result.
#' @keywords internal
#' @noRd
gllrm_item_matrix <- function(data, items) {
  out <- matrix(0L, nrow = nrow(data), ncol = nrow(items))
  for (item_index in seq_len(nrow(items))) {
    out[, item_index] <- as.integer(data[[items$name[[item_index]]]])
  }
  out
}

#' Internal gllrm background matrix helper
#'
#' Supports the gllrm context implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param data Input data for the computation.
#' @param backgrounds Internal `backgrounds` value used by this helper.
#' @return The internal `gllrm_background_matrix()` computation result.
#' @keywords internal
#' @noRd
gllrm_background_matrix <- function(data, backgrounds) {
  out <- matrix(0L, nrow = nrow(data), ncol = nrow(backgrounds))
  for (background_index in seq_len(nrow(backgrounds))) {
    out[, background_index] <- as.integer(data[[backgrounds$name[[background_index]]]])
  }
  out
}

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas and
#' source/digram_source_20260817/skunits/skbias12b.pas store observed IJ margins for each included
#' local-dependence term before calling the iterative GLLRM update. The R code
#' builds the same item-by-item score margins from the input rows.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_observed_ld()` computation result.
#' @keywords internal
#' @noRd
gllrm_observed_ld <- function(context) {
  lapply(context$ld_specs, function(spec) {
    rows <- context$item_raw_max[[spec$item1]]
    cols <- context$item_raw_max[[spec$item2]]
    out <- matrix(
      0L,
      nrow = rows,
      ncol = cols,
      dimnames = list(as.character(seq.int(0L, rows - 1L)), as.character(seq.int(0L, cols - 1L)))
    )
    if (length(context$valid_rows) > 0L) {
      score1 <- context$item_matrix[context$valid_rows, spec$item1]
      score2 <- context$item_matrix[context$valid_rows, spec$item2]
      out[] <- tabulate(score1 + score2 * rows + 1L, nbins = length(out))
    }
    out
  })
}

#' Source trace: source/digram_source_20260817/skunits/skbias22.pas and
#' source/digram_source_20260817/skunits/skbias12b.pas store observed IX margins for each included
#' DIF term. The R code materializes those target-item-by-exogenous margins
#' directly from the complete records used by the current GLLRM.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_observed_dif()` computation result.
#' @keywords internal
#' @noRd
gllrm_observed_dif <- function(context) {
  lapply(context$dif_specs, function(spec) {
    rows <- context$item_raw_max[[spec$item]]
    cols <- context$background_raw_max[[spec$background]]
    out <- matrix(
      0L,
      nrow = rows,
      ncol = cols,
      dimnames = list(as.character(seq.int(0L, rows - 1L)), as.character(seq_len(cols)))
    )
    if (length(context$valid_rows) > 0L) {
      score <- context$item_matrix[context$valid_rows, spec$item]
      value <- context$background_matrix[context$valid_rows, spec$background]
      out[] <- tabulate(score + (value - 1L) * rows + 1L, nbins = length(out))
    }
    out
  })
}

#' Source trace: source/digram_source_20260817/skunits/skbias12b.pas::Estimate_GLLRM iterates over
#' total-score and background combinations when evaluating expected margins. The
#' R helper makes those combinations explicit so the later fit code can perform
#' the same source-shaped summation without Pascal global state.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param context Prepared GLLRM computation context.
#' @param rows Rows used by the computation.
#' @return The internal `gllrm_score_exo_groups()` computation result.
#' @keywords internal
#' @noRd
gllrm_score_exo_groups <- function(context, rows = context$valid_rows) {
  if (length(rows) == 0L) {
    out <- data.frame(score = integer(), count = integer())
    for (name in context$backgrounds$name) {
      out[[name]] <- integer()
    }
    return(out)
  }

  groups <- list()
  index <- new.env(parent = emptyenv())
  for (row in rows) {
    values <- if (context$n_backgrounds > 0L) context$background_matrix[row, ] else integer()
    dif_values <- if (length(context$dif_background_indices) > 0L) {
      values[context$dif_background_indices]
    } else {
      integer()
    }
    key <- paste(c(context$score[[row]], dif_values), collapse = "\r")
    pos <- index[[key]]
    if (is.null(pos)) {
      pos <- length(groups) + 1L
      index[[key]] <- pos
      groups[[pos]] <- c(score = context$score[[row]], count = 0L, values)
    }
    groups[[pos]][["count"]] <- groups[[pos]][["count"]] + 1L
  }

  mat <- do.call(rbind, groups)
  out <- data.frame(score = as.integer(mat[, "score"]), count = as.integer(mat[, "count"]))
  if (context$n_backgrounds > 0L) {
    for (background_index in seq_len(context$n_backgrounds)) {
      out[[context$backgrounds$name[[background_index]]]] <-
        as.integer(mat[, 2L + background_index])
    }
  }
  out
}

#' Implementation guard for explicit GLLRM component enumeration.
#'
#' Source trace: source/digram_source_20260817/skunits/SKTypes.pas fixes the Pascal table and array
#' bounds used by source/digram_source_20260817/skunits/skbias22.pas::LD_Gamma_calculation. R is not
#' constrained by those exact global arrays, but this guard prevents the explicit
#' component-configuration enumeration from silently leaving the intended
#' source-shaped calculation regime.
#' Source trace: `source/digram_source_20260817/skunits/skbias12b.pas::Initialize_GLLRMinfo`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `gllrm_check_component_complexity()` computation result.
#' @keywords internal
#' @noRd
gllrm_check_component_complexity <- function(context) {
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  for (component_items in components) {
    configs <- prod(context$item_raw_max[component_items])
    if (configs > context$max_joint_configs) {
      stop(
        "GLLRM component enumeration exceeds the R implementation guard of ",
        context$max_joint_configs,
        " configurations.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}
