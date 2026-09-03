# CM2/CM3 source-order margin specifications.

#' Internal cm2 cm3 context items helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `cm2_cm3_context_items()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_context_items <- function(context) {
  items <- context$project$items %||% context$items %||% context$bundle$model$items
  if (is.null(items)) {
    stop("CM2/CM3 context does not contain item metadata.", call. = FALSE)
  }
  items
}

#' Internal cm2 cm3 context backgrounds helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `cm2_cm3_context_backgrounds()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_context_backgrounds <- function(context) {
  backgrounds <- context$project$backgrounds %||%
    context$backgrounds %||%
    context$bundle$model$backgrounds
  if (is.null(backgrounds)) {
    return(data.frame(name = character(), label_code = character(), stringsAsFactors = FALSE))
  }
  backgrounds
}

#' Internal cm2 cm3 variable labels helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param variables Internal `variables` value used by this helper.
#' @return The internal `cm2_cm3_variable_labels()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_variable_labels <- function(variables) {
  labels <- variables$label_code %||% variables$name
  as.character(labels)
}

#' Internal cm2 cm3 spec helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param kind Internal `kind` value used by this helper.
#' @param ... Additional internal arguments passed through this helper.
#' @return The internal `cm2_cm3_spec()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_spec <- function(kind, ...) {
  c(list(kind = kind), list(...))
}

#' Internal cm2 cm3 selected items helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias2.pas::SelectItems`,
#' `source/digram_source_20260817/skunits/SKbias2.pas::Select_items_on_string`,
#' `source/digram_source_20260817/skunits/BIASvars.pas::This_is_an_item`, and
#' `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' @param context Prepared GLLRM computation context.
#' @param items Item selection or item metadata.
#' @return A structured selection record containing the evaluated request,
#'   source-ordered item resolution, and automatic exogenous/score scope.
#' @keywords internal
#' @noRd
cm2_cm3_selected_items <- function(context, items = NULL) {
  # SKbias2.SelectItems turns its blank dialog response into every current item
  # label. NULL is the R representation of that default; empty vectors and
  # empty strings are not treated as a second, ambiguous spelling of the UI
  # default.
  item_metadata <- cm2_cm3_context_items(context)
  item_names <- as.character(item_metadata$name)
  n_items <- length(item_names)

  if (is.null(items)) {
    mode <- "default_all"
    requested_items <- NULL
    selected <- seq_len(n_items)
  } else if (!is.null(dim(items))) {
    stop("`items` must be a one-dimensional vector, not a matrix or array.", call. = FALSE)
  } else if (is.character(items)) {
    mode <- "explicit_names"
    requested_items <- unname(items)
    if (anyNA(items)) {
      stop("`items` must not contain missing values.", call. = FALSE)
    }
    if (!length(items) || any(!nzchar(trimws(items)))) {
      stop("Character `items` must contain nonempty fitted item names.", call. = FALSE)
    }
    selected <- match(items, item_names)
    if (anyNA(selected)) {
      unknown <- unique(items[is.na(selected)])
      background_names <- as.character(cm2_cm3_context_backgrounds(context)$name)
      invalid_scope <- unknown %in% c(background_names, "Score group", "#")
      if (any(invalid_scope)) {
        stop(
          "`items` may select fitted items only; exogenous variables and the score group are automatic.",
          call. = FALSE
        )
      }
      stop(
        "Unknown fitted item name in `items`: ",
        paste(unknown, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  } else if (is.numeric(items) && !is.logical(items)) {
    mode <- "explicit_indices"
    requested_items <- unname(items)
    whole_index <- is.finite(items) & items == floor(items)
    if (anyNA(items) || !all(whole_index)) {
      stop("Numeric `items` must be whole one-based item indices.", call. = FALSE)
    }
    # Validate the original doubles before integer coercion. Otherwise a large
    # finite whole value can overflow to NA_integer_ and escape the intended
    # source item-range error.
    if (!length(items) || any(items < 1 | items > n_items)) {
      stop("Numeric `items` are outside the item index range.", call. = FALSE)
    }
    selected <- as.integer(items)
  } else {
    stop("`items` must be NULL, item names, or one-based item indices.", call. = FALSE)
  }

  if (anyDuplicated(selected)) {
    stop("`items` must not contain duplicates.", call. = FALSE)
  }
  if (length(selected) < 2L) {
    stop("At least two selected items are required for CM2/CM3.", call. = FALSE)
  }

  # Select_items_on_string stores valid labels in a Boolean vector, discarding
  # input order, and CM3_analysis later constructs its heading by scanning item
  # indices from 1..nitems. R rejects duplicates rather than copying the legacy
  # parser's unsafe/permissive behavior, then performs the same source-order
  # normalization. CM3_analysis independently enforces the two-item minimum.
  selected <- sort(selected)
  resolved <- cm2_cm3_selected_item_table(context, selected)
  backgrounds <- cm2_cm3_context_backgrounds(context)
  background_labels <- cm2_cm3_variable_labels(backgrounds)
  exogenous <- data.frame(
    exogenous_index = seq_len(nrow(backgrounds)),
    exogenous_label = background_labels,
    exogenous_name = as.character(backgrounds$name),
    stringsAsFactors = FALSE
  )

  list(
    schema_version = 1L,
    mode = mode,
    requested_items = requested_items,
    resolved_items = resolved,
    selected_count = nrow(resolved),
    model_item_count = n_items,
    exogenous_scope = "all_fitted",
    exogenous = exogenous,
    score_group_included = TRUE,
    score_total_scope = "all_fitted_items"
  )
}

#' Internal cm2 cm3 prepare margins helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param context Prepared GLLRM computation context.
#' @param selected_items Internal `selected_items` value used by this helper.
#' @param include_three_way Whether to include three-way CM3 margins.
#' @return The internal `cm2_cm3_prepare_margins()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_prepare_margins <- function(context, selected_items, include_three_way) {
  # skbias14.Prepare_CM3tests treats UseItems as the focal diagnostic set only:
  # all loops still run in fitted/source item order, every fitted exogenous
  # variable remains automatic, and the score pseudo-variable follows the final
  # exogenous variable. It emits all CM2 rows first, then appends CM3 rows in
  # source item/exogenous/score order.
  backgrounds <- cm2_cm3_context_backgrounds(context)
  n_backgrounds <- nrow(backgrounds)
  selected_items <- sort(as.integer(selected_items))
  included_ld <- cm2_cm3_included_ld_lookup(context)
  included_dif <- cm2_cm3_included_dif_lookup(context)
  specs <- list()

  add_spec <- function(spec) {
    specs[[length(specs) + 1L]] <<- spec
  }

  # Prepare_CM3tests.Use12 requires both endpoints to be selected. An active LD
  # edge suppresses only its corresponding CM2 item-item margin; related
  # three-way margins below remain present.
  if (length(selected_items) >= 2L) {
    for (left_pos in seq_len(length(selected_items) - 1L)) {
      item1 <- selected_items[[left_pos]]
      for (right_pos in seq.int(left_pos + 1L, length(selected_items))) {
        item2 <- selected_items[[right_pos]]
        if (!isTRUE(included_ld[item1, item2])) {
          add_spec(cm2_cm3_spec("item_item", item1 = item1, item2 = item2))
        }
      }
    }
  }

  # Each selected item is crossed with every fitted exogenous variable and then
  # the score group. Active DIF suppresses only the eligible CM2 item-exogenous
  # row; item-score is never suppressed.
  for (item in selected_items) {
    if (n_backgrounds > 0L) {
      for (background in seq_len(n_backgrounds)) {
        if (!isTRUE(included_dif[item, background])) {
          add_spec(cm2_cm3_spec("item_exogenous", item = item, exogenous = background))
        }
      }
    }
    add_spec(cm2_cm3_spec("item_score_group", item = item, score_group = TRUE))
  }

  if (!isTRUE(include_three_way)) {
    return(specs)
  }

  # Prepare_CM3tests.Use123 requires all three focal items. It deliberately does
  # not consult LocalDependence or ItemBias for any CM3-only family.
  if (length(selected_items) >= 3L) {
    for (first_pos in seq_len(length(selected_items) - 2L)) {
      item1 <- selected_items[[first_pos]]
      for (second_pos in seq.int(first_pos + 1L, length(selected_items) - 1L)) {
        item2 <- selected_items[[second_pos]]
        for (third_pos in seq.int(second_pos + 1L, length(selected_items))) {
          item3 <- selected_items[[third_pos]]
          add_spec(cm2_cm3_spec("item_item_item", item1 = item1, item2 = item2, item3 = item3))
        }
      }
    }
  }

  # Every selected item pair is crossed with all fitted exogenous variables and
  # then score, even when the pair or a related item-exogenous margin is a fitted
  # LD/DIF term.
  if (length(selected_items) >= 2L) {
    for (left_pos in seq_len(length(selected_items) - 1L)) {
      item1 <- selected_items[[left_pos]]
      for (right_pos in seq.int(left_pos + 1L, length(selected_items))) {
        item2 <- selected_items[[right_pos]]
        if (n_backgrounds > 0L) {
          for (background in seq_len(n_backgrounds)) {
            add_spec(cm2_cm3_spec(
              "item_item_exogenous",
              item1 = item1,
              item2 = item2,
              exogenous = background
            ))
          }
        }
        add_spec(cm2_cm3_spec(
          "item_item_score_group",
          item1 = item1,
          item2 = item2,
          score_group = TRUE
        ))
      }
    }
  }

  # The final source family crosses each selected item with unordered pairs from
  # (all fitted exogenous variables + score). The loop spelling below preserves
  # Pascal's exogenous-exogenous rows before exogenous-score rows.
  for (item in selected_items) {
    if (n_backgrounds > 0L) {
      for (first_background in seq_len(n_backgrounds)) {
        if (first_background < n_backgrounds) {
          for (second_background in seq.int(first_background + 1L, n_backgrounds)) {
            add_spec(cm2_cm3_spec(
              "item_exogenous_exogenous",
              item = item,
              exogenous1 = first_background,
              exogenous2 = second_background
            ))
          }
        }
        add_spec(cm2_cm3_spec(
          "item_exogenous_score_group",
          item = item,
          exogenous = first_background,
          score_group = TRUE
        ))
      }
    }
  }

  specs
}

#' Internal cm2 cm3 included ld lookup helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `cm2_cm3_included_ld_lookup()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_included_ld_lookup <- function(context) {
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests
  # skips included IJ terms
  # for CM2 item-item margins only.
  n_items <- nrow(cm2_cm3_context_items(context))
  lookup <- matrix(FALSE, nrow = n_items, ncol = n_items)
  ld_specs <- context$ld_specs %||% list()
  for (spec in ld_specs) {
    item1 <- as.integer(spec$item1)
    item2 <- as.integer(spec$item2)
    if (length(item1) == 1L && length(item2) == 1L && !is.na(item1) && !is.na(item2)) {
      lookup[item1, item2] <- TRUE
      lookup[item2, item1] <- TRUE
    }
  }
  lookup
}

#' Internal cm2 cm3 included dif lookup helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param context Prepared GLLRM computation context.
#' @return The internal `cm2_cm3_included_dif_lookup()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_included_dif_lookup <- function(context) {
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::
  # Prepare_CM3tests. DIGRAM 7.04 checks ItemBias(.i1,i2.) while i1 is the
  # candidate item and i2 is the candidate exogenous variable, matching the
  # documented ItemBias(.item,exo.) coordinate order. The preceding
  # ItemBias(.i2,i1.) expression remains only as a commented historical line.
  n_items <- nrow(cm2_cm3_context_items(context))
  n_backgrounds <- nrow(cm2_cm3_context_backgrounds(context))
  lookup <- matrix(FALSE, nrow = n_items, ncol = n_backgrounds)
  dif_specs <- context$dif_specs %||% list()
  for (spec in dif_specs) {
    item <- as.integer(spec$item)
    background <- as.integer(spec$background)
    if (length(item) == 1L && length(background) == 1L && !is.na(item) && !is.na(background)) {
      if (
        item >= 1L && item <= n_items &&
          background >= 1L && background <= n_backgrounds
      ) {
        lookup[item, background] <- TRUE
      }
    }
  }
  lookup
}

#' Internal cm2 cm3 margin label helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @return The internal `cm2_cm3_margin_label()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_margin_label <- function(context, spec) {
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests
  # reports DIGRAM labels
  # for oracle/internal comparison.
  cm2_cm3_margin_string(context, spec, use_labels = TRUE)
}

#' Internal cm2 cm3 margin name helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @return The internal `cm2_cm3_margin_name()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_margin_name <- function(context, spec) {
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests
  # order is retained while
  # public margin names use R-facing variable names.
  cm2_cm3_margin_string(context, spec, use_labels = FALSE)
}

#' Internal cm2 cm3 margin public variables helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @return The internal `cm2_cm3_margin_public_variables()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_margin_public_variables <- function(context, spec) {
  # Source trace: source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests
  # includes the source
  # score-group pseudo-variable; public tables use names, not label codes.
  cm2_cm3_margin_variables(context, spec, use_labels = FALSE)
}

#' Internal cm2 cm3 margin string helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @param use_labels Internal `use_labels` value used by this helper.
#' @return The internal `cm2_cm3_margin_string()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_margin_string <- function(context, spec, use_labels) {
  paste(cm2_cm3_margin_variables(context, spec, use_labels = use_labels), collapse = ":")
}

#' Internal cm2 cm3 margin variables helper
#'
#' Supports the cm2 cm3 specs implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests`.
#' @param context Prepared GLLRM computation context.
#' @param spec GLLRM model specification.
#' @param use_labels Internal `use_labels` value used by this helper.
#' @return The internal `cm2_cm3_margin_variables()` computation result.
#' @keywords internal
#' @noRd
cm2_cm3_margin_variables <- function(context, spec, use_labels) {
  items <- cm2_cm3_context_items(context)
  backgrounds <- cm2_cm3_context_backgrounds(context)
  item_values <- if (use_labels) cm2_cm3_variable_labels(items) else as.character(items$name)
  background_values <- if (use_labels) cm2_cm3_variable_labels(backgrounds) else as.character(backgrounds$name)
  score_group <- "Score group"

  switch(
    spec$kind,
    item_item = item_values[c(spec$item1, spec$item2)],
    item_exogenous = c(item_values[[spec$item]], background_values[[spec$exogenous]]),
    item_score_group = c(item_values[[spec$item]], score_group),
    item_item_item = item_values[c(spec$item1, spec$item2, spec$item3)],
    item_item_exogenous = c(item_values[c(spec$item1, spec$item2)], background_values[[spec$exogenous]]),
    item_item_score_group = c(item_values[c(spec$item1, spec$item2)], score_group),
    item_exogenous_exogenous = c(
      item_values[[spec$item]],
      background_values[c(spec$exogenous1, spec$exogenous2)]
    ),
    item_exogenous_score_group = c(item_values[[spec$item]], background_values[[spec$exogenous]], score_group),
    stop("Unknown CM2/CM3 margin kind.", call. = FALSE)
  )
}
