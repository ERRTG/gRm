# M2/M3 source-order margin specifications.

m2_m3_context_items <- function(context) {
  items <- context$project$items %||% context$items %||% context$bundle$model$items
  if (is.null(items)) {
    stop("M2/M3 context does not contain item metadata.", call. = FALSE)
  }
  items
}

m2_m3_context_backgrounds <- function(context) {
  backgrounds <- context$project$backgrounds %||%
    context$backgrounds %||%
    context$bundle$model$backgrounds
  if (is.null(backgrounds)) {
    return(data.frame(name = character(), label_code = character(), stringsAsFactors = FALSE))
  }
  backgrounds
}

m2_m3_variable_labels <- function(variables) {
  labels <- variables$label_code %||% variables$name
  as.character(labels)
}

m2_m3_spec <- function(kind, ...) {
  c(list(kind = kind), list(...))
}

m2_m3_selected_items <- function(context, items = NULL) {
  # Source trace: source/PAS_skunits/skbias14.pas:5951-6115,
  # source/PAS_skunits/skbias14.pas::Prepare_CM3tests uses a UseItems-style
  # selected set and then loops in source item-index order.
  item_names <- as.character(m2_m3_context_items(context)$name)
  n_items <- length(item_names)

  if (is.null(items)) {
    selected <- seq_len(n_items)
  } else if (is.character(items)) {
    if (anyNA(items)) {
      stop("`items` must not contain missing values.", call. = FALSE)
    }
    selected <- match(items, item_names)
    if (anyNA(selected)) {
      stop("Unknown item name in `items`.", call. = FALSE)
    }
  } else if (is.numeric(items) && !is.logical(items)) {
    whole_index <- is.finite(items) & items == floor(items)
    if (anyNA(items) || !all(whole_index)) {
      stop("Numeric `items` must be whole one-based item indices.", call. = FALSE)
    }
    selected <- as.integer(items)
    if (any(selected < 1L | selected > n_items)) {
      stop("Numeric `items` are outside the item index range.", call. = FALSE)
    }
  } else {
    stop("`items` must be NULL, item names, or one-based item indices.", call. = FALSE)
  }

  if (anyDuplicated(selected)) {
    stop("`items` must not contain duplicates.", call. = FALSE)
  }
  if (length(selected) < 2L) {
    stop("At least two selected items are required for CM2/CM3.", call. = FALSE)
  }

  sort(selected)
}

m2_m3_prepare_margins <- function(context, selected_items, include_three_way) {
  # Source trace: source/PAS_skunits/skbias14.pas:5951-6115,
  # source/PAS_skunits/skbias14.pas::Prepare_CM3tests prepares all CM2 rows
  # first, then appends CM3 rows in source item/exogeneous/score order.
  backgrounds <- m2_m3_context_backgrounds(context)
  n_backgrounds <- nrow(backgrounds)
  selected_items <- sort(as.integer(selected_items))
  included_ld <- m2_m3_included_ld_lookup(context)
  included_dif <- m2_m3_included_dif_lookup(context)
  specs <- list()

  add_spec <- function(spec) {
    specs[[length(specs) + 1L]] <<- spec
  }

  if (length(selected_items) >= 2L) {
    for (left_pos in seq_len(length(selected_items) - 1L)) {
      item1 <- selected_items[[left_pos]]
      for (right_pos in seq.int(left_pos + 1L, length(selected_items))) {
        item2 <- selected_items[[right_pos]]
        if (!isTRUE(included_ld[item1, item2])) {
          add_spec(m2_m3_spec("item_item", item1 = item1, item2 = item2))
        }
      }
    }
  }

  for (item in selected_items) {
    if (n_backgrounds > 0L) {
      for (background in seq_len(n_backgrounds)) {
        if (!isTRUE(included_dif[item, background])) {
          add_spec(m2_m3_spec("item_exogenous", item = item, exogenous = background))
        }
      }
    }
    add_spec(m2_m3_spec("item_score_group", item = item, score_group = TRUE))
  }

  if (!isTRUE(include_three_way)) {
    return(specs)
  }

  if (length(selected_items) >= 3L) {
    for (first_pos in seq_len(length(selected_items) - 2L)) {
      item1 <- selected_items[[first_pos]]
      for (second_pos in seq.int(first_pos + 1L, length(selected_items) - 1L)) {
        item2 <- selected_items[[second_pos]]
        for (third_pos in seq.int(second_pos + 1L, length(selected_items))) {
          item3 <- selected_items[[third_pos]]
          add_spec(m2_m3_spec("item_item_item", item1 = item1, item2 = item2, item3 = item3))
        }
      }
    }
  }

  if (length(selected_items) >= 2L) {
    for (left_pos in seq_len(length(selected_items) - 1L)) {
      item1 <- selected_items[[left_pos]]
      for (right_pos in seq.int(left_pos + 1L, length(selected_items))) {
        item2 <- selected_items[[right_pos]]
        if (n_backgrounds > 0L) {
          for (background in seq_len(n_backgrounds)) {
            add_spec(m2_m3_spec(
              "item_item_exogenous",
              item1 = item1,
              item2 = item2,
              exogenous = background
            ))
          }
        }
        add_spec(m2_m3_spec(
          "item_item_score_group",
          item1 = item1,
          item2 = item2,
          score_group = TRUE
        ))
      }
    }
  }

  for (item in selected_items) {
    if (n_backgrounds > 0L) {
      for (first_background in seq_len(n_backgrounds)) {
        if (first_background < n_backgrounds) {
          for (second_background in seq.int(first_background + 1L, n_backgrounds)) {
            add_spec(m2_m3_spec(
              "item_exogenous_exogenous",
              item = item,
              exogenous1 = first_background,
              exogenous2 = second_background
            ))
          }
        }
        add_spec(m2_m3_spec(
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

m2_m3_included_ld_lookup <- function(context) {
  # Source trace: source/PAS_skunits/skbias14.pas:5951-6115,
  # source/PAS_skunits/skbias14.pas::Prepare_CM3tests skips included IJ terms
  # for CM2 item-item margins only.
  n_items <- nrow(m2_m3_context_items(context))
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

m2_m3_included_dif_lookup <- function(context) {
  # Source trace: source/PAS_skunits/skbias14.pas::Prepare_CM3tests
  # (6003-6020) checks ItemBias(.i2,i1.) while i1 is the candidate item and
  # i2 is the candidate exogenous variable. Other source routines document
  # ItemBias(.item,exo.), so this CM2 preparation branch applies the stored DIF
  # matrix in transposed index order.
  n_items <- nrow(m2_m3_context_items(context))
  n_backgrounds <- nrow(m2_m3_context_backgrounds(context))
  lookup <- matrix(FALSE, nrow = n_items, ncol = n_backgrounds)
  dif_specs <- context$dif_specs %||% list()
  for (spec in dif_specs) {
    item <- as.integer(spec$item)
    background <- as.integer(spec$background)
    if (length(item) == 1L && length(background) == 1L && !is.na(item) && !is.na(background)) {
      transposed_item <- background
      transposed_background <- item
      if (
        transposed_item >= 1L && transposed_item <= n_items &&
          transposed_background >= 1L && transposed_background <= n_backgrounds
      ) {
        lookup[transposed_item, transposed_background] <- TRUE
      }
    }
  }
  lookup
}

m2_m3_margin_label <- function(context, spec) {
  # Source trace: source/PAS_skunits/skbias14.pas:5951-6115,
  # source/PAS_skunits/skbias14.pas::Prepare_CM3tests reports DIGRAM labels
  # for oracle/internal comparison.
  m2_m3_margin_string(context, spec, use_labels = TRUE)
}

m2_m3_margin_name <- function(context, spec) {
  # Source trace: source/PAS_skunits/skbias14.pas:5951-6115,
  # source/PAS_skunits/skbias14.pas::Prepare_CM3tests order is retained while
  # public margin names use R-facing variable names.
  m2_m3_margin_string(context, spec, use_labels = FALSE)
}

m2_m3_margin_public_variables <- function(context, spec) {
  # Source trace: source/PAS_skunits/skbias14.pas:5951-6115,
  # source/PAS_skunits/skbias14.pas::Prepare_CM3tests includes the source
  # score-group pseudo-variable; public tables use names, not label codes.
  m2_m3_margin_variables(context, spec, use_labels = FALSE)
}

m2_m3_margin_string <- function(context, spec, use_labels) {
  paste(m2_m3_margin_variables(context, spec, use_labels = use_labels), collapse = ":")
}

m2_m3_margin_variables <- function(context, spec, use_labels) {
  items <- m2_m3_context_items(context)
  backgrounds <- m2_m3_context_backgrounds(context)
  item_values <- if (use_labels) m2_m3_variable_labels(items) else as.character(items$name)
  background_values <- if (use_labels) m2_m3_variable_labels(backgrounds) else as.character(backgrounds$name)
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
    stop("Unknown M2/M3 margin kind.", call. = FALSE)
  )
}
