#' Draw one source-cumulative weighted index
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param weights Non-negative unnormalized weights.
#' @param rng Private bootstrap RNG.
#' @param one_minus Whether to reproduce `x := 1 - random`; the bounded joint
#'   validation harness consumes its Park--Miller variate directly.
#' @return One-based selected index.
#' @keywords internal
m2_m3_bootstrap_weighted_index <- function(weights, rng, one_minus = TRUE) {
  weights <- as.numeric(weights)
  total <- sum(weights)
  if (!is.finite(total) || total <= 0 || any(!is.finite(weights)) || any(weights < 0)) {
    stop("The fitted GLLRM has no positive bootstrap probability for an observed score stratum.", call. = FALSE)
  }
  # SKbias8 uses x := 1 - random so the variate is never exactly zero.
  uniform <- rng$uniform()
  target <- if (isTRUE(one_minus)) (1 - uniform) * total else uniform * total
  selected <- which(target <= cumsum(weights))
  if (length(selected)) selected[[1L]] else length(weights)
}

#' Generate response patterns for one fixed score/exogenous stratum
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param context Fitted GLLRM context.
#' @param distribution Cached component distributions.
#' @param total_score Fixed total score.
#' @param count Number of response records to generate.
#' @param rng Private bootstrap RNG.
#' @return Integer response matrix with `count` rows.
#' @keywords internal
m2_m3_bootstrap_generate_group <- function(context,
                                            distribution,
                                            total_score,
                                            count,
                                            rng) {
  if (!is.null(distribution$joint)) {
    configurations <- which(distribution$joint$scores == total_score)
    responses <- matrix(0L, nrow = count, ncol = context$n_items)
    for (record in seq_len(count)) {
      selected <- m2_m3_bootstrap_weighted_index(
        distribution$joint$weights[configurations],
        rng,
        one_minus = FALSE
      )
      responses[record, ] <- distribution$joint$configurations[
        configurations[[selected]],
        ]
    }
    storage.mode(responses) <- "integer"
    return(responses)
  }

  components <- distribution$components
  n_components <- length(components)
  component_scores <- matrix(0L, nrow = count, ncol = n_components)

  # Source Generate_Compscores draws the first C-1 component scores for every
  # person and leaves the final component with the remaining total.
  for (record in seq_len(count)) {
    remaining <- as.integer(total_score)
    if (n_components > 1L) {
      for (component_index in seq_len(n_components - 1L)) {
        component <- components[[component_index]]
        possible <- seq.int(0L, min(component$maximum_score, remaining))
        rest <- distribution$suffix[[component_index + 1L]]
        weights <- component$gamma[possible + 1L] * rest[remaining - possible + 1L]
        selected <- m2_m3_bootstrap_weighted_index(weights, rng)
        component_score <- possible[[selected]]
        component_scores[record, component_index] <- component_score
        remaining <- remaining - component_score
      }
    }
    final <- components[[n_components]]
    if (remaining < 0L || remaining > final$maximum_score || final$gamma[[remaining + 1L]] <= 0) {
      stop("The fitted GLLRM cannot generate an observed total score.", call. = FALSE)
    }
    component_scores[record, n_components] <- remaining
  }

  responses <- matrix(0L, nrow = count, ncol = context$n_items)
  for (component_index in seq_along(components)) {
    component <- components[[component_index]]
    items <- component$items
    if (length(items) == 1L) {
      responses[, items] <- component_scores[, component_index]
      next
    }

    # SKbias8 fills the deterministic extreme component patterns first, then
    # loops interior component scores and draws one pattern per matching row.
    at_maximum <- component_scores[, component_index] == component$maximum_score
    if (any(at_maximum)) {
      responses[at_maximum, items] <- matrix(
        rep(context$item_raw_max[items] - 1L, each = sum(at_maximum)),
        nrow = sum(at_maximum)
      )
    }
    if (component$maximum_score > 1L) {
      for (component_score in seq.int(1L, component$maximum_score - 1L)) {
        records <- which(component_scores[, component_index] == component_score)
        if (!length(records)) {
          next
        }
        configurations <- which(component$scores == component_score)
        for (record in records) {
          selected <- m2_m3_bootstrap_weighted_index(
            component$weights[configurations],
            rng
          )
          responses[record, items] <- component$configurations[configurations[[selected]], ]
        }
      }
    }
  }
  storage.mode(responses) <- "integer"
  responses
}

#' Generate one complete CM2/CM3 parametric-bootstrap sample
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param context Fitted GLLRM context.
#' @param groups Fixed score/exogenous strata.
#' @param distribution_cache Cache from
#'   [new_m2_m3_bootstrap_distribution_cache()].
#' @param rng Private bootstrap RNG.
#' @return A list containing zero-based item scores, one-based exogenous values,
#'   total scores, RNG states, and draw counts.
#' @keywords internal
m2_m3_bootstrap_generate_sample <- function(context, groups, distribution_cache, rng) {
  total_records <- sum(groups$count)
  items <- matrix(0L, nrow = total_records, ncol = context$n_items)
  backgrounds <- matrix(0L, nrow = total_records, ncol = context$n_backgrounds)
  scores <- integer(total_records)
  start_state <- rng$state()
  start_draws <- rng$draws()
  next_row <- 1L
  dif_backgrounds <- context$dif_background_indices %||% integer()
  background_names <- as.character(context$backgrounds$name)
  generation_modes <- character()

  for (group_index in seq_len(nrow(groups))) {
    group <- groups[group_index, , drop = FALSE]
    count <- as.integer(group$count[[1L]])
    rows <- seq.int(next_row, length.out = count)
    score <- as.integer(group$score[[1L]])
    exogenous <- if (context$n_backgrounds) {
      as.integer(group[1L, background_names, drop = TRUE])
    } else {
      integer()
    }
    probability_exogenous <- exogenous
    if (context$n_backgrounds) {
      probability_exogenous[setdiff(seq_len(context$n_backgrounds), dif_backgrounds)] <- 1L
    }

    if (score == 0L) {
      generated <- matrix(0L, nrow = count, ncol = context$n_items)
    } else if (score == context$max_total_score) {
      generated <- matrix(
        rep(context$item_raw_max - 1L, each = count),
        nrow = count
      )
    } else {
      distribution <- distribution_cache(probability_exogenous)
      generation_modes <- c(generation_modes, distribution$generation_mode)
      generated <- m2_m3_bootstrap_generate_group(
        context,
        distribution,
        total_score = score,
        count = count,
        rng = rng
      )
    }
    items[rows, ] <- generated
    scores[rows] <- score
    if (context$n_backgrounds) {
      # Source SKbias8 writes Biasvalues, rather than Exovalues, for the two
      # extreme-score blocks. Non-DIF exogenous values are therefore collapsed
      # to category one in those blocks; retain this historical behavior.
      output_exogenous <- if (score %in% c(0L, context$max_total_score)) {
        probability_exogenous
      } else {
        exogenous
      }
      backgrounds[rows, ] <- matrix(
        rep(output_exogenous, each = count),
        nrow = count
      )
    }
    next_row <- next_row + count
  }
  if (nrow(items) && any(rowSums(items) != scores)) {
    stop("Internal CM2/CM3 bootstrap generation failed to preserve total scores.", call. = FALSE)
  }
  colnames(items) <- context$items$name
  colnames(backgrounds) <- background_names
  list(
    items = items,
    backgrounds = backgrounds,
    score = scores,
    start_state = start_state,
    final_state = rng$state(),
    draws = rng$draws() - start_draws,
    generation_mode = if (length(generation_modes)) {
      unique(generation_modes)[[1L]]
    } else {
      "deterministic_extreme_scores"
    }
  )
}

#' Convert a generated bootstrap sample to a source-shaped bundle
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param fit Public fitted model object.
#' @param sample Generated sample from [m2_m3_bootstrap_generate_sample()].
#' @return A bundle accepted by `build_gllrm_context()` and `fit_gllrm()`.
#' @keywords internal
m2_m3_bootstrap_bundle <- function(fit, sample) {
  project <- (fit$analysis %||% fit$spec$analysis)$project
  n_records <- nrow(sample$items)
  raw_columns <- max(
    ncol(project$raw_data),
    c(project$items$position, project$backgrounds$position, 0L)
  )
  raw <- matrix(-999L, nrow = n_records, ncol = raw_columns)
  for (item_index in seq_len(nrow(project$items))) {
    raw[, project$items$position[[item_index]]] <- sample$items[, item_index] + 1L
  }
  for (background_index in seq_len(nrow(project$backgrounds))) {
    raw[, project$backgrounds$position[[background_index]]] <-
      sample$backgrounds[, background_index]
  }
  project$raw_data <- raw
  project$source_data <- data.frame()
  build_item_parameters_bundle(project)
}
