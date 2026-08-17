# CM2/CM3 parametric-bootstrap implementation.
#
# The generator below is an independent R implementation of the component-wise
# conditional generator in source/PAS_skunits/SKbias8.pas. It deliberately does
# not call the Pascal harness. The harness remains an executable validation
# reference only.

#' Normalize public CM2/CM3 bootstrap controls
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param bootstrap Whether to run the source parametric bootstrap.
#' @param nsim Requested number of generated samples.
#' @param seed Optional positive Park--Miller stream seed.
#' @param reestimate Whether to refit the GLLRM in every generated sample.
#' @param bootstrap_max_step Maximum GLLRM steps per bootstrap refit.
#' @param keep_bootstrap_samples Whether to retain generated response matrices.
#' @param resample_score_distribution Source control whose preserved Pascal
#'   branch is empty; either value therefore keeps the observed distribution.
#' @return A normalized control list.
#' @keywords internal
normalize_m2_m3_bootstrap_control <- function(bootstrap,
                                               nsim,
                                               seed,
                                               reestimate,
                                               bootstrap_max_step,
                                               keep_bootstrap_samples,
                                               resample_score_distribution) {
  bootstrap <- normalize_public_logical(bootstrap, "bootstrap")
  reestimate <- normalize_public_logical(reestimate, "reestimate")
  keep_bootstrap_samples <- normalize_public_logical(
    keep_bootstrap_samples,
    "keep_bootstrap_samples"
  )
  resample_score_distribution <- normalize_public_logical(
    resample_score_distribution,
    "resample_score_distribution"
  )
  nsim <- normalize_public_integer_like(
    nsim,
    "`nsim` must be a single positive integer-like value.",
    scalar = TRUE,
    lower = 1L
  )
  bootstrap_max_step <- normalize_public_max_step(bootstrap_max_step)
  if (!is.null(seed)) {
    seed <- normalize_public_integer_like(
      seed,
      "`seed` must be NULL or a single integer-like value in 1..2147483646.",
      scalar = TRUE,
      lower = 1L,
      upper = 2147483646L
    )
  }
  list(
    enabled = bootstrap,
    nsim = nsim,
    seed = seed,
    reestimate = reestimate,
    max_step = bootstrap_max_step,
    keep_samples = keep_bootstrap_samples,
    resample_score_distribution = resample_score_distribution,
    acceptance_delta = 0.1
  )
}

#' Create the private seeded CM2/CM3 bootstrap stream
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param seed Optional normalized seed.
#' @return A list of closures for drawing, inspecting state, and counting draws.
#' @keywords internal
new_m2_m3_bootstrap_rng <- function(seed = NULL) {
  modulus <- 2147483647
  if (is.null(seed)) {
    # Source Randomize makes an unseeded run nondeterministic. Consume exactly
    # one draw from R's global stream, then keep every bootstrap draw private.
    seed <- sample.int(2147483646L, size = 1L)
  }
  stream <- new.env(parent = emptyenv())
  stream$seed <- as.integer(seed)
  stream$state <- as.double(seed)
  stream$draws <- 0L
  list(
    seed = function() stream$seed,
    uniform = function() {
      # A product is below 2^53, so binary64 evaluates this integer recurrence
      # exactly. This is the validation-harness Park--Miller stream, used as a
      # deterministic replacement for Pascal's runtime-specific Randomize.
      stream$state <- (48271 * stream$state) %% modulus
      if (stream$state <= 0) {
        stream$state <- 1
      }
      stream$draws <- stream$draws + 1L
      stream$state / modulus
    },
    state = function() as.integer(stream$state),
    draws = function() as.integer(stream$draws)
  )
}

#' Report whether a fitted context satisfies Pascal bootstrap bounds
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param context Fitted GLLRM context.
#' @return A list with `possible`, `reasons`, and the audited source bounds.
#' @keywords internal
m2_m3_bootstrap_capability <- function(context) {
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  component_sizes <- lengths(components)
  dif_backgrounds <- context$dif_background_indices %||% integer()
  dif_combinations <- if (length(dif_backgrounds)) {
    prod(as.double(context$background_raw_max[dif_backgrounds]))
  } else {
    1
  }
  exogenous_combinations <- if (context$n_backgrounds) {
    prod(as.double(context$background_raw_max))
  } else {
    1
  }
  n_complete <- length(source_complete_item_exogenous_rows(context))
  reasons <- character()
  if (context$n_items > 40L) {
    reasons <- c(reasons, "more than 40 items")
  }
  if (length(component_sizes) && any(component_sizes > 4L)) {
    reasons <- c(reasons, "an LD component contains more than four items")
  }
  if (length(context$item_raw_max) && any(context$item_raw_max - 1L > 7L)) {
    reasons <- c(reasons, "an item score exceeds seven")
  }
  if (dif_combinations > 64) {
    reasons <- c(reasons, "active DIF defines more than 64 probability strata")
  }
  if (exogenous_combinations > 216) {
    reasons <- c(reasons, "exogenous variables define more than 216 bootstrap strata")
  }
  if (n_complete == 0L) {
    reasons <- c(reasons, "there are no complete item/exogenous records")
  }
  list(
    possible = length(reasons) == 0L,
    reasons = reasons,
    bounds = list(
      maximum_items = 40L,
      maximum_component_items = 4L,
      maximum_item_score = 7L,
      maximum_dif_probability_strata = 64L,
      maximum_exogenous_strata = 216L
    ),
    observed = list(
      items = as.integer(context$n_items),
      largest_component = as.integer(max(component_sizes, 0L)),
      largest_item_score = as.integer(max(context$item_raw_max - 1L, 0L)),
      dif_probability_strata = as.numeric(dif_combinations),
      exogenous_strata = as.numeric(exogenous_combinations),
      complete_records = as.integer(n_complete)
    )
  )
}

#' Build Pascal-ordered score-by-exogenous bootstrap strata
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param context Fitted GLLRM context.
#' @return A data frame with fixed observed stratum counts. Exogenous variable
#'   one varies fastest, matching `LastBiasCombination1`; non-extreme scores
#'   precede the two extreme-score blocks within each stratum.
#' @keywords internal
m2_m3_bootstrap_groups <- function(context) {
  observed <- m2_m3_source_score_background_groups(context)
  background_names <- as.character(context$backgrounds$name)
  if (context$n_backgrounds) {
    grid <- expand.grid(
      lapply(context$background_raw_max, seq_len),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    names(grid) <- background_names
  } else {
    grid <- data.frame(.bootstrap_group = 1L)
  }
  observed_key <- if (nrow(observed)) {
    do.call(
      paste,
      c(observed[c("score", background_names)], sep = "\r")
    )
  } else {
    character()
  }
  observed_count <- stats::setNames(as.integer(observed$count), observed_key)
  score_order <- c(
    if (context$max_total_score > 1L) {
      seq.int(1L, context$max_total_score - 1L)
    } else {
      integer()
    },
    0L,
    context$max_total_score
  )
  rows <- vector("list", nrow(grid) * length(score_order))
  out_index <- 0L
  for (grid_row in seq_len(nrow(grid))) {
    values <- if (context$n_backgrounds) {
      as.integer(grid[grid_row, background_names, drop = TRUE])
    } else {
      integer()
    }
    for (score in score_order) {
      key <- paste(c(score, values), collapse = "\r")
      count <- unname(observed_count[key])
      if (!length(count) || is.na(count)) {
        count <- 0L
      }
      if (count <= 0L) {
        next
      }
      out_index <- out_index + 1L
      row <- data.frame(score = as.integer(score), count = as.integer(count))
      for (background_index in seq_along(background_names)) {
        row[[background_names[[background_index]]]] <- values[[background_index]]
      }
      rows[[out_index]] <- row
    }
  }
  if (out_index == 0L) {
    out <- data.frame(score = integer(), count = integer())
    for (name in background_names) {
      out[[name]] <- integer()
    }
    return(out)
  }
  do.call(rbind, rows[seq_len(out_index)])
}

#' Cache fitted component distributions for bootstrap generation
#'
#' Source trace: `source/PAS_skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param context Fitted GLLRM context.
#' @param state Current fitted GLLRM state.
#' @return A lookup closure keyed only by active-DIF background values.
#' @keywords internal
new_m2_m3_bootstrap_distribution_cache <- function(context, state) {
  cache <- new.env(parent = emptyenv())
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  joint_configuration_count <- prod(as.double(context$item_raw_max))
  function(background_values) {
    key <- gllrm_background_cache_key(context, background_values)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      details <- lapply(components, function(component_items) {
        component_key <- gllrm_component_key(component_items)
        gamma <- gllrm_component_gamma(context, state, component_items, background_values)
        list(
          items = as.integer(component_items),
          key = component_key,
          gamma = as.numeric(gamma$gamma),
          configurations = context$component_config_matrices[[component_key]],
          scores = as.integer(context$component_config_scores[[component_key]]),
          weights = as.numeric(gamma$config_weights),
          maximum_score = as.integer(sum(context$item_raw_max[component_items] - 1L))
        )
      })
      gammas <- lapply(details, `[[`, "gamma")
      n_components <- length(gammas)
      unit <- c(1, numeric(context$max_total_score))
      suffix <- vector("list", n_components + 1L)
      suffix[[n_components + 1L]] <- unit
      if (n_components) {
        for (component_index in rev(seq_len(n_components))) {
          suffix[[component_index]] <- convolve_score_vectors(
            gammas[[component_index]],
            suffix[[component_index + 1L]],
            context$max_total_score
          )
        }
      }
      joint <- NULL
      if (joint_configuration_count <= context$max_joint_configs) {
        # The bounded harness oracle enumerates item 1 outermost and the final
        # item innermost. Reverse expand.grid's columns to preserve that exact
        # cumulative-probability order for seeded parity.
        joint_frame <- expand.grid(
          rev(context$item_score_values),
          KEEP.OUT.ATTRS = FALSE
        )
        joint_matrix <- as.matrix(joint_frame[, rev(seq_len(context$n_items)), drop = FALSE])
        storage.mode(joint_matrix) <- "integer"
        joint_weights <- rep(1, nrow(joint_matrix))
        for (item in seq_len(context$n_items)) {
          item_scores <- joint_matrix[, item] + 1L
          joint_weights <- joint_weights * state$item_gamma[item, item_scores]
          dif_rows <- context$dif_by_item_matrices[[item]]
          if (nrow(dif_rows)) {
            for (dif_row in seq_len(nrow(dif_rows))) {
              background <- dif_rows[dif_row, 1L]
              dif_index <- dif_rows[dif_row, 2L]
              joint_weights <- joint_weights * state$dif_parameters[[dif_index]][
                item_scores,
                background_values[[background]]
              ]
            }
          }
        }
        for (ld_index in seq_along(context$ld_specs)) {
          ld <- context$ld_specs[[ld_index]]
          joint_weights <- joint_weights * state$ld_parameters[[ld_index]][
            cbind(joint_matrix[, ld$item1] + 1L, joint_matrix[, ld$item2] + 1L)
          ]
        }
        joint <- list(
          configurations = joint_matrix,
          scores = as.integer(rowSums(joint_matrix)),
          weights = as.numeric(joint_weights)
        )
      }
      assign(
        key,
        list(
          components = details,
          suffix = suffix,
          joint = joint,
          generation_mode = if (is.null(joint)) {
            "source_component_conditional"
          } else {
            "bounded_joint_harness_equivalent"
          }
        ),
        envir = cache
      )
    }
    get(key, envir = cache, inherits = FALSE)
  }
}
