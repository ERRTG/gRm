# CM2/CM3 parametric-bootstrap implementation.
#
# The generator below is an independent R implementation of the component-wise
# conditional generator in
# source/digram_source_20260817/skunits/SKbias8.pas. It deliberately does not
# call the Pascal harness. The harness remains an executable validation
# reference only.

#' Normalize public CM2/CM3 bootstrap controls
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' Canonical command trace:
#' `source/digram_source_20260817/scd/DIGRAM1f.pas`, command 197, which invokes
#' `skbias14.pas::CM3_analysis`.
#' The command always passes `Estimate = true`; the normalized `FALSE` branch
#' exposes the procedure's no-refit path as an explicit R extension.
#' @param bootstrap Whether to run the source parametric bootstrap.
#' @param nsim Requested number of generated samples.
#' @param seed Optional Delphi `RandSeed` bit pattern represented as an exact
#'   numeric value in `0..4294967295`.
#' @param reestimate Whether to refit the GLLRM in every generated sample. This
#'   controls the nested `CM3_analysis.DoSomething` fit, not the unconditional
#'   observed-data refit and footer in `SKbias8.Finish_random_Gllrm`.
#' @param bootstrap_max_step Maximum GLLRM steps per bootstrap refit.
#' @param bootstrap_jobs Maximum number of forked workers used for bootstrap
#'   refitting and diagnostic analysis. Sample generation remains serial.
#' @param keep_bootstrap_samples Whether to retain generated response matrices.
#' @param resample_score_distribution Source control whose preserved Pascal
#'   branch is empty; either value therefore keeps the observed distribution.
#' @return A normalized control list.
#' @keywords internal
normalize_cm2_cm3_bootstrap_control <- function(bootstrap,
                                               nsim,
                                               seed,
                                               reestimate,
                                               bootstrap_max_step,
                                               bootstrap_jobs,
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
  bootstrap_jobs <- normalize_public_jobs(bootstrap_jobs)
  if (!is.null(seed)) {
    # Delphi exposes RandSeed as a signed LongInt, but its LCG consumes all 32
    # bits. R numeric values represent every uint32 value exactly, including
    # bit patterns above R's signed-integer maximum.
    if (
      (!is.numeric(seed) && !is.integer(seed)) || length(seed) != 1L ||
        is.na(seed) || !is.finite(seed) || seed != floor(seed) ||
        seed < 0 || seed > 4294967295
    ) {
      stop(
        "`seed` must be NULL or a single integer-like value in 0..4294967295.",
        call. = FALSE
      )
    }
    seed <- as.double(seed)
  }
  list(
    enabled = bootstrap,
    nsim = nsim,
    seed = seed,
    reestimate = reestimate,
    max_step = bootstrap_max_step,
    jobs = bootstrap_jobs,
    keep_samples = keep_bootstrap_samples,
    resample_score_distribution = resample_score_distribution,
    acceptance_delta = 0.1
  )
}

#' Create the private Delphi CM2/CM3 bootstrap stream
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias14.pas::CM3_analysis`.
#' Source trace: `source/digram_source_20260817/skunits/SKbias8.pas::Generate_random_GLLRM_responses_with_exogene`.
#' @param seed Optional normalized seed.
#' @return A list of closures for drawing, inspecting state, and counting draws.
#' @keywords internal
new_cm2_cm3_bootstrap_rng <- function(seed = NULL) {
  modulus <- 4294967296
  if (is.null(seed)) {
    # DIGRAM was built with Delphi 4 (scd/scd.dof links VCL40). Its Randomize
    # initializes RandSeed from UTC milliseconds since midnight. CM3_analysis
    # calls Randomize immediately before the simulation loop, so this is the
    # effective seed; the earlier call in Start_random_Gllrm_with_exogene is
    # overwritten before any CM3 bootstrap draw is consumed.
    utc_seconds <- as.numeric(Sys.time()) %% 86400
    if (utc_seconds < 0) {
      utc_seconds <- utc_seconds + 86400
    }
    seed <- floor(utc_seconds * 1000)
  }
  stream <- new.env(parent = emptyenv())
  stream$seed <- as.double(seed)
  stream$state <- as.double(seed)
  stream$draws <- 0L

  next_state <- function(state) {
    # Delphi 4 System.Random advances the unsigned 32-bit RandSeed as
    #   RandSeed := 134775813 * RandSeed + 1 (mod 2^32).
    # The full product can exceed binary64's exact-integer range. Split both
    # operands into 16-bit limbs so every intermediate remains below 2^53 and
    # the wrapped uint32 result is bit-for-bit exact in an R double.
    state_high <- floor(state / 65536)
    state_low <- state - state_high * 65536
    low_product <- 33797 * state_low + 1
    result_low <- low_product %% 65536
    carry <- floor(low_product / 65536)
    result_high <- (2056 * state_low + 33797 * state_high + carry) %% 65536
    result_high * 65536 + result_low
  }

  list(
    seed = function() stream$seed,
    uniform = function() {
      stream$state <- next_state(stream$state)
      stream$draws <- stream$draws + 1L
      # Delphi's parameterless Random interprets the wrapped RandSeed as an
      # unsigned integer and scales it by exactly 2^-32, returning [0, 1).
      stream$state / modulus
    },
    state = function() stream$state,
    draws = function() as.integer(stream$draws)
  )
}

#' Report whether a fitted context satisfies Pascal bootstrap bounds
#'
#' Source trace: `source/digram_source_20260817/skunits/skbias1.pas::ReviewComponents`.
#' Source bounds: `source/digram_source_20260817/skunits/SKbias8.pas::MaxnumberOfMulticomp`
#' and `source/digram_source_20260817/skunits/SKTypes.pas::NsimMax`.
#' @param context Fitted GLLRM context.
#' @param state Current fitted GLLRM state.
#' @return A list with `possible`, `reasons`, and the audited source bounds.
#' @keywords internal
cm2_cm3_bootstrap_capability <- function(context, state) {
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
  largest_positive_pattern_count <- 0L
  can_count_patterns <- context$n_items <= 10L &&
    (!length(component_sizes) || all(component_sizes <= 4L)) &&
    (!length(context$item_raw_max) || all(context$item_raw_max - 1L <= 7L)) &&
    dif_combinations <= 64
  if (can_count_patterns && any(component_sizes > 1L)) {
    # SKbias1.CollectComprecords increments an Integer for every positive
    # fixed-score configuration, then assigns it to CompRecords.npatterns:
    # Byte. Count the same positive fitted weights over every active-DIF
    # probability stratum and reject before that assignment could exceed 255.
    if (length(dif_backgrounds)) {
      dif_grid <- expand.grid(
        lapply(context$background_raw_max[dif_backgrounds], seq_len),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      dif_grid <- matrix(integer(), nrow = 1L, ncol = 0L)
    }
    for (grid_row in seq_len(nrow(dif_grid))) {
      background_values <- rep(1L, context$n_backgrounds)
      if (length(dif_backgrounds)) {
        background_values[dif_backgrounds] <- as.integer(
          dif_grid[grid_row, , drop = TRUE]
        )
      }
      for (component_items in components[component_sizes > 1L]) {
        component <- gllrm_component_gamma(
          context,
          state,
          component_items,
          background_values
        )
        key <- gllrm_component_key(component_items)
        scores <- as.integer(context$component_config_scores[[key]])
        positive <- is.finite(component$config_weights) &
          component$config_weights > 0
        if (any(positive)) {
          fixed_score_counts <- tabulate(scores[positive] + 1L)
          largest_positive_pattern_count <- max(
            largest_positive_pattern_count,
            fixed_score_counts
          )
        }
      }
    }
  }
  reasons <- character()
  if (context$n_items > 10L) {
    # SKTypes.MaxItems8 is 10 in canonical DIGRAM 7.04, and
    # SKbias1.ReviewComponents rejects parametric generation above that bound.
    reasons <- c(reasons, "more than 10 items")
  }
  if (length(component_sizes) && any(component_sizes > 4L)) {
    reasons <- c(reasons, "an LD component contains more than four items")
  }
  if (sum(component_sizes > 1L) > 4L) {
    # SKbias8.TMultiPatterns is indexed by MaxnumberOfMulticomp = 4. This is a
    # bound on multi-item components, not on the total number of components.
    reasons <- c(reasons, "more than four multi-item LD components")
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
  if (n_complete > 15000L) {
    # SKTypes.NsimMax sizes Tsimdata, the generated bootstrap record store.
    reasons <- c(reasons, "more than 15000 complete item/exogenous records")
  }
  if (largest_positive_pattern_count > 255L) {
    reasons <- c(
      reasons,
      "a fixed component score has more than 255 positive response patterns"
    )
  }
  list(
    possible = length(reasons) == 0L,
    reasons = reasons,
    bounds = list(
      maximum_items = 10L,
      maximum_component_items = 4L,
      maximum_multi_item_components = 4L,
      maximum_positive_fixed_score_patterns = 255L,
      maximum_item_score = 7L,
      maximum_dif_probability_strata = 64L,
      maximum_exogenous_strata = 216L,
      maximum_complete_records = 15000L
    ),
    observed = list(
      items = as.integer(context$n_items),
      largest_component = as.integer(max(component_sizes, 0L)),
      multi_item_components = as.integer(sum(component_sizes > 1L)),
      largest_positive_fixed_score_pattern_count = as.integer(
        largest_positive_pattern_count
      ),
      largest_item_score = as.integer(max(context$item_raw_max - 1L, 0L)),
      dif_probability_strata = as.numeric(dif_combinations),
      exogenous_strata = as.numeric(exogenous_combinations),
      complete_records = as.integer(n_complete)
    )
  )
}

#' Build Pascal-ordered score-by-exogenous bootstrap strata
#'
#' Source trace: `source/digram_source_20260817/skunits/SKbias8.pas::Generate_GLLRM_boot_sample_with_exogene`.
#' @param context Fitted GLLRM context.
#' @return A data frame with fixed observed stratum counts. Exogenous variable
#'   one varies fastest, matching `LastBiasCombination1`; non-extreme scores
#'   precede the two extreme-score blocks within each stratum.
#' @keywords internal
cm2_cm3_bootstrap_groups <- function(context) {
  observed <- cm2_cm3_source_score_background_groups(context)
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
#' Source trace: `source/digram_source_20260817/skunits/SKbias8.pas::CalculateCondprobs`.
#' Pattern trace: `source/digram_source_20260817/skunits/skbias1.pas::CollectComprecords`.
#' Consumption trace:
#' `source/digram_source_20260817/skunits/SKbias8.pas::Generate_random_GLLRM_responses_with_exogene`.
#' @param context Fitted GLLRM context.
#' @param state Current fitted GLLRM state.
#' @return A lookup closure keyed only by active-DIF background values.
#' @keywords internal
new_cm2_cm3_bootstrap_distribution_cache <- function(context, state) {
  cache <- new.env(parent = emptyenv())
  components <- context$ld_components_items %||% gllrm_ld_components(context)$items
  function(background_values) {
    key <- gllrm_background_cache_key(context, background_values)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      details <- lapply(components, function(component_items) {
        component_key <- gllrm_component_key(component_items)
        gamma <- gllrm_component_gamma(context, state, component_items, background_values)
        configurations <- context$component_config_matrices[[component_key]]
        scores <- as.integer(context$component_config_scores[[component_key]])
        weights <- as.numeric(gamma$config_weights)
        component_gamma <- as.numeric(gamma$gamma)
        component_gamma_raw <- component_gamma * as.numeric(gamma$scale)

        # SKbias1.Generate_CompPatterns/CollectComprecords nests item-score
        # loops with the first component item outermost and the last innermost.
        # expand.grid() uses the reverse convention, so sort lexicographically
        # before any cumulative pattern draw. This ordering is numerically
        # material even though the conditional distribution is unchanged.
        source_order <- do.call(order, as.data.frame(configurations))
        configurations <- configurations[source_order, , drop = FALSE]
        scores <- scores[source_order]
        weights <- weights[source_order]
        maximum_score <- as.integer(sum(context$item_raw_max[component_items] - 1L))

        # SKbias1.CalculateCompProbs divides each positive configuration weight
        # by the gamma for its component score; CollectComprecords then sums
        # those already-normalized probabilities in source configuration order.
        # Cache that exact representation instead of rescaling x by a raw
        # weight total during generation.
        pattern_tables <- lapply(seq.int(0L, maximum_score), function(score) {
          indices <- which(scores == score & weights > 0)
          if (!length(indices)) {
            return(list(indices = integer(), cumulative = numeric()))
          }
          gamma_score <- component_gamma_raw[[score + 1L]]
          if (!is.finite(gamma_score) || gamma_score <= 0) {
            return(list(indices = integer(), cumulative = numeric()))
          }
          probabilities <- weights[indices] / gamma_score
          list(indices = indices, cumulative = cumsum(probabilities))
        })
        list(
          items = as.integer(component_items),
          key = component_key,
          gamma = component_gamma,
          configurations = configurations,
          scores = scores,
          weights = weights,
          maximum_score = maximum_score,
          pattern_tables = pattern_tables
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
      component_score_tables <- vector("list", max(n_components - 1L, 0L))
      if (n_components > 1L) {
        for (component_index in seq_len(n_components - 1L)) {
          component <- details[[component_index]]
          rest <- suffix[[component_index + 1L]]
          denominator <- suffix[[component_index]]
          component_score_tables[[component_index]] <- lapply(
            seq.int(0L, context$max_total_score),
            function(remaining) {
              possible <- seq.int(0L, min(component$maximum_score, remaining))
              gamma0 <- denominator[[remaining + 1L]]
              if (!is.finite(gamma0) || gamma0 <= 0) {
                return(list(scores = possible, cumulative = rep(0, length(possible))))
              }
              # SKbias8.CalculateCondprobs evaluates gamma1 * (gamma2/gamma0)
              # in this order and accumulates from component score zero upward.
              probabilities <- component$gamma[possible + 1L] *
                (rest[remaining - possible + 1L] / gamma0)
              list(scores = possible, cumulative = cumsum(probabilities))
            }
          )
        }
      }
      assign(
        key,
        list(
          components = details,
          suffix = suffix,
          component_score_tables = component_score_tables,
          generation_mode = "source_component_conditional"
        ),
        envir = cache
      )
    }
    get(key, envir = cache, inherits = FALSE)
  }
}
