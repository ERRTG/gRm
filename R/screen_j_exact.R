#' Source-shaped Monte Carlo exact partial gamma p-value
#'
#' Source trace: `SKbias3.XYZ_bias_ANALYSE` calls `SKrandom.GENTAB1` for each
#' score-conditioned slice, adds simulated `PPQ` and `PMQ` over slices, and
#' counts `abs(simulated_gamma) >= abs(observed_gamma)`.
#'
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param strata Three-way table with dimensions item/background by item/exo by
#'   conditioning score.
#' @param observed_gamma Observed partial gamma.
#' @param nsim Number of simulated tables.
#' @return Two-sided Monte Carlo exact p-value.
#' @param seed Random-stream seed.
#' @param observed_chi Internal `observed_chi` value used by this helper.
#' @param sequential Internal `sequential` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @keywords internal
#' @noRd
screen_j_exact_partial_gamma <- function(strata,
                                         observed_gamma,
                                         nsim,
                                         seed = NULL,
                                         observed_chi = NULL,
                                         sequential = FALSE,
                                         seq_limit = nsim,
                                         seq_p0 = 0.05,
                                         seq_boundary = 1.058) {
  slices <- screen_j_strata_slices(strata)
  if (!is.null(observed_chi)) {
    return(screen_j_exact_chi_gamma_slices(
      slices,
      observed_chi,
      observed_gamma,
      nsim,
      seed,
      sequential = sequential,
      seq_limit = seq_limit,
      seq_p0 = seq_p0,
      seq_boundary = seq_boundary
    )$p_gamma)
  }
  screen_j_exact_gamma_slices(slices, observed_gamma, nsim, seed)
}

#' Internal screen j exact gamma slices helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param slices Internal `slices` value used by this helper.
#' @param observed_gamma Internal `observed_gamma` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @return The internal `screen_j_exact_gamma_slices()` computation result.
#' @keywords internal
#' @noRd
screen_j_exact_gamma_slices <- function(slices, observed_gamma, nsim, seed = NULL) {
  native <- screen_j_exact_gamma_slices_native(slices, observed_gamma, nsim, seed)
  if (!is.null(native)) {
    return(as.numeric(native[[1L]]))
  }
  screen_j_exact_gamma_prepared(screen_j_prepare_exact_slices(slices), observed_gamma, nsim, seed)
}

#' Internal screen j exact gamma prepared helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param prepared_slices Internal `prepared_slices` value used by this helper.
#' @param observed_gamma Internal `observed_gamma` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @return The internal `screen_j_exact_gamma_prepared()` computation result.
#' @keywords internal
#' @noRd
screen_j_exact_gamma_prepared <- function(prepared_slices, observed_gamma, nsim, seed = NULL) {
  screen_j_exact_gamma_prepared_r(prepared_slices, observed_gamma, nsim, seed)
}

#' Internal screen j exact gamma prepared r helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param prepared_slices Internal `prepared_slices` value used by this helper.
#' @param observed_gamma Internal `observed_gamma` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @return The internal `screen_j_exact_gamma_prepared_r()` computation result.
#' @keywords internal
#' @noRd
screen_j_exact_gamma_prepared_r <- function(prepared_slices, observed_gamma, nsim, seed = NULL) {
  random_draw <- screen_j_source_random_stream(seed)
  exceed <- 0L
  observed_gamma <- as.numeric(observed_gamma)
  for (sim in seq_len(nsim)) {
    ppq_total <- 0
    pmq_total <- 0
    for (prepared in prepared_slices) {
      generated <- exo_select_gentab1_prepared(prepared, random_draw = random_draw)
      stats <- screen_rc_gamma_counts(generated)
      ppq_total <- ppq_total + stats$ppq
      pmq_total <- pmq_total + stats$pmq
    }
    simulated_gamma <- if (ppq_total > 0) pmq_total / ppq_total else 0
    if (abs(simulated_gamma) >= abs(observed_gamma)) {
      exceed <- exceed + 1L
    }
  }
  exceed / nsim
}

#' Source-shaped Monte Carlo exact partial chi-square p-value
#'
#' Source trace: `SKbias3.XYZ_bias_ANALYSE` adds simulated chi-square statistics
#' over score-conditioned slices and counts `simulated_chi >= observed_chi`.
#' When the DIGRAM exact command state is sequential/repeated, the same
#' `SEQUENTIAL`, `SEQ_P0`, `SEQ_B`, and `seq_limit` state is used while
#' evaluating simulated chi results.
#'
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param strata Three-way table.
#' @param observed_chi Observed partial chi-square.
#' @param nsim Number of simulated tables.
#' @return Monte Carlo exact p-value.
#' @param seed Random-stream seed.
#' @param sequential Internal `sequential` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @keywords internal
#' @noRd
screen_j_exact_partial_chi <- function(strata,
                                       observed_chi,
                                       nsim,
                                       seed = NULL,
                                       sequential = FALSE,
                                       seq_limit = nsim,
                                       seq_p0 = 0.05,
                                       seq_boundary = 1.058) {
  screen_j_exact_chi_slices(
    screen_j_strata_slices(strata),
    observed_chi,
    nsim,
    seed,
    sequential = sequential,
    seq_limit = seq_limit,
    seq_p0 = seq_p0,
    seq_boundary = seq_boundary
  )
}

#' Internal screen j exact chi slices helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param slices Internal `slices` value used by this helper.
#' @param observed_chi Internal `observed_chi` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param sequential Internal `sequential` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_exact_chi_slices()` computation result.
#' @keywords internal
#' @noRd
screen_j_exact_chi_slices <- function(slices,
                                      observed_chi,
                                      nsim,
                                      seed = NULL,
                                      sequential = FALSE,
                                      seq_limit = nsim,
                                      seq_p0 = 0.05,
                                      seq_boundary = 1.058) {
  use_native <- !isTRUE(sequential)
  native <- if (use_native) {
    screen_j_exact_chi_slices_native(
      slices,
      observed_chi,
      nsim,
      seed,
      sequential = sequential,
      seq_limit = seq_limit
    )
  } else {
    NULL
  }
  if (!is.null(native)) {
    return(as.numeric(native[[1L]]))
  }
  screen_j_exact_chi_prepared(
    screen_j_prepare_exact_slices(slices),
    observed_chi,
    nsim,
    seed,
    sequential = sequential,
    seq_limit = seq_limit,
    seq_p0 = seq_p0,
    seq_boundary = seq_boundary
  )
}

#' Internal screen j exact chi prepared helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param prepared_slices Internal `prepared_slices` value used by this helper.
#' @param observed_chi Internal `observed_chi` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param sequential Internal `sequential` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_exact_chi_prepared()` computation result.
#' @keywords internal
#' @noRd
screen_j_exact_chi_prepared <- function(prepared_slices,
                                        observed_chi,
                                        nsim,
                                        seed = NULL,
                                        sequential = FALSE,
                                        seq_limit = nsim,
                                        seq_p0 = 0.05,
                                        seq_boundary = 1.058) {
  screen_j_exact_chi_prepared_r(
    prepared_slices,
    observed_chi,
    nsim,
    seed,
    sequential = sequential,
    seq_limit = seq_limit,
    seq_p0 = seq_p0,
    seq_boundary = seq_boundary
  )
}

#' Internal screen j exact chi prepared r helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param prepared_slices Internal `prepared_slices` value used by this helper.
#' @param observed_chi Internal `observed_chi` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param sequential Internal `sequential` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_exact_chi_prepared_r()` computation result.
#' @keywords internal
#' @noRd
screen_j_exact_chi_prepared_r <- function(prepared_slices,
                                          observed_chi,
                                          nsim,
                                          seed = NULL,
                                          sequential = FALSE,
                                          seq_limit = nsim,
                                          seq_p0 = 0.05,
                                          seq_boundary = 1.058) {
  random_draw <- screen_j_source_random_stream(seed)
  exceed <- 0L
  seq_p0 <- as.numeric(seq_p0[[1L]])
  seq_boundary <- as.numeric(seq_boundary[[1L]])
  seq_limit <- as.integer(seq_limit[[1L]])
  observed_chi <- screen_j_source_single(observed_chi)
  for (sim in seq_len(nsim)) {
    chi_total <- 0
    for (prepared in prepared_slices) {
      generated <- exo_select_gentab1_prepared(prepared, random_draw = random_draw)
      slice_chi <- screen_j_rc_chi_square_prepared_expected(generated, prepared)
      chi_total <- chi_total + slice_chi
    }
    if (chi_total >= observed_chi) {
      exceed <- exceed + 1L
    }
    if (isTRUE(sequential) &&
        (screen_j_seq_t(exceed, sim, seq_p0) >= seq_boundary || exceed >= seq_limit)) {
      break
    }
  }
  exceed / sim
}

#' Internal screen j exact chi gamma slices helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param slices Internal `slices` value used by this helper.
#' @param observed_chi Internal `observed_chi` value used by this helper.
#' @param observed_gamma Internal `observed_gamma` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param sequential Internal `sequential` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_exact_chi_gamma_slices()` computation result.
#' @keywords internal
#' @noRd
screen_j_exact_chi_gamma_slices <- function(slices,
                                            observed_chi,
                                            observed_gamma,
                                            nsim,
                                            seed = NULL,
                                            sequential = FALSE,
                                            seq_limit = nsim,
                                            seq_p0 = 0.05,
                                            seq_boundary = 1.058) {
  use_native <- !isTRUE(sequential) ||
    (abs(as.numeric(seq_p0[[1L]]) - 0.05) < 1e-15 &&
      abs(as.numeric(seq_boundary[[1L]]) - 1.058) < 1e-15)
  native <- if (use_native) {
    screen_j_exact_chi_gamma_slices_native(
      slices,
      observed_chi,
      observed_gamma,
      nsim,
      seed,
      sequential = sequential,
      seq_limit = seq_limit
    )
  } else {
    NULL
  }
  if (!is.null(native)) {
    return(native)
  }
  prepared_slices <- screen_j_prepare_exact_slices(slices)
  screen_j_exact_chi_gamma_prepared_r(
    prepared_slices,
    observed_chi,
    observed_gamma,
    nsim,
    seed,
    sequential = sequential,
    seq_limit = seq_limit,
    seq_p0 = seq_p0,
    seq_boundary = seq_boundary
  )
}

#' Internal screen j exact chi gamma prepared r helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param prepared_slices Internal `prepared_slices` value used by this helper.
#' @param observed_chi Internal `observed_chi` value used by this helper.
#' @param observed_gamma Internal `observed_gamma` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param sequential Internal `sequential` value used by this helper.
#' @param random_draw Internal `random_draw` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @param trace Internal `trace` value used by this helper.
#' @return The internal `screen_j_exact_chi_gamma_prepared_r()` computation result.
#' @keywords internal
#' @noRd
screen_j_exact_chi_gamma_prepared_r <- function(prepared_slices,
                                                observed_chi,
                                                observed_gamma,
                                                nsim,
                                                seed = NULL,
                                                sequential = FALSE,
                                                random_draw = NULL,
                                                seq_limit = nsim,
                                                seq_p0 = 0.05,
                                                seq_boundary = 1.058,
                                                trace = FALSE) {
  source_stream <- is.null(random_draw)
  if (source_stream) {
    random_draw <- screen_j_source_random_stream(seed)
  }
  draw_count <- 0L
  current_table_draws <- numeric()
  counted_draw <- function() {
    draw_count <<- draw_count + 1L
    value <- random_draw()
    if (isTRUE(trace)) {
      current_table_draws <<- c(current_table_draws, value)
    }
    value
  }
  chi_exceed <- 0L
  gamma_exceed <- 0L
  chi_status <- FALSE
  gamma_status <- FALSE
  table_trace <- if (isTRUE(trace)) list() else NULL
  simulation_trace <- if (isTRUE(trace)) vector("list", nsim) else NULL
  seq_p0 <- as.numeric(seq_p0[[1L]])
  seq_boundary <- as.numeric(seq_boundary[[1L]])
  seq_limit <- as.integer(seq_limit[[1L]])
  # Source trace: SKbias3.XYZ_bias_ANALYSE assigns observed CHITOT to
  # RESULTS[1,1] (SKTypes.RESARRAY is SINGLE) before the exact comparison.
  # Simulated CHI/PPQ/PMQ/GAMMA and AbsGammaTot remain Pascal REAL (Double).
  observed_chi_source <- screen_j_source_single(observed_chi)
  observed_gamma_source <- as.numeric(observed_gamma)
  for (sim in seq_len(nsim)) {
    chi_total <- 0
    ppq_total <- 0
    pmq_total <- 0
    for (slice_index in seq_along(prepared_slices)) {
      prepared <- prepared_slices[[slice_index]]
      current_table_draws <- numeric()
      generated <- exo_select_gentab1_prepared(prepared, random_draw = counted_draw)
      slice_chi <- screen_j_rc_chi_square_prepared_expected(generated, prepared)
      chi_total <- chi_total + slice_chi
      gamma_counts <- screen_rc_gamma_counts(generated)
      slice_ppq <- gamma_counts$ppq
      slice_pmq <- gamma_counts$pmq
      ppq_total <- ppq_total + slice_ppq
      pmq_total <- pmq_total + slice_pmq
      if (isTRUE(trace)) {
        table_trace[[length(table_trace) + 1L]] <- list(
          sim = as.integer(sim),
          slice = as.integer(slice_index),
          table = generated,
          random_draws = current_table_draws,
          chi_square = slice_chi,
          ppq = slice_ppq,
          pmq = slice_pmq,
          draw_count = draw_count,
          final_seed = if (source_stream) {
            screen_j_source_seed_after_draws(seed, draw_count)
          } else {
            NA_real_
          }
        )
      }
    }
    chi_ge <- chi_total >= observed_chi_source
    if (chi_ge) {
      chi_exceed <- chi_exceed + 1L
    }
    simulated_gamma <- if (ppq_total > 0) pmq_total / ppq_total else 0
    gamma_ge <- abs(simulated_gamma) >= abs(observed_gamma_source)
    if (gamma_ge) {
      gamma_exceed <- gamma_exceed + 1L
    }
    if (isTRUE(sequential) &&
        (screen_j_seq_t(chi_exceed, sim, seq_p0) >= seq_boundary || chi_exceed >= seq_limit)) {
      chi_status <- TRUE
    }
    if (isTRUE(sequential) &&
        (screen_j_seq_t(gamma_exceed, sim, seq_p0) >= seq_boundary || gamma_exceed >= seq_limit)) {
      gamma_status <- TRUE
    }
    stop_now <- isTRUE(sequential) && chi_status && gamma_status
    if (isTRUE(trace)) {
      simulation_trace[[sim]] <- data.frame(
        sim = as.integer(sim),
        chi_square = chi_total,
        ppq = ppq_total,
        pmq = pmq_total,
        gamma = simulated_gamma,
        chi_ge_observed = chi_ge,
        gamma_ge_observed = gamma_ge,
        chi_exceed = chi_exceed,
        gamma_exceed = gamma_exceed,
        chi_status = chi_status,
        gamma_status = gamma_status,
        stop = stop_now,
        draw_count = draw_count,
        final_seed = if (source_stream) {
          screen_j_source_seed_after_draws(seed, draw_count)
        } else {
          NA_real_
        }
      )
    }
    if (stop_now) {
      break
    }
  }
  result <- list(
    p_chi = chi_exceed / sim,
    p_gamma = gamma_exceed / sim,
    nsim = as.integer(sim),
    chi_exceed = chi_exceed,
    gamma_exceed = gamma_exceed,
    draw_count = draw_count,
    rng_draws = draw_count,
    final_seed = if (source_stream) {
      screen_j_source_seed_after_draws(seed, draw_count)
    } else {
      NA_real_
    }
  )
  if (isTRUE(trace)) {
    result$trajectory <- list(
      tables = table_trace,
      simulations = do.call(rbind, simulation_trace[seq_len(sim)])
    )
  }
  result
}

#' Internal screen j seq t helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param exceed Internal `exceed` value used by this helper.
#' @param sim Internal `sim` value used by this helper.
#' @param p0 Internal `p0` value used by this helper.
#' @return The internal `screen_j_seq_t()` computation result.
#' @keywords internal
#' @noRd
screen_j_seq_t <- function(exceed, sim, p0) {
  if (sim < 21L) {
    return(0)
  }
  root <- sqrt(sim)
  exceed / root - root * p0
}
