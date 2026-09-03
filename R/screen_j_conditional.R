#' Source trace: source/digram_source_20260817/skunits/SKbias3.pas::XYZ_bias_ANALYSE receives the
#' conditional item/background table prepared by
#' source/digram_source_20260817/skunits/SKxyz1.PAS::MAKE_XYZ_TABLE or
#' source/digram_source_20260817/skunits/SKbigtab.pas::Transfer_BT_to_XYZ_TABLE. It then computes
#' chi-square, gamma, asymptotic p-values, and optional exact/repeated
#' random-table summaries. This R helper is the shared implementation for those
#' screen J conditional tests.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param x Object or value to process.
#' @param y Internal `y` value used by this helper.
#' @param x_dim Internal `x_dim` value used by this helper.
#' @param y_dim Internal `y_dim` value used by this helper.
#' @param condition_values Internal `condition_values` value used by this helper.
#' @param condition_dims Internal `condition_dims` value used by this helper.
#' @param valid Internal `valid` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param repeated Whether to use repeated sequential simulation.
#' @param random_draw Internal `random_draw` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param use_native Internal `use_native` value used by this helper.
#' @return The internal `screen_j_conditional_try_native()` computation result.
#' @keywords internal
#' @noRd
screen_j_conditional_try_native <- function(x,
                                            y,
                                            x_dim,
                                            y_dim,
                                            condition_values,
                                            condition_dims,
                                            valid,
                                            exact,
                                            nsim,
                                            seed,
                                            repeated,
                                            random_draw,
                                            seq_limit,
                                            use_native) {
  if (!isTRUE(use_native) || !isTRUE(exact) || !is.null(random_draw) || !screen_j_exact_native_available()) {
    return(NULL)
  }
  screen_j_conditional_bias_test_native(
    x = x,
    y = y,
    x_dim = x_dim,
    y_dim = y_dim,
    condition_values = condition_values,
    condition_dims = condition_dims,
    valid = valid,
    exact = exact,
    nsim = nsim,
    seed = seed,
    repeated = repeated,
    seq_limit = seq_limit
  )
}

#' Internal screen j conditional valid rows helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param x Object or value to process.
#' @param y Internal `y` value used by this helper.
#' @param x_dim Internal `x_dim` value used by this helper.
#' @param y_dim Internal `y_dim` value used by this helper.
#' @param condition_values Internal `condition_values` value used by this helper.
#' @param condition_dims Internal `condition_dims` value used by this helper.
#' @param valid Internal `valid` value used by this helper.
#' @return The internal `screen_j_conditional_valid_rows()` computation result.
#' @keywords internal
#' @noRd
screen_j_conditional_valid_rows <- function(x,
                                            y,
                                            x_dim,
                                            y_dim,
                                            condition_values,
                                            condition_dims,
                                            valid) {
  condition_values <- as.matrix(condition_values)
  condition_dims <- as.integer(condition_dims)
  keep <- valid & x >= 1L & x <= x_dim & y >= 1L & y <= y_dim
  if (ncol(condition_values) > 0L) {
    for (condition_index in seq_len(ncol(condition_values))) {
      keep <- keep &
        condition_values[, condition_index] >= 1L &
        condition_values[, condition_index] <= condition_dims[[condition_index]]
    }
  }
  keep
}

#' Internal screen j conditional stratum index helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param condition_values Internal `condition_values` value used by this helper.
#' @param condition_dims Internal `condition_dims` value used by this helper.
#' @return The internal `screen_j_conditional_stratum_index()` computation result.
#' @keywords internal
#' @noRd
screen_j_conditional_stratum_index <- function(condition_values, condition_dims) {
  condition_values <- as.matrix(condition_values)
  condition_dims <- as.integer(condition_dims)
  if (ncol(condition_values) == 0L) {
    return(list(values = 1, id = rep(1, nrow(condition_values))))
  }

  condition_id <- rep(1, nrow(condition_values))
  multiplier <- 1
  for (condition_column in seq_len(ncol(condition_values))) {
    condition_id <- condition_id + (condition_values[, condition_column] - 1L) * multiplier
    multiplier <- multiplier * condition_dims[[condition_column]]
  }
  list(values = sort(unique(condition_id)), id = condition_id)
}

#' Internal screen j conditional empty result helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param nsim Requested simulation count.
#' @return The internal `screen_j_conditional_empty_result()` computation result.
#' @keywords internal
#' @noRd
screen_j_conditional_empty_result <- function(exact, nsim) {
  list(
    chi_square = 0,
    df = 0L,
    p_chi = 1,
    gamma = 0,
    p_gamma = 1,
    p_chi_asymp = 1,
    p_gamma_asymp = 1,
    p_chi_exact = if (isTRUE(exact)) 1 else NA_real_,
    p_gamma_exact = if (isTRUE(exact)) 1 else NA_real_,
    exact_nsim = if (isTRUE(exact)) nsim else 0L,
    ppq = 0,
    pmq = 0,
    s = 0
  )
}

#' Internal screen j conditional slice stats helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param x Object or value to process.
#' @param y Internal `y` value used by this helper.
#' @param x_dim Internal `x_dim` value used by this helper.
#' @param y_dim Internal `y_dim` value used by this helper.
#' @param condition_values Internal `condition_values` value used by this helper.
#' @param condition_dims Internal `condition_dims` value used by this helper.
#' @param valid Internal `valid` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @return The internal `screen_j_conditional_slice_stats()` computation result.
#' @keywords internal
#' @noRd
screen_j_conditional_slice_stats <- function(x,
                                             y,
                                             x_dim,
                                             y_dim,
                                             condition_values,
                                             condition_dims,
                                             valid,
                                             exact = FALSE) {
  condition_values <- as.matrix(condition_values)
  condition_dims <- as.integer(condition_dims)
  keep <- screen_j_conditional_valid_rows(
    x,
    y,
    x_dim,
    y_dim,
    condition_values,
    condition_dims,
    valid
  )
  if (!any(keep)) {
    return(list(has_rows = FALSE))
  }

  condition_index <- screen_j_conditional_stratum_index(
    condition_values[keep, , drop = FALSE],
    condition_dims
  )
  x_keep <- x[keep]
  y_keep <- y[keep]
  chi_total <- 0
  df_total <- 0L
  ppq_total <- 0
  pmq_total <- 0
  s_total <- 0
  slices <- list()

  for (condition_key in condition_index$values) {
    in_stratum <- condition_index$id == condition_key
    tab <- matrix(0L, nrow = x_dim, ncol = y_dim)
    index <- x_keep[in_stratum] + (y_keep[in_stratum] - 1L) * x_dim
    tab[] <- tabulate(index, nbins = x_dim * y_dim)
    if (!isTRUE(exact) || screen_j_source_informative_slice(tab)) {
      # Source trace: source/digram_source_20260817/skunits/SKxyz1.PAS::MAKE_XYZ_TABLE only
      # materializes a conditioning slice when both tested variables have at
      # least two nonempty categories. Exact GENTAB1 simulations consume that
      # reduced XYZ slice set before running the fixed or sequential exact
      # tests.
      slices[[length(slices) + 1L]] <- tab
    }
    chi <- screen_j_rc_chi_source_expected(tab)
    gamma <- screen_rc_gamma(tab)
    chi_total <- chi_total + chi$chi_square
    df_total <- df_total + chi$df
    ppq_total <- ppq_total + gamma$ppq
    pmq_total <- pmq_total + gamma$pmq
    s_total <- s_total + gamma$s
  }

  gamma <- if (ppq_total > 0) pmq_total / ppq_total else 0
  s <- if (ppq_total > 0) s_total / ppq_total / ppq_total else 0
  p_gamma <- if (ppq_total <= 0) {
    2
  } else if (s <= 0) {
    2 * source_tail_norm(4, TRUE)
  } else {
    2 * source_tail_norm(abs(gamma / sqrt(s)), TRUE)
  }
  p_chi <- if (df_total > 0L) source_pfchi(df_total, chi_total) else 1
  list(
    has_rows = TRUE,
    chi_square = chi_total,
    df = df_total,
    gamma = gamma,
    p_chi_asymp = p_chi,
    p_gamma_asymp = p_gamma,
    ppq = ppq_total,
    pmq = pmq_total,
    s = s,
    slices = slices
  )
}

#' Internal screen j conditional exact results helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param slices Internal `slices` value used by this helper.
#' @param chi_total Internal `chi_total` value used by this helper.
#' @param gamma Internal `gamma` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param repeated Whether to use repeated sequential simulation.
#' @param random_draw Internal `random_draw` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @param use_native Internal `use_native` value used by this helper.
#' @return The internal `screen_j_conditional_exact_results()` computation result.
#' @keywords internal
#' @noRd
screen_j_conditional_exact_results <- function(slices,
                                               chi_total,
                                               gamma,
                                               nsim,
                                               seed,
                                               repeated,
                                               random_draw,
                                               seq_limit,
                                               seq_p0,
                                               seq_boundary,
                                               use_native) {
  # Source trace: source/digram_source_20260817/skunits/SKrandom.pas::GENTAB1 is the exact
  # random-table generator behind SKbias3.XYZ_bias_ANALYSE. Native and R
  # branches must consume the same prepared slices, seed, sequential controls,
  # and draw order.
  if (isTRUE(use_native)) {
    screen_j_exact_chi_gamma_slices(
      slices,
      chi_total,
      gamma,
      nsim,
      seed,
      sequential = repeated,
      seq_limit = seq_limit,
      seq_p0 = seq_p0,
      seq_boundary = seq_boundary
    )
  } else {
    screen_j_exact_chi_gamma_prepared_r(
      screen_j_prepare_exact_slices(slices),
      chi_total,
      gamma,
      nsim,
      seed,
      sequential = repeated,
      random_draw = random_draw,
      seq_limit = seq_limit,
      seq_p0 = seq_p0,
      seq_boundary = seq_boundary
    )
  }
}

#' Internal screen j conditional bias test helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param x Object or value to process.
#' @param y Internal `y` value used by this helper.
#' @param x_dim Internal `x_dim` value used by this helper.
#' @param y_dim Internal `y_dim` value used by this helper.
#' @param condition_values Internal `condition_values` value used by this helper.
#' @param condition_dims Internal `condition_dims` value used by this helper.
#' @param valid Internal `valid` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param repeated Whether to use repeated sequential simulation.
#' @param native Internal `native` value used by this helper.
#' @param random_draw Internal `random_draw` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_conditional_bias_test()` computation result.
#' @keywords internal
#' @noRd
screen_j_conditional_bias_test <- function(x,
                                           y,
                                           x_dim,
                                           y_dim,
                                           condition_values,
                                           condition_dims,
                                           valid,
                                           exact = FALSE,
                                           nsim = 1000L,
                                           seed = NULL,
                                           repeated = FALSE,
                                           native = TRUE,
                                           random_draw = NULL,
                                           seq_limit = nsim,
                                           seq_p0 = 0.05,
                                           seq_boundary = 1.058) {
  use_native <- isTRUE(native) &&
    screen_j_conditional_native_allowed(repeated, seq_p0, seq_boundary)
  native_result <- screen_j_conditional_try_native(
    x = x,
    y = y,
    x_dim = x_dim,
    y_dim = y_dim,
    condition_values = condition_values,
    condition_dims = condition_dims,
    valid = valid,
    exact = exact,
    nsim = nsim,
    seed = seed,
    repeated = repeated,
    random_draw = random_draw,
    seq_limit = seq_limit,
    use_native = use_native
  )
  if (!is.null(native_result)) {
    return(native_result)
  }

  stats <- screen_j_conditional_slice_stats(
    x = x,
    y = y,
    x_dim = x_dim,
    y_dim = y_dim,
    condition_values = condition_values,
    condition_dims = condition_dims,
    valid = valid,
    exact = exact
  )
  if (!isTRUE(stats$has_rows)) {
    return(screen_j_conditional_empty_result(exact, nsim))
  }
  p_chi <- stats$p_chi_asymp
  p_gamma <- stats$p_gamma_asymp
  p_chi_exact <- NA_real_
  p_gamma_exact <- NA_real_
  exact_nsim <- 0L
  if (isTRUE(exact)) {
    exact_nsim <- as.integer(nsim)
    exact_results <- screen_j_conditional_exact_results(
      slices = stats$slices,
      chi_total = stats$chi_square,
      gamma = stats$gamma,
      nsim = exact_nsim,
      seed = seed,
      repeated = repeated,
      random_draw = random_draw,
      seq_limit = seq_limit,
      seq_p0 = seq_p0,
      seq_boundary = seq_boundary,
      use_native = use_native
    )
    p_chi_exact <- exact_results$p_chi
    p_gamma_exact <- exact_results$p_gamma
    exact_nsim <- exact_results$nsim
    p_chi <- p_chi_exact
    p_gamma <- p_gamma_exact
  }
  list(
    chi_square = stats$chi_square,
    df = stats$df,
    p_chi = p_chi,
    gamma = stats$gamma,
    p_gamma = p_gamma,
    p_chi_asymp = stats$p_chi_asymp,
    p_gamma_asymp = stats$p_gamma_asymp,
    p_chi_exact = p_chi_exact,
    p_gamma_exact = p_gamma_exact,
    exact_nsim = exact_nsim,
    ppq = stats$ppq,
    pmq = stats$pmq,
    s = stats$s
  )
}

#' Internal screen j conditional bias test native helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param x Object or value to process.
#' @param y Internal `y` value used by this helper.
#' @param x_dim Internal `x_dim` value used by this helper.
#' @param y_dim Internal `y_dim` value used by this helper.
#' @param condition_values Internal `condition_values` value used by this helper.
#' @param condition_dims Internal `condition_dims` value used by this helper.
#' @param valid Internal `valid` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param repeated Whether to use repeated sequential simulation.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @return The internal `screen_j_conditional_bias_test_native()` computation result.
#' @keywords internal
#' @noRd
screen_j_conditional_bias_test_native <- function(x,
                                                  y,
                                                  x_dim,
                                                  y_dim,
                                                  condition_values,
                                                  condition_dims,
                                                  valid,
                                                  exact,
                                                  nsim,
                                                  seed,
                                                  repeated,
                                                  seq_limit = nsim) {
  raw <- .Call(
    "gRm_screen_j_conditional_bias_test",
    as.integer(x),
    as.integer(y),
    as.integer(x_dim),
    as.integer(y_dim),
    as.matrix(condition_values),
    as.integer(condition_dims),
    as.logical(valid),
    isTRUE(exact),
    as.integer(nsim),
    as.integer(screen_j_source_seed(seed)),
    isTRUE(repeated),
    as.integer(seq_limit[[1L]]),
    PACKAGE = "gRm"
  )
  gamma <- raw$gamma
  p_gamma_asymp <- if (raw$ppq <= 0) {
    2
  } else if (raw$s <= 0) {
    2 * source_tail_norm(4, TRUE)
  } else {
    2 * source_tail_norm(abs(gamma / sqrt(raw$s)), TRUE)
  }
  p_chi_asymp <- if (raw$df > 0L) source_pfchi(raw$df, raw$chi_square) else 1
  list(
    chi_square = raw$chi_square,
    df = raw$df,
    p_chi = raw$p_chi_exact,
    gamma = gamma,
    p_gamma = raw$p_gamma_exact,
    p_chi_asymp = p_chi_asymp,
    p_gamma_asymp = p_gamma_asymp,
    p_chi_exact = raw$p_chi_exact,
    p_gamma_exact = raw$p_gamma_exact,
    exact_nsim = raw$exact_nsim,
    ppq = raw$ppq,
    pmq = raw$pmq,
    s = raw$s
  )
}

#' Internal screen j hypothesis label helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param first Internal `first` value used by this helper.
#' @param second Internal `second` value used by this helper.
#' @param given Internal `given` value used by this helper.
#' @return The internal `screen_j_hypothesis_label()` computation result.
#' @keywords internal
#' @noRd
screen_j_hypothesis_label <- function(first, second, given) {
  paste0(first, "&", second, "|", paste(given, collapse = ""))
}

#' Internal screen j exa p values helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param rows Rows used by the computation.
#' @return The internal `screen_j_exa_p_values()` computation result.
#' @keywords internal
#' @noRd
screen_j_exa_p_values <- function(rows) {
  # Source trace: source/digram_source_20260817/skunits/SKexa1.pas::EXA_SUMMARY1_2 skips
  # hypotheses with RESULTS[hypnr, 2] = 0 ("No tests") before adding p-values
  # to the Benjamini-Hochberg input vector.
  tested <- rows$df > 0L
  p_values <- rows$p_chi[tested]
  if (any(rows$use_gamma[tested])) {
    p_values <- c(p_values, rows$p_gamma[tested & rows$use_gamma])
  }
  p_values
}

#' Internal screen j source stepwise p min helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param test Internal `test` value used by this helper.
#' @param use_gamma Internal `use_gamma` value used by this helper.
#' @param exact Whether to use the exact Monte Carlo branch.
#' @return The internal `screen_j_source_stepwise_p_min()` computation result.
#' @keywords internal
#' @noRd
screen_j_source_stepwise_p_min <- function(test, use_gamma = TRUE, exact = FALSE) {
  # Source trace: SKbias13.StepwiseItemBiasAnalysis and
  # AnalysisOfSpuriousItemBias remove only the candidate with max(p-min) > 0.05.
  # They do not inspect the no-test marker set by SKexa1.EXA_SUMMARY1_2; they
  # read the current RESULTS p-value columns. In the exact reports, no-test
  # hypotheses leave the exact p-value columns at zero and therefore cannot be
  # selected for removal. The asymptotic reports keep the ordinary p-value
  # placeholders and can still remove those candidates, matching DIGRAM output.
  if (isTRUE(exact) && (is.null(test$df) || is.na(test$df) || test$df <= 0L)) {
    return(0)
  }
  if (isTRUE(use_gamma)) {
    min(test$p_chi, test$p_gamma)
  } else {
    test$p_chi
  }
}
