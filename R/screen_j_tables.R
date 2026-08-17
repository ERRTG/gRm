#' Internal screen j pair table helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param x Object or value to process.
#' @param y Internal `y` value used by this helper.
#' @param x_dim Internal `x_dim` value used by this helper.
#' @param y_dim Internal `y_dim` value used by this helper.
#' @param valid Internal `valid` value used by this helper.
#' @return The internal `screen_j_pair_table()` computation result.
#' @keywords internal
#' @noRd
screen_j_pair_table <- function(x, y, x_dim, y_dim, valid) {
  keep <- valid & x >= 1L & x <= x_dim & y >= 1L & y <= y_dim
  tab <- matrix(0, nrow = x_dim, ncol = y_dim)
  if (any(keep)) {
    index <- x[keep] + (y[keep] - 1L) * x_dim
    tab[] <- tabulate(index, nbins = x_dim * y_dim)
  }
  tab
}

#' Internal screen rc chi helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param tab Internal `tab` value used by this helper.
#' @return The internal `screen_rc_chi()` computation result.
#' @keywords internal
#' @noRd
screen_rc_chi <- function(tab) {
  row_totals <- rowSums(tab)
  col_totals <- colSums(tab)
  total <- sum(tab)
  if (total <= 0) {
    return(list(chi_square = 0, df = 0L, p_value = 1))
  }

  expected <- outer(row_totals, col_totals) / total
  positive <- expected > 0
  chi <- sum((tab[positive] - expected[positive])^2 / expected[positive])
  df <- (sum(row_totals > 0) - 1L) * (sum(col_totals > 0) - 1L)
  if (df <= 0L) {
    list(chi_square = chi, df = 0L, p_value = 1)
  } else {
    # Source trace: SkStat.RCCHI/SourceRaschCore.SourcePFCHI computes the upper
    # chi-square tail for the Pearson statistic.
    list(chi_square = chi, df = df, p_value = source_pfchi(df, chi))
  }
}

#' Internal screen rc chi square helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param tab Internal `tab` value used by this helper.
#' @return The internal `screen_rc_chi_square()` computation result.
#' @keywords internal
#' @noRd
screen_rc_chi_square <- function(tab) {
  screen_rc_chi(tab)$chi_square
}

#' Internal screen rc gamma helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param tab Internal `tab` value used by this helper.
#' @return The internal `screen_rc_gamma()` computation result.
#' @keywords internal
#' @noRd
screen_rc_gamma <- function(tab) {
  stats <- source_rc_gamma_stats(tab, include_cells = TRUE)
  ppq <- stats$ppq
  pmq <- stats$pmq
  gamma_value <- stats$gamma
  n <- sum(tab)
  if (ppq <= 0) {
    return(list(gamma = 0, ppq = ppq, pmq = pmq, s = 0, p_value = 1, success = FALSE))
  }

  # Source trace: SkStat.RCGAMMA/SourceRaschCore.SourceRCGammaStats computes
  # S = 4 * (-PMQ * PMQ / N + sum(n_ij * (AIJ - DIJ)^2)) and tests
  # abs(gamma / (sqrt(S) / PPQ)) against the upper normal tail.
  s <- if (n > 0) -pmq * (pmq / n) else 0
  m <- stats$aij - stats$dij
  s <- s + sum(tab * m * m)
  s <- 4 * s
  p_value <- if (s > 0) {
    source_tail_norm(abs(gamma_value / (sqrt(s) / ppq)), TRUE)
  } else {
    1
  }

  list(gamma = gamma_value, ppq = ppq, pmq = pmq, s = s, p_value = p_value, success = TRUE)
}

#' Internal screen rc gamma counts helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param tab Internal `tab` value used by this helper.
#' @return The internal `screen_rc_gamma_counts()` computation result.
#' @keywords internal
#' @noRd
screen_rc_gamma_counts <- function(tab) {
  source_rc_gamma_counts(tab)
}

#' Internal screen j source seed helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param seed Random-stream seed.
#' @return The internal `screen_j_source_seed()` computation result.
#' @keywords internal
#' @noRd
screen_j_source_seed <- function(seed) {
  seed <- as.integer(seed[[1L]])
  if (is.na(seed)) {
    return(9L)
  }
  max(0L, min(255L, seed))
}

#' Internal screen j source random stream helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param seed Random-stream seed.
#' @return The internal `screen_j_source_random_stream()` computation result.
#' @keywords internal
#' @noRd
screen_j_source_random_stream <- function(seed) {
  state <- as.numeric(screen_j_source_seed(seed))
  base <- 65536
  modulus <- base * base
  multiplier_hi <- 2056
  multiplier_lo <- 33797

  function() {
    # Delphi Random advances RandSeed as a 32-bit LCG. The 16-bit limb update
    # keeps the low bits exact in R's double-only arithmetic.
    state_lo <- state %% base
    state_hi <- floor(state / base)
    low_product <- multiplier_lo * state_lo + 1
    next_lo <- low_product %% base
    carry <- floor(low_product / base)
    next_hi <- (multiplier_hi * state_lo + multiplier_lo * state_hi + carry) %% base
    state <<- next_hi * base + next_lo
    state / modulus
  }
}

#' Internal screen j source seed after draws helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param seed Random-stream seed.
#' @param draw_count Internal `draw_count` value used by this helper.
#' @return The internal `screen_j_source_seed_after_draws()` computation result.
#' @keywords internal
#' @noRd
screen_j_source_seed_after_draws <- function(seed, draw_count) {
  state <- as.numeric(screen_j_source_seed(seed))
  base <- 65536
  multiplier_hi <- 2056
  multiplier_lo <- 33797
  draw_count <- as.integer(draw_count)
  if (draw_count <= 0L) {
    return(state)
  }
  for (draw in seq_len(draw_count)) {
    state_lo <- state %% base
    state_hi <- floor(state / base)
    low_product <- multiplier_lo * state_lo + 1
    next_lo <- low_product %% base
    carry <- floor(low_product / base)
    next_hi <- (multiplier_hi * state_lo + multiplier_lo * state_hi + carry) %% base
    state <- next_hi * base + next_lo
  }
  state
}

#' Internal screen j repeated seq limit helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param project Encoded gRm project.
#' @param repeated Whether to use repeated sequential simulation.
#' @param nsim Requested simulation count.
#' @param exact_state Internal `exact_state` value used by this helper.
#' @return The internal `screen_j_repeated_seq_limit()` computation result.
#' @keywords internal
#' @noRd
screen_j_repeated_seq_limit <- function(project, repeated, nsim, exact_state = NULL) {
  nsim <- as.integer(nsim[[1L]])
  if (!is.null(exact_state)) {
    return(as.integer(exact_state$seq_limit[[1L]]))
  }
  if (!isTRUE(repeated)) {
    return(nsim)
  }
  as.integer(gRm_exact_command_state("repeated", nsim = nsim)$seq_limit)
}

#' Internal screen j strata table helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param x Object or value to process.
#' @param y Internal `y` value used by this helper.
#' @param z Internal `z` value used by this helper.
#' @param x_dim Internal `x_dim` value used by this helper.
#' @param y_dim Internal `y_dim` value used by this helper.
#' @param z_dim Internal `z_dim` value used by this helper.
#' @param valid Internal `valid` value used by this helper.
#' @return The internal `screen_j_strata_table()` computation result.
#' @keywords internal
#' @noRd
screen_j_strata_table <- function(x, y, z, x_dim, y_dim, z_dim, valid) {
  keep <- valid & x >= 1L & x <= x_dim & y >= 1L & y <= y_dim & z >= 1L & z <= z_dim
  tab <- array(0, dim = c(x_dim, y_dim, z_dim))
  if (any(keep)) {
    index <- x[keep] + (y[keep] - 1L) * x_dim + (z[keep] - 1L) * x_dim * y_dim
    tab[] <- tabulate(index, nbins = x_dim * y_dim * z_dim)
  }
  tab
}

#' Internal screen j strata slices helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param strata Internal `strata` value used by this helper.
#' @return The internal `screen_j_strata_slices()` computation result.
#' @keywords internal
#' @noRd
screen_j_strata_slices <- function(strata) {
  lapply(seq_len(dim(strata)[[3L]]), function(level) {
    strata[, , level, drop = FALSE][, , 1L]
  })
}

#' Internal screen j source informative slice helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param slice Internal `slice` value used by this helper.
#' @return The internal `screen_j_source_informative_slice()` computation result.
#' @keywords internal
#' @noRd
screen_j_source_informative_slice <- function(slice) {
  sum(rowSums(slice) > 0L) >= 2L && sum(colSums(slice) > 0L) >= 2L
}

#' Internal screen j partial gamma helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param strata Internal `strata` value used by this helper.
#' @return The internal `screen_j_partial_gamma()` computation result.
#' @keywords internal
#' @noRd
screen_j_partial_gamma <- function(strata) {
  ppq_total <- 0
  pmq_total <- 0
  s_total <- 0
  for (level in seq_len(dim(strata)[[3L]])) {
    slice <- strata[, , level, drop = FALSE][, , 1L]
    stats <- screen_rc_gamma(slice)
    ppq_total <- ppq_total + stats$ppq
    pmq_total <- pmq_total + stats$pmq
    s_total <- s_total + stats$s
  }
  if (ppq_total <= 0) {
    return(list(gamma = 0, p_value = 2, ppq = 0, pmq = 0, s = 0))
  }
  gamma <- pmq_total / ppq_total
  s_total <- s_total / ppq_total
  s_total <- s_total / ppq_total
  u <- if (s_total <= 0) 4 else abs(gamma / sqrt(s_total))
  p <- 2 * if (u > 4) 0.000001 else source_tail_norm(u, TRUE)
  list(gamma = gamma, p_value = p, ppq = ppq_total, pmq = pmq_total, s = s_total)
}

#' Internal screen j source single helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param value Value to validate or transform.
#' @return The internal `screen_j_source_single()` computation result.
#' @keywords internal
#' @noRd
screen_j_source_single <- function(value) {
  storage.mode(value) <- "double"
  readBin(writeBin(value, raw(), size = 4L), "numeric", n = length(value), size = 4L)
}

#' Internal screen j partial chi helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param strata Internal `strata` value used by this helper.
#' @return The internal `screen_j_partial_chi()` computation result.
#' @keywords internal
#' @noRd
screen_j_partial_chi <- function(strata) {
  chi_total <- 0
  df_total <- 0L
  for (level in seq_len(dim(strata)[[3L]])) {
    slice <- strata[, , level, drop = FALSE][, , 1L]
    chi <- screen_j_rc_chi_source_expected(slice)
    chi_total <- chi_total + chi$chi_square
    df_total <- df_total + chi$df
  }
  if (df_total <= 0L) {
    list(stat = 0, df = 0L, p_value = 1)
  } else {
    list(stat = chi_total, df = df_total, p_value = source_pfchi(df_total, chi_total))
  }
}

#' Internal screen j rc chi source expected helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param tab Internal `tab` value used by this helper.
#' @return The internal `screen_j_rc_chi_source_expected()` computation result.
#' @keywords internal
#' @noRd
screen_j_rc_chi_source_expected <- function(tab) {
  chi_square <- screen_j_rc_chi_square_source_expected(tab)
  row_totals <- rowSums(tab)
  col_totals <- colSums(tab)
  df <- (sum(row_totals > 0) - 1L) * (sum(col_totals > 0) - 1L)
  if (sum(tab) <= 0 || df <= 0L) {
    list(chi_square = chi_square, df = 0L, p_value = 1)
  } else {
    list(chi_square = chi_square, df = df, p_value = source_pfchi(df, chi_square))
  }
}

#' Internal screen j rc chi square source expected helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/PAS_skunits/SKbias7.pas::Item_Screening`.
#' @param tab Internal `tab` value used by this helper.
#' @return The internal `screen_j_rc_chi_square_source_expected()` computation result.
#' @keywords internal
#' @noRd
screen_j_rc_chi_square_source_expected <- function(tab) {
  row_totals <- rowSums(tab)
  col_totals <- colSums(tab)
  total <- sum(tab)
  if (total <= 0) {
    return(0)
  }

  chi_square <- 0
  # Source trace: SkStat.RCCHI loops the first Pascal table index before the
  # second (`I := 1..C`, then `J := 1..R`). R stores those indices as matrix
  # rows and columns, so retain row-outer/column-inner accumulation order.
  for (row in seq_len(nrow(tab))) {
    for (col in seq_len(ncol(tab))) {
      col_share <- col_totals[[col]] / total
      # Source trace: source/PAS_skunits/SKxyz1.PAS::MAKE_XYZ_TABLE and
      # source/PAS_skunits/SKbigtab.pas::Transfer_BT_to_XYZ_TABLE store
      # expected cells as row_margin * (col_margin / total). SKrandom.GENTAB1
      # passes that stored table to SkStat.RCCHI for exact chi-square
      # comparisons.
      expected <- row_totals[[row]] * col_share
      if (expected > 0) {
        residual <- tab[row, col] - expected
        chi_square <- chi_square + residual * (residual / expected)
      }
    }
  }
  chi_square
}
