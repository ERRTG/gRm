# Derive DIGRAM SCREEN J item-screening values
#
# Computes the source-shaped values printed by DIGRAM's `SCREEN J` command.
# This is the item-screening branch from `DIGRAM1f.pas` command 3 and
# `SKbias7.pas::Item_screening`; it is distinct from the graphical
# `SCREEN` command, which is outside the current package scope.
#
# The validation runtime target uses the asymptotic, two-sided branch by default.
# Exact mode follows the source `GENTAB1` Monte Carlo branch and computes
# partial p-values directly from the project data.
#
# @param project A parsed DIGRAM project.
# @param exact Logical; when `TRUE`, replace the partial item-item and partial
#   DIF p-values with the source-shaped Monte Carlo exact p-values used by
#   `SKbias3.XYZ_bias_ANALYSE`.
# @param nsim Number of random tables for exact p-values.
# @param seed Random seed for the source-shaped exact branch.
# @return A `gRm_screen_j_values` object.
screen_j_values <- function(project,
                            exact = FALSE,
                            repeated = FALSE,
                            nsim = 1000L,
                            seed = 9L,
                            exact_state = NULL) {
  exact <- isTRUE(exact)
  repeated <- isTRUE(repeated)
  nsim <- as.integer(nsim)
  seed <- as.integer(seed)
  if (is.null(exact_state)) {
    exact_state <- gRm_exact_state_from_flags(exact, repeated, nsim = nsim, seed = seed)
  }
  exact <- isTRUE(exact_state$exact)
  repeated <- isTRUE(exact_state$sequential)
  nsim <- as.integer(exact_state$nsim)
  seed <- as.integer(exact_state$seed)
  if (exact && (length(nsim) != 1L || is.na(nsim) || nsim < 1L)) {
    stop("`nsim` must be a positive integer for exact SCREEN J.", call. = FALSE)
  }
  if (exact && (length(seed) != 1L || is.na(seed))) {
    stop("`seed` must be a single integer for exact SCREEN J.", call. = FALSE)
  }
  seq_limit <- screen_j_repeated_seq_limit(project, repeated, nsim, exact_state = exact_state)
  seq_p0 <- if (isTRUE(repeated)) as.numeric(exact_state$seq_p0[[1L]]) else 0.05
  seq_boundary <- if (isTRUE(repeated)) as.numeric(exact_state$seq_b[[1L]]) else 1.058
  old_seed <- NULL
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (exact) {
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  items <- project$items
  backgrounds <- project$backgrounds
  n_items <- nrow(items)
  n_exo <- nrow(backgrounds)
  largest_possible_score <- sum(items$raw_max - 1L)

  item_matrix <- matrix(0L, nrow = nrow(project$raw_data), ncol = n_items)
  for (item_index in seq_len(n_items)) {
    item_matrix[, item_index] <- as.integer(project$raw_data[, items$position[[item_index]]])
  }
  complete_items <- rep(TRUE, nrow(item_matrix))
  for (item_index in seq_len(n_items)) {
    complete_items <- complete_items &
      item_matrix[, item_index] >= 1L &
      item_matrix[, item_index] <= items$raw_max[[item_index]]
  }
  item_score <- rowSums(sweep(item_matrix, 2L, 1L, "-"))
  max_score <- if (any(complete_items)) max(item_score[complete_items]) else largest_possible_score
  good_items <- complete_items & item_score >= 1L &
    item_score <= max_score &
    item_score <= (largest_possible_score - 1L)

  marginal_gamma <- matrix(0, nrow = n_items, ncol = n_items)
  marginal_p <- matrix(2, nrow = n_items, ncol = n_items)
  partial_gamma <- matrix(0, nrow = n_items, ncol = n_items)
  partial_p <- matrix(2, nrow = n_items, ncol = n_items)
  partial_ppq <- matrix(0, nrow = n_items, ncol = n_items)
  partial_pmq <- matrix(0, nrow = n_items, ncol = n_items)
  rest_gamma <- rep(0, n_items)
  rest_p <- rep(2, n_items)

  dimnames(marginal_gamma) <- dimnames(marginal_p) <-
    dimnames(partial_gamma) <- dimnames(partial_p) <-
    dimnames(partial_ppq) <- dimnames(partial_pmq) <-
    list(items$label_code, items$label_code)
  names(rest_gamma) <- names(rest_p) <- items$label_code

  for (row_item in seq_len(n_items)) {
    # Source trace: SKbias7.inexpensive_rosenberg sets THE_ITEMS(.X.) := 1
    # before calling SCORE(THE_ITEMS), so the score column is Score \ X, not
    # Score minus the raw item code.
    rest_score <- item_score - (item_matrix[, row_item] - 1L)
    rest_table <- screen_j_pair_table(
      item_matrix[, row_item],
      rest_score + 1L,
      items$raw_max[[row_item]],
      max_score + 1L,
      complete_items
    )
    rest_stats <- screen_rc_gamma(rest_table)
    rest_gamma[[row_item]] <- rest_stats$gamma
    rest_p[[row_item]] <- rest_stats$p_value

    for (col_item in seq_len(n_items)) {
      if (row_item == col_item) {
        next
      }

      marginal_table <- screen_j_pair_table(
        item_matrix[, row_item],
        item_matrix[, col_item],
        items$raw_max[[row_item]],
        items$raw_max[[col_item]],
        complete_items
      )
      marginal_stats <- screen_rc_gamma(marginal_table)
      marginal_gamma[row_item, col_item] <- marginal_stats$gamma
      marginal_p[row_item, col_item] <- marginal_stats$p_value

      if (row_item > col_item) {
        x_values <- item_matrix[, col_item]
        y_values <- item_matrix[, row_item]
        x_dim <- items$raw_max[[col_item]]
        y_dim <- items$raw_max[[row_item]]
        partial_valid <- good_items
      } else {
        x_values <- item_matrix[, row_item]
        y_values <- item_matrix[, col_item]
        x_dim <- items$raw_max[[row_item]]
        y_dim <- items$raw_max[[col_item]]
        partial_valid <- complete_items
      }
      strata <- screen_j_strata_table(
        x_values,
        y_values,
        rest_score + 1L,
        x_dim,
        y_dim,
        max_score + 1L,
        partial_valid
      )
      partial_stats <- screen_j_partial_gamma(strata)
      partial_gamma[row_item, col_item] <- partial_stats$gamma
      partial_p[row_item, col_item] <- if (exact) {
        native_pair <- if (screen_j_item_pair_native_allowed(repeated, seq_p0, seq_boundary)) {
          screen_j_item_pair_conditional_exact_native(
            x_values,
            y_values,
            x_dim,
            y_dim,
            rest_score + 1L,
            max_score + 1L,
            partial_valid,
            nsim,
            seed,
            sequential = repeated,
            seq_limit = seq_limit
          )
        } else {
          NULL
        }
        if (!is.null(native_pair)) {
          native_pair$p_gamma_exact
        } else {
          partial_chi <- screen_j_partial_chi(strata)
          screen_j_exact_chi_gamma_prepared_r(
            screen_j_prepare_exact_slices(screen_j_strata_slices(strata)),
            partial_chi$stat,
            partial_stats$gamma,
            nsim,
            seed,
            sequential = repeated,
            seq_limit = seq_limit,
            seq_p0 = seq_p0,
            seq_boundary = seq_boundary
          )$p_gamma
        }
      } else {
        partial_stats$p_value
      }
      partial_ppq[row_item, col_item] <- partial_stats$ppq
      partial_pmq[row_item, col_item] <- partial_stats$pmq
    }
  }

  attr(project, "screen_j_exact") <- exact
  attr(project, "screen_j_nsim") <- nsim
  attr(project, "screen_j_seed") <- seed
  attr(project, "screen_j_repeated") <- repeated
  attr(project, "screen_j_seq_limit") <- seq_limit
  attr(project, "screen_j_seq_p0") <- seq_p0
  attr(project, "screen_j_seq_boundary") <- seq_boundary
  attr(project, "screen_j_exact_state") <- exact_state

  marginal_exo <- screen_j_exo_values(
    project = project,
    item_matrix = item_matrix,
    item_score = item_score,
    complete_items = complete_items,
    items = items,
    backgrounds = backgrounds,
    max_score = max_score,
    partial = FALSE
  )
  partial_exo <- screen_j_exo_values(
    project = project,
    item_matrix = item_matrix,
    item_score = item_score,
    complete_items = complete_items,
    items = items,
    backgrounds = backgrounds,
    max_score = max_score,
    partial = TRUE
  )

  bh_item_values <- c(as.vector(partial_p[row(partial_p) != col(partial_p)]))
  bh_values <- bh_item_values
  if (n_exo > 0L) {
    bh_values <- c(bh_values, as.vector(partial_exo$p))
  }
  weighted_gamma <- screen_j_weighted_partial_gamma(partial_ppq, partial_pmq)
  initial_item_bias <- partial_exo$p <= source_bh_critical(bh_values, 0.05)
  post_screen_dif <- screen_j_post_screen_dif(
    item_matrix = item_matrix,
    exo_values = partial_exo$exo_values,
    item_score = item_score,
    complete_items = complete_items,
    items = items,
    backgrounds = backgrounds,
    # Source trace: SKbias13.AnalysisOfSpuriousItemBias and
    # StepwiseItemBiasAnalysis dimension their score variable by
    # Largest_possible_score, not by the largest observed score.
    max_score = largest_possible_score,
    partial_exo = partial_exo,
    initial_item_bias = initial_item_bias,
    fdr_01 = source_bh_critical(bh_values, 0.01),
    exact = exact,
    nsim = nsim,
    seed = seed,
    repeated = repeated,
    seq_limit = seq_limit,
    seq_p0 = seq_p0,
    seq_boundary = seq_boundary
  )
  local_dependence <- screen_j_stepwise_local_dependence(
    partial_gamma,
    partial_p,
    weighted_gamma,
    source_bh_critical(bh_values, 0.05)
  )
  score_effects <- screen_j_score_effect_values(
    project = project,
    item_matrix = item_matrix,
    item_score = item_score,
    complete_items = complete_items,
    backgrounds = backgrounds,
    max_score = max_score,
    exact = exact,
    nsim = nsim,
    seed = seed,
    repeated = repeated,
    seq_limit = seq_limit,
    seq_p0 = seq_p0,
    seq_boundary = seq_boundary
  )

  values <- structure(
    list(
      items = items,
      backgrounds = backgrounds,
      marginal = list(
        item_gamma = marginal_gamma,
        item_p = marginal_p,
        rest_gamma = rest_gamma,
        rest_p = rest_p,
        exo_stat = marginal_exo$stat,
        exo_p = marginal_exo$p,
        exo_kind = marginal_exo$kind,
        score_stat = marginal_exo$score_stat,
        score_p = marginal_exo$score_p
      ),
      partial = list(
        item_gamma = partial_gamma,
        item_p = partial_p,
        item_ppq = partial_ppq,
        item_pmq = partial_pmq,
        weighted_gamma = weighted_gamma,
        # Source trace: SKbias7.Item_screening prints
        # sum(abs(part_g[i,j] + part_g[j,i])) / (nitems * (nitems - 1)).
        average_abs_gamma = screen_j_average_abs_partial_gamma(partial_gamma),
        exo_stat = partial_exo$stat,
        exo_p = partial_exo$p,
        exo_ppq = partial_exo$ppq,
        exo_pmq = partial_exo$pmq,
        exo_kind = partial_exo$kind
      ),
      model = list(
        local_dependence = local_dependence,
        initial_item_bias = initial_item_bias,
        item_bias = post_screen_dif$item_bias,
        item_bias_status = post_screen_dif$item_bias_status,
        spurious_dif = post_screen_dif$spurious_dif,
        multiple_dif = post_screen_dif$multiple_dif,
        score_effects = score_effects
      ),
      bh = list(
        fdr_05 = source_bh_critical(bh_values, 0.05),
        fdr_01 = source_bh_critical(bh_values, 0.01),
        fdr_001 = source_bh_critical(bh_values, 0.001),
        n_tests = length(bh_values)
      ),
      source_status = c(
        command = "DIGRAM1f_command_3_SCREEN_parameter_J",
        values = if (exact) {
          "SKbias7_Item_screening_inexpensive_rosenberg_inexpensive_itembias1_exact_GENTAB1"
        } else {
          "SKbias7_Item_screening_inexpensive_rosenberg_inexpensive_itembias1_asymptotic"
        },
        rendering = "SKbias7_Item_screening_WRITE_MESSAGE_cases_1_2_3"
      ),
      exact = exact,
      nsim = if (exact) nsim else 0L,
      seed = if (exact) seed else NA_integer_
    ),
    class = "gRm_screen_j_values"
  )
  values
}

screen_j_pair_table <- function(x, y, x_dim, y_dim, valid) {
  keep <- valid & x >= 1L & x <= x_dim & y >= 1L & y <= y_dim
  tab <- matrix(0, nrow = x_dim, ncol = y_dim)
  if (any(keep)) {
    index <- x[keep] + (y[keep] - 1L) * x_dim
    tab[] <- tabulate(index, nbins = x_dim * y_dim)
  }
  tab
}

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

screen_rc_chi_square <- function(tab) {
  screen_rc_chi(tab)$chi_square
}

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

screen_rc_gamma_counts <- function(tab) {
  source_rc_gamma_counts(tab)
}

screen_j_source_seed <- function(seed) {
  seed <- as.integer(seed[[1L]])
  if (is.na(seed)) {
    return(9L)
  }
  max(0L, min(255L, seed))
}

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

screen_j_strata_table <- function(x, y, z, x_dim, y_dim, z_dim, valid) {
  keep <- valid & x >= 1L & x <= x_dim & y >= 1L & y <= y_dim & z >= 1L & z <= z_dim
  tab <- array(0, dim = c(x_dim, y_dim, z_dim))
  if (any(keep)) {
    index <- x[keep] + (y[keep] - 1L) * x_dim + (z[keep] - 1L) * x_dim * y_dim
    tab[] <- tabulate(index, nbins = x_dim * y_dim * z_dim)
  }
  tab
}

screen_j_strata_slices <- function(strata) {
  lapply(seq_len(dim(strata)[[3L]]), function(level) {
    strata[, , level, drop = FALSE][, , 1L]
  })
}

screen_j_source_informative_slice <- function(slice) {
  sum(rowSums(slice) > 0L) >= 2L && sum(colSums(slice) > 0L) >= 2L
}

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

screen_j_source_single <- function(value) {
  storage.mode(value) <- "double"
  readBin(writeBin(value, raw(), size = 4L), "numeric", n = length(value), size = 4L)
}

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
    list(stat = 0, p_value = 1)
  } else {
    list(stat = chi_total, p_value = source_pfchi(df_total, chi_total))
  }
}

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

screen_j_rc_chi_square_source_expected <- function(tab) {
  row_totals <- rowSums(tab)
  col_totals <- colSums(tab)
  total <- sum(tab)
  if (total <= 0) {
    return(0)
  }

  chi_square <- 0
  for (col in seq_len(ncol(tab))) {
    col_share <- col_totals[[col]] / total
    for (row in seq_len(nrow(tab))) {
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

# Source-shaped Monte Carlo exact partial gamma p-value
#
# Source trace: `SKbias3.XYZ_bias_ANALYSE` calls `SKrandom.GENTAB1` for each
# score-conditioned slice, adds simulated `PPQ` and `PMQ` over slices, and
# counts `abs(simulated_gamma) >= abs(observed_gamma)`.
#
# @param strata Three-way table with dimensions item/background by item/exo by
#   conditioning score.
# @param observed_gamma Observed partial gamma.
# @param nsim Number of simulated tables.
# @return Two-sided Monte Carlo exact p-value.
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

screen_j_exact_gamma_slices <- function(slices, observed_gamma, nsim, seed = NULL) {
  native <- screen_j_exact_gamma_slices_native(slices, observed_gamma, nsim, seed)
  if (!is.null(native)) {
    return(as.numeric(native[[1L]]))
  }
  screen_j_exact_gamma_prepared(screen_j_prepare_exact_slices(slices), observed_gamma, nsim, seed)
}

screen_j_exact_gamma_prepared <- function(prepared_slices, observed_gamma, nsim, seed = NULL) {
  screen_j_exact_gamma_prepared_r(prepared_slices, observed_gamma, nsim, seed)
}

screen_j_exact_gamma_prepared_r <- function(prepared_slices, observed_gamma, nsim, seed = NULL) {
  random_draw <- screen_j_source_random_stream(seed)
  exceed <- 0L
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

# Source-shaped Monte Carlo exact partial chi-square p-value
#
# Source trace: `SKbias3.XYZ_bias_ANALYSE` adds simulated chi-square statistics
# over score-conditioned slices and counts `simulated_chi >= observed_chi`.
# When the DIGRAM exact command state is sequential/repeated, the same
# `SEQUENTIAL`, `SEQ_P0`, `SEQ_B`, and `seq_limit` state is used while
# evaluating simulated chi results.
#
# @param strata Three-way table.
# @param observed_chi Observed partial chi-square.
# @param nsim Number of simulated tables.
# @return Monte Carlo exact p-value.
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
  for (sim in seq_len(nsim)) {
    chi_total <- 0
    for (prepared in prepared_slices) {
      generated <- exo_select_gentab1_prepared(prepared, random_draw = random_draw)
      chi_total <- chi_total + screen_j_rc_chi_square_prepared_expected(generated, prepared)
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
    return(list(p_chi = native$p_chi, p_gamma = native$p_gamma, nsim = native$nsim))
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

screen_j_exact_chi_gamma_prepared_r <- function(prepared_slices,
                                                observed_chi,
                                                observed_gamma,
                                                nsim,
                                                seed = NULL,
                                                sequential = FALSE,
                                                random_draw = NULL,
                                                seq_limit = nsim,
                                                seq_p0 = 0.05,
                                                seq_boundary = 1.058) {
  if (is.null(random_draw)) {
    random_draw <- screen_j_source_random_stream(seed)
  }
  chi_exceed <- 0L
  gamma_exceed <- 0L
  chi_status <- FALSE
  gamma_status <- FALSE
  seq_p0 <- as.numeric(seq_p0[[1L]])
  seq_boundary <- as.numeric(seq_boundary[[1L]])
  seq_limit <- as.integer(seq_limit[[1L]])
  for (sim in seq_len(nsim)) {
    chi_total <- 0
    ppq_total <- 0
    pmq_total <- 0
    for (prepared in prepared_slices) {
      generated <- exo_select_gentab1_prepared(prepared, random_draw = random_draw)
      chi_total <- chi_total + screen_j_rc_chi_square_prepared_expected(generated, prepared)
      gamma_counts <- screen_rc_gamma_counts(generated)
      ppq_total <- ppq_total + gamma_counts$ppq
      pmq_total <- pmq_total + gamma_counts$pmq
    }
    if (chi_total >= observed_chi) {
      chi_exceed <- chi_exceed + 1L
    }
    simulated_gamma <- if (ppq_total > 0) pmq_total / ppq_total else 0
    if (abs(simulated_gamma) >= abs(observed_gamma)) {
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
    if (isTRUE(sequential) && chi_status && gamma_status) {
      break
    }
  }
  list(p_chi = chi_exceed / sim, p_gamma = gamma_exceed / sim, nsim = sim)
}

screen_j_seq_t <- function(exceed, sim, p0) {
  if (sim < 21L) {
    return(0)
  }
  root <- sqrt(sim)
  exceed / root - root * p0
}

screen_j_exact_native_available <- function() {
  disabled <- tolower(Sys.getenv("RDIGRAM_SCREEN_J_EXACT_CPP", unset = "")) %in%
    c("0", "false", "no", "off")
  if (disabled) {
    return(FALSE)
  }
  tryCatch(
    is.loaded("gRm_screen_j_exact_chi_gamma_slices", PACKAGE = "gRm"),
    error = function(e) FALSE
  )
}

screen_j_conditional_native_allowed <- function(repeated, seq_p0, seq_boundary) {
  screen_j_conditional_native_controls_allowed(repeated, seq_p0, seq_boundary) &&
    screen_j_conditional_native_source_faithful()
}

screen_j_conditional_native_controls_allowed <- function(repeated, seq_p0, seq_boundary) {
  !isTRUE(repeated) ||
    (abs(as.numeric(seq_p0[[1L]]) - 0.05) < 1e-15 &&
      abs(as.numeric(seq_boundary[[1L]]) - 1.058) < 1e-15)
}

screen_j_conditional_native_probe_fixtures <- function() {
  # Source trace: the first fixture includes condition level 2 with only one
  # nonempty x category, so parity requires the native path to follow the
  # SKxyz1.MAKE_XYZ_TABLE/source R informative-slice filter before GENTAB1.
  list(
    list(
      args = list(
        x = c(1L, 2L, 1L, 2L, 3L, 3L, 1L, 1L, 2L, 3L, 2L, 3L),
        y = c(1L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L),
        x_dim = 3L,
        y_dim = 2L,
        condition_values = matrix(c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 3L, 3L, 3L, 3L), ncol = 1L),
        condition_dims = 3L,
        valid = rep(TRUE, 12L)
      ),
      fixed = list(nsim = 200L, seed = 123L, seq_limit = 200L),
      repeated = list(nsim = 200L, seed = 123L, seq_limit = 200L)
    ),
    list(
      # Regression fixture from the optimized parity tests, retained as parity
      # coverage for the optimized native route after the p_gamma mismatch fix.
      args = list(
        x = c(1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
        y = c(1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 1L),
        x_dim = 3L,
        y_dim = 2L,
        condition_values = matrix(c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L), ncol = 1L),
        condition_dims = 4L,
        valid = rep(TRUE, 12L)
      ),
      fixed = list(nsim = 37L, seed = 9L, seq_limit = 37L),
      repeated = list(nsim = 80L, seed = 17L, seq_limit = 3L)
    )
  )
}

screen_j_conditional_native_probe_matches <- function(fixture, controls, repeated) {
  controls <- list(
    exact = TRUE,
    nsim = as.integer(controls$nsim),
    seed = as.integer(controls$seed),
    repeated = isTRUE(repeated),
    seq_limit = as.integer(controls$seq_limit)
  )
  reference <- do.call(
    screen_j_conditional_bias_test,
    c(fixture$args, controls, list(native = FALSE, seq_p0 = 0.05, seq_boundary = 1.058))
  )
  native <- do.call(screen_j_conditional_bias_test_native, c(fixture$args, controls))
  fields <- c(
    "chi_square", "df", "gamma", "p_chi", "p_gamma",
    "p_chi_exact", "p_gamma_exact", "exact_nsim"
  )
  all(vapply(fields, function(field) {
    isTRUE(all.equal(reference[[field]], native[[field]], tolerance = 1e-12, check.attributes = FALSE))
  }, logical(1L)))
}

screen_j_conditional_native_source_faithful <- local({
  cached <- NULL
  function() {
    if (!screen_j_exact_native_available()) {
      return(FALSE)
    }
    if (!is.loaded("gRm_screen_j_conditional_bias_test", PACKAGE = "gRm")) {
      return(FALSE)
    }
    if (!is.null(cached)) {
      return(cached)
    }
    cached <<- tryCatch({
      fixtures <- screen_j_conditional_native_probe_fixtures()
      all(vapply(fixtures, function(fixture) {
        screen_j_conditional_native_probe_matches(fixture, fixture$fixed, repeated = FALSE) &&
          screen_j_conditional_native_probe_matches(fixture, fixture$repeated, repeated = TRUE)
      }, logical(1L)))
    }, error = function(e) FALSE)
    cached
  }
})

screen_j_item_pair_native_allowed <- function(repeated, seq_p0, seq_boundary) {
  screen_j_conditional_native_controls_allowed(repeated, seq_p0, seq_boundary) &&
    screen_j_item_pair_native_source_faithful()
}

screen_j_item_pair_native_source_faithful <- local({
  cached <- NULL
  function() {
    if (!screen_j_exact_native_available()) {
      return(FALSE)
    }
    if (!is.null(cached)) {
      return(cached)
    }
    cached <<- tryCatch({
      available <- is.loaded("gRm_screen_j_item_pair_conditional_exact", PACKAGE = "gRm")
      if (!available) {
        return(FALSE)
      }
      item_matrix <- matrix(
        c(
          1L, 1L, 2L,
          1L, 2L, 1L,
          2L, 1L, 2L,
          2L, 2L, 1L,
          3L, 1L, 2L,
          3L, 2L, 1L,
          1L, 2L, 2L,
          2L, 1L, 1L,
          3L, 2L, 2L
        ),
        ncol = 3L,
        byrow = TRUE
      )
      item_score <- rowSums(item_matrix - 1L)
      rest_score <- item_score - (item_matrix[, 1L] - 1L)
      valid <- rep(TRUE, nrow(item_matrix))
      strata <- screen_j_strata_table(
        item_matrix[, 1L],
        item_matrix[, 2L],
        rest_score + 1L,
        3L,
        2L,
        max(rest_score) + 1L,
        valid
      )
      partial_stats <- screen_j_partial_gamma(strata)
      partial_chi <- screen_j_partial_chi(strata)
      reference <- screen_j_exact_chi_gamma_prepared_r(
        screen_j_prepare_exact_slices(screen_j_strata_slices(strata)),
        partial_chi$stat,
        partial_stats$gamma,
        41L,
        9L,
        sequential = FALSE,
        seq_limit = 41L
      )
      native <- screen_j_item_pair_conditional_exact_native(
        item_matrix[, 1L],
        item_matrix[, 2L],
        3L,
        2L,
        rest_score + 1L,
        max(rest_score) + 1L,
        valid,
        nsim = 41L,
        seed = 9L,
        sequential = FALSE,
        seq_limit = 41L
      )
      !is.null(native) &&
        identical(native$exact_nsim, reference$nsim) &&
        isTRUE(all.equal(native$p_gamma_exact, reference$p_gamma, tolerance = 1e-7, check.attributes = FALSE))
    }, error = function(e) FALSE)
    cached
  }
})

screen_j_exact_chi_gamma_slices_native <- function(slices,
                                                   observed_chi,
                                                   observed_gamma,
                                                   nsim,
                                                   seed = NULL,
                                                   sequential = FALSE,
                                                   seq_limit = nsim) {
  if (!screen_j_exact_native_available()) {
    return(NULL)
  }
  .Call(
    "gRm_screen_j_exact_chi_gamma_slices",
    slices,
    as.numeric(observed_chi),
    as.numeric(observed_gamma),
    as.integer(nsim),
    as.integer(screen_j_source_seed(seed)),
    isTRUE(sequential),
    as.integer(seq_limit[[1L]]),
    PACKAGE = "gRm"
  )
}

screen_j_exact_chi_slices_native <- function(slices,
                                             observed_chi,
                                             nsim,
                                             seed = NULL,
                                             sequential = FALSE,
                                             seq_limit = nsim) {
  if (!screen_j_exact_native_available()) {
    return(NULL)
  }
  .Call(
    "gRm_screen_j_exact_chi_slices",
    slices,
    as.numeric(observed_chi),
    as.integer(nsim),
    as.integer(screen_j_source_seed(seed)),
    isTRUE(sequential),
    as.integer(seq_limit[[1L]]),
    PACKAGE = "gRm"
  )
}

screen_j_exact_gamma_slices_native <- function(slices, observed_gamma, nsim, seed = NULL) {
  if (!screen_j_exact_native_available()) {
    return(NULL)
  }
  .Call(
    "gRm_screen_j_exact_gamma_slices",
    slices,
    as.numeric(observed_gamma),
    as.integer(nsim),
    as.integer(screen_j_source_seed(seed)),
    PACKAGE = "gRm"
  )
}

screen_j_item_pair_conditional_exact_native <- function(x,
                                                        y,
                                                        x_dim,
                                                        y_dim,
                                                        condition_values,
                                                        condition_dim,
                                                        valid,
                                                        nsim,
                                                        seed,
                                                        sequential,
                                                        seq_limit) {
  if (!screen_j_exact_native_available()) {
    return(NULL)
  }
  available <- tryCatch(
    is.loaded("gRm_screen_j_item_pair_conditional_exact", PACKAGE = "gRm"),
    error = function(e) FALSE
  )
  if (!available) {
    return(NULL)
  }
  raw <- .Call(
    "gRm_screen_j_item_pair_conditional_exact",
    as.integer(x),
    as.integer(y),
    as.integer(x_dim),
    as.integer(y_dim),
    as.integer(condition_values),
    as.integer(condition_dim),
    as.logical(valid),
    as.integer(nsim),
    as.integer(screen_j_source_seed(seed)),
    isTRUE(sequential),
    as.integer(seq_limit[[1L]]),
    PACKAGE = "gRm"
  )
  if (is.null(raw$p_chi_exact) && !is.null(raw$p_chi)) {
    raw$p_chi_exact <- raw$p_chi
  }
  if (is.null(raw$p_gamma_exact) && !is.null(raw$p_gamma)) {
    raw$p_gamma_exact <- raw$p_gamma
  }
  if (is.null(raw$exact_nsim) && !is.null(raw$nsim)) {
    raw$exact_nsim <- raw$nsim
  }
  raw
}

screen_j_prepare_exact_slices <- function(slices) {
  non_empty <- vapply(slices, function(slice) sum(slice) > 0, logical(1L))
  lapply(slices[non_empty], screen_j_prepare_exact_slice)
}

screen_j_prepare_exact_slice <- function(slice) {
  prepared <- exo_select_prepare_gentab1(slice)
  prepared$expected <- screen_j_expected_from_margins(
    prepared$row_total,
    prepared$col_total,
    prepared$grand_total
  )
  prepared
}

screen_j_expected_from_margins <- function(row_total, col_total, grand_total) {
  expected <- matrix(0, nrow = length(row_total), ncol = length(col_total))
  if (grand_total <= 0) {
    return(expected)
  }
  for (col in seq_along(col_total)) {
    col_share <- col_total[[col]] / grand_total
    for (row in seq_along(row_total)) {
      # Source trace: SKbigtab.Transfer_BT_to_XYZ_TABLE fills RTAB2 once as
      # row margin * (column margin / total). SKrandom.GENTAB1 then passes
      # that stored expected table to SkStat.RCCHI for every generated table.
      expected[row, col] <- row_total[[row]] * col_share
    }
  }
  expected
}

screen_j_rc_chi_square_prepared_expected <- function(tab, prepared) {
  expected <- prepared$expected
  if (is.null(expected)) {
    expected <- screen_j_expected_from_margins(
      prepared$row_total,
      prepared$col_total,
      prepared$grand_total
    )
  }
  positive <- expected > 0
  if (!any(positive)) {
    return(0)
  }
  residual <- tab[positive] - expected[positive]
  sum(residual * (residual / expected[positive]))
}

screen_j_exo_values <- function(project,
                                item_matrix,
                                item_score,
                                complete_items,
                                items,
                                backgrounds,
                                max_score,
                                partial) {
  n_items <- nrow(items)
  n_exo <- nrow(backgrounds)
  stat <- matrix(0, nrow = n_items, ncol = n_exo)
  p <- matrix(2, nrow = n_items, ncol = n_exo)
  ppq <- matrix(0, nrow = n_items, ncol = n_exo)
  pmq <- matrix(0, nrow = n_items, ncol = n_exo)
  score_stat <- rep(0, n_exo)
  score_p <- rep(2, n_exo)
  kind <- rep("Gamma", n_exo)
  exo_matrix <- matrix(0L, nrow = nrow(project$raw_data), ncol = n_exo)
  if (n_exo == 0L) {
    return(list(stat = stat, p = p, ppq = ppq, pmq = pmq, score_stat = score_stat, score_p = score_p, kind = kind, exo_values = exo_matrix))
  }
  dimnames(stat) <- dimnames(p) <- dimnames(ppq) <- dimnames(pmq) <- list(items$label_code, backgrounds$label_code)
  names(score_stat) <- names(score_p) <- names(kind) <- backgrounds$label_code

  for (exo_index in seq_len(n_exo)) {
    exo_values <- as.integer(project$raw_data[, backgrounds$position[[exo_index]]])
    exo_matrix[, exo_index] <- exo_values
    exo_dim <- backgrounds$raw_max[[exo_index]]
    exo_valid <- exo_values >= 1L & exo_values <= exo_dim
    use_gamma <- backgrounds$vtype[[exo_index]] > 2L || exo_dim == 2L
    kind[[exo_index]] <- if (use_gamma) "Gamma" else "CHI*2"

    score_tab <- screen_j_pair_table(item_score + 1L, exo_values, max_score + 1L, exo_dim, complete_items & exo_valid)
    if (use_gamma) {
      score_stats <- screen_rc_gamma(score_tab)
      score_stat[[exo_index]] <- score_stats$gamma
      score_p[[exo_index]] <- score_stats$p_value
    } else {
      score_stats <- screen_rc_chi(score_tab)
      score_stat[[exo_index]] <- score_stats$chi_square
      score_p[[exo_index]] <- score_stats$p_value
    }

    for (item_index in seq_len(n_items)) {
      item_values <- item_matrix[, item_index]
      rest_score <- item_score - item_values
      if (!partial) {
        valid_item <- item_values >= 1L & item_values <= items$raw_max[[item_index]]
        tab <- screen_j_pair_table(
          item_values, exo_values, items$raw_max[[item_index]], exo_dim,
          valid_item & exo_valid
        )
        if (use_gamma) {
          stats <- screen_rc_gamma(tab)
          stat[item_index, exo_index] <- stats$gamma
          p[item_index, exo_index] <- stats$p_value
        } else {
          stats <- screen_rc_chi(tab)
          stat[item_index, exo_index] <- stats$chi_square
          p[item_index, exo_index] <- stats$p_value
        }
      } else {
        strata <- screen_j_strata_table(
          item_values, exo_values, item_score + 1L,
          items$raw_max[[item_index]], exo_dim, max_score + 1L,
          complete_items & exo_valid & item_score > 0L & item_score < max_score
        )
        if (use_gamma) {
          stats <- screen_j_partial_gamma(strata)
          stat[item_index, exo_index] <- stats$gamma
          ppq[item_index, exo_index] <- stats$ppq
          pmq[item_index, exo_index] <- stats$pmq
          p[item_index, exo_index] <- if (isTRUE(attr(project, "screen_j_exact"))) {
            chi_stats <- screen_j_partial_chi(strata)
            screen_j_exact_partial_gamma(
              strata,
              stats$gamma,
              attr(project, "screen_j_nsim"),
              attr(project, "screen_j_seed"),
              chi_stats$stat,
              sequential = isTRUE(attr(project, "screen_j_repeated")),
              seq_limit = attr(project, "screen_j_seq_limit"),
              seq_p0 = attr(project, "screen_j_seq_p0"),
              seq_boundary = attr(project, "screen_j_seq_boundary")
            )
          } else {
            stats$p_value
          }
        } else {
          stats <- screen_j_partial_chi(strata)
          stat[item_index, exo_index] <- stats$stat
          p[item_index, exo_index] <- if (isTRUE(attr(project, "screen_j_exact"))) {
            screen_j_exact_partial_chi(
              strata,
              stats$stat,
              attr(project, "screen_j_nsim"),
              attr(project, "screen_j_seed"),
              sequential = isTRUE(attr(project, "screen_j_repeated")),
              seq_limit = attr(project, "screen_j_seq_limit"),
              seq_p0 = attr(project, "screen_j_seq_p0"),
              seq_boundary = attr(project, "screen_j_seq_boundary")
            )
          } else {
            stats$p_value
          }
        }
      }
    }
  }
  list(stat = stat, p = p, ppq = ppq, pmq = pmq, score_stat = score_stat, score_p = score_p, kind = kind, exo_values = exo_matrix)
}

screen_j_weighted_partial_gamma <- function(ppq, pmq) {
  n_items <- nrow(ppq)
  weighted <- matrix(0, nrow = n_items, ncol = n_items, dimnames = dimnames(ppq))
  diag(weighted) <- 1
  for (i in seq_len(n_items - 1L)) {
    for (j in seq.int(i + 1L, n_items)) {
      denom <- ppq[i, j] + ppq[j, i]
      value <- if (denom > 0) (pmq[i, j] + pmq[j, i]) / denom else 0
      weighted[i, j] <- value
      weighted[j, i] <- value
    }
  }
  weighted
}

screen_j_average_abs_partial_gamma <- function(gamma) {
  n_items <- nrow(gamma)
  if (n_items < 2L) {
    return(0)
  }
  total <- 0
  for (i in seq_len(n_items - 1L)) {
    for (j in seq.int(i + 1L, n_items)) {
      total <- total + abs(gamma[i, j] + gamma[j, i])
    }
  }
  total / (n_items * (n_items - 1L))
}

screen_j_problem_counts <- function(gamma, p, fdr05, i, j, temp_p = p) {
  positive <- 0L
  negative <- 0L
  if (temp_p[i, j] <= fdr05) {
    if (gamma[i, j] > 0) positive <- positive + 1L else negative <- negative + 1L
  }
  if (temp_p[j, i] <= fdr05) {
    if (gamma[j, i] > 0) positive <- positive + 1L else negative <- negative + 1L
  }
  c(positive = positive, negative = negative)
}

#' Select and finalize SCREEN J local-dependence evidence
#'
#' Reproduce the four-stage greedy selection in
#' `SKbias7.stepwise_elimination`, then apply the separate command-level
#' positive-local-dependence filter in `DIGRAM1f.pas`. The provisional result
#' is retained for diagnostics, whereas `matrix` and `rows` contain only the
#' terms that DIGRAM permits in the screen model.
#'
#' @param gamma Directed partial-gamma matrix.
#' @param p Directed partial-gamma p-value matrix.
#' @param weighted_gamma Symmetric weighted-partial-gamma matrix.
#' @param fdr05 Global SCREEN J Benjamini-Hochberg boundary at FDR 0.05.
#' @return A list containing final model `matrix` and `rows`, provisional
#'   `stepwise_matrix` and `stepwise_rows`, and excluded `negative_matrix` and
#'   `negative_rows`.
#' @keywords internal
#' @noRd
screen_j_stepwise_local_dependence <- function(gamma, p, weighted_gamma, fdr05) {
  n_items <- nrow(gamma)
  temp_p <- p
  selected <- matrix(FALSE, nrow = n_items, ncol = n_items, dimnames = dimnames(gamma))
  rows <- list()

  add_selected <- function(i, j, value, stage) {
    selected[i, j] <<- TRUE
    selected[j, i] <<- TRUE
    temp_p[i, ] <<- 2
    temp_p[j, ] <<- 2
    rows[[length(rows) + 1L]] <<- data.frame(
      row = i,
      col = j,
      pair = paste0(rownames(gamma)[[i]], colnames(gamma)[[j]]),
      value = value,
      stage = stage,
      stringsAsFactors = FALSE
    )
  }

  choose_pair <- function(target, sign) {
    best_i <- 0L
    best_j <- 0L
    best_value <- 0
    for (i in seq_len(n_items - 1L)) {
      for (j in seq.int(i + 1L, n_items)) {
        counts <- screen_j_problem_counts(gamma, p, fdr05, i, j, p)
        adjusted <- screen_j_problem_counts(gamma, p, fdr05, i, j, temp_p)
        if (counts[[target]] == sign$count && adjusted[[target]] > 0L) {
          value <- weighted_gamma[i, j]
          if ((sign$direction == "positive" && value > best_value) ||
              (sign$direction == "negative" && value < best_value)) {
            best_i <- i
            best_j <- j
            best_value <- value
          }
        }
      }
    }
    if (best_i > 0L) {
      add_selected(best_i, best_j, best_value, paste0(sign$direction, sign$count))
      TRUE
    } else {
      FALSE
    }
  }

  while (choose_pair("positive", list(count = 2L, direction = "positive"))) {}
  while (choose_pair("positive", list(count = 1L, direction = "positive"))) {}
  while (choose_pair("negative", list(count = 2L, direction = "negative"))) {}
  while (choose_pair("negative", list(count = 1L, direction = "negative"))) {}

  if (length(rows) == 0L) {
    rows_df <- data.frame(
      row = integer(),
      col = integer(),
      pair = character(),
      value = numeric(),
      stage = character(),
      stringsAsFactors = FALSE
    )
  } else {
    rows_df <- do.call(rbind, rows)
  }

  screen_j_finalize_local_dependence(gamma, selected, rows_df)
}

#' Apply DIGRAM's final SCREEN J local-dependence model rule
#'
#' The greedy SCREEN J procedure reports both positive and negative evidence.
#' Since DIGRAM version 3.37, a provisional pair becomes a screen-model term
#' only when the sum of its two directed partial gammas is strictly positive.
#'
#' @param gamma Directed partial-gamma matrix.
#' @param stepwise_matrix Symmetric logical matrix of provisional selections.
#' @param stepwise_rows Ordered provisional-selection table.
#' @return A local-dependence selection list with separate provisional, final,
#'   and excluded-negative representations.
#' @keywords internal
#' @noRd
screen_j_finalize_local_dependence <- function(gamma, stepwise_matrix, stepwise_rows) {
  rows <- stepwise_rows
  if (!is.data.frame(rows)) {
    rows <- data.frame()
  }

  # Source trace: DIGRAM1f command 3, parameters I/J, removes a provisional
  # LocalDependence[i,j] unless IJXgamma[i,j] + IJXgamma[j,i] > 0. This is an
  # unweighted sum of the two directed partial gammas; it is not the WPG sign
  # and it is not implied by the provisional stage label.
  if (nrow(rows)) {
    index_forward <- cbind(rows$row, rows$col)
    index_reverse <- cbind(rows$col, rows$row)
    rows$directed_gamma_sum <- gamma[index_forward] + gamma[index_reverse]
    rows$included <- !is.na(rows$directed_gamma_sum) & rows$directed_gamma_sum > 0
  } else {
    rows$directed_gamma_sum <- numeric()
    rows$included <- logical()
  }

  final_matrix <- matrix(
    FALSE,
    nrow = nrow(stepwise_matrix),
    ncol = ncol(stepwise_matrix),
    dimnames = dimnames(stepwise_matrix)
  )
  if (nrow(rows) && any(rows$included)) {
    included_rows <- rows[rows$included, , drop = FALSE]
    final_matrix[cbind(included_rows$row, included_rows$col)] <- TRUE
    final_matrix[cbind(included_rows$col, included_rows$row)] <- TRUE
  }

  negative_matrix <- stepwise_matrix & !final_matrix
  included_rows <- rows[rows$included, , drop = FALSE]
  negative_rows <- rows[!rows$included, , drop = FALSE]
  rownames(rows) <- NULL
  rownames(included_rows) <- NULL
  rownames(negative_rows) <- NULL

  list(
    matrix = final_matrix,
    rows = included_rows,
    stepwise_matrix = stepwise_matrix,
    stepwise_rows = rows,
    negative_matrix = negative_matrix,
    negative_rows = negative_rows
  )
}

screen_j_post_screen_dif <- function(item_matrix,
                                     exo_values,
                                     item_score,
                                     complete_items,
                                     items,
                                     backgrounds,
                                     max_score,
                                     partial_exo,
                                     initial_item_bias,
                                     fdr_01,
                                     exact = FALSE,
                                     nsim = 1000L,
                                     seed = NULL,
                                     repeated = FALSE,
                                     seq_limit = nsim,
                                     seq_p0 = 0.05,
                                     seq_boundary = 1.058) {
  n_items <- nrow(items)
  n_exo <- nrow(backgrounds)
  item_bias <- initial_item_bias
  item_bias_status <- matrix(
    0L,
    nrow = n_items,
    ncol = n_exo,
    dimnames = dimnames(initial_item_bias)
  )
  item_bias_status[initial_item_bias] <- 1L
  spurious_rows <- list()
  multiple_rows <- list()

  if (n_exo == 0L) {
    return(list(
      item_bias = item_bias,
      item_bias_status = item_bias_status,
      spurious_dif = list(rows = spurious_rows),
      multiple_dif = list(rows = multiple_rows)
    ))
  }

  for (exo_index in seq_len(n_exo)) {
    biased_items <- which(item_bias[, exo_index])
    entry <- NULL

    if (length(biased_items) == n_items) {
      retained <- which(partial_exo$p[, exo_index] <= fdr_01)
      item_bias[, exo_index] <- FALSE
      item_bias[retained, exo_index] <- TRUE
      biased_items <- retained
      entry <- list(
        background = exo_index,
        all_items_warning = TRUE,
        fdr01_reduction = TRUE,
        still_all_warning = length(biased_items) == n_items,
        iterations = list(),
        remaining = biased_items
      )
    }

    if (length(biased_items) > 1L && length(biased_items) < n_items) {
      analysis <- screen_j_spurious_item_bias_analysis(
        item_matrix = item_matrix,
        exo_values = exo_values[, exo_index],
        item_score = item_score,
        complete_items = complete_items,
        items = items,
        background = backgrounds[exo_index, , drop = FALSE],
        exo_index = exo_index,
        biased_items = biased_items,
        max_score = max_score,
        exact = exact,
        nsim = nsim,
        seed = seed,
        repeated = repeated,
        seq_limit = seq_limit,
        seq_p0 = seq_p0,
        seq_boundary = seq_boundary
      )
      if (!is.null(entry)) {
        analysis$all_items_warning <- entry$all_items_warning
        analysis$fdr01_reduction <- entry$fdr01_reduction
        analysis$still_all_warning <- entry$still_all_warning
      }
      spurious_rows[[length(spurious_rows) + 1L]] <- analysis
      item_bias[, exo_index] <- FALSE
      item_bias[analysis$remaining, exo_index] <- TRUE
      for (item_index in analysis$remaining) {
        item_bias_status[item_index, exo_index] <- 2L
      }
    } else {
      if (!is.null(entry)) {
        spurious_rows[[length(spurious_rows) + 1L]] <- entry
      }
      if (length(biased_items) == 1L) {
        item_bias_status[biased_items, exo_index] <- 2L
      }
    }
  }

  for (item_index in seq_len(n_items)) {
    biased_exo <- which(item_bias[item_index, ])
    if (length(biased_exo) > 1L) {
      analysis <- screen_j_multiple_item_bias_analysis(
        item_matrix = item_matrix,
        exo_values = exo_values,
        item_score = item_score,
        complete_items = complete_items,
        items = items,
        backgrounds = backgrounds,
        item_index = item_index,
        biased_exo = biased_exo,
        max_score = max_score,
        exact = exact,
        nsim = nsim,
        seed = seed,
        repeated = repeated,
        seq_limit = seq_limit,
        seq_p0 = seq_p0,
        seq_boundary = seq_boundary
      )
      multiple_rows[[length(multiple_rows) + 1L]] <- analysis
      for (exo_index in biased_exo) {
        if (exo_index %in% analysis$remaining) {
          item_bias_status[item_index, exo_index] <- 2L
        } else if (item_bias[item_index, exo_index]) {
          item_bias_status[item_index, exo_index] <- 1L
        }
      }
      item_bias[item_index, ] <- FALSE
      item_bias[item_index, analysis$remaining] <- TRUE
    } else if (length(biased_exo) == 1L) {
      item_bias_status[item_index, biased_exo] <- 2L
    }
  }

  list(
    item_bias = item_bias,
    item_bias_status = item_bias_status,
    spurious_dif = list(rows = spurious_rows),
    multiple_dif = list(rows = multiple_rows)
  )
}

screen_j_spurious_item_bias_analysis <- function(item_matrix,
                                                exo_values,
                                                item_score,
                                                complete_items,
                                                items,
                                                background,
                                                exo_index,
                                                biased_items,
                                                max_score,
                                                exact = FALSE,
                                                nsim = 1000L,
                                                seed = NULL,
                                                repeated = FALSE,
                                                seq_limit = nsim,
                                                seq_p0 = 0.05,
                                                seq_boundary = 1.058) {
  current <- biased_items
  iterations <- list()
  use_gamma <- background$vtype[[1L]] > 2L || background$raw_max[[1L]] == 2L

  repeat {
    rows <- vector("list", length(current))
    p_min <- rep(0, length(current))
    for (candidate_index in seq_along(current)) {
      item_index <- current[[candidate_index]]
      other_items <- setdiff(current, item_index)
      condition_values <- cbind(
        if (length(other_items) > 0L) item_matrix[, other_items, drop = FALSE] else NULL,
        item_score
      )
      condition_dims <- c(items$raw_max[other_items], max_score - 1L)
      test <- screen_j_conditional_bias_test(
        x = exo_values,
        y = item_matrix[, item_index],
        x_dim = background$raw_max[[1L]],
        y_dim = items$raw_max[[item_index]],
        condition_values = condition_values,
        condition_dims = condition_dims,
        valid = complete_items &
          exo_values >= 1L & exo_values <= background$raw_max[[1L]] &
          item_score > 0L & item_score < max_score,
        exact = exact,
        nsim = nsim,
        seed = seed,
        repeated = repeated,
        native = screen_j_conditional_native_allowed(repeated, seq_p0, seq_boundary),
        seq_limit = seq_limit,
        seq_p0 = seq_p0,
        seq_boundary = seq_boundary
      )
      p_min[[candidate_index]] <- screen_j_source_stepwise_p_min(test, use_gamma, exact)
      rows[[candidate_index]] <- data.frame(
        index = candidate_index,
        hypothesis = screen_j_hypothesis_label(
          background$label_code[[1L]],
          items$label_code[[item_index]],
          c(items$label_code[other_items], "#")
        ),
        candidate = item_index,
        candidate_label = items$label_code[[item_index]],
        candidate_name = items$name[[item_index]],
        chi_square = test$chi_square,
        df = test$df,
        p_chi = test$p_chi,
        p_chi_asymp = test$p_chi_asymp,
        p_chi_exact = test$p_chi_exact,
        gamma = test$gamma,
        p_gamma = test$p_gamma,
        p_gamma_asymp = test$p_gamma_asymp,
        p_gamma_exact = test$p_gamma_exact,
        p_min = p_min[[candidate_index]],
        use_gamma = use_gamma,
        n = test$exact_nsim,
        stringsAsFactors = FALSE
      )
    }

    rows_df <- do.call(rbind, rows)
    rows_df$bh_05 <- source_bh_critical(screen_j_exa_p_values(rows_df), 0.05)
    rows_df$bh_01 <- source_bh_critical(screen_j_exa_p_values(rows_df), 0.01)
    excluded <- integer()
    excluded_label <- NA_character_
    excluded_name <- NA_character_
    if (length(p_min) > 0L && max(p_min) > 0.05) {
      excluded_index <- which.max(p_min)
      excluded <- current[[excluded_index]]
      excluded_label <- items$label_code[[excluded]]
      excluded_name <- items$name[[excluded]]
    }
    iterations[[length(iterations) + 1L]] <- list(
      rows = rows_df,
      excluded = excluded,
      excluded_label = excluded_label,
      excluded_name = excluded_name
    )
    if (length(excluded) == 0L) {
      break
    }
    current <- current[current != excluded]
    if (length(current) == 0L) {
      break
    }
  }

  list(
    background = exo_index,
    all_items_warning = FALSE,
    fdr01_reduction = FALSE,
    still_all_warning = FALSE,
    iterations = iterations,
    remaining = current
  )
}

screen_j_multiple_item_bias_analysis <- function(item_matrix,
                                                exo_values,
                                                item_score,
                                                complete_items,
                                                items,
                                                backgrounds,
                                                item_index,
                                                biased_exo,
                                                max_score,
                                                exact = FALSE,
                                                nsim = 1000L,
                                                seed = NULL,
                                                repeated = FALSE,
                                                seq_limit = nsim,
                                                seq_p0 = 0.05,
                                                seq_boundary = 1.058) {
  current <- biased_exo
  iterations <- list()

  repeat {
    rows <- vector("list", length(current))
    p_min <- rep(0, length(current))
    for (candidate_index in seq_along(current)) {
      exo_index <- current[[candidate_index]]
      other_exo <- setdiff(current, exo_index)
      condition_values <- cbind(
        if (length(other_exo) > 0L) exo_values[, other_exo, drop = FALSE] else NULL,
        item_score
      )
      condition_dims <- c(backgrounds$raw_max[other_exo], max_score - 1L)
      valid <- complete_items & item_score > 0L & item_score < max_score
      for (valid_exo in current) {
        valid <- valid &
          exo_values[, valid_exo] >= 1L &
          exo_values[, valid_exo] <= backgrounds$raw_max[[valid_exo]]
      }
      test <- screen_j_conditional_bias_test(
        x = item_matrix[, item_index],
        y = exo_values[, exo_index],
        x_dim = items$raw_max[[item_index]],
        y_dim = backgrounds$raw_max[[exo_index]],
        condition_values = condition_values,
        condition_dims = condition_dims,
        valid = valid,
        exact = exact,
        nsim = nsim,
        seed = seed,
        repeated = repeated,
        native = screen_j_conditional_native_allowed(repeated, seq_p0, seq_boundary),
        seq_limit = seq_limit,
        seq_p0 = seq_p0,
        seq_boundary = seq_boundary
      )
      use_gamma <- backgrounds$vtype[[exo_index]] > 2L || backgrounds$raw_max[[exo_index]] == 2L
      p_min[[candidate_index]] <- screen_j_source_stepwise_p_min(test, use_gamma, exact)
      rows[[candidate_index]] <- data.frame(
        index = candidate_index,
        hypothesis = screen_j_hypothesis_label(
          items$label_code[[item_index]],
          backgrounds$label_code[[exo_index]],
          c(backgrounds$label_code[other_exo], "#")
        ),
        candidate = exo_index,
        candidate_label = backgrounds$label_code[[exo_index]],
        candidate_name = backgrounds$name[[exo_index]],
        chi_square = test$chi_square,
        df = test$df,
        p_chi = test$p_chi,
        p_chi_asymp = test$p_chi_asymp,
        p_chi_exact = test$p_chi_exact,
        gamma = test$gamma,
        p_gamma = test$p_gamma,
        p_gamma_asymp = test$p_gamma_asymp,
        p_gamma_exact = test$p_gamma_exact,
        p_min = p_min[[candidate_index]],
        use_gamma = use_gamma,
        n = test$exact_nsim,
        stringsAsFactors = FALSE
      )
    }

    rows_df <- do.call(rbind, rows)
    rows_df$bh_05 <- source_bh_critical(screen_j_exa_p_values(rows_df), 0.05)
    rows_df$bh_01 <- source_bh_critical(screen_j_exa_p_values(rows_df), 0.01)
    excluded <- integer()
    excluded_label <- NA_character_
    excluded_name <- NA_character_
    if (length(p_min) > 0L && max(p_min) > 0.05) {
      excluded_index <- which.max(p_min)
      excluded <- current[[excluded_index]]
      excluded_label <- backgrounds$label_code[[excluded]]
      excluded_name <- backgrounds$name[[excluded]]
    }
    iterations[[length(iterations) + 1L]] <- list(
      rows = rows_df,
      excluded = excluded,
      excluded_label = excluded_label,
      excluded_name = excluded_name
    )
    if (length(excluded) == 0L) {
      break
    }
    current <- current[current != excluded]
    if (length(current) == 0L) {
      break
    }
  }

  list(item = item_index, iterations = iterations, remaining = current)
}

# Source trace: source/PAS_skunits/SKbias3.pas::XYZ_bias_ANALYSE receives the
# conditional item/background table prepared by
# source/PAS_skunits/SKxyz1.PAS::MAKE_XYZ_TABLE or
# source/PAS_skunits/SKbigtab.pas::Transfer_BT_to_XYZ_TABLE. It then computes
# chi-square, gamma, asymptotic p-values, and optional exact/repeated
# random-table summaries. This R helper is the shared implementation for those
# screen J conditional tests.
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
      # Source trace: source/PAS_skunits/SKxyz1.PAS::MAKE_XYZ_TABLE only
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
  # Source trace: source/PAS_skunits/SKrandom.pas::GENTAB1 is the exact
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

screen_j_hypothesis_label <- function(first, second, given) {
  paste0(first, "&", second, "|", paste(given, collapse = ""))
}

screen_j_exa_p_values <- function(rows) {
  # Source trace: source/PAS_skunits/SKexa1.pas::EXA_SUMMARY1_2 skips
  # hypotheses with RESULTS[hypnr, 2] = 0 ("No tests") before adding p-values
  # to the Benjamini-Hochberg input vector.
  tested <- rows$df > 0L
  p_values <- rows$p_chi[tested]
  if (any(rows$use_gamma[tested])) {
    p_values <- c(p_values, rows$p_gamma[tested & rows$use_gamma])
  }
  p_values
}

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

# Source trace: source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening builds
# score-effect screening rows from the same conditional bias-test machinery
# used by screen J. The R helper returns those rows as numeric values rather
# than writing DIGRAM text.
screen_j_score_effect_values <- function(project,
                                         item_matrix,
                                         item_score,
                                         complete_items,
                                         backgrounds,
                                         max_score,
                                         exact = FALSE,
                                         nsim = 1000L,
                                         seed = NULL,
                                         repeated = FALSE,
                                         seq_limit = nsim,
                                         seq_p0 = 0.05,
                                         seq_boundary = 1.058) {
  n_exo <- nrow(backgrounds)
  maxdim <- 57L
  score_top <- maxdim - 1L
  score_category <- pmin(item_score, score_top) + 1L
  score_valid <- complete_items
  rows <- vector("list", n_exo)
  if (n_exo == 0L) {
    return(list(
      rows = data.frame(),
      fdr_05 = 0,
      fdr_01 = 0,
      selected = character(),
      score_top = score_top,
      exact = exact,
      screening = list(rows = data.frame(), iterations = list())
    ))
  }
  p_values <- numeric(0)
  nsim <- as.integer(nsim)
  marginal_p_min <- rep(1, n_exo)
  exo_matrix <- matrix(0L, nrow = nrow(project$raw_data), ncol = n_exo)
  for (exo_index in seq_len(n_exo)) {
    exo_values <- as.integer(project$raw_data[, backgrounds$position[[exo_index]]])
    exo_matrix[, exo_index] <- exo_values
    exo_dim <- backgrounds$raw_max[[exo_index]]
    exo_valid <- exo_values >= 1L & exo_values <= exo_dim
    test <- screen_j_conditional_bias_test(
      x = score_category,
      y = exo_values,
      x_dim = score_top + 1L,
      y_dim = exo_dim,
      condition_values = matrix(integer(), nrow = nrow(item_matrix), ncol = 0L),
      condition_dims = integer(),
      valid = score_valid & exo_valid,
      exact = exact,
      nsim = nsim,
      seed = seed,
      repeated = repeated,
      native = screen_j_conditional_native_allowed(repeated, seq_p0, seq_boundary),
      seq_limit = seq_limit,
      seq_p0 = seq_p0,
      seq_boundary = seq_boundary
    )
    p_chi <- test$p_chi
    p_gamma <- test$p_gamma
    p_chi_exact <- test$p_chi_exact
    p_gamma_exact <- test$p_gamma_exact
    exact_nsim <- test$exact_nsim
    marginal_p_min[[exo_index]] <- min(p_chi, p_gamma)
    p_values <- c(p_values, p_chi, p_gamma)
    rows[[exo_index]] <- data.frame(
      label = backgrounds$label_code[[exo_index]],
      name = backgrounds$name[[exo_index]],
      chi_square = test$chi_square,
      df = test$df,
      p_chi = p_chi,
      p_chi_asymp = test$p_chi_asymp,
      p_chi_exact = p_chi_exact,
      p_chi_exact_low = if (isTRUE(exact)) screen_j_conflimit99_field(exact_nsim, p_chi_exact, "low") else NA_real_,
      p_chi_exact_high = if (isTRUE(exact)) screen_j_conflimit99_field(exact_nsim, p_chi_exact, "high") else NA_real_,
      gamma = test$gamma,
      p_gamma = p_gamma,
      p_gamma_asymp = test$p_gamma_asymp,
      p_gamma_exact = p_gamma_exact,
      p_gamma_exact_low = if (isTRUE(exact)) screen_j_conflimit99_field(exact_nsim, p_gamma_exact, "low") else NA_real_,
      p_gamma_exact_high = if (isTRUE(exact)) screen_j_conflimit99_field(exact_nsim, p_gamma_exact, "high") else NA_real_,
      exact_nsim = exact_nsim,
      selected = FALSE,
      stringsAsFactors = FALSE
    )
  }
  rows_df <- do.call(rbind, rows)
  screening <- screen_j_stepwise_score_screening(
    score_category = score_category,
    score_valid = score_valid,
    score_in_range = item_score <= score_top,
    exo_values = exo_matrix,
    backgrounds = backgrounds,
    marginal_p_min = marginal_p_min,
    score_top = score_top,
    exact = exact,
    nsim = nsim,
    seed = seed,
    repeated = repeated,
    seq_limit = seq_limit,
    seq_p0 = seq_p0,
    seq_boundary = seq_boundary
  )
  rows_df$selected <- seq_len(n_exo) %in% screening$remaining
  report_rows <- screen_j_score_report_rows(rows_df, screening$rows, exact)
  p_values <- c(p_values, screen_j_score_screening_p_values(screening$rows))
  list(
    rows = report_rows,
    fdr_05 = source_bh_critical(p_values, 0.05),
    fdr_01 = source_bh_critical(p_values, 0.01),
    marginal_selected = rows_df$label[marginal_p_min <= 0.05],
    selected = rows_df$label[rows_df$selected],
    score_top = score_top,
    exact = exact,
    screening = screening
  )
}

# Source trace: source/PAS_skunits/SKbias13.pas::StepwiseScoreScreening runs
# the marginal and conditional score-screening loop. The R function keeps the
# same candidate order and source random-table convention when exact or
# repeated inference is requested.
screen_j_stepwise_score_screening <- function(score_category,
                                              score_valid,
                                              score_in_range,
                                              exo_values,
                                              backgrounds,
                                              marginal_p_min,
                                              score_top,
                                              exact = FALSE,
                                              nsim = 1000L,
                                              seed = NULL,
                                              repeated = FALSE,
                                              seq_limit = nsim,
                                              seq_p0 = 0.05,
                                              seq_boundary = 1.058) {
  n_exo <- nrow(backgrounds)
  current <- which(marginal_p_min <= 0.05)
  iterations <- list()
  if (length(current) <= 1L) {
    return(list(rows = data.frame(), iterations = iterations, remaining = current))
  }

  while (length(current) > 1L) {
    rows <- vector("list", length(current))
    p_min <- rep(0, length(current))
    for (candidate_index in seq_along(current)) {
      exo_index <- current[[candidate_index]]
      other_exo <- setdiff(current, exo_index)
      condition_values <- if (length(other_exo) > 0L) {
        exo_values[, other_exo, drop = FALSE]
      } else {
        matrix(integer(), nrow = nrow(exo_values), ncol = 0L)
      }
      valid <- score_valid & score_in_range
      for (valid_exo in current) {
        valid <- valid &
          exo_values[, valid_exo] >= 1L &
          exo_values[, valid_exo] <= backgrounds$raw_max[[valid_exo]]
      }
      test <- screen_j_conditional_bias_test(
        x = score_category,
        y = exo_values[, exo_index],
        x_dim = score_top + 1L,
        y_dim = backgrounds$raw_max[[exo_index]],
        condition_values = condition_values,
        condition_dims = backgrounds$raw_max[other_exo],
        valid = valid,
        exact = exact,
        nsim = nsim,
        seed = seed,
        repeated = repeated,
        native = screen_j_conditional_native_allowed(repeated, seq_p0, seq_boundary),
        seq_limit = seq_limit,
        seq_p0 = seq_p0,
        seq_boundary = seq_boundary
      )
      p_min[[candidate_index]] <- screen_j_source_stepwise_p_min(test, use_gamma = TRUE, exact = exact)
      rows[[candidate_index]] <- data.frame(
        index = candidate_index,
        hypothesis = screen_j_hypothesis_label(
          "#",
          backgrounds$label_code[[exo_index]],
          backgrounds$label_code[other_exo]
        ),
        candidate = exo_index,
        candidate_label = backgrounds$label_code[[exo_index]],
        candidate_name = backgrounds$name[[exo_index]],
        chi_square = test$chi_square,
        df = test$df,
        p_chi = test$p_chi,
        p_chi_asymp = test$p_chi_asymp,
        p_chi_exact = test$p_chi_exact,
        gamma = test$gamma,
        p_gamma = test$p_gamma,
        p_gamma_asymp = test$p_gamma_asymp,
        p_gamma_exact = test$p_gamma_exact,
        p_min = p_min[[candidate_index]],
        n = test$exact_nsim,
        stringsAsFactors = FALSE
      )
    }

    rows_df <- do.call(rbind, rows)
    excluded <- integer()
    excluded_label <- NA_character_
    excluded_name <- NA_character_
    if (length(p_min) > 0L && max(p_min) > 0.05) {
      excluded_index <- which.max(p_min)
      excluded <- current[[excluded_index]]
      excluded_label <- backgrounds$label_code[[excluded]]
      excluded_name <- backgrounds$name[[excluded]]
    }
    iterations[[length(iterations) + 1L]] <- list(
      rows = rows_df,
      excluded = excluded,
      excluded_label = excluded_label,
      excluded_name = excluded_name
    )
    if (length(excluded) == 0L) {
      break
    }
    current <- current[current != excluded]
    if (length(current) <= 1L) {
      break
    }
  }

  rows_df <- if (length(iterations) == 0L) data.frame() else do.call(rbind, lapply(iterations, `[[`, "rows"))
  list(rows = rows_df, iterations = iterations, remaining = current)
}

screen_j_score_screening_p_values <- function(rows) {
  if (is.null(rows) || nrow(rows) == 0L) {
    return(numeric(0))
  }
  c(rows$p_chi, rows$p_gamma)
}

screen_j_score_report_rows <- function(marginal_rows, screening_rows, exact) {
  marginal_rows$hypothesis <- paste0("#&", marginal_rows$label)
  if (is.null(screening_rows) || nrow(screening_rows) == 0L) {
    return(marginal_rows)
  }
  screening <- data.frame(
    label = screening_rows$candidate_label,
    name = screening_rows$candidate_name,
    chi_square = screening_rows$chi_square,
    df = screening_rows$df,
    p_chi = screening_rows$p_chi,
    p_chi_asymp = screening_rows$p_chi_asymp,
    p_chi_exact = screening_rows$p_chi_exact,
    p_chi_exact_low = if (isTRUE(exact)) screen_j_conflimit99_field(screening_rows$n, screening_rows$p_chi_exact, "low") else NA_real_,
    p_chi_exact_high = if (isTRUE(exact)) screen_j_conflimit99_field(screening_rows$n, screening_rows$p_chi_exact, "high") else NA_real_,
    gamma = screening_rows$gamma,
    p_gamma = screening_rows$p_gamma,
    p_gamma_asymp = screening_rows$p_gamma_asymp,
    p_gamma_exact = screening_rows$p_gamma_exact,
    p_gamma_exact_low = if (isTRUE(exact)) screen_j_conflimit99_field(screening_rows$n, screening_rows$p_gamma_exact, "low") else NA_real_,
    p_gamma_exact_high = if (isTRUE(exact)) screen_j_conflimit99_field(screening_rows$n, screening_rows$p_gamma_exact, "high") else NA_real_,
    exact_nsim = screening_rows$n,
    selected = FALSE,
    hypothesis = screening_rows$hypothesis,
    stringsAsFactors = FALSE
  )
  rbind(marginal_rows, screening[, names(marginal_rows)])
}

screen_j_conflimit99_field <- function(n, p, field) {
  n <- rep_len(as.integer(n), length(p))
  p <- as.numeric(p)
  vapply(seq_along(p), function(i) {
    source_conflimit99(n[[i]], p[[i]])[[field]]
  }, numeric(1L))
}
