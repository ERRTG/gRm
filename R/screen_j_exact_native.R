#' Internal screen j exact native available helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @return The internal `screen_j_exact_native_available()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal screen j conditional native allowed helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param repeated Whether to use repeated sequential simulation.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_conditional_native_allowed()` computation result.
#' @keywords internal
#' @noRd
screen_j_conditional_native_allowed <- function(repeated, seq_p0, seq_boundary) {
  screen_j_conditional_native_controls_allowed(repeated, seq_p0, seq_boundary) &&
    screen_j_conditional_native_source_faithful()
}

#' Internal screen j conditional native controls allowed helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param repeated Whether to use repeated sequential simulation.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_conditional_native_controls_allowed()` computation result.
#' @keywords internal
#' @noRd
screen_j_conditional_native_controls_allowed <- function(repeated, seq_p0, seq_boundary) {
  !isTRUE(repeated) ||
    (abs(as.numeric(seq_p0[[1L]]) - 0.05) < 1e-15 &&
      abs(as.numeric(seq_boundary[[1L]]) - 1.058) < 1e-15)
}

#' Internal screen j conditional native probe fixtures helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @return The internal `screen_j_conditional_native_probe_fixtures()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal screen j conditional native probe matches helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param fixture Internal `fixture` value used by this helper.
#' @param controls Internal `controls` value used by this helper.
#' @param repeated Whether to use repeated sequential simulation.
#' @return The internal `screen_j_conditional_native_probe_matches()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal screen j item pair native allowed helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param repeated Whether to use repeated sequential simulation.
#' @param seq_p0 Internal `seq_p0` value used by this helper.
#' @param seq_boundary Internal `seq_boundary` value used by this helper.
#' @return The internal `screen_j_item_pair_native_allowed()` computation result.
#' @keywords internal
#' @noRd
screen_j_item_pair_native_allowed <- function(repeated, seq_p0, seq_boundary) {
  screen_j_conditional_native_controls_allowed(repeated, seq_p0, seq_boundary) &&
    screen_j_item_pair_native_source_faithful()
}

#' Internal screen j item pair native probe fixtures helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @return The internal `screen_j_item_pair_native_probe_fixtures()` computation result.
#' @keywords internal
#' @noRd
screen_j_item_pair_native_probe_fixtures <- function() {
  tie_heavy <- list(
    x = c(rep(1L, 6L), 2L, 1L),
    y = c(rep(1L, 7L), 2L),
    x_dim = 2L,
    y_dim = 2L,
    condition_values = rep(1L, 8L),
    condition_dim = 1L,
    valid = rep(TRUE, 8L)
  )
  ordinary <- list(
    x = c(1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L),
    y = c(1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L),
    x_dim = 3L,
    y_dim = 2L,
    condition_values = c(1L, 2L, 2L, 3L, 3L, 4L, 3L, 3L, 4L),
    condition_dim = 4L,
    valid = rep(TRUE, 9L)
  )
  list(
    list(args = tie_heavy, fixed = c(40L, 9L, 40L), repeated = c(40L, 9L, 1L)),
    list(args = ordinary, fixed = c(41L, 9L, 41L), repeated = c(80L, 17L, 3L))
  )
}

#' Internal screen j item pair probe reference helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param args Internal `args` value used by this helper.
#' @param controls Internal `controls` value used by this helper.
#' @param sequential Internal `sequential` value used by this helper.
#' @return The internal `screen_j_item_pair_probe_reference()` computation result.
#' @keywords internal
#' @noRd
screen_j_item_pair_probe_reference <- function(args, controls, sequential) {
  strata <- screen_j_strata_table(
    args$x,
    args$y,
    args$condition_values,
    args$x_dim,
    args$y_dim,
    args$condition_dim,
    args$valid
  )
  chi <- screen_j_partial_chi(strata)
  gamma <- screen_j_partial_gamma(strata)
  exact <- screen_j_exact_chi_gamma_prepared_r(
    screen_j_prepare_exact_slices(screen_j_strata_slices(strata)),
    chi$stat,
    gamma$gamma,
    controls[[1L]],
    controls[[2L]],
    sequential = sequential,
    seq_limit = controls[[3L]]
  )
  list(
    chi_square = chi$stat,
    df = chi$df,
    gamma = gamma$gamma,
    ppq = gamma$ppq,
    pmq = gamma$pmq,
    s = gamma$s,
    p_chi_asymp = chi$p_value,
    p_gamma_asymp = gamma$p_value,
    p_chi_exact = exact$p_chi,
    p_gamma_exact = exact$p_gamma,
    exact_nsim = exact$nsim,
    chi_exceed = exact$chi_exceed,
    gamma_exceed = exact$gamma_exceed,
    draw_count = exact$draw_count,
    final_seed = exact$final_seed
  )
}

#' Internal screen j item pair probe matches helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param fixture Internal `fixture` value used by this helper.
#' @param controls Internal `controls` value used by this helper.
#' @param sequential Internal `sequential` value used by this helper.
#' @return The internal `screen_j_item_pair_probe_matches()` computation result.
#' @keywords internal
#' @noRd
screen_j_item_pair_probe_matches <- function(fixture, controls, sequential) {
  reference <- screen_j_item_pair_probe_reference(
    fixture$args,
    controls,
    sequential
  )
  native <- do.call(
    screen_j_item_pair_conditional_exact_native,
    c(
      fixture$args,
      list(
        nsim = controls[[1L]],
        seed = controls[[2L]],
        sequential = sequential,
        seq_limit = controls[[3L]]
      )
    )
  )
  fields <- c(
    "chi_square", "df", "gamma", "ppq", "pmq", "s",
    "p_chi_asymp", "p_gamma_asymp", "p_chi_exact", "p_gamma_exact",
    "exact_nsim", "chi_exceed", "gamma_exceed", "draw_count", "final_seed"
  )
  !is.null(native) && all(vapply(fields, function(field) {
    isTRUE(all.equal(
      native[[field]],
      reference[[field]],
      tolerance = 1e-12,
      check.attributes = FALSE
    ))
  }, logical(1L)))
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
      fixtures <- screen_j_item_pair_native_probe_fixtures()
      all(vapply(fixtures, function(fixture) {
        screen_j_item_pair_probe_matches(fixture, fixture$fixed, FALSE) &&
          screen_j_item_pair_probe_matches(fixture, fixture$repeated, TRUE)
      }, logical(1L)))
    }, error = function(e) FALSE)
    cached
  }
})

#' Internal screen j exact chi gamma slices native helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param slices Internal `slices` value used by this helper.
#' @param observed_chi Internal `observed_chi` value used by this helper.
#' @param observed_gamma Internal `observed_gamma` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param sequential Internal `sequential` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @return The internal `screen_j_exact_chi_gamma_slices_native()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal screen j exact chi gamma trace slices native helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param slices Internal `slices` value used by this helper.
#' @param observed_chi Internal `observed_chi` value used by this helper.
#' @param observed_gamma Internal `observed_gamma` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param sequential Internal `sequential` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @return The internal `screen_j_exact_chi_gamma_trace_slices_native()` computation result.
#' @keywords internal
#' @noRd
screen_j_exact_chi_gamma_trace_slices_native <- function(slices,
                                                         observed_chi,
                                                         observed_gamma,
                                                         nsim,
                                                         seed = NULL,
                                                         sequential = FALSE,
                                                         seq_limit = nsim) {
  if (!screen_j_exact_native_available()) {
    return(NULL)
  }
  available <- tryCatch(
    is.loaded("gRm_screen_j_exact_chi_gamma_trace_slices", PACKAGE = "gRm"),
    error = function(e) FALSE
  )
  if (!available) {
    return(NULL)
  }
  raw <- .Call(
    "gRm_screen_j_exact_chi_gamma_trace_slices",
    slices,
    as.numeric(observed_chi),
    as.numeric(observed_gamma),
    as.integer(nsim),
    as.integer(screen_j_source_seed(seed)),
    isTRUE(sequential),
    as.integer(seq_limit[[1L]]),
    PACKAGE = "gRm"
  )
  table_state <- as.data.frame(raw$trajectory$table_state, stringsAsFactors = FALSE)
  tables <- Map(function(table, random_draws, index) {
    list(
      sim = table_state$sim[[index]],
      slice = table_state$slice[[index]],
      table = table,
      random_draws = random_draws,
      chi_square = table_state$chi_square[[index]],
      ppq = table_state$ppq[[index]],
      pmq = table_state$pmq[[index]],
      draw_count = table_state$draw_count[[index]],
      final_seed = table_state$final_seed[[index]]
    )
  }, raw$trajectory$tables, raw$trajectory$random_draws, seq_along(raw$trajectory$tables))
  raw$trajectory <- list(
    tables = tables,
    simulations = as.data.frame(
      raw$trajectory$simulations,
      stringsAsFactors = FALSE
    )
  )
  raw
}

#' Internal screen j exact chi slices native helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param slices Internal `slices` value used by this helper.
#' @param observed_chi Internal `observed_chi` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param sequential Internal `sequential` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @return The internal `screen_j_exact_chi_slices_native()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal screen j exact gamma slices native helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param slices Internal `slices` value used by this helper.
#' @param observed_gamma Internal `observed_gamma` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @return The internal `screen_j_exact_gamma_slices_native()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal screen j item pair conditional exact native helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param x Object or value to process.
#' @param y Internal `y` value used by this helper.
#' @param x_dim Internal `x_dim` value used by this helper.
#' @param y_dim Internal `y_dim` value used by this helper.
#' @param condition_values Internal `condition_values` value used by this helper.
#' @param condition_dim Internal `condition_dim` value used by this helper.
#' @param valid Internal `valid` value used by this helper.
#' @param nsim Requested simulation count.
#' @param seed Random-stream seed.
#' @param sequential Internal `sequential` value used by this helper.
#' @param seq_limit Internal `seq_limit` value used by this helper.
#' @return The internal `screen_j_item_pair_conditional_exact_native()` computation result.
#' @keywords internal
#' @noRd
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
  raw$chi_square <- raw$chi
  raw$p_chi_asymp <- if (raw$df > 0L) source_pfchi(raw$df, raw$chi) else 1
  raw$p_gamma_asymp <- if (raw$ppq <= 0) {
    2
  } else if (raw$s <= 0) {
    2 * source_tail_norm(4, TRUE)
  } else {
    2 * source_tail_norm(abs(raw$gamma / sqrt(raw$s)), TRUE)
  }
  raw$draw_count <- raw$rng_draws
  raw
}

#' Internal screen j prepare exact slices helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param slices Internal `slices` value used by this helper.
#' @return The internal `screen_j_prepare_exact_slices()` computation result.
#' @keywords internal
#' @noRd
screen_j_prepare_exact_slices <- function(slices) {
  non_empty <- vapply(slices, function(slice) sum(slice) > 0, logical(1L))
  lapply(slices[non_empty], screen_j_prepare_exact_slice)
}

#' Internal screen j prepare exact slice helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param slice Internal `slice` value used by this helper.
#' @return The internal `screen_j_prepare_exact_slice()` computation result.
#' @keywords internal
#' @noRd
screen_j_prepare_exact_slice <- function(slice) {
  prepared <- exo_select_prepare_gentab1(slice)
  prepared$expected <- screen_j_expected_from_margins(
    prepared$row_total,
    prepared$col_total,
    prepared$grand_total
  )
  prepared
}

#' Internal screen j expected from margins helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param row_total Internal `row_total` value used by this helper.
#' @param col_total Internal `col_total` value used by this helper.
#' @param grand_total Internal `grand_total` value used by this helper.
#' @return The internal `screen_j_expected_from_margins()` computation result.
#' @keywords internal
#' @noRd
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

#' Internal screen j rc chi square prepared expected helper
#'
#' Supports the screen j values implementation while preserving its internal contract.
#' Source trace: `source/digram_source_20260817/skunits/SKbias7.pas::Item_Screening`.
#' @param tab Internal `tab` value used by this helper.
#' @param prepared Internal `prepared` value used by this helper.
#' @return The internal `screen_j_rc_chi_square_prepared_expected()` computation result.
#' @keywords internal
#' @noRd
screen_j_rc_chi_square_prepared_expected <- function(tab, prepared) {
  expected <- prepared$expected
  if (is.null(expected)) {
    expected <- screen_j_expected_from_margins(
      prepared$row_total,
      prepared$col_total,
      prepared$grand_total
    )
  }
  chi_square <- 0
  # Source trace: source/digram_source_20260817/skunits/SkStat.pas::RCCHI accumulates Pascal REAL
  # cells in first-index/second-index order. Do not use R's column-major logical
  # extraction here: exact comparisons are deliberately sensitive to this sum.
  for (row in seq_len(nrow(expected))) {
    for (col in seq_len(ncol(expected))) {
      if (expected[row, col] > 0) {
        residual <- tab[row, col] - expected[row, col]
        chi_square <- chi_square + residual * (residual / expected[row, col])
      }
    }
  }
  chi_square
}
