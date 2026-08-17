screen_j_tie_heavy_trace_fixture <- function() {
  strata <- array(as.integer(c(6L, 1L, 1L, 0L)), dim = c(2L, 2L, 1L))
  slices <- screen_j_strata_slices(strata)
  list(
    strata = strata,
    slices = slices,
    prepared = screen_j_prepare_exact_slices(slices),
    chi = screen_j_partial_chi(strata),
    gamma = screen_j_partial_gamma(strata)
  )
}

screen_j_tie_heavy_r_trace <- function(sequential) {
  fixture <- screen_j_tie_heavy_trace_fixture()
  screen_j_exact_chi_gamma_prepared_r(
    fixture$prepared,
    fixture$chi$stat,
    fixture$gamma$gamma,
    nsim = 40L,
    seed = 9L,
    sequential = sequential,
    seq_limit = 1L,
    seq_p0 = 0.05,
    seq_boundary = 1.058,
    trace = TRUE
  )
}

test_that("SCREEN J exact ties use the source REAL-to-SINGLE comparison boundary", {
  fixture <- screen_j_tie_heavy_trace_fixture()
  fixed <- screen_j_tie_heavy_r_trace(FALSE)
  repeated <- screen_j_tie_heavy_r_trace(TRUE)

  # Source trace: SKbias3.XYZ_bias_ANALYSE stores observed CHITOT in the
  # SINGLE-valued RESARRAY before Evaluate_simulated_biasresults compares the
  # simulated REAL CHI with >=. The 35 equal tables therefore all exceed.
  expect_gt(fixture$chi$stat, screen_j_source_single(fixture$chi$stat))
  expect_equal(fixed$p_chi, 1)
  expect_equal(fixed$chi_exceed, 40L)
  expect_equal(fixed$nsim, 40L)
  expect_equal(fixed$draw_count, 40L)
  expect_equal(fixed$final_seed, 460123585)
  expect_equal(sum(vapply(
    fixed$trajectory$tables,
    function(x) identical(as.numeric(x$table), c(6, 1, 1, 0)),
    logical(1L)
  )), 35L)
  expect_true(all(fixed$trajectory$simulations$chi_ge_observed))

  expect_equal(repeated$p_chi, 1)
  expect_equal(repeated$chi_exceed, 1L)
  expect_equal(repeated$nsim, 1L)
  expect_equal(repeated$draw_count, 1L)
  expect_equal(repeated$final_seed, 1212982318)
  expect_true(repeated$trajectory$simulations$chi_status[[1L]])
  expect_true(repeated$trajectory$simulations$gamma_status[[1L]])
  expect_true(repeated$trajectory$simulations$stop[[1L]])
})

test_that("native SCREEN J trace matches every R table, draw, statistic, and state", {
  fixture <- screen_j_tie_heavy_trace_fixture()
  for (sequential in c(FALSE, TRUE)) {
    reference <- screen_j_tie_heavy_r_trace(sequential)
    native <- screen_j_exact_chi_gamma_trace_slices_native(
      fixture$slices,
      fixture$chi$stat,
      fixture$gamma$gamma,
      nsim = 40L,
      seed = 9L,
      sequential = sequential,
      seq_limit = 1L
    )
    skip_if(is.null(native), "Native SCREEN J trace endpoint is unavailable.")

    aggregate_fields <- c(
      "p_chi", "p_gamma", "chi_exceed", "gamma_exceed", "nsim",
      "draw_count", "rng_draws", "final_seed"
    )
    expect_equal(native[aggregate_fields], reference[aggregate_fields], tolerance = 0)
    expect_equal(
      lapply(native$trajectory$tables, `[[`, "table"),
      lapply(reference$trajectory$tables, `[[`, "table"),
      tolerance = 0,
      ignore_attr = TRUE
    )
    expect_equal(
      lapply(native$trajectory$tables, `[[`, "random_draws"),
      lapply(reference$trajectory$tables, `[[`, "random_draws"),
      tolerance = 0
    )
    table_fields <- c(
      "sim", "slice", "chi_square", "ppq", "pmq", "draw_count", "final_seed"
    )
    for (field in table_fields) {
      expect_equal(
        vapply(native$trajectory$tables, `[[`, native$trajectory$tables[[1L]][[field]], field),
        vapply(reference$trajectory$tables, `[[`, reference$trajectory$tables[[1L]][[field]], field),
        tolerance = 0,
        info = field
      )
    }
    expect_equal(
      native$trajectory$simulations,
      reference$trajectory$simulations,
      tolerance = 0,
      ignore_attr = TRUE
    )
  }
})

test_that("SCREEN J repeated public chi path preserves source stopping", {
  fixture <- screen_j_tie_heavy_trace_fixture()
  reference <- screen_j_tie_heavy_r_trace(TRUE)
  actual <- screen_j_exact_partial_chi(
    fixture$strata,
    fixture$chi$stat,
    nsim = 40L,
    seed = 9L,
    sequential = TRUE,
    seq_limit = 1L,
    seq_p0 = 0.05,
    seq_boundary = 1.058
  )
  expect_equal(actual, reference$p_chi, tolerance = 0)
})
