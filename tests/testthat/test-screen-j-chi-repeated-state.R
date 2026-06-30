test_that("SCREEN J exact chi path honors repeated command state", {
  # Source trace: DIGRAM1f REP sets SEQUENTIAL, SEQ_P0, SEQ_B, and seq_limit.
  # SKbias3.XYZ_bias_ANALYSE then evaluates simulated chi values through
  # Evaluate_simulated_biasresults() and Stop_exact_biastests().
  old_cpp <- Sys.getenv("RDIGRAM_SCREEN_J_EXACT_CPP", unset = NA_character_)
  Sys.setenv(RDIGRAM_SCREEN_J_EXACT_CPP = "false")
  on.exit({
    if (is.na(old_cpp)) {
      Sys.unsetenv("RDIGRAM_SCREEN_J_EXACT_CPP")
    } else {
      Sys.setenv(RDIGRAM_SCREEN_J_EXACT_CPP = old_cpp)
    }
  }, add = TRUE)

  strata <- array(as.integer(c(6L, 1L, 1L, 0L)), dim = c(2L, 2L, 1L))
  observed_chi <- screen_j_partial_chi(strata)$stat

  fixed <- screen_j_exact_partial_chi(
    strata,
    observed_chi,
    nsim = 40L,
    seed = 9L
  )
  repeated <- screen_j_exact_partial_chi(
    strata,
    observed_chi,
    nsim = 40L,
    seed = 9L,
    sequential = TRUE,
    seq_limit = 1L,
    seq_p0 = 0.05,
    seq_boundary = 1.058
  )

  expect_equal(fixed, 0.125, tolerance = 1e-12)
  expect_equal(repeated, 1 / 3, tolerance = 1e-12)
})

screen_j_repeated_chi_reference <- function(strata, observed_chi, nsim, seed, seq_limit) {
  prepared_slices <- screen_j_prepare_exact_slices(screen_j_strata_slices(strata))
  random_draw <- screen_j_source_random_stream(seed)
  draw_count <- 0L
  counted_draw <- function() {
    draw_count <<- draw_count + 1L
    random_draw()
  }
  exceed <- 0L
  completed <- nsim
  for (sim in seq_len(nsim)) {
    chi_total <- 0
    for (prepared in prepared_slices) {
      generated <- exo_select_gentab1_prepared(prepared, random_draw = counted_draw)
      chi_total <- chi_total + screen_j_rc_chi_square_prepared_expected(generated, prepared)
    }
    if (chi_total >= observed_chi) {
      exceed <- exceed + 1L
    }
    if (screen_j_seq_t(exceed, sim, 0.05) >= 1.058 || exceed >= seq_limit) {
      completed <- sim
      break
    }
  }
  list(
    p_chi = exceed / completed,
    nsim = completed,
    draw_count = draw_count,
    final_seed = {
      state <- as.numeric(screen_j_source_seed(seed))
      base <- 65536
      multiplier_hi <- 2056
      multiplier_lo <- 33797
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
  )
}

test_that("SCREEN J exact chi repeated public path preserves source stopping with native available", {
  strata <- array(as.integer(c(6L, 1L, 1L, 0L)), dim = c(2L, 2L, 1L))
  observed_chi <- screen_j_partial_chi(strata)$stat
  actual <- screen_j_exact_partial_chi(
    strata,
    observed_chi,
    nsim = 40L,
    seed = 9L,
    sequential = TRUE,
    seq_limit = 1L,
    seq_p0 = 0.05,
    seq_boundary = 1.058
  )
  reference <- screen_j_repeated_chi_reference(
    strata,
    observed_chi,
    nsim = 40L,
    seed = 9L,
    seq_limit = 1L
  )

  expect_equal(actual, reference$p_chi, tolerance = 0)
})
