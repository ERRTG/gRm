screen_j_opt_internal <- function(name, required = TRUE) {
  if (exists(name, envir = asNamespace("gRm"), mode = "function", inherits = FALSE)) {
    return(get(name, envir = asNamespace("gRm"), mode = "function"))
  }
  if (isTRUE(required)) {
    stop("Internal function is not available: ", name, call. = FALSE)
  }
  NULL
}

screen_j_opt_slices <- function() {
  list(
    matrix(c(7L, 2L, 1L, 0L, 3L, 4L), nrow = 3L),
    matrix(c(0L, 5L, 2L, 1L, 0L, 6L), nrow = 3L),
    matrix(c(3L, 0L, 1L, 2L, 8L, 1L), nrow = 3L),
    matrix(c(4L, 1L, 0L, 0L, 2L, 5L), nrow = 3L)
  )
}

screen_j_opt_advance_source_seed <- function(seed, draws) {
  state <- as.numeric(screen_j_opt_internal("screen_j_source_seed")(seed))
  base <- 65536
  multiplier_hi <- 2056
  multiplier_lo <- 33797
  for (draw in seq_len(as.integer(draws))) {
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

screen_j_opt_native <- function(slices, observed_chi, observed_gamma, nsim, seed, sequential, seq_limit) {
  screen_j_opt_internal("screen_j_exact_chi_gamma_slices_native")(
    slices,
    observed_chi,
    observed_gamma,
    nsim,
    seed,
    sequential = sequential,
    seq_limit = seq_limit
  )
}

screen_j_opt_reference <- function(slices, observed_chi, observed_gamma, nsim, seed, sequential, seq_limit) {
  prepared <- screen_j_opt_internal("screen_j_prepare_exact_slices")(slices)
  draw_count <- 0L
  random_draw <- screen_j_opt_internal("screen_j_source_random_stream")(seed)
  counted_draw <- function() {
    draw_count <<- draw_count + 1L
    random_draw()
  }
  out <- screen_j_opt_internal("screen_j_exact_chi_gamma_prepared_r")(
    prepared,
    observed_chi,
    observed_gamma,
    nsim,
    seed,
    sequential = sequential,
    random_draw = counted_draw,
    seq_limit = seq_limit,
    seq_p0 = 0.05,
    seq_boundary = 1.058
  )
  # Source trace: SKexa2.XYZ_TEST assigns NCHI / NSIM and NGAMMA / NSIM
  # directly to Real result cells; do not single-round the exact p-values here.
  out$draw_count <- draw_count
  out$final_seed <- screen_j_opt_advance_source_seed(seed, draw_count)
  out
}

screen_j_opt_expect_exact_equal <- function(native, reference) {
  expect_equal(native$p_chi, reference$p_chi, tolerance = 0)
  expect_equal(native$p_gamma, reference$p_gamma, tolerance = 0)
  expect_equal(native$nsim, reference$nsim)
  expect_equal(native$draw_count, reference$draw_count)
  expect_equal(native$rng_draws, reference$draw_count)
  expect_equal(native$final_seed, reference$final_seed)
}

screen_j_opt_project_data <- function() {
  data.frame(
    ID = seq_len(12L),
    I1 = c(1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
    I2 = c(1L, 2L, 1L, 3L, 2L, 3L, 3L, 1L, 2L, 2L, 3L, 1L),
    I3 = c(2L, 1L, 3L, 1L, 2L, 3L, 2L, 1L, 3L, 3L, 2L, 1L),
    site = c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 2L)
  )
}

screen_j_opt_with_native_env <- function(value, code) {
  old_cpp <- Sys.getenv("RDIGRAM_SCREEN_J_EXACT_CPP", unset = NA_character_)
  if (is.na(value)) {
    Sys.unsetenv("RDIGRAM_SCREEN_J_EXACT_CPP")
  } else {
    Sys.setenv(RDIGRAM_SCREEN_J_EXACT_CPP = value)
  }
  on.exit({
    if (is.na(old_cpp)) {
      Sys.unsetenv("RDIGRAM_SCREEN_J_EXACT_CPP")
    } else {
      Sys.setenv(RDIGRAM_SCREEN_J_EXACT_CPP = old_cpp)
    }
  }, add = TRUE)
  force(code)
}

test_that("optimized exact kernel preserves fixed exact p-values, draw count, and final seed", {
  slices <- screen_j_opt_slices()
  observed_chi <- 4.25
  observed_gamma <- 0.18
  native <- screen_j_opt_native(
    slices,
    observed_chi,
    observed_gamma,
    nsim = 53L,
    seed = 9L,
    sequential = FALSE,
    seq_limit = 53L
  )
  if (is.null(native)) {
    skip("No native SCREEN J exact chi-gamma wrapper is available.")
  }
  reference <- screen_j_opt_reference(
    slices,
    observed_chi,
    observed_gamma,
    nsim = 53L,
    seed = 9L,
    sequential = FALSE,
    seq_limit = 53L
  )

  screen_j_opt_expect_exact_equal(native, reference)
})

test_that("optimized exact kernel preserves repeated stopping, draw count, and final seed", {
  slices <- screen_j_opt_slices()
  observed_chi <- 0.25
  observed_gamma <- 0.01
  native <- screen_j_opt_native(
    slices,
    observed_chi,
    observed_gamma,
    nsim = 100L,
    seed = 17L,
    sequential = TRUE,
    seq_limit = 3L
  )
  if (is.null(native)) {
    skip("No native SCREEN J exact chi-gamma wrapper is available.")
  }
  reference <- screen_j_opt_reference(
    slices,
    observed_chi,
    observed_gamma,
    nsim = 100L,
    seed = 17L,
    sequential = TRUE,
    seq_limit = 3L
  )

  screen_j_opt_expect_exact_equal(native, reference)
})

test_that("two consecutive native exact calls are independent", {
  slices <- screen_j_opt_slices()
  first <- screen_j_opt_native(slices, 4.25, 0.18, 67L, 9L, FALSE, 67L)
  second <- screen_j_opt_native(slices, 4.25, 0.18, 67L, 9L, FALSE, 67L)
  if (is.null(first) || is.null(second)) {
    skip("No native SCREEN J exact chi-gamma wrapper is available.")
  }

  expect_equal(second$p_chi, first$p_chi, tolerance = 0)
  expect_equal(second$p_gamma, first$p_gamma, tolerance = 0)
  expect_equal(second$nsim, first$nsim)
  expect_equal(second$draw_count, first$draw_count)
  expect_equal(second$final_seed, first$final_seed)
})

test_that("conditional bias native and R fallback agree for exact and repeated", {
  helper <- screen_j_opt_internal("screen_j_conditional_bias_test")

  x <- c(1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 2L, 3L)
  y <- c(1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 1L)
  condition_values <- data.frame(score = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L))
  condition_dims <- 4L
  valid <- rep(TRUE, length(x))

  native_fixed <- helper(
    x, y, 3L, 2L, condition_values, condition_dims, valid,
    exact = TRUE, nsim = 37L, seed = 9L, repeated = FALSE,
    seq_limit = 37L, native = TRUE
  )
  r_fixed <- helper(
    x, y, 3L, 2L, condition_values, condition_dims, valid,
    exact = TRUE, nsim = 37L, seed = 9L, repeated = FALSE,
    seq_limit = 37L, native = FALSE
  )
  native_repeated <- helper(
    x, y, 3L, 2L, condition_values, condition_dims, valid,
    exact = TRUE, nsim = 80L, seed = 17L, repeated = TRUE,
    seq_limit = 3L, native = TRUE
  )
  r_repeated <- helper(
    x, y, 3L, 2L, condition_values, condition_dims, valid,
    exact = TRUE, nsim = 80L, seed = 17L, repeated = TRUE,
    seq_limit = 3L, native = FALSE
  )

  expect_equal(native_fixed$p_chi_exact, r_fixed$p_chi_exact, tolerance = 0)
  expect_equal(native_fixed$p_gamma_exact, r_fixed$p_gamma_exact, tolerance = 0)
  expect_equal(native_fixed$exact_nsim, r_fixed$exact_nsim)
  expect_equal(native_repeated$p_chi_exact, r_repeated$p_chi_exact, tolerance = 0)
  expect_equal(native_repeated$p_gamma_exact, r_repeated$p_gamma_exact, tolerance = 0)
  expect_equal(native_repeated$exact_nsim, r_repeated$exact_nsim)
})

test_that("item-pair native gate covers every statistic and deterministic state field", {
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
  strata <- screen_j_opt_internal("screen_j_strata_table")(
    item_matrix[, 1L],
    item_matrix[, 2L],
    rest_score + 1L,
    3L,
    2L,
    max(rest_score) + 1L,
    valid
  )
  partial_stats <- screen_j_opt_internal("screen_j_partial_gamma")(strata)
  partial_chi <- screen_j_opt_internal("screen_j_partial_chi")(strata)
  expect_true(is.finite(partial_stats$gamma))

  native_helper <- screen_j_opt_internal("screen_j_item_pair_conditional_exact_native", required = FALSE)
  if (is.null(native_helper)) {
    skip("Native item-pair conditional exact wrapper is not available.")
  }

  native <- native_helper(
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
  if (is.null(native)) {
    skip("Native item-pair conditional exact wrapper is disabled.")
  }
  reference <- screen_j_opt_internal("screen_j_item_pair_probe_reference")(
    list(
      x = item_matrix[, 1L],
      y = item_matrix[, 2L],
      x_dim = 3L,
      y_dim = 2L,
      condition_values = rest_score + 1L,
      condition_dim = max(rest_score) + 1L,
      valid = valid
    ),
    c(41L, 9L, 41L),
    FALSE
  )
  fields <- c(
    "chi_square", "df", "gamma", "ppq", "pmq", "s",
    "p_chi_asymp", "p_gamma_asymp", "p_chi_exact", "p_gamma_exact",
    "exact_nsim", "chi_exceed", "gamma_exceed", "draw_count", "final_seed"
  )

  expect_equal(native[fields], reference[fields], tolerance = 1e-12)
  expect_equal(native$chi_square, partial_chi$stat, tolerance = 0)
  expect_equal(native$gamma, partial_stats$gamma, tolerance = 0)
  expect_true(screen_j_opt_internal("screen_j_item_pair_native_source_faithful")())
})

test_that("native exact routing preserves selected SCREEN J model terms", {
  analysis <- gRm(
    screen_j_opt_project_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID"
  )
  values_native <- screen_j_opt_with_native_env(NA_character_, {
    screen_j_opt_internal("screen_j_values")(analysis$project, exact = TRUE, nsim = 23L, seed = 9L)
  })
  values_fallback <- screen_j_opt_with_native_env("false", {
    screen_j_opt_internal("screen_j_values")(analysis$project, exact = TRUE, nsim = 23L, seed = 9L)
  })

  expect_equal(values_native$partial$item_p, values_fallback$partial$item_p, tolerance = 1e-7)
  expect_equal(values_native$partial$exo_p, values_fallback$partial$exo_p, tolerance = 1e-7)
  expect_equal(values_native$model$local_dependence$matrix, values_fallback$model$local_dependence$matrix)
  expect_equal(values_native$model$local_dependence$rows, values_fallback$model$local_dependence$rows)
})
