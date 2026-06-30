screen_j_internal <- function(name) {
  get(name, envir = asNamespace("gRm"), mode = "function")
}

screen_j_prepare_exact_slices <- screen_j_internal("screen_j_prepare_exact_slices")
screen_j_source_random_stream <- screen_j_internal("screen_j_source_random_stream")
screen_j_source_seed <- screen_j_internal("screen_j_source_seed")
screen_j_exact_native_available <- screen_j_internal("screen_j_exact_native_available")
screen_j_conditional_native_allowed <- screen_j_internal("screen_j_conditional_native_allowed")
screen_j_conditional_native_source_faithful <- screen_j_internal("screen_j_conditional_native_source_faithful")
screen_j_conditional_try_native <- screen_j_internal("screen_j_conditional_try_native")
screen_j_conditional_bias_test <- screen_j_internal("screen_j_conditional_bias_test")
screen_j_conditional_bias_test_native <- screen_j_internal("screen_j_conditional_bias_test_native")
screen_j_rc_chi_square_prepared_expected <- screen_j_internal("screen_j_rc_chi_square_prepared_expected")
screen_rc_gamma_counts <- screen_j_internal("screen_rc_gamma_counts")
exo_select_gentab1_prepared <- screen_j_internal("exo_select_gentab1_prepared")

screen_j_native_exact_get <- function(names) {
  for (name in names) {
    if (exists(name, mode = "function", inherits = TRUE)) {
      return(list(name = name, fun = get(name, mode = "function", inherits = TRUE)))
    }
    if ("gRm" %in% loadedNamespaces() &&
        exists(name, envir = asNamespace("gRm"), mode = "function", inherits = FALSE)) {
      return(list(name = name, fun = get(name, envir = asNamespace("gRm"), mode = "function")))
    }
  }
  NULL
}

screen_j_native_exact_fused <- function() {
  screen_j_native_exact_get(c(
    "screen_j_exact_chi_gamma_slices_cpp",
    "screen_j_exact_chi_gamma_slices_native",
    "screen_j_exact_chi_gamma_kernel_cpp",
    "screen_j_exact_chi_gamma_prepared_cpp"
  ))
}

screen_j_native_exact_gamma <- function() {
  screen_j_native_exact_get(c(
    "screen_j_exact_gamma_slices_cpp",
    "screen_j_exact_gamma_slices_native",
    "screen_j_exact_gamma_kernel_cpp",
    "screen_j_exact_gamma_prepared_cpp"
  ))
}

screen_j_native_exact_call <- function(candidate, slices, observed_chi, observed_gamma, nsim, seed, kind) {
  use_prepared <- grepl("prepared|kernel", candidate$name)
  first_arg <- if (use_prepared) screen_j_prepare_exact_slices(slices) else slices
  args <- if (identical(kind, "gamma")) {
    list(first_arg, observed_gamma, nsim, seed)
  } else {
    list(first_arg, observed_chi, observed_gamma, nsim, seed)
  }
  fn_formals <- names(formals(candidate$fun))
  if (identical(kind, "fused") && "sequential" %in% fn_formals) {
    args$sequential <- TRUE
  }
  draw_flags <- c("return_draw_count", "return_rng_draws", "include_draw_count", "count_draws")
  for (flag in intersect(draw_flags, fn_formals)) {
    args[[flag]] <- TRUE
  }
  result <- do.call(candidate$fun, args)
  screen_j_native_exact_normalize(result, kind = kind)
}

screen_j_native_exact_first_name <- function(candidates, available) {
  matched <- intersect(candidates, available)
  if (length(matched)) matched[[1L]] else character()
}

screen_j_native_exact_normalize <- function(result, kind) {
  draw_names <- c("rng_draws", "draw_count", "n_draws", "draws", "random_draws")
  seed_names <- c("final_seed", "rng_seed", "source_seed")
  nsim_names <- c("nsim", "n_sim", "completed", "completed_nsim")
  p_gamma_names <- c("p_gamma", "gamma_p", "p_gamma_exact", "gamma", "p")
  p_chi_names <- c("p_chi", "chi_p", "p_chi_exact", "chi")

  if (is.list(result)) {
    draw_name <- screen_j_native_exact_first_name(draw_names, names(result))
    seed_name <- screen_j_native_exact_first_name(seed_names, names(result))
    nsim_name <- screen_j_native_exact_first_name(nsim_names, names(result))
    p_gamma_name <- screen_j_native_exact_first_name(p_gamma_names, names(result))
    p_chi_name <- screen_j_native_exact_first_name(p_chi_names, names(result))
    return(list(
      p_chi = if (length(p_chi_name)) unname(result[[p_chi_name]]) else NA_real_,
      p_gamma = if (length(p_gamma_name)) unname(result[[p_gamma_name]]) else NA_real_,
      rng_draws = if (length(draw_name)) as.integer(result[[draw_name]]) else NA_integer_,
      final_seed = if (length(seed_name)) as.numeric(result[[seed_name]]) else NA_real_,
      nsim = if (length(nsim_name)) as.integer(result[[nsim_name]]) else NA_integer_
    ))
  }

  numeric_result <- unname(result)
  result_names <- names(result)
  draw <- attr(result, "rng_draws", exact = TRUE)
  if (is.null(draw)) draw <- attr(result, "draw_count", exact = TRUE)
  final_seed <- attr(result, "final_seed", exact = TRUE)
  nsim <- attr(result, "nsim", exact = TRUE)
  if (is.null(draw) && !is.null(result_names)) {
    draw_index <- match(TRUE, result_names %in% draw_names)
    draw <- if (!is.na(draw_index)) numeric_result[[draw_index]] else NA_integer_
  }
  if (is.null(final_seed) && !is.null(result_names)) {
    seed_index <- match(TRUE, result_names %in% seed_names)
    final_seed <- if (!is.na(seed_index)) numeric_result[[seed_index]] else NA_real_
  }
  if (is.null(nsim) && !is.null(result_names)) {
    nsim_index <- match(TRUE, result_names %in% nsim_names)
    nsim <- if (!is.na(nsim_index)) numeric_result[[nsim_index]] else NA_integer_
  }

  p_chi <- NA_real_
  p_gamma <- NA_real_
  if (!is.null(result_names)) {
    chi_index <- match(TRUE, result_names %in% p_chi_names)
    gamma_index <- match(TRUE, result_names %in% p_gamma_names)
    if (!is.na(chi_index)) p_chi <- numeric_result[[chi_index]]
    if (!is.na(gamma_index)) p_gamma <- numeric_result[[gamma_index]]
  } else if (identical(kind, "gamma") && length(numeric_result) >= 1L) {
    p_gamma <- numeric_result[[1L]]
  }
  list(
    p_chi = p_chi,
    p_gamma = p_gamma,
    rng_draws = as.integer(draw),
    final_seed = as.numeric(final_seed),
    nsim = as.integer(nsim)
  )
}

screen_j_native_exact_final_seed <- function(seed, draws) {
  state <- as.numeric(screen_j_source_seed(seed))
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

screen_j_reference_exact <- function(slices, observed_chi, observed_gamma, nsim, seed, kind) {
  prepared_slices <- screen_j_prepare_exact_slices(slices)
  random_draw <- screen_j_source_random_stream(seed)
  draw_count <- 0L
  counted_draw <- function() {
    draw_count <<- draw_count + 1L
    random_draw()
  }
  chi_exceed <- 0L
  gamma_exceed <- 0L
  chi_status <- !identical(kind, "fused")
  gamma_status <- FALSE
  completed <- nsim

  for (sim in seq_len(nsim)) {
    chi_total <- 0
    ppq_total <- 0
    pmq_total <- 0
    for (prepared in prepared_slices) {
      generated <- exo_select_gentab1_prepared(prepared, random_draw = counted_draw)
      if (identical(kind, "fused")) {
        chi_total <- chi_total + screen_j_rc_chi_square_prepared_expected(generated, prepared)
      }
      gamma_counts <- screen_rc_gamma_counts(generated)
      ppq_total <- ppq_total + gamma_counts$ppq
      pmq_total <- pmq_total + gamma_counts$pmq
    }
    if (identical(kind, "fused") && chi_total >= observed_chi) {
      chi_exceed <- chi_exceed + 1L
    }
    simulated_gamma <- if (ppq_total > 0) pmq_total / ppq_total else 0
    if (abs(simulated_gamma) >= abs(observed_gamma)) {
      gamma_exceed <- gamma_exceed + 1L
    }
    if (identical(kind, "fused")) {
      if (sim >= 21L) {
        root <- sqrt(sim)
        if (chi_exceed / root - root * 0.05 >= 1.058) {
          chi_status <- TRUE
        }
        if (gamma_exceed / root - root * 0.05 >= 1.058) {
          gamma_status <- TRUE
        }
      }
      if (chi_status && gamma_status) {
        completed <- sim
        break
      }
    }
  }

  list(
    p_chi = if (identical(kind, "fused")) chi_exceed / completed else NA_real_,
    p_gamma = gamma_exceed / completed,
    rng_draws = draw_count,
    final_seed = screen_j_native_exact_final_seed(seed, draw_count),
    nsim = completed
  )
}

screen_j_representative_exact_slices <- function() {
  list(
    matrix(c(4L, 1L, 2L, 3L), nrow = 2L),
    matrix(c(0L, 2L, 1L, 5L), nrow = 2L),
    matrix(c(3L, 0L, 1L, 2L), nrow = 2L)
  )
}

screen_j_conditional_native_fixture <- function() {
  list(
    x = c(1L, 2L, 1L, 2L, 3L, 3L, 1L, 1L, 2L, 3L, 2L, 3L),
    y = c(1L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L),
    x_dim = 3L,
    y_dim = 2L,
    condition_values = matrix(c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 3L, 3L, 3L, 3L), ncol = 1L),
    condition_dims = 3L,
    valid = rep(TRUE, 12L)
  )
}

expect_screen_j_conditional_native_parity <- function(reference, native) {
  fields <- c("chi_square", "df", "gamma", "p_chi", "p_gamma", "exact_nsim")
  expect_equal(native[fields], reference[fields], tolerance = 0)
  expect_equal(native$p_chi_exact, reference$p_chi_exact, tolerance = 0)
  expect_equal(native$p_gamma_exact, reference$p_gamma_exact, tolerance = 0)
}

test_that("SCREEN J native exact gamma kernel matches R simulation and draw count", {
  candidate <- screen_j_native_exact_gamma()
  if (is.null(candidate)) {
    skip("No native SCREEN J exact gamma wrapper is available.")
  }

  slices <- screen_j_representative_exact_slices()
  observed_gamma <- 0.2
  nsim <- 29L
  seed <- 9L
  native <- screen_j_native_exact_call(candidate, slices, NA_real_, observed_gamma, nsim, seed, "gamma")
  reference <- screen_j_reference_exact(slices, NA_real_, observed_gamma, nsim, seed, "gamma")

  expect_false(is.na(native$rng_draws), info = paste(candidate$name, "must report RNG draw count"))
  expect_false(is.na(native$final_seed), info = paste(candidate$name, "must report final source seed"))
  expect_false(is.na(native$nsim), info = paste(candidate$name, "must report completed simulation count"))
  expect_equal(native$rng_draws, reference$rng_draws)
  expect_equal(native$final_seed, reference$final_seed)
  expect_equal(native$nsim, reference$nsim)
  expect_equal(native$p_gamma, reference$p_gamma, tolerance = 1e-7)
})

test_that("SCREEN J native exact chi-gamma kernel matches R simulation and draw count", {
  candidate <- screen_j_native_exact_fused()
  if (is.null(candidate)) {
    skip("No native SCREEN J exact chi-gamma wrapper is available.")
  }

  slices <- screen_j_representative_exact_slices()
  observed_chi <- 1.25
  observed_gamma <- 0.2
  nsim <- 31L
  seed <- 17L
  native <- screen_j_native_exact_call(candidate, slices, observed_chi, observed_gamma, nsim, seed, "fused")
  reference <- screen_j_reference_exact(slices, observed_chi, observed_gamma, nsim, seed, "fused")

  expect_false(is.na(native$rng_draws), info = paste(candidate$name, "must report RNG draw count"))
  expect_false(is.na(native$final_seed), info = paste(candidate$name, "must report final source seed"))
  expect_false(is.na(native$nsim), info = paste(candidate$name, "must report completed simulation count"))
  expect_equal(native$rng_draws, reference$rng_draws)
  expect_equal(native$final_seed, reference$final_seed)
  expect_equal(native$nsim, reference$nsim)
  expect_equal(native$p_chi, reference$p_chi, tolerance = 1e-7)
  expect_equal(native$p_gamma, reference$p_gamma, tolerance = 1e-7)
})

test_that("SCREEN J native scalar exact attributes are explicitly protected", {
  path <- package_path("src", "screen_j_exact.cpp")
  skip_if(
    !file.exists(path),
    "source-tree C++ file is not installed"
  )
  cpp <- readLines(path, warn = FALSE)

  expect_true(any(grepl("SEXP exceed_attr = PROTECT", cpp, fixed = TRUE)))
  expect_true(any(grepl("SEXP final_seed_attr = PROTECT", cpp, fixed = TRUE)))
  expect_true(any(grepl("UNPROTECT(6)", cpp, fixed = TRUE)))
})

test_that("SCREEN J conditional native path matches R exact reference", {
  if (!screen_j_exact_native_available()) {
    skip("No native SCREEN J conditional wrapper is available.")
  }
  fixture <- screen_j_conditional_native_fixture()

  reference <- do.call(
    screen_j_conditional_bias_test,
    c(fixture, list(exact = TRUE, nsim = 200L, seed = 123L, repeated = FALSE, native = FALSE))
  )
  native <- do.call(
    screen_j_conditional_bias_test_native,
    c(fixture, list(exact = TRUE, nsim = 200L, seed = 123L, repeated = FALSE))
  )

  expect_screen_j_conditional_native_parity(reference, native)
})

test_that("SCREEN J conditional native path matches R repeated exact reference", {
  if (!screen_j_exact_native_available()) {
    skip("No native SCREEN J conditional wrapper is available.")
  }
  fixture <- screen_j_conditional_native_fixture()

  reference <- do.call(
    screen_j_conditional_bias_test,
    c(
      fixture,
      list(
        exact = TRUE, nsim = 200L, seed = 123L, repeated = TRUE,
        native = FALSE, seq_limit = 200L, seq_p0 = 0.05, seq_boundary = 1.058
      )
    )
  )
  native <- do.call(
    screen_j_conditional_bias_test_native,
    c(fixture, list(exact = TRUE, nsim = 200L, seed = 123L, repeated = TRUE, seq_limit = 200L))
  )

  expect_screen_j_conditional_native_parity(reference, native)
})

test_that("SCREEN J conditional native route is parity gated", {
  if (!screen_j_exact_native_available()) {
    skip("No native SCREEN J conditional wrapper is available.")
  }
  fixture <- screen_j_conditional_native_fixture()
  source_faithful <- screen_j_conditional_native_source_faithful()

  result <- do.call(
    screen_j_conditional_try_native,
    c(
      fixture,
      list(
        exact = TRUE, nsim = 200L, seed = 123L, repeated = FALSE,
        random_draw = NULL, seq_limit = 200L,
        use_native = source_faithful
      )
    )
  )

  if (isTRUE(source_faithful)) {
    expect_false(is.null(result))
  } else {
    expect_null(result)
  }
})

test_that("SCREEN J conditional native source-faithfulness gate opens", {
  if (!screen_j_exact_native_available()) {
    skip("No native SCREEN J conditional wrapper is available.")
  }

  expect_true(screen_j_conditional_native_source_faithful())
})

test_that("SCREEN J conditional native route honors the native disable flag", {
  withr::local_envvar(RDIGRAM_SCREEN_J_EXACT_CPP = "false")

  expect_false(screen_j_exact_native_available())
  expect_false(screen_j_conditional_native_allowed(FALSE, 0.05, 1.058))
})

test_that("SCREEN J exact score screening confidence limits accept vector rows", {
  screen_j_score_report_rows <- screen_j_internal("screen_j_score_report_rows")

  marginal_rows <- data.frame(
    label = "H",
    name = "Age_PF",
    chi_square = 1,
    df = 1L,
    p_chi = 0.2,
    p_chi_asymp = 0.2,
    p_chi_exact = 0.2,
    p_chi_exact_low = 0.1,
    p_chi_exact_high = 0.3,
    gamma = 0.1,
    p_gamma = 0.4,
    p_gamma_asymp = 0.4,
    p_gamma_exact = 0.4,
    p_gamma_exact_low = 0.3,
    p_gamma_exact_high = 0.5,
    exact_nsim = 1000L,
    selected = FALSE,
    stringsAsFactors = FALSE
  )
  screening_rows <- data.frame(
    candidate_label = c("I", "J"),
    candidate_name = c("Education", "Etiology"),
    chi_square = c(2, 3),
    df = c(1L, 1L),
    p_chi = c(0.05, 0.06),
    p_chi_asymp = c(0.05, 0.06),
    p_chi_exact = c(0.04, 0.07),
    gamma = c(0.2, -0.3),
    p_gamma = c(0.08, 0.09),
    p_gamma_asymp = c(0.08, 0.09),
    p_gamma_exact = c(0.09, 0.1),
    n = c(200L, 250L),
    hypothesis = c("#&H|I", "#&H|J"),
    stringsAsFactors = FALSE
  )

  rows <- screen_j_score_report_rows(marginal_rows, screening_rows, exact = TRUE)

  expect_equal(nrow(rows), 3L)
  expect_false(anyNA(rows$p_chi_exact_low))
  expect_false(anyNA(rows$p_gamma_exact_high))
})
