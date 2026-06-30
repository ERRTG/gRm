test_that("exact command state follows DIGRAM source command defaults", {
  asymptotic <- gRm_exact_command_state("asymptotic")
  expect_equal(asymptotic$command_no, 17L)
  expect_false(asymptotic$exact)
  expect_false(asymptotic$sequential)
  expect_equal(asymptotic$nsim, 0L)
  expect_equal(asymptotic$seed, 9L)
  expect_equal(asymptotic$seq_limit, 0L)

  exact <- gRm_exact_command_state("exact")
  expect_equal(exact$command_no, 2L)
  expect_true(exact$exact)
  expect_false(exact$sequential)
  expect_equal(exact$nsim, 1000L)
  expect_equal(exact$seed, 9L)
  expect_equal(exact$seq_limit, 1000L)

  repeated <- gRm_exact_command_state("repeated")
  expect_equal(repeated$command_no, 74L)
  expect_true(repeated$exact)
  expect_true(repeated$sequential)
  expect_equal(repeated$nsim, 1000L)
  expect_equal(repeated$seed, 9L)
  expect_equal(repeated$seq_p0, 0.05)
  expect_equal(repeated$seq_alpha, 0.001)
  expect_equal(repeated$seq_limit, 20L)

  sequential <- gRm_exact_command_state("sequential")
  expect_equal(sequential$command_no, 57L)
  expect_true(sequential$exact)
  expect_true(sequential$sequential)
  expect_equal(sequential$nsim, 1000L)
  expect_equal(sequential$seed, 9L)
  expect_equal(sequential$seq_p0, 0.05)
  expect_equal(sequential$seq_alpha, 0)
  expect_equal(sequential$seq_b, 1000)
  expect_equal(sequential$seq_limit, 20L)
})

test_that("exact command state applies DIGRAM source parameter semantics", {
  exact <- gRm_exact_command_state("exact", nsim = 400L, seed = 11L)
  expect_equal(exact$nsim, 400L)
  expect_equal(exact$seed, 11L)
  expect_equal(exact$seq_limit, 400L)

  repeated <- gRm_exact_command_state(
    "repeated",
    nsim = 400L,
    seed = 11L,
    critlevel = 25L,
    risk = 2L
  )
  expect_equal(repeated$nsim, 400L)
  expect_equal(repeated$seed, 11L)
  expect_equal(repeated$seq_p0, 0.05)
  expect_equal(repeated$seq_alpha, 0.001)
  expect_equal(repeated$seq_limit, 20L)

  sequential <- gRm_exact_command_state(
    "sequential",
    nsim = 400L,
    seed = 11L,
    critlevel = 25L
  )
  expect_equal(sequential$nsim, 400L)
  expect_equal(sequential$seed, 11L)
  expect_equal(sequential$seq_p0, 0.05)
  expect_equal(sequential$seq_alpha, 0)
  expect_equal(sequential$seq_b, 1000)
  expect_equal(sequential$seq_limit, 25L)
})

test_that("exact command state rejects malformed scalar parameters", {
  bad_nonnegative <- list(1.9, c(25L, 999L), NA_integer_, Inf, -1L)
  for (x in bad_nonnegative) {
    expect_error(
      gRm_exact_command_state("exact", nsim = x),
      "`nsim` must be a single non-negative integer-like value"
    )
  }

  bad_seed <- list(9.9, c(9L, 10L), NA_integer_, Inf)
  for (x in bad_seed) {
    expect_error(
      gRm_exact_command_state("exact", seed = x),
      "`seed` must be a single integer-like value"
    )
  }

  bad_per1000 <- list(25.5, c(25L, 50L), NA_integer_, Inf, -1L)
  for (x in bad_per1000) {
    expect_error(
      gRm_exact_command_state("repeated", critlevel = x),
      "`critlevel` must be a single non-negative integer-like value"
    )
    expect_error(
      gRm_exact_command_state("repeated", risk = x),
      "`risk` must be a single non-negative integer-like value"
    )
  }

  expect_equal(gRm_exact_command_state("exact", nsim = 0L)$nsim, 1000L)
})

test_that("package does not expose a DIGRAM command-file parser", {
  expect_false(exists("parse_gRm_exact_command", envir = asNamespace("gRm"), inherits = FALSE))
})

test_that("source command state replaces artifact-specific branch cues", {
  repeated <- gRm_exact_command_state("repeated", nsim = 1000L, seed = 9L)
  sequential <- gRm_exact_command_state("sequential", nsim = 1000L, seed = 9L)

  project_a <- list(paths = list(input_dir = "/tmp/project_a"))
  project_b <- list(paths = list(input_dir = "/tmp/project_b"))
  expect_equal(screen_j_repeated_seq_limit(project_a, TRUE, 1000L, exact_state = repeated), 20L)
  expect_equal(screen_j_repeated_seq_limit(project_b, TRUE, 1000L, exact_state = repeated), 20L)
  expect_equal(screen_j_repeated_seq_limit(project_b, TRUE, 1000L, exact_state = sequential), 20L)

  y_labels <- list(backgrounds = data.frame(label_code = c("y", "z", "A", "B", "C", "D")))
  other_labels <- list(backgrounds = data.frame(label_code = c("I", "J", "K", "L", "M", "N")))
  expect_equal(exo_select_repeated_count_cutoff(y_labels, exact_state = repeated), 20L)
  expect_equal(exo_select_repeated_count_cutoff(other_labels, exact_state = repeated), 20L)
})
