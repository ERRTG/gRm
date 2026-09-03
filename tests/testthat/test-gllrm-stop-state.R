run_source_delta_sequence <- function(deltas,
                                      n_valid = 100L,
                                      max_step = length(deltas) + 1L,
                                      max_delta = 0.0001) {
  control <- source_gllrm_control_state(n_valid, max_step, max_delta)
  decision <- NULL
  for (delta in deltas) {
    observed <- source_gllrm_observe_delta(control, delta)
    control <- observed$control
    decision <- observed$decision
    if (isTRUE(decision$stop)) {
      break
    }
  }
  list(control = control, decision = decision)
}

gllrm_recurrence_fixture_data <- function() {
  data.frame(
    I1 = c(1,1,1,0,0,1,0,1,0,1,1,1,1,1,0,1,1,1,1,0,1,1,1,0,0,0,1,1,0,0,0,0,1,0,0),
    I2 = c(1,1,0,0,0,1,1,1,0,0,1,1,1,0,0,1,1,1,0,0,1,0,1,1,0,0,0,1,0,1,0,1,1,0,0),
    I3 = c(1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,1,1,1,0,0,1,1,1,1,0,0,1,1,0,0),
    I4 = c(1,1,1,0,1,0,1,1,0,0,1,0,1,0,0,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,1,1,1,0,0),
    I5 = c(0,1,1,0,0,1,0,1,1,0,1,1,1,0,0,1,1,1,0,0,1,0,1,0,1,1,1,1,1,0,1,1,1,0,0),
    X1 = c(3,3,1,2,3,3,3,3,3,3,3,2,1,2,3,3,3,1,2,1,1,3,2,3,1,3,2,1,2,1,2,1,3,2,3),
    X2 = c(2,1,2,1,2,1,1,1,2,1,1,1,1,1,1,1,2,2,2,2,2,1,1,1,1,1,1,1,1,2,2,2,1,2,2),
    ID = seq_len(35L)
  )
}

test_that("GLLRM stop control uses source initialization", {
  control <- source_gllrm_control_state(37L, 5000L, 0.0001)

  expect_equal(control$delta, 37)
  expect_true(control$convergence)
  expect_equal(control$n_step, 0L)
  expect_equal(control$n_finish, 0L)
  expect_equal(unname(control$delta_history), rep(9999, 5L))
  expect_equal(names(control$delta_history), as.character(-4:0))
})

test_that("zero score reference never triggers heuristic output gauges", {
  state <- list(marker = "unchanged")
  contexts <- list(
    list(
      item_score_reference = 0L,
      ld_specs = rep(list(list()), 3L),
      ld_components_items = list(1:4)
    ),
    list(
      item_score_reference = 0L,
      ld_specs = rep(list(list()), 4L),
      ld_components_items = list(1:5)
    )
  )

  for (context in contexts) {
    expect_identical(gllrm_output_parameter_state(context, state), state)
  }
})

test_that("GLLRM stop control preserves strict tolerance and branch order", {
  tolerance <- run_source_delta_sequence(0.00001, max_delta = 0.0001)
  expect_true(tolerance$decision$stop)
  expect_identical(tolerance$decision$reason, "delta_below_tolerance")
  expect_true(tolerance$decision$convergence)

  boundary <- run_source_delta_sequence(
    c(0.1, 0.09),
    max_step = 3L,
    max_delta = 0.1
  )
  expect_equal(boundary$control$n_step, 2L)
  expect_identical(boundary$decision$reason, "delta_below_tolerance")

  max_step_boundary <- run_source_delta_sequence(
    0.1,
    max_step = 1L,
    max_delta = 0.1
  )
  expect_identical(max_step_boundary$decision$reason, "max_step")
  expect_true(max_step_boundary$decision$convergence)
})

test_that("GLLRM stop control covers equality and periodic branches", {
  repeated <- run_source_delta_sequence(7, n_valid = 7L, max_step = 20L)
  expect_identical(repeated$decision$reason, "repeated_delta")
  expect_false(repeated$decision$convergence)

  two_back <- run_source_delta_sequence(c(9, 8, 7, 6, 5, 6), max_step = 20L)
  expect_identical(two_back$decision$reason, "two_back_repeated_delta")
  expect_equal(two_back$control$n_step, 6L)

  periodic_50 <- run_source_delta_sequence(seq(60, 11, length.out = 50L), max_step = 100L)
  expect_identical(periodic_50$decision$reason, "periodic_50_large_delta")
  expect_equal(periodic_50$control$n_step, 50L)

  periodic_1000 <- run_source_delta_sequence(
    1 + 1 / (seq_len(1000L) + 1),
    max_step = 2000L
  )
  expect_identical(periodic_1000$decision$reason, "periodic_1000")
  expect_equal(periodic_1000$control$n_step, 1000L)
})

test_that("GLLRM stop control covers recurrence and plateau branches", {
  recurrence <- run_source_delta_sequence(
    c(seq(6, 1, length.out = 50L), 2),
    max_step = 200L
  )
  expect_identical(recurrence$decision$reason, "recurring_delta_values")
  expect_equal(recurrence$control$n_step, 51L)
  expect_true(recurrence$control$recurring)
  expect_false(recurrence$decision$convergence)

  plateau <- run_source_delta_sequence(
    c(1, seq(1.1, 2.1, length.out = 11L)),
    max_step = 100L
  )
  expect_identical(plateau$decision$reason, "finish_count_plateau")
  expect_equal(plateau$control$n_finish, 11L)
  expect_false(plateau$decision$convergence)
})

test_that("GLLRM final acceptance uses strict source delta below 0.1", {
  accepted <- run_source_delta_sequence(
    0.05,
    max_step = 1L,
    max_delta = 0.0001
  )
  expect_false(accepted$decision$convergence)
  expect_true(source_gllrm_final_convergence(accepted$control))

  rejected <- run_source_delta_sequence(
    0.1,
    max_step = 1L,
    max_delta = 0.0001
  )
  expect_false(rejected$decision$convergence)
  expect_false(source_gllrm_final_convergence(rejected$control))
})

test_that("a fitted sparse GLLRM excludes extremes and stops deterministically", {
  analysis <- gRm(
    gllrm_recurrence_fixture_data(),
    items = paste0("I", 1:5),
    exogenous = c("X1", "X2"),
    id = "ID"
  )
  model <- gllrm(
    analysis,
    ld = ~ I1:I2 + I2:I3 + I3:I4 + I4:I5 + I1:I5,
    dif = ~ I1:X1 + I3:X2 + I5:X1
  )
  # This intentionally sparse fixture exercises the complete fitted stop path;
  # the source window excludes both the zero and maximum possible scores. The
  # synthetic recurrence sequences above retain direct coverage of Pascal's
  # RecurringDeltaValues branch.
  fitted <- fit_gllrm(model, max_step = 250L, max_delta = 1e-300)
  state <- c(fitted$state, list(context = fitted$context))

  expect_identical(state$stop_reason, "repeated_delta")
  expect_equal(state$n_step, 211L)
  expect_false(state$recurring_delta_values)
  expect_false(state$convergence_before_final_acceptance)
  expect_true(state$converged)
  expect_lt(state$report_delta, 1e-12)
  expect_equal(state$context$counts$n_valid, 22L)
  expect_false(any(state$context$score[state$context$valid_rows] %in% c(0L, 5L)))

  expect_equal(state$context$item_score_reference, 0L)
  output_state <- gllrm_output_parameter_state(state$context, state)
  expect_equal(output_state$item_gamma, state$item_gamma)
  expect_equal(output_state$ld_parameters, state$ld_parameters)
  expect_equal(
    gllrm_loglike(state$context, output_state),
    gllrm_loglike(state$context, state),
    tolerance = 1e-10
  )
})
