convergence_contract_analysis <- function() {
  data <- expand.grid(
    I1 = 0:1,
    I2 = 0:1,
    I3 = 0:1,
    X = 0:1,
    KEEP.OUT.ATTRS = FALSE
  )
  data$ID <- seq_len(nrow(data))
  gRm(
    data[c("ID", "I1", "I2", "I3", "X")],
    items = c("I1", "I2", "I3"),
    exogenous = "X",
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1),
    exogenous_levels = list(X = 0:1),
    score_cuts = c(1L, 3L)
  )
}

test_that("Rasch, GLLRM, and candidate fits share one convergence schema", {
  analysis <- convergence_contract_analysis()
  rasch <- fit(gllrm(analysis), max_step = 50L, max_delta = 1e-5)
  general <- fit(
    gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:X),
    max_step = 50L,
    max_delta = 1e-5
  )
  candidate <- fit_gllrm_candidate_ld(
    rasch,
    item1 = 1L,
    item2 = 2L,
    max_step = 50L,
    max_delta = 1e-5
  )

  expected_names <- c(
    "schema", "converged", "source_converged",
    "source_converged_before_post_acceptance", "iterations", "report_delta",
    "final_delta", "tolerance", "max_step", "stop_reason",
    "post_stop_accepted", "attempted", "report_value_source", "delta",
    "max_delta"
  )
  for (state in list(rasch$convergence, general$convergence, candidate$convergence)) {
    expect_s3_class(state, "gRm_convergence_state")
    expect_named(state, expected_names)
    expect_identical(state$schema, "gRm-convergence-state-v1")
    expect_identical(state$source_converged, state$converged)
    expect_identical(state$delta, state$report_delta)
    expect_identical(state$max_delta, state$tolerance)
    expect_identical(state$report_value_source, "attempted_fit")
    expect_identical(state$attempted$iterations, state$iterations)
    expect_identical(state$attempted$report_delta, state$report_delta)
    expect_identical(state$attempted$final_delta, state$final_delta)
    expect_identical(state$attempted$converged, state$converged)
    expect_identical(state$attempted$stop_reason, state$stop_reason)
  }

  expect_identical(rasch$convergence$report_delta, rasch$fit$report_delta)
  expect_identical(rasch$convergence$final_delta, rasch$fit$delta)
  expect_identical(general$convergence$report_delta, general$fit$report_delta)
  expect_identical(general$convergence$final_delta, general$fit$delta)
  expect_identical(candidate$convergence$report_delta, candidate$state$report_delta)
  expect_identical(candidate$convergence$final_delta, candidate$state$delta)
})

test_that("GLLRM post-stop acceptance is explicit in convergence metadata", {
  state <- list(
    converged = TRUE,
    convergence_before_final_acceptance = FALSE,
    n_step = 17L,
    report_delta = 0.05,
    delta = 0.04,
    stop_reason = "max_step"
  )
  convergence <- gRm_gllrm_convergence_state(state, max_step = 17L, max_delta = 1e-4)

  expect_true(convergence$converged)
  expect_false(convergence$source_converged_before_post_acceptance)
  expect_true(convergence$post_stop_accepted)
  expect_equal(convergence$tolerance, 1e-4)
})
