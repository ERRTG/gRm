gllrm_context_data <- function() {
  data.frame(
    ID = 1:6,
    I1 = c(1L, 2L, 1L, 2L, 1L, 2L),
    I2 = c(1L, 1L, 2L, 2L, 1L, 2L),
    I3 = c(1L, 1L, 2L, 1L, 2L, 2L),
    X1 = c(1L, 1L, 1L, 2L, 2L, 2L)
  )
}

test_that("GLLRM context stores terms in source position order", {
  ia <- gRm(gllrm_context_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  spec <- gllrm(ia, ld = ~ I2:I1, dif = ~ I3:X1)
  bundle <- build_item_parameters_bundle(ia$project)

  context <- build_gllrm_context(spec, bundle)

  expect_equal(context$ld_specs[[1L]]$item1, 1L)
  expect_equal(context$ld_specs[[1L]]$item2, 2L)
  expect_equal(context$dif_specs[[1L]]$item, 3L)
  expect_equal(context$dif_specs[[1L]]$background, 1L)
  expect_equal(context$item_score_values[[1L]], 0:1)
  expect_equal(context$background_values[[1L]], 1:2)
})

test_that("GLLRM context counts observed LD, DIF, and first-seen score/exogenous groups", {
  ia <- gRm(gllrm_context_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  spec <- gllrm(ia, ld = ~ I1:I2, dif = ~ I3:X1)
  context <- build_gllrm_context(spec, build_item_parameters_bundle(ia$project))

  expect_equal(context$observed_ld[[1L]][1L, 1L], 1L)
  expect_equal(context$observed_dif[[1L]][1L, 1L], 1L)
  expect_equal(context$score_exo_groups$score, c(1L, 2L, 2L, 1L, 3L))
  expect_equal(context$score_exo_groups$count, c(1L, 1L, 1L, 1L, 1L))
  expect_equal(context$score_exo_groups$X1, c(1L, 1L, 2L, 2L, 2L))
})

test_that("score/exogenous groups ignore inincluded DIF backgrounds in the grouping key", {
  context <- list(
    valid_rows = 1:4,
    n_backgrounds = 2L,
    backgrounds = data.frame(name = c("X_dif", "X_other")),
    background_matrix = matrix(
      c(
        1L, 1L,
        1L, 2L,
        2L, 1L,
        2L, 2L
      ),
      ncol = 2L,
      byrow = TRUE
    ),
    score = c(2L, 2L, 2L, 2L),
    dif_background_indices = 1L
  )

  groups <- gllrm_score_exo_groups(context)

  expect_equal(groups$score, c(2L, 2L))
  expect_equal(groups$count, c(2L, 2L))
  expect_equal(groups$X_dif, c(1L, 2L))
})

test_that("GLLRM components preserve source item order and initialize state", {
  ia <- gRm(gllrm_context_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  spec <- gllrm(ia, ld = ~ I1:I2 + I2:I3, dif = ~ I1:X1)
  context <- build_gllrm_context(spec, build_item_parameters_bundle(ia$project))

  components <- gllrm_ld_components(context)
  state <- initialize_gllrm_state(context)

  expect_equal(components$items, list(1:3))
  expect_equal(length(state$ld_parameters), 2L)
  expect_equal(length(state$dif_parameters), 1L)
  expect_true(all(state$item_gamma[, 1:2] == 1))
})

test_that("GLLRM context stores matrix fast paths matching data-frame metadata", {
  ia <- gRm(gllrm_context_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  spec <- gllrm(ia, ld = ~ I1:I2 + I2:I3, dif = ~ I1:X1 + I3:X1)
  context <- build_gllrm_context(spec, build_item_parameters_bundle(ia$project))

  expect_equal(context$dif_by_item_matrices[[1L]], unname(as.matrix(context$dif_by_item[[1L]])))
  expect_equal(context$dif_by_item_matrices[[2L]], matrix(integer(), ncol = 2L))
  expect_equal(context$dif_by_item_matrices[[3L]], unname(as.matrix(context$dif_by_item[[3L]])))

  key <- gllrm_component_key(context$ld_components_items[[1L]])
  expect_equal(
    context$component_ld_local_matrices[[key]],
    unname(as.matrix(context$component_ld_local_indices[[key]]))
  )
})

test_that("GLLRM component complexity limit is described as an R implementation guard", {
  analysis <- gRm(gllrm_context_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  spec <- gllrm(analysis, ld = ~ I1:I2 + I2:I3)
  bundle <- build_item_parameters_bundle(analysis$project)

  expect_error(
    build_gllrm_context(spec, bundle, max_joint_configs = 1L),
    "R implementation guard",
    fixed = TRUE
  )
})

test_that("superseded non-context candidate diagnostic helpers are not retained", {
  namespace <- asNamespace("gRm")
  superseded_helpers <- c(
    "candidate_ld_counts",
    "calculate_candidate_ld_expected",
    "candidate_ld_loglike",
    "candidate_dif_counts",
    "calculate_candidate_dif_expected",
    "candidate_dif_loglike"
  )

  expect_false(any(vapply(
    superseded_helpers,
    exists,
    logical(1L),
    envir = namespace,
    inherits = FALSE
  )))
})
