active_fit_data <- function() {
  data.frame(
    ID = seq_len(12L),
    I1 = c(1L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L),
    I2 = c(1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L),
    I3 = c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 2L),
    X1 = c(1L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L)
  )
}

active_lid_data <- function() {
  rows <- expand.grid(
    I1 = 1:2,
    I2 = 1:2,
    I3 = 1:2,
    I4 = 1:2,
    X1 = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  data.frame(ID = seq_len(nrow(rows)), rows)
}

scalar_active_item_probabilities_reference <- function(context, state, total_score, background_values, components = NULL) {
  components <- components %||% context$ld_components_items %||% gllrm_ld_components(context)$items
  out <- lapply(context$item_raw_max, numeric)
  if (total_score < 0L || total_score > context$max_total_score) {
    return(out)
  }

  component_gamma <- lapply(components, function(component_items) {
    one <- gllrm_component_gamma(context, state, component_items, background_values)
    one$gamma * one$scale
  })
  convolutions <- gllrm_component_convolutions(component_gamma, context$max_total_score)
  denominator <- convolutions$full[[total_score + 1L]]
  if (denominator <= 0) {
    return(out)
  }

  for (component_index in seq_along(components)) {
    component_items <- components[[component_index]]
    rest_gamma <- convolutions$rest[[component_index]]
    configs <- gllrm_component_configurations(context, component_items)
    for (config_index in seq_len(nrow(configs))) {
      item_values <- as.integer(configs[config_index, ])
      names(item_values) <- names(configs)
      config_score <- sum(item_values)
      if (config_score > total_score) {
        next
      }
      rest_weight <- rest_gamma[[total_score - config_score + 1L]]
      if (rest_weight <= 0) {
        next
      }
      weight <- gllrm_component_config_weight(context, state, component_items, item_values, background_values)
      if (weight <= 0) {
        next
      }
      probability <- weight * rest_weight / denominator
      for (item_index in component_items) {
        item_score <- item_values[[as.character(item_index)]]
        out[[item_index]][[item_score + 1L]] <- out[[item_index]][[item_score + 1L]] + probability
      }
    }
  }

  out
}

scalar_active_ld_probabilities_reference <- function(context, state, total_score, background_values, components = NULL) {
  components <- components %||% context$ld_components_items %||% gllrm_ld_components(context)$items
  out <- lapply(context$ld_specs, function(spec) {
    matrix(0, nrow = context$item_raw_max[[spec$item1]], ncol = context$item_raw_max[[spec$item2]])
  })
  if (total_score < 0L || total_score > context$max_total_score) {
    return(out)
  }

  component_gamma <- lapply(components, function(component_items) {
    one <- gllrm_component_gamma(context, state, component_items, background_values)
    one$gamma * one$scale
  })
  convolutions <- gllrm_component_convolutions(component_gamma, context$max_total_score)
  denominator <- convolutions$full[[total_score + 1L]]
  if (denominator <= 0) {
    return(out)
  }

  for (component_index in seq_along(components)) {
    component_items <- components[[component_index]]
    key <- gllrm_component_key(component_items)
    ld_local <- context$component_ld_local_indices[[key]]
    if (is.null(ld_local) || nrow(ld_local) == 0L) {
      next
    }
    rest_gamma <- convolutions$rest[[component_index]]
    configs <- context$component_config_matrices[[key]]
    config_scores <- context$component_config_scores[[key]]
    for (config_index in seq_len(nrow(configs))) {
      item_values <- configs[config_index, ]
      config_score <- config_scores[[config_index]]
      if (config_score > total_score) {
        next
      }
      rest_weight <- rest_gamma[[total_score - config_score + 1L]]
      if (rest_weight <= 0) {
        next
      }
      weight <- gllrm_component_config_weight_fast(context, state, component_items, item_values, background_values, key = key)
      if (weight <= 0) {
        next
      }
      probability <- weight * rest_weight / denominator
      for (ld_row in seq_len(nrow(ld_local))) {
        ld_index <- ld_local$ld_index[[ld_row]]
        score1 <- item_values[[ld_local$item1_pos[[ld_row]]]] + 1L
        score2 <- item_values[[ld_local$item2_pos[[ld_row]]]] + 1L
        out[[ld_index]][score1, score2] <- out[[ld_index]][score1, score2] + probability
      }
    }
  }

  out
}

scalar_active_uniform_ld_scoregroup_tables_reference <- function(context, state, groups, spec, ld_index, components = NULL) {
  components <- components %||% context$ld_components_items %||% gllrm_ld_components(context)$items
  dimensions <- c(context$item_raw_max[[spec$item1]], context$item_raw_max[[spec$item2]], nrow(groups))
  observed <- array(0, dim = dimensions)
  expected <- array(0, dim = dimensions)

  score_group_index <- function(score) {
    hit <- which(groups$from_score <= score & score <= groups$to_score)
    if (length(hit) > 0L) hit[[1L]] else NA_integer_
  }

  for (row in context$valid_rows) {
    group_index <- score_group_index(context$score[[row]])
    if (is.na(group_index)) {
      next
    }
    score1 <- context$item_matrix[row, spec$item1] + 1L
    score2 <- context$item_matrix[row, spec$item2] + 1L
    observed[score1, score2, group_index] <- observed[score1, score2, group_index] + 1
  }

  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    score <- group$score[[1L]]
    homogeneity_group <- score_group_index(score)
    if (is.na(homogeneity_group)) {
      next
    }
    background_values <- gllrm_group_background_values(context, group)
    probabilities <- scalar_active_ld_probabilities_reference(
      context,
      state,
      total_score = score,
      background_values = background_values,
      components = components
    )[[ld_index]]
    expected[, , homogeneity_group] <- expected[, , homogeneity_group] +
      group$count[[1L]] * probabilities
  }

  list(observed = observed, expected = expected)
}

test_that("fit accepts active LD and DIF terms", {
  ia <- gRm(active_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  spec <- gllrm(ia, ld = ~ I1:I2, dif = ~ I3:X1)

  fit <- fit(spec, max_step = 200L, max_delta = 1e-6)

  expect_s3_class(fit, "gRm_fit")
  expect_equal(nrow(fit$spec$ld), 1L)
  expect_equal(nrow(fit$spec$dif), 1L)
  expect_true(is.list(fit$fit$ld_parameters))
  expect_true(is.list(fit$fit$dif_parameters))
  expect_true(isTRUE(fit$convergence$converged) || isFALSE(fit$convergence$converged))
  expect_false(any(grepl("not implemented", fit$unmodeled, fixed = TRUE)))
})

test_that("active expected-margin cache keys use only active DIF backgrounds", {
  context <- list(active_background_indices = 1L)

  expect_identical(
    gllrm_background_cache_key(context, c(X1 = 1L, X2 = 1L)),
    gllrm_background_cache_key(context, c(X1 = 1L, X2 = 2L))
  )
  expect_false(identical(
    gllrm_background_cache_key(context, c(X1 = 1L, X2 = 1L)),
    gllrm_background_cache_key(context, c(X1 = 2L, X2 = 1L))
  ))

  context$active_background_indices <- integer()
  expect_identical(gllrm_background_cache_key(context, c(X1 = 1L, X2 = 1L)), ".")
  expect_identical(gllrm_background_cache_key(context, c(X1 = 2L, X2 = 2L)), ".")
})

test_that("global homogeneity score group lookup matches source interval scan", {
  groups <- data.frame(
    group = 1:3,
    from_score = c(1L, 3L, 5L),
    to_score = c(2L, 4L, 7L)
  )

  lookup <- global_homogeneity_score_group_lookup(groups, max_score = 8L)
  scan_index <- function(score) {
    hit <- which(groups$from_score <= score & score <= groups$to_score)
    if (length(hit) > 0L) hit[[1L]] else NA_integer_
  }

  expect_equal(lookup[seq.int(0L, 8L) + 1L], vapply(seq.int(0L, 8L), scan_index, integer(1L)))
})

test_that("active fit summaries expose source-shaped parameter and expected-margin tables", {
  ia <- gRm(active_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I3:X1), max_step = 200L, max_delta = 1e-6)

  tables <- gllrm_active_detail_tables(fit$values)

  expect_true(all(c(
    "ld_parameters",
    "dif_parameters",
    "expected_items",
    "expected_ld",
    "expected_dif",
    "update_items",
    "update_ld",
    "update_dif"
  ) %in% names(tables)))
  expect_true(all(c("term", "score1", "score2", "gamma") %in% names(tables$ld_parameters)))
  expect_true(all(c("item", "exogenous", "score", "value", "gamma") %in% names(tables$dif_parameters)))
})

test_that("single LD general fit agrees with existing local-independence candidate fitter", {
  ia <- gRm(active_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  bundle <- build_item_parameters_bundle(ia$project)
  base_counts <- rasch_counts(bundle)
  candidate <- fit_active_ld_candidate(bundle, base_counts, item1 = 1L, item2 = 2L, max_step = 200L, max_delta = 1e-6)
  general <- fit(gllrm(ia, ld = ~ I1:I2), max_step = 200L, max_delta = 1e-6)

  expect_equal(general$fit$ld_parameters[[1L]], candidate$ld_gamma, tolerance = 1e-6)
})

test_that("single-pass active update matches source double-pass update", {
  ia <- gRm(active_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  spec <- gllrm(ia, ld = ~ I1:I2, dif = ~ I3:X1)
  context <- build_gllrm_active_context(spec, build_item_parameters_bundle(ia$project))
  state <- initialize_gllrm_active_state(context)
  state <- calculate_gllrm_joint_expected_margins(context, state)

  double_pass <- state
  double_pass$delta <- 0
  double_pass <- update_gllrm_active_parameters(context, double_pass, apply_update = FALSE, track_delta = TRUE)
  double_pass_report_delta <- double_pass$delta
  double_pass <- update_gllrm_active_parameters(context, double_pass, apply_update = TRUE, track_delta = TRUE)

  single_pass <- state
  single_pass$delta <- 0
  single_pass <- update_gllrm_active_parameters_once(context, single_pass)

  expect_equal(single_pass$report_delta, double_pass_report_delta, tolerance = 1e-12)
  expect_equal(single_pass$delta, double_pass$delta, tolerance = 1e-12)
  expect_equal(single_pass$update_items, double_pass$update_items, tolerance = 1e-12)
  expect_equal(single_pass$update_ld, double_pass$update_ld, tolerance = 1e-12)
  expect_equal(single_pass$update_dif, double_pass$update_dif, tolerance = 1e-12)
  expect_equal(single_pass$item_gamma, double_pass$item_gamma, tolerance = 1e-12)
  expect_equal(single_pass$ld_parameters, double_pass$ld_parameters, tolerance = 1e-12)
  expect_equal(single_pass$dif_parameters, double_pass$dif_parameters, tolerance = 1e-12)
})

test_that("active LD reference adjustment reselects source densest ties per table", {
  observed <- matrix(
    c(
      1, 1, 0,
      1, 1, 0,
      0, 0, 1
    ),
    nrow = 3L,
    byrow = TRUE
  )
  gamma <- matrix(
    c(
      2, 3, 0,
      5, 7, 0,
      0, 0, 11
    ),
    nrow = 3L,
    byrow = TRUE
  )

  carried <- adjust_ld_gamma_source_reference_details(
    observed,
    gamma,
    i_ref = 2L,
    j_ref = 2L,
    preserve_current_ties = TRUE
  )
  source_table <- adjust_ld_gamma_source_reference_details(
    observed,
    gamma,
    i_ref = 2L,
    j_ref = 2L,
    preserve_current_ties = FALSE
  )

  expect_equal(carried$i_ref, 2L)
  expect_equal(carried$j_ref, 2L)
  expect_equal(source_table$i_ref, 1L)
  expect_equal(source_table$j_ref, 1L)
  expect_false(isTRUE(all.equal(carried$adjusted, source_table$adjusted)))
})

test_that("active dependency adjustment applies source tie reselection after carried references", {
  context <- list(
    ld_specs = list(
      list(item1 = 1L, item2 = 2L),
      list(item1 = 2L, item2 = 3L)
    ),
    dif_specs = list(),
    observed_ld = list(
      matrix(
        c(
          1, 0, 0,
          1, 1, 0,
          0, 0, 1
        ),
        nrow = 3L,
        byrow = TRUE
      ),
      matrix(
        c(
          1, 1, 0,
          1, 1, 0,
          0, 0, 1
        ),
        nrow = 3L,
        byrow = TRUE
      )
    ),
    bundle = list(
      model = list(
        items = data.frame(
          name = c("I1", "I2", "I3"),
          raw_max = c(3L, 3L, 3L),
          stringsAsFactors = FALSE
        )
      )
    )
  )
  state <- list(
    item_gamma = matrix(
      1,
      nrow = 3L,
      ncol = 3L,
      dimnames = list(c("I1", "I2", "I3"), as.character(0:2))
    ),
    ld_parameters = list(
      matrix(c(2, 3, 0, 5, 7, 0, 0, 0, 11), nrow = 3L, byrow = TRUE),
      matrix(c(13, 17, 0, 19, 23, 0, 0, 0, 29), nrow = 3L, byrow = TRUE)
    ),
    dif_parameters = list()
  )

  carried <- adjust_gllrm_dependency_parameters(context, state, preserve_current_ties = TRUE)
  source_table <- adjust_gllrm_dependency_parameters(context, state, preserve_current_ties = FALSE)

  expect_false(isTRUE(all.equal(carried$ld_parameters[[2L]], source_table$ld_parameters[[2L]])))
  expect_equal(
    source_table$ld_parameters[[2L]],
    adjust_ld_gamma_source_reference_details(
      context$observed_ld[[2L]],
      state$ld_parameters[[2L]],
      i_ref = 2L,
      j_ref = 2L,
      preserve_current_ties = FALSE
    )$adjusted
  )
})

test_that("active component-restscore tables seed complete-item endpoints", {
  context <- list(
    item_raw_max = c(4L, 3L),
    max_total_score = 5L,
    item_matrix = matrix(
      c(
        0L, 0L,
        0L, 1L,
        1L, 0L,
        2L, 2L,
        2L, -1L
      ),
      ncol = 2L,
      byrow = TRUE
    ),
    score = c(0L, 1L, 1L, 4L, -1L),
    valid_rows = c(2L, 3L),
    ld_components_items = list(2L),
    score_exo_groups = data.frame(score = integer(), count = integer())
  )

  tables <- active_gllrm_component_restscore_tables(context, state = list(), component_items = 1L)

  expect_equal(tables$observed[1L, 1L], 1)
  expect_equal(tables$observed[4L, 3L], 0)
  expect_equal(tables$expected[1L, 1L], 1)
  expect_equal(tables$expected[4L, 3L], 0)
  expect_equal(tables$observed[1L, 2L], 1)
  expect_equal(tables$observed[2L, 1L], 1)
})

test_that("active component-restscore tables keep source interior score loop", {
  context <- list(
    item_raw_max = c(4L, 3L),
    max_total_score = 5L,
    item_matrix = matrix(
      c(
        0L, 0L,
        0L, 1L,
        3L, 2L
      ),
      ncol = 2L,
      byrow = TRUE
    ),
    score = c(0L, 1L, 5L),
    valid_rows = 1:3,
    ld_components_items = list(2L),
    score_exo_groups = data.frame(score = integer(), count = integer())
  )

  tables <- active_gllrm_component_restscore_tables(context, state = list(), component_items = 1L)

  expect_equal(tables$observed[1L, 1L], 1)
  expect_equal(tables$observed[4L, 3L], 1)
  expect_equal(tables$observed[1L, 2L], 1)
})

test_that("single DIF general fit agrees with existing missing-DIF candidate fitter", {
  ia <- gRm(active_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  bundle <- build_item_parameters_bundle(ia$project)
  base_counts <- rasch_counts(bundle)
  candidate <- fit_active_dif_candidate(bundle, base_counts, target_item = 3L, background_index = 1L, max_step = 200L, max_delta = 1e-6)
  general <- fit(gllrm(ia, dif = ~ I3:X1), max_step = 200L, max_delta = 1e-6)

  expect_equal(general$fit$dif_parameters[[1L]], candidate$dif_gamma, tolerance = 1e-6)
})

test_that("active GLLRM item-parameter details are available from the fit object", {
  ia <- gRm(active_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I3:X1), max_step = 200L, max_delta = 1e-6)
  terms <- model_terms(fit)
  tables <- gllrm_active_detail_tables(fit$values)

  expect_equal(terms$ld$item1, "I1")
  expect_equal(terms$ld$item2, "I2")
  expect_equal(terms$dif$item, "I3")
  expect_equal(terms$dif$exogenous, "X1")
  expect_true(all(c("ld_parameters", "dif_parameters") %in% names(tables)))
  expect_gt(nrow(tables$ld_parameters), 0L)
  expect_gt(nrow(tables$dif_parameters), 0L)
  expect_equal(fit$values$n_parameters, 4L)
})

test_that("active GLLRM local-independence checks condition on the current model", {
  ia <- gRm(active_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)

  values <- local_independence_values(fit, max_step = 80L, max_delta = 1e-5, jobs = 2L)

  expect_s3_class(values, "gRm_local_independence_values")
  expect_equal(sort(values$tests$pair_label), c("ac", "ad", "bc", "bd", "cd"))
  expect_false(any(values$tests$pair_label == "ab"))
  expect_equal(sort(values$tests$item1_name), c("I1", "I1", "I2", "I2", "I3"))
  expect_equal(sort(values$tests$item2_name), c("I3", "I3", "I4", "I4", "I4"))
})

test_that("active GLLRM DIF checks use the current model", {
  ia <- gRm(active_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)

  values <- dif_tests_values(fit, max_step = 80L, max_delta = 1e-5, jobs = 2L)

  expect_s3_class(values, "gRm_dif_tests_values")
  expect_equal(sort(paste0(values$tests$item_label, values$tests$background_label)), c("ae", "be", "ce"))
  expect_equal(paste0(values$active_tests$item_label, values$active_tests$background_label), "de")
  expect_true("active_tests" %in% names(details(values)$tables))
  expect_false(any(values$tests$item_name == "I4" & values$tests$background_name == "X1"))
})

test_that("active GLLRM global homogeneity refits the active model in score groups", {
  ia <- gRm(active_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)

  values <- global_homogeneity_values(fit, score_cuts = c(2L, 4L), max_step = 80L, max_delta = 1e-5)

  expect_s3_class(values, "gRm_global_homogeneity_values")
  expect_equal(values$summary$n_parameters, fit$values$n_parameters)
  expect_equal(values$summary$full_log_likelihood, fit$values$log_likelihood, tolerance = 1e-8)
  expect_equal(nrow(values$score_groups), 2L)
  expect_true(all(vapply(values$group_values, function(x) inherits(x$fit$values, "gRm_active_gllrm_values"), logical(1L))))
  expect_true(all(is.na(values$items$residual)))
  expect_false(any(values$items$residual_runtime_source_backed))
  expect_false(any(values$items$marker_runtime_source_backed))
  expect_true(all(c("summary", "groups", "items") %in% detail_names(values)))
})

test_that("batched uniform LD scoregroup tables match scalar source-faithful references", {
  ia <- gRm(active_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2 + I2:I3, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)
  context <- fit$fit$context
  state <- fit$fit
  groups <- global_homogeneity_score_groups(context$bundle, c(2L, 4L))
  components <- context$ld_components_items

  batched <- active_uniform_ld_scoregroup_tables_all(
    context,
    groups,
    new_active_gllrm_probability_cache(context, state, components = components)
  )

  for (ld_index in seq_along(context$ld_specs)) {
    expected <- scalar_active_uniform_ld_scoregroup_tables_reference(
      context,
      state,
      groups,
      context$ld_specs[[ld_index]],
      ld_index,
      components
    )
    expect_equal(batched[[ld_index]], expected, tolerance = 1e-12)
  }
})

test_that("active GLLRM item-fit values use active conditional probabilities", {
  ia <- gRm(active_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)

  values <- item_fits_values(fit, include_extended = TRUE)

  expect_s3_class(values, "gRm_item_fits_values")
  expect_equal(values$fit$model, "active_gllrm")
  expect_true(all(c("outfit", "infit", "observed_gamma", "expected_gamma") %in% names(values$items)))
  expect_false(is.null(values$extended))
  expect_true(all(c("statistics", "score_level_fit", "item_fit_summaries") %in% detail_names(values)))
})

test_that("active probability helpers match scalar source-faithful references", {
  ia <- gRm(active_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2 + I2:I3, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)
  context <- fit$fit$context
  state <- fit$fit
  components <- context$ld_components_items

  checked <- 0L
  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    total_score <- group$score[[1L]]
    background_values <- gllrm_group_background_values(context, group)

    item_expected <- scalar_active_item_probabilities_reference(context, state, total_score, background_values, components)
    item_actual <- active_gllrm_group_item_probabilities(context, state, total_score, background_values, components)
    expect_equal(item_actual, item_expected, tolerance = 1e-12)

    ld_expected <- scalar_active_ld_probabilities_reference(context, state, total_score, background_values, components)
    ld_actual <- active_gllrm_group_ld_probabilities(context, state, total_score, background_values, components)
    expect_equal(ld_actual, ld_expected, tolerance = 1e-12)

    checked <- checked + 1L
    if (checked >= 4L) {
      break
    }
  }
  expect_equal(checked, 4L)

  item_fit_values <- item_fits_values(fit, include_extended = TRUE)
  expect_s3_class(item_fit_values, "gRm_item_fits_values")
  expect_true(all(is.finite(item_fit_values$items$expected_gamma)))

  gh_values <- global_homogeneity_values(fit, score_cuts = c(2L, 4L), max_step = 80L, max_delta = 1e-5)
  expect_s3_class(gh_values, "gRm_global_homogeneity_values")
  expect_true(is.finite(gh_values$summary$clr))
})
