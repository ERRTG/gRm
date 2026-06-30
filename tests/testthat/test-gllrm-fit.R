gllrm_fit_data <- function() {
  data.frame(
    ID = seq_len(12L),
    I1 = c(1L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L),
    I2 = c(1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L),
    I3 = c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 2L),
    X1 = c(1L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L)
  )
}

non_syntactic_gllrm_fit_data <- function() {
  rows <- expand.grid(
    `item:one` = 0:1,
    `item two` = 0:1,
    `item-three` = 0:1,
    `site one` = 0:1,
    KEEP.OUT.ATTRS = FALSE
  )
  rows <- data.frame(rows, check.names = FALSE)
  rows$ID <- seq_len(nrow(rows))
  rows[, c("ID", "item:one", "item two", "item-three", "site one")]
}

gllrm_lid_data <- function() {
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

scalar_gllrm_item_probabilities_reference <- function(context, state, total_score, background_values, components = NULL) {
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

scalar_gllrm_ld_probabilities_reference <- function(context, state, total_score, background_values, components = NULL) {
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

scalar_gllrm_uniform_ld_scoregroup_tables_reference <- function(context, state, groups, spec, ld_index, components = NULL) {
  components <- components %||% context$ld_components_items %||% gllrm_ld_components(context)$items
  dimensions <- c(context$item_raw_max[[spec$item1]], context$item_raw_max[[spec$item2]], nrow(groups))
  observed <- array(0, dim = dimensions)
  expected <- array(0, dim = dimensions)

  score_group_lookup <- global_homogeneity_uniform_score_group_lookup(groups, context$max_total_score)
  uniform_rows <- gllrm_uniform_complete_rows(context, score_group_lookup)

  for (row in uniform_rows) {
    group_index <- global_homogeneity_lookup_score(score_group_lookup, context$score[[row]])
    if (is.na(group_index)) {
      next
    }
    score1 <- context$item_matrix[row, spec$item1] + 1L
    score2 <- context$item_matrix[row, spec$item2] + 1L
    observed[score1, score2, group_index] <- observed[score1, score2, group_index] + 1
  }

  score_exo_groups <- gllrm_score_exo_groups(context, rows = uniform_rows)
  for (group_index in seq_len(nrow(score_exo_groups))) {
    group <- score_exo_groups[group_index, , drop = FALSE]
    score <- group$score[[1L]]
    homogeneity_group <- global_homogeneity_lookup_score(score_group_lookup, score)
    if (is.na(homogeneity_group)) {
      next
    }
    background_values <- gllrm_group_background_values(context, group)
    probabilities <- scalar_gllrm_ld_probabilities_reference(
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

scalar_gllrm_uniform_dif_scoregroup_tables_reference <- function(context, state, groups, spec, components = NULL) {
  components <- components %||% context$ld_components_items %||% gllrm_ld_components(context)$items
  item_max <- context$item_raw_max[[spec$item]] - 1L
  background_max <- context$background_raw_max[[spec$background]]
  dimensions <- c(item_max + 1L, background_max, nrow(groups))
  observed <- array(0, dim = dimensions)
  expected <- array(0, dim = dimensions)

  score_group_lookup <- global_homogeneity_uniform_score_group_lookup(groups, context$max_total_score)
  uniform_rows <- gllrm_uniform_complete_rows(context, score_group_lookup)

  for (row in uniform_rows) {
    group_index <- global_homogeneity_lookup_score(score_group_lookup, context$score[[row]])
    if (is.na(group_index)) {
      next
    }
    item_score <- context$item_matrix[row, spec$item] + 1L
    background_value <- context$background_matrix[row, spec$background]
    observed[item_score, background_value, group_index] <-
      observed[item_score, background_value, group_index] + 1
  }

  score_exo_groups <- gllrm_score_exo_groups(context, rows = uniform_rows)
  for (group_index in seq_len(nrow(score_exo_groups))) {
    group <- score_exo_groups[group_index, , drop = FALSE]
    score <- group$score[[1L]]
    homogeneity_group <- global_homogeneity_lookup_score(score_group_lookup, score)
    if (is.na(homogeneity_group)) {
      next
    }
    background_values <- gllrm_group_background_values(context, group)
    background_value <- background_values[[spec$background]]
    probabilities <- scalar_gllrm_item_probabilities_reference(
      context,
      state,
      total_score = score,
      background_values = background_values,
      components = components
    )[[spec$item]]
    expected[seq_along(probabilities), background_value, homogeneity_group] <-
      expected[seq_along(probabilities), background_value, homogeneity_group] +
        group$count[[1L]] * probabilities
  }

  list(observed = observed, expected = expected)
}

test_that("fit accepts included LD and DIF terms", {
  ia <- gRm(gllrm_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  spec <- gllrm(ia, ld = ~ I1:I2, dif = ~ I3:X1)

  fit <- fit(spec, max_step = 200L, max_delta = 1e-6)

  expect_s3_class(fit, "gRm_fit")
  expect_equal(nrow(fit$spec$ld), 1L)
  expect_equal(nrow(fit$spec$dif), 1L)
  expect_true(is.list(fit$fit$ld_parameters))
  expect_true(is.list(fit$fit$dif_parameters))
  expect_true(isTRUE(fit$convergence$converged) || isFALSE(fit$convergence$converged))
  expect_false(any(grepl("not implemented", fit$unmodeled, fixed = TRUE)))

  limited <- fit_gllrm(spec, max_step = 1L, max_delta = 0)
  expect_equal(limited$state$stop_reason, "max_step")
  expect_equal(limited$state$n_step, 1L)
})

test_that("GLLRM diagnostics preserve non-syntactic LD and DIF names without formula roundtrips", {
  data <- non_syntactic_gllrm_fit_data()
  analysis <- gRm(
    data,
    items = c("item:one", "item two", "item-three"),
    exogenous = "site one",
    id = "ID",
    item_levels = list(
      `item:one` = 0:1,
      `item two` = 0:1,
      `item-three` = 0:1
    ),
    exogenous_levels = list(`site one` = 0:1),
    score_cuts = c(1L, 3L)
  )
  model <- gllrm(
    analysis,
    ld = ~ `item:one`:`item two`,
    dif = ~ `item-three`:`site one`
  )
  fitted <- fit(model, max_step = 50L, max_delta = 1e-4)

  ld <- local_dependence(fitted, max_step = 25L, max_delta = 1e-4, jobs = 1L)
  dif_result <- dif(fitted, max_step = 25L, max_delta = 1e-4, jobs = 1L)

  expect_s3_class(ld, "gRm_local_dependence")
  expect_s3_class(dif_result, "gRm_dif")
  expect_setequal(unlist(fitted$model$ld[1L, c("item1", "item2")]), c("item:one", "item two"))
  expect_equal(fitted$model$dif$item[[1L]], "item-three")
  expect_equal(fitted$model$dif$exogenous[[1L]], "site one")
})

test_that("GLLRM stores source negative log likelihood and logLik returns R sign", {
  analysis <- gRm(
    gllrm_fit_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "X1",
    id = "ID"
  )
  fitted <- fit(gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:X1), max_step = 100L, max_delta = 1e-4)

  expect_true(is.numeric(fitted$values$log_likelihood))
  expect_equal(as.numeric(logLik(fitted)), -fitted$values$log_likelihood)
  expect_equal(as.numeric(logLik(fitted)), -fitted$fit$log_likelihood)
})

test_that("GLLRM expected-margin cache keys use only included DIF backgrounds", {
  context <- list(dif_background_indices = 1L)

  expect_identical(
    gllrm_background_cache_key(context, c(X1 = 1L, X2 = 1L)),
    gllrm_background_cache_key(context, c(X1 = 1L, X2 = 2L))
  )
  expect_false(identical(
    gllrm_background_cache_key(context, c(X1 = 1L, X2 = 1L)),
    gllrm_background_cache_key(context, c(X1 = 2L, X2 = 1L))
  ))

  context$dif_background_indices <- integer()
  expect_identical(gllrm_background_cache_key(context, c(X1 = 1L, X2 = 1L)), ".")
  expect_identical(gllrm_background_cache_key(context, c(X1 = 2L, X2 = 2L)), ".")
})

test_that("native GLLRM expected-margin backend is registered", {
  expect_true(gllrm_expected_native_available())
})

test_that("native GLLRM expected margins match R reference at initial state", {
  ia <- gRm(gllrm_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  spec <- gllrm(ia, ld = ~ I1:I2 + I2:I3, dif = ~ I3:X1 + I4:X1)
  context <- build_gllrm_context(spec, build_item_parameters_bundle(ia$project))
  state <- initialize_gllrm_state(context)

  expect_gllrm_expected_matches_reference(context, state)
})

test_that("native GLLRM expected margins match R reference after updates", {
  ia <- gRm(gllrm_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  spec <- gllrm(ia, ld = ~ I1:I2 + I2:I3, dif = ~ I3:X1 + I4:X1)
  context <- build_gllrm_context(spec, build_item_parameters_bundle(ia$project))
  state <- initialize_gllrm_state(context)
  state <- calculate_gllrm_joint_expected_margins_r(context, state)
  state <- update_gllrm_parameters_once(context, state)

  expect_gllrm_expected_matches_reference(context, state)
})

test_that("GLLRM expected-margin dispatcher matches R reference", {
  ia <- gRm(gllrm_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  spec <- gllrm(ia, ld = ~ I1:I2, dif = ~ I4:X1)
  context <- build_gllrm_context(spec, build_item_parameters_bundle(ia$project))
  state <- initialize_gllrm_state(context)

  reference <- calculate_gllrm_joint_expected_margins_r(context, state)
  dispatched <- calculate_gllrm_joint_expected_margins(context, state)

  expect_equal(dispatched$expected_items, reference$expected_items, tolerance = 1e-10)
  expect_equal(dispatched$expected_ld, reference$expected_ld, tolerance = 1e-10)
  expect_equal(dispatched$expected_dif, reference$expected_dif, tolerance = 1e-10)
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

test_that("GLLRM fit summaries expose source-shaped parameter and expected-margin tables", {
  ia <- gRm(gllrm_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I3:X1), max_step = 200L, max_delta = 1e-6)

  tables <- gllrm_detail_tables(fit$values)

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

test_that("LD reference adjustment scans source-style strict improvements", {
  observed <- matrix(0, nrow = 4L, ncol = 4L)
  observed[1:3, 1L] <- 1
  observed[1:4, 2L] <- 1
  observed[1:4, 3L] <- 1
  observed[1:2, 4L] <- 1

  adjustment <- adjust_ld_gamma_source_reference_details(
    observed,
    matrix(1, nrow = 4L, ncol = 4L),
    i_ref = 4L,
    j_ref = 4L,
    preserve_current_ties = TRUE
  )

  expect_equal(adjustment$i_ref, 1L)
  expect_equal(adjustment$j_ref, 2L)
})

test_that("single LD general fit agrees with existing local-independence candidate fitter", {
  ia <- gRm(gllrm_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  bundle <- build_item_parameters_bundle(ia$project)
  base_counts <- rasch_counts(bundle)
  candidate <- fit_ld_candidate(bundle, base_counts, item1 = 1L, item2 = 2L, max_step = 200L, max_delta = 1e-6)
  general <- fit(gllrm(ia, ld = ~ I1:I2), max_step = 200L, max_delta = 1e-6)

  expect_equal(general$fit$ld_parameters[[1L]], candidate$ld_gamma, tolerance = 1e-6)
})

test_that("single-pass GLLRM update matches source double-pass update", {
  ia <- gRm(gllrm_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  spec <- gllrm(ia, ld = ~ I1:I2, dif = ~ I3:X1)
  context <- build_gllrm_context(spec, build_item_parameters_bundle(ia$project))
  state <- initialize_gllrm_state(context)
  state <- calculate_gllrm_joint_expected_margins(context, state)

  double_pass <- state
  double_pass$delta <- 0
  double_pass <- update_gllrm_parameters(context, double_pass, apply_update = FALSE, track_delta = TRUE)
  double_pass_report_delta <- double_pass$delta
  double_pass <- update_gllrm_parameters(context, double_pass, apply_update = TRUE, track_delta = TRUE)

  single_pass <- state
  single_pass$delta <- 0
  single_pass <- update_gllrm_parameters_once(context, single_pass)

  expect_equal(single_pass$report_delta, double_pass_report_delta, tolerance = 1e-12)
  expect_equal(single_pass$delta, double_pass$delta, tolerance = 1e-12)
  expect_equal(single_pass$update_items, double_pass$update_items, tolerance = 1e-12)
  expect_equal(single_pass$update_ld, double_pass$update_ld, tolerance = 1e-12)
  expect_equal(single_pass$update_dif, double_pass$update_dif, tolerance = 1e-12)
  expect_equal(single_pass$item_gamma, double_pass$item_gamma, tolerance = 1e-12)
  expect_equal(single_pass$ld_parameters, double_pass$ld_parameters, tolerance = 1e-12)
  expect_equal(single_pass$dif_parameters, double_pass$dif_parameters, tolerance = 1e-12)
})

test_that("included LD reference adjustment reselects source densest ties per table", {
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

test_that("GLLRM dependency adjustment applies source tie reselection after carried references", {
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

test_that("GLLRM component-restscore tables seed complete-item endpoints", {
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

  tables <- gllrm_component_restscore_tables(context, state = list(), component_items = 1L)

  expect_equal(tables$observed[1L, 1L], 1)
  expect_equal(tables$observed[4L, 3L], 0)
  expect_equal(tables$expected[1L, 1L], 1)
  expect_equal(tables$expected[4L, 3L], 0)
  expect_equal(tables$observed[1L, 2L], 1)
  expect_equal(tables$observed[2L, 1L], 1)
})

test_that("GLLRM component-restscore tables keep source interior score loop", {
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

  tables <- gllrm_component_restscore_tables(context, state = list(), component_items = 1L)

  expect_equal(tables$observed[1L, 1L], 1)
  expect_equal(tables$observed[4L, 3L], 1)
  expect_equal(tables$observed[1L, 2L], 1)
})

test_that("single DIF general fit agrees with existing missing-DIF candidate fitter", {
  ia <- gRm(gllrm_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  bundle <- build_item_parameters_bundle(ia$project)
  base_counts <- rasch_counts(bundle)
  candidate <- fit_dif_candidate(bundle, base_counts, target_item = 3L, background_index = 1L, max_step = 200L, max_delta = 1e-6)
  general <- fit(gllrm(ia, dif = ~ I3:X1), max_step = 200L, max_delta = 1e-6)

  expect_equal(general$fit$dif_parameters[[1L]], candidate$dif_gamma, tolerance = 1e-6)
})

test_that("GLLRM item-parameter value tables are available from the fit object", {
  ia <- gRm(gllrm_fit_data(), items = c("I1", "I2", "I3"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I3:X1), max_step = 200L, max_delta = 1e-6)
  terms <- model_terms(fit)
  tables <- gllrm_detail_tables(fit$values)

  expect_equal(terms$ld$item1, "I1")
  expect_equal(terms$ld$item2, "I2")
  expect_equal(terms$dif$item, "I3")
  expect_equal(terms$dif$exogenous, "X1")
  expect_true(all(c("ld_parameters", "dif_parameters") %in% names(tables)))
  expect_gt(nrow(tables$ld_parameters), 0L)
  expect_gt(nrow(tables$dif_parameters), 0L)
  expect_equal(fit$values$n_parameters, 4L)
})

test_that("GLLRM local-independence checks condition on the current model", {
  ia <- gRm(gllrm_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)

  values <- local_independence_values(fit, max_step = 80L, max_delta = 1e-5, jobs = 2L)

  expect_s3_class(values, "gRm_local_independence_values")
  expect_equal(sort(values$tests$pair_label), c("ac", "ad", "bc", "bd", "cd"))
  expect_false(any(values$tests$pair_label == "ab"))
  expect_equal(sort(values$tests$item1_name), c("I1", "I1", "I2", "I2", "I3"))
  expect_equal(sort(values$tests$item2_name), c("I3", "I3", "I4", "I4", "I4"))
  expect_true(all(c("wpg_gamma", "wpg_gamma_source") %in% names(values$tests)))
  expect_equal(unique(values$tests$wpg_gamma_source), "item_screening")
  expect_true(all(c(
    "stop_reason", "attempted_n_step", "attempted_delta",
    "attempted_converged", "attempted_stop_reason",
    "reported_checkpoint_step", "report_value_source"
  ) %in% names(values$tests)))
  expect_true(all(values$tests$report_value_source %in% c("attempted_fit", "source_first_post_50_checkpoint")))

  gamma_context <- build_local_independence_gamma_context(fit$project)
  for (row in seq_len(nrow(values$tests))) {
    item1 <- match(values$tests$item1_name[[row]], fit$project$items$name)
    item2 <- match(values$tests$item2_name[[row]], fit$project$items$name)
    expected <- local_independence_wpg_gamma(fit$project, gamma_context, item1, item2)
    expect_equal(values$tests$wpg_gamma[[row]], expected)
  }

  checkpoint_values <- local_independence_values(fit, max_step = 80L, max_delta = 0, jobs = 1L)
  checkpoint_rows <- checkpoint_values$tests$report_value_source == "source_first_post_50_checkpoint"
  expect_true(any(checkpoint_rows))
  expect_true(all(!is.na(checkpoint_values$tests$attempted_n_step[checkpoint_rows])))
  expect_true(all(!is.na(checkpoint_values$tests$attempted_stop_reason[checkpoint_rows])))
  expect_true(all(checkpoint_values$tests$reported_checkpoint_step[checkpoint_rows] == 51L))
})

test_that("GLLRM DIF checks use the current model", {
  ia <- gRm(gllrm_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)

  values <- dif_tests_values(fit, max_step = 80L, max_delta = 1e-5, jobs = 2L)

  expect_s3_class(values, "gRm_dif_tests_values")
  expect_equal(sort(paste0(values$tests$item_label, values$tests$background_label)), c("ae", "be", "ce"))
  expect_equal(paste0(values$included_tests$item_label, values$included_tests$background_label), "de")
  expect_true(is.data.frame(values$included_tests))
  expect_gt(nrow(values$included_tests), 0L)
  expect_false(any(values$tests$item_name == "I4" & values$tests$background_name == "X1"))
  expect_true(all(c("gamma", "p_gamma", "gamma_source") %in% names(values$tests)))
  expect_equal(unique(values$tests$gamma_source), "item_screening")
  expect_true(any(is.finite(values$tests$gamma)))
  expect_false(anyNA(values$tests$p_gamma))
  expect_true(all(c(
    "stop_reason", "attempted_n_step", "attempted_delta",
    "attempted_converged", "attempted_stop_reason",
    "reported_checkpoint_step", "report_value_source"
  ) %in% names(values$tests)))
  expect_true(all(is.na(values$tests$reported_checkpoint_step)))
  expect_equal(unique(values$tests$report_value_source), "attempted_fit")

  gamma_context <- build_dif_gamma_context(fit$project)
  for (row in seq_len(nrow(values$tests))) {
    item <- match(values$tests$item_name[[row]], fit$project$items$name)
    background <- match(values$tests$background_name[[row]], fit$project$backgrounds$name)
    expected <- dif_tests_partial_gamma_stats(fit$project, gamma_context, item, background)
    expect_equal(values$tests$gamma[[row]], expected$gamma)
    expect_equal(values$tests$p_gamma[[row]], expected$p_value)
  }
})

test_that("GLLRM global homogeneity refits the current GLLRM in score groups", {
  ia <- gRm(gllrm_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)

  values <- global_homogeneity_values(fit, score_cuts = c(2L, 4L), max_step = 80L, max_delta = 1e-5)

  expect_s3_class(values, "gRm_global_homogeneity_values")
  expect_equal(values$summary$n_parameters, fit$values$n_parameters)
  expect_equal(values$summary$full_log_likelihood, fit$values$log_likelihood, tolerance = 1e-8)
  expect_equal(nrow(values$score_groups), 2L)
  expect_true(all(vapply(values$group_values, function(x) inherits(x$fit$values, "gRm_gllrm_values"), logical(1L))))
  expect_true(all(is.na(values$items$residual)))
  expect_false(any(values$items$residual_runtime_source_backed))
  expect_false(any(values$items$marker_runtime_source_backed))
  expect_true(is.list(values$summary))
  expect_true(is.data.frame(values$score_groups))
  expect_true(is.data.frame(values$items))
})

test_that("GLLRM uniform DIF score-group tables follow source minscore convention", {
  ia <- gRm(gllrm_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)
  context <- fit$fit$context
  groups <- global_homogeneity_score_groups(context$bundle, c(2L, 4L))
  cache <- new_gllrm_probability_cache(context, fit$fit)

  tables <- gllrm_uniform_dif_scoregroup_tables(context, groups, context$dif_specs[[1L]], cache)
  complete_rows <- which(rowSums(context$item_matrix < 0L) == 0L & context$background_matrix[, 1L] >= 1L)
  expected_first_group_n <- sum(context$score[complete_rows] >= 0L & context$score[complete_rows] <= 2L)
  valid_first_group_n <- sum(context$score[context$valid_rows] >= 1L & context$score[context$valid_rows] <= 2L)

  expect_gt(expected_first_group_n, valid_first_group_n)
  expect_equal(sum(tables$observed[, , 1L]), expected_first_group_n)
})

test_that("batched uniform LD scoregroup tables match scalar source-faithful references", {
  ia <- gRm(gllrm_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2 + I2:I3, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)
  context <- fit$fit$context
  state <- fit$fit
  groups <- global_homogeneity_score_groups(context$bundle, c(2L, 4L))
  components <- context$ld_components_items

  batched <- gllrm_uniform_ld_scoregroup_tables_all(
    context,
    groups,
    new_gllrm_probability_cache(context, state, components = components)
  )

  for (ld_index in seq_along(context$ld_specs)) {
    expected <- scalar_gllrm_uniform_ld_scoregroup_tables_reference(
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

test_that("batched uniform DIF scoregroup tables match scalar source-faithful references", {
  ia <- gRm(gllrm_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2 + I2:I3, dif = ~ I3:X1 + I4:X1), max_step = 200L, max_delta = 1e-6)
  context <- fit$fit$context
  state <- fit$fit
  groups <- global_homogeneity_score_groups(context$bundle, c(2L, 4L))
  components <- context$ld_components_items

  batched <- gllrm_uniform_dif_scoregroup_tables_all(
    context,
    groups,
    new_gllrm_probability_cache(context, state, components = components)
  )

  expect_equal(length(batched), length(context$dif_specs))
  for (dif_index in seq_along(context$dif_specs)) {
    expected <- scalar_gllrm_uniform_dif_scoregroup_tables_reference(
      context,
      state,
      groups,
      context$dif_specs[[dif_index]],
      components
    )
    expect_equal(batched[[dif_index]], expected, tolerance = 1e-12)
  }
})

test_that("GLLRM item fit values use GLLRM conditional probabilities", {
  ia <- gRm(gllrm_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)

  values <- item_fits_values(fit, include_extended = TRUE)

  expect_s3_class(values, "gRm_item_fits_values")
  expect_equal(values$fit$model, "gllrm")
  expect_false(values$incomplete_records_used)
  expect_equal(values$incomplete_records_status, "not_source_backed")
  expect_true(all(c("outfit", "infit", "observed_gamma", "expected_gamma") %in% names(values$items)))
  expect_false(is.null(values$extended))
  expect_true(is.data.frame(values$items))
  expect_true(is.data.frame(values$extended$scores))
  expect_true(is.data.frame(values$extended$summaries))
})

test_that("GLLRM probability helpers match scalar source-faithful references", {
  ia <- gRm(gllrm_lid_data(), items = c("I1", "I2", "I3", "I4"), exogenous = "X1", id = "ID")
  fit <- fit(gllrm(ia, ld = ~ I1:I2 + I2:I3, dif = ~ I4:X1), max_step = 200L, max_delta = 1e-6)
  context <- fit$fit$context
  state <- fit$fit
  components <- context$ld_components_items

  checked <- 0L
  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    total_score <- group$score[[1L]]
    background_values <- gllrm_group_background_values(context, group)

    item_expected <- scalar_gllrm_item_probabilities_reference(context, state, total_score, background_values, components)
    item_actual <- gllrm_group_item_probabilities(context, state, total_score, background_values, components)
    expect_equal(item_actual, item_expected, tolerance = 1e-12)

    ld_expected <- scalar_gllrm_ld_probabilities_reference(context, state, total_score, background_values, components)
    ld_actual <- gllrm_group_ld_probabilities(context, state, total_score, background_values, components)
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

test_that("native GLLRM expected margins validate matrices before duplicating outputs", {
  path <- package_path("src", "gllrm_expected.cpp")
  skip_if(
    !file.exists(path),
    "source-tree C++ file is not installed"
  )
  cpp <- readLines(path, warn = FALSE)

  validate_line <- grep('require_numeric_matrix\\(item_gamma, "item_gamma"\\);', cpp)
  duplicate_line <- grep("Rf_duplicate\\(item_gamma\\)", cpp)

  expect_length(validate_line, 1L)
  expect_true(validate_line < duplicate_line[[1L]])
})
