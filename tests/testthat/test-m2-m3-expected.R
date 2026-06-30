m2_m3_expected_data <- function() {
  rows <- expand.grid(
    I1 = 1:2,
    I2 = 1:2,
    I3 = 1:2,
    I4 = 1:2,
    X1 = 1:2,
    X2 = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  weights <- 1L +
    as.integer(rows$I1 == rows$I2) +
    as.integer(rows$I2 == rows$I3 & rows$X1 == 2L) +
    as.integer(rows$I4 == 2L & rows$X2 == 1L)
  rows <- rows[rep(seq_len(nrow(rows)), weights), , drop = FALSE]
  data.frame(ID = seq_len(nrow(rows)), rows, row.names = NULL)
}

m2_m3_expected_fixture <- local({
  fixtures <- list()

  function(kind = c("ld", "no_ld", "rasch")) {
    kind <- match.arg(kind)
    if (!is.null(fixtures[[kind]])) {
      return(fixtures[[kind]])
    }

    analysis <- gRm(
      m2_m3_expected_data(),
      items = c("I1", "I2", "I3", "I4"),
      exogenous = c("X1", "X2"),
      id = "ID"
    )
    model <- if (identical(kind, "ld")) {
      gllrm(analysis, ld = ~ I1:I2 + I2:I3, dif = ~ I4:X1 + I3:X2)
    } else if (identical(kind, "no_ld")) {
      gllrm(analysis, dif = ~ I4:X1 + I3:X2)
    } else {
      gllrm(analysis)
    }
    fitted <- fit(model, max_step = 80L, max_delta = 1e-5)
    context_state <- m2_m3_fit_context_state(fitted)
    fixtures[[kind]] <<- list(
      context = context_state$context,
      state = context_state$state
    )
    fixtures[[kind]]
  }
})

m2_m3_test_score_group_lookup <- function(context) {
  m2_m3_score_group_lookup(context, c(1L, 3L))
}

m2_m3_scalar_convolve <- function(gammas, max_total_score) {
  out <- c(1, numeric(max_total_score))
  for (gamma in gammas) {
    out <- convolve_score_vectors(out, gamma, max_total_score)
  }
  out
}

m2_m3_scalar_focus_probabilities <- function(context, state, focus_items, total_score, background_values) {
  focus_items <- as.integer(focus_items)
  dims <- as.integer(context$item_raw_max[focus_items])
  out <- array(0, dim = dims)
  if (total_score < 0L || total_score > context$max_total_score) {
    return(out)
  }

  components <- context$ld_components_items
  focus_component_indices <- unique(context$ld_component_of[focus_items])
  component_gammas <- lapply(components, function(component_items) {
    gamma <- gllrm_component_gamma(context, state, component_items, background_values)
    gamma$gamma * gamma$scale
  })
  denominator <- m2_m3_scalar_convolve(component_gammas, context$max_total_score)[[total_score + 1L]]
  if (denominator <= 0) {
    return(out)
  }

  rest_component_indices <- setdiff(seq_along(components), focus_component_indices)
  rest_gamma <- m2_m3_scalar_convolve(component_gammas[rest_component_indices], context$max_total_score)
  focus_blocks <- lapply(focus_component_indices, function(component_index) {
    component_items <- components[[component_index]]
    key <- gllrm_component_key(component_items)
    config_matrix <- context$component_config_matrices[[key]]
    weights <- gllrm_component_config_weights_fast(
      context,
      state,
      component_items,
      config_matrix,
      background_values,
      key = key
    )
    data.frame(
      row = seq_len(nrow(config_matrix)),
      score = context$component_config_scores[[key]],
      weight = weights
    )
  })
  focus_blocks <- lapply(focus_blocks, function(block) block[block$weight > 0, , drop = FALSE])
  if (any(vapply(focus_blocks, nrow, integer(1L)) == 0L)) {
    return(out)
  }

  row_grid <- expand.grid(lapply(focus_blocks, function(block) seq_len(nrow(block))), KEEP.OUT.ATTRS = FALSE)
  for (grid_row in seq_len(nrow(row_grid))) {
    block_score <- 0L
    block_weight <- 1
    focus_values <- integer(length(focus_items))
    for (block_index in seq_along(focus_blocks)) {
      block <- focus_blocks[[block_index]]
      block_row <- block[row_grid[grid_row, block_index], , drop = FALSE]
      block_score <- block_score + block_row$score[[1L]]
      block_weight <- block_weight * block_row$weight[[1L]]

      component_items <- components[[focus_component_indices[[block_index]]]]
      key <- gllrm_component_key(component_items)
      config <- context$component_config_matrices[[key]][block_row$row[[1L]], ]
      matched <- match(focus_items, component_items)
      present <- !is.na(matched)
      focus_values[present] <- config[matched[present]]
    }
    if (block_score > total_score) {
      next
    }
    rest_value <- rest_gamma[[total_score - block_score + 1L]]
    if (rest_value <= 0) {
      next
    }
    out[matrix(focus_values + 1L, nrow = 1L)] <-
      out[matrix(focus_values + 1L, nrow = 1L)] + block_weight * rest_value / denominator
  }
  out
}

m2_m3_reference_item_exogenous <- function(context, state, spec) {
  out <- matrix(0, nrow = context$item_raw_max[[spec$item]], ncol = context$background_raw_max[[spec$exogenous]])
  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    background_values <- gllrm_group_background_values(context, group)
    probabilities <- gllrm_group_item_probabilities(
      context,
      state,
      total_score = group$score[[1L]],
      background_values = background_values
    )[[spec$item]]
    background_value <- background_values[[spec$exogenous]]
    out[seq_along(probabilities), background_value] <-
      out[seq_along(probabilities), background_value] + group$count[[1L]] * probabilities
  }
  out
}

m2_m3_reference_item_score_group <- function(context, state, spec, score_group_lookup) {
  out <- matrix(0, nrow = context$item_raw_max[[spec$item]], ncol = m2_m3_score_group_count(score_group_lookup))
  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    score_group <- global_homogeneity_lookup_score(score_group_lookup, group$score[[1L]])
    if (is.na(score_group)) {
      next
    }
    background_values <- gllrm_group_background_values(context, group)
    probabilities <- gllrm_group_item_probabilities(
      context,
      state,
      total_score = group$score[[1L]],
      background_values = background_values
    )[[spec$item]]
    out[seq_along(probabilities), score_group] <-
      out[seq_along(probabilities), score_group] + group$count[[1L]] * probabilities
  }
  out
}

m2_m3_reference_focus_expected <- function(context, state, spec, focus_items, extra_index = NULL) {
  dims <- c(as.integer(context$item_raw_max[focus_items]), extra_index$size %||% integer())
  out <- array(0, dim = dims)
  for (group_index in seq_len(nrow(context$score_exo_groups))) {
    group <- context$score_exo_groups[group_index, , drop = FALSE]
    background_values <- gllrm_group_background_values(context, group)
    suffix <- if (is.null(extra_index)) {
      integer()
    } else {
      extra_index$value(group, background_values)
    }
    if (length(suffix) > 0L && any(is.na(suffix))) {
      next
    }
    probabilities <- m2_m3_scalar_focus_probabilities(
      context,
      state,
      focus_items = focus_items,
      total_score = group$score[[1L]],
      background_values = background_values
    )
    out <- out + group$count[[1L]] * probabilities
    if (length(suffix) > 0L) {
      # Rebuild for suffix-bearing tables so R does not recycle over all slices.
      out <- out - group$count[[1L]] * probabilities
      prefix_grid <- arrayInd(seq_along(probabilities), dim(probabilities))
      for (cell_index in seq_len(nrow(prefix_grid))) {
        out[matrix(c(prefix_grid[cell_index, ], suffix), nrow = 1L)] <-
          out[matrix(c(prefix_grid[cell_index, ], suffix), nrow = 1L)] +
          group$count[[1L]] * probabilities[prefix_grid[cell_index, , drop = FALSE]]
      }
    }
  }
  out
}

test_that("arbitrary focus probability cache matches scalar references", {
  ld_fixture <- m2_m3_expected_fixture("ld")
  no_ld_fixture <- m2_m3_expected_fixture("no_ld")

  cases <- list(
    list(fixture = ld_fixture, focus = 1L),
    list(fixture = ld_fixture, focus = c(1L, 4L)),
    list(fixture = ld_fixture, focus = c(1L, 3L)),
    list(fixture = no_ld_fixture, focus = c(1L, 3L, 4L)),
    list(fixture = ld_fixture, focus = c(1L, 3L, 4L))
  )

  for (case in cases) {
    context <- case$fixture$context
    state <- case$fixture$state
    cache <- new_m2_m3_focus_probability_cache(context, state)
    for (group_index in seq_len(nrow(context$score_exo_groups))) {
      group <- context$score_exo_groups[group_index, , drop = FALSE]
      background_values <- gllrm_group_background_values(context, group)
      expected <- m2_m3_scalar_focus_probabilities(
        context,
        state,
        focus_items = case$focus,
        total_score = group$score[[1L]],
        background_values = background_values
      )
      actual <- m2_m3_focus_probabilities(
        cache,
        focus_items = case$focus,
        total_score = group$score[[1L]],
        background_values = background_values
      )
      expect_equal(actual, expected, tolerance = 1e-10)
    }
  }
})

test_that("focus probability cache preserves caller focus item order", {
  fixture <- m2_m3_expected_fixture("ld")
  context <- fixture$context
  state <- fixture$state
  group <- context$score_exo_groups[1L, , drop = FALSE]
  background_values <- gllrm_group_background_values(context, group)
  cache <- new_m2_m3_focus_probability_cache(context, state)

  ordered <- m2_m3_focus_probabilities(cache, c(1L, 3L), group$score[[1L]], background_values)
  reversed <- m2_m3_focus_probabilities(cache, c(3L, 1L), group$score[[1L]], background_values)

  expect_equal(dim(reversed), as.integer(context$item_raw_max[c(3L, 1L)]))
  expect_equal(reversed, aperm(ordered, c(2L, 1L)), tolerance = 1e-10)
})

test_that("two-way expected tables are current-GLLRM model margins", {
  fixture <- m2_m3_expected_fixture("ld")
  context <- fixture$context
  state <- fixture$state
  cache <- new_m2_m3_focus_probability_cache(context, state)
  score_group_lookup <- m2_m3_test_score_group_lookup(context)

  exo_spec <- m2_m3_spec("item_exogenous", item = 2L, exogenous = 1L)
  score_spec <- m2_m3_spec("item_score_group", item = 2L, score_group = TRUE)
  item_item_spec <- m2_m3_spec("item_item", item1 = 3L, item2 = 1L)

  expect_equal(
    m2_m3_expected_item_exogenous(context, exo_spec, NULL, cache),
    m2_m3_reference_item_exogenous(context, state, exo_spec),
    tolerance = 1e-10
  )
  expect_equal(
    m2_m3_expected_item_score_group(context, score_spec, score_group_lookup, cache),
    m2_m3_reference_item_score_group(context, state, score_spec, score_group_lookup),
    tolerance = 1e-10
  )
  expect_equal(
    m2_m3_expected_item_item(context, item_item_spec, NULL, cache),
    m2_m3_reference_focus_expected(context, state, item_item_spec, c(3L, 1L)),
    tolerance = 1e-10
  )
})

test_that("item-exogenous expected tables keep non-DIF exogenous level totals", {
  fixture <- m2_m3_expected_fixture("rasch")
  context <- fixture$context
  state <- fixture$state
  cache <- new_m2_m3_focus_probability_cache(context, state)
  exo_spec <- m2_m3_spec("item_exogenous", item = 1L, exogenous = 1L)

  actual <- m2_m3_expected_item_exogenous(context, exo_spec, NULL, cache)
  observed_exogenous_counts <- tabulate(
    context$background_matrix[context$valid_rows, exo_spec$exogenous],
    nbins = context$background_raw_max[[exo_spec$exogenous]]
  )

  expect_equal(colSums(actual), as.numeric(observed_exogenous_counts), tolerance = 1e-8)
})

test_that("three-way expected tables have source-shaped dimensions and totals", {
  fixture <- m2_m3_expected_fixture("ld")
  context <- fixture$context
  state <- fixture$state
  cache <- new_m2_m3_focus_probability_cache(context, state)
  score_group_lookup <- m2_m3_test_score_group_lookup(context)
  total_n <- sum(context$score_exo_groups$count)
  score_group_n <- sum(context$score_exo_groups$count[
    !is.na(vapply(
      context$score_exo_groups$score,
      function(score) global_homogeneity_lookup_score(score_group_lookup, score),
      integer(1L)
    ))
  ])

  specs <- list(
    item_item_item = m2_m3_spec("item_item_item", item1 = 1L, item2 = 3L, item3 = 4L),
    item_item_exogenous = m2_m3_spec("item_item_exogenous", item1 = 1L, item2 = 4L, exogenous = 2L),
    item_item_score_group = m2_m3_spec("item_item_score_group", item1 = 3L, item2 = 1L, score_group = TRUE),
    # Validation-sensitive: source/PAS_skunits/SKbias2.pas::calculate_EIXZ_table
    # has a preserved Cprob1/i2 branch wrinkle for item-exogenous-exogenous rows.
    item_exogenous_exogenous = m2_m3_spec("item_exogenous_exogenous", item = 2L, exogenous1 = 1L, exogenous2 = 2L),
    # Validation-sensitive: the same calculate_EIXZ_table branch covers
    # item-exogenous-score-group rows until oracle coverage is added.
    item_exogenous_score_group = m2_m3_spec("item_exogenous_score_group", item = 2L, exogenous = 1L, score_group = TRUE)
  )

  actual <- list(
    item_item_item = m2_m3_expected_item_item_item(context, specs$item_item_item, NULL, cache),
    item_item_exogenous = m2_m3_expected_item_item_exogenous(context, specs$item_item_exogenous, NULL, cache),
    item_item_score_group = m2_m3_expected_item_item_score_group(context, specs$item_item_score_group, score_group_lookup, cache),
    item_exogenous_exogenous = m2_m3_expected_item_exogenous_exogenous(context, specs$item_exogenous_exogenous, NULL, cache),
    item_exogenous_score_group = m2_m3_expected_item_exogenous_score_group(context, specs$item_exogenous_score_group, score_group_lookup, cache)
  )

  expect_equal(dim(actual$item_item_item), as.integer(context$item_raw_max[c(1L, 3L, 4L)]))
  expect_equal(dim(actual$item_item_exogenous), c(as.integer(context$item_raw_max[c(1L, 4L)]), context$background_raw_max[[2L]]))
  expect_equal(dim(actual$item_item_score_group), c(as.integer(context$item_raw_max[c(3L, 1L)]), m2_m3_score_group_count(score_group_lookup)))
  expect_equal(dim(actual$item_exogenous_exogenous), c(context$item_raw_max[[2L]], context$background_raw_max[[1L]], context$background_raw_max[[2L]]))
  expect_equal(dim(actual$item_exogenous_score_group), c(context$item_raw_max[[2L]], context$background_raw_max[[1L]], m2_m3_score_group_count(score_group_lookup)))

  expected_totals <- c(
    item_item_item = total_n,
    item_item_exogenous = total_n,
    item_item_score_group = score_group_n,
    item_exogenous_exogenous = total_n,
    item_exogenous_score_group = score_group_n
  )
  expect_equal(vapply(actual, sum, numeric(1L)), expected_totals, tolerance = 1e-8)
})

test_that("expected table dispatcher covers every Task 5 margin family", {
  fixture <- m2_m3_expected_fixture("ld")
  context <- fixture$context
  state <- fixture$state
  score_group_lookup <- m2_m3_test_score_group_lookup(context)
  cache <- new_m2_m3_focus_probability_cache(context, state)
  specs <- list(
    m2_m3_spec("item_item", item1 = 1L, item2 = 4L),
    m2_m3_spec("item_exogenous", item = 1L, exogenous = 1L),
    m2_m3_spec("item_score_group", item = 1L, score_group = TRUE),
    m2_m3_spec("item_item_item", item1 = 1L, item2 = 3L, item3 = 4L),
    m2_m3_spec("item_item_exogenous", item1 = 1L, item2 = 4L, exogenous = 1L),
    m2_m3_spec("item_item_score_group", item1 = 1L, item2 = 4L, score_group = TRUE),
    m2_m3_spec("item_exogenous_exogenous", item = 1L, exogenous1 = 1L, exogenous2 = 2L),
    m2_m3_spec("item_exogenous_score_group", item = 1L, exogenous = 1L, score_group = TRUE)
  )

  for (spec in specs) {
    expect_gt(sum(m2_m3_expected_table(context, state, spec, score_group_lookup, cache)), 0)
  }
})
