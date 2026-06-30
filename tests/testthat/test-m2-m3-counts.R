m2_m3_counts_fixture_context <- function() {
  item_matrix <- rbind(
    c(0L, 0L, 0L),
    c(0L, 1L, 1L),
    c(1L, 0L, 1L),
    c(1L, 1L, 0L),
    c(1L, 1L, 1L),
    c(0L, 1L, 0L),
    c(1L, 0L, 0L),
    c(0L, 0L, 1L)
  )
  background_matrix <- rbind(
    c(1L, 1L),
    c(1L, 2L),
    c(2L, 2L),
    c(2L, 3L),
    c(1L, 3L),
    c(2L, 1L),
    c(1L, 2L),
    c(2L, 3L)
  )
  score <- rowSums(item_matrix)

  list(
    project = list(
      items = data.frame(
        name = c("I1", "I2", "I3"),
        label_code = c("A", "B", "C"),
        raw_max = c(2L, 2L, 2L),
        stringsAsFactors = FALSE
      ),
      backgrounds = data.frame(
        name = c("X1", "X2"),
        label_code = c("X", "Y"),
        raw_max = c(2L, 3L),
        stringsAsFactors = FALSE
      )
    ),
    bundle = list(
      model = list(
        least_score = 0L,
        largest_score = 3L,
        max_total_score = 3L
      )
    ),
    analysis = list(score_groups = c(1L, 3L)),
    n_items = 3L,
    n_backgrounds = 2L,
    item_raw_max = c(2L, 2L, 2L),
    background_raw_max = c(2L, 3L),
    item_matrix = item_matrix,
    background_matrix = background_matrix,
    score = score,
    valid_rows = 1:6
  )
}

expect_count_array <- function(observed, dimensions, cells) {
  expected <- array(0L, dim = dimensions)
  for (cell in cells) {
    expected[matrix(as.integer(cell), nrow = 1L)] <- expected[matrix(as.integer(cell), nrow = 1L)] + 1L
  }
  expect_equal(unname(observed), expected)
}

test_that("M2 observed count helpers tabulate all two-way families from valid rows", {
  context <- m2_m3_counts_fixture_context()
  score_group_lookup <- m2_m3_score_group_lookup(context, context$analysis$score_groups)

  expect_equal(
    unname(m2_m3_count_item_item(context, 1L, 2L)),
    matrix(c(1L, 1L, 2L, 2L), nrow = 2L)
  )
  expect_equal(
    unname(m2_m3_count_item_exogenous(context, 1L, 1L)),
    matrix(c(2L, 1L, 1L, 2L), nrow = 2L)
  )
  expect_equal(
    unname(m2_m3_count_item_score_group(context, 1L, score_group_lookup)),
    matrix(c(2L, 0L, 1L, 3L), nrow = 2L)
  )

  expect_equal(
    unname(m2_m3_count_observed(
      context,
      m2_m3_spec("item_score_group", item = 1L, score_group = TRUE),
      score_group_lookup
    )),
    matrix(c(2L, 0L, 1L, 3L), nrow = 2L)
  )
})

test_that("M2/M3 score group lookup honors explicit diagnostic score cuts", {
  context <- m2_m3_counts_fixture_context()
  score_group_lookup <- m2_m3_score_group_lookup(context, c(1L, 2L))

  expect_equal(
    attr(score_group_lookup, "score_groups")$to_score,
    c(1L, 2L)
  )
  expect_equal(
    unname(m2_m3_count_item_score_group(context, 1L, score_group_lookup)),
    matrix(c(2L, 0L, 1L, 2L), nrow = 2L)
  )
})

test_that("M3 observed count helpers tabulate all three-way families from valid rows", {
  context <- m2_m3_counts_fixture_context()
  score_group_lookup <- m2_m3_score_group_lookup(context, context$analysis$score_groups)

  expect_count_array(
    m2_m3_count_item_item_item(context, 1L, 2L, 3L),
    c(2L, 2L, 2L),
    list(
      c(1L, 1L, 1L),
      c(1L, 2L, 2L),
      c(2L, 1L, 2L),
      c(2L, 2L, 1L),
      c(2L, 2L, 2L),
      c(1L, 2L, 1L)
    )
  )
  expect_count_array(
    m2_m3_count_item_item_exogenous(context, 1L, 2L, 1L),
    c(2L, 2L, 2L),
    list(
      c(1L, 1L, 1L),
      c(1L, 2L, 1L),
      c(2L, 1L, 2L),
      c(2L, 2L, 2L),
      c(2L, 2L, 1L),
      c(1L, 2L, 2L)
    )
  )
  expect_count_array(
    m2_m3_count_item_item_score_group(context, 1L, 2L, score_group_lookup),
    c(2L, 2L, 2L),
    list(
      c(1L, 1L, 1L),
      c(1L, 2L, 2L),
      c(2L, 1L, 2L),
      c(2L, 2L, 2L),
      c(2L, 2L, 2L),
      c(1L, 2L, 1L)
    )
  )
  expect_count_array(
    m2_m3_count_item_exogenous_exogenous(context, 1L, 1L, 2L),
    c(2L, 2L, 3L),
    list(
      c(1L, 1L, 1L),
      c(1L, 1L, 2L),
      c(2L, 2L, 2L),
      c(2L, 2L, 3L),
      c(2L, 1L, 3L),
      c(1L, 2L, 1L)
    )
  )
  expect_count_array(
    m2_m3_count_item_exogenous_score_group(context, 1L, 1L, score_group_lookup),
    c(2L, 2L, 2L),
    list(
      c(1L, 1L, 1L),
      c(1L, 1L, 2L),
      c(2L, 2L, 2L),
      c(2L, 2L, 2L),
      c(2L, 1L, 2L),
      c(1L, 2L, 1L)
    )
  )

  expect_equal(
    unname(m2_m3_count_observed(
      context,
      m2_m3_spec("item_exogenous_score_group", item = 1L, exogenous = 1L, score_group = TRUE),
      score_group_lookup
    )),
    unname(m2_m3_count_item_exogenous_score_group(context, 1L, 1L, score_group_lookup))
  )
})

test_that("observed count helpers ignore rows outside context valid_rows", {
  context <- m2_m3_counts_fixture_context()
  score_group_lookup <- m2_m3_score_group_lookup(context, context$analysis$score_groups)
  valid_context <- context
  valid_context$item_matrix <- context$item_matrix[context$valid_rows, , drop = FALSE]
  valid_context$background_matrix <- context$background_matrix[context$valid_rows, , drop = FALSE]
  valid_context$score <- context$score[context$valid_rows]
  valid_context$valid_rows <- seq_along(valid_context$score)

  expect_equal(
    m2_m3_count_item_item(context, 1L, 2L),
    m2_m3_count_item_item(valid_context, 1L, 2L)
  )
  expect_equal(
    m2_m3_count_item_exogenous_score_group(context, 1L, 1L, score_group_lookup),
    m2_m3_count_item_exogenous_score_group(valid_context, 1L, 1L, score_group_lookup)
  )
})
