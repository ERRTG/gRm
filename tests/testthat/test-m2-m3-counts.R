m2_m3_counts_fixture_context <- function() {
  item_matrix <- rbind(
    c(0L, 0L, 0L),
    c(0L, 1L, 1L),
    c(1L, 0L, 1L),
    c(1L, 1L, 0L),
    c(1L, 1L, 1L),
    c(0L, 1L, 0L),
    c(1L, 0L, 0L),
    c(0L, 0L, 1L),
    c(0L, 1L, 1L)
  )
  background_matrix <- rbind(
    c(1L, 1L),
    c(1L, 2L),
    c(2L, 2L),
    c(2L, 3L),
    c(1L, 3L),
    c(2L, 1L),
    c(1L, 2L),
    c(2L, 3L),
    c(-1L, 2L)
  )
  score <- rowSums(item_matrix)
  score[[9L]] <- -1L

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
    estimation_rows = 2:8,
    valid_rows = 2:8
  )
}

expect_count_array <- function(observed, dimensions, cells) {
  expected <- array(0L, dim = dimensions)
  for (cell in cells) {
    expected[matrix(as.integer(cell), nrow = 1L)] <- expected[matrix(as.integer(cell), nrow = 1L)] + 1L
  }
  expect_equal(unname(observed), expected)
}

test_that("M2 observed counters include complete score-zero records outside estimation", {
  context <- m2_m3_counts_fixture_context()
  score_group_lookup <- m2_m3_score_group_lookup(context, context$analysis$score_groups)

  expect_equal(
    unname(m2_m3_count_item_item(context, 1L, 2L)),
    matrix(c(2L, 2L, 2L, 2L), nrow = 2L)
  )
  expect_equal(
    unname(m2_m3_count_item_exogenous(context, 1L, 1L)),
    matrix(c(2L, 2L, 2L, 2L), nrow = 2L)
  )
  expect_equal(
    unname(m2_m3_count_item_score_group(context, 1L, score_group_lookup)),
    matrix(c(3L, 1L, 1L, 3L), nrow = 2L)
  )

  expect_equal(
    unname(m2_m3_count_observed(
      context,
      m2_m3_spec("item_score_group", item = 1L, score_group = TRUE),
      score_group_lookup
    )),
    matrix(c(3L, 1L, 1L, 3L), nrow = 2L)
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
    matrix(c(3L, 1L, 1L, 2L), nrow = 2L)
  )
})

test_that("M3 observed counters follow routine-specific source record policies", {
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
      c(1L, 2L, 1L),
      c(2L, 1L, 1L),
      c(1L, 1L, 2L),
      c(1L, 2L, 2L)
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
      c(1L, 2L, 2L),
      c(2L, 1L, 1L),
      c(1L, 1L, 2L)
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
      c(1L, 2L, 1L),
      c(2L, 1L, 1L),
      c(1L, 1L, 1L)
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
      c(1L, 2L, 1L),
      c(2L, 1L, 2L),
      c(1L, 2L, 3L)
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
      c(1L, 2L, 1L),
      c(2L, 1L, 1L),
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

test_that("CM3 row policies preserve the missing-exogenous IJK asymmetry", {
  context <- m2_m3_counts_fixture_context()

  expect_equal(source_estimation_rows(context), 2:8)
  expect_equal(source_complete_item_exogenous_rows(context), 1:8)
  expect_equal(source_cm3_observed_ijk_rows(context), 1:9)
  expect_equal(sum(m2_m3_count_item_item(context, 1L, 2L)), 8L)
  expect_equal(sum(m2_m3_count_item_item_exogenous(context, 1L, 2L, 1L)), 8L)
  expect_equal(sum(m2_m3_count_item_item_item(context, 1L, 2L, 3L)), 9L)
})
