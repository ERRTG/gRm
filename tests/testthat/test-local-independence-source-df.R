test_that("base local-independence df follows observed nonzero source margins", {
  data <- data.frame(
    ID = seq_len(12L),
    I1 = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 1L),
    I2 = c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 1L)
  )
  analysis <- gRm(
    data,
    items = c("I1", "I2"),
    id = "ID",
    item_levels = list(I1 = 0:2, I2 = 0:2),
    score_cuts = c(1L, 2L)
  )

  values <- local_independence_values(
    analysis$project,
    max_step = 80L,
    max_delta = 0.0001,
    jobs = 1L
  )

  expect_equal(analysis$project$items$raw_max, c(3L, 3L))
  expect_equal(values$tests$degrees_of_freedom, 1L)
  expect_equal(
    values$tests$p_value,
    source_pfchi(values$tests$degrees_of_freedom, values$tests$chi_square)
  )
})
