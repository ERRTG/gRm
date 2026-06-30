auto_score_cut_extreme_data <- function() {
  counts <- c(`0` = 100L, `1` = 1L, `2` = 2L, `3` = 1L, `4` = 100L)
  rows <- lapply(names(counts), function(score_name) {
    score <- as.integer(score_name)
    pattern <- rep.int(0L, 4L)
    if (score > 0L) {
      pattern[seq_len(score)] <- 1L
    }
    data.frame(
      I1 = rep.int(pattern[[1L]], counts[[score_name]]),
      I2 = rep.int(pattern[[2L]], counts[[score_name]]),
      I3 = rep.int(pattern[[3L]], counts[[score_name]]),
      I4 = rep.int(pattern[[4L]], counts[[score_name]])
    )
  })
  data <- do.call(rbind, rows)
  data.frame(ID = seq_len(nrow(data)), data)
}

test_that("automatic score cuts use the source non-extreme denominator", {
  score_counts <- c(100L, 1L, 2L, 1L, 100L)
  expect_equal(item_selection_default_cut(score_counts), 2L)

  data <- auto_score_cut_extreme_data()
  analysis <- gRm(
    data,
    items = c("I1", "I2", "I3", "I4"),
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1, I4 = 0:1),
    score_cuts = "auto"
  )

  expect_equal(analysis$score_groups, c(2L, 4L))
})

test_that("automatic score cuts fail when source groups cannot be defined", {
  invalid_two_item_cases <- list(
    all_missing = data.frame(
      ID = seq_len(3L),
      I1 = c(NA_integer_, NA_integer_, NA_integer_),
      I2 = c(NA_integer_, NA_integer_, NA_integer_)
    ),
    no_complete_item_rows = data.frame(
      ID = seq_len(3L),
      I1 = c(0L, NA_integer_, 1L),
      I2 = c(NA_integer_, 1L, NA_integer_)
    ),
    lower_boundary_only = data.frame(
      ID = seq_len(3L),
      I1 = c(0L, 0L, 0L),
      I2 = c(0L, 0L, 0L)
    )
  )

  for (case_data in invalid_two_item_cases) {
    expect_error(
      gRm(
        case_data,
        items = c("I1", "I2"),
        id = "ID",
        item_levels = list(I1 = 0:1, I2 = 0:1),
        score_cuts = "auto"
      ),
      "Automatic `score_cuts`"
    )
  }

  expect_error(
    gRm(
      data.frame(ID = seq_len(4L), I1 = c(0L, 1L, 0L, 1L)),
      items = "I1",
      id = "ID",
      item_levels = list(I1 = 0:1),
      score_cuts = "auto"
    ),
    "Automatic `score_cuts`"
  )
})
