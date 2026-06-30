test_that("source statistics helpers live in their neutral source file", {
  source_file <- repo_path("gRm", "R", "source_statistics.R")
  expect_true(file.exists(source_file))
  if (!file.exists(source_file)) {
    return(invisible())
  }

  source_text <- readLines(source_file, warn = FALSE)
  dif_text <- readLines(repo_path("gRm", "R", "dif_tests_values.R"), warn = FALSE)

  expect_true(any(grepl("^source_tail_norm <- function", source_text)))
  expect_true(any(grepl("^source_pfchi <- function", source_text)))
  expect_true(any(grepl("^source_bh_critical <- function", source_text)))
  expect_false(any(grepl("^source_tail_norm <- function", dif_text)))
  expect_false(any(grepl("^source_pfchi <- function", dif_text)))
  expect_false(any(grepl("^source_bh_critical <- function", dif_text)))
})

test_that("source statistics helpers retain current DIGRAM conventions", {
  expect_equal(source_tail_norm(0, TRUE), 0.5, tolerance = 0)
  expect_equal(source_tail_norm(2.5, TRUE), 0.0062096653257761323, tolerance = 1e-15)
  expect_equal(source_tail_norm(-2.5, TRUE), 0.99379033467422395, tolerance = 1e-15)

  expect_equal(source_pfchi(0L, 10), 1, tolerance = 0)
  expect_equal(source_pfchi(1L, 3.84), 0.050043521248704592, tolerance = 1e-15)
  expect_equal(source_pfchi(4L, 9.5), 0.049747247417943646, tolerance = 1e-15)
  expect_equal(source_pfchi(120L, 150), 0.033075278818784526, tolerance = 1e-15)

  expect_equal(source_bh_critical(c(0.001, 0.02, 0.2), 0.05), 0.03333333333333333, tolerance = 0)
  expect_equal(source_bh_critical(numeric(), 0.05), 0, tolerance = 0)
})

test_that("item parameter Nresponses counts non-missing incomplete item cells by default", {
  bundle <- list(
    data = data.frame(
      item1 = c(0L, 2L),
      item2 = c(1L, 0L),
      item3 = c(1L, 1L),
      item4 = c(-1L, -1L),
      background = c(1L, 1L),
      status = c(0L, 0L),
      stringsAsFactors = FALSE
    ),
    model = list(
      items = data.frame(
        name = c("item1", "item2", "item3", "item4"),
        raw_max = c(4L, 4L, 4L, 4L),
        stringsAsFactors = FALSE
      ),
      backgrounds = data.frame(
        name = "background",
        raw_max = 2L,
        stringsAsFactors = FALSE
      ),
      least_score = 1L,
      largest_score = 7L
    )
  )

  expect_equal(item_parameters_input_stats(bundle)$n_responses, 3L)
})

test_that("item parameter Nresponses counts positive scores when focal item is missing", {
  bundle <- list(
    data = data.frame(
      item1 = c(0L, 2L),
      item2 = c(1L, 0L),
      item3 = c(1L, 1L),
      item4 = c(-1L, -1L),
      background = c(1L, 1L),
      status = c(0L, 0L),
      stringsAsFactors = FALSE
    ),
    model = list(
      items = data.frame(
        name = c("item1", "item2", "item3", "item4"),
        raw_max = c(4L, 4L, 4L, 4L),
        stringsAsFactors = FALSE
      ),
      backgrounds = data.frame(
        name = "background",
        raw_max = 2L,
        stringsAsFactors = FALSE
      ),
      least_score = 1L,
      largest_score = 7L
    )
  )
  gllrm_context <- list(
    ld_specs = list(list(item1 = 1L, item2 = 4L)),
    dif_specs = list()
  )

  expect_equal(
    item_parameters_input_stats(bundle, gllrm_context = gllrm_context)$n_responses,
    2L
  )
})
