test_that("read_digram_csv builds internal representation and writes import files", {
  csv_path <- tempfile(fileext = ".csv")
  data <- data.frame(
    ID = 1:3,
    item_a = c(1L, 2L, 3L),
    item_b = c(2L, 1L, 3L),
    group = c(1L, 2L, 1L)
  )
  utils::write.csv(data, csv_path, row.names = FALSE, quote = FALSE)
  output_dir <- tempfile("digr")

  project <- read_digram_csv(
    csv_path = csv_path,
    items = c("item_a", "item_b"),
    exo = "group",
    idvar = "ID",
    save_digram_files = TRUE,
    output_dir = output_dir
  )

  expect_s3_class(project, "gRm_data")
  expect_equal(project$items$name, c("item_a", "item_b"))
  expect_equal(project$backgrounds$name, "group")
  expect_true(file.exists(file.path(output_dir, "DIGRAM.csv")))
  expect_true(file.exists(file.path(output_dir, "DIGRAM.imp")))
  expect_true(file.exists(file.path(output_dir, "DIGRAM.imv")))
  expect_false(file.exists(file.path(output_dir, "DIGRAM.cmd")))
})

test_that("read_digram_files reconstructs a saved import bundle", {
  csv_path <- tempfile(fileext = ".csv")
  data <- data.frame(
    ID = 1:3,
    item_a = c(1L, 2L, 3L),
    item_b = c(2L, 1L, 3L),
    group = c(1L, 2L, 1L)
  )
  utils::write.csv(data, csv_path, row.names = FALSE, quote = FALSE)
  output_dir <- tempfile("digr")

  original <- read_digram_csv(
    csv_path = csv_path,
    items = c("item_a", "item_b"),
    exo = "group",
    idvar = "ID",
    save_digram_files = TRUE,
    output_dir = output_dir
  )
  restored <- read_digram_files(
    input_dir = output_dir,
    items = c("item_a", "item_b"),
    exo = "group",
    idvar = NULL
  )

  expect_s3_class(restored, "gRm_data")
  expect_equal(restored$variables, original$variables)
  expect_equal(restored$raw_data, original$raw_data)
  expect_equal(restored$import$loader, "read_digram_files")
})

test_that("gRm maps zero-based observed levels to DIGRAM raw categories", {
  data <- data.frame(
    ID = 1:8,
    item_a = c(0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L),
    item_b = c(1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L),
    group = c(0L, 0L, 1L, 1L, 0L, 1L, 0L, 1L)
  )

  analysis <- gRm(
    data = data,
    items = c("item_a", "item_b"),
    exogenous = "group",
    id = "ID"
  )
  bundle <- build_item_parameters_bundle(analysis$project)

  expect_equal(analysis$project$items$raw_max, c(2L, 2L))
  expect_equal(analysis$project$backgrounds$raw_max, 2L)
  expect_equal(analysis$project$raw_data[, 1L], data$item_a + 1L)
  expect_equal(bundle$data$item_a, data$item_a)
  expect_equal(bundle$data$item_b, data$item_b)
  expect_gt(bundle$manifest$nvalid, 0L)
  expect_no_error(fit(gllrm(analysis)))
})

test_that("gRm accepts explicit levels and rejects uncovered observed values", {
  data <- data.frame(
    ID = 1:4,
    item_a = c(0L, 1L, 0L, 1L),
    item_b = c(0L, 1L, 2L, 1L),
    group = c("control", "case", "control", "case")
  )

  analysis <- gRm(
    data = data,
    items = c("item_a", "item_b"),
    exogenous = "group",
    id = "ID",
    item_levels = list(item_a = 0:1, item_b = 0:2),
    exogenous_levels = list(group = c("control", "case"))
  )

  expect_equal(analysis$project$items$raw_max, c(2L, 3L))
  expect_equal(analysis$project$backgrounds$raw_max, 2L)
  expect_equal(analysis$project$raw_data[, 2L], data$item_b + 1L)
  expect_equal(analysis$project$raw_data[, 3L], c(1L, 2L, 1L, 2L))

  expect_error(
    gRm(
      data = data,
      items = c("item_a", "item_b"),
      exogenous = "group",
      id = "ID",
      item_levels = 0:1,
      exogenous_levels = list(group = c("control", "case"))
    ),
    "item_levels for item_b do not cover observed values"
  )
})

test_that("new data entry points do not depend on old constructors", {
  source_lines <- readLines(repo_path("gRm", "R", "data_entry_points.R"), warn = FALSE)
  read_digram_csv_body <- paste(deparse(read_digram_csv), collapse = "\n")
  read_digram_files_body <- paste(deparse(read_digram_files), collapse = "\n")
  entry_bodies <- paste(read_digram_csv_body, read_digram_files_body, collapse = "\n")

  expect_false(grepl("read_raw_csv_project\\(", entry_bodies))
  expect_false(grepl("read_digram_project\\(", entry_bodies))
  expect_false(grepl("gRm_project\\(", entry_bodies))
  expect_false(grepl("gRm_project_from_data_frame\\(", entry_bodies))
  expect_true(any(grepl("build_gRm_internal_project", source_lines, fixed = TRUE)))
})

test_that("public data entry point defaults are limited to file conventions", {
  csv_formals <- formals(read_digram_csv)
  files_formals <- formals(read_digram_files)

  expect_equal(csv_formals$name, "DIGRAM")
  expect_equal(csv_formals$digram_folder, ".")
  expect_equal(csv_formals$na.strings, "NA")
  expect_false("sep" %in% names(csv_formals))
  expect_false("item_labels" %in% names(csv_formals))
  expect_false("exo_labels" %in% names(csv_formals))
  expect_false("line_ending" %in% names(csv_formals))
  expect_equal(files_formals$name, "DIGRAM")
  expect_equal(files_formals$na.strings, "NA")
})
