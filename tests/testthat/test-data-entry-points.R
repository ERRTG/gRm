test_that("read_digram_csv builds internal representation and writes import files", {
  csv_path <- tempfile(fileext = ".csv")
  data <- data.frame(
    ID = 1:3,
    item_a = c(1L, 2L, 3L),
    item_b = c(2L, 1L, 3L),
    site = c(1L, 2L, 1L)
  )
  utils::write.csv(data, csv_path, row.names = FALSE, quote = FALSE)
  output_dir <- tempfile("digr")

  project <- read_digram_csv(
    csv_path = csv_path,
    items = c("item_a", "item_b"),
    exo = "site",
    idvar = "ID",
    save_digram_files = TRUE,
    output_dir = output_dir
  )

  expect_s3_class(project, "gRm_data")
  expect_equal(project$items$name, c("item_a", "item_b"))
  expect_equal(project$backgrounds$name, "site")
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
    site = c(1L, 2L, 1L)
  )
  utils::write.csv(data, csv_path, row.names = FALSE, quote = FALSE)
  output_dir <- tempfile("digr")

  original <- read_digram_csv(
    csv_path = csv_path,
    items = c("item_a", "item_b"),
    exo = "site",
    idvar = "ID",
    save_digram_files = TRUE,
    output_dir = output_dir
  )
  restored <- read_digram_files(
    input_dir = output_dir,
    items = c("item_a", "item_b"),
    exo = "site",
    idvar = NULL
  )

  expect_s3_class(restored, "gRm_data")
  expect_equal(restored$variables, original$variables)
  expect_equal(restored$raw_data, original$raw_data)
  expect_equal(restored$import$loader, "read_digram_files")
})

test_that("saved DIGRAM import files use encoded raw_data categories", {
  input_dir <- tempfile("digram_csv_source")
  dir.create(input_dir)
  csv_path <- file.path(input_dir, "input.csv")
  data <- data.frame(
    id = seq_len(5L),
    I1 = c(0L, 1L, 0L, 1L, NA),
    I2 = c(1L, 1L, 0L, 0L, 1L),
    X = c("clinic_a", "clinic_b", "clinic_a", NA, "clinic_b"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(data, csv_path, row.names = FALSE, na = "")
  output_dir <- file.path(input_dir, "digram")

  project <- read_digram_csv(
    csv_path = csv_path,
    items = c("I1", "I2"),
    exo = "X",
    idvar = "id",
    output_dir = output_dir,
    save_digram_files = TRUE,
    na.strings = ""
  )

  saved <- utils::read.csv(
    file.path(output_dir, "DIGRAM.csv"),
    na.strings = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_equal(names(saved), c("id", "I1", "I2", "X"))
  expect_equal(saved$I1[1:4], c(1L, 2L, 1L, 2L))
  expect_equal(saved$I2, c(2L, 2L, 1L, 1L, 2L))
  expect_true(all(na.omit(saved$X) %in% c(1L, 2L)))
  expect_false(any(na.omit(saved$X) %in% c("clinic_a", "clinic_b")))
  expect_true(is.na(saved$I1[[5L]]))
  expect_true(is.na(saved$X[[4L]]))

  restored <- read_digram_files(
    input_dir = output_dir,
    items = c("I1", "I2"),
    exo = "X",
    idvar = "id",
    name = "DIGRAM",
    na.strings = ""
  )
  expect_equal(restored$raw_data, project$raw_data)
})

test_that("read_digram_files follows DIGRAM.imp project name and imv structure", {
  input_dir <- tempfile("digr_imp")
  project_dir <- file.path(input_dir, "legacy_source")
  dir.create(project_dir, recursive = TRUE)

  data <- data.frame(
    ID = 1:3,
    item_a = c(1L, 2L, 1L),
    item_b = c(2L, 1L, 2L),
    site = c(1L, 2L, 1L)
  )
  utils::write.csv(
    data,
    file.path(project_dir, "legacy_project.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  writeLines(
    c(project_dir, "legacy_project", "-", "legacy_project.cmd"),
    file.path(input_dir, "DIGRAM.imp"),
    useBytes = TRUE
  )
  writeLines(
    c(
      "a,item_a,1,one,2,two",
      "< recursive level marker",
      "b,item_b,1,one,2,two",
      "c,site,1,one,2,two"
    ),
    file.path(project_dir, "legacy_project.imv"),
    useBytes = TRUE
  )

  restored <- read_digram_files(
    input_dir = input_dir,
    items = c("item_a", "item_b"),
    exo = "site",
    idvar = NULL
  )

  expect_equal(basename(restored$paths$csv), "legacy_project.csv")
  expect_equal(basename(restored$paths$imv), "legacy_project.imv")
  expect_equal(restored$import$project_name, "legacy_project")
  expect_false("command_file" %in% names(restored$import))
  expect_equal(restored$import$imv_level_markers, "< recursive level marker")
  expect_equal(restored$variables$name, c("item_a", "item_b", "site"))
  expect_equal(restored$raw_data, unname(as.matrix(data[c("item_a", "item_b", "site")])))

  analysis <- read_digram_project(
    input_dir,
    items = c("item_a", "item_b"),
    exogenous = "site",
    id = NULL,
    score_cuts = c(1L, 2L)
  )
  expect_equal(analysis$name, "legacy_project")
  expect_equal(basename(analysis$project$paths$csv), "legacy_project.csv")
})

test_that("read_digram_imv enforces contiguous one-based category declarations", {
  expect_invalid_imv <- function(line) {
    imv_path <- tempfile(fileext = ".imv")
    writeLines(line, imv_path, useBytes = TRUE)
    expect_error(read_digram_imv(imv_path), "contiguous one-based")
  }

  expect_invalid_imv("a,item_a,1,one,3,three")
  expect_invalid_imv("a,item_a,1,one,1,uno")
  expect_invalid_imv("a,item_a,0,zero,1,one")
  expect_invalid_imv("a,item_a,-1,minus,1,one")
  expect_invalid_imv("a,item_a,1.5,one-point-five,2,two")
  expect_invalid_imv("a,item_a,x,unknown,1,one")

  imv_path <- tempfile(fileext = ".imv")
  writeLines("a,item_a,1,one,2,two,3,three", imv_path, useBytes = TRUE)
  parsed <- read_digram_imv(imv_path)
  expect_equal(parsed$raw_max, 3L)
})

test_that("DIGRAM imv writer labels supported categories as one-based", {
  labels <- digram_category_labels()
  expect_false("0" %in% names(labels))
  expect_false("zero" %in% unname(labels))
  expect_equal(labels[[as.character(1L)]], "one")
  expect_equal(labels[[as.character(2L)]], "two")

  io_text <- paste(readLines(repo_path("gRm", "R", "digram_project_io.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("seq_along(labels) - 1L", io_text, fixed = TRUE))
  expect_true(grepl("seq_along(labels)", io_text, fixed = TRUE))
})

test_that("gRm maps zero-based observed levels to DIGRAM raw categories", {
  data <- data.frame(
    ID = 1:8,
    item_a = c(0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L),
    item_b = c(1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L),
    site = c(0L, 0L, 1L, 1L, 0L, 1L, 0L, 1L)
  )

  analysis <- gRm(
    data = data,
    items = c("item_a", "item_b"),
    exogenous = "site",
    id = "ID",
    score_cuts = c(1L, 2L)
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
    site = c("clinic_a", "clinic_b", "clinic_a", "clinic_b")
  )

  analysis <- gRm(
    data = data,
    items = c("item_a", "item_b"),
    exogenous = "site",
    id = "ID",
    item_levels = list(item_a = 0:1, item_b = 0:2),
    exogenous_levels = list(site = c("clinic_a", "clinic_b")),
    score_cuts = c(1L, 3L)
  )

  expect_equal(analysis$project$items$raw_max, c(2L, 3L))
  expect_equal(analysis$project$backgrounds$raw_max, 2L)
  expect_equal(analysis$project$raw_data[, 2L], data$item_b + 1L)
  expect_equal(analysis$project$raw_data[, 3L], c(1L, 2L, 1L, 2L))

  expect_error(
    gRm(
      data = data,
      items = c("item_a", "item_b"),
      exogenous = "site",
      id = "ID",
      item_levels = 0:1,
      exogenous_levels = list(site = c("clinic_a", "clinic_b"))
    ),
    "item_levels for item_b do not cover observed values"
  )
})

test_that("new data entry points do not depend on old constructors", {
  source_lines <- readLines(repo_path("gRm", "R", "project_input.R"), warn = FALSE)
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

test_that("source bundle construction avoids growing objects in row loops", {
  source_lines <- readLines(repo_path("gRm", "R", "source_bundle.R"), warn = FALSE)
  source_text <- paste(source_lines, collapse = "\n")

  expect_false(grepl("complete_item_scores\\s*<-\\s*c\\s*\\(", source_text))
  expect_false(grepl("model_rows\\s*<-\\s*rbind\\s*\\(", source_text))
})

test_that("write_source_bundle preserves source model row order", {
  data <- data.frame(
    ID = 1:4,
    item_a = c(0L, 1L, 0L, 1L),
    item_b = c(1L, 0L, 1L, 0L),
    clinic = c(0L, 0L, 1L, 1L)
  )
  analysis <- gRm(
    data = data,
    items = c("item_a", "item_b"),
    exogenous = "clinic",
    id = "ID",
    score_cuts = c(1L, 2L)
  )
  bundle <- build_item_parameters_bundle(analysis$project)
  output_dir <- tempfile("source_bundle")

  write_source_bundle(bundle, output_dir)
  model_rows <- utils::read.delim(file.path(output_dir, "model.tsv"), stringsAsFactors = FALSE)

  expect_equal(model_rows$record_type, c("scenario", "item", "item", "background"))
  expect_equal(model_rows$position, c(0L, 1L, 2L, 1L))
  expect_equal(model_rows$name, c("DIGRAM", "item_a", "item_b", "clinic"))
  expect_equal(model_rows$max_score, c(2L, 1L, 1L, 2L))
})

test_that("source bundle manifest preserves DIGRAM incomplete and useless counters", {
  data <- data.frame(
    ID = 1:5,
    item_a = c(1L, NA, 1L, NA, 0L),
    item_b = c(0L, 1L, 0L, 1L, 0L),
    clinic = c(0L, 0L, NA, NA, NA)
  )
  analysis <- gRm(
    data = data,
    items = c("item_a", "item_b"),
    exogenous = "clinic",
    id = "ID",
    item_levels = list(item_a = 0:1, item_b = 0:1),
    exogenous_levels = list(clinic = 0:1),
    score_cuts = c(1L, 2L)
  )

  bundle <- build_item_parameters_bundle(analysis$project)

  expect_equal(bundle$manifest$nmissing_items, 1L)
  expect_equal(bundle$manifest$nmissing_backgrounds, 2L)
  expect_equal(bundle$data$missing_items, c(0L, 1L, 0L, 1L, 0L))
  expect_equal(bundle$data$missing_backgrounds, c(0L, 0L, 1L, 1L, 1L))
  expect_equal(bundle$model$least_score, 1L)
  expect_equal(bundle$data$score[[5L]], -1L)
})
