source_package_root <- function() {
  root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = FALSE)
  if (
    file.exists(file.path(root, "DESCRIPTION")) &&
      dir.exists(file.path(root, "R")) &&
      file.exists(file.path(root, "R", "api-summary.R"))
  ) {
    return(root)
  }
  NA_character_
}

test_that("namespace exports only the statistical modeling API", {
  exports <- getNamespaceExports("gRm")
  expected <- c(
    "gRm",
    "read_digram_project",
    "gllrm",
    "model_graph",
    "fit",
    "ari",
    "screen",
    "score_effects",
    "item_fit",
    "local_dependence",
    "dif",
    "global_homogeneity",
    "m2",
    "m3"
  )
  removed <- c(
    "sum_score",
    "score_groups_auto",
    "score_groups_cut",
    "check",
    "model_terms",
    "tidy",
    "glance",
    "details",
    "detail_names",
    "status",
    "source_trace",
    "validation",
    "unmodeled",
    "gRm_warnings"
  )

  expect_setequal(exports, expected)
  expect_false(any(removed %in% exports))
})

test_that("old output accessor generics are absent from exports and namespace", {
  old_accessors <- c("tidy", "glance", "details", "detail_names")

  expect_false(any(old_accessors %in% getNamespaceExports("gRm")))

  ns <- asNamespace("gRm")
  for (name in old_accessors) {
    expect_false(exists(name, envir = ns, inherits = FALSE), info = name)
  }
})

test_that("internal null-coalescing helper is package-owned and not exported", {
  expect_true(exists("%||%", envir = asNamespace("gRm"), inherits = FALSE))
  expect_false("%||%" %in% getNamespaceExports("gRm"))
})

test_that("public help topics coexist with explicitly internal documentation", {
  root <- source_package_root()
  skip_if(is.na(root), "source-tree help files are not installed")

  man_dir <- file.path(root, "man")
  help_pages <- sub("[.]Rd$", "", basename(list.files(man_dir, pattern = "[.]Rd$")))
  expected <- c(
    "dif",
    "fit",
    "gRm",
    "gRm-likelihood",
    "gRm-package",
    "gRm-print",
    "gRm-summary",
    "ari",
    "gllrm",
    "global_homogeneity",
    "item_fit",
    "local_dependence",
    "m2",
    "m3",
    "model_graph",
    "plot.gRm_ari",
    "read_digram_project",
    "score_effects",
    "screen",
    "update.gRm_model"
  )

  expect_true(all(expected %in% help_pages))
  internal_pages <- setdiff(help_pages, expected)
  for (page in internal_pages) {
    rd <- readLines(file.path(man_dir, paste0(page, ".Rd")), warn = FALSE)
    expect_true(
      any(grepl("\\\\keyword\\{internal\\}", rd)),
      info = paste(page, "must be marked as internal")
    )
  }
})

test_that("every registered S3 method has an exact Rd alias", {
  root <- source_package_root()
  skip_if(is.na(root), "source-tree help files are not installed")

  namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  registrations <- grep("^S3method\\(", namespace, value = TRUE)
  registered_methods <- sub(
    "^S3method\\(([^,]+),([^\\)]+)\\).*$",
    "\\1.\\2",
    registrations
  )
  rd <- unlist(lapply(
    list.files(file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE),
    readLines,
    warn = FALSE
  ), use.names = FALSE)
  aliases <- sub("^\\\\alias\\{([^}]+)\\}$", "\\1", grep("^\\\\alias\\{", rd, value = TRUE))

  expect_setequal(intersect(registered_methods, aliases), registered_methods)
})

test_that("contributor object and internal model-term contracts are generated", {
  root <- source_package_root()
  skip_if(is.na(root), "source-tree help files are not installed")

  for (topic in c("gRm-object-shapes.Rd", "gRm-model-terms.Rd")) {
    path <- file.path(root, "man", topic)
    expect_true(file.exists(path), info = topic)
    if (file.exists(path)) {
      expect_true(
        any(grepl("\\\\keyword\\{internal\\}", readLines(path, warn = FALSE))),
        info = topic
      )
    }
  }
})

test_that("workflow objects and namespace use only current public class names", {
  data <- data.frame(
    ID = 1:8,
    I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
    I2 = c(1, 0, 1, 0, 1, 0, 1, 0)
  )
  analysis <- gRm(data, items = c("I1", "I2"), id = "ID", score_cuts = c(1L, 2L))
  model <- gllrm(analysis)
  fitted <- fit(model, max_step = 20L)

  expect_identical(class(analysis), c("gRm_analysis", "list"))
  expect_identical(class(model), c("gRm_model", "list"))
  expect_identical(class(fitted), c("gRm_fit", "list"))

  root <- source_package_root()
  if (!is.na(root)) {
    namespace_file <- file.path(root, "NAMESPACE")
    namespace_lines <- readLines(namespace_file, warn = FALSE)
    stale_classes <- c("gRm_item_analysis", "gRm_gllrm_spec", "gRm_gllrm_fit")

    expect_false(any(grepl(paste(stale_classes, collapse = "|"), namespace_lines)))
  }
})

test_that("package algorithm docs do not describe retired API files or helpers", {
  doc_paths <- c(
    test_path("..", "..", "ALGORITHMS.md"),
    test_path("..", "..", "..", "docs", "GRM_R_IMPLEMENTATION_FILES.md")
  )
  docs <- unlist(lapply(doc_paths[file.exists(doc_paths)], readLines, warn = FALSE))
  stale_patterns <- c(
    "api-check[.]R",
    "api-provenance[.]R",
    "api-output[.]R",
    "screen_values[.]R",
    "global_invariance_values[.]R",
    "Graphical SCREEN Report Status",
    "algorithm check[(]",
    "status[(]x[)]",
    "source_trace[(]x[)]",
    "validation[(]x[)]",
    "unmodeled[(]x[)]",
    "gRm_warnings[(]x[)]",
    "vignettes/gRm-api-example-emo"
  )

  for (pattern in stale_patterns) {
    expect_false(any(grepl(pattern, docs)), info = pattern)
  }
})

test_that("R source directory contains only R source files", {
  root <- source_package_root()
  skip_if(is.na(root), "source-tree R directory is not installed")

  r_dir <- file.path(root, "R")
  markdown_in_r <- list.files(r_dir, pattern = "[.](md|Rmd)$", ignore.case = TRUE)

  expect_equal(markdown_in_r, character(), info = paste(markdown_in_r, collapse = ", "))
})
