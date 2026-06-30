test_that("global invariance implementation is outside the current package scope", {
  r_dir <- repo_path("gRm", "R")
  man_dir <- repo_path("gRm", "man")

  expect_false(file.exists(file.path(r_dir, "global_invariance_values.R")))
  expect_false(file.exists(file.path(man_dir, "global_invariance_values.Rd")))
  expect_false(file.exists(file.path(man_dir, "subset_bundle_to_background_value.Rd")))
  expect_false(file.exists(file.path(man_dir, "source_exo_value_name.Rd")))

  package_files <- c(
    list.files(r_dir, pattern = "[.]R$", full.names = TRUE),
    list.files(man_dir, pattern = "[.]Rd$", full.names = TRUE)
  )
  package_text <- unlist(lapply(package_files, readLines, warn = FALSE), use.names = FALSE)

  expect_false(any(grepl("global_invariance|global invariance|Global Invariance", package_text)))
  expect_false(any(grepl("global[_ -]?DIF|global dif", package_text, ignore.case = TRUE)))
})

test_that("local invariance is outside the current package scope", {
  expect_false("local_invariance" %in% getNamespaceExports("gRm"))
  expect_false("local_invariance_values" %in% ls(asNamespace("gRm"), all.names = TRUE))
})
