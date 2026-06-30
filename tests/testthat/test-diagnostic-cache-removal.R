test_that("LD and DIF diagnostics are not short-circuited by stale cache entries", {
  data <- data.frame(
    ID = seq_len(12L),
    I1 = c(1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
    I2 = c(1L, 2L, 1L, 3L, 2L, 3L, 3L, 1L, 2L, 2L, 3L, 1L),
    I3 = c(2L, 1L, 3L, 1L, 2L, 3L, 2L, 1L, 3L, 3L, 2L, 1L),
    site = c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 2L)
  )
  analysis <- gRm(data, items = c("I1", "I2", "I3"), exogenous = "site", id = "ID")
  project <- analysis$project

  if (
    exists("local_independence_values_cache", envir = globalenv(), inherits = FALSE) &&
      exists("local_independence_cache_key", envir = globalenv(), inherits = FALSE)
  ) {
    ld_key <- local_independence_cache_key(project, 5L, 0.1, 1L)
    assign(
      ld_key,
      structure(
        list(
          tests = data.frame(pair_label = "STALE", stringsAsFactors = FALSE),
          bh_critical_p = 0,
          suggested_ld = "STALE",
          max_step = 5L,
          max_delta = 0.1
        ),
        class = "gRm_local_independence_values"
      ),
      envir = local_independence_values_cache
    )
    on.exit(rm(list = ld_key, envir = local_independence_values_cache), add = TRUE)
  }

  ld_values <- local_independence_values(project, max_step = 5L, max_delta = 0.1, jobs = 1L)
  expect_false(identical(ld_values$tests$pair_label, "STALE"))
  expect_true(nrow(ld_values$tests) > 1L)

  if (
    exists("dif_tests_values_cache", envir = globalenv(), inherits = FALSE) &&
      exists("dif_tests_cache_key", envir = globalenv(), inherits = FALSE)
  ) {
    dif_key <- dif_tests_cache_key(project, 5L, 0.1, 1L)
    assign(
      dif_key,
      structure(
        list(
          tests = data.frame(item_label = "STALE", stringsAsFactors = FALSE),
          bh_critical_p = 0,
          max_step = 5L,
          max_delta = 0.1
        ),
        class = "gRm_dif_tests_values"
      ),
      envir = dif_tests_values_cache
    )
    on.exit(rm(list = dif_key, envir = dif_tests_values_cache), add = TRUE)
  }

  dif_values <- dif_tests_values(project, max_step = 5L, max_delta = 0.1, jobs = 1L)
  expect_false(identical(dif_values$tests$item_label, "STALE"))
  expect_true(nrow(dif_values$tests) > 1L)
})
