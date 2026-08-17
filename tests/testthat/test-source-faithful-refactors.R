test_that("GLLRM update logic has a single shared core", {
  gllrm_fit_text <- readLines(repo_path("gRm", "R", "gllrm_fit.R"), warn = FALSE)

  expect_true(any(grepl("^update_gllrm_parameters_core <- function", gllrm_fit_text)))
  expect_true(any(grepl("^update_gllrm_parameters <- function", gllrm_fit_text)))
  expect_true(any(grepl("^update_gllrm_parameters_once <- function", gllrm_fit_text)))
})

test_that("candidate dispatch helper preserves row order", {
  expect_true(exists("source_candidate_map", mode = "function"))
  if (!exists("source_candidate_map", mode = "function")) {
    return(invisible())
  }

  rows <- source_candidate_map(
    n_candidates = 4L,
    jobs = 1L,
    fit_one = function(index) {
      data.frame(index = index, value = index * 10L)
    }
  )

  expect_equal(rows$index, 1:4)
  expect_equal(rows$value, c(10L, 20L, 30L, 40L))
  expect_null(source_candidate_map(0L, jobs = 128L, fit_one = function(index) data.frame(index = index)))
})

test_that("item fit result assembly is shared", {
  item_fit_text <- readLines(repo_path("gRm", "R", "item_fits_values.R"), warn = FALSE)

  expect_true(any(grepl("^assemble_item_fits_values <- function", item_fit_text)))
  expect_true(exists("assemble_item_fits_values", mode = "function"))
})

test_that("GLLRM global homogeneity uniform summaries share an accumulator", {
  gh_text <- readLines(repo_path("gRm", "R", "global_homogeneity_ld.R"), warn = FALSE)

  expect_true(any(grepl("^gllrm_uniform_summary_stats <- function", gh_text)))
  expect_true(exists("gllrm_uniform_summary_stats", mode = "function"))
})

test_that("gamma cell-table helper lives with source gamma primitives", {
  gamma_text <- readLines(repo_path("gRm", "R", "source_gamma_stats.R"), warn = FALSE)
  item_fit_text <- readLines(repo_path("gRm", "R", "item_fits_values.R"), warn = FALSE)

  expect_true(any(grepl("^gamma_cell_tables <- function", gamma_text)))
  expect_false(any(grepl("^gamma_cell_tables <- function", item_fit_text)))
})

test_that("summary construction has a shared internal builder", {
  summary_text <- readLines(repo_path("gRm", "R", "api-summary.R"), warn = FALSE)

  expect_true(any(grepl("^make_gRm_summary <- function", summary_text)))
  expect_true(exists("make_gRm_summary", mode = "function"))
})

test_that("summary table normalization no longer refers to retired tidy API", {
  namespace <- asNamespace("gRm")
  retired_name <- paste0("transform", "_for", "_tidy")

  expect_true(exists("normalize_summary_table", envir = namespace, inherits = FALSE))
  expect_false(exists(retired_name, envir = namespace, inherits = FALSE))
})

test_that("base global homogeneity does not expose unsupported variance helpers", {
  namespace <- asNamespace("gRm")

  expect_false(exists("global_homogeneity_expected_score_variances", envir = namespace, inherits = FALSE))
  expect_false(exists("global_homogeneity_expected_summary_cpp", envir = namespace, inherits = FALSE))
  expect_false(exists("global_homogeneity_score_item_n", envir = namespace, inherits = FALSE))
  expect_false(exists("global_homogeneity_component_score_gamma", envir = namespace, inherits = FALSE))
  expect_false(exists("global_homogeneity_full_score_gamma", envir = namespace, inherits = FALSE))
})

test_that("included DIF diagnostics do not expose orphan conditional-test helpers", {
  namespace <- asNamespace("gRm")
  helper_names <- c(
    "gllrm_no_dif_pair_test",
    "gllrm_no_dif_conditioning",
    "gllrm_item_dif_backgrounds",
    "gllrm_model_dif_backgrounds",
    "gllrm_background_dif_items"
  )

  for (helper_name in helper_names) {
    expect_false(
      exists(helper_name, envir = namespace, inherits = FALSE),
      info = paste(helper_name, "should not be available as an alternate included-DIF path")
    )
  }
})

test_that("included LD diagnostics do not expose orphan conditional-test helpers", {
  namespace <- asNamespace("gRm")
  helper_names <- c(
    "gllrm_li_pair_test",
    "gllrm_li_bias_backgrounds"
  )

  for (helper_name in helper_names) {
    expect_false(
      exists(helper_name, envir = namespace, inherits = FALSE),
      info = paste(helper_name, "should not be available as an alternate included-LD path")
    )
  }
})

test_that("screen J conditional bias calculations have source-faithful helper boundaries", {
  helper_names <- c(
    "screen_j_conditional_try_native",
    "screen_j_conditional_valid_rows",
    "screen_j_conditional_slice_stats",
    "screen_j_conditional_exact_results"
  )
  for (helper_name in helper_names) {
    expect_true(
      exists(helper_name, mode = "function"),
      info = paste(helper_name, "should be an explicit helper")
    )
  }

  x <- c(1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 2L, 3L)
  y <- c(1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 1L)
  condition_values <- data.frame(score = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L))
  condition_dims <- 4L
  valid <- rep(TRUE, length(x))

  slice_stats <- screen_j_conditional_slice_stats(
    x = x,
    y = y,
    x_dim = 3L,
    y_dim = 2L,
    condition_values = condition_values,
    condition_dims = condition_dims,
    valid = valid,
    exact = TRUE
  )
  full_stats <- screen_j_conditional_bias_test(
    x = x,
    y = y,
    x_dim = 3L,
    y_dim = 2L,
    condition_values = condition_values,
    condition_dims = condition_dims,
    valid = valid,
    exact = FALSE,
    native = FALSE
  )

  expect_true(slice_stats$has_rows)
  expect_equal(slice_stats$chi_square, full_stats$chi_square, tolerance = 0)
  expect_equal(slice_stats$df, full_stats$df)
  expect_equal(slice_stats$gamma, full_stats$gamma, tolerance = 0)
  expect_equal(slice_stats$p_chi_asymp, full_stats$p_chi_asymp, tolerance = 0)
  expect_equal(slice_stats$p_gamma_asymp, full_stats$p_gamma_asymp, tolerance = 0)
  expect_true(length(slice_stats$slices) > 0L)
  expect_true(all(vapply(slice_stats$slices, screen_j_source_informative_slice, logical(1L))))
})

test_that("SCREEN J native source-random aliases are registered centrally", {
  conditional_text <- paste(readLines(repo_path("gRm", "R", "screen_j_conditional.R"), warn = FALSE), collapse = "\n")
  native_text <- paste(readLines(repo_path("gRm", "src", "screen_j_exact.cpp"), warn = FALSE), collapse = "\n")
  registration_text <- paste(readLines(repo_path("gRm", "src", "init.c"), warn = FALSE), collapse = "\n")

  expect_match(conditional_text, "screen_j_conditional_try_native", fixed = TRUE)
  expect_match(conditional_text, "screen_j_conditional_exact_results", fixed = TRUE)
  expect_false(grepl("screen_j_source_random_draws_native", conditional_text, fixed = TRUE))
  expect_false(grepl("gRm_screen_j_source_random_draws", conditional_text, fixed = TRUE))
  expect_match(native_text, "gRm_screen_j_source_random_draws", fixed = TRUE)
  expect_match(registration_text, '"gRm_screen_j_source_random_draws"', fixed = TRUE)
  expect_match(registration_text, '"_gRm_screen_j_source_random_draws"', fixed = TRUE)
})

test_that("SCREEN J conditional native comments describe parity coverage", {
  screen_j_text <- paste(readLines(repo_path("gRm", "R", "screen_j_exact_native.R"), warn = FALSE), collapse = "\n")

  expect_false(grepl("keeps the gate closed until the C++ semantics are fixed", screen_j_text, fixed = TRUE))
  expect_match(screen_j_text, "retained as parity", fixed = TRUE)
  expect_match(screen_j_text, "coverage for the optimized native route", fixed = TRUE)
})

test_that("native registration omits unused validation/debug probes", {
  screen_j_native <- paste(readLines(repo_path("gRm", "src", "screen_j_exact.cpp"), warn = FALSE), collapse = "\n")
  item_parameter_native <- paste(readLines(repo_path("gRm", "src", "item_parameters_extended.cpp"), warn = FALSE), collapse = "\n")
  removed_symbols <- c(
    "_gRm_screen_j_exact_kernel",
    "gRm_screen_j_exact_kernel",
    "gRm_item_parameters_extended_ice_fields",
    "build_gamma_from_matrix",
    "gRm_global_homogeneity_expected_summary",
    "gRm_global_homogeneity_expected_summary_from_fit",
    "HomogeneityItemSummary",
    "global_homogeneity_item_summary_cpp"
  )

  for (symbol in removed_symbols) {
    expect_false(
      grepl(symbol, screen_j_native, fixed = TRUE),
      info = paste(symbol, "should not be registered or forward-declared")
    )
    expect_false(
      grepl(symbol, item_parameter_native, fixed = TRUE),
      info = paste(symbol, "should not remain as a compiled validation/debug probe")
    )
  }
})

test_that("top ICE diagnostic wrapper delegates source cancellation formatting", {
  item_parameter_native <- paste(readLines(repo_path("gRm", "src", "item_parameters_extended.cpp"), warn = FALSE), collapse = "\n")
  count_fixed <- function(text, pattern) {
    matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
    if (identical(matches, -1L)) 0L else length(matches)
  }

  expect_equal(count_fixed(item_parameter_native, "volatile long double z = source_ln"), 1L)
  expect_match(item_parameter_native, "return TopIceDiagnostic{value, field9_3_long_double(value)};", fixed = TRUE)
  expect_match(item_parameter_native, "TopIceDiagnostic diagnostic = top_ice_diagnostic(", fixed = TRUE)
})

test_that("score-cut API omits retired score-spec scaffolding", {
  namespace <- asNamespace("gRm")
  retired_names <- c(
    "sum_score",
    "score_groups_auto",
    "score_groups_cut",
    "validate_score_spec",
    "validate_score_group_spec",
    "resolve_gRm_score_groups"
  )

  for (name in retired_names) {
    expect_false(
      exists(name, envir = namespace, inherits = FALSE),
      info = paste(name, "should not remain as score-spec API scaffolding")
    )
  }

  constructor_text <- paste(readLines(repo_path("gRm", "R", "api-constructors.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("gRm_score_group_spec", constructor_text, fixed = TRUE))
  expect_false(grepl("gRm_score_spec", constructor_text, fixed = TRUE))
})

test_that("read_digram_project does not hide explicit roles behind a pass-through resolver", {
  namespace <- asNamespace("gRm")
  constructor_text <- paste(readLines(repo_path("gRm", "R", "api-constructors.R"), warn = FALSE), collapse = "\n")

  expect_false(exists("resolve_gRm_project_roles", envir = namespace, inherits = FALSE))
  expect_false(grepl("resolve_gRm_project_roles", constructor_text, fixed = TRUE))
  expect_true(grepl("items <- as.character(items)", constructor_text, fixed = TRUE))
  expect_true(grepl("exogenous <- as.character(exogenous %||% character())", constructor_text, fixed = TRUE))
})

test_that("project input does not expose unused variable-name validator", {
  namespace <- asNamespace("gRm")
  project_input_text <- paste(readLines(repo_path("gRm", "R", "project_input.R"), warn = FALSE), collapse = "\n")

  expect_false(exists("validate_gRm_variable_names", envir = namespace, inherits = FALSE))
  expect_false(grepl("validate_gRm_variable_names", project_input_text, fixed = TRUE))
})

test_that("summary tables do not expose dead item-analysis helper cluster", {
  namespace <- asNamespace("gRm")
  summary_table_text <- paste(readLines(repo_path("gRm", "R", "api-summary-tables.R"), warn = FALSE), collapse = "\n")
  dead_helpers <- c(
    "bind_result_rows",
    "item_analysis_score_data",
    "item_analysis_exogenous_complete",
    "score_distribution_table",
    "item_analysis_estimation_score_distribution"
  )

  for (helper in dead_helpers) {
    expect_false(exists(helper, envir = namespace, inherits = FALSE), info = helper)
    expect_false(grepl(helper, summary_table_text, fixed = TRUE), info = helper)
  }
})

test_that("source bundle emits constant invalid fields without dead local storage", {
  source_bundle_lines <- readLines(repo_path("gRm", "R", "source_bundle.R"), warn = FALSE)
  source_bundle_text <- paste(source_bundle_lines, collapse = "\n")

  expect_false(any(grepl("^\\s*invalid_items\\s*<-\\s*integer\\(nrow\\(raw\\)\\)", source_bundle_lines)))
  expect_false(any(grepl("^\\s*invalid_backgrounds\\s*<-\\s*integer\\(nrow\\(raw\\)\\)", source_bundle_lines)))
  expect_true(grepl("invalid_items = integer(nrow(raw))", source_bundle_text, fixed = TRUE))
  expect_true(grepl("invalid_backgrounds = integer(nrow(raw))", source_bundle_text, fixed = TRUE))
  expect_true(grepl("ninvalid_items = 0L", source_bundle_text, fixed = TRUE))
  expect_true(grepl("ninvalid_backgrounds = 0L", source_bundle_text, fixed = TRUE))
})

test_that("base DIF tests use candidate dispatch without dead preallocation", {
  dif_text <- paste(readLines(repo_path("gRm", "R", "dif_tests_values.R"), warn = FALSE), collapse = "\n")

  expect_false(grepl("tests <- data.frame(", dif_text, fixed = TRUE))
  expect_true(grepl("tests <- source_candidate_map(nrow(candidates), jobs, fit_one)", dif_text, fixed = TRUE))
})

test_that("native GLLRM expected native boundary omits dead metadata fields", {
  gllrm_fit_text <- paste(readLines(repo_path("gRm", "R", "gllrm_fit.R"), warn = FALSE), collapse = "\n")
  native_text <- paste(readLines(repo_path("gRm", "src", "gllrm_expected.cpp"), warn = FALSE), collapse = "\n")

  expect_true(grepl("keys <- vapply(components, gllrm_component_key, character(1L))", gllrm_fit_text, fixed = TRUE))
  expect_false(grepl("component_keys =", gllrm_fit_text, fixed = TRUE))
  expect_false(grepl("int item_gamma_rows = matrix_nrows(item_gamma);", native_text, fixed = TRUE))
  expect_true(grepl("matrix_nrows(item_gamma) != n_items", native_text, fixed = TRUE))
})

test_that("rbind_fill preserves existing column types", {
  first <- data.frame(
    id = 1L,
    score = 1.5,
    flag = TRUE,
    label = "a",
    stringsAsFactors = FALSE
  )
  second <- data.frame(
    id = 2L,
    other_score = 2.5,
    other_flag = FALSE,
    other_label = "b",
    stringsAsFactors = FALSE
  )

  out <- rbind_fill(first, second)

  expect_type(out$id, "integer")
  expect_type(out$score, "double")
  expect_type(out$flag, "logical")
  expect_type(out$label, "character")
  expect_type(out$other_score, "double")
  expect_type(out$other_flag, "logical")
  expect_type(out$other_label, "character")
  expect_true(is.na(out$score[[2L]]))
  expect_true(is.na(out$flag[[2L]]))
  expect_true(is.na(out$label[[2L]]))
})
