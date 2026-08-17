summary_print_data <- function() {
  rows <- expand.grid(
    I1 = 1:3,
    I2 = 1:3,
    I3 = 1:3,
    site = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  data.frame(ID = seq_len(nrow(rows)), rows)
}

summary_print_analysis <- function() {
  source_data <- summary_print_data()
  gRm(
    source_data,
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID",
    score_cuts = c(2L, 6L)
  )
}

expect_summary_surface <- function(object, expected_class, which = NULL, title_pattern = "gRm:") {
  out <- if (is.null(which)) summary(object) else summary(object, which = which)

  expect_s3_class(out, expected_class)
  expect_true(any(c("tables", "tests", "parameters") %in% names(out)))
  expect_output(print(out), title_pattern)
  invisible(out)
}

public_print_summary_method_assignments <- function() {
  r_files <- list.files(repo_path("gRm", "R"), pattern = "[.]R$", full.names = TRUE)
  rows <- list()
  for (file in r_files) {
    expressions <- parse(file)
    for (expression in as.list(expressions)) {
      if (!is.call(expression) || !identical(expression[[1L]], as.name("<-"))) {
        next
      }
      lhs <- expression[[2L]]
      if (!is.symbol(lhs)) {
        next
      }
      method <- as.character(lhs)
      if (!grepl("^(print|summary)[.]gRm_", method)) {
        next
      }
      rows[[length(rows) + 1L]] <- data.frame(
        method = method,
        file = basename(file),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

test_that("public print and summary S3 methods have one implementation", {
  assignments <- public_print_summary_method_assignments()
  duplicated_methods <- unique(assignments$method[duplicated(assignments$method)])

  expect_equal(
    duplicated_methods,
    character(),
    info = paste(
      utils::capture.output(assignments[assignments$method %in% duplicated_methods, ]),
      collapse = "\n"
    )
  )
})

test_that("diagnostic print methods show compact status overviews", {
  ld <- structure(
    list(
      values = structure(
        list(
          tests = data.frame(
            item1_name = c("item01", "item02", "item03"),
            item2_name = c("item02", "item03", "item04"),
            p_value = c(0.003, 0.04, 0.20),
            converged = c(TRUE, FALSE, TRUE),
            stringsAsFactors = FALSE
          ),
          bh_critical_p = 0.01
        ),
        class = "gRm_local_independence_values"
      ),
      analysis = list(items = paste0("item", sprintf("%02d", 1:4)), exogenous = character())
    ),
    class = "gRm_local_dependence"
  )
  dif_result <- structure(
    list(
      values = structure(
        list(
          tests = data.frame(
            item_name = c("item01", "item02", "item03"),
            background_name = c("site", "site", "site"),
            p_value = c(0.003, 0.04, 0.20),
            converged = c(TRUE, FALSE, TRUE),
            output_stable = c(TRUE, FALSE, TRUE),
            stringsAsFactors = FALSE
          ),
          included_tests = data.frame(
            item_name = "item04",
            background_name = "site",
            stringsAsFactors = FALSE
          ),
          bh_critical_p = 0.01
        ),
        class = "gRm_dif_tests_values"
      ),
      analysis = list(items = paste0("item", sprintf("%02d", 1:4)), exogenous = "site")
    ),
    class = "gRm_dif"
  )

  expect_equal(capture.output(print(ld)), c(
    "gRm: Local dependence tests",
    "",
    "  Candidate pairs: 3",
    "  Selected by BH FDR 0.05: 1",
    "  Non-converged fits: 1",
    "  BH threshold: 0.01",
    "",
    "Use summary(x) to show the complete test table."
  ))
  expect_equal(capture.output(print(dif_result)), c(
    "gRm: DIF tests",
    "",
    "  Candidate tests: 3",
    "  Included DIF terms: 1",
    "  Selected by BH FDR 0.05: 1",
    "  Non-converged fits: 1",
    "  Output-stable fits: 2",
    "  BH threshold: 0.01",
    "",
    "Use summary(x) to show the complete test table."
  ))
})

test_that("fit print method shows a compact status overview", {
  fitted <- structure(
    list(
      spec = list(
        model_type = "gllrm",
        ld = data.frame(item1 = "I1", item2 = "I2"),
        dif = data.frame(item = "I3", exogenous = "site")
      ),
      convergence = list(
        converged = TRUE,
        iterations = 184L,
        delta = 0.0000023
      ),
      values = list(
        log_likelihood = 1234.56
      )
    ),
    class = "gRm_fit"
  )

  expect_equal(capture.output(print(fitted)), c(
    "gRm: Graphical Log-Linear Rasch Model fit",
    "",
    "  Converged: yes",
    "  Iterations: 184",
    "  Delta: 2.3e-06",
    "  Negative log likelihood (DIGRAM): 1234.56",
    "",
    "Use summary(x) to show the fitted-model details."
  ))
})

test_that("analysis summary prints a Graphical Log-Linear Rasch Model header", {
  analysis <- summary_print_analysis()
  printed <- capture.output(print(summary(analysis)))

  expect_equal(printed[seq_len(18L)], c(
    "gRm: Graphical Log-Linear Rasch Model analysis",
    "",
    "Data",
    "  Source: source_data",
    "  Rows: 54",
    "  ID: ID",
    "",
    "Variables",
    "  Items: 3 (I1, I2, I3)",
    "  Item levels: I1=3, I2=3, I3=3",
    "  Exogenous: 1 (site)",
    "  Exogenous levels: site=2",
    "",
    "Score groups",
    "  Groups: 2 (0-2, 3-6)",
    "",
    "Score group distribution",
    "Observed score range: 0-6"
  ))
  expect_true("Score group distribution" %in% printed)
  expect_true("Score-group cases: 54" %in% printed)
  expect_true("Missing item-score rows: 0" %in% printed)
  expect_false("data" %in% printed)

  summary_object <- summary(analysis)
  expect_equal(summary_object$data$items, "I1, I2, I3")
  expect_equal(summary_object$score_groups$score, c("0-2", "3-6", "Total"))
  expect_equal(summary_object$score_groups$count, c(20L, 34L, 54L))
})

test_that("analysis print shows a compact status overview", {
  analysis <- summary_print_analysis()

  expect_equal(capture.output(print(analysis)), c(
    "gRm: Graphical Log-Linear Rasch Model analysis",
    "",
    "  Source: source_data",
    "  Rows: 54",
    "  ID: ID",
    "  Items: 3",
    "  Exogenous: 1",
    "  Score groups: 2",
    "",
    "Use summary(x) to show variable names, level counts, and score-group distribution."
  ))
})

test_that("analysis summary has one output surface without which", {
  analysis <- summary_print_analysis()

  default <- summary(analysis)

  expect_equal(names(default$tables), c("data", "score_groups"))
  expect_equal(default$which, c("data", "score_groups"))
  expect_error(summary(analysis, which = "data"), "does not accept `which`")
})

test_that("model and fit summaries use statistical labels", {
  analysis <- summary_print_analysis()
  rasch_model <- gllrm(analysis)
  gllrm_model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:site)
  rasch_fit <- fit(rasch_model, max_step = 50L)
  gllrm_fit <- fit(gllrm_model, max_step = 50L)

  expect_output(print(rasch_model), "gRm: Rasch model specification", fixed = TRUE)
  expect_output(print(gllrm_model), "gRm: Graphical Log-Linear Rasch Model specification", fixed = TRUE)
  expect_output(print(rasch_fit), "gRm: Rasch model fit", fixed = TRUE)
  expect_output(print(gllrm_fit), "gRm: Graphical Log-Linear Rasch Model fit", fixed = TRUE)

  expect_summary_surface(
    rasch_model,
    "summary.gRm_model",
    title_pattern = "gRm: Rasch model specification"
  )
  expect_summary_surface(
    gllrm_model,
    "summary.gRm_model",
    title_pattern = "gRm: Graphical Log-Linear Rasch Model specification"
  )
  expect_summary_surface(
    rasch_fit,
    "summary.gRm_fit",
    title_pattern = "gRm: Rasch model fit"
  )
  expect_summary_surface(
    gllrm_fit,
    "summary.gRm_fit",
    title_pattern = "gRm: Graphical Log-Linear Rasch Model fit"
  )
})

test_that("model summary has one complete model surface without which", {
  model <- gllrm(summary_print_analysis(), ld = ~ I1:I2, dif = ~ I3:site)

  default <- summary(model)

  expect_equal(default$which, "model")
  expect_equal(names(default$tables), c("model", "ld", "dif"))
  expect_equal(names(default$model), c("model_type", "n_items", "n_exogenous", "n_ld", "n_dif"))
  printed <- capture.output(print(default))
  expect_true("  Model type: Graphical Log-Linear Rasch Model" %in% printed)
  expect_true("  Items: 3" %in% printed)
  expect_true("  Exogenous variables: 1" %in% printed)
  expect_true("  Local dependence terms: 1" %in% printed)
  expect_true("  DIF terms: 1" %in% printed)
  expect_false(any(grepl("model_type|n_items|n_exogenous|n_ld|n_dif", printed)))
  expect_true("Model" %in% printed)
  expect_false("model" %in% printed)
  expect_true("Local dependence terms" %in% printed)
  expect_true("  I1 -- I2" %in% printed)
  expect_true("DIF terms" %in% printed)
  expect_true("  I3 by site" %in% printed)
  expect_false("ld" %in% printed)
  expect_false("dif" %in% printed)
  expect_false(any(grepl("type item1 item2 source|type item exogenous source", printed)))
  expect_error(summary(model, which = "model"), "does not accept `which`")
  expect_error(summary(model, which = "ld"), "does not accept `which`")
  expect_error(summary(model, which = "dif"), "does not accept `which`")

  rasch_printed <- capture.output(print(summary(gllrm(summary_print_analysis()))))
  expect_equal(rasch_printed[rasch_printed %in% c("Local dependence terms", "DIF terms")], c(
    "Local dependence terms",
    "DIF terms"
  ))
  expect_equal(sum(rasch_printed == "  None"), 2L)
})

test_that("model print uses the complete model summary text", {
  model <- gllrm(summary_print_analysis(), ld = ~ I1:I2, dif = ~ I3:site)

  expect_equal(
    capture.output(print(model)),
    capture.output(print(summary(model)))
  )
  returned <- NULL
  capture.output(returned <- print(model))
  expect_identical(returned, model)
})

test_that("fit summary exposes fitted contents without fit status as a summary section", {
  fitted <- fit(gllrm(summary_print_analysis(), ld = ~ I1:I2, dif = ~ I3:site), max_step = 50L)

  default <- summary(fitted)
  parameters_requested <- summary(fitted, which = "parameters")
  thresholds_requested <- summary(fitted, which = "thresholds")

  expect_equal(default$which, c("parameters", "thresholds"))
  expect_equal(parameters_requested$which, "parameters")
  expect_equal(thresholds_requested$which, "thresholds")
  expect_equal(names(default$tables), c("parameters", "thresholds", "terms"))
  expect_equal(names(parameters_requested$tables), "parameters")
  expect_equal(names(thresholds_requested$tables), "thresholds")
  expect_equal(names(default$fit), c(
    "model_type", "log_likelihood", "n_parameters", "likelihood_n",
    "converged", "iterations", "delta"
  ))

  printed <- capture.output(print(default))
  expect_false("Fit" %in% printed)
  expect_false(any(grepl("^  Model type: ", printed)))
  expect_false(any(grepl("^  Converged: ", printed)))
  expect_false(any(grepl("^  Iterations: ", printed)))
  expect_false(any(grepl("^  Delta: ", printed)))
  expect_false(any(grepl("^  Log likelihood: ", printed)))
  expect_false(any(grepl("^  Likelihood rows: ", printed)))
  expect_true("Item parameters" %in% printed)
  expect_equal(names(default$parameters), c(
    "item", "location", "midpoint", "target", "info_at_target", "info_per_step"
  ))
  expect_equal(names(default$thresholds), c("item", "score", "threshold"))
  expect_true(any(grepl("Item +Location +Midpoint +Target score +Info at target +Info per step", printed)))
  expect_true(any(grepl("Item +Step +Score +Threshold", printed)))
  expect_false(any(grepl("info_at_target|info_per_step", printed)))
  expect_true(is.numeric(default$parameters$location))
  expect_true(is.numeric(default$thresholds$threshold))
  expect_true(is.integer(default$thresholds$score))
  expect_true("Model terms" %in% printed)
  expect_true("Local dependence terms" %in% printed)
  expect_true("  I1 -- I2" %in% printed)
  expect_true("DIF terms" %in% printed)
  expect_true("  I3 by site" %in% printed)
  expect_false(any(grepl("model_type|log_likelihood|n_parameters|likelihood_n", printed)))
  expect_false(any(grepl("type item1 item2 source|type item exogenous source", printed)))
  expect_false("Model terms" %in% capture.output(print(parameters_requested)))
  expect_false("Model terms" %in% capture.output(print(thresholds_requested)))
  expect_error(summary(fitted, which = "fit"), "which")
  expect_error(summary(fitted, which = "terms"), "which")
  expect_error(summary(fitted, which = "all"), "which")
})

test_that("summary which selects designed tables and rejects old accessor names", {
  fitted <- fit(gllrm(summary_print_analysis()), max_step = 50L)

  expect_false("item_parameters" %in% getNamespaceExports("gRm"))
  expect_summary_surface(
    local_dependence(fitted),
    "summary.gRm_local_dependence",
    title_pattern = "gRm: Local dependence tests"
  )
  expect_summary_surface(
    dif(fitted),
    "summary.gRm_dif",
    title_pattern = "gRm: DIF tests"
  )
  expect_summary_surface(
    global_homogeneity(fitted),
    "summary.gRm_global_homogeneity",
    "test",
    title_pattern = "gRm: Global homogeneity test"
  )

  expect_error(summary(fitted, which = "details"), "which")
  expect_error(summary(fitted, which = "tidy"), "which")
  expect_error(summary(fitted, which = "glance"), "which")
})

test_that("local dependence summary prints a compact BH-marked test table", {
  ld <- structure(
    list(
      values = structure(
        list(
          tests = data.frame(
            pair_label = c("I1I2", "I2I3"),
            item1_label = c("I1", "I2"),
            item2_label = c("I2", "I3"),
            item1_name = c("item01", "item02"),
            item2_name = c("item02", "item03"),
            chi_square = c(8.41234, 2.106),
            degrees_of_freedom = c(1L, 1L),
            p_value = c(0.003704, 0.1467),
            wpg_gamma = c(0.182345, 0.071234),
            converged = c(TRUE, TRUE),
            stop_reason = c(NA_character_, NA_character_),
            delta = c(0.00001234, 0.000009),
            stringsAsFactors = FALSE
          ),
          bh_critical_p = 0.012345
        ),
        class = "gRm_local_independence_values"
      )
    ),
    class = "gRm_local_dependence"
  )

  printed <- capture.output(print(summary(ld)))

  expect_equal(printed[[1L]], "gRm: Local dependence tests")
  expect_error(summary(ld, which = "tests"), "does not accept `which`")
  expect_false("tests" %in% printed)
  expect_true(any(grepl("Item 1 +Item 2 +Chisq +Df +Pr\\(>Chisq\\) +WPG +Converged +delta", printed)))
  expect_true(any(grepl("item01 +item02 +8\\.41 +1 +0\\.0037 +0\\.18230 +yes +1\\.234e-05 +\\*", printed)))
  expect_true(any(grepl("item02 +item03 +2\\.11 +1 +0\\.1467 +0\\.07123 +yes +9\\.000e-06", printed)))
  expect_true("---" %in% printed)
  expect_true("*: p <= Benjamini-Hochberg threshold for FDR = 0.05 (threshold = 0.0123)" %in% printed)
  expect_false(any(grepl("pair_label|item1_label|item2_label|TRUE|FALSE|Benjamini-Hochberg$", printed)))
})

test_that("local dependence summary prints all rows regardless of max.print", {
  n_rows <- 80L
  ld <- structure(
    list(
      values = structure(
        list(
          tests = data.frame(
            pair_label = paste0("I", seq_len(n_rows), "I", seq_len(n_rows) + 1L),
            item1_label = paste0("I", seq_len(n_rows)),
            item2_label = paste0("I", seq_len(n_rows) + 1L),
            item1_name = sprintf("item%03d", seq_len(n_rows)),
            item2_name = sprintf("item%03d", seq_len(n_rows) + 1L),
            chi_square = seq_len(n_rows) / 10,
            degrees_of_freedom = rep(1L, n_rows),
            p_value = seq_len(n_rows) / 1000,
            wpg_gamma = seq_len(n_rows) / 100,
            converged = rep(TRUE, n_rows),
            stop_reason = rep(NA_character_, n_rows),
            delta = seq_len(n_rows) / 100000,
            stringsAsFactors = FALSE
          ),
          bh_critical_p = 0.01
        ),
        class = "gRm_local_independence_values"
      )
    ),
    class = "gRm_local_dependence"
  )

  old_max <- getOption("max.print")
  options(max.print = 50L)
  on.exit(options(max.print = old_max), add = TRUE)

  printed <- capture.output(print(summary(ld)))

  expect_false(any(grepl("max.print|omitted", printed)))
  expect_true(any(grepl("item080 +item081", printed)))
})

test_that("DIF summary prints a compact BH-marked test table", {
  dif_result <- structure(
    list(
      values = structure(
        list(
          tests = data.frame(
            item_label = c("I1", "I2"),
            background_label = c("X1", "X1"),
            item_name = c("item01", "item02"),
            background_name = c("site", "site"),
            chi_square = c(8.41234, 2.106),
            degrees_of_freedom = c(1L, 1L),
            p_value = c(0.003704, 0.1467),
            p_chi = c(0.003704, 0.1467),
            gamma = c(0.182345, 0.071234),
            p_gamma = c(0.11, 0.22),
            gamma_source = c("item_screening", "item_screening"),
            test_type = c("no_dif", "no_dif"),
            status = c("tested", "tested"),
            converged = c(TRUE, FALSE),
            output_stable = c(TRUE, FALSE),
            delta = c(0.00001234, 0.000009),
            n_step = c(300L, 5000L),
            stop_reason = c(NA_character_, "source_periodic_checkpoint"),
            stringsAsFactors = FALSE
          ),
          included_tests = data.frame(),
          bh_critical_p = 0.012345
        ),
        class = "gRm_dif_tests_values"
      )
    ),
    class = "gRm_dif"
  )

  printed <- capture.output(print(summary(dif_result)))

  expect_equal(printed[[1L]], "gRm: DIF tests")
  expect_error(summary(dif_result, which = "tests"), "does not accept `which`")
  expect_false("tests" %in% printed)
  expect_true(any(grepl("Item +Exogenous +Chisq +Df +Pr\\(>Chisq\\) +Gamma +Converged +Stable +delta", printed)))
  expect_true(any(grepl("item01 +site +8\\.41 +1 +0\\.0037 +0\\.18230 +yes +yes +1\\.234e-05 +\\*", printed)))
  expect_true(any(grepl("item02 +site +2\\.11 +1 +0\\.1467 +0\\.07123 +no +no +9\\.000e-06", printed)))
  expect_true("---" %in% printed)
  expect_true("*: p <= Benjamini-Hochberg threshold for FDR = 0.05 (threshold = 0.0123)" %in% printed)
  expect_false(any(grepl("item_label|background_label|output_stable|n_step|p_chi|p_gamma|gamma_source|test_type|status|Benjamini-Hochberg$", printed)))
})

test_that("DIF summary prints all rows regardless of max.print", {
  n_rows <- 80L
  dif_result <- structure(
    list(
      values = structure(
        list(
          tests = data.frame(
            item_label = paste0("I", seq_len(n_rows)),
            background_label = rep("X1", n_rows),
            item_name = sprintf("item%03d", seq_len(n_rows)),
            background_name = rep("site", n_rows),
            chi_square = seq_len(n_rows) / 10,
            degrees_of_freedom = rep(1L, n_rows),
            p_value = seq_len(n_rows) / 1000,
            p_chi = seq_len(n_rows) / 1000,
            gamma = seq_len(n_rows) / 100,
            p_gamma = seq_len(n_rows) / 200,
            gamma_source = rep("item_screening", n_rows),
            test_type = rep("no_dif", n_rows),
            status = rep("tested", n_rows),
            converged = rep(TRUE, n_rows),
            output_stable = rep(FALSE, n_rows),
            delta = seq_len(n_rows) / 100000,
            n_step = seq_len(n_rows),
            stop_reason = rep(NA_character_, n_rows),
            stringsAsFactors = FALSE
          ),
          included_tests = data.frame(),
          bh_critical_p = 0.01
        ),
        class = "gRm_dif_tests_values"
      )
    ),
    class = "gRm_dif"
  )

  old_max <- getOption("max.print")
  options(max.print = 50L)
  on.exit(options(max.print = old_max), add = TRUE)

  printed <- capture.output(print(summary(dif_result)))

  expect_false(any(grepl("max.print|omitted", printed)))
  expect_true(any(grepl("item080 +site", printed)))
})

test_that("global homogeneity summary explains unbacked residual cells", {
  fitted <- fit(gllrm(summary_print_analysis()), max_step = 50L)
  gh <- global_homogeneity(fitted, max_step = 50L)

  gh_summary <- summary(gh)
  source_status <- attr(gh_summary, "source_status", exact = TRUE)
  gh_items <- gh$values$items

  expect_true(all(is.na(gh_items$residual)))
  expect_true(all(is.na(gh_items$marker)))
  expect_false(any(gh_items$residual_runtime_source_backed))
  expect_false(any(gh_items$marker_runtime_source_backed))
  expect_equal(source_status$item_residual_status, "not_source_backed")
  expect_equal(source_status$item_marker_status, "not_source_backed")
  expect_false(source_status$item_residual_runtime_source_backed)
  expect_false(source_status$item_marker_runtime_source_backed)

  printed <- capture.output(print(gh_summary))
  expect_true(any(grepl("Item means", printed, fixed = TRUE)))
  expect_true(any(grepl("not source-backed in gRm and are not printed", printed, fixed = TRUE)))
  expect_false(any(grepl("residual_runtime_source_backed", printed, fixed = TRUE)))
  expect_false(any(grepl("marker_runtime_source_backed", printed, fixed = TRUE)))
})

test_that("global homogeneity print shows compact status", {
  fitted <- fit(gllrm(summary_print_analysis()), max_step = 50L)
  gh <- global_homogeneity(fitted, max_step = 50L)
  printed <- capture.output(print(gh))

  expect_equal(printed[seq_len(2L)], c(
    "gRm: Global homogeneity test",
    ""
  ))
  expect_true(any(grepl("^  Score groups: 2$", printed)))
  expect_true(any(grepl("^  Parameters: 5$", printed)))
  expect_true(any(grepl("^  CLR: ", printed)))
  expect_true(any(grepl("^  Df: 5$", printed)))
  expect_true(any(grepl("^  Pr\\(>CLR\\): ", printed)))
  expect_true(any(grepl("^  Non-converged group fits: ", printed)))
  expect_true(any(grepl("Use summary\\(x\\) to show the global test, score groups, and item means.", printed)))
})

test_that("global homogeneity print and summary expose uniform interaction sections", {
  gh <- structure(
    list(
      values = structure(
        list(
          summary = list(
            n_groups = 2L,
            n_parameters = 38L,
            full_log_likelihood = -1234.5,
            subgroup_log_likelihood_sum = -1012.7,
            clr = 443.6,
            df = 38L,
            p_value = 1e-8
          ),
          score_groups = data.frame(
            group = c(1L, 2L),
            from_score = c(0L, 19L),
            to_score = c(18L, 53L),
            n = c(150L, 155L),
            log_likelihood = c(-500.1, -512.2),
            converged = c(TRUE, TRUE),
            delta = c(0.000091, 0.000087),
            stringsAsFactors = FALSE
          ),
          items = data.frame(),
          fit = list(
            context = list(
              items = data.frame(
                label_code = c("b", "c", "d"),
                name = c("Item02", "Item03", "Item04"),
                stringsAsFactors = FALSE
              ),
              backgrounds = data.frame(
                label_code = c("f", "g"),
                name = c("site", "cohort"),
                stringsAsFactors = FALSE
              )
            )
          ),
          uniform_ld = data.frame(
            item1_label = "b",
            item2_label = "c",
            observed_gamma = I(list(c(0.25, 0.24))),
            expected_gamma = I(list(c(0.36, 0.34))),
            chi_square = 97.3,
            df = 9L,
            p_value = 1e-8,
            stringsAsFactors = FALSE
          ),
          uniform_dif = data.frame(
            item_label = c("b", "d"),
            background_label = c("f", "g"),
            observed_gamma = I(list(c(0.31, 0.00), c(0.34, 0.05))),
            expected_gamma = I(list(c(0.19, 0.05), c(0.22, 0.21))),
            chi_square = c(100.9, 63.7),
            df = c(6L, 3L),
            p_value = c(1e-8, 2e-5),
            stringsAsFactors = FALSE
          )
        ),
        class = "gRm_global_homogeneity_values"
      )
    ),
    class = "gRm_global_homogeneity"
  )

  printed_result <- capture.output(print(gh))
  expect_true(any(grepl("^  Uniform LD tests: 1$", printed_result)))
  expect_true(any(grepl("^  Uniform DIF tests: 2$", printed_result)))
  expect_true(any(grepl("uniform interaction tests", printed_result, fixed = TRUE)))

  gh_summary <- summary(gh)
  printed <- capture.output(print(gh_summary))

  expect_equal(gh_summary$which, c("test", "score_groups", "item_means", "uniform_ld", "uniform_dif"))
  expect_true("Uniform local dependence" %in% printed)
  expect_true("Uniform DIF" %in% printed)
  expect_true(any(grepl("Obs gamma 0-18", printed, fixed = TRUE)))
  expect_true(any(grepl("Exp gamma 19-53", printed, fixed = TRUE)))
  expect_true(any(grepl("Item02 +Item03 +0\\.25 +0\\.36 +0\\.24 +0\\.34", printed)))
  expect_true(any(grepl("97\\.3 +9 +1e-08", printed)))
  expect_true(any(grepl("Item02 +site +0\\.31 +0\\.19 +0\\.00 +0\\.05", printed)))
  expect_true(any(grepl("101\\.0 +6 +1e-08", printed)))
  expect_true(any(grepl("Item04 +cohort +0\\.34 +0\\.22 +0\\.05 +0\\.21", printed)))
  expect_true(any(grepl("63\\.7 +3 +2e-05", printed)))
  expect_false(any(grepl("Benjamini-Hochberg|\\*", printed)))
})

test_that("global homogeneity default summary omits empty uniform sections", {
  gh <- structure(
    list(
      values = structure(
        list(
          summary = list(n_groups = 2L, n_parameters = 5L, clr = 4.2, df = 5L, p_value = 0.5),
          score_groups = data.frame(
            group = c(1L, 2L),
            from_score = c(0L, 3L),
            to_score = c(2L, 6L),
            n = c(4L, 8L),
            log_likelihood = c(3.3, 12.1),
            converged = c(TRUE, TRUE),
            delta = c(0.00001, 0.00002),
            stringsAsFactors = FALSE
          ),
          items = data.frame(),
          uniform_ld = data.frame(),
          uniform_dif = data.frame()
        ),
        class = "gRm_global_homogeneity_values"
      )
    ),
    class = "gRm_global_homogeneity"
  )

  gh_summary <- summary(gh)
  printed <- capture.output(print(gh_summary))

  expect_equal(gh_summary$which, c("test", "score_groups", "item_means"))
  expect_false("uniform_ld" %in% names(gh_summary$tables))
  expect_false("uniform_dif" %in% names(gh_summary$tables))
  expect_false("Uniform local dependence" %in% printed)
  expect_false("Uniform DIF" %in% printed)
  expect_s3_class(summary(gh, which = "uniform_ld"), "summary.gRm_global_homogeneity")
  expect_s3_class(summary(gh, which = "uniform_dif"), "summary.gRm_global_homogeneity")
})

test_that("global homogeneity residual restriction is documented", {
  api_docs <- readLines(repo_path("gRm", "R", "api-results.R"), warn = FALSE)
  package_docs <- readLines(repo_path("gRm", "R", "gRm-package.R"), warn = FALSE)

  expect_true(any(grepl("Global homogeneity residuals", api_docs, fixed = TRUE)))
  expect_true(any(grepl("not_source_backed", api_docs, fixed = TRUE)))
  expect_true(any(grepl("Global homogeneity residuals", package_docs, fixed = TRUE)))
  expect_true(any(grepl("not source-backed", package_docs, fixed = TRUE)))
})

test_that("item fit tests and items expose distinct public tables and backend values", {
  fitted <- fit(gllrm(summary_print_analysis()), max_step = 50L)
  tests <- item_fit(fitted, include_extended = TRUE)
  items <- item_fit(fitted, which = "items", include_extended = TRUE)
  backend <- attr(tests, "values", exact = TRUE)$items
  item_backend <- attr(items, "values", exact = TRUE)$extended$summaries
  bh <- attr(tests, "bh", exact = TRUE)

  expect_s3_class(tests, "gRm_item_fit")
  expect_s3_class(tests, "gRm_direct_table")
  expect_s3_class(items, "gRm_item_fit")
  expect_s3_class(items, "gRm_direct_table")
  expect_true(is.data.frame(tests))
  expect_true(is.data.frame(items))
  expect_equal(nrow(items), nrow(tests))
  expect_equal(attr(tests, "which", exact = TRUE), "tests")
  expect_equal(attr(items, "which", exact = TRUE), "items")
  expect_s3_class(attr(tests, "values", exact = TRUE), "gRm_item_fits_values")
  expect_null(attr(tests, "analysis", exact = TRUE))
  expect_null(attr(tests, "fit", exact = TRUE))
  expect_null(attr(items, "analysis", exact = TRUE))
  expect_null(attr(items, "fit", exact = TRUE))
  expect_true(is.data.frame(bh))
  expect_true(all(c("fdr_5", "fdr_1") %in% bh$threshold))
  expect_false(" " %in% names(items))

  expect_equal(names(tests), c(
    "Item",
    "Outfit",
    "Outfit SE",
    "Pr(>Outfit)",
    "Outfit FDR",
    "Infit",
    "Infit SE",
    "Pr(>Infit)",
    "Infit FDR",
    "Observed gamma",
    "Expected gamma",
    "Gamma SE",
    "Pr(>Gamma)",
    "Gamma FDR",
    "Gamma direction"
  ))
  expect_false(any(c(
    "item_label", "item_name", "outfit_fdr", "infit_fdr", "gamma_fdr",
    "Selected", " "
  ) %in% names(tests)))
  expect_true(all(c(
    "item_label", "item_name", "outfit_fdr", "infit_fdr", "gamma_fdr"
  ) %in% names(backend)))
  expect_true(is.numeric(tests$Outfit))
  expect_true(is.numeric(tests[["Pr(>Outfit)"]]))
  expect_true(is.numeric(tests[["Observed gamma"]]))
  expect_true(all(tests[["Outfit FDR"]] %in% c("", "*", "**", "***")))
  expect_true(all(tests[["Infit FDR"]] %in% c("", "*", "**", "***")))
  expect_true(all(tests[["Gamma FDR"]] %in% c("", "*", "**", "***")))

  expect_equal(names(items), c(
    "Item",
    "Outfit N",
    "Outfit observed",
    "Outfit expected",
    "Outfit total",
    "Infit observed",
    "Infit expected",
    "Infit variance",
    "Infit ratio"
  ))
  expect_false(any(c("item_label", "item_name") %in% names(items)))
  expect_true(all(c(
    "item_label", "item_name",
    "outfit_total_n", "outfit_total_observed",
    "outfit_total_expected", "outfit_total_value",
    "infit_observed", "infit_expected",
    "infit_variance", "infit_value"
  ) %in% names(item_backend)))
  expect_true(is.numeric(items[["Outfit observed"]]))
  expect_true(is.numeric(items[["Infit ratio"]]))

  expect_false("p_gamma" %in% names(items))
  expect_false("outfit_total_value" %in% names(tests))
  expect_error(item_fit(fitted, which = "bh"), "which")
})

test_that("item fit output prints the selected table directly", {
  fitted <- fit(gllrm(summary_print_analysis()), max_step = 50L)

  tests <- item_fit(fitted, include_extended = TRUE)
  items <- item_fit(fitted, which = "items", include_extended = TRUE)
  tests_printed <- capture.output(print(tests))
  items_printed <- capture.output(print(items))

  expect_equal(tests_printed[[1L]], "gRm: Item fit tests")
  expect_equal(items_printed[[1L]], "gRm: Item fit item diagnostics")
  expect_true(any(grepl("Item +Outfit +Outfit SE +Pr\\(>Outfit\\) +Infit +Infit SE +Pr\\(>Infit\\)", tests_printed)))
  expect_true(any(grepl("Observed gamma", tests_printed)))
  expect_true(any(grepl("Expected gamma +Gamma SE +Pr\\(>Gamma\\) +Gamma direction", tests_printed)))
  expect_true(any(grepl("Item +Outfit N +Outfit observed +Outfit expected +Outfit total", items_printed)))
  expect_true(any(grepl("Infit observed", items_printed)))
  expect_true(any(grepl("Infit expected +Infit variance +Infit ratio", items_printed)))
  expect_true(any(grepl("Benjamini-Hochberg thresholds: FDR 0\\.05 = .*FDR 0\\.01 =", tests_printed)))
  expect_false(any(grepl("item_label|item_name|outfit_fdr|infit_fdr|gamma_fdr|Selected", tests_printed)))
  expect_false(any(grepl("item_label|item_name|infit_value|outfit_total_value", items_printed)))
  expect_true(any(grepl("Stars are based on separate Benjamini-Hochberg adjustments across items", tests_printed, fixed = TRUE)))
  expect_true(any(grepl("* = selected at 5% FDR", tests_printed, fixed = TRUE)))
  expect_true(any(grepl("*** = selected at 0.1% FDR", tests_printed, fixed = TRUE)))
  expect_true(any(grepl("The thresholds below are based on all item-fit p-values pooled together", tests_printed, fixed = TRUE)))
  expect_false(any(grepl("Benjamini-Hochberg threshold", items_printed)))
  expect_false(any(grepl("gRm_item_fit result", tests_printed, fixed = TRUE)))
})

test_that("item fit exposes each source FDR grade beside its own p-value", {
  values <- list(
    items = data.frame(
      item_name = paste0("I", 1:4),
      outfit = rep(1, 4),
      outfit_sd = rep(0.1, 4),
      p_outfit = c(0.4, 0.04, 0.004, 0.0004),
      outfit_fdr = 0:3,
      infit = rep(1, 4),
      infit_sd = rep(0.1, 4),
      p_infit = c(0.04, 0.4, 0.0004, 0.004),
      infit_fdr = c(1L, 0L, 3L, 2L),
      observed_gamma = c(-0.2, 0.2, -0.2, 0.2),
      expected_gamma = rep(0, 4),
      gamma_sd = rep(0.1, 4),
      p_gamma = c(0.0004, 0.004, 0.04, 0.4),
      gamma_fdr = c(3L, 2L, 1L, 0L),
      direction = c("low", "high", "low", "high"),
      stringsAsFactors = FALSE
    ),
    bh_limits = c(fdr_5 = 0.004, fdr_1 = 0.0004)
  )

  table <- gRm:::item_fit_tests_table(values)

  expect_equal(table[["Outfit FDR"]], c("", "*", "**", "***"))
  expect_equal(table[["Infit FDR"]], c("*", "", "***", "**"))
  expect_equal(table[["Gamma FDR"]], c("***", "**", "*", ""))
  expect_equal(table$`Gamma direction`[[1L]], "low")
  expect_equal(table$`Outfit FDR`[[1L]], "")

  printed <- capture.output(gRm:::print_item_fit_tests_table(table))
  expect_true(any(grepl("4e-04***", printed, fixed = TRUE)))
  expect_false(any(grepl("Outfit FDR|Infit FDR|Gamma FDR", printed)))
  item_lines <- printed[grepl("^\\s+I[1-4]\\s", printed)]
  outfit_p_positions <- vapply(
    item_lines,
    function(line) unname(regexpr("4e-0[1-4]", line, perl = TRUE)[[1L]]),
    integer(1L)
  )
  expect_length(item_lines, 4L)
  expect_equal(length(unique(outfit_p_positions)), 1L)
})

test_that("screen summary and score-effect table print BH selection context", {
  analysis <- summary_print_analysis()
  screen_result <- screen(analysis, inference = "asymptotic")
  effects <- score_effects(analysis)

  screen_summary <- expect_summary_surface(
    screen_result,
    "summary.gRm_screen",
    title_pattern = "gRm: SCREEN J tests"
  )
  screen_printed <- capture.output(print(screen_result))
  summary_printed <- capture.output(print(screen_summary))
  expect_s3_class(effects, "gRm_score_effects")
  expect_s3_class(effects, "gRm_direct_table")
  expect_true(is.data.frame(effects))
  expect_false(exists("summary.gRm_score_effects", envir = asNamespace("gRm"), inherits = FALSE))
  expect_false(exists("print.gRm_score_effects", envir = asNamespace("gRm"), inherits = FALSE))
  expect_equal(screen_printed[[1L]], "gRm: SCREEN J tests")
  expect_true(any(grepl("Tested relations:", screen_printed, fixed = TRUE)))
  expect_true(any(grepl("Use summary\\(x\\) to show the SCREEN J test tables\\.", screen_printed)))
  expect_true(any(grepl("Local dependence", summary_printed, fixed = TRUE)))
  expect_true(any(grepl("DIF", summary_printed, fixed = TRUE)))
  expect_true(any(grepl("Score effects", summary_printed, fixed = TRUE)))
  expect_true(any(grepl("Pr\\(>Chisq\\)", summary_printed)))
  expect_true(any(grepl("selected by the SCREEN J source decision path at the 5% level", summary_printed, fixed = TRUE)))
  expect_true(is.data.frame(screen_summary$selected))
  expect_true(is.data.frame(attr(screen_summary, "bh", exact = TRUE)))
  expect_true(is.data.frame(attr(screen_summary, "model_terms", exact = TRUE)))
  expect_error(summary(screen_result, which = "selected"), "does not accept `which`")
  expect_output(print(effects), "Benjamini-Hochberg thresholds")
  expect_output(print(effects), "gRm: Score-effect tests", fixed = TRUE)
  expect_equal(names(effects), c(
    "Exogenous", "Hypothesis", "Chisq", "Df", "Pr(>Chisq)", "Gamma",
    "Pr(Gamma+)", "Pr(|Gamma|)"
  ))
  expect_false(any(grepl("exo_label|exo_name|selected|chi_marker|gamma_marker", capture.output(print(effects)))))
  expect_true(is.data.frame(attr(effects, "selected", exact = TRUE)))
  expect_true(is.data.frame(attr(effects, "bh", exact = TRUE)))
  expect_s3_class(attr(effects, "values", exact = TRUE), "gRm_exo_select_values")
})
