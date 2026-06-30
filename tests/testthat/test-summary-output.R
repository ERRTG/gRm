result_output_data <- function() {
  data.frame(
    ID = seq_len(12L),
    I1 = c(1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
    I2 = c(1L, 2L, 1L, 3L, 2L, 3L, 3L, 1L, 2L, 2L, 3L, 1L),
    I3 = c(2L, 1L, 3L, 1L, 2L, 3L, 2L, 1L, 3L, 3L, 2L, 1L),
    site = c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 2L)
  )
}

test_that("analysis summaries expose input metadata through summary", {
  ia <- gRm(
    result_output_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID"
  )

  out <- summary(ia)

  expect_s3_class(out, "summary.gRm_analysis")
  expect_true("data" %in% names(out$tables))
  expect_equal(out$data$n_items, 3L)
  expect_equal(out$data$n_exogenous, 1L)
  expect_error(summary(ia, which = "not_a_table"), "which")
})

test_that("screen summaries expose selected model terms and BH metadata", {
  ia <- gRm(result_output_data(), items = c("I1", "I2"), exogenous = "site", id = "ID")
  scr <- screen(ia, inference = "asymptotic")

  out <- summary(scr)
  bh <- attr(out, "bh", exact = TRUE)
  model_terms <- attr(out, "model_terms", exact = TRUE)

  expect_true(is.data.frame(out$local_dependence))
  expect_true(is.data.frame(out$dif))
  expect_true(is.data.frame(out$score_effects))
  expect_true(is.data.frame(out$selected))
  expect_true(is.data.frame(model_terms))
  expect_true(is.data.frame(bh))
  expect_equal(bh$fdr, c("0.05", "0.01", "0.001"))
  expect_output(print(out), "selected by the SCREEN J source decision path at the 5% level")
  expect_output(print(out), "global Benjamini-Hochberg threshold for FDR = 0.05")
  expect_error(summary(scr, which = "selected"), "does not accept `which`")
})

test_that("diagnostic summaries expose result tables through summary", {
  ia <- gRm(result_output_data(), items = c("I1", "I2"), exogenous = "site", id = "ID")
  fitted <- fit(gllrm(ia), max_step = 50L)

  ld <- local_dependence(fitted, jobs = 1L)
  dif_result <- dif(fitted, jobs = 1L)

  expect_true(is.data.frame(summary(ld)$tests))
  expect_true(is.data.frame(summary(dif_result)$tests))
  expect_error(summary(ld, which = "details"), "does not accept `which`")
})

test_that("diagnostic summaries remark on non-converged candidate fits", {
  ld <- structure(
    list(
      values = structure(
        list(
          tests = data.frame(
            pair_label = c("I1I2", "I1I3"),
            item1_label = c("I1", "I1"),
            item2_label = c("I2", "I3"),
            item1_name = c("I1", "I1"),
            item2_name = c("I2", "I3"),
            chi_square = c(1, 2),
            degrees_of_freedom = c(1L, 1L),
            p_value = c(0.5, 0.25),
            wpg_gamma = c(0.12, -0.34),
            converged = c(TRUE, FALSE),
            stop_reason = c(NA_character_, "source_periodic_checkpoint"),
            delta = c(0.00001, 3.5),
            stringsAsFactors = FALSE
          ),
          bh_critical_p = 0
        ),
        class = "gRm_local_independence_values"
      )
    ),
    class = "gRm_local_dependence"
  )

  out <- summary(ld)

  expect_equal(out$which, "tests")
  expect_equal(names(out$tables), "tests")
  expect_true(is.data.frame(out$selected))
  expect_true(is.data.frame(attr(out, "bh")))
  expect_equal(attr(out, "bh")$fdr, "0.05")
  expect_equal(attr(out, "bh")$p_value, 0)
  expect_equal(names(out$tests), c(
    "Item 1", "Item 2", "Chisq", "Df", "Pr(>Chisq)", "WPG",
    "Converged", "delta", " "
  ))
  expect_false(any(c("pair_label", "item1_label", "item2_label") %in% names(out$tests)))
  expect_equal(out$tests[["Item 1"]], c("I1", "I1"))
  expect_equal(out$tests[["Item 2"]], c("I2", "I3"))
  expect_equal(out$tests$Converged, c("yes", "no"))
  expect_equal(out$tests[[" "]], c("", ""))
  expect_true(is.data.frame(out$remarks))
  expect_equal(out$remarks$n_nonconverged, 1L)
  expect_match(out$remarks$message, "1 candidate fit did not converge")
  expect_match(out$remarks$message, "stop reasons: source_periodic_checkpoint=1", fixed = TRUE)
  expect_output(print(out), "Non-converged candidate fits: 1")
  expect_output(print(out), "source_periodic_checkpoint=1")
  expect_output(print(out), "\\*: p <= Benjamini-Hochberg threshold for FDR = 0.05")
  expect_false(any(grepl("^tests$", capture.output(print(out)))))
  expect_error(summary(ld, which = "selected"), "does not accept `which`")
  expect_error(summary(ld, which = "bh"), "does not accept `which`")
})

test_that("DIF summaries expose a compact public table with programmatic metadata", {
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
          included_tests = data.frame(
            item_name = "item03",
            background_name = "site",
            chi_square = 1.5,
            degrees_of_freedom = 1L,
            p_value = 0.22,
            stringsAsFactors = FALSE
          ),
          bh_critical_p = 0.012345
        ),
        class = "gRm_dif_tests_values"
      )
    ),
    class = "gRm_dif"
  )

  out <- summary(dif_result)

  expect_equal(out$which, "tests")
  expect_equal(names(out$tables), "tests")
  expect_true(is.data.frame(out$selected))
  expect_true(is.data.frame(out$included))
  expect_true(is.data.frame(attr(out, "bh")))
  expect_equal(attr(out, "bh")$fdr, "0.05")
  expect_equal(attr(out, "bh")$p_value, 0.012345)
  expect_equal(names(out$tests), c(
    "Item", "Exogenous", "Chisq", "Df", "Pr(>Chisq)", "Gamma",
    "Converged", "Stable", "delta", " "
  ))
  expect_false(any(c("output_stable", "n_step", "p_chi", "p_gamma", "gamma_source", "test_type", "status") %in% names(out$tests)))
  expect_equal(out$tests$Item, c("item01", "item02"))
  expect_equal(out$tests$Exogenous, c("site", "site"))
  expect_equal(out$tests$Converged, c("yes", "no"))
  expect_equal(out$tests$Stable, c("yes", "no"))
  expect_equal(out$tests[[" "]], c("*", ""))
  expect_output(print(out), "Non-converged candidate fits: 1")
  expect_output(print(out), "source_periodic_checkpoint=1")
  expect_output(print(out), "\\*: p <= Benjamini-Hochberg threshold for FDR = 0.05")
  expect_false(any(grepl("^tests$", capture.output(print(out)))))
  expect_error(summary(dif_result, which = "selected"), "does not accept `which`")
  expect_error(summary(dif_result, which = "included"), "does not accept `which`")
  expect_error(summary(dif_result, which = "bh"), "does not accept `which`")
})

test_that("fitted summaries expose item parameters without a separate item_parameters accessor", {
  analysis <- gRm(result_output_data(), items = c("I1", "I2", "I3"), exogenous = "site", id = "ID")
  fitted <- fit(gllrm(analysis), max_step = 50L)
  parameters <- summary(fitted, which = "parameters")
  thresholds <- summary(fitted, which = "thresholds")

  expect_false("item_parameters" %in% getNamespaceExports("gRm"))
  expect_false(exists("item_parameters", envir = asNamespace("gRm"), inherits = FALSE))
  expect_true(is.data.frame(parameters$parameters))
  expect_true(is.data.frame(thresholds$thresholds))
  expect_equal(names(parameters$parameters), c(
    "item", "location", "midpoint", "target", "info_at_target", "info_per_step"
  ))
  expect_equal(names(thresholds$thresholds), c("item", "score", "threshold"))
  expect_error(summary(fitted, which = "coefficients"), "which")
  expect_error(summary(fitted, which = "fit"), "which")
  expect_error(summary(fitted, which = "tests"), "which")
  expect_error(summary(fitted, which = "items"), "which")
})

test_that("global homogeneity summaries expose one test table plus groups and items", {
  analysis <- gRm(
    result_output_data(),
    items = c("I1", "I2", "I3"),
    exogenous = "site",
    id = "ID",
    score_cuts = c(1L, 6L)
  )
  fitted <- fit(gllrm(analysis), max_step = 50L)
  gh <- global_homogeneity(fitted, score_cuts = c(1L, 6L))

  tests <- summary(gh, which = "test")
  groups <- summary(gh, which = "score_groups")
  items <- summary(gh, which = "item_means")

  expect_equal(names(tests$test), c(
    "Score groups", "Parameters", "LogLik full", "LogLik groups", "CLR", "Df", "Pr(>CLR)"
  ))
  expect_equal(names(groups$score_groups), c("Score group", "Cases", "LogLik", "Converged", "delta"))
  expect_equal(names(items$item_means), c(
    "Score group", "Item", "Cases", "Observed mean", "Expected mean"
  ))
  expect_true(is.data.frame(tests$test))
  expect_true(is.data.frame(groups$score_groups))
  expect_true(is.data.frame(items$item_means))
  printed_tests <- capture.output(print(tests))
  expect_false(any(grepl("Benjamini-Hochberg", printed_tests, fixed = TRUE)))
  expect_error(summary(gh, which = "tests"), "which")
  expect_error(summary(gh, which = "groups"), "which")
  expect_error(summary(gh, which = "items"), "which")
  expect_error(summary(gh, which = "summary"), "which")
})

test_that("global homogeneity summaries expose uniform LD and DIF tables when present", {
  gh <- structure(
    list(
      values = structure(
        list(
          summary = list(n_groups = 2L, n_parameters = 38L, clr = 443.6, df = 38L, p_value = 1e-8),
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
                label_code = c("b", "c"),
                name = c("Item02", "Item03"),
                stringsAsFactors = FALSE
              ),
              backgrounds = data.frame(
                label_code = "f",
                name = "site",
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
            item_label = "b",
            background_label = "f",
            observed_gamma = I(list(c(0.31, 0.00))),
            expected_gamma = I(list(c(0.19, 0.05))),
            chi_square = 100.9,
            df = 6L,
            p_value = 1e-8,
            stringsAsFactors = FALSE
          )
        ),
        class = "gRm_global_homogeneity_values"
      )
    ),
    class = "gRm_global_homogeneity"
  )

  ld <- summary(gh, which = "uniform_ld")
  dif <- summary(gh, which = "uniform_dif")

  expect_equal(names(ld$uniform_ld), c(
    "Item 1", "Item 2",
    "Obs gamma 0-18", "Exp gamma 0-18",
    "Obs gamma 19-53", "Exp gamma 19-53",
    "Chisq", "Df", "Pr(>Chisq)"
  ))
  expect_equal(names(dif$uniform_dif), c(
    "Item", "Exogenous",
    "Obs gamma 0-18", "Exp gamma 0-18",
    "Obs gamma 19-53", "Exp gamma 19-53",
    "Chisq", "Df", "Pr(>Chisq)"
  ))
  expect_equal(ld$uniform_ld[["Item 1"]], "Item02")
  expect_equal(ld$uniform_ld[["Item 2"]], "Item03")
  expect_equal(dif$uniform_dif$Item, "Item02")
  expect_equal(dif$uniform_dif$Exogenous, "site")
  expect_equal(ld$uniform_ld[["Obs gamma 0-18"]], gh$values$uniform_ld$observed_gamma[[1L]][[1L]])
  expect_equal(ld$uniform_ld[["Exp gamma 19-53"]], gh$values$uniform_ld$expected_gamma[[1L]][[2L]])
  expect_equal(dif$uniform_dif[["Obs gamma 19-53"]], gh$values$uniform_dif$observed_gamma[[1L]][[2L]])
  expect_equal(dif$uniform_dif[["Pr(>Chisq)"]], gh$values$uniform_dif$p_value)
})
