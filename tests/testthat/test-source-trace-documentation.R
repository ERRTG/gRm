test_that("calculation-heavy GLLRM files carry Pascal source traces", {
  required_patterns <- list(
    "R/gllrm_context.R" = c(
      "Initialize_GLLRMinfo",
      "Estimate_GLLRM",
      "Implementation guard for explicit GLLRM component enumeration",
      "source/PAS_skunits/SKTypes.pas",
      "source/PAS_skunits/skbias22.pas::LD_Gamma_calculation"
    ),
    "R/gllrm_components.R" = c(
      "LD_Gamma_calculation",
      "Gamma_calculation",
      "CalculateBiasedGammaValues2",
      "source/PAS_skunits/skbias12b.pas::InitializeParameters"
    ),
    "R/gllrm_fit.R" = c(
      "CalculateBiasedGammaValues2",
      "Find_new_IJparameters",
      "Find_new_IXparameters",
      "Adjust_IJparameters",
      "Adjust_IXparameters",
      "GLLRM_estim"
    ),
    "R/gllrm_values.R" = c(
      "GLLRM_output",
      "PREPARE_REAL_GAMMA_STATISTICS",
      "source/PAS_skunits/skfit2.pas::Standardize_tab4"
    ),
    "R/source_gamma_stats.R" = c(
      "source/PAS_skunits/skfit2.pas::Standardize_tab4",
      "first source table index",
      "30 fixed row/column scaling passes"
    ),
    "R/gllrm_probability_cache.R" = c(
      "Implementation-only",
      "does not change the source algorithm"
    ),
    "src/gllrm_expected.cpp" = c(
      "source/GLLRM_ESTIM.txt::",
      "CalculateBiasedGammaValues2",
      "source/PAS_skunits/skbias12b.pas::",
      "Estimate_GLLRM",
      "source/PAS_skunits/skbias22.pas::Gamma_calculation",
      "source/PAS_skunits/skbias22.pas::LD_Gamma_calculation"
    ),
    "R/m2_m3_specs.R" = c(
      "source/PAS_skunits/skbias14.pas::Prepare_CM3tests",
      "all CM2 rows",
      "appends CM3 rows"
    ),
    "R/m2_m3_counts.R" = c(
      "source/PAS_skunits/skbias14.pas::Count_IJtable",
      "Count_IXtable",
      "Count_IStable",
      "Count_IJK",
      "Count_IJXtable",
      "Count_IXZtable",
      "Count_IJStable",
      "Count_IXStable"
    ),
    "R/m2_m3_expected.R" = c(
      "source/PAS_skunits/SKbias2.pas::calculate_expected_IJ_table",
      "calculate_expected_IX_table",
      "calculate_expected_IJK_table",
      "calculate_EIJX_table",
      "calculate_EIXZ_table",
      "Calculate_Cprob1",
      "Calculate_Cprob2",
      "Calculate_Cprob3"
    ),
    "R/m2_m3_values.R" = c(
      "source/PAS_skunits/skbias14.pas::Twoway_analysis",
      "Threeway_analysis",
      "CM2 aggregate",
      "CM3 aggregate",
      "PFCHI"
    )
  )

  for (relative_path in names(required_patterns)) {
    path <- repo_path("gRm", relative_path)
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    for (pattern in required_patterns[[relative_path]]) {
      expect_match(
        text,
        pattern,
        fixed = TRUE,
        info = paste(relative_path, "should document", pattern)
      )
    }
  }
})

test_that("included DIF, included LD, and screen J helpers carry local source traces", {
  required_patterns <- list(
    "R/dif_tests_values.R" = c(
      "CHECK D",
      "Find_new_IXparameters",
      "Adjust_IXparameters"
    ),
    "R/local_independence_values.R" = c(
      "MissingLD",
      "Find_new_IJparameters",
      "Adjust_IJparameters"
    ),
    "R/screen_j_values.R" = "XYZ_bias_ANALYSE",
    "R/screen_j_conditional.R" = c(
      "source/PAS_skunits/SKxyz1.PAS::MAKE_XYZ_TABLE",
      "source/PAS_skunits/SKbigtab.pas::Transfer_BT_to_XYZ_TABLE"
    ),
    "R/screen_j_score_effects.R" = "StepwiseScoreScreening",
    "R/global_homogeneity_ld.R" = c(
      "gllrm_uniform_summary_stats",
      "degrees-of-freedom floor",
      "requires source or validator evidence before changing"
    ),
    "R/global_homogeneity_dif.R" = c(
      "source/PAS_skunits/skfit2.pas::Standardize_ETAB2_to_TAB2_margins",
      "source/PAS_skunits/skfit2.pas::Standardize_tab4"
    )
  )

  for (relative_path in names(required_patterns)) {
    path <- repo_path("gRm", relative_path)
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    for (pattern in required_patterns[[relative_path]]) {
      expect_match(
        text,
        pattern,
        fixed = TRUE,
        info = paste(relative_path, "should document", pattern)
      )
    }
  }
})

test_that("source bundle manifest counters document Pascal counter names", {
  source_bundle_text <- paste(readLines(repo_path("gRm", "R", "source_bundle.R"), warn = FALSE), collapse = "\n")

  expect_match(source_bundle_text, "Nincomplete", fixed = TRUE)
  expect_match(source_bundle_text, "Nuseless", fixed = TRUE)
  expect_match(source_bundle_text, "nmissing_items = classified$n_incomplete", fixed = TRUE)
  expect_match(source_bundle_text, "nmissing_backgrounds = classified$n_useless", fixed = TRUE)
})
