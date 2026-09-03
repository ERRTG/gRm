test_that("calculation-heavy GLLRM files carry Pascal source traces", {
  required_patterns <- list(
    "R/gllrm_context.R" = c(
      "Initialize_GLLRMinfo",
      "Estimate_GLLRM",
      "Implementation guard for explicit GLLRM component enumeration",
      "source/digram_source_20260817/skunits/SKTypes.pas",
      "source/digram_source_20260817/skunits/skbias22.pas::LD_Gamma_calculation"
    ),
    "R/gllrm_components.R" = c(
      "LD_Gamma_calculation",
      "Gamma_calculation",
      "CalculateBiasedGammaValues2",
      "source/digram_source_20260817/skunits/skbias12b.pas::InitializeParameters"
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
      "source/digram_source_20260817/skunits/skfit2.pas::Standardize_tab4"
    ),
    "R/source_gamma_stats.R" = c(
      "source/digram_source_20260817/skunits/skfit2.pas::Standardize_tab4",
      "first source table index",
      "30 fixed row/column scaling passes"
    ),
    "R/gllrm_probability_cache.R" = c(
      "Implementation-only",
      "does not change the source algorithm"
    ),
    "src/gllrm_expected.cpp" = c(
      "source/digram_source_20260817/skunits/skbias22.pas::",
      "CalculateBiasedGammaValues2",
      "source/digram_source_20260817/skunits/skbias12b.pas::",
      "Estimate_GLLRM",
      "source/digram_source_20260817/skunits/skbias22.pas::Gamma_calculation",
      "source/digram_source_20260817/skunits/skbias22.pas::LD_Gamma_calculation"
    ),
    "R/cm2_cm3_specs.R" = c(
      "source/digram_source_20260817/skunits/skbias14.pas::Prepare_CM3tests",
      "all CM2 rows",
      "appends CM3 rows"
    ),
    "R/cm2_cm3_counts.R" = c(
      "source/digram_source_20260817/skunits/skbias14.pas::Count_IJtable",
      "Count_IXtable",
      "Count_IStable",
      "Count_IJK",
      "Count_IJXtable",
      "Count_IXZtable",
      "Count_IJStable",
      "Count_IXStable"
    ),
    "R/cm2_cm3_expected.R" = c(
      "source/digram_source_20260817/skunits/SKbias2.pas::calculate_expected_IJ_table",
      "calculate_expected_IX_table",
      "calculate_expected_IJK_table",
      "calculate_EIJX_table",
      "calculate_EIXZ_table",
      "Calculate_Cprob1",
      "Calculate_Cprob2",
      "Calculate_Cprob3"
    ),
    "R/cm2_cm3_values.R" = c(
      "source/digram_source_20260817/skunits/skbias14.pas::Twoway_analysis",
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
      "source/digram_source_20260817/skunits/SKxyz1.PAS::MAKE_XYZ_TABLE",
      "source/digram_source_20260817/skunits/SKbigtab.pas::Transfer_BT_to_XYZ_TABLE"
    ),
    "R/screen_j_score_effects.R" = "StepwiseScoreScreening",
    "R/global_homogeneity_ld.R" = c(
      "gllrm_uniform_summary_stats",
      "degrees-of-freedom floor",
      "requires source or validator evidence before changing"
    ),
    "R/global_homogeneity_dif.R" = c(
      "source/digram_source_20260817/skunits/skfit2.pas::Standardize_ETAB2_to_TAB2_margins",
      "source/digram_source_20260817/skunits/skfit2.pas::Standardize_tab4"
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

test_that("repository source traces follow the exact scd.dpr unit paths", {
  dpr_path <- repo_path("source", "digram_source_20260817", "scd", "scd.dpr")
  if (!file.exists(dpr_path)) {
    skip("Canonical DIGRAM source tree is not present in this package check.")
  }

  dpr <- readLines(dpr_path, warn = FALSE, encoding = "latin1")
  unit_lines <- grep("[[:space:]]in[[:space:]]+'[^']+[.]pas'", dpr,
    value = TRUE, ignore.case = TRUE
  )
  unit_paths <- sub(
    ".*[[:space:]]in[[:space:]]+'([^']+[.]pas)'.*",
    "\\1",
    unit_lines,
    ignore.case = TRUE
  )
  unit_paths <- gsub("\\\\", "/", unit_paths)
  unit_dirs <- ifelse(
    startsWith(tolower(unit_paths), "../skunits/"),
    "skunits",
    ifelse(!grepl("/", unit_paths, fixed = TRUE), "scd", NA_character_)
  )
  keep <- !is.na(unit_dirs)
  selected_by_name <- setNames(unit_dirs[keep], tolower(basename(unit_paths[keep])))
  selected_path_by_name <- setNames(
    file.path(
      "source",
      "digram_source_20260817",
      unit_dirs[keep],
      basename(unit_paths[keep])
    ),
    tolower(basename(unit_paths[keep]))
  )

  expect_identical(unname(selected_by_name[["digram1f.pas"]]), "scd")
  expect_identical(unname(selected_by_name[["skbias14.pas"]]), "skunits")
  expect_identical(unname(selected_by_name[["skstat.pas"]]), "skunits")

  validation_trace_files <- list.files(
    repo_path("validation"),
    pattern = "[.](R|md)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  # Durable/generated validation artifacts can be very large and merely copy
  # source-trace text already checked in code and maintained documentation.
  validation_trace_files <- validation_trace_files[
    !grepl("/validation-output/", validation_trace_files, fixed = TRUE) &
      !grepl("/validation/runs/", validation_trace_files, fixed = TRUE)
  ]
  trace_files <- c(
    list.files(repo_path("gRm", "R"), pattern = "[.]R$", full.names = TRUE),
    list.files(repo_path("gRm", "src"), pattern = "[.](c|cpp)$", full.names = TRUE),
    list.files(repo_path("gRm", "man"), pattern = "[.]Rd$", full.names = TRUE),
    repo_path("gRm", "ALGORITHMS.md"),
    repo_path("gRm", "tools", "numerical-source-trace.csv"),
    repo_path("README.md"),
    repo_path("AGENTS.md"),
    list.files(repo_path("docs"), pattern = "[.]md$", recursive = TRUE, full.names = TRUE),
    list.files(
      repo_path("pascal_harness"),
      pattern = "[.](pas|md)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    ),
    validation_trace_files
  )
  trace_text <- paste(
    unlist(lapply(trace_files, readLines, warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )
  citation_pattern <- paste0(
    "source/digram_source_20260817/(scd|skunits)/",
    "[^[:space:]`\\\"'():,;]+[.]pas"
  )
  citations <- unique(unlist(
    regmatches(trace_text, gregexpr(citation_pattern, trace_text,
      perl = TRUE, ignore.case = TRUE
    )),
    use.names = FALSE
  ))
  expect_gt(length(citations), 0L)

  cited_dirs <- sub(
    "^source/digram_source_20260817/(scd|skunits)/.*$",
    "\\1",
    citations,
    ignore.case = TRUE
  )
  cited_names <- tolower(basename(citations))
  selected_dirs <- unname(selected_by_name[cited_names])
  mapped <- !is.na(selected_dirs)
  wrong <- citations[mapped & tolower(cited_dirs) != selected_dirs]
  # These two references deliberately identify copies that scd.dpr does not
  # select; their surrounding documentation labels them as unselected.
  allowed_unselected <- c(
    "source/digram_source_20260817/scd/skbias14.pas",
    "source/digram_source_20260817/skunits/DIGRAM1f.pas"
  )
  wrong <- setdiff(wrong, allowed_unselected)
  expect_equal(
    length(wrong),
    0L,
    info = paste("Trace(s) disagree with scd.dpr:", paste(wrong, collapse = ", "))
  )

  missing <- citations[!file.exists(file.path(repo_root(), citations))]
  expect_equal(
    length(missing),
    0L,
    info = paste("Canonical trace path(s) do not exist:", paste(missing, collapse = ", "))
  )

  pascal_source_cache <- new.env(parent = emptyenv())
  read_pascal_source <- function(relative_path) {
    source_path <- file.path(repo_root(), relative_path)
    if (!file.exists(source_path)) {
      candidates <- list.files(dirname(source_path), full.names = TRUE)
      matches <- candidates[
        tolower(basename(candidates)) == tolower(basename(source_path))
      ]
      if (length(matches) == 1L) {
        source_path <- matches[[1L]]
      }
    }
    cache_key <- normalizePath(source_path, mustWork = FALSE)
    if (!exists(cache_key, envir = pascal_source_cache, inherits = FALSE)) {
      source_text <- paste(
        readLines(source_path, warn = FALSE, encoding = "latin1"),
        collapse = "\n"
      )
      assign(cache_key, source_text, envir = pascal_source_cache)
    }
    get(cache_key, envir = pascal_source_cache, inherits = FALSE)
  }

  routine_pattern <- paste0(citation_pattern, "::[A-Za-z_][A-Za-z0-9_]*")
  routine_citations <- unique(unlist(
    regmatches(trace_text, gregexpr(routine_pattern, trace_text,
      perl = TRUE, ignore.case = TRUE
    )),
    use.names = FALSE
  ))
  missing_routines <- routine_citations[!vapply(routine_citations, function(citation) {
    parts <- strsplit(citation, "::", fixed = TRUE)[[1L]]
    source_text <- read_pascal_source(parts[[1L]])
    grepl(tolower(parts[[2L]]), tolower(source_text), fixed = TRUE)
  }, logical(1L))]
  expect_equal(
    length(missing_routines),
    0L,
    info = paste(
      "Routine trace(s) are absent from the cited canonical unit:",
      paste(missing_routines, collapse = ", ")
    )
  )

  relative_routine_pattern <- paste0(
    "(?<![/A-Za-z0-9_])",
    "[A-Za-z_][A-Za-z0-9_]*[.]pas::[A-Za-z_][A-Za-z0-9_]*"
  )
  relative_routine_citations <- unique(unlist(
    regmatches(trace_text, gregexpr(relative_routine_pattern, trace_text,
      perl = TRUE, ignore.case = TRUE
    )),
    use.names = FALSE
  ))
  missing_relative_routines <- relative_routine_citations[!vapply(
    relative_routine_citations,
    function(citation) {
      parts <- strsplit(citation, "::", fixed = TRUE)[[1L]]
      selected_path <- unname(selected_path_by_name[tolower(parts[[1L]])])
      if (length(selected_path) != 1L || is.na(selected_path)) {
        return(TRUE)
      }
      source_text <- read_pascal_source(selected_path)
      grepl(tolower(parts[[2L]]), tolower(source_text), fixed = TRUE)
    },
    logical(1L)
  )]
  expect_equal(
    length(missing_relative_routines),
    0L,
    info = paste(
      "Relative routine trace(s) are absent from the scd.dpr-selected unit:",
      paste(missing_relative_routines, collapse = ", ")
    )
  )
})
