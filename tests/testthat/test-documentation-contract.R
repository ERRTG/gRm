documentation_top_level_functions <- function(path) {
  parsed <- parse(path, keep.source = TRUE)
  references <- attr(parsed, "srcref")
  rows <- list()
  for (index in seq_along(parsed)) {
    expression <- parsed[[index]]
    if (
      is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && is.call(expression[[3L]]) &&
      identical(expression[[3L]][[1L]], as.name("function"))
    ) {
      formal_names <- names(as.list(expression[[3L]][[2L]]))
      if (is.null(formal_names)) {
        formal_names <- character()
      }
      rows[[length(rows) + 1L]] <- list(
        name = as.character(expression[[2L]]),
        line = as.integer(references[[index]][1L]),
        formals = formal_names
      )
    }
  }
  rows
}

documentation_roxygen_block <- function(lines, function_line) {
  if (function_line <= 1L || !grepl("^#'", lines[[function_line - 1L]])) {
    return(character())
  }
  first <- function_line - 1L
  while (first > 1L && grepl("^#'", lines[[first - 1L]])) {
    first <- first - 1L
  }
  lines[first:(function_line - 1L)]
}

documentation_param_names <- function(block) {
  rows <- grep("@param\\s+", block, value = TRUE)
  if (!length(rows)) {
    return(character())
  }
  names <- sub("^.*@param\\s+(`[^`]+`|[^ ,]+).*$", "\\1", rows)
  gsub("`", "", names, fixed = TRUE)
}

test_that("every top-level function has an adjacent role-appropriate roxygen block", {
  namespace <- readLines(package_path("NAMESPACE"), warn = FALSE)
  exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", namespace, value = TRUE))
  methods <- sub(
    "^S3method\\(([^,]+),([^\\)]+)\\)$",
    "\\1.\\2",
    grep("^S3method\\(", namespace, value = TRUE)
  )
  public <- unique(c(exports, methods))
  r_files <- sort(list.files(package_path("R"), pattern = "[.]R$", full.names = TRUE))

  observed <- list()
  for (path in r_files) {
    lines <- readLines(path, warn = FALSE)
    functions <- documentation_top_level_functions(path)
    for (entry in functions) {
      block <- documentation_roxygen_block(lines, entry$line)
      expect_true(
        length(block) > 0L,
        info = paste(basename(path), entry$name, "must have immediately adjacent roxygen")
      )
      is_public <- entry$name %in% public
      if (!is_public) {
        expect_true(
          any(grepl("@keywords\\s+internal(?:\\s|$)", block, perl = TRUE)),
          info = paste(basename(path), entry$name, "must be marked internal")
        )
        expect_true(
          any(grepl("@return(?:\\s|$)", block, perl = TRUE)),
          info = paste(basename(path), entry$name, "must document its return value")
        )
        expect_setequal(
          documentation_param_names(block),
          entry$formals
        )
      }
      observed[[length(observed) + 1L]] <- data.frame(
        r_function = entry$name,
        r_file = file.path("R", basename(path)),
        line = entry$line,
        exported = entry$name %in% exports,
        s3_method = entry$name %in% methods,
        roxygen_attached = length(block) > 0L,
        keywords_internal = any(grepl("@keywords\\s+internal", block)),
        formals = paste(entry$formals, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
  observed <- do.call(rbind, observed)
  inventory <- utils::read.csv(
    package_path("tools", "r-function-inventory.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  key <- c("r_function", "r_file", "line")
  inventory_key <- inventory[order(inventory$r_file, inventory$line), key]
  observed_key <- observed[order(observed$r_file, observed$line), key]
  rownames(inventory_key) <- NULL
  rownames(observed_key) <- NULL
  expect_equal(
    inventory_key,
    observed_key
  )
  expect_true(all(inventory$roxygen_attached))
  expect_true(all(inventory$keywords_internal[!inventory$exported & !inventory$s3_method]))
})

test_that("the numerical source-trace manifest is complete and locally anchored", {
  inventory <- utils::read.csv(
    package_path("tools", "r-function-inventory.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  manifest <- utils::read.csv(
    package_path("tools", "numerical-source-trace.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  numerical <- inventory[inventory$numerical, c("r_function", "r_file")]
  manifest_key <- data.frame(
    r_function = manifest$r_function,
    r_file = sub("^gRm/", "", manifest$r_file),
    stringsAsFactors = FALSE
  )
  numerical <- numerical[order(numerical$r_file, numerical$r_function), ]
  manifest_key <- manifest_key[order(manifest_key$r_file, manifest_key$r_function), ]
  rownames(numerical) <- NULL
  rownames(manifest_key) <- NULL
  expect_equal(numerical, manifest_key)
  expect_true(all(nzchar(manifest$original_pascal_file)))
  expect_true(all(nzchar(manifest$pascal_routine)))
  expect_true(all(nzchar(manifest$algorithms_section)))
  expect_true(all(nzchar(manifest$unit_tests)))
  expect_true(all(nzchar(manifest$oracle_targets)))

  for (index in seq_len(nrow(manifest))) {
    row <- manifest[index, , drop = FALSE]
    r_path <- package_path(sub("^gRm/R/", "R/", row$r_file))
    lines <- readLines(r_path, warn = FALSE)
    functions <- documentation_top_level_functions(r_path)
    hit <- which(vapply(functions, `[[`, character(1L), "name") == row$r_function)
    expect_length(hit, 1L)
    block <- documentation_roxygen_block(lines, functions[[hit]]$line)
    expect_true(
      any(grepl(
        paste0(row$original_pascal_file, "::", row$pascal_routine),
        block,
        fixed = TRUE
      )),
      info = paste(row$r_file, row$r_function, "must carry its exact local trace")
    )
  }

  # Original Pascal, harness, algorithm, and oracle paths belong to the full
  # development repository and are deliberately absent from an installed
  # package. Validate them whenever that repository surface is available.
  if (!dir.exists(repo_path("source"))) {
    skip("Full development source tree is not present in this package check.")
  }
  algorithms <- readLines(repo_path("gRm", "ALGORITHMS.md"), warn = FALSE)
  source_cache <- new.env(parent = emptyenv())
  source_text <- function(path) {
    if (!exists(path, envir = source_cache, inherits = FALSE)) {
      assign(
        path,
        tolower(paste(readLines(path, warn = FALSE, encoding = "latin1"), collapse = "\n")),
        envir = source_cache
      )
    }
    get(path, envir = source_cache, inherits = FALSE)
  }
  for (index in seq_len(nrow(manifest))) {
    row <- manifest[index, , drop = FALSE]
    original <- repo_path(row$original_pascal_file)
    expect_true(file.exists(original), info = row$original_pascal_file)
    expect_true(
      grepl(tolower(row$pascal_routine), source_text(original), fixed = TRUE),
      info = paste(row$original_pascal_file, row$pascal_routine)
    )
    expect_true(file.exists(repo_path(row$unit_tests)), info = row$unit_tests)
    expect_true(file.exists(repo_path(row$oracle_targets)), info = row$oracle_targets)
    expect_true(
      any(grepl(paste0("#", "#", "? ", row$algorithms_section), algorithms)),
      info = paste("Missing ALGORITHMS section", row$algorithms_section)
    )
    if (nzchar(row$harness_file)) {
      harness <- repo_path(row$harness_file)
      expect_true(file.exists(harness), info = row$harness_file)
      expect_true(
        grepl(tolower(row$harness_routine), source_text(harness), fixed = TRUE),
        info = paste(row$harness_file, row$harness_routine)
      )
    }
  }
})

test_that("native numerical kernels carry exact local original-source paths", {
  required <- list(
    "item_parameters_extended.cpp" = c(
      "source/GLLRM.txt::write_iteminformation1",
      "source/PAS_skunits/skbias12.pas::CalculateICEandMICE"
    ),
    "gllrm_expected.cpp" = c(
      "source/PAS_skunits/skbias12b.pas::Initialize_GLLRMinfo",
      "source/PAS_skunits/skbias22.pas::Gamma_calculation",
      "source/GLLRM_ESTIM.txt::CalculateBiasedGammaValues2"
    ),
    "screen_j_exact.cpp" = c(
      "source/PAS_skunits/SKrandom.pas::GENTAB1",
      "source/PAS_skunits/SKbigtab.pas::Transfer_BT_to_XYZ_TABLE",
      "source/PAS_skunits/SkStat.pas::RCCHI",
      "source/PAS_skunits/SKbias3.pas::XYZ_bias_ANALYSE"
    )
  )
  for (file in names(required)) {
    text <- paste(readLines(package_path("src", file), warn = FALSE), collapse = "\n")
    for (trace in required[[file]]) {
      expect_match(text, trace, fixed = TRUE, info = paste(file, trace))
    }
  }
})

test_that("the active algorithm inventory matches files, exports, and S3 registrations", {
  algorithms_path <- repo_path("gRm", "ALGORITHMS.md")
  if (!file.exists(algorithms_path)) {
    skip("Contributor algorithm specification is not installed with the package.")
  }
  algorithms <- readLines(algorithms_path, warn = FALSE)
  source_trace_end <- match("Primary source traces:", algorithms)
  expect_false(is.na(source_trace_end))
  inventory <- algorithms[seq_len(source_trace_end - 1L)]

  documented_r <- sub(
    "^- `([^`]+)`.*$", "\\1",
    grep("^- `R/[^`]+[.]R`", inventory, value = TRUE)
  )
  active_r <- file.path(
    "R",
    sort(list.files(repo_path("gRm", "R"), pattern = "[.]R$"))
  )
  expect_setequal(documented_r, active_r)

  documented_native <- sub(
    "^- `([^`]+)`.*$", "\\1",
    grep("^- `src/[^`]+[.](?:c|cpp)`", inventory, value = TRUE, perl = TRUE)
  )
  active_native <- file.path(
    "src",
    sort(list.files(repo_path("gRm", "src"), pattern = "[.](?:c|cpp)$"))
  )
  expect_setequal(documented_native, active_native)

  between_markers <- function(start, end) {
    first <- match(start, algorithms)
    last <- match(end, algorithms)
    expect_false(is.na(first))
    expect_false(is.na(last))
    expect_gt(last, first)
    algorithms[seq.int(first + 1L, last - 1L)]
  }
  signature_lines <- between_markers(
    "<!-- BEGIN EXPORTED SIGNATURES -->",
    "<!-- END EXPORTED SIGNATURES -->"
  )
  documented_exports <- sub(
    "^([[:alnum:]_.]+)\\(.*$", "\\1",
    grep("^[[:alnum:]_.]+\\(", signature_lines, value = TRUE)
  )
  namespace <- readLines(package_path("NAMESPACE"), warn = FALSE)
  exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", namespace, value = TRUE))
  expect_setequal(documented_exports, exports)

  method_lines <- between_markers(
    "<!-- BEGIN REGISTERED S3 METHODS -->",
    "<!-- END REGISTERED S3 METHODS -->"
  )
  documented_methods <- grep(
    "^[[:alnum:]_.]+[.][[:alnum:]_.]+$",
    method_lines,
    value = TRUE
  )
  methods <- sub(
    "^S3method\\(([^,]+),([^\\)]+)\\)$",
    "\\1.\\2",
    grep("^S3method\\(", namespace, value = TRUE)
  )
  expect_setequal(documented_methods, methods)
})
