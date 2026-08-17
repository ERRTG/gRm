gllrm_native_boundary_internal <- function(name) {
  get(name, envir = asNamespace("gRm"), mode = "function")
}

gllrm_native_boundary_fixture <- function() {
  data <- data.frame(
    ID = 1:8,
    I1 = rep(1:2, each = 4L),
    I2 = rep(rep(1:2, each = 2L), 2L),
    I3 = rep(1:2, 4L),
    X1 = rep(1:2, 4L)
  )
  analysis <- gRm(
    data,
    items = c("I1", "I2", "I3"),
    exogenous = "X1",
    id = "ID",
    score_cuts = c(1L, 3L)
  )
  model <- gllrm(
    analysis,
    ld = ~ I1:I2 + I2:I3,
    dif = ~ I1:X1 + I3:X1
  )
  bundle <- gllrm_native_boundary_internal("build_item_parameters_bundle")(
    analysis$project
  )
  context <- gllrm_native_boundary_internal("build_gllrm_context")(model, bundle)
  state <- gllrm_native_boundary_internal("initialize_gllrm_state")(context)
  list(
    input = gllrm_native_boundary_internal("gllrm_expected_native_input")(context),
    item_gamma = state$item_gamma,
    ld = state$ld_parameters,
    dif = state$dif_parameters
  )
}

gllrm_native_boundary_call <- function(fixture) {
  .Call(
    "gRm_gllrm_expected_margins",
    fixture$input,
    fixture$item_gamma,
    fixture$ld,
    fixture$dif,
    PACKAGE = "gRm"
  )
}

gllrm_native_boundary_mutate <- function(fixture, mutate) {
  copy <- unserialize(serialize(fixture, NULL, version = 2L))
  mutate(copy)
}

test_that("native GLLRM boundary accepts the complete package-constructed graph", {
  skip_if_not(
    is.loaded("gRm_gllrm_expected_margins", PACKAGE = "gRm"),
    "native GLLRM expected-margin routine is unavailable"
  )
  result <- gllrm_native_boundary_call(gllrm_native_boundary_fixture())

  expect_named(result, c("expected_items", "expected_ld", "expected_dif"))
  expect_true(all(is.finite(result$expected_items)))
})

test_that("native GLLRM boundary rejects invalid item and configuration indices", {
  fixture <- gllrm_native_boundary_fixture()
  mutations <- list(
    component_zero = function(x) {
      x$input$components[[1L]][1L] <- 0L
      x
    },
    component_upper = function(x) {
      x$input$components[[1L]][1L] <- nrow(x$item_gamma) + 1L
      x
    },
    component_na = function(x) {
      x$input$components[[1L]][1L] <- NA_integer_
      x
    },
    component_fractional = function(x) {
      x$input$components[[1L]] <- as.numeric(x$input$components[[1L]])
      x$input$components[[1L]][1L] <- 1.5
      x
    },
    config_negative = function(x) {
      x$input$config_matrices[[1L]][1L, 1L] <- -1L
      x
    },
    config_upper = function(x) {
      item <- x$input$components[[1L]][1L]
      x$input$config_matrices[[1L]][1L, 1L] <- x$input$item_raw_max[[item]]
      x
    },
    config_na = function(x) {
      x$input$config_matrices[[1L]][1L, 1L] <- NA_integer_
      x
    },
    config_fractional = function(x) {
      storage.mode(x$input$config_matrices[[1L]]) <- "double"
      x$input$config_matrices[[1L]][1L, 1L] <- 0.5
      x
    },
    score_negative = function(x) {
      x$input$config_scores[[1L]][1L] <- -1L
      x
    },
    score_upper = function(x) {
      x$input$config_scores[[1L]][1L] <- x$input$max_total_score + 1L
      x
    },
    score_mismatch = function(x) {
      x$input$config_scores[[1L]][1L] <- x$input$config_scores[[1L]][1L] + 1L
      x
    }
  )

  for (name in names(mutations)) {
    malformed <- gllrm_native_boundary_mutate(fixture, mutations[[name]])
    expect_error(gllrm_native_boundary_call(malformed), info = name)
  }
})

test_that("native GLLRM boundary rejects invalid LD and DIF index graphs", {
  fixture <- gllrm_native_boundary_fixture()
  mutations <- list(
    ld_index_zero = function(x) {
      x$input$ld_local_matrices[[1L]][1L, 1L] <- 0L
      x
    },
    ld_index_upper = function(x) {
      x$input$ld_local_matrices[[1L]][1L, 1L] <- length(x$ld) + 1L
      x
    },
    ld_local_zero = function(x) {
      x$input$ld_local_matrices[[1L]][1L, 2L] <- 0L
      x
    },
    ld_local_upper = function(x) {
      x$input$ld_local_matrices[[1L]][1L, 2L] <-
        length(x$input$components[[1L]]) + 1L
      x
    },
    ld_local_same = function(x) {
      x$input$ld_local_matrices[[1L]][1L, 3L] <-
        x$input$ld_local_matrices[[1L]][1L, 2L]
      x
    },
    ld_na = function(x) {
      x$input$ld_local_matrices[[1L]][1L, 1L] <- NA_integer_
      x
    },
    ld_fractional = function(x) {
      storage.mode(x$input$ld_local_matrices[[1L]]) <- "double"
      x$input$ld_local_matrices[[1L]][1L, 1L] <- 1.5
      x
    },
    dif_background_zero = function(x) {
      x$input$dif_by_item[[1L]][1L, 1L] <- 0L
      x
    },
    dif_background_upper = function(x) {
      x$input$dif_by_item[[1L]][1L, 1L] <- ncol(x$input$group_backgrounds) + 1L
      x
    },
    dif_index_zero = function(x) {
      x$input$dif_by_item[[1L]][1L, 2L] <- 0L
      x
    },
    dif_index_upper = function(x) {
      x$input$dif_by_item[[1L]][1L, 2L] <- length(x$dif) + 1L
      x
    },
    dif_na = function(x) {
      x$input$dif_by_item[[1L]][1L, 2L] <- NA_integer_
      x
    },
    dif_fractional = function(x) {
      storage.mode(x$input$dif_by_item[[1L]]) <- "double"
      x$input$dif_by_item[[1L]][1L, 2L] <- 1.5
      x
    },
    dif_cache_missing = function(x) {
      x$input$dif_backgrounds <- integer()
      x
    },
    dif_cache_duplicate = function(x) {
      x$input$dif_backgrounds <- rep(x$input$dif_backgrounds, 2L)
      x
    }
  )

  for (name in names(mutations)) {
    malformed <- gllrm_native_boundary_mutate(fixture, mutations[[name]])
    expect_error(gllrm_native_boundary_call(malformed), info = name)
  }
})

test_that("native GLLRM boundary rejects invalid score, background, and parameter cells", {
  fixture <- gllrm_native_boundary_fixture()
  mutations <- list(
    group_score_negative = function(x) {
      x$input$group_scores[[1L]] <- -1L
      x
    },
    group_score_upper = function(x) {
      x$input$group_scores[[1L]] <- x$input$max_total_score + 1L
      x
    },
    group_score_fractional = function(x) {
      x$input$group_scores <- as.numeric(x$input$group_scores)
      x$input$group_scores[[1L]] <- 1.5
      x
    },
    group_score_na = function(x) {
      x$input$group_scores[[1L]] <- NA_integer_
      x
    },
    group_count_negative = function(x) {
      x$input$group_counts[[1L]] <- -1
      x
    },
    group_count_na = function(x) {
      x$input$group_counts[[1L]] <- NA_real_
      x
    },
    background_zero = function(x) {
      x$input$group_backgrounds[1L, 1L] <- 0L
      x
    },
    background_upper = function(x) {
      x$input$group_backgrounds[1L, 1L] <- ncol(x$dif[[1L]]) + 1L
      x
    },
    background_fractional = function(x) {
      storage.mode(x$input$group_backgrounds) <- "double"
      x$input$group_backgrounds[1L, 1L] <- 1.5
      x
    },
    background_na = function(x) {
      x$input$group_backgrounds[1L, 1L] <- NA_integer_
      x
    },
    item_gamma_negative = function(x) {
      x$item_gamma[1L, 1L] <- -1
      x
    },
    item_gamma_na = function(x) {
      x$item_gamma[1L, 1L] <- NA_real_
      x
    },
    ld_negative = function(x) {
      x$ld[[1L]][1L, 1L] <- -1
      x
    },
    ld_na = function(x) {
      x$ld[[1L]][1L, 1L] <- NA_real_
      x
    },
    dif_negative = function(x) {
      x$dif[[1L]][1L, 1L] <- -1
      x
    },
    dif_na = function(x) {
      x$dif[[1L]][1L, 1L] <- NA_real_
      x
    }
  )

  for (name in names(mutations)) {
    malformed <- gllrm_native_boundary_mutate(fixture, mutations[[name]])
    expect_error(gllrm_native_boundary_call(malformed), info = name)
  }
})

test_that("native GLLRM boundary rejects a deterministic fuzz corpus of malformed shapes", {
  fixture <- gllrm_native_boundary_fixture()
  mutations <- list(
    config_rows = function(x) {
      x$input$config_matrices[[1L]] <- x$input$config_matrices[[1L]][-1L, , drop = FALSE]
      x
    },
    config_columns = function(x) {
      x$input$config_matrices[[1L]] <- x$input$config_matrices[[1L]][, -1L, drop = FALSE]
      x
    },
    ld_columns = function(x) {
      x$input$ld_local_matrices[[1L]] <- x$input$ld_local_matrices[[1L]][, -1L, drop = FALSE]
      x
    },
    dif_columns = function(x) {
      x$input$dif_by_item[[1L]] <- x$input$dif_by_item[[1L]][, 1L, drop = FALSE]
      x
    },
    background_rows = function(x) {
      x$input$group_backgrounds <- x$input$group_backgrounds[-1L, , drop = FALSE]
      x
    },
    item_gamma_columns = function(x) {
      x$item_gamma <- x$item_gamma[, 1L, drop = FALSE]
      x
    },
    ld_parameter_rows = function(x) {
      x$ld[[1L]] <- x$ld[[1L]][-1L, , drop = FALSE]
      x
    },
    dif_parameter_columns = function(x) {
      x$dif[[1L]] <- x$dif[[1L]][, 1L, drop = FALSE]
      x
    }
  )
  set.seed(20260816)

  for (name in sample(names(mutations))) {
    malformed <- gllrm_native_boundary_mutate(fixture, mutations[[name]])
    expect_error(gllrm_native_boundary_call(malformed), info = name)
  }
})
