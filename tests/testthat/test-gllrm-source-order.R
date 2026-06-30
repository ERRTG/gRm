source_order_analysis <- function() {
  data <- expand.grid(
    I1 = 1:2,
    I2 = 1:2,
    I3 = 1:2,
    I4 = 1:2,
    X1 = 1:2,
    X2 = 1:2,
    X3 = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  data$ID <- seq_len(nrow(data))
  data <- data[, c("ID", "I1", "I2", "I3", "I4", "X1", "X2", "X3")]
  gRm(
    data,
    items = c("I1", "I2", "I3", "I4"),
    exogenous = c("X1", "X2", "X3"),
    id = "ID"
  )
}

source_order_term_labels <- function(terms, family) {
  if (identical(family, "ld")) {
    return(paste(terms$item1, terms$item2, sep = ":"))
  }
  paste(terms$item, terms$exogenous, sep = ":")
}

test_that("GLLRM LD formula terms are stored in DIGRAM source order", {
  analysis <- source_order_analysis()
  model <- gllrm(
    analysis,
    ld = ~ I4:I2 + I1:I3 + I3:I2 + I4:I1
  )

  expect_equal(model$ld$item1, c("I1", "I1", "I2", "I2"))
  expect_equal(model$ld$item2, c("I3", "I4", "I3", "I4"))
  expect_equal(
    source_order_term_labels(model$ld, "ld"),
    c("I1:I3", "I1:I4", "I2:I3", "I2:I4")
  )
})

test_that("GLLRM DIF formula terms are stored in DIGRAM source order", {
  analysis <- source_order_analysis()
  model <- gllrm(
    analysis,
    dif = ~ X3:I4 + I2:X1 + X2:I1 + I1:X1 + I4:X1
  )

  expect_equal(model$dif$item, c("I1", "I1", "I2", "I4", "I4"))
  expect_equal(model$dif$exogenous, c("X1", "X2", "X1", "X1", "X3"))
  expect_equal(
    source_order_term_labels(model$dif, "dif"),
    c("I1:X1", "I1:X2", "I2:X1", "I4:X1", "I4:X3")
  )
})

test_that("GLLRM reversed single terms are accepted but canonicalized", {
  analysis <- source_order_analysis()
  model <- gllrm(
    analysis,
    ld = ~ I3:I1,
    dif = ~ X2:I4
  )

  expect_equal(model$ld$item1, "I1")
  expect_equal(model$ld$item2, "I3")
  expect_equal(model$dif$item, "I4")
  expect_equal(model$dif$exogenous, "X2")
})

test_that("GLLRM duplicate terms still error after canonicalization", {
  analysis <- source_order_analysis()

  expect_error(
    gllrm(analysis, ld = ~ I1:I2 + I2:I1),
    "Duplicate LD terms are not allowed",
    fixed = TRUE
  )
  expect_error(
    gllrm(analysis, dif = ~ I1:X1 + X1:I1),
    "Duplicate DIF terms are not allowed",
    fixed = TRUE
  )
})

test_that("screen-to-model conversion uses DIGRAM model order", {
  analysis <- source_order_analysis()
  item_bias <- matrix(FALSE, nrow = 4L, ncol = 3L)
  item_bias[4L, 3L] <- TRUE
  item_bias[1L, 2L] <- TRUE
  item_bias[2L, 1L] <- TRUE
  values <- list(
    items = data.frame(name = analysis$items, stringsAsFactors = FALSE),
    backgrounds = data.frame(name = analysis$exogenous, stringsAsFactors = FALSE),
    model = list(
      local_dependence = list(
        rows = data.frame(
          row = c(4L, 1L, 2L),
          col = c(3L, 4L, 3L),
          stage = c("selected", "selected", "selected")
        )
      ),
      item_bias = item_bias,
      score_effects = list(rows = data.frame())
    )
  )
  screen_obj <- list(
    analysis = analysis,
    values = values,
    source_trace = character()
  )
  class(screen_obj) <- c("gRm_screen", "list")

  model <- gllrm(screen_obj)

  expect_equal(
    source_order_term_labels(model$ld, "ld"),
    c("I1:I4", "I2:I3", "I3:I4")
  )
  expect_equal(
    source_order_term_labels(model$dif, "dif"),
    c("I1:X2", "I2:X1", "I4:X3")
  )
})

test_that("update.gRm_model reorders replacement terms", {
  analysis <- source_order_analysis()
  model <- gllrm(analysis)

  updated <- update(
    model,
    ld = ~ I4:I2 + I1:I3,
    dif = ~ X3:I4 + I1:X2
  )

  expect_equal(source_order_term_labels(updated$ld, "ld"), c("I1:I3", "I2:I4"))
  expect_equal(source_order_term_labels(updated$dif, "dif"), c("I1:X2", "I4:X3"))
})

test_that("fit objects expose source-ordered term tables", {
  analysis <- source_order_analysis()
  model_a <- gllrm(
    analysis,
    ld = ~ I1:I2 + I2:I3,
    dif = ~ I1:X1 + I3:X2
  )
  model_b <- gllrm(
    analysis,
    ld = ~ I2:I3 + I1:I2,
    dif = ~ I3:X2 + I1:X1
  )

  expect_equal(model_a$ld, model_b$ld)
  expect_equal(model_a$dif, model_b$dif)

  fit_a <- fit(model_a, max_step = 100L, max_delta = 1e-6)
  fit_b <- fit(model_b, max_step = 100L, max_delta = 1e-6)

  expect_equal(fit_a$values$log_likelihood, fit_b$values$log_likelihood, tolerance = 1e-8)
  expect_equal(fit_a$values$n_parameters, fit_b$values$n_parameters)
  expect_equal(
    vapply(fit_a$values$ld_parameter_tables, function(x) x$spec$term, character(1L)),
    c("I1:I2", "I2:I3")
  )
  expect_equal(
    vapply(fit_a$values$dif_parameter_tables, function(x) x$spec$term, character(1L)),
    c("I1:X1", "I3:X2")
  )
  expect_equal(
    vapply(fit_a$values$ld_parameter_tables, function(x) x$spec$term, character(1L)),
    vapply(fit_b$values$ld_parameter_tables, function(x) x$spec$term, character(1L))
  )
  expect_equal(
    vapply(fit_a$values$dif_parameter_tables, function(x) x$spec$term, character(1L)),
    vapply(fit_b$values$dif_parameter_tables, function(x) x$spec$term, character(1L))
  )
})

test_that("model summary prints source-ordered terms", {
  analysis <- source_order_analysis()
  model <- gllrm(
    analysis,
    ld = ~ I4:I2 + I1:I3,
    dif = ~ X3:I4 + I1:X2
  )

  printed <- capture.output(summary(model))
  ld_i1_i3 <- grep("  I1 -- I3", printed, fixed = TRUE)
  ld_i2_i4 <- grep("  I2 -- I4", printed, fixed = TRUE)
  dif_i1_x2 <- grep("  I1 by X2", printed, fixed = TRUE)
  dif_i4_x3 <- grep("  I4 by X3", printed, fixed = TRUE)

  expect_equal(length(ld_i1_i3), 1L)
  expect_equal(length(ld_i2_i4), 1L)
  expect_true(ld_i1_i3 < ld_i2_i4)
  expect_equal(length(dif_i1_x2), 1L)
  expect_equal(length(dif_i4_x3), 1L)
  expect_true(dif_i1_x2 < dif_i4_x3)
})
