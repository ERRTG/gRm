screen_j_ld_fixture <- function() {
  item_names <- paste0("I", 1:4)
  gamma <- matrix(0, nrow = 4L, ncol = 4L, dimnames = list(item_names, item_names))
  p_value <- matrix(2, nrow = 4L, ncol = 4L, dimnames = list(item_names, item_names))
  wpg <- matrix(0, nrow = 4L, ncol = 4L, dimnames = list(item_names, item_names))

  gamma[1L, 2L] <- 0.40
  gamma[2L, 1L] <- 0.30
  p_value[1L, 2L] <- 0.001
  p_value[2L, 1L] <- 0.002
  wpg[1L, 2L] <- wpg[2L, 1L] <- 0.35

  gamma[3L, 4L] <- -0.30
  gamma[4L, 3L] <- -0.20
  p_value[3L, 4L] <- 0.003
  p_value[4L, 3L] <- 0.004
  wpg[3L, 4L] <- wpg[4L, 3L] <- -0.25

  list(
    names = item_names,
    gamma = gamma,
    p_value = p_value,
    wpg = wpg,
    selection = screen_j_stepwise_local_dependence(
      gamma,
      p_value,
      wpg,
      fdr05 = 0.01
    )
  )
}

test_that("SCREEN J separates provisional LD evidence from final model terms", {
  fixture <- screen_j_ld_fixture()
  ld <- fixture$selection

  expect_true(ld$stepwise_matrix[1L, 2L])
  expect_true(ld$stepwise_matrix[3L, 4L])
  expect_equal(ld$stepwise_rows$stage, c("positive2", "negative2"))
  expect_equal(ld$stepwise_rows$directed_gamma_sum, c(0.7, -0.5))
  expect_equal(ld$stepwise_rows$included, c(TRUE, FALSE))

  expect_true(ld$matrix[1L, 2L])
  expect_false(ld$matrix[3L, 4L])
  expect_equal(cbind(ld$rows$row, ld$rows$col), matrix(c(1L, 2L), nrow = 1L))
  expect_equal(
    cbind(ld$negative_rows$row, ld$negative_rows$col),
    matrix(c(3L, 4L), nrow = 1L)
  )
  expect_true(ld$negative_matrix[3L, 4L])
  expect_false(ld$negative_matrix[1L, 2L])
})

test_that("the final LD rule uses the directed-gamma sum rather than WPG or stage", {
  gamma <- matrix(0, nrow = 2L, ncol = 2L)
  gamma[1L, 2L] <- 0.20
  gamma[2L, 1L] <- -0.30
  provisional <- matrix(c(FALSE, TRUE, TRUE, FALSE), nrow = 2L)
  rows <- data.frame(
    row = 1L,
    col = 2L,
    pair = "I1I2",
    value = 0.80,
    stage = "positive1",
    stringsAsFactors = FALSE
  )

  ld <- screen_j_finalize_local_dependence(gamma, provisional, rows)

  expect_equal(ld$stepwise_rows$directed_gamma_sum, -0.10)
  expect_false(ld$stepwise_rows$included)
  expect_false(any(ld$matrix))
  expect_true(ld$negative_matrix[1L, 2L])
})

test_that("screen summaries and screen-derived models use final LD membership", {
  fixture <- screen_j_ld_fixture()
  analysis <- gRm(
    data.frame(
      ID = 1:8,
      I1 = c(0L, 0L, 1L, 1L, 0L, 1L, 0L, 1L),
      I2 = c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L),
      I3 = c(0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L),
      I4 = c(1L, 0L, 1L, 0L, 1L, 0L, 0L, 1L)
    ),
    items = fixture$names,
    id = "ID"
  )
  values <- list(
    items = data.frame(name = fixture$names, stringsAsFactors = FALSE),
    backgrounds = data.frame(name = character(), stringsAsFactors = FALSE),
    partial = list(
      item_gamma = fixture$gamma,
      item_p = fixture$p_value,
      weighted_gamma = fixture$wpg,
      exo_p = matrix(numeric(), nrow = 4L, ncol = 0L),
      exo_stat = matrix(numeric(), nrow = 4L, ncol = 0L),
      exo_kind = character()
    ),
    model = list(
      local_dependence = fixture$selection,
      item_bias = matrix(FALSE, nrow = 4L, ncol = 0L),
      score_effects = list(rows = data.frame())
    ),
    bh = list(fdr_05 = 0.01, fdr_01 = 0.002, fdr_001 = 0.001, n_tests = 12L)
  )
  screen_object <- structure(
    list(
      analysis = analysis,
      values = values,
      inference = "asymptotic",
      exact_state = list(exact = FALSE),
      source_trace = character(),
      terms = NULL
    ),
    class = c("gRm_screen", "list")
  )
  screen_object$terms <- model_terms(screen_object)

  screen_summary <- summary(screen_object)
  ld_table <- screen_summary$local_dependence
  positive_row <- ld_table$`Item 1` == "I1" & ld_table$`Item 2` == "I2"
  negative_row <- ld_table$`Item 1` == "I3" & ld_table$`Item 2` == "I4"

  expect_equal(screen_summary$selected$type, "ld")
  expect_equal(screen_summary$selected$item1, "I1")
  expect_equal(screen_summary$selected$item2, "I2")
  expect_equal(screen_summary$selected$status, "selected")
  expect_equal(
    names(screen_summary$tables),
    c("local_dependence", "dif", "score_effects")
  )

  expect_equal(ld_table$`Gamma 1->2`[positive_row], 0.40)
  expect_equal(ld_table$`Pr(>|Gamma 1->2|)`[positive_row], 0.001)
  expect_equal(ld_table$`Gamma 2->1`[positive_row], 0.30)
  expect_equal(ld_table$`Pr(>|Gamma 2->1|)`[positive_row], 0.002)
  expect_equal(ld_table$`Gamma sum`[positive_row], 0.70)
  expect_equal(ld_table$Decision[positive_row], "included")
  expect_equal(ld_table$` `[positive_row], "*")
  expect_equal(ld_table$Decision[negative_row], "negative LD; not included")
  expect_equal(ld_table$` `[negative_row], "")

  expect_output(print(screen_object), "Selected LD terms: 1", fixed = TRUE)
  expect_output(print(screen_summary), "Selected model terms", fixed = TRUE)
  expect_output(print(screen_summary), "I1 -- I2", fixed = TRUE)
  expect_output(print(screen_summary), "negative LD; not included", fixed = TRUE)

  model <- gllrm(screen_object)
  expect_equal(nrow(model$ld), 1L)
  expect_equal(model$ld$item1, "I1")
  expect_equal(model$ld$item2, "I2")
  expect_false(any(c(model$ld$item1, model$ld$item2) %in% c("I3", "I4")))
})
