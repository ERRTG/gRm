glm_like_gRm_fit <- function(model, name) {
  ll <- logLik(model)
  n <- as.integer(attr(ll, "nobs"))
  analysis_fingerprint <- paste0("test-analysis:", name)
  likelihood_sample <- list(
    schema = "gRm-likelihood-sample-v1",
    fingerprint = paste0("test-sample:", name),
    row_mask = rep(TRUE, n),
    row_indices = seq_len(n),
    n = n
  )
  structure(
    list(
      values = list(
        log_likelihood = -as.numeric(ll),
        n_parameters = attr(ll, "df"),
        likelihood_n = attr(ll, "nobs")
      ),
      convergence = list(converged = TRUE),
      analysis_fingerprint = analysis_fingerprint,
      likelihood_sample = likelihood_sample,
      model = list(
        ld = data.frame(),
        dif = data.frame(),
        analysis = list(name = name)
      ),
      spec = list(
        ld = data.frame(),
        dif = data.frame(),
        analysis = list(name = name)
      ),
      call = model$call
    ),
    class = c("gRm_fit", "list")
  )
}

test_that("anova.gRm_fit matches R glm likelihood-ratio format and p-values", {
  glm0 <- glm(am ~ 1, data = mtcars, family = binomial())
  glm1 <- glm(am ~ mpg, data = mtcars, family = binomial())
  expected <- anova(glm0, glm1, test = "Chisq")

  fit0 <- glm_like_gRm_fit(glm0, "same-analysis")
  fit1 <- glm_like_gRm_fit(glm1, "same-analysis")
  actual <- anova(fit0, fit1, test = "Chisq")

  expect_s3_class(actual, "anova")
  expect_named(actual, c("Model Df", "-2 logLik", "Df", "Chisq", "Pr(>Chisq)"))
  expect_equal(actual[["Model Df"]], c(attr(logLik(glm0), "df"), attr(logLik(glm1), "df")))
  expect_equal(actual[["-2 logLik"]], expected[["Resid. Dev"]])
  expect_equal(actual[["Df"]], expected[["Df"]])
  expect_equal(actual[["Chisq"]], expected[["Deviance"]])
  expect_equal(actual[["Pr(>Chisq)"]], expected[["Pr(>Chi)"]])
  expect_equal(row.names(actual), row.names(expected))
})

test_that("anova.gRm_fit suppresses p-values for invalid comparison order", {
  glm0 <- glm(am ~ 1, data = mtcars, family = binomial())
  glm1 <- glm(am ~ mpg, data = mtcars, family = binomial())
  fit0 <- glm_like_gRm_fit(glm0, "same-analysis")
  fit1 <- glm_like_gRm_fit(glm1, "same-analysis")

  actual <- NULL
  expect_warning(
    actual <- anova(fit1, fit0, test = "Chisq"),
    "p-values set to NA"
  )

  expect_lt(actual[["Df"]][[2L]], 0)
  expect_lt(actual[["Chisq"]][[2L]], 0)
  expect_true(is.na(actual[["Pr(>Chisq)"]][[2L]]))
})

test_that("logLik.gRm_fit exposes standard R logLik sign and attributes", {
  glm0 <- glm(am ~ 1, data = mtcars, family = binomial())
  fit0 <- glm_like_gRm_fit(glm0, "same-analysis")

  actual <- logLik(fit0)
  expected <- logLik(glm0)

  expect_s3_class(actual, "logLik")
  expect_equal(as.numeric(actual), as.numeric(expected))
  expect_equal(attr(actual, "df"), attr(expected, "df"))
  expect_equal(attr(actual, "nobs"), attr(expected, "nobs"))
})

likelihood_identity_analysis <- function(data, score_cuts = c(1L, 3L)) {
  gRm(
    data,
    items = c("I1", "I2", "I3", "I4"),
    id = "ID",
    item_levels = list(I1 = 0:1, I2 = 0:1, I3 = 0:1, I4 = 0:1),
    score_cuts = score_cuts
  )
}

likelihood_identity_data <- function() {
  data <- expand.grid(
    I1 = 0:1,
    I2 = 0:1,
    I3 = 0:1,
    I4 = 0:1,
    KEEP.OUT.ATTRS = FALSE
  )
  data$ID <- seq_len(nrow(data))
  data[c("ID", "I1", "I2", "I3", "I4")]
}

test_that("analysis identity is carried unchanged through models, fits, and results", {
  analysis <- likelihood_identity_analysis(likelihood_identity_data())
  model <- gllrm(analysis)
  fitted <- fit(model, max_step = 50L)
  result <- local_dependence(fitted, max_step = 1L)

  expect_match(analysis$analysis_fingerprint, "^gRm-analysis-identity-v1:")
  expect_identical(model$analysis_fingerprint, analysis$analysis_fingerprint)
  expect_identical(fitted$analysis_fingerprint, analysis$analysis_fingerprint)
  expect_identical(result$analysis_fingerprint, analysis$analysis_fingerprint)
  expect_identical(fitted$likelihood_sample, analysis$likelihood_sample)
  expect_identical(result$likelihood_sample, analysis$likelihood_sample)
})

test_that("anova requires exact analysis data, row order, missingness, and cuts", {
  data <- likelihood_identity_data()
  reference <- fit(gllrm(likelihood_identity_analysis(data)), max_step = 50L)

  changed <- data
  changed$I1[c(2L, 3L)] <- rev(changed$I1[c(2L, 3L)])
  different_data <- fit(gllrm(likelihood_identity_analysis(changed)), max_step = 50L)
  expect_error(
    anova(reference, different_data),
    "same analysis data and encoding"
  )

  permuted <- data[rev(seq_len(nrow(data))), , drop = FALSE]
  rownames(permuted) <- NULL
  different_order <- fit(gllrm(likelihood_identity_analysis(permuted)), max_step = 50L)
  expect_error(anova(reference, different_order), "same analysis data and encoding")

  missing <- data
  missing$I1[[2L]] <- NA_integer_
  different_missingness <- fit(gllrm(likelihood_identity_analysis(missing)), max_step = 50L)
  expect_error(anova(reference, different_missingness), "same analysis data and encoding")

  different_cuts <- fit(
    gllrm(likelihood_identity_analysis(data, score_cuts = c(2L, 4L))),
    max_step = 50L
  )
  expect_error(anova(reference, different_cuts), "same analysis data and encoding")
})

test_that("anova accepts independently fitted models from a copied analysis", {
  analysis <- likelihood_identity_analysis(likelihood_identity_data())
  analysis_copy <- analysis
  fit0 <- fit(gllrm(analysis), max_step = 50L)
  fit1 <- fit(gllrm(analysis_copy, ld = ~ I1:I2), max_step = 50L)

  comparison <- anova(fit0, fit1)
  expect_s3_class(comparison, "anova")
  expect_equal(nrow(comparison), 2L)
})

test_that("anova rejects changed likelihood masks and non-nested terms", {
  analysis <- likelihood_identity_analysis(likelihood_identity_data())
  fit0 <- fit(gllrm(analysis), max_step = 50L)
  fit12 <- fit(gllrm(analysis, ld = ~ I1:I2), max_step = 50L)
  fit34 <- fit(gllrm(analysis, ld = ~ I3:I4), max_step = 50L)

  changed_sample <- fit12
  changed_sample$likelihood_sample$row_mask[c(1L, 2L)] <-
    rev(changed_sample$likelihood_sample$row_mask[c(1L, 2L)])
  changed_sample$likelihood_sample$row_indices <-
    which(changed_sample$likelihood_sample$row_mask)
  expect_error(anova(fit0, changed_sample), "same likelihood rows")

  expect_error(anova(fit12, fit34), "non-nested")
})

test_that("likelihood boundary rejects invalid parameter counts", {
  glm0 <- glm(am ~ 1, data = mtcars, family = binomial())
  fitted <- glm_like_gRm_fit(glm0, "same-analysis")

  fitted$values$n_parameters <- -1
  expect_error(logLik(fitted), "non-negative integer-like parameter count")

  fitted$values$n_parameters <- 1.5
  expect_error(logLik(fitted), "non-negative integer-like parameter count")
})
