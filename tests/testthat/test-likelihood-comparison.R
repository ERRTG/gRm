glm_like_gRm_fit <- function(model, name) {
  ll <- logLik(model)
  structure(
    list(
      values = list(
        log_likelihood = -as.numeric(ll),
        n_parameters = attr(ll, "df"),
        likelihood_n = attr(ll, "nobs")
      ),
      convergence = list(converged = TRUE),
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
