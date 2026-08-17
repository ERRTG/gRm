extended_native_internal <- function(name) {
  get(name, envir = asNamespace("gRm"), mode = "function")
}

test_that("native Extended capability is explicit and internally consistent", {
  capability <- extended_native_internal("source_extended_native_capability")()

  expect_identical(capability$schema, "gRm-pascal-extended-capability-v1")
  expect_type(capability$mantissa_bits, "integer")
  expect_type(capability$max_exponent, "integer")
  expect_type(capability$storage_bytes, "integer")
  expect_type(capability$pascal_extended_fixed_field_supported, "logical")
  expect_length(capability$pascal_extended_fixed_field_supported, 1L)
  expect_identical(
    capability$pascal_extended_fixed_field_supported,
    isTRUE(capability$x87_available) &&
      capability$mantissa_bits >= 64L && capability$max_exponent >= 16384L
  )
})

test_that("top-category ICE cancellation derives positive, negative, and zero signs", {
  skip_if_not(
    is.loaded("gRm_item_parameters_top_ice_field", PACKAGE = "gRm"),
    "native item-parameter diagnostic helper is not loaded"
  )
  capability <- extended_native_internal("source_extended_native_capability")()
  skip_if_not(
    isTRUE(capability$pascal_extended_fixed_field_supported),
    "this platform does not provide the audited Pascal Extended arithmetic path"
  )
  field <- extended_native_internal("source_extended_top_ice_cancellation_field")

  negative <- field(11 / 12345, 3L)
  positive <- field(3 / 7, 3L)
  exact_zero <- field(1, 3L)

  expect_identical(unname(as.character(negative)), "   -0.000")
  expect_identical(attr(negative, "cancellation_sign"), "negative")
  expect_lt(attr(negative, "cancellation_value"), 0)

  expect_identical(unname(as.character(positive)), "    0.000")
  expect_identical(attr(positive, "cancellation_sign"), "positive")
  expect_gt(attr(positive, "cancellation_value"), 0)

  expect_identical(unname(as.character(exact_zero)), "    0.000")
  expect_identical(attr(exact_zero, "cancellation_sign"), "zero")
  expect_identical(attr(exact_zero, "cancellation_value"), 0)
})

test_that("source ICE fixed fields cover both sides of the three-decimal boundary", {
  skip_if_not(
    is.loaded("gRm_item_parameters_ice_fields_from_gamma", PACKAGE = "gRm"),
    "native item-parameter formatter is not loaded"
  )
  capability <- extended_native_internal("source_extended_native_capability")()
  skip_if_not(
    isTRUE(capability$pascal_extended_fixed_field_supported),
    "exact Pascal Extended field boundaries require the audited x87 path"
  )
  fields <- extended_native_internal("source_item_parameters_extended_ice_fields")
  bundle <- list(model = list(items = data.frame(raw_max = rep(3L, 4L))))
  fit <- list(item_gamma = rbind(
    below_positive = c(1, exp(0.00049), 1),
    above_positive = c(1, exp(0.00051), 1),
    below_negative = c(1, exp(-0.00049), 1),
    above_negative = c(1, exp(-0.00051), 1)
  ))

  actual <- fields(fit, bundle)[, 2L]
  expect_identical(
    unname(actual),
    c("    0.000", "    0.001", "   -0.000", "   -0.001")
  )
})

test_that("unused native Extended replay solvers are absent", {
  expect_false(is.loaded("gRm_item_parameters_extended_top_ice_fields", PACKAGE = "gRm"))
  expect_false(is.loaded("gRm_item_parameters_extended_gamma", PACKAGE = "gRm"))

  namespace <- asNamespace("gRm")
  expect_false(exists("source_item_parameters_top_ice_fields", namespace, inherits = FALSE))
  expect_false(exists("source_item_parameters_extended_gamma", namespace, inherits = FALSE))
})
