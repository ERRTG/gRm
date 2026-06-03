test_that("top-category ICE cancellation uses Extended expression order", {
  skip_if_not(
    is.loaded("gRm_item_parameters_top_ice_field", PACKAGE = "gRm"),
    "native item-parameter diagnostic helper is not loaded"
  )
  field <- source_extended_top_ice_cancellation_field(11 / 12345, 3L)

  expect_identical(field, "   -0.000")
})
