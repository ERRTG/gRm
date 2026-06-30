gllrm_expected_state_diff <- function(a, b) {
  c(
    expected_items = max(abs(a$expected_items - b$expected_items), 0, na.rm = TRUE),
    expected_ld = max(c(0, unlist(Map(function(x, y) {
      max(abs(x - y), 0, na.rm = TRUE)
    }, a$expected_ld, b$expected_ld)))),
    expected_dif = max(c(0, unlist(Map(function(x, y) {
      max(abs(x - y), 0, na.rm = TRUE)
    }, a$expected_dif, b$expected_dif))))
  )
}

expect_gllrm_expected_matches_reference <- function(context, state, tolerance = 1e-10) {
  reference <- calculate_gllrm_joint_expected_margins_r(context, state)
  native <- calculate_gllrm_joint_expected_margins_cpp(context, state)
  expect_false(is.null(native), info = "native GLLRM expected-margin backend should be registered")
  diff <- gllrm_expected_state_diff(reference, native)
  expect_lte(diff[["expected_items"]], tolerance)
  expect_lte(diff[["expected_ld"]], tolerance)
  expect_lte(diff[["expected_dif"]], tolerance)
}

with_gllrm_expected_margin_backend <- function(fun, expr) {
  ns <- asNamespace("gRm")
  name <- "calculate_gllrm_joint_expected_margins"
  old <- get(name, envir = ns, inherits = FALSE)
  unlockBinding(name, ns)
  assign(name, fun, envir = ns)
  lockBinding(name, ns)
  on.exit({
    unlockBinding(name, ns)
    assign(name, old, envir = ns)
    lockBinding(name, ns)
  }, add = TRUE)
  force(expr)
}
