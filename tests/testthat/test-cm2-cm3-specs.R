cm2_cm3_fixture_context <- function(
    item_names = c("I1", "I2", "I3", "I4"),
    item_label_codes = c("A", "B", "C", "D"),
    background_names = c("X1", "X2"),
    background_label_codes = c("X", "Y"),
    ld_specs = list(),
    dif_specs = list()) {
  list(
    project = list(
      items = data.frame(
        name = item_names,
        label_code = item_label_codes,
        stringsAsFactors = FALSE
      ),
      backgrounds = data.frame(
        name = background_names,
        label_code = background_label_codes,
        stringsAsFactors = FALSE
      )
    ),
    n_items = length(item_names),
    n_backgrounds = length(background_names),
    ld_specs = ld_specs,
    dif_specs = dif_specs
  )
}

test_that("CM2/CM3 selected items follow DIGRAM UseItems source-order semantics", {
  context <- cm2_cm3_fixture_context()

  default <- cm2_cm3_selected_items(context, NULL)
  expect_equal(default$schema_version, 1L)
  expect_equal(default$mode, "default_all")
  expect_null(default$requested_items)
  expect_equal(default$resolved_items$item_index, 1:4)
  expect_equal(default$resolved_items$item_label, c("A", "B", "C", "D"))
  expect_equal(default$resolved_items$item_name, c("I1", "I2", "I3", "I4"))
  expect_equal(default$selected_count, 4L)
  expect_equal(default$model_item_count, 4L)
  expect_equal(default$exogenous_scope, "all_fitted")
  expect_equal(default$exogenous$exogenous_index, 1:2)
  expect_equal(default$exogenous$exogenous_label, c("X", "Y"))
  expect_equal(default$exogenous$exogenous_name, c("X1", "X2"))
  expect_true(default$score_group_included)
  expect_equal(default$score_total_scope, "all_fitted_items")

  selected <- cm2_cm3_selected_items(context, c("I3", "I1"))
  expect_equal(selected$mode, "explicit_names")
  expect_equal(selected$requested_items, c("I3", "I1"))
  expect_equal(selected$resolved_items$item_index, c(1L, 3L))

  indexed <- cm2_cm3_selected_items(context, c(4L, 2L))
  expect_equal(indexed$mode, "explicit_indices")
  expect_equal(indexed$requested_items, c(4L, 2L))
  expect_equal(indexed$resolved_items$item_index, c(2L, 4L))

  expect_error(
    cm2_cm3_selected_items(context, c("I1", "I1")),
    "duplicate",
    ignore.case = TRUE
  )
  expect_error(
    cm2_cm3_selected_items(context, "unknown"),
    "unknown",
    ignore.case = TRUE
  )
  expect_error(
    cm2_cm3_selected_items(context, "A"),
    "unknown",
    ignore.case = TRUE
  )
  expect_error(
    cm2_cm3_selected_items(context, c(1L, 1L)),
    "duplicate",
    ignore.case = TRUE
  )
  expect_error(
    cm2_cm3_selected_items(context, c(1.5, 2)),
    "whole|integer",
    ignore.case = TRUE
  )
  expect_error(
    cm2_cm3_selected_items(context, c(0L, 2L)),
    "range|positive",
    ignore.case = TRUE
  )
  expect_error(
    cm2_cm3_selected_items(context, "I1"),
    "at least two",
    ignore.case = TRUE
  )
})

test_that("CM2/CM3 item selectors reject every ambiguous or invalid shape", {
  context <- cm2_cm3_fixture_context()

  invalid <- list(
    character(), numeric(), "", "  ",
    c("I1", NA_character_),
    TRUE, c(TRUE, FALSE), factor(c("I1", "I2")), list("I1", "I2"),
    matrix(c("I1", "I2"), nrow = 1L), array(c(1L, 2L), dim = c(2L, 1L, 1L)),
    c(1, NA_real_), c(1, NaN), c(1, Inf), c(1, -Inf), c(1, 2.5),
    c(0L, 2L), c(-1L, 2L), c(1L, 5L), c(1, 3e9),
    c("I1", "X1"), c("I1", "Score group"), c("I1", "#")
  )
  for (selector in invalid) {
    expect_error(cm2_cm3_selected_items(context, selector))
  }
})

test_that("the string all is an ordinary exact fitted item name", {
  context <- cm2_cm3_fixture_context(
    item_names = c("all", "other"),
    item_label_codes = c("A", "B")
  )

  selected <- cm2_cm3_selected_items(context, c("other", "all"))
  expect_equal(selected$mode, "explicit_names")
  expect_equal(selected$resolved_items$item_name, c("all", "other"))
  expect_error(cm2_cm3_selected_items(context, "all"), "at least two", ignore.case = TRUE)
})

test_that("CM2 margin specs are prepared in source order with source LD/DIF skips", {
  context <- cm2_cm3_fixture_context(
    ld_specs = list(list(item1 = 1L, item2 = 2L)),
    dif_specs = list(list(item = 2L, background = 1L))
  )

  specs <- cm2_cm3_prepare_margins(context, 1:3, include_three_way = FALSE)

  expect_equal(
    vapply(specs, `[[`, character(1L), "kind"),
    c(
      "item_item", "item_item",
      "item_exogenous", "item_exogenous", "item_score_group",
      "item_exogenous", "item_score_group",
      "item_exogenous", "item_exogenous", "item_score_group"
    )
  )
  expect_equal(
    vapply(specs, function(spec) cm2_cm3_margin_name(context, spec), character(1L)),
    c(
      "I1:I3", "I2:I3",
      "I1:X1", "I1:X2", "I1:Score group",
      "I2:X2", "I2:Score group",
      "I3:X1", "I3:X2", "I3:Score group"
    )
  )
})

test_that("CM2 item-exogenous source skip follows corrected DIGRAM 7.04 ItemBias order", {
  context <- cm2_cm3_fixture_context(
    item_names = sprintf("Item%02d", 1:5),
    item_label_codes = letters[1:5],
    background_names = c("Exo3", "Group2", "Site4"),
    background_label_codes = letters[6:8],
    ld_specs = list(list(item1 = 2L, item2 = 3L)),
    dif_specs = list(
      list(item = 2L, background = 1L),
      list(item = 3L, background = 1L),
      list(item = 4L, background = 2L)
    )
  )

  specs <- cm2_cm3_prepare_margins(context, 1:5, include_three_way = FALSE)
  names <- vapply(specs, function(spec) cm2_cm3_margin_name(context, spec), character(1L))

  expect_length(specs, 26L)
  expect_true("Item01:Group2" %in% names)
  expect_true("Item01:Site4" %in% names)
  expect_false("Item02:Exo3" %in% names)
  expect_false("Item03:Exo3" %in% names)
  expect_false("Item04:Group2" %in% names)
})

test_that("CM3 margin specs append all three-way rows after CM2 rows", {
  context <- cm2_cm3_fixture_context(
    ld_specs = list(list(item1 = 1L, item2 = 2L)),
    dif_specs = list(list(item = 2L, background = 1L))
  )

  specs <- cm2_cm3_prepare_margins(context, 1:3, include_three_way = TRUE)
  names <- vapply(specs, function(spec) cm2_cm3_margin_name(context, spec), character(1L))
  kinds <- vapply(specs, `[[`, character(1L), "kind")

  expect_equal(length(specs), 29L)
  expect_equal(kinds[1:10], vapply(
    cm2_cm3_prepare_margins(context, 1:3, include_three_way = FALSE),
    `[[`,
    character(1L),
    "kind"
  ))
  expect_true(all(kinds[seq_len(10L)] != "item_item_item"))
  expect_equal(names[11L], "I1:I2:I3")
  expect_equal(
    names[12:20],
    c(
      "I1:I2:X1", "I1:I2:X2", "I1:I2:Score group",
      "I1:I3:X1", "I1:I3:X2", "I1:I3:Score group",
      "I2:I3:X1", "I2:I3:X2", "I2:I3:Score group"
    )
  )
  expect_equal(
    names[21:29],
    c(
      "I1:X1:X2", "I1:X1:Score group", "I1:X2:Score group",
      "I2:X1:X2", "I2:X1:Score group", "I2:X2:Score group",
      "I3:X1:X2", "I3:X1:Score group", "I3:X2:Score group"
    )
  )

  expect_true("I1:I2:X1" %in% names)
  expect_true("I3:X1:X2" %in% names)
})

test_that("CM2/CM3 no-exogenous combinatorics retain the automatic score group", {
  context <- cm2_cm3_fixture_context(
    background_names = character(),
    background_label_codes = character()
  )

  cm2_specs <- cm2_cm3_prepare_margins(context, 1:3, include_three_way = FALSE)
  cm3_specs <- cm2_cm3_prepare_margins(context, 1:3, include_three_way = TRUE)
  expect_length(cm2_specs, 6L)
  expect_length(cm3_specs, 10L)
  expect_equal(
    vapply(cm2_specs, `[[`, character(1L), "kind"),
    c(rep("item_item", 3L), rep("item_score_group", 3L))
  )
  expect_equal(
    vapply(cm3_specs[7:10], `[[`, character(1L), "kind"),
    c("item_item_item", rep("item_item_score_group", 3L))
  )
})

test_that("fitted LD and DIF suppress only eligible selected CM2 rows", {
  context <- cm2_cm3_fixture_context(
    ld_specs = list(list(item1 = 1L, item2 = 2L)),
    dif_specs = list(list(item = 1L, background = 1L))
  )
  names_for <- function(selected, three_way = FALSE) {
    specs <- cm2_cm3_prepare_margins(context, selected, include_three_way = three_way)
    vapply(specs, function(spec) cm2_cm3_margin_name(context, spec), character(1L))
  }

  both <- names_for(c(1L, 2L))
  one <- names_for(c(1L, 3L))
  neither <- names_for(c(3L, 4L))
  expect_false("I1:I2" %in% both)
  expect_false("I1:X1" %in% both)
  expect_true("I1:I3" %in% one)
  expect_false("I1:X1" %in% one)
  expect_true("I3:I4" %in% neither)
  expect_true("I3:X1" %in% neither)
  expect_false(any(grepl("^I1", neither)))

  three_way <- names_for(c(1L, 2L), three_way = TRUE)
  expect_true("I1:I2:X1" %in% three_way)
  expect_true("I1:I2:Score group" %in% three_way)
})

test_that("CM2/CM3 lookup and margin naming helpers separate labels from public names", {
  context <- cm2_cm3_fixture_context(
    ld_specs = list(list(item1 = 2L, item2 = 1L)),
    dif_specs = list(list(item = 2L, background = 1L))
  )
  specs <- cm2_cm3_prepare_margins(context, 1:3, include_three_way = TRUE)
  included_ld <- cm2_cm3_included_ld_lookup(context)
  included_dif <- cm2_cm3_included_dif_lookup(context)

  expect_true(included_ld[1L, 2L])
  expect_true(included_ld[2L, 1L])
  expect_false(included_dif[1L, 2L])
  expect_true(included_dif[2L, 1L])

  margin_names <- vapply(specs, function(spec) cm2_cm3_margin_name(context, spec), character(1L))
  item_score <- specs[[which(margin_names == "I1:Score group")]]
  three_way <- specs[[which(margin_names == "I1:I2:X1")]]

  expect_equal(cm2_cm3_margin_label(context, item_score), "A:Score group")
  expect_equal(cm2_cm3_margin_name(context, item_score), "I1:Score group")
  expect_equal(cm2_cm3_margin_public_variables(context, item_score), c("I1", "Score group"))

  expect_equal(cm2_cm3_margin_label(context, three_way), "A:B:X")
  expect_equal(cm2_cm3_margin_name(context, three_way), "I1:I2:X1")
  expect_equal(cm2_cm3_margin_public_variables(context, three_way), c("I1", "I2", "X1"))
})
