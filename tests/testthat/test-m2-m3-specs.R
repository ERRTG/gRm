m2_m3_fixture_context <- function(
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

test_that("M2/M3 selected items follow DIGRAM UseItems source-order semantics", {
  context <- m2_m3_fixture_context()

  expect_equal(m2_m3_selected_items(context, NULL), 1:4)

  selected <- m2_m3_selected_items(context, c("I3", "I1"))
  expect_equal(selected, c(1L, 3L))

  expect_equal(m2_m3_selected_items(context, c(4L, 2L)), c(2L, 4L))

  expect_error(
    m2_m3_selected_items(context, c("I1", "I1")),
    "duplicate",
    ignore.case = TRUE
  )
  expect_error(
    m2_m3_selected_items(context, "unknown"),
    "unknown",
    ignore.case = TRUE
  )
  expect_error(
    m2_m3_selected_items(context, "A"),
    "unknown",
    ignore.case = TRUE
  )
  expect_error(
    m2_m3_selected_items(context, c(1L, 1L)),
    "duplicate",
    ignore.case = TRUE
  )
  expect_error(
    m2_m3_selected_items(context, c(1.5, 2)),
    "whole|integer",
    ignore.case = TRUE
  )
  expect_error(
    m2_m3_selected_items(context, c(0L, 2L)),
    "range|positive",
    ignore.case = TRUE
  )
  expect_error(
    m2_m3_selected_items(context, "I1"),
    "at least two",
    ignore.case = TRUE
  )
})

test_that("M2 margin specs are prepared in source order with source LD/DIF skips", {
  context <- m2_m3_fixture_context(
    ld_specs = list(list(item1 = 1L, item2 = 2L)),
    dif_specs = list(list(item = 2L, background = 1L))
  )

  specs <- m2_m3_prepare_margins(context, 1:3, include_three_way = FALSE)

  expect_equal(
    vapply(specs, `[[`, character(1L), "kind"),
    c(
      "item_item", "item_item",
      "item_exogenous", "item_score_group",
      "item_exogenous", "item_exogenous", "item_score_group",
      "item_exogenous", "item_exogenous", "item_score_group"
    )
  )
  expect_equal(
    vapply(specs, function(spec) m2_m3_margin_name(context, spec), character(1L)),
    c(
      "I1:I3", "I2:I3",
      "I1:X1", "I1:Score group",
      "I2:X1", "I2:X2", "I2:Score group",
      "I3:X1", "I3:X2", "I3:Score group"
    )
  )
})

test_that("M2 item-exogenous source skip follows DIGRAM ItemBias index order", {
  context <- m2_m3_fixture_context(
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

  specs <- m2_m3_prepare_margins(context, 1:5, include_three_way = FALSE)
  names <- vapply(specs, function(spec) m2_m3_margin_name(context, spec), character(1L))

  expect_length(specs, 27L)
  expect_false("Item01:Group2" %in% names)
  expect_false("Item01:Site4" %in% names)
  expect_true("Item02:Exo3" %in% names)
  expect_true("Item03:Exo3" %in% names)
  expect_true("Item04:Group2" %in% names)
})

test_that("M3 margin specs append all three-way rows after M2 rows", {
  context <- m2_m3_fixture_context(
    ld_specs = list(list(item1 = 1L, item2 = 2L)),
    dif_specs = list(list(item = 2L, background = 1L))
  )

  specs <- m2_m3_prepare_margins(context, 1:3, include_three_way = TRUE)
  names <- vapply(specs, function(spec) m2_m3_margin_name(context, spec), character(1L))
  kinds <- vapply(specs, `[[`, character(1L), "kind")

  expect_equal(length(specs), 29L)
  expect_equal(kinds[1:10], vapply(
    m2_m3_prepare_margins(context, 1:3, include_three_way = FALSE),
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

test_that("M2/M3 lookup and margin naming helpers separate labels from public names", {
  context <- m2_m3_fixture_context(
    ld_specs = list(list(item1 = 2L, item2 = 1L)),
    dif_specs = list(list(item = 2L, background = 1L))
  )
  specs <- m2_m3_prepare_margins(context, 1:3, include_three_way = TRUE)
  included_ld <- m2_m3_included_ld_lookup(context)
  included_dif <- m2_m3_included_dif_lookup(context)

  expect_true(included_ld[1L, 2L])
  expect_true(included_ld[2L, 1L])
  expect_true(included_dif[1L, 2L])
  expect_false(included_dif[2L, 1L])

  margin_names <- vapply(specs, function(spec) m2_m3_margin_name(context, spec), character(1L))
  item_score <- specs[[which(margin_names == "I1:Score group")]]
  three_way <- specs[[which(margin_names == "I1:I2:X1")]]

  expect_equal(m2_m3_margin_label(context, item_score), "A:Score group")
  expect_equal(m2_m3_margin_name(context, item_score), "I1:Score group")
  expect_equal(m2_m3_margin_public_variables(context, item_score), c("I1", "Score group"))

  expect_equal(m2_m3_margin_label(context, three_way), "A:B:X")
  expect_equal(m2_m3_margin_name(context, three_way), "I1:I2:X1")
  expect_equal(m2_m3_margin_public_variables(context, three_way), c("I1", "I2", "X1"))
})
