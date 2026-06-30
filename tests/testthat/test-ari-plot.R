ari_plot_fixture <- function() {
  out <- data.frame(
    ItemNo = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
    Item = c(rep("I1", 4L), rep("I2", 4L)),
    Score = c(3L, 1L, 2L, 4L, 3L, 1L, 2L, 4L),
    n = c(4, 3, 2, 1, 4, 3, 2, 1),
    Obs0 = 0,
    Obs1 = 0,
    ObsMean = c(0.30, 0.10, 0.20, 0.40, 0.60, 0.20, 0.40, 0.80),
    ObsVar = 0,
    Exp0 = 0,
    Exp1 = 0,
    ExpMean = c(0.25, 0.15, 0.25, 0.35, 0.55, 0.25, 0.45, 0.75),
    ExpVar = c(0.09, 0.04, 0.16, 0.25, 0.36, 0.09, 0.25, 0.49),
    z = 0,
    stringsAsFactors = FALSE
  )
  class(out) <- c("gRm_ari", "data.frame")
  out
}

ari_score_intervals_internal <- function(...) {
  getFromNamespace("ari_score_intervals", "gRm")(...)
}

ari_plot_data_internal <- function(...) {
  getFromNamespace("ari_plot_data", "gRm")(...)
}

test_that("ari_score_intervals follows the SAS class-size algorithm in score order", {
  intervals <- ari_score_intervals_internal(ari_plot_fixture(), class_size = 5L)

  expect_equal(
    intervals,
    data.frame(
      Score = c(1L, 2L, 3L, 4L),
      interval = c(1L, 1L, 2L, 2L),
      stringsAsFactors = FALSE
    )
  )
})

test_that("ari_score_intervals merges a final short interval into the previous one", {
  x <- ari_plot_fixture()
  x$n[x$ItemNo == 1L] <- c(1, 5, 4, 1)

  intervals <- ari_score_intervals_internal(x, class_size = 5L)

  expect_equal(intervals$Score, c(1L, 2L, 3L, 4L))
  expect_equal(intervals$interval, c(1L, 2L, 2L, 2L))
})

test_that("ari_score_intervals keeps a single degenerate interval as one", {
  x <- ari_plot_fixture()
  x <- x[x$Score %in% c(1L, 2L), , drop = FALSE]
  x$n[x$ItemNo == 1L] <- c(2, 3)

  intervals <- ari_score_intervals_internal(x, class_size = 40L)

  expect_equal(unique(intervals$interval), 1L)
})

test_that("ari_score_intervals rejects invalid class sizes and duplicate first-item scores", {
  expect_error(ari_score_intervals_internal(ari_plot_fixture(), class_size = 0L), "class_size")
  expect_error(ari_score_intervals_internal(ari_plot_fixture(), class_size = 1.5), "class_size")

  x <- rbind(ari_plot_fixture(), ari_plot_fixture()[1L, , drop = FALSE])
  expect_error(ari_score_intervals_internal(x, class_size = 5L), "duplicate")
})

test_that("ari_plot_data computes SAS-compatible weighted summaries and qnorm bands", {
  plot_data <- ari_plot_data_internal(ari_plot_fixture(), class_size = 5L)

  expect_equal(names(plot_data), c("ItemNo", "Item", "interval", "O", "N", "E", "V", "lower", "upper"))
  expect_equal(nrow(plot_data), 4L)

  first <- plot_data[plot_data$Item == "I1" & plot_data$interval == 1L, , drop = FALSE]
  expect_equal(first$N, 5)
  expect_equal(first$O, (3 * 0.10 + 2 * 0.20) / 5)
  expect_equal(first$E, (3 * 0.15 + 2 * 0.25) / 5)
  expect_equal(first$V, (3 * 0.04 + 2 * 0.16) / 5)
  expect_equal(first$lower, first$E - stats::qnorm(0.975) * sqrt(first$V / first$N))
  expect_equal(first$upper, first$E + stats::qnorm(0.975) * sqrt(first$V / first$N))

  narrower <- ari_plot_data_internal(ari_plot_fixture(), class_size = 5L, confidence = 0.90)
  first_narrower <- narrower[narrower$Item == "I1" & narrower$interval == 1L, , drop = FALSE]
  expect_equal(first_narrower$lower, first$E - stats::qnorm(0.95) * sqrt(first$V / first$N))
})

test_that("ari_plot_data rejects invalid inputs", {
  expect_error(ari_plot_data_internal(data.frame()), "gRm_ari")
  expect_error(ari_plot_data_internal(ari_plot_fixture()[, -which(names(ari_plot_fixture()) == "ExpVar")]), "ExpVar")
  expect_error(ari_plot_data_internal(ari_plot_fixture(), confidence = 1), "confidence")

  x <- ari_plot_fixture()
  x$n[1L] <- -1
  expect_error(ari_plot_data_internal(x), "negative")
})

test_that("plot.gRm_ari filters items and returns a ggplot object without modifying input", {
  x <- ari_plot_fixture()
  before <- x

  p <- plot(x, class_size = 5L, items = "I2", rows = 1L, columns = 1L)

  expect_s3_class(p, "ggplot")
  expect_identical(x, before)
  expect_equal(unique(p$data$Item), "I2")
  expect_equal(unique(p$data$ItemNo), 2L)
})

test_that("plot.gRm_ari validates item selectors, layout, and reserved dots", {
  x <- ari_plot_fixture()

  expect_error(plot(x, class_size = 5L, items = "missing"), "items")
  expect_error(plot(x, class_size = 5L, items = list("I1")), "items")
  expect_error(plot(x, class_size = 5L, rows = 1L, columns = 1L), "rows.*columns|layout")
  expect_error(plot(x, class_size = 5L, colour = "red"), "unused|reserved|\\.\\.\\.")
})

test_that("plot.gRm_ari includes SAS-compatible default layers and optional expected line", {
  x <- ari_plot_fixture()

  default_plot <- plot(x, class_size = 5L)
  expected_plot <- plot(x, class_size = 5L, show_expected = TRUE)

  default_geoms <- vapply(default_plot$layers, function(layer) class(layer$geom)[[1L]], character(1L))
  expected_geoms <- vapply(expected_plot$layers, function(layer) class(layer$geom)[[1L]], character(1L))

  expect_equal(length(default_plot$layers), 2L)
  expect_equal(length(expected_plot$layers), 3L)
  expect_false("GeomText" %in% default_geoms)
  expect_false("GeomText" %in% expected_geoms)
  expect_s3_class(default_plot$facet$params$labeller, "labeller")
  expect_s3_class(default_plot$theme$strip.text, "element_text")
  expect_s3_class(default_plot$theme$panel.background, "element_blank")
  expect_s3_class(default_plot$theme$plot.background, "element_blank")
  expect_equal(default_plot$labels$x, "class interval")
  expect_equal(default_plot$labels$y, "Mean item score")
})
