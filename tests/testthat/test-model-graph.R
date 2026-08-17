graph_api_data <- function() {
  rows <- expand.grid(
    I1 = 0:1,
    I2 = 0:1,
    I3 = 0:1,
    site = 0:1,
    age = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  data.frame(ID = seq_len(nrow(rows)), rows)
}

graph_api_analysis <- function() {
  gRm(
    graph_api_data(),
    items = c("I1", "I2", "I3"),
    exogenous = c("site", "age"),
    id = "ID"
  )
}

edge_keys <- function(graph) {
  edges <- igraph::as_data_frame(graph, what = "edges")
  paste(edges$type, edges$from_label, edges$to_label, sep = ":")
}

test_that("model_graph creates the score-inclusive Rasch graph", {
  analysis <- graph_api_analysis()
  model <- gllrm(analysis)

  graph <- model_graph(model)
  vertices <- igraph::as_data_frame(graph, what = "vertices")
  edges <- igraph::as_data_frame(graph, what = "edges")

  expect_s3_class(graph, "igraph")
  expect_equal(vertices$label, c("Score", "I1", "I2", "I3"))
  expect_equal(vertices$type, c("score", "item", "item", "item"))
  expect_equal(vertices$x, c(0, 1, 1, 1))
  expect_equal(sum(edges$type == "score"), 3L)
  expect_equal(sum(edges$type == "ld"), 0L)
  expect_equal(sum(edges$type == "dif"), 0L)
  expect_setequal(edge_keys(graph), c(
    "score:Score:I1",
    "score:Score:I2",
    "score:Score:I3"
  ))
})

test_that("model_graph adds LD and DIF edges to the score graph", {
  analysis <- graph_api_analysis()
  model <- gllrm(
    analysis,
    ld = ~ I1:I2 + I2:I3,
    dif = ~ I1:site + I3:age
  )

  graph <- model_graph(model)
  vertices <- igraph::as_data_frame(graph, what = "vertices")

  expect_equal(vertices$label, c("Score", "I1", "I2", "I3", "site", "age"))
  expect_equal(vertices$type, c("score", "item", "item", "item", "exogenous", "exogenous"))
  expect_equal(vertices$x, c(0, 1, 1, 1, 2, 2))
  expect_setequal(edge_keys(graph), c(
    "score:Score:I1",
    "score:Score:I2",
    "score:Score:I3",
    "ld:I1:I2",
    "ld:I2:I3",
    "dif:I1:site",
    "dif:I3:age"
  ))
})

test_that("model_graph works for fitted models", {
  analysis <- graph_api_analysis()
  model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I3:site)
  fit_obj <- fit(model, max_step = 50L)

  model_graph_obj <- model_graph(model)
  fit_graph_obj <- model_graph(fit_obj)

  expect_equal(
    igraph::as_data_frame(fit_graph_obj, what = "vertices"),
    igraph::as_data_frame(model_graph_obj, what = "vertices")
  )
  expect_equal(edge_keys(fit_graph_obj), edge_keys(model_graph_obj))
})

test_that("plot methods draw and invisibly return the graph", {
  analysis <- graph_api_analysis()
  model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I1:site)

  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit(grDevices::dev.off(), add = TRUE)

  returned <- plot(model)

  expect_s3_class(returned, "igraph")
  expect_equal(edge_keys(returned), edge_keys(model_graph(model)))
})

test_that("plot supports spacing controls for the score-left layout", {
  analysis <- graph_api_analysis()
  model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I1:site)
  graph <- model_graph(model)

  layout <- gRm_graph_plot_layout(graph, layout = "score_left", x_spacing = 3, y_spacing = 2)

  expect_equal(layout[, 1], igraph::V(graph)$x * 3)
  expect_equal(layout[, 2], igraph::V(graph)$y * 2)
})

test_that("score-left layout can separate the item layer", {
  analysis <- graph_api_analysis()
  model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I1:site)
  graph <- model_graph(model)
  vertices <- igraph::as_data_frame(graph, what = "vertices")
  item_index <- which(vertices$type == "item")

  compact <- gRm_graph_plot_layout(
    graph,
    layout = "score_left",
    x_spacing = 1,
    y_spacing = 1,
    item_spacing = 1
  )
  separated <- gRm_graph_plot_layout(
    graph,
    layout = "score_left",
    x_spacing = 1,
    y_spacing = 1,
    item_spacing = 2
  )

  expect_equal(diff(sort(separated[item_index, 2])), diff(sort(compact[item_index, 2])) * 2)
  expect_equal(separated[-item_index, 2], compact[-item_index, 2])
})

test_that("score-left item spacing changes the rendered plot", {
  skip_if_not(capabilities("png"))
  analysis <- graph_api_analysis()
  model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I1:site)

  render_graph <- function(item_spacing) {
    path <- tempfile(fileext = ".png")
    grDevices::png(path, width = 700, height = 500)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot(model, item_spacing = item_spacing)
    path
  }

  compact <- render_graph(1.5)
  separated <- render_graph(3.5)

  expect_false(identical(
    unname(tools::md5sum(compact)),
    unname(tools::md5sum(separated))
  ))
})

test_that("score-left renderer places labels beside node markers", {
  analysis <- graph_api_analysis()
  model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I1:site)
  graph <- model_graph(model)

  layers <- gRm_graph_score_left_plot_layers(graph)
  vertices <- layers$vertices

  score <- vertices[vertices$type == "score", ]
  items <- vertices[vertices$type == "item", ]

  expect_true(score$label_x < score$x)
  expect_true(all(items$label_x > items$x))
  expect_true(all(abs(items$label_x - items$x) > 0))
  expect_true(all(vertices$node_cex < 1.2))
})

test_that("score-left renderer routes LD edges away from the item column", {
  item_names <- sprintf("I%02d", 1:8)
  data <- as.data.frame(matrix(rep(c(0L, 1L), length.out = 80L), nrow = 10L))
  names(data) <- item_names
  data$ID <- seq_len(nrow(data))
  data$site <- rep(0:1, length.out = nrow(data))
  analysis <- gRm(data, items = item_names, exogenous = "site", id = "ID")
  model <- gllrm(analysis, ld = ~ I01:I02 + I01:I08, dif = ~ I01:site)
  graph <- model_graph(model)

  layers <- gRm_graph_score_left_plot_layers(graph)
  vertices <- layers$vertices
  ld_edges <- layers$edges[layers$edges$type == "ld", ]
  item_x <- unique(vertices$x[vertices$type == "item"])
  adjacent <- ld_edges[ld_edges$from_label == "I01" & ld_edges$to_label == "I02", ]
  distant <- ld_edges[ld_edges$from_label == "I01" & ld_edges$to_label == "I08", ]

  expect_equal(nrow(ld_edges), 2L)
  expect_equal(adjacent$edge_shape, "segment")
  expect_true(is.na(adjacent$control_x))
  expect_equal(distant$edge_shape, "curve")
  expect_false(is.na(distant$control_x))
  expect_true(all(distant$control_x != item_x))
})

test_that("score-left renderer determines LD adjacency from plotted item rows", {
  item_names <- sprintf("I%02d", 1:8)
  data <- as.data.frame(matrix(rep(c(0L, 1L), length.out = 80L), nrow = 10L))
  names(data) <- item_names
  data$ID <- seq_len(nrow(data))
  data$site <- rep(0:1, length.out = nrow(data))
  analysis <- gRm(data, items = item_names, exogenous = "site", id = "ID")
  model <- gllrm(analysis, ld = ~ I01:I08)
  graph <- model_graph(model)

  item_vertices <- which(igraph::V(graph)$type == "item")
  igraph::V(graph)$y[item_vertices] <- c(4, -3, -4, -5, -6, -7, -8, 3)

  layers <- gRm_graph_score_left_plot_layers(graph)
  ld_edges <- layers$edges[layers$edges$type == "ld", ]

  expect_equal(ld_edges$edge_shape, "segment")
  expect_true(is.na(ld_edges$control_x))
})

test_that("plot supports automatic igraph layouts and plot sizing arguments", {
  analysis <- graph_api_analysis()
  model <- gllrm(analysis, ld = ~ I1:I2 + I2:I3, dif = ~ I1:site + I3:age)
  graph <- model_graph(model)

  for (layout_name in c("fr", "kk", "nicely")) {
    layout <- gRm_graph_plot_layout(graph, layout = layout_name)
    expect_equal(dim(layout), c(igraph::vcount(graph), 2L))
    expect_true(all(is.finite(layout)))
  }

  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_no_error(
    returned <- plot(
      model,
      layout = "fr",
      vertex.size = 10,
      vertex.label.cex = 0.45,
      margin = 0.6
    )
  )
  expect_s3_class(returned, "igraph")
})
