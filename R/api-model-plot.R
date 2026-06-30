# Plot score-inclusive gRm model graphs
#
# Graph data construction lives in api-model-graph.R. This file contains the
# S3 plot methods and base-graphics drawing helpers for that graph.

#' @rdname model_graph
#' @param x A `gRm_model` or `gRm_fit` object.
#' @param layout Graph layout used by `plot()`. `"score_left"` keeps the
#'   source-faithful deterministic layout with `Score` on the left and uses a
#'   gRm-specific layered renderer. `"fr"`, `"kk"`, and `"nicely"` use
#'   automatic `igraph` layouts and the `igraph` renderer.
#' @param x_spacing Horizontal spacing multiplier for `layout = "score_left"`.
#' @param y_spacing Vertical spacing multiplier for `layout = "score_left"`.
#' @param item_spacing Additional vertical spacing multiplier for item nodes
#'   when `layout = "score_left"`. Increase this when item nodes are visually
#'   too close together.
#' @param vertex.label Optional vertex labels.
#' @param vertex.color Optional vertex colors.
#' @param vertex.shape Optional vertex shapes. For `layout = "score_left"`,
#'   `"circle"` and `"rectangle"` are supported.
#' @param vertex.size Optional vertex sizes.
#' @param vertex.label.cex Optional vertex-label scaling.
#' @param vertex.label.dist Optional vertex-label distance. For
#'   `layout = "score_left"`, this scales the distance between the node marker
#'   and its label.
#' @param vertex.label.degree Optional vertex-label angle passed to
#'   `plot.igraph` for automatic `igraph` layouts.
#' @param edge.color Optional edge colors.
#' @param edge.lty Optional edge line types.
#' @param edge.width Optional edge widths.
#' @param margin Optional plot margin.
#' @param rescale Logical value passed to `plot.igraph`. The default is
#'   `FALSE`. This is used only for automatic `igraph` layouts.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#' @return `plot()` methods draw the graph and return the `igraph` object
#'   invisibly.
#' @export
plot.gRm_model <- function(x,
                           ...,
                           layout = c("score_left", "fr", "kk", "nicely"),
                           x_spacing = 2.5,
                           y_spacing = 1.4,
                           item_spacing = 1.6,
                           vertex.label = NULL,
                           vertex.color = NULL,
                           vertex.shape = NULL,
                           vertex.size = NULL,
                           vertex.label.cex = NULL,
                           vertex.label.dist = NULL,
                           vertex.label.degree = NULL,
                           edge.color = NULL,
                           edge.lty = NULL,
                           edge.width = NULL,
                           margin = NULL,
                           rescale = FALSE,
                           xlim = NULL,
                           ylim = NULL) {
  graph <- model_graph(x)
  plot_gRm_model_graph(
    graph,
    ...,
    layout = layout,
    x_spacing = x_spacing,
    y_spacing = y_spacing,
    item_spacing = item_spacing,
    vertex.label = vertex.label,
    vertex.color = vertex.color,
    vertex.shape = vertex.shape,
    vertex.size = vertex.size,
    vertex.label.cex = vertex.label.cex,
    vertex.label.dist = vertex.label.dist,
    vertex.label.degree = vertex.label.degree,
    edge.color = edge.color,
    edge.lty = edge.lty,
    edge.width = edge.width,
    margin = margin,
    rescale = rescale,
    xlim = xlim,
    ylim = ylim
  )
  invisible(graph)
}

#' @export
plot.gRm_fit <- function(x,
                         ...,
                         layout = c("score_left", "fr", "kk", "nicely"),
                         x_spacing = 2.5,
                         y_spacing = 1.4,
                         item_spacing = 1.6,
                         vertex.label = NULL,
                         vertex.color = NULL,
                         vertex.shape = NULL,
                         vertex.size = NULL,
                         vertex.label.cex = NULL,
                         vertex.label.dist = NULL,
                         vertex.label.degree = NULL,
                         edge.color = NULL,
                         edge.lty = NULL,
                         edge.width = NULL,
                         margin = NULL,
                         rescale = FALSE,
                         xlim = NULL,
                         ylim = NULL) {
  graph <- model_graph(x)
  plot_gRm_model_graph(
    graph,
    ...,
    layout = layout,
    x_spacing = x_spacing,
    y_spacing = y_spacing,
    item_spacing = item_spacing,
    vertex.label = vertex.label,
    vertex.color = vertex.color,
    vertex.shape = vertex.shape,
    vertex.size = vertex.size,
    vertex.label.cex = vertex.label.cex,
    vertex.label.dist = vertex.label.dist,
    vertex.label.degree = vertex.label.degree,
    edge.color = edge.color,
    edge.lty = edge.lty,
    edge.width = edge.width,
    margin = margin,
    rescale = rescale,
    xlim = xlim,
    ylim = ylim
  )
  invisible(graph)
}

plot_gRm_model_graph <- function(graph,
                                 ...,
                                 layout = c("score_left", "fr", "kk", "nicely"),
                                 x_spacing = 2.5,
                                 y_spacing = 1.4,
                                 item_spacing = 1.6,
                                 vertex.label = NULL,
                                 vertex.color = NULL,
                                 vertex.shape = NULL,
                                 vertex.size = NULL,
                                 vertex.label.cex = NULL,
                                 vertex.label.dist = NULL,
                                 vertex.label.degree = NULL,
                                 edge.color = NULL,
                                 edge.lty = NULL,
                                 edge.width = NULL,
                                 margin = NULL,
                                 rescale = FALSE,
                                 xlim = NULL,
                                 ylim = NULL) {
  layout <- match.arg(layout)
  if (identical(layout, "score_left")) {
    layers <- gRm_graph_score_left_plot_layers(
      graph,
      x_spacing = x_spacing,
      y_spacing = y_spacing,
      item_spacing = item_spacing,
      vertex.label = vertex.label,
      vertex.color = vertex.color,
      vertex.shape = vertex.shape,
      vertex.size = vertex.size,
      vertex.label.cex = vertex.label.cex,
      vertex.label.dist = vertex.label.dist,
      edge.color = edge.color,
      edge.lty = edge.lty,
      edge.width = edge.width,
      margin = margin,
      xlim = xlim,
      ylim = ylim
    )
    gRm_draw_score_left_graph(layers, ...)
    return(invisible(layers))
  }

  vertex_type <- igraph::V(graph)$type
  edge_type <- igraph::E(graph)$type
  graph_layout <- gRm_graph_plot_layout(
    graph,
    layout = layout,
    x_spacing = x_spacing,
    y_spacing = y_spacing,
    item_spacing = item_spacing
  )

  if (is.null(vertex.label)) {
    vertex.label <- igraph::V(graph)$label
  }
  if (is.null(vertex.color)) {
    vertex.color <- c(
      score = "grey85",
      item = "white",
      exogenous = "grey95"
    )[vertex_type]
  }
  if (is.null(vertex.shape)) {
    vertex.shape <- c(
      score = "circle",
      item = "circle",
      exogenous = "rectangle"
    )[vertex_type]
  }
  if (is.null(edge.color)) {
    edge.color <- c(
      score = "grey70",
      ld = "black",
      dif = "grey35"
    )[edge_type]
  }
  if (is.null(edge.lty)) {
    edge.lty <- c(
      score = 1,
      ld = 1,
      dif = 2
    )[edge_type]
  }
  if (is.null(edge.width)) {
    edge.width <- c(
      score = 1,
      ld = 2,
      dif = 1.5
    )[edge_type]
  }
  if (is.null(vertex.size)) {
    vertex.size <- 9
  }
  if (is.null(vertex.label.cex)) {
    vertex.label.cex <- 0.65
  }
  if (is.null(vertex.label.dist)) {
    vertex.label.dist <- 0.75
  }
  if (is.null(vertex.label.degree)) {
    vertex.label.degree <- 0
  }
  if (is.null(margin)) {
    margin <- 0.35
  }
  if (isTRUE(rescale)) {
    if (is.null(xlim)) {
      xlim <- c(-1, 1)
    }
    if (is.null(ylim)) {
      ylim <- c(-1, 1)
    }
  } else if (is.null(xlim) || is.null(ylim)) {
    plot_limits <- gRm_graph_plot_limits(
      graph_layout,
      vertex.size = vertex.size,
      vertex.label.dist = vertex.label.dist,
      margin = margin
    )
    if (is.null(xlim)) {
      xlim <- plot_limits$xlim
    }
    if (is.null(ylim)) {
      ylim <- plot_limits$ylim
    }
  }

  graphics::plot(
    graph,
    layout = graph_layout,
    vertex.label = vertex.label,
    vertex.color = vertex.color,
    vertex.shape = vertex.shape,
    edge.color = edge.color,
    edge.lty = edge.lty,
    edge.width = edge.width,
    vertex.size = vertex.size,
    vertex.label.cex = vertex.label.cex,
    vertex.label.dist = vertex.label.dist,
    vertex.label.degree = vertex.label.degree,
    margin = margin,
    rescale = rescale,
    xlim = xlim,
    ylim = ylim,
    ...
  )
}

gRm_graph_score_left_plot_layers <- function(graph,
                                             x_spacing = 2.5,
                                             y_spacing = 1.4,
                                             item_spacing = 1.6,
                                             vertex.label = NULL,
                                             vertex.color = NULL,
                                             vertex.shape = NULL,
                                             vertex.size = NULL,
                                             vertex.label.cex = NULL,
                                             vertex.label.dist = NULL,
                                             edge.color = NULL,
                                             edge.lty = NULL,
                                             edge.width = NULL,
                                             margin = NULL,
                                             xlim = NULL,
                                             ylim = NULL) {
  vertices <- igraph::as_data_frame(graph, what = "vertices")
  edges <- igraph::as_data_frame(graph, what = "edges")
  graph_layout <- gRm_graph_plot_layout(
    graph,
    layout = "score_left",
    x_spacing = x_spacing,
    y_spacing = y_spacing,
    item_spacing = item_spacing
  )

  if (is.null(vertex.label)) {
    vertex.label <- vertices$label
  }
  if (is.null(vertex.color)) {
    vertex.color <- c(
      score = "#D7DEE7",
      item = "#FFFFFF",
      exogenous = "#F2F4F7"
    )[vertices$type]
  }
  if (is.null(vertex.shape)) {
    vertex.shape <- c(
      score = "circle",
      item = "circle",
      exogenous = "rectangle"
    )[vertices$type]
  }
  if (is.null(vertex.size)) {
    vertex.size <- 8
  }
  if (is.null(vertex.label.cex)) {
    vertex.label.cex <- 0.78
  }
  if (is.null(vertex.label.dist)) {
    vertex.label.dist <- 1
  }
  if (is.null(edge.color)) {
    edge.color <- c(
      score = "#BFC5CF",
      ld = "#1F2933",
      dif = "#737B86"
    )[edges$type]
  }
  if (is.null(edge.lty)) {
    edge.lty <- c(
      score = 1,
      ld = 1,
      dif = 2
    )[edges$type]
  }
  if (is.null(edge.width)) {
    edge.width <- c(
      score = 1,
      ld = 2,
      dif = 1.3
    )[edges$type]
  }
  if (is.null(margin)) {
    margin <- 0.16
  }

  vertex_count <- nrow(vertices)
  edge_count <- nrow(edges)
  vertices$x <- graph_layout[, 1]
  vertices$y <- graph_layout[, 2]
  vertices$plot_label <- gRm_recycle_plot_value(vertex.label, vertex_count)
  vertices$node_fill <- gRm_recycle_plot_value(vertex.color, vertex_count)
  vertices$node_shape <- gRm_graph_base_pch(gRm_recycle_plot_value(vertex.shape, vertex_count))
  vertices$node_cex <- gRm_recycle_plot_value(vertex.size, vertex_count) / 10
  vertices$label_cex <- gRm_recycle_plot_value(vertex.label.cex, vertex_count)
  vertices$label_dist <- gRm_recycle_plot_value(vertex.label.dist, vertex_count)

  label_offset <- max(0.18, x_spacing * 0.075)
  vertices$label_x <- vertices$x +
    ifelse(vertices$type == "score", -label_offset, label_offset) * vertices$label_dist
  vertices$label_y <- vertices$y
  vertices$label_hjust <- ifelse(vertices$type == "score", 1, 0)
  vertices$plot_row <- gRm_graph_item_plot_rows(vertices)

  if (edge_count > 0L) {
    from <- match(edges$from, vertices$name)
    to <- match(edges$to, vertices$name)
    edges$x0 <- vertices$x[from]
    edges$y0 <- vertices$y[from]
    edges$x1 <- vertices$x[to]
    edges$y1 <- vertices$y[to]
    edges$from_order <- vertices$order[from]
    edges$to_order <- vertices$order[to]
    edges$from_plot_row <- vertices$plot_row[from]
    edges$to_plot_row <- vertices$plot_row[to]
    edges$edge_color <- gRm_recycle_plot_value(edge.color, edge_count)
    edges$edge_lty <- gRm_recycle_plot_value(edge.lty, edge_count)
    edges$edge_width <- gRm_recycle_plot_value(edge.width, edge_count)
    adjacent_ld <- edges$type == "ld" & abs(edges$from_plot_row - edges$to_plot_row) == 1L
    curved_ld <- edges$type == "ld" & !adjacent_ld
    edges$edge_shape <- ifelse(curved_ld, "curve", "segment")
    ld_offset <- max(0.35, x_spacing * 0.22)
    edges$control_x <- ifelse(curved_ld, pmin(edges$x0, edges$x1) - ld_offset, NA_real_)
    edges$control_y <- ifelse(curved_ld, (edges$y0 + edges$y1) / 2, NA_real_)
    edges$draw_order <- match(edges$type, c("score", "dif", "ld"))
    edges <- edges[order(edges$draw_order), , drop = FALSE]
  } else {
    edges$edge_color <- character()
    edges$edge_lty <- numeric()
    edges$edge_width <- numeric()
    edges$from_order <- integer()
    edges$to_order <- integer()
    edges$from_plot_row <- integer()
    edges$to_plot_row <- integer()
    edges$edge_shape <- character()
    edges$control_x <- numeric()
    edges$control_y <- numeric()
    edges$draw_order <- integer()
  }

  limits <- gRm_score_left_plot_limits(
    vertices,
    x_spacing = x_spacing,
    margin = margin,
    xlim = xlim,
    ylim = ylim
  )

  list(
    vertices = vertices,
    edges = edges,
    xlim = limits$xlim,
    ylim = limits$ylim
  )
}

gRm_graph_item_plot_rows <- function(vertices) {
  plot_row <- rep(NA_integer_, nrow(vertices))
  item_index <- which(vertices$type == "item")
  if (!length(item_index)) {
    return(plot_row)
  }
  item_order <- order(-vertices$y[item_index], vertices$order[item_index])
  plot_row[item_index[item_order]] <- seq_along(item_index)
  plot_row
}

gRm_draw_score_left_graph <- function(layers, ...) {
  dots <- list(...)
  plot_args <- dots[names(dots) %in% c("main", "sub", "xlab", "ylab", "axes", "frame.plot")]
  plot_args$x <- NA_real_
  plot_args$y <- NA_real_
  plot_args$type <- "n"
  plot_args$xlim <- layers$xlim
  plot_args$ylim <- layers$ylim
  if (is.null(plot_args$xlab)) {
    plot_args$xlab <- ""
  }
  if (is.null(plot_args$ylab)) {
    plot_args$ylab <- ""
  }
  if (is.null(plot_args$axes)) {
    plot_args$axes <- FALSE
  }
  if (is.null(plot_args$frame.plot)) {
    plot_args$frame.plot <- FALSE
  }

  old_par <- graphics::par(xaxs = "i", yaxs = "i")
  on.exit(graphics::par(old_par), add = TRUE)
  do.call(graphics::plot, plot_args)

  edges <- layers$edges
  if (nrow(edges) > 0L) {
    segment_edges <- edges[edges$edge_shape == "segment", , drop = FALSE]
    curve_edges <- edges[edges$edge_shape == "curve", , drop = FALSE]
  } else {
    segment_edges <- edges
    curve_edges <- edges
  }
  if (nrow(segment_edges) > 0L) {
    graphics::segments(
      segment_edges$x0,
      segment_edges$y0,
      segment_edges$x1,
      segment_edges$y1,
      col = segment_edges$edge_color,
      lty = segment_edges$edge_lty,
      lwd = segment_edges$edge_width,
      lend = "round"
    )
  }
  if (nrow(curve_edges) > 0L) {
    for (edge_index in seq_len(nrow(curve_edges))) {
      edge <- curve_edges[edge_index, , drop = FALSE]
      curve <- gRm_quadratic_edge_path(
        x0 = edge$x0,
        y0 = edge$y0,
        control_x = edge$control_x,
        control_y = edge$control_y,
        x1 = edge$x1,
        y1 = edge$y1
      )
      graphics::lines(
        curve$x,
        curve$y,
        col = edge$edge_color,
        lty = edge$edge_lty,
        lwd = edge$edge_width,
        lend = "round"
      )
    }
  }

  vertices <- layers$vertices
  graphics::points(
    vertices$x,
    vertices$y,
    pch = vertices$node_shape,
    cex = vertices$node_cex,
    col = "#2E3440",
    bg = vertices$node_fill
  )

  visible_labels <- !is.na(vertices$plot_label) & nzchar(vertices$plot_label)
  label_groups <- split(which(visible_labels), vertices$label_hjust[visible_labels])
  for (group in label_groups) {
    graphics::text(
      vertices$label_x[group],
      vertices$label_y[group],
      labels = vertices$plot_label[group],
      adj = c(vertices$label_hjust[group][1], 0.5),
      cex = vertices$label_cex[group],
      col = "#111827"
    )
  }
}

gRm_quadratic_edge_path <- function(x0,
                                    y0,
                                    control_x,
                                    control_y,
                                    x1,
                                    y1,
                                    n = 80L) {
  t <- seq(0, 1, length.out = n)
  data.frame(
    x = (1 - t)^2 * x0 + 2 * (1 - t) * t * control_x + t^2 * x1,
    y = (1 - t)^2 * y0 + 2 * (1 - t) * t * control_y + t^2 * y1
  )
}

gRm_graph_base_pch <- function(shape) {
  shape <- as.character(shape)
  pch <- rep(21L, length(shape))
  pch[shape %in% c("rectangle", "square")] <- 22L
  numeric_shape <- suppressWarnings(as.integer(shape))
  pch[!is.na(numeric_shape)] <- numeric_shape[!is.na(numeric_shape)]
  pch
}

gRm_score_left_plot_limits <- function(vertices,
                                       x_spacing,
                                       margin,
                                       xlim = NULL,
                                       ylim = NULL) {
  margin <- rep(margin, length.out = 4L)
  x_range <- range(vertices$x, finite = TRUE)
  y_range <- range(vertices$y, finite = TRUE)
  label_widths <- nchar(vertices$plot_label)
  label_widths <- label_widths[is.finite(label_widths)]
  max_label_width <- if (length(label_widths)) max(label_widths) else 0
  left_pad <- max(0.75, x_spacing * (0.35 + margin[2]))
  right_pad <- max(1.1, x_spacing * (0.45 + margin[4]), max_label_width * 0.08)
  y_pad <- max(0.6, x_spacing * 0.18, margin[1])

  if (is.null(xlim)) {
    xlim <- c(x_range[1] - left_pad, x_range[2] + right_pad)
  }
  if (is.null(ylim)) {
    ylim <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  }

  list(xlim = xlim, ylim = ylim)
}

gRm_recycle_plot_value <- function(value, n) {
  if (n == 0L) {
    return(value[FALSE])
  }
  rep(value, length.out = n)
}

gRm_graph_plot_limits <- function(graph_layout,
                                  vertex.size,
                                  vertex.label.dist,
                                  margin) {
  margin <- rep(margin, length.out = 4L)
  x_range <- range(graph_layout[, 1], finite = TRUE)
  y_range <- range(graph_layout[, 2], finite = TRUE)
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  if (!is.finite(x_span) || x_span <= 0) {
    x_span <- 1
  }
  if (!is.finite(y_span) || y_span <= 0) {
    y_span <- 1
  }
  vertex_pad <- max(vertex.size, na.rm = TRUE) / 200
  label_pad <- vertex_pad * max(vertex.label.dist, na.rm = TRUE)
  base_pad <- max(0.5, vertex_pad + label_pad)
  x_pad <- base_pad + x_span * c(margin[2], margin[4])
  y_pad <- base_pad + max(x_span, y_span) * c(margin[1], margin[3])

  list(
    xlim = c(x_range[1] - x_pad[1], x_range[2] + x_pad[2]),
    ylim = c(y_range[1] - y_pad[1], y_range[2] + y_pad[2])
  )
}

gRm_graph_plot_layout <- function(graph,
                                  layout = c("score_left", "fr", "kk", "nicely"),
                                  x_spacing = 2.5,
                                  y_spacing = 1.4,
                                  item_spacing = 1) {
  layout <- match.arg(layout)
  if (identical(layout, "score_left")) {
    vertex_type <- igraph::V(graph)$type
    y <- igraph::V(graph)$y * y_spacing
    y[vertex_type == "item"] <- y[vertex_type == "item"] * item_spacing
    return(cbind(
      igraph::V(graph)$x * x_spacing,
      y
    ))
  }
  switch(
    layout,
    fr = igraph::layout_with_fr(graph),
    kk = igraph::layout_with_kk(graph),
    nicely = igraph::layout_nicely(graph)
  )
}
