#' Build or plot the score-inclusive GLLRM graph
#'
#' `model_graph()` returns the score-inclusive IRT graph for a gRm model or
#' fitted model. The graph always contains the score node. Items are connected
#' to score by baseline score edges, LD terms add item-item edges, and DIF terms
#' add item-exogenous edges.
#'
#' @param object A `gRm_model` or `gRm_fit` object.
#' @param ... Reserved for S3 dispatch compatibility; ignored by
#'   `model_graph()` methods.
#' @return An undirected `igraph` object with vertex attributes `label`, `type`,
#'   `order`, `x`, and `y`, and edge attributes `type`, `source`, `status`,
#'   `from_label`, and `to_label`.
#' @details
#' The official gRm graph is the IRT graph with score. The package deliberately
#' does not expose a score-free graph because removing score changes the graph's
#' interpretation. A Rasch model is therefore a star centered on `Score`; GLLRM
#' local-dependence and DIF terms are added to that baseline graph.
#'
#' The default score-left plot is a layered rendering of the same graph. LD
#' edges between adjacent plotted item rows are drawn as straight item-item
#' links. LD edges between non-adjacent plotted item rows are curved away from
#' the item column so they are not mistaken for dependence among all intervening
#' items.
#' @export
#' @examples
#' data <- data.frame(
#'   ID = 1:8,
#'   I1 = c(0, 1, 0, 1, 0, 1, 1, 0),
#'   I2 = c(1, 0, 1, 0, 1, 0, 1, 0),
#'   site = c(0, 0, 1, 1, 0, 1, 0, 1)
#' )
#' analysis <- gRm(
#'   data,
#'   items = c("I1", "I2"),
#'   exogenous = "site",
#'   id = "ID",
#'   score_cuts = c(1L, 2L)
#' )
#' model <- gllrm(analysis, ld = ~ I1:I2, dif = ~ I1:site)
#' graph <- model_graph(model)
model_graph <- function(object, ...) {
  UseMethod("model_graph")
}

#' @export
model_graph.gRm_model <- function(object, ...) {
  build_gRm_model_graph(object)
}

#' @export
model_graph.gRm_fit <- function(object, ...) {
  model <- object$spec %||% object$model
  if (is.null(model)) {
    stop("Could not find the model specification attached to this fit.", call. = FALSE)
  }
  model_graph(model, ...)
}

build_gRm_model_graph <- function(model) {
  analysis <- model$analysis
  items <- analysis$items
  ld <- model$ld %||% empty_ld_terms()
  dif <- model$dif %||% empty_dif_terms()
  used_exogenous <- analysis$exogenous[analysis$exogenous %in% unique(dif$exogenous)]

  vertices <- gRm_graph_vertices(items, used_exogenous)
  edges <- rbind_fill(
    gRm_graph_score_edges(items),
    gRm_graph_ld_edges(ld),
    gRm_graph_dif_edges(dif)
  )

  graph <- igraph::graph_from_data_frame(
    d = edges,
    directed = FALSE,
    vertices = vertices
  )
  igraph::graph_attr(graph, "gRm_layout") <- "score_left"
  graph
}

gRm_graph_vertices <- function(items, exogenous) {
  items <- as.character(items %||% character())
  exogenous <- as.character(exogenous %||% character())
  item_y <- gRm_centered_y(length(items), top_to_bottom = TRUE)
  exo_y <- gRm_centered_y(length(exogenous), top_to_bottom = TRUE)

  data.frame(
    name = c(
      "score:Score",
      gRm_graph_ids("item", items),
      gRm_graph_ids("exogenous", exogenous)
    ),
    label = c("Score", items, exogenous),
    type = c("score", rep("item", length(items)), rep("exogenous", length(exogenous))),
    order = c(1L, seq_along(items), seq_along(exogenous)),
    x = c(0, rep(1, length(items)), rep(2, length(exogenous))),
    y = c(0, item_y, exo_y),
    stringsAsFactors = FALSE
  )
}

gRm_graph_ids <- function(type, values) {
  values <- as.character(values %||% character())
  if (!length(values)) {
    return(character())
  }
  paste0(type, ":", values)
}

gRm_centered_y <- function(n, top_to_bottom = TRUE) {
  if (n <= 0L) {
    return(numeric())
  }
  y <- seq_len(n) - mean(seq_len(n))
  if (isTRUE(top_to_bottom)) {
    y <- rev(y)
  }
  as.numeric(y)
}

gRm_graph_edge_frame <- function() {
  data.frame(
    from = character(),
    to = character(),
    type = character(),
    source = character(),
    status = character(),
    from_label = character(),
    to_label = character(),
    stringsAsFactors = FALSE
  )
}

gRm_graph_score_edges <- function(items) {
  data.frame(
    from = "score:Score",
    to = gRm_graph_ids("item", items),
    type = "score",
    source = "baseline",
    status = "baseline",
    from_label = "Score",
    to_label = items,
    stringsAsFactors = FALSE
  )
}

gRm_graph_ld_edges <- function(ld) {
  if (!is.data.frame(ld) || nrow(ld) == 0L) {
    return(gRm_graph_edge_frame())
  }
  data.frame(
    from = gRm_graph_ids("item", ld$item1),
    to = gRm_graph_ids("item", ld$item2),
    type = "ld",
    source = ld$source %||% rep(NA_character_, nrow(ld)),
    status = ld$status %||% rep(NA_character_, nrow(ld)),
    from_label = ld$item1,
    to_label = ld$item2,
    stringsAsFactors = FALSE
  )
}

gRm_graph_dif_edges <- function(dif) {
  if (!is.data.frame(dif) || nrow(dif) == 0L) {
    return(gRm_graph_edge_frame())
  }
  data.frame(
    from = gRm_graph_ids("item", dif$item),
    to = gRm_graph_ids("exogenous", dif$exogenous),
    type = "dif",
    source = dif$source %||% rep(NA_character_, nrow(dif)),
    status = dif$status %||% rep(NA_character_, nrow(dif)),
    from_label = dif$item,
    to_label = dif$exogenous,
    stringsAsFactors = FALSE
  )
}
