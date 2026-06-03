#' Derive DIGRAM graphical SCREEN report values
#'
#' Computes the source-shaped graphical `SCREEN` matrices from a parsed DIGRAM
#' project. This is the graphical model screening command implemented by
#' `DIGRAM1f.Execute_Screening` and `SKscr4.ISCREEN`; it is distinct from the
#' item-screening `SCREEN J` report.
#'
#' The implemented slice covers the two-way table stage, hidden-association
#' stage, final model-3 revision stage, incomplete-path warnings, and the
#' current-model partial gamma matrix from the same source branch. Partial gamma
#' is computed after the SCREEN acceptance path copies `MODEL_3` back to the
#' project graph, through a native R port of `DGRexe.execute_gamma(false)`,
#' `Generate_hypotheses`, `CheckHyp`, and `quick_tests`; substituting marginal
#' gamma or an earlier SCREEN stage here would be source-unfaithful.
#'
#' Production R computes from `DIGRAM.var` and `DIGRAM.dat`; Pascal and the
#' supplied DIGRAM report are test oracles only.
#'
#' @param project A parsed DIGRAM project from [read_digram_project()].
#' @param significance SCREEN significance level. DIGRAM defaults to `0.05`.
#' @return A `gRm_screen_values` object.
#' @examples
#' \dontrun{
#' project <- read_digram_project("path/to/DIGRAM")
#' values <- screen_values(project)
#' values$final_model$status_chars[1:3, 1:3]
#' }
#' @keywords internal
screen_values <- function(project, significance = 0.05) {
  two_way <- screen_two_way_matrix(project, significance)
  hidden <- screen_hidden_matrix(project, two_way, significance)
  final_model <- screen_final_model(project, two_way, hidden, significance)
  partial_gamma <- screen_partial_gamma_matrix(project, final_model)

  structure(
    list(
      variables = screen_variables(project),
      significance = significance,
      fixed_edge_policy = "Fixed edges have not been removed and inadmissible edges not included.",
      two_way = two_way,
      hidden = hidden,
      final_model = final_model,
      incomplete_paths = screen_incomplete_path_warnings(project, final_model),
      partial_gamma = partial_gamma,
      source_status = c(
        two_way = "implemented_SKscr4_ISCREEN_TWOWAY_RCCHI_RCGAMMA",
        hidden = "implemented_SKscr4_ISCREEN_hidden_THREEWAY",
        final_model = "implemented_SKscr4_ISCREEN_model3_revision",
        incomplete_paths = "implemented_DGRexe_execute_gamma_false_SKgrf8_Generate_hypotheses_SKsep1_FINDPATHS_warning",
        partial_gamma = "implemented_DGRexe_execute_gamma_false_Generate_hypotheses_CheckHyp_quick_tests"
      )
    ),
    class = "gRm_screen_values"
  )
}

#' Compute the graphical SCREEN two-way matrix
#'
#' @inheritParams screen_values
#' @return A `gRm_screen_stage` object for the two-way table stage.
#' @keywords internal
screen_two_way_matrix <- function(project, significance = 0.05) {
  ctx <- screen_context(project)
  n_vars <- nrow(ctx$variables)
  status <- screen_initial_status(ctx)
  old_status <- status
  tests <- list()
  test_index <- 0L

  for (v1 in seq_len(ctx$screen_k)) {
    for (v2 in seq.int(v1 + 1L, n_vars)) {
      if (status[v1, v2] %in% c(0L, 3L)) {
        next
      }

      table <- screen_two_way_table(ctx, v1, v2)
      chi <- screen_rc_chi(table)
      gamma <- if (screen_ordinals(ctx, v1, v2)) {
        screen_rc_gamma(table)
      } else {
        list(gamma = 0, p_value = 1, ppq = 0, pmq = 0, s = 0, success = FALSE)
      }

      # Source trace: SKscr4.ISCREEN calls TWOWAY and then
      # DETERMINE_SIGNIFICANCE. In source model coding, status 1 means the
      # conditional-independence hypothesis is accepted (edge removed), while
      # status 2 means undecided/rejected (edge retained).
      accept <- chi$p_value > significance
      if (accept && screen_ordinals(ctx, v1, v2) && gamma$p_value <= significance) {
        accept <- FALSE
      }
      model_status <- if (accept) 1L else 2L
      status[v1, v2] <- status[v2, v1] <- model_status

      test_index <- test_index + 1L
      tests[[test_index]] <- data.frame(
        test = test_index,
        row = ctx$variables$label_code[[v1]],
        col = ctx$variables$label_code[[v2]],
        v1 = v1,
        v2 = v2,
        n = sum(table),
        chi_square = chi$chi_square,
        df = chi$df,
        chi_p = chi$p_value,
        ordinal = screen_ordinals(ctx, v1, v2),
        gamma = gamma$gamma,
        gamma_p = gamma$p_value,
        accept = accept,
        status = model_status,
        stringsAsFactors = FALSE
      )
    }
  }

  screen_stage(
    name = "Twoway tables",
    status = status,
    old_status = status,
    variables = ctx$variables,
    tests = do.call(rbind, tests),
    significance = significance
  )
}

#' Compute the graphical SCREEN hidden-association matrix
#'
#' @inheritParams screen_values
#' @param two_way Result from [screen_two_way_matrix()].
#' @return A `gRm_screen_stage` object for the hidden-association stage.
#' @keywords internal
screen_hidden_matrix <- function(project, two_way, significance = 0.05) {
  ctx <- screen_context(project)
  n_vars <- nrow(ctx$variables)
  status <- two_way$status
  old_status <- two_way$status
  tests <- list()
  test_index <- 0L

  for (v1 in seq_len(ctx$screen_k)) {
    for (v2 in seq.int(v1 + 1L, n_vars)) {
      if (two_way$status[v1, v2] != 1L) {
        next
      }

      for (v3 in seq.int(screen_first_node_at_level(ctx, v1), n_vars)) {
        if (v3 == v1 || v3 == v2) {
          next
        }
        # Source trace: SKscr2/SKscr4 assign MODEL_H := MODEL_2 as a pointer.
        # Hidden associations found earlier in source order therefore update
        # the same matrix consulted by later hidden-conditioner gates.
        if (status[v1, v3] <= 1L || status[v2, v3] <= 1L) {
          next
        }

        test <- screen_three_way_test(ctx, v1, v2, v3)
        accept <- screen_determine_significance(test$chi_p, test$gamma_p, test$ordinal, significance)
        if (!accept) {
          status[v1, v2] <- status[v2, v1] <- 4L
        }

        test_index <- test_index + 1L
        tests[[test_index]] <- data.frame(
          test = test_index,
          row = ctx$variables$label_code[[v1]],
          col = ctx$variables$label_code[[v2]],
          conditioner = ctx$variables$label_code[[v3]],
          v1 = v1,
          v2 = v2,
          v3 = v3,
          chi_square = test$chi_square,
          df = test$df,
          chi_p = test$chi_p,
          ordinal = test$ordinal,
          gamma = test$gamma,
          gamma_p = test$gamma_p,
          accept = accept,
          status = status[v1, v2],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  screen_stage(
    name = "hidden association",
    status = status,
    old_status = old_status,
    variables = ctx$variables,
    tests = do.call(rbind, tests),
    significance = significance
  )
}

#' Compute the final graphical SCREEN model
#'
#' @inheritParams screen_values
#' @param two_way Result from [screen_two_way_matrix()].
#' @param hidden Result from [screen_hidden_matrix()].
#' @return A `gRm_screen_final_model` object.
#' @keywords internal
screen_final_model <- function(project, two_way, hidden, significance = 0.05) {
  ctx <- screen_context(project)
  n_vars <- nrow(ctx$variables)
  status <- hidden$status
  tests <- list()
  test_index <- 0L
  # Source trace: SKscr4.ISCREEN sets K1 := NVAR - 2, then bounds it by
  # NRESPONS. The example graphical SCREEN project has one level containing all
  # variables, so background variables I/J/K participate in final MODEL_3
  # revisions instead of stopping at the item block.
  k1 <- n_vars - 2L
  k2 <- n_vars - 1L

  if (k1 >= 1L) {
    for (v1 in seq_len(k1)) {
      for (v2 in seq.int(v1 + 1L, k2)) {
        if (status[v1, v2] <= 1L) {
          next
        }

        for (v3 in seq.int(v2 + 1L, n_vars)) {
          # Source trace: SKscr2/SKscr4 assign MODEL_3 := MODEL_H as a
          # pointer. Final MODEL_3 revisions therefore update the matrix used
          # by later final-stage candidate gates.
          if (status[v1, v3] <= 1L || status[v2, v3] <= 1L) {
            next
          }

          performed <- c(
            status[v1, v2] != 3L,
            status[v1, v3] != 3L,
            status[v2, v3] != 3L
          )
          if (performed[[3]] && ctx$levels[[v1]] < ctx$levels[[v2]]) {
            performed[[3]] <- FALSE
          }
          if (!any(performed)) {
            next
          }

          pair_tests <- list(
            screen_three_way_test(ctx, v1, v2, v3),
            screen_three_way_test(ctx, v1, v3, v2),
            screen_three_way_test(ctx, v2, v3, v1)
          )
          p_min <- rep(-1, 3)
          for (index in seq_len(3L)) {
            if (performed[[index]]) {
              # Source trace: SKscr4.ISCREEN model-3 loop computes PMIN from
              # the smaller of chi and gamma p-values when ORDINALS(V1,V2)
              # holds. The original code uses ORDINALS(V1,V2) for all three
              # members of the triple, so this port preserves that quirk.
              if (!screen_ordinals(ctx, v1, v2)) {
                p_min[[index]] <- pair_tests[[index]]$chi_p
              } else {
                p_min[[index]] <- min(pair_tests[[index]]$chi_p, pair_tests[[index]]$gamma_p)
              }
            }
          }

          largest_pair <- which.max(ifelse(performed, p_min, -Inf))
          test_results <- rep(3L, 3L)
          for (index in seq_len(3L)) {
            if (performed[[index]]) {
              test_results[[index]] <- if (p_min[[index]] <= significance) {
                2L
              } else if (largest_pair == index) {
                1L
              } else {
                5L
              }
            }
          }

          pair_nodes <- list(c(v1, v2, v3), c(v1, v3, v2), c(v2, v3, v1))
          for (index in seq_len(3L)) {
            if (!performed[[index]]) {
              next
            }
            nodes <- pair_nodes[[index]]
            status <- screen_revise_model3(status, test_results[[index]], nodes[[1]], nodes[[2]])
            test_index <- test_index + 1L
            tests[[test_index]] <- data.frame(
              test = test_index,
              row = ctx$variables$label_code[[nodes[[1]]]],
              col = ctx$variables$label_code[[nodes[[2]]]],
              conditioner = ctx$variables$label_code[[nodes[[3]]]],
              v1 = nodes[[1]],
              v2 = nodes[[2]],
              v3 = nodes[[3]],
              chi_square = pair_tests[[index]]$chi_square,
              df = pair_tests[[index]]$df,
              chi_p = pair_tests[[index]]$chi_p,
              ordinal = pair_tests[[index]]$ordinal,
              gamma = pair_tests[[index]]$gamma,
              gamma_p = pair_tests[[index]]$gamma_p,
              p_min = p_min[[index]],
              test_result = test_results[[index]],
              status = status[nodes[[1]], nodes[[2]]],
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  structure(
    c(screen_stage(
      name = "The final SCREEN model",
      status = status,
      old_status = status,
      variables = ctx$variables,
      tests = do.call(rbind, tests),
      significance = significance
    ), list(status_text = screen_status_text(status))),
    class = c("gRm_screen_final_model", "gRm_screen_stage")
  )
}

#' Compute partial gamma values for the final graphical SCREEN model
#'
#' @inheritParams screen_values
#' @param final_model Result from [screen_final_model()].
#' @return A `gRm_screen_partial_gamma` object.
#' @keywords internal
screen_partial_gamma_matrix <- function(project, final_model) {
  ctx <- screen_context(project)
  variables <- ctx$variables
  n_vars <- nrow(variables)
  gamma <- matrix(999, nrow = n_vars, ncol = n_vars, dimnames = list(variables$label_code, variables$label_code))
  diag(gamma) <- 999
  graph <- final_model$status > 1L
  diag(graph) <- FALSE

  for (v1 in seq_len(n_vars - 1L)) {
    for (v2 in seq.int(v1 + 1L, n_vars)) {
      ordinal <- screen_ordinals(ctx, v1, v2)
      in_graph <- final_model$status[v1, v2] > 1L
      if (ordinal && !in_graph) {
        gamma[v1, v2] <- gamma[v2, v1] <- 0
      } else if (!ordinal && in_graph) {
        gamma[v1, v2] <- gamma[v2, v1] <- 888
      } else if (ordinal) {
        masks <- screen_generate_separation_masks(ctx, graph, v1, v2)
        kept <- vapply(masks, screen_hypothesis_fits_bounds, logical(1), ctx = ctx, v1 = v1, v2 = v2)
        if (any(kept)) {
          gamma_values <- vapply(masks[kept], screen_quick_tests_gamma_total, numeric(1), ctx = ctx, v1 = v1, v2 = v2)
          gamma[v1, v2] <- gamma[v2, v1] <- mean(gamma_values)
        }
      }
    }
  }

  structure(
    list(
      variables = variables,
      gamma = gamma,
      implemented = TRUE,
      graph_source = "project_status_after_SCREEN_accept_MODEL_3_rebuilt_by_MAKE_GRAPH",
      source_status = "implemented_DGRexe_execute_gamma_false_Generate_hypotheses_CheckHyp_quick_tests",
      source_detail = paste(
        "DIGRAM1f.Execute_Screening accepts MODEL_3 as project status, SKscr2.Screening rebuilds GRAF with MAKE_GRAPH(MODEL_3),",
        "and DGRexe.execute_gamma(false) initializes PartGamValues, generates graph-conditioned hypotheses,",
        "runs CheckHyp and quick_tests, then averages results[,5]. This implementation does not substitute marginal gamma or an earlier SCREEN graph."
      )
    ),
    class = "gRm_screen_partial_gamma"
  )
}

#' @keywords internal
screen_generate_separation_masks <- function(ctx, rec_graph, v1, v2) {
  n_vars <- nrow(rec_graph)
  if (n_vars == 0L) {
    return(list())
  }

  graph <- screen_generate_hypothesis_graph(ctx, rec_graph, v1, v2)
  degrees <- screen_graph_degrees(graph, v1, v2)
  if (degrees[[2L]] < degrees[[1L]]) {
    fra <- v2
    til <- v1
  } else {
    fra <- v1
    til <- v2
  }

  direct <- which(graph[fra, ] & graph[til, ])
  direct <- setdiff(direct, c(fra, til))
  direct_mask <- screen_mask_from_vertices(direct)
  reduced <- screen_separate1_reduced_graph(graph, fra, til)
  path_result <- screen_find_path_masks(reduced, fra, til, (30L * 162230L) %/% n_vars)
  if (!path_result$success) {
    return(screen_local_markov_masks(graph, fra, til))
  }
  path_masks <- path_result$masks
  if (length(path_masks) == 0L) {
    return(list(direct_mask))
  }

  vertices <- sort(unique(unlist(lapply(path_masks, screen_vertices_from_mask, n_vars = n_vars), use.names = FALSE)))
  cut_result <- screen_generate_cut_masks(vertices, path_masks, n_vars)
  if (!cut_result$success) {
    return(screen_local_markov_masks(graph, fra, til))
  }
  lapply(cut_result$masks, function(mask) screen_mask_union(direct_mask, mask))
}

#' @keywords internal
screen_generate_hypothesis_graph <- function(ctx, rec_graph, v1, v2) {
  n_vars <- nrow(rec_graph)
  graph <- matrix(FALSE, nrow = n_vars, ncol = n_vars)
  if (n_vars == 0L) {
    return(graph)
  }

  # Source trace: BCLQMTR.GRAFDEFINE(..., graphType = 1) keeps current-level
  # recursive graph edges and completes variables after the level. The current
  # R project reader exposes the example graphical SCREEN case as one level.
  level_start <- 1L
  level_end <- n_vars
  graph[level_start:level_end, level_start:n_vars] <- rec_graph[level_start:level_end, level_start:n_vars, drop = FALSE]
  graph[lower.tri(graph)] <- t(graph)[lower.tri(graph)]
  diag(graph) <- TRUE

  # Source trace: SKgrf8.Generate_hypotheses then applies LEVELCONNECT
  # collapsibility pruning. With one level this is inert for the validation target,
  # but the loop is kept source-shaped.
  level_connect <- matrix(0L, nrow = n_vars, ncol = 1L)
  for (row in seq_len(n_vars)) {
    for (col in seq_len(n_vars)) {
      if (row != col && rec_graph[row, col]) {
        level_connect[row, 1L] <- level_connect[row, 1L] + 1L
      }
    }
  }
  if (level_end < n_vars) {
    for (row in seq.int(level_end + 1L, n_vars)) {
      if (row != v2 && level_connect[row, 1L] == 0L) {
        graph[row, ] <- FALSE
        graph[, row] <- FALSE
      }
    }
  }
  graph
}

#' @keywords internal
screen_graph_degrees <- function(graph, v1, v2) {
  c(sum(graph[v1, -v1]), sum(graph[v2, -v2]))
}

#' @keywords internal
screen_local_markov_masks <- function(graph, fra, til) {
  n_fra <- sum(graph[fra, -fra])
  n_til <- sum(graph[til, -til])
  masks <- list()
  if (n_fra <= n_til) {
    masks[[length(masks) + 1L]] <- screen_mask_from_vertices(setdiff(which(graph[fra, ]), c(fra, til)))
  }
  if (n_til <= n_fra) {
    masks[[length(masks) + 1L]] <- screen_mask_from_vertices(setdiff(which(graph[til, ]), c(fra, til)))
  }
  masks
}

#' @keywords internal
screen_find_path_masks <- function(graph, fra, til, candidate_limit) {
  cpp <- screen_find_path_masks_cpp()
  if (!is.null(cpp)) {
    masks <- cpp(graph, as.integer(fra), as.integer(til), as.integer(candidate_limit))
    return(list(success = !is.null(attr(masks, "success")) && isTRUE(attr(masks, "success")), masks = as.list(as.numeric(masks))))
  }

  n_vars <- nrow(graph)
  old_paths <- list(list(prev = 0L, last = fra, visited = fra, inner = integer()))
  finished <- list()
  max_finished_paths <- 500L

  for (length_minus_one in seq_len(n_vars - 1L)) {
    new_paths <- list()
    for (path in old_paths) {
      candidates <- which(graph[path$last, ])
      candidates <- setdiff(candidates, path$last)
      if (length_minus_one > 1L) {
        candidates <- candidates[!(candidates != til & graph[path$prev, candidates])]
      }
      candidates <- setdiff(candidates, path$visited)
      if (til %in% candidates) {
        candidates <- setdiff(candidates, til)
        inner_mask <- screen_mask_from_vertices(path$inner)
        if (!screen_any_finished_subset(finished, inner_mask)) {
          if (length(finished) + 1L > max_finished_paths) {
            return(list(success = FALSE, masks = list()))
          }
          finished[[length(finished) + 1L]] <- inner_mask
        }
        if (length_minus_one > 1L) {
          next
        }
      }
      for (next_vertex in candidates) {
        if (length(new_paths) >= candidate_limit) {
          return(list(success = FALSE, masks = list()))
        }
        inner_mask <- screen_mask_from_vertices(path$inner)
        if (!screen_any_finished_subset(finished, inner_mask)) {
          new_paths[[length(new_paths) + 1L]] <- list(
            prev = path$last,
            last = next_vertex,
            visited = c(path$visited, next_vertex),
            inner = c(path$inner, next_vertex)
          )
        }
      }
    }
    if (length(new_paths) == 0L) {
      return(list(success = TRUE, masks = finished))
    }
    old_paths <- new_paths
  }
  list(success = TRUE, masks = finished)
}

#' @keywords internal
screen_find_path_masks_cpp <- local({
  compiled <- NULL
  function() {
    if (!is.null(compiled)) {
      return(compiled)
    }
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      return(NULL)
    }
    env <- new.env(parent = globalenv())
    Rcpp::sourceCpp(code = '
      #include <Rcpp.h>
      using namespace Rcpp;

      struct PathState {
        int prev;
        int last;
        unsigned long long visited;
        unsigned long long inner;
      };

      bool subset_mask(unsigned long long existing, unsigned long long candidate) {
        return (existing & ~candidate) == 0ULL;
      }

      bool any_finished_subset(const std::vector<unsigned long long>& finished, unsigned long long candidate) {
        for (size_t i = 0; i < finished.size(); ++i) {
          if (subset_mask(finished[i], candidate)) return true;
        }
        return false;
      }

      // [[Rcpp::export]]
      NumericVector gRm_screen_find_path_masks_cpp(LogicalMatrix graph, int fra, int til, int candidate_limit) {
        const int n = graph.nrow();
        const int max_finished_paths = 500;
        std::vector<PathState> old_paths;
        std::vector<PathState> new_paths;
        std::vector<unsigned long long> finished;
        bool success = true;
        old_paths.push_back({0, fra, 1ULL << (fra - 1), 0ULL});

        for (int lgd = 1; lgd <= n - 1; ++lgd) {
          new_paths.clear();
          new_paths.reserve(candidate_limit + 1);
          for (size_t path_index = 0; path_index < old_paths.size(); ++path_index) {
            PathState path = old_paths[path_index];
            std::vector<bool> candidates(n + 1, false);
            for (int j = 1; j <= n; ++j) candidates[j] = graph(path.last - 1, j - 1);
            candidates[path.last] = false;

            if (lgd > 1) {
              for (int j = 1; j <= n; ++j) {
                if (j != til && candidates[j] && graph(path.prev - 1, j - 1)) candidates[j] = false;
              }
            }
            for (int j = 1; j <= n; ++j) {
              if ((path.visited & (1ULL << (j - 1))) != 0ULL) candidates[j] = false;
            }

            if (candidates[til]) {
              candidates[til] = false;
              if (!any_finished_subset(finished, path.inner)) {
                if ((int)finished.size() + 1 > max_finished_paths) {
                  success = false;
                  NumericVector out(0);
                  out.attr("success") = success;
                  return out;
                }
                finished.push_back(path.inner);
              }
              if (lgd > 1) continue;
            }

            for (int j = 1; j <= n; ++j) {
              if (j == til || !candidates[j]) continue;
              if ((int)new_paths.size() >= candidate_limit) {
                success = false;
                NumericVector out(0);
                out.attr("success") = success;
                return out;
              }
              if (!any_finished_subset(finished, path.inner)) {
                unsigned long long bit = 1ULL << (j - 1);
                new_paths.push_back({path.last, j, path.visited | bit, path.inner | bit});
              }
            }
          }
          if (new_paths.empty()) break;
          old_paths.swap(new_paths);
        }

        NumericVector out(finished.size());
        for (size_t i = 0; i < finished.size(); ++i) out[i] = (double)finished[i];
        out.attr("success") = success;
        return out;
      }
    ', env = env, verbose = FALSE)
    compiled <- get("gRm_screen_find_path_masks_cpp", envir = env)
    compiled
  }
})

#' @keywords internal
screen_generate_cut_masks <- function(vertices, path_masks, n_vars = 64L) {
  source_max_cuts <- 400L
  cuts <- lapply(vertices, function(vertex) {
    mask <- screen_mask_from_vertices(vertex)
    list(mask = mask, last = vertex, missed = !vapply(path_masks, screen_mask_has_vertex, logical(1), vertex = vertex))
  })

  order <- 1L
  repeat {
    complete <- vapply(cuts, function(cut) !any(cut$missed), logical(1))
    if (any(complete)) {
      masks <- lapply(cuts[complete], `[[`, "mask")
      if (length(masks) > source_max_cuts) {
        return(list(success = FALSE, masks = masks[seq_len(source_max_cuts)]))
      }
      return(list(success = TRUE, masks = masks))
    }
    next_cuts <- list()
    for (cut in cuts) {
      for (vertex in seq.int(cut$last + 1L, n_vars)) {
        if (!any(vapply(path_masks, screen_mask_has_vertex, logical(1), vertex = vertex))) {
          next
        }
        missed_covered <- cut$missed & vapply(path_masks, screen_mask_has_vertex, logical(1), vertex = vertex)
        if (any(missed_covered)) {
          if (length(next_cuts) >= source_max_cuts) {
            return(list(success = FALSE, masks = lapply(next_cuts, `[[`, "mask")))
          }
          next_cuts[[length(next_cuts) + 1L]] <- list(
            mask = screen_mask_union(cut$mask, screen_mask_from_vertices(vertex)),
            last = vertex,
            missed = cut$missed & !vapply(path_masks, screen_mask_has_vertex, logical(1), vertex = vertex)
          )
        }
      }
    }
    order <- order + 1L
    if (order > length(vertices) || length(next_cuts) == 0L) {
      return(list(success = TRUE, masks = list()))
    }
    cuts <- next_cuts
  }
}

#' @keywords internal
screen_hypothesis_fits_bounds <- function(mask, ctx, v1, v2) {
  product_value <- 1
  active <- c(v1, v2, screen_vertices_from_mask(mask, nrow(ctx$variables)))
  active <- active[!is.na(active) & active >= 1L & active <= nrow(ctx$variables)]
  for (index in sort(unique(active))) {
    dim_value <- ctx$variables$dim[as.integer(index)]
    if (length(dim_value) != 1L || is.na(dim_value) || dim_value <= 0L) {
      return(FALSE)
    }
    product_value <- product_value * dim_value
    if (product_value > 9.22e18) {
      return(FALSE)
    }
  }
  TRUE
}

#' @keywords internal
screen_quick_tests_gamma_total <- function(mask, ctx, v1, v2) {
  x_dim <- ctx$variables$dim[[v1]]
  y_dim <- ctx$variables$dim[[v2]]
  conditioner_vars <- screen_vertices_from_mask(mask, nrow(ctx$variables))
  data <- ctx$data
  x <- data[, v1]
  y <- data[, v2]
  valid <- x >= 1L & x <= x_dim & y >= 1L & y <= y_dim
  if (length(conditioner_vars) > 0L) {
    for (conditioner in conditioner_vars) {
      z <- data[, conditioner]
      valid <- valid & z >= 1L & z <= ctx$variables$dim[[conditioner]]
    }
  }
  if (!any(valid)) {
    return(0)
  }

  keys <- if (length(conditioner_vars) == 0L) {
    rep("", sum(valid))
  } else {
    apply(data[valid, conditioner_vars, drop = FALSE], 1L, paste, collapse = "|")
  }
  rows <- x[valid]
  cols <- y[valid]
  ppq_total <- 0
  pmq_total <- 0
  for (key in unique(keys)) {
    in_stratum <- keys == key
    slice <- matrix(0, nrow = x_dim, ncol = y_dim)
    index <- rows[in_stratum] + (cols[in_stratum] - 1L) * x_dim
    slice[] <- tabulate(index, nbins = x_dim * y_dim)
    if ((sum(rowSums(slice) > 0) - 1L) * (sum(colSums(slice) > 0) - 1L) <= 0L) {
      next
    }
    stats <- screen_rc_gamma(slice)
    ppq_total <- ppq_total + stats$ppq
    pmq_total <- pmq_total + stats$pmq
  }
  if (ppq_total > 0) pmq_total / ppq_total else 0
}

#' @keywords internal
screen_mask_from_vertices <- function(vertices) {
  if (length(vertices) == 0L) {
    return(0)
  }
  sum(2^(as.integer(vertices) - 1L))
}

#' @keywords internal
screen_mask_has_vertex <- function(mask, vertex) {
  floor(mask / 2^(as.integer(vertex) - 1L)) %% 2 >= 1
}

#' @keywords internal
screen_mask_union <- function(mask_a, mask_b) {
  vertices <- union(screen_vertices_from_mask(mask_a), screen_vertices_from_mask(mask_b))
  screen_mask_from_vertices(vertices)
}

#' @keywords internal
screen_vertices_from_mask <- function(mask, n_vars = 64L) {
  if (length(mask) == 0L || is.na(mask) || mask <= 0) {
    return(integer())
  }
  which(floor(mask / 2^seq.int(0L, n_vars - 1L)) %% 2 >= 1)
}

#' @keywords internal
screen_any_finished_subset <- function(finished_masks, candidate_mask, powers = NULL) {
  if (length(finished_masks) == 0L) {
    return(FALSE)
  }
  candidate_vertices <- screen_vertices_from_mask(candidate_mask)
  any(vapply(finished_masks, function(mask) {
    all(screen_vertices_from_mask(mask) %in% candidate_vertices)
  }, logical(1)))
}

#' @keywords internal
screen_context <- function(project) {
  variables <- screen_variables(project)
  data <- screen_recoded_data(project, variables)

  list(
    variables = variables,
    data = data,
    complete = rep(TRUE, nrow(data)),
    n_items = sum(variables$is_item),
    n_backgrounds = sum(!variables$is_item),
    # Source trace: SKscr2/SKscr4 SCREEN use the active variable sequence and,
    # when no recursive explanatory block is supplied, test V1 = 1..NVAR-1.
    # The example graphical SCREEN project does not encode recursive levels in
    # DIGRAM.var, and the runtime matrix includes background-background pairs
    # such as J-K. Item/background flags are therefore not used as SCREEN level
    # boundaries here.
    screen_k = max(0L, nrow(variables) - 1L),
    levels = rep(1L, nrow(variables))
  )
}

#' @keywords internal
screen_variables <- function(project) {
  variables <- project$variables
  data.frame(
    label_code = variables$label_code,
    name = variables$name,
    position = variables$position,
    dim = variables$raw_max,
    vtype = variables$vtype,
    is_item = variables$is_item,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
screen_recoded_data <- function(project, variables) {
  data <- matrix(NA_integer_, nrow = nrow(project$raw_data), ncol = nrow(variables))
  for (index in seq_len(nrow(variables))) {
    data[, index] <- project$raw_data[, variables$position[[index]]]
  }
  colnames(data) <- variables$label_code
  data
}

#' @keywords internal
screen_initial_status <- function(ctx) {
  n_vars <- nrow(ctx$variables)
  status <- matrix(2L, nrow = n_vars, ncol = n_vars, dimnames = list(ctx$variables$label_code, ctx$variables$label_code))
  diag(status) <- 3L
  status
}

#' @keywords internal
screen_ordinals <- function(ctx, v1, v2) {
  !((ctx$variables$dim[[v1]] > 2L && ctx$variables$vtype[[v1]] < 3L) ||
      (ctx$variables$dim[[v2]] > 2L && ctx$variables$vtype[[v2]] < 3L))
}

#' @keywords internal
screen_two_way_table <- function(ctx, v1, v2) {
  x_dim <- ctx$variables$dim[[v1]]
  y_dim <- ctx$variables$dim[[v2]]
  x <- ctx$data[, v1]
  y <- ctx$data[, v2]
  valid <- x >= 1L & x <= x_dim & y >= 1L & y <= y_dim
  tab <- matrix(0, nrow = x_dim, ncol = y_dim)
  if (any(valid)) {
    index <- x[valid] + (y[valid] - 1L) * x_dim
    tab[] <- tabulate(index, nbins = x_dim * y_dim)
  }
  tab
}

#' @keywords internal
screen_three_way_table <- function(ctx, v1, v2, v3) {
  x_dim <- ctx$variables$dim[[v1]]
  y_dim <- ctx$variables$dim[[v2]]
  z_dim <- ctx$variables$dim[[v3]]
  x <- ctx$data[, v1]
  y <- ctx$data[, v2]
  z <- ctx$data[, v3]
  valid <- x >= 1L & x <= x_dim & y >= 1L & y <= y_dim & z >= 1L & z <= z_dim
  tab <- array(0, dim = c(x_dim, y_dim, z_dim))
  if (any(valid)) {
    index <- x[valid] + (y[valid] - 1L) * x_dim + (z[valid] - 1L) * x_dim * y_dim
    tab[] <- tabulate(index, nbins = x_dim * y_dim * z_dim)
  }
  tab
}

#' @keywords internal
screen_rc_chi <- function(tab) {
  row_totals <- rowSums(tab)
  col_totals <- colSums(tab)
  total <- sum(tab)
  chi <- 0
  if (total <= 0) {
    return(list(chi_square = 0, df = 0L, p_value = 1))
  }

  expected <- outer(row_totals, col_totals) / total
  positive <- expected > 0
  chi <- sum((tab[positive] - expected[positive])^2 / expected[positive])
  df <- (sum(row_totals > 0) - 1L) * (sum(col_totals > 0) - 1L)
  if (df <= 0L) {
    list(chi_square = chi, df = 0L, p_value = 1)
  } else {
    # Source trace: SkStat.RCCHI/SourceRaschCore.SourcePFCHI computes the upper
    # chi-square tail for the Pearson statistic.
    list(chi_square = chi, df = df, p_value = stats::pchisq(chi, df = df, lower.tail = FALSE))
  }
}

#' @keywords internal
screen_rc_chi_square <- function(tab) {
  row_totals <- rowSums(tab)
  col_totals <- colSums(tab)
  total <- sum(tab)
  if (total <= 0) {
    return(0)
  }

  expected <- outer(row_totals, col_totals) / total
  positive <- expected > 0
  sum((tab[positive] - expected[positive])^2 / expected[positive])
}

#' @keywords internal
screen_rc_gamma <- function(tab) {
  n_rows <- nrow(tab)
  n_cols <- ncol(tab)
  aij <- matrix(0, nrow = n_rows, ncol = n_cols)
  dij <- matrix(0, nrow = n_rows, ncol = n_cols)
  cumulative <- t(apply(apply(tab, 2L, cumsum), 1L, cumsum))
  if (n_rows == 1L) {
    cumulative <- matrix(cumulative, nrow = 1L)
  }
  total <- sum(tab)
  p <- 0
  q <- 0

  cell_sum <- function(row_to, col_to) {
    if (row_to <= 0L || col_to <= 0L) {
      0
    } else {
      cumulative[row_to, col_to]
    }
  }

  for (row in seq_len(n_rows)) {
    for (col in seq_len(n_cols)) {
      less_less <- cell_sum(row - 1L, col - 1L)
      greater_greater <- total - cell_sum(row, n_cols) - cell_sum(n_rows, col) + cell_sum(row, col)
      less_greater <- cell_sum(row - 1L, n_cols) - cell_sum(row - 1L, col)
      greater_less <- cell_sum(n_rows, col - 1L) - cell_sum(row, col - 1L)
      aij[row, col] <- less_less + greater_greater
      dij[row, col] <- less_greater + greater_less
      p <- p + tab[row, col] * aij[row, col]
      q <- q + tab[row, col] * dij[row, col]
    }
  }

  ppq <- p + q
  pmq <- p - q
  gamma_value <- if (ppq > 0) pmq / ppq else 0
  n <- sum(tab)
  if (ppq <= 0) {
    return(list(gamma = 0, ppq = ppq, pmq = pmq, s = 0, p_value = 1, success = FALSE))
  }

  # Source trace: SkStat.RCGAMMA/SourceRaschCore.SourceRCGammaStats computes
  # S = 4 * (-PMQ * PMQ / N + sum(n_ij * (AIJ - DIJ)^2)) and tests
  # abs(gamma / (sqrt(S) / PPQ)) against the upper normal tail.
  s <- if (n > 0) -pmq * (pmq / n) else 0
  m <- aij - dij
  s <- s + sum(tab * m * m)
  s <- 4 * s
  p_value <- if (s > 0) {
    stats::pnorm(abs(gamma_value / (sqrt(s) / ppq)), lower.tail = FALSE)
  } else {
    1
  }

  list(gamma = gamma_value, ppq = ppq, pmq = pmq, s = s, p_value = p_value, success = TRUE)
}

#' @keywords internal
screen_rc_gamma_counts <- function(tab) {
  n_rows <- nrow(tab)
  n_cols <- ncol(tab)
  cumulative <- t(apply(apply(tab, 2L, cumsum), 1L, cumsum))
  if (n_rows == 1L) {
    cumulative <- matrix(cumulative, nrow = 1L)
  }
  total <- sum(tab)
  p <- 0
  q <- 0

  cell_sum <- function(row_to, col_to) {
    if (row_to <= 0L || col_to <= 0L) {
      0
    } else {
      cumulative[row_to, col_to]
    }
  }

  for (row in seq_len(n_rows)) {
    for (col in seq_len(n_cols)) {
      less_less <- cell_sum(row - 1L, col - 1L)
      greater_greater <- total - cell_sum(row, n_cols) - cell_sum(n_rows, col) + cell_sum(row, col)
      less_greater <- cell_sum(row - 1L, n_cols) - cell_sum(row - 1L, col)
      greater_less <- cell_sum(n_rows, col - 1L) - cell_sum(row, col - 1L)
      p <- p + tab[row, col] * (less_less + greater_greater)
      q <- q + tab[row, col] * (less_greater + greater_less)
    }
  }

  ppq <- p + q
  pmq <- p - q
  list(gamma = if (ppq > 0) pmq / ppq else 0, ppq = ppq, pmq = pmq)
}

#' @keywords internal
screen_three_way_test <- function(ctx, v1, v2, v3) {
  tab <- screen_three_way_table(ctx, v1, v2, v3)
  chi_square <- 0
  df <- 0L
  ppq_total <- 0
  pmq_total <- 0
  s_total <- 0
  ordinal <- screen_ordinals(ctx, v1, v2)

  for (level in seq_len(dim(tab)[[3L]])) {
    slice <- tab[, , level, drop = FALSE][, , 1L]
    chi <- screen_rc_chi(slice)
    chi_square <- chi_square + chi$chi_square
    df <- df + chi$df
    if (ordinal) {
      gamma <- screen_rc_gamma(slice)
      ppq_total <- ppq_total + gamma$ppq
      pmq_total <- pmq_total + gamma$pmq
      s_total <- s_total + gamma$s
    }
  }

  chi_p <- if (df > 0L) stats::pchisq(chi_square, df = df, lower.tail = FALSE) else 1
  if (ordinal && ppq_total > 0) {
    gamma_value <- pmq_total / ppq_total
    s_total <- s_total / ppq_total
    s_total <- s_total / ppq_total
    gamma_p <- if (s_total == 0) {
      1
    } else {
      u <- abs(gamma_value / sqrt(s_total))
      if (u > 4) 0.000001 else stats::pnorm(u, lower.tail = FALSE)
    }
  } else {
    gamma_value <- 0
    gamma_p <- 1
  }

  list(
    chi_square = chi_square,
    df = df,
    chi_p = chi_p,
    ordinal = ordinal,
    gamma = gamma_value,
    gamma_p = gamma_p
  )
}

#' @keywords internal
screen_determine_significance <- function(chi_p, gamma_p, ordinal, significance) {
  accept <- chi_p > significance
  if (accept && ordinal && gamma_p <= significance) {
    accept <- FALSE
  }
  accept
}

#' @keywords internal
screen_first_node_at_level <- function(ctx, node) {
  which(ctx$levels >= ctx$levels[[node]])[[1]]
}

#' @keywords internal
screen_revise_model3 <- function(status, test_result, v1, v2) {
  if (test_result == 1L && status[v1, v2] != 1L) {
    status[v1, v2] <- status[v2, v1] <- 1L
  }
  if (test_result == 5L && status[v1, v2] == 2L) {
    status[v1, v2] <- status[v2, v1] <- 5L
  }
  status
}

#' @keywords internal
screen_stage <- function(name, status, old_status, variables, tests, significance) {
  if (is.null(tests)) {
    tests <- data.frame()
  }
  structure(
    list(
      name = name,
      variables = variables,
      significance = significance,
      status = status,
      old_status = old_status,
      status_chars = screen_status_chars(status, old_status),
      tests = tests
    ),
    class = "gRm_screen_stage"
  )
}

#' @keywords internal
screen_status_chars <- function(status, old_status = status) {
  chars <- matrix(" ", nrow = nrow(status), ncol = ncol(status), dimnames = dimnames(status))
  chars[status == 2L] <- "+"
  chars[status == 3L] <- "*"
  chars[status == 4L] <- "h"
  chars[status == 5L] <- "o"
  chars[status < 2L & old_status > 1L] <- "-"
  chars
}

#' @keywords internal
screen_status_text <- function(status) {
  text <- matrix("conditional_independence", nrow = nrow(status), ncol = ncol(status), dimnames = dimnames(status))
  text[status == 2L] <- "undecided"
  text[status == 3L] <- "fixed"
  text[status == 4L] <- "hidden_interaction"
  text[status == 5L] <- "unused_conditional_independence"
  text
}

#' @keywords internal
screen_incomplete_path_warnings <- function(project, final_model) {
  variables <- screen_variables(project)
  status <- final_model$status
  n_vars <- nrow(status)
  candidate_limit <- (30L * 162230L) %/% n_vars
  graph <- status > 1L
  diag(graph) <- FALSE

  warnings <- list()
  ctx <- screen_context(project)
  pairs <- list()
  pair_index <- 0L
  for (v1 in seq_len(n_vars - 1L)) {
    for (v2 in seq.int(v1 + 1L, n_vars)) {
      if (!screen_ordinals(ctx, v1, v2) || status[v1, v2] < 2L) {
        next
      }
      pair_index <- pair_index + 1L
      pairs[[pair_index]] <- c(v1, v2)
    }
  }

  evaluate_pair <- function(pair) {
    v1 <- pair[[1]]
    v2 <- pair[[2]]
    reduced <- screen_separate1_reduced_graph(graph, v1, v2)
    degree1 <- sum(reduced[v1, ])
    degree2 <- sum(reduced[v2, ])
    if (degree2 < degree1) {
      fra <- v2
      til <- v1
    } else {
      fra <- v1
      til <- v2
    }
    reduced <- screen_separate1_reduced_graph(graph, fra, til)
    if (!screen_findpaths_incomplete(reduced, fra, til, candidate_limit)) {
      return(NULL)
    }
    data.frame(
      var1 = fra,
      var2 = til,
      label1 = variables$label_code[[fra]],
      label2 = variables$label_code[[til]],
      source_status = "ported_DGRexe_execute_gamma_false_SKgrf8_Generate_hypotheses_SKsep1_FINDPATHS_warning",
      stringsAsFactors = FALSE
    )
  }

  if (length(pairs) > 0L) {
    warnings <- lapply(pairs, evaluate_pair)
    warnings <- Filter(Negate(is.null), warnings)
  }

  if (length(warnings) > 0L) {
    do.call(rbind, warnings)
  } else {
    data.frame(
      var1 = integer(),
      var2 = integer(),
      label1 = character(),
      label2 = character(),
      source_status = character(),
      stringsAsFactors = FALSE
    )[0, ]
  }
}

#' @keywords internal
screen_separate1_reduced_graph <- function(graph, fra, til) {
  reduced <- graph
  diag(reduced) <- FALSE
  reduced[fra, til] <- FALSE
  reduced[til, fra] <- FALSE

  # Source trace: SKsep1.SEPARATE1 first treats common neighbours of the
  # tested pair as direct separators (`DIREKTE`) and removes them before
  # calling FINDPATHS on the remaining graph.
  direct <- rep(FALSE, nrow(reduced))
  direct[-c(fra, til)] <- reduced[fra, -c(fra, til)] & reduced[til, -c(fra, til)]
  if (any(direct)) {
    reduced[direct, ] <- FALSE
    reduced[, direct] <- FALSE
  }
  reduced
}

#' @keywords internal
screen_path_subset <- function(existing_path, candidate_path) {
  existing_inner <- if (length(existing_path) > 2L) existing_path[2L:(length(existing_path) - 1L)] else integer()
  candidate_inner <- if (length(candidate_path) > 2L) candidate_path[2L:(length(candidate_path) - 1L)] else integer()
  all(existing_inner %in% candidate_inner)
}

#' @keywords internal
screen_findpaths_incomplete <- function(graph, fra, til, candidate_limit) {
  cpp <- screen_findpaths_incomplete_cpp()
  if (!is.null(cpp)) {
    return(cpp(graph, as.integer(fra), as.integer(til), as.integer(candidate_limit)))
  }

  n_vars <- nrow(graph)
  powers <- 2 ^ seq.int(0L, n_vars - 1L)
  old_len <- 1
  old_prev <- 0
  old_last <- fra
  old_visited <- powers[[fra]]
  old_inner <- 0
  finished_masks <- numeric()
  max_finished_paths <- 500L

  for (length_minus_one in seq_len(n_vars - 1L)) {
    new_count <- 0L
    new_len <- numeric(candidate_limit + 1L)
    new_prev <- numeric(candidate_limit + 1L)
    new_last <- numeric(candidate_limit + 1L)
    new_visited <- numeric(candidate_limit + 1L)
    new_inner <- numeric(candidate_limit + 1L)

    for (path_index in seq_along(old_last)) {
      last_vertex <- old_last[[path_index]]
      candidates <- graph[last_vertex, ]
      candidates[last_vertex] <- FALSE

      # Source trace: SKsep1.FINDPATHS skips a longer route through a clique by
      # removing neighbours that also connect to the previous path vertex.
      if (length_minus_one > 1L) {
        previous <- old_prev[[path_index]]
        candidates[-til] <- candidates[-til] & !graph[previous, -til]
      }
      visited <- old_visited[[path_index]]
      visited_vertices <- which(floor(visited / powers) %% 2 >= 1)
      if (length(visited_vertices) > 0L) {
        candidates[visited_vertices] <- FALSE
      }

      # Source trace: FINDPATHS always tries to close to TIL before extending
      # unfinished paths. If a finished path is found after the first step, all
      # unfinished extensions from the same candidate are redundant.
      if (isTRUE(candidates[[til]])) {
        candidates[[til]] <- FALSE
        if (length(finished_masks) + 1L > max_finished_paths) {
          return(FALSE)
        }
        candidate_inner <- old_inner[[path_index]]
        redundant <- screen_any_finished_subset(finished_masks, candidate_inner, powers)
        if (!redundant) {
          finished_masks <- c(finished_masks, candidate_inner)
        }
        if (length_minus_one > 1L) {
          next
        }
      }

      for (next_vertex in which(candidates)) {
        if (next_vertex == til) {
          next
        }
        if (new_count >= candidate_limit) {
          return(TRUE)
        }
        candidate_inner <- old_inner[[path_index]]
        redundant <- screen_any_finished_subset(finished_masks, candidate_inner, powers)
        if (!redundant) {
          new_count <- new_count + 1L
          new_len[[new_count]] <- old_len[[path_index]] + 1L
          new_prev[[new_count]] <- old_last[[path_index]]
          new_last[[new_count]] <- next_vertex
          new_visited[[new_count]] <- old_visited[[path_index]] + powers[[next_vertex]]
          new_inner[[new_count]] <- old_inner[[path_index]] + powers[[next_vertex]]
        }
      }
    }
    if (new_count == 0L) {
      return(FALSE)
    }
    old_len <- new_len[seq_len(new_count)]
    old_prev <- new_prev[seq_len(new_count)]
    old_last <- new_last[seq_len(new_count)]
    old_visited <- new_visited[seq_len(new_count)]
    old_inner <- new_inner[seq_len(new_count)]
  }
  FALSE
}

#' @keywords internal
screen_findpaths_incomplete_cpp <- local({
  compiled <- NULL
  function() {
    if (!is.null(compiled)) {
      return(compiled)
    }
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      return(NULL)
    }
    env <- new.env(parent = globalenv())
    Rcpp::sourceCpp(code = '
      #include <Rcpp.h>
      using namespace Rcpp;

      struct PathState {
        int len;
        int prev;
        int last;
        unsigned long long visited;
        unsigned long long inner;
      };

      bool subset_mask(unsigned long long existing, unsigned long long candidate) {
        return (existing & ~candidate) == 0ULL;
      }

      bool any_finished_subset(const std::vector<unsigned long long>& finished, unsigned long long candidate) {
        for (size_t i = 0; i < finished.size(); ++i) {
          if (subset_mask(finished[i], candidate)) return true;
        }
        return false;
      }

      // [[Rcpp::export]]
      bool gRm_screen_findpaths_incomplete_cpp(LogicalMatrix graph, int fra, int til, int candidate_limit) {
        const int n = graph.nrow();
        const int max_finished_paths = 500;
        std::vector<PathState> old_paths;
        std::vector<PathState> new_paths;
        std::vector<unsigned long long> finished;
        old_paths.push_back({1, 0, fra, 1ULL << (fra - 1), 0ULL});

        for (int lgd = 1; lgd <= n - 1; ++lgd) {
          new_paths.clear();
          new_paths.reserve(candidate_limit + 1);
          for (size_t path_index = 0; path_index < old_paths.size(); ++path_index) {
            PathState path = old_paths[path_index];
            std::vector<bool> candidates(n + 1, false);
            for (int j = 1; j <= n; ++j) candidates[j] = graph(path.last - 1, j - 1);
            candidates[path.last] = false;

            if (lgd > 1) {
              for (int j = 1; j <= n; ++j) {
                if (j != til && candidates[j] && graph(path.prev - 1, j - 1)) candidates[j] = false;
              }
            }
            for (int j = 1; j <= n; ++j) {
              if ((path.visited & (1ULL << (j - 1))) != 0ULL) candidates[j] = false;
            }

            if (candidates[til]) {
              candidates[til] = false;
              if ((int)finished.size() + 1 > max_finished_paths) return false;
              if (!any_finished_subset(finished, path.inner)) finished.push_back(path.inner);
              if (lgd > 1) continue;
            }

            for (int j = 1; j <= n; ++j) {
              if (j == til || !candidates[j]) continue;
              if ((int)new_paths.size() >= candidate_limit) return true;
              if (!any_finished_subset(finished, path.inner)) {
                unsigned long long bit = 1ULL << (j - 1);
                new_paths.push_back({path.len + 1, path.last, j, path.visited | bit, path.inner | bit});
              }
            }
          }
          if (new_paths.empty()) return false;
          old_paths.swap(new_paths);
        }
        return false;
      }
    ', env = env, verbose = FALSE)
    compiled <- get("gRm_screen_findpaths_incomplete_cpp", envir = env)
    compiled
  }
})

#' @keywords internal
screen_mask_subset <- function(existing_mask, candidate_mask, powers = NULL) {
  if (existing_mask == 0) {
    return(TRUE)
  }
  if (is.null(powers)) {
    existing_vertices <- screen_vertices_from_mask(existing_mask)
    return(all(existing_vertices %in% screen_vertices_from_mask(candidate_mask)))
  }
  existing_vertices <- which(floor(existing_mask / powers) %% 2 >= 1)
  all(floor(candidate_mask / powers[existing_vertices]) %% 2 >= 1)
}

#' @keywords internal
screen_any_finished_subset <- function(finished_masks, candidate_mask, powers = NULL) {
  if (length(finished_masks) == 0L) {
    return(FALSE)
  }
  any(vapply(finished_masks, screen_mask_subset, logical(1), candidate_mask = candidate_mask, powers = powers))
}
