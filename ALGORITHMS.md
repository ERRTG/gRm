# gRm Algorithm Pseudocode

This file describes the algorithms currently implemented in `gRm/R/` in a
rigorous pseudocode style. It is intentionally equivalent to the current R
implementation, not a cleaned-up replacement model. Where the notation is
clearer than the code, the equivalence is stated explicitly.

Production R code must continue to implement these algorithms directly from the
Pascal source or the source-faithful Pascal harness. Pascal output may be used
as a test oracle, but not as a production input.

## Source Trace

Primary implemented R files:

- `R/api-constructors.R`
- `R/api-ari.R`
- `R/api-model-spec.R`
- `R/api-fit.R`
- `R/api-likelihood-comparison.R`
- `R/api-screen.R`
- `R/api-results.R`
- `R/api-summary.R`
- `R/api-summary-print-docs.R`
- `R/api-summary-tables.R`
- `R/api-table-helpers.R`
- `R/api-model-graph.R`
- `R/api-model-plot.R`
- `R/project_input.R`
- `R/digram_project_io.R`
- `R/source_bundle.R`
- `R/ari_values.R`
- `R/ari_plot.R`
- `R/rasch_base_fit.R`
- `R/item_parameters_values.R`
- `R/item_fits_values.R`
- `R/screen_j_values.R`
- `R/dif_tests_values.R`
- `R/exact_command_state.R`
- `R/local_independence_values.R`
- `R/global_homogeneity_values.R`
- `R/source_score_groups.R`
- `R/exo_select_values.R`
- `R/gamma_values.R`
- `R/gRm-package.R`
- `R/gllrm_candidate_fit.R`
- `R/gllrm_context.R`
- `R/gllrm_components.R`
- `R/gllrm_fit.R`
- `R/gllrm_probability_cache.R`
- `R/gllrm_values.R`
- `R/internal-utils.R`
- `R/source_candidate_map.R`
- `R/source_gamma_stats.R`
- `R/source_statistics.R`

Primary source traces:

- `pascal_harness/SourceRaschCore.pas::LoadCounts`
- `pascal_harness/SourceRaschCore.pas::InitializeRaschFit`
- `pascal_harness/SourceRaschCore.pas::BuildGammaExcludingItem`
- `pascal_harness/SourceRaschCore.pas::CalculateRaschExpectedItems`
- `pascal_harness/SourceRaschCore.pas::CalculateRaschUpdateRatiosFromScore`
- `pascal_harness/SourceRaschCore.pas::AdjustItemGammasSourceScale`
- `pascal_harness/SourceRaschCore.pas::EstimateBaseRaschParameters`
- `pascal_harness/SourceRaschCore.pas::SourceItemThresholdFromGamma`
- `pascal_harness/SourceRaschCore.pas::SourceItemLocationFromGamma`
- `pascal_harness/SourceRaschCore.pas::CalculateTrueScoreFromGamma`
- `pascal_harness/SourceRaschCore.pas::InitializePersonParameter2`
- `pascal_harness/SourceRaschCore.pas::EstimatePersonParameter`
- `pascal_harness/SourceRaschCore.pas::SourceItemDifficultyFromGamma`
- `pascal_harness/SourceRaschCore.pas::ScoreProbabilitiesForLogTheta`
- `pascal_harness/SourceRaschCore.pas::CalculateSourceItemTarget`
- `pascal_harness/SourceRaschCore.pas::SourceItemInformation`
- `pascal_harness/SourceRaschCore.pas::EmitGLLRMOutputRows` output(4)
- `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits`
- `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits.Count_Observed`
- `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits.CalculateMeans`
- `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits.CalculateInAndOutfits`
- `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits.CalculateInAndOutfits.CalculateOutfit`
- `source/PAS_skunits/skbias15.pas::Calculate_residuals_and_item_fits.CalculateInAndOutfits.CalculateInfit`
- `source/PAS_skunits/skbias12a.pas::IncompleteItemfits`
- `source/PAS_scd/DGRirtD.pas::TargetSlut`
- `source/PAS_skunits/skbias15.pas::Calculate_Ari`
- `source/PAS_skunits/skbias15.pas::Calculate_Ari.Count_Observed`
- `source/PAS_skunits/skbias15.pas::Calculate_Ari.CalculateExpectedValues`
- `source/PAS_skunits/skbias15.pas::Calculate_Ari.CalculateMeans`
- `source/PAS_skunits/skbias15.pas::Calculate_Ari.SaveARI`
- `pascal_harness/item_parameters_report/example_ITEM_PARAMETERS_REPORT.pas::CalculateObservedScoreRange`
- `pascal_harness/item_parameters_report/example_ITEM_PARAMETERS_REPORT.pas::CalculateSourceNParameters`

## ARI Item Score Curves

Source trace:

- `source/PAS_scd/DGRirtD.pas::TargetSlut`, ITA 18 branch;
- `source/PAS_skunits/skbias15.pas::Calculate_Ari`;
- `Calculate_Ari.Count_Observed`;
- `Calculate_Ari.CalculateExpectedValues`;
- `Calculate_Ari.CalculateMeans`;
- `Calculate_Ari.SaveARI`.

```text
algorithm ari_values(fit)
  bundle <- fit source-shaped estimation bundle
  items <- bundle item metadata
  max_score <- sum(item_dimension[item] - 1 over all items)
  first_score <- 1
  last_score <- max_score - 1
  global_item_max <- max(item_dimension[item]) - 1

  initialize observed[score, item, item_score] to zero for all scores,
    items, and global item-score columns

  for each source-valid row in bundle data:
    if row total score is outside first_score..last_score, skip it
    if any item score is invalid for its item, skip it
    row_count <- source row count if present, otherwise 1
    for each item:
      observed[score, item, item_score] += row_count
    end for
  end for

  if fit is a base Rasch fit:
    for each item and score:
      compute expected item-score probabilities from
        item_gamma[item, item_score] * gamma_without_item[score - item_score]
        normalized by the full conditional gamma denominator
    end for
  else:
    for each score and DIF exogeneous cell:
      compute GLLRM item probabilities using the fitted LD/DIF component
        probability machinery
      weight each cell by the source score/exogeneous cell count
    end for
    normalize weighted GLLRM probabilities within each total score
  end if

  for item in source item order:
    for score from first_score to last_score:
      n <- sum observed[score, item, ]
      if n is zero, do not emit a row

      Obs0..ObsK <- observed[score, item, ] / n
      ObsMean <- sum item_score * ObsProbability over supported item scores
      raw_ObsVar <- sum item_score^2 * ObsProbability - ObsMean^2
      if n > 1:
        ObsVar <- ((n - 1) / n) * raw_ObsVar
      else:
        ObsVar <- 0

      Exp0..ExpK <- fitted expected item-score probabilities, with unsupported
        global category columns filled by zero
      ExpMean <- fitted expected mean over supported item scores
      ExpVar <- fitted expected variance over supported item scores
      if ExpVar > 0:
        z <- sqrt(n) * (ObsMean - ExpMean) / sqrt(ExpVar)
      else:
        z <- 0

      emit ItemNo, Item, Score, n, Obs0..ObsK, ObsMean, ObsVar,
        Exp0..ExpK, ExpMean, ExpVar, z
    end for
  end for
end algorithm
```

The returned `ari()` object is an R data frame. Unlike DIGRAM, the package does
not write `Ari_dot.csv` or `Ari_comma.csv`; those files are validation oracles
only.

## ARI Plot Data

Source reference:

- `ARIplot.sas`, which reads `Ari_dot.csv`, collapses raw total-score rows into
  class intervals, aggregates item means by item and interval, and draws one
  panel per item.

The package implementation does not read `Ari_dot.csv`. It consumes the
source-shaped `gRm_ari` table returned by `ari()`.

```text
algorithm ari_score_intervals(ari_table, class_size)
  require ari_table is a non-empty gRm_ari table with ItemNo, Score, and n
  require class_size is a positive whole number

  first_item <- minimum ItemNo in ari_table
  score_distribution <- rows with ItemNo == first_item, columns Score and n
  require score_distribution has no duplicated Score values
  sort score_distribution by Score ascending

  interval <- 1
  running_frequency <- 0
  for each score row in sorted score_distribution:
    if running_frequency >= class_size:
      interval <- interval + 1
      running_frequency <- 0
    end if

    running_frequency <- running_frequency + n
    assign current score to interval
  end for

  if there is more than one interval and the final interval total n is less
     than class_size:
    merge the final interval into the previous interval
  end if

  if there is only one interval, keep it as interval 1 even when total n is
    below class_size. This is the R-side degenerate guard for a usable plot.

  return Score, interval
end algorithm

algorithm ari_plot_data(ari_table, class_size, confidence)
  require ari_table has ItemNo, Item, Score, n, ObsMean, ExpMean, ExpVar
  require confidence is strictly between 0 and 1

  intervals <- ari_score_intervals(ari_table, class_size)
  join intervals back to all ari_table rows by Score
  multiplier <- qnorm((1 + confidence) / 2)

  for each ItemNo, Item, interval in source item and interval order:
    rows <- ARI rows for this item and interval
    N <- sum(rows.n)
    O <- sum(rows.n * rows.ObsMean) / N
    E <- sum(rows.n * rows.ExpMean) / N
    V <- sum(rows.n * rows.ExpVar) / N
    lower <- E - multiplier * sqrt(V / N)
    upper <- E + multiplier * sqrt(V / N)
    emit ItemNo, Item, interval, O, N, E, V, lower, upper
  end for
end algorithm

algorithm plot.gRm_ari(ari_table)
  plot_data <- ari_plot_data(ari_table)
  optionally filter plot_data by ItemNo or Item after interval construction
  validate requested facet rows and columns

  draw a ggplot with:
    one expected confidence-band ribbon,
    one observed mean line,
    optional dashed expected mean line,
    facet panels by ItemNo with item names as strip titles,
    integer x-axis breaks 1..maximum interval,
    x label "class interval",
    y label "Mean item score",
    ggplot2::theme_minimal()
end algorithm
```

## Retired Non-J SCREEN Scope

The package no longer implements DIGRAM's older dialog/non-J SCREEN
workflow. In the current package, `screen()` means the validated SCREEN J item
screening workflow described below. Historical repository notes about the
retired workflow may remain outside the package, but they are not part of this
R implementation or its installed API.

## SCREEN J Item Screening

Source trace:

- `source/PAS_scd/DIGRAM1f.pas`, command 3 `SCREEN`, parameter `J`;
- `source/PAS_scd/SKbias7.pas::Item_screening`;
- `source/PAS_scd/SKbias7.pas::inexpensive_rosenberg`;
- `source/PAS_scd/SKbias7.pas::inexpensive_itembias1`;
- `source/PAS_scd/SKbias3.pas::XYZ_bias_ANALYSE`;
- `source/PAS_skunits/SKbias7.pas::number_of_Tjur_problems`,
  `Adjusted_number_of_Tjur_problems`, and `stepwise_elimination`;
- `source/PAS_skunits/DIGRAM1f.pas`, command 3 `SCREEN`, parameters `I` and
  `J`, for the final positive-local-dependence model filter;
- `source/PAS_skunits/SKVars.pas`, version 3.37, for the rule that negative
  local dependence is reported but not included in the screen model;
- `SkStat.RCGAMMA`, `SkStat.RCCHI`, `PNORMAL`, `PFCHI`, and
  `BenjaminiHochberg`.

```text
algorithm screen_j_values(project, exact = false, nsim = 1000, seed = 9)
  read source-coded item responses, exogeneous responses, item dimensions, and
  source variable types from the parsed DIGRAM project

  if exact is true:
    seed the R random number generator with seed while preserving the caller's
      previous RNG state

  complete_items[row] <- all item responses are valid source categories
  total_score[row] <- sum(item_response[row, item] - 1 over all items)
  max_score <- sum(item_dimension[item] - 1 over all items)

  for each row item X:
    rest_score_X[row] <- total_score[row] - (item_response[row, X] - 1)

    # Source: inexpensive_rosenberg, item-restscore LongRCgamma.
    build table item_response[X] by rest_score_X + 1 over complete_items rows
    store marginal item-restscore gamma and asymptotic p from RCGAMMA

    for each other item Y:
      # Source: inexpensive_rosenberg, marginal item-item RCGAMMA.
      build table item_response[X] by item_response[Y] over complete_items rows
      store marginal gamma and asymptotic p from RCGAMMA

      # Source: inexpensive_rosenberg and XYZ_bias_ANALYSE.
      build strata tables item_response[X] by item_response[Y] within
        rest_score_X + 1 over complete_items rows
      for each stratum:
        compute RCGAMMA, accumulating PPQ, PMQ, and S
      partial_gamma <- sum(PMQ) / sum(PPQ)
      partial_variance <- sum(S) / sum(PPQ)^2
      if exact is false:
        partial_p <- 2 * upper_normal_tail(abs(partial_gamma /
          sqrt(partial_variance)))
      else:
        repeat nsim times:
          for each score stratum:
            generate one random two-way table with the same margins using the
              source SKrandom.GENTAB1 conditional hypergeometric traversal
            compute RCGAMMA on the generated table
          accumulate simulated PPQ and PMQ over strata
          simulated_gamma <- simulated_PMQ / simulated_PPQ, or 0 if
            simulated_PPQ is zero
          count abs(simulated_gamma) >= abs(partial_gamma)
        partial_p <- exceedance_count / nsim
      store directed partial item-item result X -> Y
    end for
  end for

  for each exogeneous variable E:
    use_gamma <- source type is ordinal or E is dichotomous

    # Source: inexpensive_itembias1, score-exogeneous LongRCgamma/LongRCchi.
    build table total_score + 1 by exogeneous_response[E] over complete item
      rows with valid E
    store score-exogeneous gamma or chi-square result

    for each item X:
      # Source detail: this marginal branch reads the item and exogeneous
      # variables directly from TABDATA. It does not require other items to be
      # complete.
      build table item_response[X] by exogeneous_response[E] over rows where
        X and E are valid
      store marginal DIF gamma or chi-square result

      # Source detail: partial DIF uses Z := SCORES from the precomputed total
      # score array and excludes only Z = 0 and Z = maxscore. It does not use
      # Score\X, despite the local-independence heading above the DIF block.
      build strata tables item_response[X] by exogeneous_response[E] within
        total_score + 1 over complete item rows with valid E and
        0 < total_score < max_score
      if use_gamma:
        accumulate partial gamma as in XYZ_bias_ANALYSE
        if exact is true:
          compute the Monte Carlo exact gamma p-value by the same GENTAB1
          repeated-table procedure used for directed partial item-item tests
      else:
        accumulate Pearson chi-square and degrees of freedom over strata and
        compute PFCHI(df_total, chi_total)
        if exact is true:
          repeat nsim times:
            for each score stratum, generate one GENTAB1 random table with
              fixed margins and compute its Pearson chi-square statistic
            count simulated_chi_total >= observed_chi_total
          partial_p <- exceedance_count / nsim
      store directed item-exogeneous partial DIF result
    end for
  end for

  p_values <- all directed partial item-item p-values plus all partial DIF
    p-values
  for alpha in {0.05, 0.01, 0.001}:
    sort p_values ascending
    return the source Benjamini-Hochberg critical boundary alpha * i / m for
      the largest accepted index i
  end for
end algorithm
```

The 5% Benjamini-Hochberg boundary is not itself the final local-dependence
model rule. DIGRAM first classifies the two directed tests for every unordered
item pair, combines their concordance counts into a weighted partial gamma,
and applies a staged greedy procedure. The resulting rows are provisional
SCREEN evidence. A separate command-level filter then determines which of
those rows are eligible to become GLLRM terms.

For an unordered pair \(\{i,j\}\), write \(\gamma_{ij}\) and \(p_{ij}\) for
the partial gamma and p-value obtained when item \(i\) is the row item and item
\(j\) is the other item. Let \(c_{.05}\) be the global 5% boundary returned by
the preceding Benjamini-Hochberg calculation. The source treats a zero gamma
as non-positive.

```text
algorithm screen_j_directed_evidence_counts(gamma, p, i, j, c_05)
  positive <- 0
  negative <- 0

  for direction (a, b) in {(i, j), (j, i)}:
    if p[a, b] <= c_05:
      if gamma[a, b] > 0:
        positive <- positive + 1
      else:
        negative <- negative + 1
      end if
    end if
  end for

  return (positive, negative)
end algorithm
```

Thus `positive = 2` means that both directed tests cross the common 5% FDR
boundary and both directed gammas are positive. `positive = 1` means exactly
one directed test supplies significant positive evidence; it does not imply
that the other direction is non-significant, because that direction may
instead supply significant negative evidence. The analogous interpretation
holds for `negative = 2` and `negative = 1`.

The pairwise ordering statistic is DIGRAM's weighted partial gamma (WPG), not
the unweighted mean of the two directed gammas. If `PMQ[a,b]` and `PPQ[a,b]`
are the directed concordance-minus-discordance and comparable-pair totals
accumulated by item screening, then

\[
  \operatorname{WPG}_{ij}
  =
  \frac{\operatorname{PMQ}_{ij}+\operatorname{PMQ}_{ji}}
       {\operatorname{PPQ}_{ij}+\operatorname{PPQ}_{ji}}.
\]

The same symmetric value is stored at positions \((i,j)\) and \((j,i)\). A
pair with zero combined `PPQ` has no usable WPG and cannot improve the
zero-valued positive or negative search sentinel, so it is not chosen by the
source stepwise procedure.

```text
algorithm screen_j_greedy_ld_evidence(gamma, p, PMQ, PPQ, c_05)
  temp_p <- p
  provisional_matrix <- symmetric all-false item-by-item matrix
  provisional_rows <- empty ordered table

  for each unordered pair {i,j}, i < j:
    original_count[i,j] <-
      screen_j_directed_evidence_counts(gamma, p, i, j, c_05)
    if PPQ[i,j] + PPQ[j,i] > 0:
      WPG[i,j] <- (PMQ[i,j] + PMQ[j,i]) /
                  (PPQ[i,j] + PPQ[j,i])
      WPG[j,i] <- WPG[i,j]
    else:
      WPG[i,j] and WPG[j,i] <- unavailable
    end if
  end for

  stages <- [(positive, 2), (positive, 1),
             (negative, 2), (negative, 1)]

  for (sign, required_count) in stages, in this exact order:
    repeat:
      best_pair <- none
      best_value <- 0

      scan pairs {i,j} in increasing i and then increasing j:
        adjusted_count <-
          screen_j_directed_evidence_counts(gamma, temp_p, i, j, c_05)

        eligible <-
          original_count[i,j][sign] = required_count and
          adjusted_count[sign] > 0 and
          WPG[i,j] is available

        if sign is positive and eligible and WPG[i,j] > best_value:
          best_pair <- {i,j}
          best_value <- WPG[i,j]
        else if sign is negative and eligible and WPG[i,j] < best_value:
          best_pair <- {i,j}
          best_value <- WPG[i,j]
        end if
      end scan

      if best_pair is none:
        leave this stage
      end if

      let best_pair be {i0,j0}
      provisional_matrix[i0,j0] <- true
      provisional_matrix[j0,i0] <- true
      append (i0, j0, WPG[i0,j0],
              stage = sign concatenated with required_count)
        to provisional_rows

      # Source-faithful directed suppression. The original gamma and p-value
      # matrices remain unchanged; only these two outgoing temporary rows are
      # made ineligible for later choices.
      temp_p[i0, every item] <- 2.0
      temp_p[j0, every item] <- 2.0
    end repeat
  end for

  return (provisional_matrix, provisional_rows, WPG)
end algorithm
```

Comparisons with `best_value` are strict. Consequently, a positive stage only
selects a pair with WPG greater than zero, a negative stage only selects a pair
with WPG less than zero, and an exact WPG tie is retained for the first pair in
source scan order. Suppression is directional because `temp_p` stores directed
tests: after selecting \(\{i_0,j_0\}\), evidence directed *from* either selected
item is disabled. A competing pair that shares one of those items can still be
eligible if its remaining direction, from the unselected item, crosses the FDR
boundary. This is the behavior of `stepwise_elimination`; replacing it with
undirected whole-item deletion would be a different algorithm.

The stage label describes why a row was chosen during the greedy evidence
procedure. It is not the final model-inclusion decision. After item screening
returns, the command-level code in `DIGRAM1f.pas` applies the following filter
to every provisional pair.

```text
algorithm screen_j_final_ld_model_filter(provisional_rows, gamma)
  included_rows <- empty ordered table
  excluded_negative_rows <- empty ordered table

  for row in provisional_rows:
    i <- row.row
    j <- row.col

    if gamma[i,j] + gamma[j,i] > 0:
      row.included <- true
      append row to included_rows
    else:
      row.included <- false
      append row to excluded_negative_rows
    end if
  end for

  final_model_matrix <- symmetric matrix containing only included_rows
  return (included_rows, excluded_negative_rows, final_model_matrix)
end algorithm
```

This last test deliberately uses the unweighted sum of the two directed
partial gammas, not the sign of WPG and not the stage label. DIGRAM reports the
excluded rows as negative local dependence followed by "Not included in the
model". The version 3.37 source history records the same invariant. Negative
evidence must therefore remain available for diagnostics while being absent
from the screen-derived GLLRM.

At the public API boundary, SCREEN evidence and model terms are distinct data
products. The following contract prevents a stage name from being mistaken for
an exclusion status and prevents a reported negative pair from silently
entering a fitted model.

```text
algorithm screen_j_public_summary_and_model_contract(screen_values)
  diagnostics.local_dependence <- all tested unordered item pairs with their
    directed partial results, WPG, directed-gamma sum, final inclusion marker,
    and a decision that retains provisional negative evidence

  evidence_rows <- every row chosen by screen_j_greedy_ld_evidence
  require each evidence row retains its stage in
    {positive2, positive1, negative2, negative1}
  require each evidence row has an explicit included decision from
    screen_j_final_ld_model_filter

  summary(screen).local_dependence <- all pair diagnostics, including a
    decision that distinguishes final included pairs from provisional
    negative evidence

  summary(screen).selected_ld <- evidence_rows where included is true
  summary(screen).selected <- selected_ld together with final selected DIF rows
  summary(screen).model_terms.ld <- selected_ld
  summary(screen).model_terms.dif <- final selected DIF rows

  model_terms(screen).ld <- summary(screen).model_terms.ld
  model_terms(screen).dif <- summary(screen).model_terms.dif
  gllrm(screen), when ld and dif are not explicitly overridden, uses exactly
    model_terms(screen).ld and model_terms(screen).dif

  require no row with included = false occurs in model_terms(screen).ld
  require excluded negative evidence remains visible in the SCREEN summary
end algorithm
```

The current R package exposes these values through `screen()` and
`summary.gRm_screen`; it no longer implements a DIGRAM text renderer for the
SCREEN J report.

Implementation note: exact and repeated SCREEN J may use C++ slice entry points
for the source-shaped conditional chi-square/gamma simulations
(`gRm_screen_j_exact_chi_gamma_slices`, `gRm_screen_j_exact_chi_slices`,
`gRm_screen_j_exact_gamma_slices`, `gRm_screen_j_conditional_bias_test`, and
`gRm_screen_j_item_pair_conditional_exact`). Those native routes are optional
accelerators for the pseudocode above. The R implementation keeps parity probes
for the conditional native route, and falls back to the R source-shaped path when
native use is disabled or a probe does not pass. The older generic exact-kernel
entry point is not part of the current package boundary.

```text
algorithm screen_two_way_matrix(project, significance)
  start from the current status graph
  set diagonals to fixed status 3
  set each non-fixed admissible candidate to status 1

  for V1 in source order up to K
    for V2 in source order from V1 + 1 to NVAR
      if the pair is fixed or inadmissible, skip it
      build the observed two-way table from complete positive values
      compute expected counts from row and column marginals
      p_chi <- source RCCHI p-value
      if both variables are ordinal
        p_gamma <- source RCGAMMA p-value
      else
        p_gamma <- 1
      accepted <- p_chi > significance
      if accepted and ordinal and p_gamma <= significance
        accepted <- false
      if not accepted
        set symmetric status to 2
      end if
    end for
  end for

  return status matrix plus chi-square/gamma diagnostics
end algorithm
```

```text
algorithm screen_hidden_matrix(project, two_way, significance)
  hidden <- two_way status matrix
  gate <- hidden

  Note: in the Pascal source, MODEL_H := MODEL_2 is a pointer assignment. The
  hidden matrix and the matrix used by later candidate-gate checks are the same
  mutable object. This aliasing is intentional source behavior and creates
  cascade hidden links.

  for V1 in source order up to K
    for V2 in source order from V1 + 1 to NVAR
      if gate[V1,V2] != 1, skip the pair
      for V3 from LEVEL[LEVNR[V1] - 1] + 1 to NVAR
        require V3 different from V1 and V2
        require gate[V1,V3] > 1 and gate[V2,V3] > 1
        run the source THREEWAY conditional test for pair V1,V2 | V3
        apply DETERMINE_SIGNIFICANCE to PCHI3[1] and PGAMMA3[1]
        if not accepted
          set hidden[V1,V2] and hidden[V2,V1] to status 4
          gate is updated by the same assignment because hidden and gate alias
        end if
      end for
    end for
  end for

  return hidden matrix and conditioner diagnostics
end algorithm
```

```text
algorithm screen_final_model(project, two_way, hidden, significance)
  final <- hidden status matrix

  for each source-ordered triple satisfying the MODEL_H candidate rules
    mark fixed/suppressed pair tests as not performed
    run source THREEWAY for the performed pair tests
    for each performed test
      p_min <- min(PCHI3, PGAMMA3) for ordinal pairs, otherwise PCHI3
    end for
    largest_non_significant <- performed test with largest p_min
    for each performed test
      if p_min <= significance
        result <- 2  // reject, keep edge
      else if test is largest_non_significant
        result <- 1  // accept, remove edge
      else
        result <- 5  // not used, keep as unused conditional independence
      end if
      if result == 1 and final[pair] != 1, set final[pair] to 1
      if result == 5 and final[pair] == 2, set final[pair] to 5
    end for
  end for

  return final status matrix for SHOW_REVISED_STATUSMATRIX
end algorithm
```

The Pascal and R implementations of this source-shaped loop agree and
reproduce the supplied example final SCREEN matrix. This depends on preserving
the source `MODEL_3 := MODEL_H` pointer assignment: final revisions update the
same matrix consulted by later final-stage gates.

```text
algorithm screen_incomplete_path_warnings(project, final)
  graph <- MAKE_GRAPH(final status), where graph[i,j] = 1 when status[i,j] > 1
  candidate_limit <- (30 * 162230) div NVAR
  warnings <- empty list

  for i in 1..NVAR-1
    for j in i+1..NVAR
      if pair is not ordinal or final status[i,j] < 2
        skip
      end if

      reduced <- separate1_reduced_graph(graph, i, j)
      if degree(j, reduced) < degree(i, reduced)
        fra <- j; til <- i
      else
        fra <- i; til <- j
      end if

      reduced <- separate1_reduced_graph(graph, fra, til)
      if findpaths_incomplete(reduced, fra, til, candidate_limit)
        append warning (fra, til)
      end if
    end for
  end for

  return warnings in source traversal order
end algorithm
```

```text
algorithm separate1_reduced_graph(graph, fra, til)
  reduced <- graph
  remove diagonal
  remove tested edge fra--til
  direct <- vertices connected to both fra and til, excluding fra and til
  remove every direct vertex from reduced
  return reduced
end algorithm
```

```text
algorithm findpaths_incomplete(graph, fra, til, candidate_limit)
  old_paths <- one unfinished path with last = fra, visited = {fra}, inner = {}
  finished_inner_sets <- empty list
  max_finished_paths <- 500

  for length_minus_one in 1..NVAR-1
    new_paths <- empty list
    for each path in old_paths
      candidates <- neighbours(path.last)
      remove path.last and all visited vertices from candidates
      if length_minus_one > 1
        remove candidate j != til when j is also adjacent to path.previous
      end if

      if til is a candidate
        remove til from candidates
        if adding this finished path would exceed max_finished_paths
          return false  // source exits without the incomplete-path warning
        end if
        if no existing finished inner set is a subset of path.inner
          add path.inner to finished_inner_sets
        end if
        if length_minus_one > 1
          continue to next old path
        end if
      end if

      for candidate j in source vertex order
        if count(new_paths) >= candidate_limit
          return true
        end if
        if no existing finished inner set is a subset of path.inner
          append unfinished path with previous = path.last,
                 last = j, visited = path.visited union {j},
                 inner = path.inner union {j}
        end if
      end for
    end for
    if new_paths is empty
      return false
    end if
    old_paths <- new_paths
  end for

  return false
end algorithm
```

Implementation note: the R production function owns the value computation and
does not read Pascal TSV output or `screen.txt`. For speed, the inner
`findpaths_incomplete` state expansion uses an Rcpp helper that is a direct
state-record translation of the pseudocode above; the R fallback remains in
the file for auditability.

```text
algorithm screen_partial_gamma_matrix(project, final)
  initialize every off-diagonal value to source sentinel 999

  for each unordered pair i,j
    if both variables are ordinal and final status < 2
      value[i,j] <- value[j,i] <- 0
    else if not both ordinal and final status > 1
      value[i,j] <- value[j,i] <- 888
    else if both variables are ordinal
      hypotheses <- Generate_hypotheses("test <label_i><label_j>")
      reduced <- CheckHyp(hypotheses)
      if hypotheses exist
        results <- quick_tests(hypotheses)
        value[i,j] <- value[j,i] <- mean(results[, 5])
      end if
    end if
  end for

  keep source sentinels compatible with DGRVARS.Show_gammavalues semantics:
  999 for diagonal/unavailable ordinal cells, 888 for nonordinal retained
  graph cells, 0 for accepted ordinal nonedges, and numeric partial gamma for
  tested graph-conditioned ordinal cells
end algorithm
```

Implementation note: native R now ports the
`Generate_hypotheses`/`CheckHyp`/`quick_tests` route for the final graph and
tests it against the Pascal structured oracle. The supplied example partial-gamma
values must still not be used as target constants; the documented `b-d`
contradiction is explicitly OPEN / UNRESOLVED until a runtime command log or
saved graph/status state explains the appended block.

## Notation

The pseudocode uses zero-based score categories, matching the R representation:

- `i = 1..I` indexes items.
- `x = 0..m_i` indexes the score category of item `i`.
- `s = 0..M` indexes the total score.
- `m_i = raw_max_i - 1`.
- `M = sum_i m_i`.
- `Y[r,i]` is the zero-based item score for row `r`.
- `B[r,b]` is the one-based background value for row `r`, or `-1` if invalid.
- `score[r] = sum_i Y[r,i]` for complete valid rows, otherwise `-1`.
- `status[r] = 1` means row `r` is complete for all items and backgrounds and
  lies inside the source estimation score window.
- `g[i,x]` is the fitted multiplicative item score parameter, called
  `item_gamma` in R.
- `N_s[s]` is the observed count of valid rows with total score `s`.
- `N_ix[i,x]` is the observed count of valid rows where item `i` has score `x`.
- `E_ix[i,x]` is the expected count under the current Rasch fit.
- `G_without_i[s]` is the score-gamma convolution for all items except `i`.

All vectors that are indexed by scores are conceptually indexed from zero even
when R stores them at position `score + 1`.

## Data Entry and Project Construction

The current R package builds its internal DIGRAM project representation from
either a user data frame, a raw CSV file, or a DIGRAM import-file directory. The
old `DIGRAM.var`/`DIGRAM.dat` reader is no longer the package data-entry
surface.

```text
algorithm gRm(data, items, exogenous, id, item_levels,
              exogenous_levels, score_cuts, name)
  require data is a data frame
  require at least one declared item column
  if id is supplied, require it names a column in data

  project <- build_gRm_internal_project(
    data,
    items = declared item columns,
    exo = declared exogenous columns,
    item labels = automatic DIGRAM aliases,
    exogenous labels = automatic DIGRAM aliases,
    item levels = item_levels,
    exogenous levels = exogenous_levels
  )

  attach source_data, import metadata, and source_trace
  score_groups <- normalize_gRm_score_cuts(score_cuts, project)
  return gRm_analysis(project, data, id, score_groups, name, call)
end algorithm
```

```text
algorithm read_digram_project(path, items, exogenous, id, score_cuts, name)
  normalize path and require it exists
  require explicit item names

  project <- read_digram_files(path, items, exogenous, id, name)
  score_groups <- normalize_gRm_score_cuts(score_cuts, project)

  return gRm_analysis(project, project.source_data, project.import.idvar,
                         score_groups, basename(path), call)
end algorithm
```

```text
algorithm read_digram_files(input_dir, items, exo, idvar, name)
  require input_dir contains name.imp
  imp <- read_digram_imp(name.imp)
  # DIGRAM.imp supplies only the source directory and project prefix here.
  # Later command-file lines are intentionally ignored; the R package does
  # not parse or execute .cmd files.
  require imp.path contains imp.project_name.csv and imp.project_name.imv
  metadata <- read_digram_imv(imp.project_name.imv)
  require every declared item/exogenous variable appears exactly once in metadata
  reorder metadata to declared item order followed by declared exogenous order

  data <- read imp.project_name.csv with literal column names
  if idvar is NULL, use first CSV column
  require idvar appears in data

  project <- build_gRm_internal_project(
    data,
    items = declared item columns,
    exo = declared exogenous columns,
    item labels/maxima = metadata rows for items,
    exogenous labels/maxima = metadata rows for exogenous variables,
    item levels = NULL,
    exogenous levels = NULL,
    paths = import-file paths
  )

  attach source_data, import metadata, and source_trace
  return project
end algorithm
```

```text
algorithm read_digram_csv(csv_path, items, exo, idvar, output_dir,
                          save_digram_files)
  read csv_path with literal column names
  if idvar is NULL, use first CSV column
  require idvar appears in data

  project <- build_gRm_internal_project with observed source levels and
             automatic DIGRAM aliases
  attach source_data, import metadata, and source_trace

  if save_digram_files:
    require output_dir
    write DIGRAM.csv with user id column and one-based encoded
      project.raw_data item/exogenous categories; write internal missing
      sentinel -999 as blank cells
    write DIGRAM.imp with folder, project name, and two source "-" placeholder
      lines; no command file is written or parsed by the package
    write DIGRAM.imv with label code, variable name, and one-based category
      labels for the supported categories
    record written file paths in project.import
  end if

  return project
end algorithm
```

```text
algorithm build_gRm_internal_project(data, items, exo, labels, maxima,
                                     item_levels, exogenous_levels, paths)
  require item and exogenous declarations are unique and present in data
  require no more than the 50 currently defined DIGRAM aliases

  assign label codes from supplied metadata or automatic aliases
  if maxima are supplied from DIGRAM.imv:
    resolve levels as 1..maxima for each variable
    validate supplied maxima cover all observed non-missing values
  else:
    resolve item/exogenous levels from explicit level vectors or sorted
      observed non-missing values
  end if

  raw_data <- integer matrix over item columns followed by exogenous columns
  encode each observed value as its one-based DIGRAM category index
  replace missing values in raw_data by -999

  variables <- table with label_code, position, raw_max, vtype = 3, name,
               and is_item
  return gRm_project with variables, items, backgrounds, raw_data, and paths
end algorithm
```

```text
algorithm normalize_gRm_score_cuts(score_cuts, project)
  if score_cuts is "auto":
    cuts <- upper score cuts from items_select_values(project).score_groups
    require cuts define at least two usable source score groups
    return cuts

  if score_cuts is numeric:
    require non-missing integer-like values
    require strictly increasing values
    require every cut is within 0..sum(raw_max_i - 1)
    return integer cuts

  otherwise error
end algorithm
```

## Public Modeling API

The exported package surface is intentionally the statistical modeling API:
`gRm`, `read_digram_project`, `gllrm`, `fit`, `screen`,
`score_effects`, `item_parameters`, `item_fit`, `local_dependence`, `dif`, and
`global_homogeneity`. The package no longer exports report-generation,
validation, `tidy()`, `glance()`, `details()`, `detail_names()`, or provenance
helper APIs.

```text
algorithm gllrm(project, ld, dif)
  if project is a gRm_screen:
    selected <- model_terms(project)
    if ld is NULL, use selected local-dependence terms; otherwise parse ld
    if dif is NULL, use selected DIF terms; otherwise parse dif
    return gRm_model with screen provenance
  end if

  analysis <- as_gRm_analysis(project)
  ld_terms <- parse_ld_formula(ld, analysis)
  dif_terms <- parse_dif_formula(dif, analysis)
  model_type <- "gllrm" if either term table is nonempty, otherwise "rasch"
  return gRm_model(analysis, ld_terms, dif_terms, model_type, call)
end algorithm
```

```text
algorithm parse_gRm_interaction_formula(formula, items, exogenous, type)
  if formula is NULL or one-sided formula ~0:
    return an empty term table of the requested type
  require formula is one-sided
  require every formula term is exactly a two-variable ":" interaction
  reject formula operators other than "+" between ":" terms

  for each term label:
    if type is "ld":
      require both variables are item names
      require the two item names differ
      order the item pair by source item order
      append canonical row (type, item1, item2, source, status)
    else if type is "dif":
      require exactly one item name and exactly one exogenous name
      append canonical row (type, item, exogenous, source, status)
    end if
  end for

  require no duplicate canonical terms
  return canonical term table
end algorithm
```

```text
algorithm screen(analysis, inference, nsim, seed, jobs)
  inference <- asymptotic, exact, or repeated
  exact <- inference is exact or repeated
  repeated <- inference is repeated

  values <- screen_j_values(analysis.project, exact, repeated, nsim, seed)
  out <- gRm_screen containing analysis, values, inference controls,
         source_trace, and call
  out.terms <- model_terms(out)
  return out
end algorithm
```

```text
algorithm fit(model, max_step, max_delta)
  require model is a gRm_model

  if model has any LD or DIF terms:
    gllrm_fit <- fit_gllrm(model, max_step, max_delta)
    values <- gllrm_values(gllrm_fit, model)
    return gRm_fit(model, gllrm_fit.bundle, gllrm_fit.state/context,
                      values, GLLRM convergence metadata)
  else:
    bundle <- build_item_parameters_bundle(model.analysis.project)
    base_fit <- fit_rasch_base(bundle, max_step, max_delta)
    values <- item_parameters_values(base_fit, bundle)
    return gRm_fit(model, bundle, base_fit, values,
                      base convergence metadata)
  end if
end algorithm
```

```text
algorithm score_effects(analysis, inference, nsim, seed, jobs)
  analysis <- as_public_gRm_analysis(analysis)
  exact <- inference is exact or repeated
  repeated <- inference is repeated
  values <- exo_select_values(analysis.project, exact, repeated, nsim, seed)
  return gRm_score_effects result with inference metadata
end algorithm
```

```text
algorithm item_parameters(fit)
  fit <- as_public_gRm_fit(fit)
  return gRm_item_parameters result using fit.values
end algorithm

algorithm item_fit(fit, include_extended)
  fit <- as_public_gRm_fit(fit)
  if fit is an GLLRM fit:
    values <- item_fits_values(fit, include_extended)
  else:
    values <- item_fits_values(fit.project, include_extended)
  end if
  return gRm_item_fit result
end algorithm

algorithm local_dependence(fit, max_step, max_delta, jobs)
  fit <- as_public_gRm_fit(fit)
  if fit is an GLLRM fit:
    values <- gllrm_local_independence_values(fit, max_step, max_delta, jobs)
  else:
    values <- local_independence_values(fit.project, max_step, max_delta, jobs)
  end if
  return gRm_local_dependence result
end algorithm

algorithm dif(fit, max_step, max_delta, jobs)
  fit <- as_public_gRm_fit(fit)
  if fit is an GLLRM fit:
    values <- gllrm_dif_tests_values(fit, max_step, max_delta, jobs)
  else:
    values <- dif_tests_values(fit.project, max_step, max_delta, jobs)
  end if
  return gRm_dif result
end algorithm

algorithm global_homogeneity(fit, score_cuts, max_step, max_delta, jobs)
  fit <- as_public_gRm_fit(fit)
  score_cuts <- supplied score_cuts, or analysis score groups, normalized to
                source-compatible increasing integer cuts
  values <- global_homogeneity_values(fit or fit.project, score_cuts,
                                      max_step, max_delta)
  return gRm_global_homogeneity result
end algorithm
```

## Public Summary Output

Public output is exposed through `summary()` methods and compact print methods.
The package does not provide an internal or exported `details()`,
`detail_names()`, `tidy()`, or `glance()` output API.

```text
algorithm summary.gRm_analysis(object, which)
  validate which is "data"
  return summary.gRm with one data table:
    row count, item count, exogenous count, item names, exogenous names,
    and score-group cuts
end algorithm

algorithm summary.gRm_model(object, which)
  validate which is drawn from model, ld, dif
  return summary.gRm with requested tables:
    model counts/type, canonical LD term table, canonical DIF term table
end algorithm

algorithm summary.gRm_fit(object, which)
  validate which is drawn from fit, parameters, terms
  return summary.gRm with requested tables:
    fit likelihood/convergence/counts,
    item-parameter item statistics,
    canonical current GLLRM terms
end algorithm

algorithm summary for result objects
  score_effects: expose selected score effects, all score-effect tests, and BH
  item_parameters: expose coefficients/items, thresholds, and fit summary
  item_fit: expose item fit tests/items and BH thresholds
  local_dependence: expose selected LD, all LD tests, and BH threshold
  dif: expose selected DIF, all no-DIF tests, included-DIF tests, and BH threshold
  global_homogeneity: expose summary, score groups, and item rows

  for every exposed table with a converged column:
    count rows where converged is FALSE
    if any rows did not converge:
      add a remarks table row naming the table, the number of non-converged
      candidate fits, and stop-reason counts when stop_reason is present
    end if
end algorithm
```

## Result Object Metadata Fields

The package does not expose separate check-orchestration, provenance, warning,
or validation helper APIs. Each public workflow returns one typed object, and
users inspect numeric values through `summary()` plus compact `print()`
methods. Result objects may still carry internal list fields such as
`source_trace`, `warnings`, and `unmodeled` for package tests and source
maintenance, but those fields are not a public accessor layer.

## Global Homogeneity Report

The global homogeneity implementation is a source-shaped R port of the DIGRAM
path in `DGRirtD.pas` that loops over score cuts, refits the model inside each
score group, calls the `skbias15.pas` item-margin branch, and prints a
conditional likelihood-ratio summary. Production R computes from the parsed
DIGRAM project and fitted Rasch/GLLRM state; Pascal and supplied runtime reports
remain test oracles only.

Inputs:

- Parsed DIGRAM project `P`
- Score cuts `C[1..G]`, such as `30, 87` for the example runtime example. Public
  report/workflow wrappers use explicit cuts when supplied; otherwise they
  derive the source item-selection default score groups and pass their upper
  bounds as `C`.
- Base Rasch fit on all valid complete rows
- Source controls `LeastScore`, `LargestScore`, `MaxTotalScore`

```text
algorithm global_homogeneity_score_groups(bundle, score_cuts)
  previous_cut <- LeastScore - 1

  for each cut c in score_cuts
    fra <- previous_cut + 1
    if fra < LeastScore then fra <- LeastScore

    til <- c
    if til > LargestScore then til <- LargestScore

    if fra <= til and fra < MaxTotalScore and til > 0
      emit score group with lower bound fra and upper bound til
    end if

    previous_cut <- c
  end for
end algorithm
```

```text
algorithm global_homogeneity_values(project, score_cuts)
  bundle <- build_item_parameters_bundle(project)

  fit full Rasch model on all valid complete rows
  L_full <- source negative log likelihood of the full model
  N_parameters <- source item-parameter count from full observed item margins

  groups <- global_homogeneity_score_groups(bundle, score_cuts)
  L_groups <- 0

  for each score group g with bounds fra_g..til_g
    group_bundle <- copy bundle
    mark rows invalid unless status = 1 and fra_g <= total score <= til_g

    fit subgroup Rasch model on group_bundle
    L_g <- source negative log likelihood of subgroup model
    L_groups <- L_groups + L_g

    for each item i
      observed_counts[k] <- subgroup observed count for item score k
      expected_counts[k] <- expected item count under the full-model item
                            gammas and the subgroup score distribution

      n_i <- sum_k observed_counts[k]
      observed_mean_i <- sum_k k * observed_counts[k] / n_i
      expected_mean_i <- sum_k k * expected_counts[k] / sum_k expected_counts[k]

      full_gamma <- skbias12-style full score gamma array
      without_item_gamma <- skbias12-style score gamma array excluding item i
      ExpectedItemMargTables <- zero score-by-item-score table
      for each total score s in the subgroup score distribution
        for each item score k
          e <- gamma_i[k] * without_item_gamma[s-k] * N_s / full_gamma[s]
          ExpectedItemMargTables[s, k] <- ExpectedItemMargTables[s, k] + e
        end for
      end for

      expected_variance_i <- NA
      residual_i <- NA
      marker_i <- NA
    end for
  end for

  CLR <- 2 * abs(L_full - L_groups)
  df <- (number of groups - 1) * N_parameters
  p <- PFCHI(df, CLR)

  return score groups, item mean rows, and summary
end algorithm
```

The Pascal source oracle used in tests follows the same pseudo-code and emits
the intermediate denominator terms as structured rows. It is deliberately a
test/reference surface: production R calls `global_homogeneity_values()` and
does not call Pascal or read DIGRAM report text. The historical executable's
printed residual variance materialization has not been recovered
source-faithfully, so the base Rasch global-homogeneity rows report the
expected-variance, residual, and marker cells as `NA` and validate only the
source-backed counts, means, CLR, df, and p values.

Residual boundary note:

- `global_homogeneity_item_mean_rows()` follows the source report branch that
  calls `skbias15.pas::Calculate_residuals_and_item_fits` for global
  homogeneity.
- The helper computes expected item-margin tables needed for source-backed item
  means, but it deliberately does not expose the unsupported residual variance
  path.
- Reporting expected variance, residual, or marker values for base global
  homogeneity would require a fresh source trace of the hidden runtime residual
  materialization. The package must not fill these cells from hard-coded report
  values, empirical correction factors, or the removed native expected-summary
  helper.

For GLLRM fits, `gllrm_global_homogeneity_values()` repeats the
same score-group refit structure with included LD/DIF terms present. GLLRM item
mean rows compute expected item-margin tables from the GLLRM probability cache
and expose source-shaped expected variances from those tables. The residual and
marker cells remain `NA` with `*_runtime_source_backed = FALSE`, because the
runtime residual denominator is still not source-backed.

The current R package returns structured homogeneity rows through
`global_homogeneity()` and `summary.gRm_global_homogeneity`; it does not
assemble the historical DIGRAM text report.

## Compact ItemFits Report

The compact ItemFits implementation is a source-shaped R port of the report family
rooted in `skbias15.pas::Calculate_item_fits` and its compact `output(123)`
branch, with shared algorithms inherited from the corresponding `skbias14.pas`
and `skbias12a.pas` routines. Production R computes values from the parsed
DIGRAM project and the fitted base Rasch state; Pascal is only a test oracle.
The example item-restscore gamma anchors and the compact Motiva01/Motiva32
outfit/infit anchors match the supplied compact report.

The outfit/infit calculation intentionally preserves a mathematically strange
DIGRAM source convention. When exogenous variables are selected,
`skbias15.pas::Count_Observed` applies `GET_EXOGENE` before filling
`ItemMargTables`, so observed item-score proportions can be based on
valid-background rows. Later, `CalculateInAndOutfits` weights those proportions
with the global complete-item `ScoreDistribution`. If a complete score group has
no valid-background observed item margin, the source leaves that score's outfit
contribution at zero but still includes the complete score weight in the final
outfit denominator; infit receives no variance weight for that score.

The supplied example `ItemFits-extended.txt` does not contain the source line
stating that incomplete item responses were included. Therefore this report
does not synthesize `NincompleteRecs` merely from missing item rows in the bundle.
The `IncompleteItemFits` branch remains a source-traced branch for future runtime
examples where DIGRAM explicitly populates it.

Visible implementation gap: the R ItemFits implementation is not complete until
the source `IncompleteItemFits` branch has been implemented and validated
against a DIGRAM runtime example that actually exercises it. See
`gRm/INCOMPLETE_ITEM_FITS_GAP.md`.

Inputs:

- Parsed DIGRAM project `P`
- Base Rasch fit with item gammas `gamma[i,k]`
- Complete item-response rows `Y[p,i]` and item total scores `S[p]`
- Valid-background complete rows used by `Count_Observed` when exogenous
  variables are selected
- Source score window `LeastScore..LargestScore`
- Runtime-provided incomplete item-fit records, used by
  `skbias12a.pas::IncompleteItemfits` only when the source runtime has
  populated `NincompleteRecs`

```text
algorithm item_fits_values(project)
  bundle <- build_item_parameters_bundle(project)
  fit <- fit_rasch_base(bundle)

  for each item i
    build G_without_i, the score polynomial excluding item i

    for each total score s
      denom <- sum_k gamma[i,k] * G_without_i[s-k]
      Pr_i,s(k) <- gamma[i,k] * G_without_i[s-k] / denom
      mu_i,s <- sum_k k * Pr_i,s(k)
      var_i,s <- sum_k (k - mu_i,s)^2 * Pr_i,s(k)
      fourth_i,s <- sum_k (k - mu_i,s)^4 * Pr_i,s(k)
    end for

    initialize outfit and infit accumulators
  end for

  score_weight[s] <- count complete item-response rows with total score s
  observed_rows <- complete item-response rows
  if exogenous variables are selected
    observed_rows <- observed_rows with source-valid exogenous values
  end if

  for each score s in the source score window
    for each item i
      observed_count[k] <- count observed_rows with S[p] = s and Y[p,i] = k
      observed_total <- sum_k observed_count[k]

      if score_weight[s] = 0
        continue
      end if

      if observed_total = 0
        n_used[i] <- n_used[i] + score_weight[s]
        continue
      end if

      observed_frequency[k] <- observed_count[k] / observed_total
      z2[k] <- (k - mu_i,s)^2 / var_i,s
      r2[k] <- (k - mu_i,s)^2

      outfit_score <- sum_k z2[k] * observed_frequency[k]
      outfit_expected <- sum_k z2[k] * Pr_i,s(k)
      outfit_variance <- sum_k z2[k]^2 * Pr_i,s(k) - outfit_expected^2

      infit_score <- sum_k r2[k] * observed_frequency[k]
      infit_expected <- sum_k r2[k] * Pr_i,s(k)
      infit_variance <- sum_k r2[k]^2 * Pr_i,s(k) - infit_expected^2

      outfit_sum[i] <- outfit_sum[i] + score_weight[s] * outfit_score
      outfit_mean[i] <- outfit_mean[i] + score_weight[s] * outfit_expected
      outfit_var[i] <- outfit_var[i] + score_weight[s] * outfit_variance

      infit_sum[i] <- infit_sum[i] + score_weight[s] * infit_score
      infit_mean[i] <- infit_mean[i] + score_weight[s] * infit_expected
      infit_var[i] <- infit_var[i] + score_weight[s] * infit_variance
      infit_weight[i] <- infit_weight[i] + score_weight[s] * var_i,s
      n_used[i] <- n_used[i] + score_weight[s]
    end for
  end for

  for each runtime incomplete person p
    for each observed item i in p
      use the observed item subset and partial observed score
      compute Pr_i,p(k) from the subset score polynomial
      add the source `IncompleteItemfits` outfit/infit contribution
    end for
  end for

  for each item i
    outfit <- outfit_sum[i] / n_used[i]
    outfit_sd <- sqrt(outfit_var[i] / n_used[i]^2)
    p_outfit <- 2 * PNORMAL(abs((outfit - outfit_mean[i] / n_used[i]) / outfit_sd))

    infit <- infit_sum[i] / infit_weight[i]
    infit_sd <- sqrt(infit_var[i] / infit_weight[i]^2)
    p_infit <- 2 * PNORMAL(abs((infit - infit_mean[i] / infit_weight[i]) / infit_sd))
  end for
end algorithm
```

```text
algorithm item_restscore_gamma(bundle, fit)
  complete_items[p] <- all item responses are valid
  total_score[p] <- sum_i Y[p,i] over complete-item rows
  low_extreme <- count complete_items with total_score < LeastScore
  high_extreme <- count complete_items with total_score > LargestScore
                  or total_score = highest_possible_score

  observed_rows <- complete_items
  if exogenous variables are selected
    observed_rows <- observed_rows with source-valid exogenous values
  end if
  item_margin_n[s] <- count observed_rows with total score s

  for each item i
    create observed and expected item-score by restscore tables
    add low_extreme to cell (item score 0, restscore 0)
    add high_extreme to cell (max item score, max restscore)

    for each score s in the source score window
      observed_count[k] <- count observed_rows with S[p] = s and Y[p,i] = k
      observed_frequency[k] <- observed_count[k] / sum_k observed_count[k]

      for each possible item score k
        rest <- s - k
        observed[k, rest] += observed_frequency[k] * item_margin_n[s]
        expected[k, rest] += Pr_i,s(k) * item_margin_n[s]
      end for
    end for

    do not add the legacy incomplete item-restscore gamma branch here.
    The supplied example compact gamma reports do not print the source
    line announcing that incomplete item-restscore records were included. The
    original incomplete branch uses an inherited `iscore`/`itemscore`
    expected-table indexing convention, so the R port keeps that branch absent
    unless a runtime report explicitly exercises it.

    ObsGamma <- GoodmanKruskalGamma(observed)
    ExpGamma <- GoodmanKruskalGamma(expected)
    sdGamma <- sqrt(source CalculateFittedGAMMA variance)
    p_gamma <- 2 * PNORMAL(abs((ObsGamma - ExpGamma) / sdGamma))
  end for
end algorithm
```

```text
algorithm compact_item_fits_bh(rows)
  infit risk marks <- BenjaminiHochberg1(p_infit for all items)
  outfit risk marks <- BenjaminiHochberg1(p_outfit for all items)
  gamma risk marks <- BenjaminiHochberg1(p_gamma for all items)

  for each diagnostic-specific risk mark:
    display 0 as "", 1 as "*", 2 as "**", and 3 as "***" in a
    fixed three-character suffix field after the corresponding p-value
  end for

  all_p <- concatenate p_infit, p_outfit, p_gamma in source order
  fdr_5_limit <- BenjaminiHochberg(all_p, alpha = 0.05)
  fdr_1_limit <- BenjaminiHochberg(all_p, alpha = 0.01)

  direction is "low" when any item risk mark is present and ObsGamma < ExpGamma
  direction is "high" when any item risk mark is present and ObsGamma > ExpGamma
end algorithm
```

The three marker vectors answer three different questions: which outfits are
selected when outfits are adjusted across items, which infits are selected
when infits are adjusted across items, and which item-restscore gamma tests are
selected when gamma tests are adjusted across items. The combined limits below
the table answer a fourth question involving the pooled collection of all
three kinds of p-value. Consequently, a p-value can have a source marker even
when it exceeds the displayed combined limit, or lack a marker even when it is
below that limit.

The `low` or `high` annotation is item-level rather than diagnostic-specific.
It is shown when any one of the item's three risk marks is nonzero and records
whether the observed item-restscore gamma is below or above its fitted value.
It therefore can appear on a row whose outfit or infit marker is blank; the
gamma test, or the other fit statistic, may be the diagnostic that selected the
item.

The current R package exposes compact item-fit values through `item_fit()`.
The returned test table keeps p-values numeric and places the corresponding
`Outfit FDR`, `Infit FDR`, and `Gamma FDR` marker columns immediately after
them. The source risk grades and the separately pooled BH limits remain
available in the backend attributes; the package does not assemble the full
DIGRAM fixed-width report text.

## Extended ItemFits Report

The extended ItemFits design is a source trace for the report generated by
`skbias15.pas::Calculate_residuals_and_item_fits`, especially the nested
`CalculateInAndOutfits` branch. The example runtime call comes from `DGRirtD.pas`
with `LeastScore=1`, `LargestScore=highest_possible_score-1`,
`printOutput=true`, `printICC=true`, `CalculateIRgamma=true`, and
`ExtendedICCoutput=ScaleDescription`.

Implementation status: the Pascal non-GUI path now reproduces the supplied example
`ItemFits-extended.txt` with focused text parity. The independent R path has the
main score-level value surface and source-shaped anchors, but it must not be
called complete until the full extended value tests pass without calling Pascal
or reading the supplied DIGRAM report in production code. The
`IncompleteItemFits` branch remains outside the validated example slice because the
supplied example extended report does not print the source line saying incomplete
responses were included.

```text
algorithm extended_item_fits_values(project, fit)
  compute Pr_i,s(k), mu_i,s, var_i,s, fourth_i,s as in item_fits_values

  score_weight[s] <- count complete item-response rows with total score s
  observed_rows <- complete item-response rows
  if exogenous variables are selected
    observed_rows <- observed_rows with source-valid exogenous values
  end if

  for each item i
    initialize item-level outfit and infit summary accumulators

    for each score s in LeastScore..LargestScore
      if score_weight[s] = 0
        continue
      end if

      observed_count[k] <- count observed_rows with S[p] = s and Y[p,i] = k
      observed_total <- sum_k observed_count[k]

      if observed_total > 0
        observed_frequency[k] <- observed_count[k] / observed_total
        observed_source_value[k] <- observed_frequency[k] * score_weight[s]
        printed_observed_frequency[k] <- source_print_round(observed_source_value[k])
      else
        observed_frequency[k] <- 0
        observed_source_value[k] <- 0
        printed_observed_frequency[k] <- 0
      end if

      residual[k] <- k - mu_i,s
      squared_residual[k] <- residual[k]^2
      standardized[k] <- residual[k] / sqrt(var_i,s)
      squared_standardized[k] <- standardized[k]^2

      if observed_total > 0 and var_i,s > 0
        outfit_contribution <- sum_k observed_frequency[k] * squared_standardized[k]
        outfit_expected <- sum_k Pr_i,s(k) * squared_standardized[k]
        outfit_variance <- sum_k Pr_i,s(k) * squared_standardized[k]^2
                           - outfit_expected^2
        outfit_standard_error <- sqrt(outfit_variance / score_weight[s])
        outfit_z <- (outfit_contribution - 1) / outfit_standard_error
        outfit_p <- 2 * PNORMAL(abs(outfit_z))
      else
        outfit_contribution <- 0
        outfit_expected <- 0
        outfit_variance <- 0
        outfit_standard_error <- 0
        outfit_z <- missing
        outfit_p <- missing
      end if

      infit_average <- sum_k observed_frequency[k] * squared_residual[k]
      infit_expected <- var_i,s
      infit_ratio <- infit_average / infit_expected
      infit_variance <- fourth_i,s - infit_expected^2

      append one extended score row containing:
        n = score_weight[s]
        observed source values = observed_source_value
        printed observed frequencies = printed_observed_frequency
        probabilities = Pr_i,s(k)
        residuals, squared residuals, standardized residuals
        outfit contribution, standard error, z, p
        infit average, expected, ratio, variance

      update outfit summary accumulators using score_weight[s]
      update infit summary accumulators using score_weight[s] * var_i,s
    end for

    append per-score `Outfit summary` rows for item i
    append one `Infit Summary` row for item i
  end for
end algorithm
```

```text
algorithm extended_item_restscore_tables(project, fit)
  # Source trace: skbias15.pas::Calculate_item_restscore_gamma and
  # Calculate_item_restscore_gamma1, printed with PrintCrosstab.

  build source score/item arrays:
    Observed[score,item,item_score] from complete item rows that pass any
      background/status range check
    ObservedGamma[score,item,item_score] from valid-background effective scores
    Expected[score,item,item_score] from the complete extended score
      distribution times Pr_i,s(item_score)

  for each item i
    rest_max <- max_total_score - max_item_score[i]

    initialize global observed and expected tables over item score x restscore
    set low-score and high-score boundary cells from base Counts.ScoreCounts

    for score s from 1 to max_total_score - 1
      n_base <- Counts.ScoreCounts[s]
      for item score k from 0 to max_item_score[i]
        r <- s - k
        if 0 <= r <= rest_max
          observed[k,r] += ObservedGamma[s,i,k]
          expected[k,r] += n_base * Pr_i,s(k)
        end if
      end for
    end for

    emit global observed and expected item-restscore tables
    compute observed gamma with CalculateGamma
    compute fitted gamma, sd, and p with CalculateFittedGamma
  end for

  for each adjacent item score a in 0..max_item_score-1
    for each item i
      initialize adjacent observed and expected tables
      add low boundary only for a = 0
      add high boundary only for a = max_item_score[i] - 1

      for score s from 1 to max_total_score - 1
        obs_total <- sum_k Observed[s,i,k]
        exp_total <- sum_k Expected[s,i,k]
        for item score k in {a, a + 1}
          r <- s - k
          if 0 <= r <= rest_max
            observed[k,r] += Observed[s,i,k]
            if exp_total > 0
              expected[k,r] += (Expected[s,i,k] / exp_total) * obs_total
            end if
          end if
        end for
      end for

      emit adjacent observed and expected tables and their n totals
      compute observed gamma, fitted gamma, sd, and p
    end for

    apply BenjaminiHochberg1 across item p-values for this adjacent score pair
    emit the local item-restscore gamma summary for the adjacent score pair
  end for
end algorithm
```

Extended item-fit rows are retained as structured detail tables under the
item-fit value object. They are not routed through a package-owned text
renderer.

## Item-Parameters Bundle Construction

The bundle builder converts DIGRAM raw categories into the zero-based item score
surface used by the R estimator. Items are recoded from `1..raw_max` to
`0..raw_max-1`; background variables remain one-based.

```text
algorithm build_item_parameters_bundle(project)
  items <- project.items
  backgrounds <- project.backgrounds
  raw <- project.raw_data

  initialize item_data, background_data, score, status, missing flags
  initialize manifest counters to zero

  for each row r in raw
    row_score <- 0
    complete_items <- true
    complete_backgrounds <- true

    for each item i
      value <- raw[r, items[i].position]

      if value < 1 or value > items[i].raw_max
        item_data[r,i] <- -1
        complete_items <- false
        missing_items[r] <- 1
      else
        item_data[r,i] <- value - 1
        row_score <- row_score + item_data[r,i]
      end if
    end for

    for each background b
      value <- raw[r, backgrounds[b].position]

      if value < 1 or value > backgrounds[b].raw_max
        background_data[r,b] <- -1
        complete_backgrounds <- false
        missing_backgrounds[r] <- 1
      else
        background_data[r,b] <- value
      end if
    end for

    if complete_items
      n_complete_items <- n_complete_items + 1
      complete_item_scores[n_complete_items] <- row_score
    end if

    if complete_backgrounds
      n_complete_backgrounds <- n_complete_backgrounds + 1
    end if

    if complete_items and complete_backgrounds
      score[r] <- row_score
    else
      status[r] <- 0
      score[r] <- -1
    end if

    complete_item_flags[r] <- complete_items
    complete_background_flags[r] <- complete_backgrounds
    row_scores[r] <- row_score
  end for

  least_score <- 1
  largest_score <- maximum complete item score, or 0 when there are no
                   complete item rows

  for each row r
    if score[r] is in least_score..largest_score
      status[r] <- 1
      n_valid <- n_valid + 1
    end if
  end for

  for each row r
    if complete_item_flags[r]
      # Source Count_Margins checks the complete item score window before
      # reading/checking exogeneous values. Boundary complete-item rows do not
      # become Nuseless merely because their backgrounds are missing.
      if row_scores[r] is outside least_score..largest_score
        continue
      end if
      if not complete_background_flags[r]
        n_useless <- n_useless + 1
      end if
    else if complete_background_flags[r]
      n_incomplete <- n_incomplete + 1
    else
      n_useless <- n_useless + 1
    end if
  end for

  return bundle with:
    model.items
    model.backgrounds
    model.max_total_score = sum_i(items[i].raw_max - 1)
    model.least_score = least_score
    model.largest_score = largest_score
    manifest counters
    data = item_data, background_data, score, status, flags
end algorithm
```

## Base Rasch Estimation

### Observed Margins

Source trace: `LoadCounts`.

Only rows with `status == 1` contribute to the estimator.

```text
algorithm rasch_counts(bundle)
  valid_rows <- rows where bundle.data.status == 1

  for each item i
    for x in 0..m_i
      N_ix[i,x] <- count valid rows r where Y[r,i] == x
    end for
  end for

  for s in 0..M
    N_s[s] <- count valid rows r where score[r] == s
  end for

  return { n_valid, N_ix, N_s }
end algorithm
```

### Initial Gamma Parameters

Source trace: `InitializeRaschFit`.

```text
algorithm initial_item_gamma(bundle)
  for each item i
    for x in 0..m_i
      g[i,x] <- 1
    end for
  end for

  for any unused rectangular matrix cells
    g[cell] <- 0
  end for

  return g
end algorithm
```

### Score-Gamma Convolution Excluding One Item

Source trace: `BuildGammaExcludingItem`.

This is the coefficient table of the score-generating polynomial for all items
except one. In polynomial notation it is equivalent to:

```text
G_without_i(t) = product over j != i of (sum_{x=0..m_j} g[j,x] * t^x)
```

The implementation uses dynamic programming rather than constructing a
polynomial object.

```text
algorithm build_gamma_excluding_item(bundle, g, excluded_item)
  G[0] <- 1
  G[s > 0] <- 0
  current_max <- 0

  for item j in 1..I
    if j == excluded_item
      continue
    end if

    Next[s] <- 0 for all s
    next_max <- current_max + m_j

    for s in 0..current_max
      if G[s] != 0
        for x in 0..m_j
          Next[s + x] <- Next[s + x] + G[s] * g[j,x]
        end for
      end if
    end for

    G[0..next_max] <- Next[0..next_max]
    current_max <- next_max
  end for

  return G
end algorithm
```

### Expected Item Margins

Source trace: `CalculateRaschExpectedItems`.

For a fixed total score `s`, the conditional probability that item `i` has score
`x` is:

```text
P(X_i = x | S = s)
  = g[i,x] * G_without_i[s - x]
    / sum_{z=0..m_i, z<=s} g[i,z] * G_without_i[s - z]
```

The expected item margin is the observed total-score count multiplied by this
conditional probability and summed over total scores.

```text
algorithm calculate_rasch_expected_items(bundle, counts, g)
  E_ix[i,x] <- 0 for all i,x

  for item i in 1..I
    G_without_i <- build_gamma_excluding_item(bundle, g, i)

    for total score s in 0..M
      if N_s[s] == 0
        continue
      end if

      denominator <- 0
      for x in 0..m_i
        if s >= x
          denominator <- denominator + g[i,x] * G_without_i[s - x]
        end if
      end for

      if denominator <= 0
        continue
      end if

      for x in 0..m_i
        if s >= x
          numerator <- g[i,x] * G_without_i[s - x]
          E_ix[i,x] <- E_ix[i,x] + N_s[s] * numerator / denominator
        end if
      end for
    end for
  end for

  return E_ix
end algorithm
```

### IPF Update Ratios

Source trace: `CalculateRaschUpdateRatiosFromScore` and
`CalculateRaschUpdateRatios`.

DIGRAM's convergence delta is the maximum absolute count discrepancy, not the
maximum gamma-ratio movement.

```text
algorithm calculate_rasch_update_ratios(counts, E_ix, g, apply_update)
  delta <- 0
  update[i,x] <- 1 for all i,x
  next_g <- g

  for item i in 1..I
    for score x in 0..m_i
      observed <- N_ix[i,x]
      fitted <- E_ix[i,x]

      if observed > 0
        if fitted > 0
          ratio <- observed / fitted
        else
          ratio <- 1
        end if

        delta <- max(delta, abs(fitted - observed))
      else
        ratio <- 0
      end if

      update[i,x] <- ratio

      if apply_update
        next_g[i,x] <- next_g[i,x] * ratio
      end if
    end for
  end for

  return { item_gamma = next_g, update, delta }
end algorithm
```

### Source Gamma Normalization

Source trace: `AdjustItemGammasSourceScale`.

This normalization is not a statistical re-estimation step; it keeps the fitted
parameters on the same source reporting scale as DIGRAM.

```text
algorithm adjust_item_gammas_source_scale(bundle, g)
  last_sgamma <- 1
  s_min <- 0
  s_max <- 0

  for item i in 1..I
    fra <- 0
    til <- m_i
    s_min <- s_min + fra
    s_max <- s_max + til

    alpha <- g[i,fra]
    if alpha > 0
      for x in fra..til
        g[i,x] <- g[i,x] / alpha
      end for
    end if

    if g[i,til] > 0
      last_sgamma <- last_sgamma * g[i,til]
    end if
  end for

  alpha <- 0
  if s_max - s_min > 0
    alpha <- -log(last_sgamma) / (s_max - s_min)
  end if

  for item i in 1..I
    for x in 0..m_i
      g[i,x] <- exp(x * alpha) * g[i,x]
    end for
  end for

  return g
end algorithm
```

### Base Rasch Iteration Loop

Source trace: `EstimateBaseRaschParameters`, base `base_Rasch_only` branch.

The R implementation fits the base Rasch item-parameter path currently needed
by the source item-parameters report.

```text
algorithm fit_rasch_base(bundle, max_step, max_delta)
  counts <- rasch_counts(bundle)
  g <- initial_item_gamma(bundle)
  n_step <- 0
  report_delta <- 0
  converged <- false

  if counts.n_valid == 0 or number of items == 0
    return trivial converged fit
  end if

  while n_step < max_step
    n_step <- n_step + 1

    E_ix <- calculate_rasch_expected_items(bundle, counts, g)

    ratio_state <- calculate_rasch_update_ratios(counts, E_ix, g, false)
    delta <- ratio_state.delta
    update <- ratio_state.update

    ratio_state <- calculate_rasch_update_ratios(counts, E_ix, g, true)
    g <- adjust_item_gammas_source_scale(bundle, ratio_state.item_gamma)

    if delta <= max_delta
      converged <- true
      report_delta <- delta
      break
    end if
  end while

  if not converged or report_delta == 0
    report_delta <- delta
  end if

  E_ix <- calculate_rasch_expected_items(bundle, counts, g)
  ratio_state <- calculate_rasch_update_ratios(counts, E_ix, g, false)
  delta <- ratio_state.delta
  update <- ratio_state.update
  converged <- delta <= max_delta

  return {
    n_step,
    delta,
    report_delta,
    converged,
    item_gamma = g,
    expected_items = E_ix,
    update_items = update,
    counts
  }
end algorithm
```

## Derived Item-Parameter Values

### PCM Thresholds

Source trace: `SourceItemThresholdFromGamma`.

For score `x >= 1`:

```text
threshold[i,x] = log(g[i,x-1] / g[i,x])
```

The source sentinel convention is preserved:

```text
algorithm source_threshold_from_gamma(gamma_values, x)
  epsilon <- 1e-12

  if gamma_values[x] <= epsilon
    return 9999999
  else if gamma_values[x - 1] <= epsilon
    return -9999999
  else
    return log(gamma_values[x - 1] / gamma_values[x])
  end if
end algorithm
```

### PCM Location

Source trace: `SourceItemLocationFromGamma`.

The item location is the mean of finite adjacent-category thresholds.

```text
algorithm source_location_from_gamma(gamma_values, max_score)
  total <- 0
  n_finite <- 0

  for x in 1..max_score
    threshold <- source_threshold_from_gamma(gamma_values, x)

    if abs(threshold) < 9999999
      total <- total + threshold
      n_finite <- n_finite + 1
    end if
  end for

  if n_finite > 0
    return total / n_finite
  else
    return 0
  end if
end algorithm
```

### Expected Score at Theta

Source trace: `CalculateTrueScoreFromGamma`.

`theta` is on the multiplicative scale.

```text
algorithm true_score_from_gamma(theta, max_score, gamma_values)
  if theta <= 0
    return 0
  end if

  numerator <- 0
  denominator <- 0
  power <- 1

  for x in 0..max_score
    if x == 0
      power <- 1
    else
      power <- power * theta
    end if

    numerator <- numerator + x * power * gamma_values[x]
    denominator <- denominator + power * gamma_values[x]
  end for

  if denominator > 0
    return numerator / denominator
  else
    return 0
  end if
end algorithm
```

### Initial Theta for Person-Parameter Search

Source trace: `InitializePersonParameter2`.

```text
algorithm initialize_person_parameter2(target_score, max_score, gamma_values)
  max_diff <- 10 * max_score
  theta <- 1

  for k in 0..100
    trial_theta <- exp(-5 + k / 10)
    fitted_score <- true_score_from_gamma(trial_theta, max_score, gamma_values)
    diff <- abs(fitted_score - target_score)

    if diff < max_diff
      theta <- trial_theta
      max_diff <- diff
    end if
  end for

  return theta
end algorithm
```

### Newton Estimate of Theta

Source trace: `EstimatePersonParameter`.

The R implementation keeps the Pascal update formula. The target equation is:

```text
E(X | theta) = target_score
```

Boundary scores use source sentinel values.

```text
algorithm estimate_person_parameter(target_score, max_score, gamma_values,
                                    max_iterations, max_delta)
  theta <- initialize_person_parameter2(target_score, max_score, gamma_values)

  if target_score == 0
    return 0
  end if

  if target_score == max_score
    return 999999
  end if

  if theta <= 0
    theta <- max_delta
  end if

  for iteration in 1..max_iterations
    H0 <- 0
    H1 <- 0
    H2 <- 0
    factor <- 1 / theta

    for x in 0..max_score
      factor <- factor * theta
      product0 <- gamma_values[x] * factor
      product1 <- product0 * x
      product2 <- product1 * x

      H0 <- H0 + product0
      H1 <- H1 + product1
      H2 <- H2 + product2
    end for

    step <- -theta * (target_score * H0 - H1)

    if H0 > 0
      if H1 == H0
        ratio <- 1
      else
        ratio <- H1 / H0
      end if

      denominator <- ratio * H1 - H2
      if denominator != 0
        step <- step / denominator
      end if

      theta <- theta + step
    end if

    if theta <= 0
      theta <- max_delta
    end if

    if abs(step) < 0.001 or theta <= 0 or theta > 999999
      break
    end if
  end for

  return theta
end algorithm
```

### Midpoint Difficulty

Source trace: `SourceItemDifficultyFromGamma`.

For dichotomous items, DIGRAM uses `-log(gamma[1])` with the source zero
sentinel. For polytomous items, it estimates the multiplicative theta whose
expected item score is half the maximum score, and returns `log(theta)`.

```text
algorithm source_difficulty_from_gamma(gamma_values, max_score)
  if max_score == 1
    return -ln_zero(gamma_values[1])
  end if

  theta <- estimate_person_parameter(
    target_score = max_score / 2,
    max_score = max_score,
    gamma_values = gamma_values,
    max_iterations = 1000
  )

  return ln_zero(theta)
end algorithm
```

### Score Probabilities at Log Theta

Source trace: `ScoreProbabilitiesForLogTheta`.

The unnormalized score mass is:

```text
mass[x] = exp(log_theta * x) * gamma_values[x]
```

```text
algorithm score_probabilities_for_log_theta(log_theta, max_score, gamma_values)
  theta <- exp(log_theta)
  power <- 1
  total <- 0

  for x in 0..max_score
    if x == 0
      power <- 1
    else
      power <- power * theta
    end if

    prob[x] <- power * gamma_values[x]
    total <- total + prob[x]
  end for

  mean <- 0
  sum_sq <- 0

  if total > 0
    for x in 0..max_score
      prob[x] <- prob[x] / total
      mean <- mean + x * prob[x]
      sum_sq <- sum_sq + x * x * prob[x]
    end for
  end if

  variance <- sum_sq - mean * mean
  return { probabilities = prob, mean_score = mean, score_variance = variance }
end algorithm
```

### Item Target and Information at Target

Source trace: `CalculateSourceItemTarget`.

In this implementation, item information at a given `log_theta` is the score
variance returned by `score_probabilities_for_log_theta`. The target is the
`log_theta` where that variance is maximized by the source grid search.

```text
algorithm source_item_target(gamma_values, max_score)
  target <- 0
  max_info <- 0
  best_target <- target
  degenerate <- false
  d <- 0.01

  for ii in -5000..5000
    candidate <- target + ii * d / 10
    props <- score_probabilities_for_log_theta(candidate, max_score, gamma_values)

    if props.score_variance == 0
      best_target <- candidate
      max_info <- 0
      degenerate <- true
      break
    end if

    if props.score_variance > max_info
      max_info <- props.score_variance
      best_target <- candidate
    end if
  end for

  target <- best_target
  d <- d / 10

  if not degenerate
    for ii in -25..25
      candidate <- target + ii * d / 25
      props <- score_probabilities_for_log_theta(candidate, max_score, gamma_values)

      if props.score_variance > max_info
        max_info <- props.score_variance
        best_target <- candidate
      end if
    end for
  end if

  return { target = best_target, info = max_info }
end algorithm
```

### Item Information Function

Source trace: `SourceItemInformation`.

This function is implemented for source parity and future use, although the
current item-statistics target values use the score-variance search above.

```text
algorithm source_item_information_from_gamma(gamma_values, max_score, log_theta)
  G <- 0
  G1 <- 0
  G2 <- 0

  for x in 0..max_score
    FX[x] <- exp(log_theta * x) * gamma_values[x]
    F1X[x] <- x * FX[x]
    F2X[x] <- x * F1X[x]

    G <- G + FX[x]
    G1 <- G1 + F1X[x]
    G2 <- G2 + F2X[x]
  end for

  if G <= 0
    return 0
  end if

  G_square <- G * G
  G_cube <- G_square * G
  information <- 0

  for x in 0..max_score
    PX <- FX[x] / G
    P1X <- F1X[x] / G - FX[x] * G1 / G_square
    P2X <- F2X[x] / G
           - 2 * F1X[x] * G1 / G_square
           - FX[x] * G2 / G_cube
           + 2 * FX[x] * G1 * G1 / G_cube

    if PX > 0
      information <- information + P1X * P1X / PX - P2X
    end if
  end for

  return information
end algorithm
```

### Observed Score Range

Source trace: `example_ITEM_PARAMETERS_REPORT.pas::CalculateObservedScoreRange`.

This is not the empirical range of total scores. It sums the lowest and highest
observed supported category for each item.

```text
algorithm calculate_observed_score_range(N_ix)
  score_min <- 0
  score_max <- 0

  for item i in 1..I
    supported_scores <- { x : N_ix[i,x] > 0 }

    if supported_scores is not empty
      score_min <- score_min + min(supported_scores)
      score_max <- score_max + max(supported_scores)
    end if
  end for

  return (score_min, score_max)
end algorithm
```

### Source Parameter Count

Source trace: `example_ITEM_PARAMETERS_REPORT.pas::CalculateSourceNParameters`.

```text
algorithm calculate_source_n_parameters(N_ix)
  result <- -1

  for item i in 1..I
    result <- result + count{x : N_ix[i,x] > 0} - 1
  end for

  return result
end algorithm
```

### Report-Ready Item Parameter Values

Source trace: `EmitGLLRMOutputRows` output(4) plus the example report footer
helpers.

The ICE/MICE rows in this report use the output(4) emission formula:

```text
z_i = log(g[i,m_i]) / m_i
ICE[i,x] = log(g[i,x]) - x * z_i
MICE[i,x] = exp(ICE[i,x])
ICE_item_effect[i] = z_i
MICE_item_effect[i] = exp(z_i)
```

This is distinct from the lower-level
`CalculateSourceICEAndMICEFromGamma` helper, which centers by PCM location and
is used by other source paths. The package stores numeric ICE/MICE values from
the output(4) formula and separately stores fixed-width source ICE fields from
the lower-level source path for report parity checks.

```text
algorithm item_parameters_input_stats(bundle)
  valid_scores <- empty list
  n_incomplete <- 0
  n_responses <- 0
  n_useful <- 0
  n_useless <- 0

  for each record r
    complete_items <- all item fields are valid
    complete_backgrounds <- all background fields are valid
    score <- sum valid item scores in r
    responses <- count valid item responses in r

    if complete_items
      if score is outside the estimated score range
        skip record before checking background validity
      end if

      if not complete_backgrounds
        n_useless <- n_useless + 1
      end if
    else if complete_backgrounds
      n_incomplete <- n_incomplete + 1

      if responses > 1 and score > 0 and score < maximum possible score
        n_useful <- n_useful + 1
        n_responses <- responses
      end if
    else
      n_useless <- n_useless + 1
    end if
  end for

  return valid score min/max, nvalid, nincomplete, nresponses, nuseful, nuseless
end algorithm
```

Note that `Nresponses` follows `example_ITEM_PARAMETERS_REPORT.pas::ReadDataStats`:
it is overwritten by each useful incomplete row and therefore records the
latest useful-row response count, not a sum.

```text
algorithm item_parameters_values(fit, bundle)
  for each item i
    gamma_values <- fit.item_gamma[i, 0..m_i]

    for x in 1..m_i
      thresholds[i,x] <- source_threshold_from_gamma(gamma_values, x)
    end for

    locations[i] <- source_location_from_gamma(gamma_values, m_i)
    midpoints[i] <- source_difficulty_from_gamma(gamma_values, m_i)

    target_result <- source_item_target(gamma_values, m_i)
    targets[i] <- target_result.target
    info_at_target[i] <- target_result.info

    if m_i >= 2
      info_per_step[i] <- info_at_target[i] / m_i
    else
      info_per_step[i] <- NA
    end if

    if m_i > 0 and gamma_values[m_i] > 0
      z <- log(gamma_values[m_i]) / m_i

      for x in 0..m_i
        if gamma_values[x] > 0
          ICE[i,x] <- log(gamma_values[x]) - x * z
          if x == m_i
            ICE[i,x] <- source signed-zero cancellation of
                        log(gamma_values[m_i]) - m_i * z
          end if
          MICE[i,x] <- exp(ICE[i,x])
        else
          ICE[i,x] <- 0
          MICE[i,x] <- 0
        end if
      end for

      ICE_item_effect[i] <- z
      MICE_item_effect[i] <- exp(z)
    end if
  end for

  observed_score_range <- calculate_observed_score_range(fit.counts.N_ix)
  n_parameters <- calculate_source_n_parameters(fit.counts.N_ix)
  log_likelihood <- base_rasch_loglike(bundle, fit.item_gamma)
  likelihood_n <- fit.counts.n_valid
  input_stats <- item_parameters_input_stats(bundle)
  ice_fields <- source_item_parameters_extended_ice_fields(fit, bundle)

  return report-ready values object
end algorithm
```

Implementation note: fixed-width `ice_fields` are derived from the already
fitted gamma matrix through the native `gRm_item_parameters_ice_fields_from_gamma`
boundary. The package no longer has a separate native entry point that refits
counts just to emit ICE fields. The still-registered extended-gamma entry point
is a diagnostic/parity surface for the source-shaped extended arithmetic, not
the production source of item-parameter ICE fields.

## Item Parameter Tables

Implemented in:

- `gRm/R/item_parameters_values.R`
- `gRm/R/api-results.R`
- `gRm/R/api-summary.R`
- `gRm/R/api-summary-tables.R`
- `gRm/R/api-table-helpers.R`

The current package keeps item-parameter output as structured value objects,
summary tables, and package-owned helper tables. It does not expose a DIGRAM
report renderer or an internal output accessor API.

```text
algorithm item_parameter_result_tables(values)
  build item_gamma table with one row for each item and score category
  build threshold table with finite source PCM thresholds
  build ICE and MICE tables from values.ice, values.mice, and item effects
  build item_statistics from location, midpoint, target, information, and
    per-step information
  build fit_summary with n_step, delta, observed score range, n_parameters,
    log_likelihood, AIC, BIC, and likelihood_n
  return named structured tables for summary() and package tests
end algorithm
```

## GLLRM Estimation

Implemented in:

- `gRm/R/gllrm_context.R`
- `gRm/R/gllrm_components.R`
- `gRm/R/gllrm_candidate_fit.R`
- `gRm/R/gllrm_fit.R`
- `gRm/R/gllrm_probability_cache.R`
- `gRm/R/gllrm_values.R`

The GLLRM path is used when a model specification contains local
dependence or DIF terms. It extends the source-faithful Rasch IPF machinery with
enumerated local-dependence components and background-specific DIF factors.

```text
algorithm fit_gllrm(spec, max_step, max_delta, max_joint_configs, bundle)
  bundle <- supplied bundle or build_item_parameters_bundle(spec.project)
  context <- build_gllrm_context(spec, bundle, max_joint_configs)
  state <- initialize_gllrm_state(context)

  if no valid rows or no items:
    mark converged and return context, state, and bundle

  repeat until max_step or a source-shaped stop condition:
    increment n_step
    calculate_gllrm_joint_expected_margins(context, state)
    update item, LD, and DIF parameters from observed / expected margins
    adjust LD/DIF reference cells and source-scale item gammas
    record delta and report_delta
    if delta < max_delta:
      mark converged and stop with stop_reason = "delta_below_tolerance"
    else if max_step is reached:
      stop with stop_reason = "max_step"
    else if source periodic checkpoint condition fires:
      stop with stop_reason = "source_periodic_checkpoint"
    else if repeated deltas are detected:
      stop with stop_reason = "repeated_delta" or "two_back_repeated_delta"
    else if non-improvement persists:
      stop with stop_reason = "finish_count_plateau"
    end if
  end repeat

  if the fit did not converge, or no report delta has been recorded:
    report_delta <- final delta
  end if
  recalculate expected margins without applying another update
  recompute final observed/expected deltas for detail tables
  state.log_likelihood <- gllrm_loglike(context, state)
  return context, state, bundle, converged, report_delta, and stop_reason
end algorithm
```

```text
algorithm build_gllrm_context(spec, bundle, max_joint_configs)
  extract item and background metadata plus zero-based item/background matrices
  valid_rows <- rows with source status 1
  translate LD and DIF formula terms to item/background indices
  build connected item components induced by LD terms
  precompute component configurations, scores, lookup maps, and local LD
    metadata for every component
  compute observed item counts from rasch_counts(bundle)
  compute observed LD margins over valid rows for each GLLRM item pair
  compute observed DIF margins over valid rows for each item/background term
  aggregate valid rows by total score and included DIF background values
  stop if any LD component has more than max_joint_configs score
    configurations
  return context
end algorithm
```

```text
algorithm gllrm_component_gamma(context, state, component_items, background_values)
  enumerate the component's item-score configurations
  for each configuration:
    weight <- product of item gammas for the configured scores
    multiply any included DIF gamma for each item and background value
    multiply any included LD gamma whose two items are inside the component
    add weight to gamma[sum(configuration scores)]
  normalize gamma by its maximum positive entry and keep the scale separately
  return normalized component gamma, scale, and configuration weights
end algorithm
```

```text
algorithm calculate_gllrm_joint_expected_margins(context, state)
  zero expected item, LD, and DIF margins
  cache component gammas and prefix/suffix convolutions by DIF background
    value pattern

  for each score/background group:
    get cached component gammas and leave-one-component convolutions
    denominator <- full score gamma at the group's total score
    for each component:
      for each component configuration compatible with the total score:
        expected mass <- group count * component weight * rest gamma /
                         denominator, accounting for component scale
        add expected mass to the corresponding item-score margins
        add expected mass to matching included DIF margins
        add expected mass to matching included LD pair margins
      end for
    end for
  end for

  return state with expected margins
end algorithm
```

Implementation note: `calculate_gllrm_joint_expected_margins()` may dispatch to
the native `gRm_gllrm_expected_margins` backend for speed. The native payload is
only a compact representation of the R context/state boundary: item gamma must
have one row per item, LD/DIF parameter arrays come from the GLLRM state, and
component identity is reconstructed from the structured component metadata. Dead
debug fields such as `component_keys` are not part of the current boundary.

```text
algorithm update_gllrm_parameters(context, state, apply_update, track_delta)
  update item gammas through the same Rasch observed/expected ratio helper used
    by the base model
  for each included LD table:
    update positive observed/expected cells by observed / expected
    set structural-zero observed cells to zero
    optionally track maximum absolute margin discrepancy
  for each included DIF table:
    update positive expected cells by observed / expected
    optionally track maximum absolute margin discrepancy
  optionally multiply parameters by the update ratios
  return state with update tables and delta
end algorithm
```

```text
algorithm adjust_gllrm_dependency_parameters(context, state)
  for each LD parameter table:
    apply source reference-cell normalization for observed item-pair cells
    optionally absorb row/column factors into item gammas
  for each DIF parameter table:
    apply source reference-cell normalization across item scores and
      background categories
  rescale item gammas on the source Rasch item scale using log-space
    normalization
  return adjusted state
end algorithm
```

```text
algorithm gllrm_loglike(context, state)
  cache conditional score gamma vectors by DIF background value pattern
  for each valid row:
    numerator <- product of fitted item gammas for the observed item scores
    multiply included DIF gammas matching observed item scores/backgrounds
    multiply included LD gammas matching observed item-pair scores
    probability <- numerator / score_gamma[observed_total_score]
    add -log(probability) when probability is positive
  return conditional negative log likelihood
end algorithm
```

```text
algorithm gllrm_values(gllrm_fit, spec)
  derive item parameter values from GLLRM item gammas using the same threshold,
    location, ICE, MICE, and item-statistic formulas as the base Rasch fit
  attach GLLRM context, item/LD/DIF expected margins, update tables, and raw
    LD/DIF parameter matrices
  build standardized LD and DIF parameter tables by iterated marginal scaling
    and convert source gamma to odds ratios
  expose GLLRM detail tables for item parameters, LD/DIF parameters, expected
    margins, and update ratios
end algorithm
```

## DIF Tests Report

Source trace:

- `source/PAS_scd/DGRirtD.pas`: CHECK D branch headed `Check assumptions of no DIF`
- `source/PAS_scd/DGRirtD.pas`: `lr:=2*abs(Raschloglike-Raschloglike1)`
- `source/PAS_scd/DGRirtD.pas`: `p := pfchi(df,lr)`
- `source/PAS_scd/DGRirtD.pas`: significant MissingDIF rows append
  `IJXgamma(.i,nitems+j.)`
- `source/PAS_skunits/SKbias7.pas`: `inexpensive_itembias1` and
  `Item_Screening` populate `IJXgamma`/`IJXgamma_pvalues`
- `source/PAS_scd/SkStat.pas`: `PFCHI` and `BenjaminiHochberg`
- `pascal_harness/SourceRaschCore.pas`: source-faithful estimator subset and
  included DIF likelihood rows

The current R slice implements the included-DIF IPF estimator directly in R.
Pascal commands remain test/reference oracles only; `dif_tests_values()` does
not call the Pascal harness for production computation.

Inputs:

- Parsed DIGRAM project `P`
- Item responses `Y[p,i]`
- Exogenous variables `X[p,e]`
- Base Rasch/GLLRM model without included DIF
- Candidate model for each item-by-exogenous pair `(i,e)`

```text
algorithm source_dif_tests_values(project)
  fit base model without candidate DIF
  L0 <- source negative log likelihood of the base model

  for each exogenous variable e in source order
    for each item i in source order
      create candidate model M_ie with item_bias[i,e] included
      fit M_ie with the same max_step and max_delta controls
      L1 <- source negative log likelihood of M_ie

      CLR[i,e] <- 2 * abs(L0 - L1)
      df[i,e] <- (number of item non-reference scores)
                 * (number of background non-reference values)
      p[i,e] <- PFCHI(df[i,e], CLR[i,e])
      gamma[i,e] <- source item-screening partial gamma IJXgamma[i,nitems+e]
      p_gamma[i,e] <- source item-screening IJXgamma_pvalues[i,nitems+e]

      store item label, exogenous label, item name, exogenous name,
            CLR[i,e], df[i,e], p[i,e], gamma[i,e], p_gamma[i,e],
            and gamma_source = "item_screening"
    end for
  end for

  BHcrit <- BenjaminiHochberg(p[all item/background tests], alpha = 0.05)
  return tests and BHcrit
end algorithm
```

For an already fitted GLLRM, `gllrm_dif_tests_values()` uses the
same source-facing candidate-refit idea for no-DIF candidates: start from the
current GLLRM, add one missing item-by-exogeneous DIF term, refit with the
requested controls, and compare the GLLRM log likelihoods. The gamma and
p-gamma columns for those no-DIF candidate rows still come from the
item-screening `IJXgamma`/`IJXgamma_pvalues` matrix and are marked
`gamma_source = "item_screening"`; they are not recomputed from the candidate
refit. Already-included DIF terms are reported separately by the conditional
score-stratified included-term test. The removed orphan conditional no-DIF helper
is not a production path.

```text
algorithm fit_dif_candidate(bundle, target item i, exogenous variable e)
  initialize item gammas from the already-fitted base Rasch model
  initialize included DIF gammas D[k,x] <- 1 for each score k and exogenous value x

  count observed score-by-exogenous margins:
    N[s,x] <- number of persons with total score s and X_e = x

  count observed included DIF margins:
    O[k,x] <- number of persons with item score Y_i = k and X_e = x

  repeat until max_delta or max_step:
    for each exogenous value x:
      build item score polynomials W_j,z:
        W_j,z[k] <- item_gamma[j,k]
        if j = i then W_j,z[k] <- W_j,z[k] * D[k,x]

      compute full conditional score polynomial:
        G_x <- convolution over all item polynomials W_j,z

      compute leave-one-item polynomials:
        G_x_without_j <- convolution over all item polynomials except item j

      for each total score s with N[s,x] > 0:
        for each item j and item score k <= s:
          expected contribution <- N[s,x] * W_j,z[k]
                                   * G_x_without_j[s-k] / G_x[s]
          add contribution to expected item margin E_item[j,k]
          if j = i:
            add contribution to expected DIF margin E_dif[k,x]

    for each item margin cell (j,k):
      if observed item count C[j,k] > 0 and E_item[j,k] > 0:
        item_gamma[j,k] <- item_gamma[j,k] * C[j,k] / E_item[j,k]
      else if C[j,k] = 0:
        item_gamma[j,k] <- 0

    for each included DIF margin cell (k,x):
      if O[k,x] > 0 and E_dif[k,x] > 0:
        D[k,x] <- D[k,x] * O[k,x] / E_dif[k,x]
      else if O[k,x] = 0:
        D[k,x] <- 0

    rescale item gammas with the same reference adjustment as the base Rasch fit
    delta <- maximum absolute fitted-observed margin difference in this step
  end repeat

  return item gammas, included DIF gammas, convergence flag, and final delta
end algorithm
```

```text
algorithm candidate_dif_loglike(bundle, fit, target item i, exogenous variable e)
  for each exogenous value x:
    compute G_x, the included-DIF conditional score polynomial

  L <- 0
  for each valid person p:
    s <- total score of p
    x <- exogenous value X_p,e
    numerator <- product over items j of item_gamma[j, Y_p,j]
    numerator <- numerator * D[Y_p,i, x]
    probability <- numerator / G_x[s]
    L <- L - log(probability)
  end for

  return L
end algorithm
```

The current R package exposes DIF checks through `dif()` and
`summary.gRm_dif`, with selected rows determined from the structured
Benjamini-Hochberg threshold. It no longer assembles the historical DIGRAM text
report.

## Local Independence Candidate Report

The local-independence report implementation is a source-shaped R port of the
candidate-addition likelihood-ratio path used by DIGRAM for
`check-local-independence.txt` and
`check-local-independence-extended.txt`. This is the `DGRirtD.pas` `MissingLD`
report path, not the related `SKbias13.PrintLDandDIFmatrix` gamma-matrix
diagnostic output.

Source trace:

- `source/PAS_scd/DGRirtD.pas`: `MissingLD`
- `source/PAS_scd/SkStat.pas`: `PFCHI`
- `source/PAS_skunits/SKmca.pas`: `BenjaminiHochberg`
- `pascal_harness/local_independence_report/example_LOCAL_INDEPENDENCE_REPORT.pas`
- `docs/example_LOCAL_INDEPENDENCE_SOURCE_TRACE.md`

```text
algorithm local_independence_values(project, max_step, max_delta)
  bundle <- build source-shaped item/background bundle from project
  fit base Rasch model M0 by the source IPF loop
  LL0 <- negative conditional log likelihood of M0

  tests <- empty table
  for item i from 1 to nitems - 1:
    for item j from i + 1 to nitems:
      fit candidate model Mij with exactly LD(i,j) included
      LLij <- included-LD negative conditional log likelihood of Mij

      CLRij <- 2 * abs(LL0 - LLij)
      dfij <- candidate LD parameter increment
              = (max_score_i) * (max_score_j) for complete example items
      pij <- PFCHI(dfij, CLRij)
      WPGgammaij <- item-screening weighted partial gamma for item pair (i,j);
        this is a MissingLD report diagnostic, not a candidate-refit statistic

      append item labels, item names, CLRij, dfij, pij, WPGgammaij,
        convergence, delta, stop_reason, attempted-fit metadata, and
        report_value_source
    end for
  end for

  BHcrit <- BenjaminiHochberg(p over all candidate item pairs, alpha = 0.05)
  suggested LD <- pair labels where p <= BHcrit, in tested order
  return all test rows, BHcrit, and suggested LD list
end algorithm
```

For an already fitted GLLRM, `gllrm_local_independence_values()` tests
missing LD candidates by adding one LD term to the current GLLRM and
refitting. If that attempted candidate fit does not converge and `max_step > 51`,
the package follows the CHECK LID report convention by reporting the first
post-50 checkpoint fit (`reported_checkpoint_step = 51`) while preserving
`attempted_*` metadata for the full attempted fit. `report_value_source`
distinguishes checkpoint-reported rows from ordinary attempted-fit rows. The
WPG gamma column remains the item-screening WPG gamma diagnostic, not a
candidate-refit statistic.

```text
algorithm fit_ld_candidate(bundle, item pair i,j)
  initialize item gammas from the base fit
  initialize LD gamma D[x,y] <- 1 for all score cells of items i and j

  count observed included LD margin:
    O[x,y] <- number of valid persons with item i score x and item j score y

  repeat until max_delta or max_step:
    build item components:
      all items except i and j remain single-item score polynomials
      items i and j form one paired component with weights
        W_pair[t] <- sum over x+y=t of gamma_i[x] * gamma_j[y] * D[x,y]

    compute full conditional score polynomial:
      G <- convolution over all single-item components and W_pair

    compute leave-one-component polynomials for every component

    for each total score s with observed score count N[s] > 0:
      for each non-pair item k and item score x <= s:
        expected contribution <- N[s] * gamma_k[x]
                                 * G_without_k[s-x] / G[s]
        add contribution to expected item margin E_item[k,x]

      for each pair score cell (x,y) with x+y <= s:
        expected contribution <- N[s] * gamma_i[x] * gamma_j[y] * D[x,y]
                                 * G_without_pair[s-x-y] / G[s]
        add contribution to E_item[i,x]
        add contribution to E_item[j,y]
        add contribution to E_ld[x,y]

    for each item margin cell (k,x):
      if observed item count C[k,x] > 0 and E_item[k,x] > 0:
        gamma_k[x] <- gamma_k[x] * C[k,x] / E_item[k,x]
      else if C[k,x] = 0:
        gamma_k[x] <- 0

    for each included LD margin cell (x,y):
      if O[x,y] > 0 and E_ld[x,y] > 0:
        D[x,y] <- D[x,y] * O[x,y] / E_ld[x,y]
      else if O[x,y] = 0:
        D[x,y] <- 0

    rescale item gammas with the source reference adjustment
    delta <- maximum absolute fitted-observed margin difference in this step
  end repeat

  return item gammas, LD gammas, convergence flag, and final delta
end algorithm
```

```text
algorithm candidate_ld_loglike(bundle, fit, item pair i,j)
  compute G, the included-LD conditional score polynomial

  L <- 0
  for each valid person p:
    s <- total score of p
    numerator <- product over items k of gamma_k[Y_p,k]
    numerator <- numerator * D[Y_p,i, Y_p,j]
    probability <- numerator / G[s]
    L <- L - log(probability)
  end for

  return L
end algorithm
```

The current R package exposes local-dependence candidates through
`local_dependence()` and `summary.gRm_local_dependence`. Compact versus
extended report text is no longer a package output surface; the full candidate
rows remain available as structured values.

## Retired Global Invariance Scope

Global invariance/global DIF is outside the current package scope. The
implemented public diagnostics are `global_homogeneity()`, `local_dependence()`,
`dif()`, `score_effects()`, `item_fit()`, and `item_parameters()`. Historical
source notes about global invariance remain repository history, but no current
package algorithm or installed API is defined for that report family.

## Items-Select Report

Implemented in:

- `gRm/R/source_score_groups.R`

Source traces:

- `source/PAS_skunits/DGRexe.pas`: `Execute_select_items`
- `source/PAS_skunits/SKbias2.pas`: `SHOW_ITEMS`, `Sort_items`,
  `Calculate_scores`, `Cut_scores`
- `source/PAS_skunits/SKbias7.pas`: `Calculate_ChronbachsAlpha`

```text
algorithm items_select_values(project)
  read selected item metadata from project.items in source order
  read source-coded records from project.raw_data

  for each item i:
    valid_i <- records with raw item response in 1..raw_max_i
    score_i <- raw response - 1
    n_i <- count(valid_i)
    mean_i <- mean(score_i over valid_i)
  end for

  complete <- records with valid responses for all selected items
  for each complete record r:
    item scores <- raw responses - 1
    total_score_r <- sum item scores
  end for

  for each item i:
    complete_mean_i <- mean(score_i over complete records)
    item observed range_i <- min/max score_i over complete records
  end for

  score_distribution[s] <- count complete records with total_score = s
  mean_score <- sum_s s * count_s / n_complete
  score_variance <- (sum_s s^2 * count_s / n_complete - mean_score^2)
                    * n_complete / (n_complete - 1)
  score_sd <- sqrt(score_variance)
  skewness <- source sample-skewness correction from SKbias2.pas

  for each item i:
    item_var_i <- E(score_i^2 over complete records) - E(score_i)^2
  end for
  score_var_population <- E(total_score^2) - E(total_score)^2
  alpha <- k / (k - 1) * (1 - sum_i item_var_i / score_var_population)

  cut <- source median-like score cut:
    accumulate non-extreme score counts until cumulative >= n_complete / 2
    if previous cumulative is closer to half, move cut down by one score

  score groups <- 0..cut and cut+1..observed_max
  return item rows, score distribution, alpha, and score groups
end algorithm
```

## Exo-Select Report

Source trace:

- `source/PAS_scd/SKbias2.pas`: exogenous overview, complete-case counts,
  missing-exogenous diagnostics, and score distribution.
- `source/PAS_skunits/SKbias13.pas`: `StepwiseScoreScreening`.
- `source/PAS_scd/DIGRAM1f.pas`: final selected exogenous-variable list.
- `source/PAS_scd/SKexa1.pas`: Benjamini-Hochberg footer and markers.
- `source/PAS_scd/SKexa2.pas`: exact-test simulation counters in
  `Inexpensive_bt_tests`.
- `source/PAS_scd/SKrandom.pas::GENTAB1`: conditional random table generation
  with fixed margins.
- `source/PAS_scd/SkStat.pas`: normal and chi-square tail probabilities.

```text
algorithm exo_select_values(project, score_cap = 56)
  read item and exogenous metadata from the project tables
  read source-coded records from project.raw_data

  for each record r:
    item_complete_r <- all item responses are in 1..item_raw_max
    if item_complete_r:
      score_r <- sum_i(raw_item_ri - 1)
    end if
    exo_complete_r <- all exogenous responses are in 1..exo_raw_max
  end for

  complete_item_cases <- count(item_complete)
  complete_item_exo_cases <- count(item_complete and exo_complete)

  for each exogenous variable e:
    missing_e <- item_complete records where e is invalid/missing
    known_e <- item_complete records where e is valid
    missing_mean_e <- mean(score over missing_e)
    known_mean_e <- mean(score over known_e)
    sd_missing_e <- sqrt((mean(score^2 over missing_e) - missing_mean_e^2)
                         / count(missing_e))
    sd_known_e <- sqrt((mean(score^2 over known_e) - known_mean_e^2)
                       / count(known_e))
    t_e <- (-known_mean_e + missing_mean_e) / (sd_missing_e + sd_known_e)
    p_e <- 2 * source_normal_upper_tail(abs(t_e))
  end for

  distribution_scores <- scores for records with item_complete and exo_complete
  score_count_s <- count(distribution_scores = s)
  mean <- sum_s s * score_count_s / n
  variance <- (sum_s s^2 * score_count_s / n - mean^2) * n / (n - 1)
  sd <- sqrt(variance)
  skewness <- source sample-skewness correction from SKbias2.pas

  for each exogenous variable e:
    valid_e <- item_complete records where e is valid
    collapsed_score <- min(score, score_cap)
    table <- counts of collapsed_score by exogenous category e
    remove empty score rows

    chi_square <- sum_cells (observed - expected)^2 / expected
    df <- (number_of_nonempty_score_rows - 1) * (exo_categories_e - 1)
    chi_p <- source_pfchi(df, chi_square)

    gamma <- (concordant - discordant) / (concordant + discordant)
    s <- source RC gamma variance numerator from SKbias13/Inexpensive_bt_tests
    gamma_p_one_sided <- source_normal_upper_tail(abs(gamma / (sqrt(s) / ppq)))
    gamma_p_two_sided <- 2 * gamma_p_one_sided

    if exact:
      set the R implementation seed requested by the validation run
      n_chi <- 0
      n_gamma_directional <- 0
      n_gamma_absolute <- 0
      for sim in 1..nsim:
        generated <- GENTAB1_source_shape(table)
        generated_chi <- chi-square statistic for generated with fixed margins
        generated_gamma <- RC gamma for generated
        if generated_chi >= chi_square:
          n_chi <- n_chi + 1
        if gamma > 0 and generated_gamma >= gamma:
          n_gamma_directional <- n_gamma_directional + 1
        else if gamma < 0 and generated_gamma <= gamma:
          n_gamma_directional <- n_gamma_directional + 1
        else if gamma == 0:
          n_gamma_directional <- n_gamma_directional + 1
        if abs(generated_gamma) >= abs(gamma):
          n_gamma_absolute <- n_gamma_absolute + 1
      exact_chi_p <- n_chi / nsim
      exact_gamma_p_one_sided <- n_gamma_directional / nsim
      exact_gamma_p_two_sided <- n_gamma_absolute / nsim
    end if

    selected_e <- min(chi_p, gamma_p_two_sided) <= 0.05
  end for

  if exact:
    bh_p_values <- concatenate all exact_chi_p and exact_gamma_p_two_sided values
  else:
    bh_p_values <- concatenate all chi_p and gamma_p_two_sided values
	  bh_05 <- source BenjaminiHochberg(bh_p_values, alpha = 0.05)
	  bh_01 <- source BenjaminiHochberg(bh_p_values, alpha = 0.01)
	  chi marker_e <- compare asymptotic or exact chi p-value to bh_01 and bh_05
	  gamma marker_e <- compare the source marker p-value to bh_01 and bh_05
	                    and use the sign of gamma. In the asymptotic example
	                    runtime this is the undoubled one-sided gamma p-value;
	                    in the exact EXA_SUMMARY1_2 branch it is the absolute
	                    two-sided simulated gamma p-value from Results[*,9].
	  selected_labels <- labels where selected_e

  return overview, missing diagnostics, score distribution,
         screening rows, BH thresholds, recursive line, selected labels
end algorithm
```

```text
algorithm GENTAB1_source_shape(observed_table)
  row_total <- observed row totals
  column_total <- observed column totals
  new_table <- all zeroes
  generated_column_totals <- zero for each column
  generated_total <- 0

  for each row i except the last:
    generated_row_total <- 0
    free_n <- grand_total - generated_total
    for each column j except the last:
      if remaining column j total is zero or free_n is zero:
        t11 <- 0
      else:
        column1 <- row_total[i] - generated_row_total
        row1 <- column_total[j] - generated_column_totals[j]
        tmin <- max(column1 + row1 - free_n, 0)
        tmax <- min(column1, row1)
        expected <- round(column1 * row1 / free_n)
        draw <- uniform(0, 1)

        visit possible t11 values in the source order:
          expected,
          expected + 1, expected - 1,
          expected + 2, expected - 2, ...,
          then whichever tail remains
        accumulate hypergeometric probabilities until cumulative >= draw
      end if

      new_table[i,j] <- t11
      update generated row, column, and total counters
    end for
    set the last column of row i to the remaining row total
  end for

  set the final row to the remaining column totals
  return new_table
end algorithm
```

## Marginal Gamma Report

Source trace:

- `source/PAS_scd/DGRexe.pas::execute_marginal_gamma`
- `source/PAS_scd/skbig3.pas::quick_tests`
- `source/PAS_scd/SKxyz1.PAS::Transfer_results`
- `source/PAS_scd/SkStat.pas::PREPARE_GAMMA_STATISTICS`
- `source/PAS_scd/SkStat.pas::RCGAMMA`
- `source/PAS_scd/DGRVARS.PAS::Show_gammavalues`

The `gamma.txt` report is the DIGRAM `GAMMA` command's pairwise marginal
Goodman-Kruskal association-gamma matrix. It is not the Rasch score-gamma or
item-gamma array used by item parameter estimation.

```text
algorithm gamma_values(project)
  variables <- project variables in source order
  raw <- project.raw_data source-coded integer matrix
  ordinal_v <- variables[v].raw_max > 1 and variables[v].vtype = 3

  initialize gamma[v,w] <- 999 for all v,w
  initialize ppq[v,w], pmq[v,w], n[v,w] <- 0
  set gamma[v,v] <- 0 for all diagonal cells

  for row variable v from 1 to V - 1:
    for column variable w from v + 1 to V:
      if not ordinal_v or not ordinal_w:
        continue with sentinel gamma = 999
      end if

      table <- zero matrix variables[v].raw_max by variables[w].raw_max
      for each record r:
        x <- raw[r, variables[v].position]
        y <- raw[r, variables[w].position]
        if 1 <= x <= variables[v].raw_max and
           1 <= y <= variables[w].raw_max:
          table[x,y] <- table[x,y] + 1
        end if
      end for

      (gamma_vw, ppq_vw, pmq_vw) <- gRm_goodman_kruskal_gamma(table)
      gamma[v,w] <- gamma_vw
      gamma[w,v] <- gamma_vw
      ppq[v,w] <- ppq[w,v] <- ppq_vw
      pmq[v,w] <- pmq[w,v] <- pmq_vw
      n[v,w] <- n[w,v] <- sum(table)
    end for
  end for

  return variables, gamma, ppq, pmq, n
end algorithm

algorithm gRm_goodman_kruskal_gamma(table)
  P <- 0
  Q <- 0
  for each cell (i,j):
    AIJ[i,j] <- sum of table[k,l] where
                (i > k and j > l) or (i < k and j < l)
    DIJ[i,j] <- sum of table[k,l] where
                (i < k and j > l) or (i > k and j < l)
    P <- P + table[i,j] * AIJ[i,j]
    Q <- Q + table[i,j] * DIJ[i,j]
  end for

  PPQ <- P + Q
  PMQ <- P - Q
  if PPQ > 0:
    Gamma <- PMQ / PPQ
  else:
    Gamma <- 0
  end if

  return Gamma, PPQ, PMQ
end algorithm
```

The current R package retains marginal gamma calculations as structured values
and does not emit the historical blocked text matrix.

## Extended ItemFits item-restscore tables

```text
algorithm extended_item_restscore_source_tables(project, fit)
  complete_score_counts <- score distribution among complete item-response rows
  valid_score_counts <- score distribution after DIGRAM background/exogenous
                        validity filters

  # Source trace:
  # skbias15.pas::Count_Observed increments Nlowscore/Nhighscore before
  # GET_EXOGENE can reject a row, but Add_count_to_tables fills interior
  # ItemMargTables only after exogenous validation succeeds.
  source_score_counts <- valid_score_counts
  source_score_counts[lowest score] <- complete_score_counts[lowest score]
  source_score_counts[highest score] <- complete_score_counts[highest score]

  for each item:
    seed observed and expected global item-restscore tables with the low/high
    endpoint counts from source_score_counts

    for each interior score:
      n_score <- source_score_counts[score]
      for each item score compatible with the total score:
        rest_score <- score - item_score
        observed[item_score, rest_score] +=
          observed item-score relative frequency from ItemMargTables * n_score
        expected[item_score, rest_score] +=
          conditional item-score probability from ExpectedItemMargTables * n_score
      end for
    end for
  end for
end algorithm

algorithm extended_local_item_restscore_tables(project, fit, adjacent_score)
  for each item:
    if adjacent_score = 0:
      seed visible low endpoint cell with Nlowscore
    else if adjacent_score is the top adjacent pair for this item:
      seed visible high endpoint cell with Nhighscore
    end if

    Isum <- 0
    Esum <- 0
    for each interior score:
      for item_score in adjacent_score, adjacent_score + 1:
        rest_score <- score - item_score
        observed_cell <- ItemMargTables relative frequency * interior n_score
        expected_cell <- ExpectedItemMargTables probability * interior n_score
        table observed[item_score, rest_score] += observed_cell
        table expected[item_score, rest_score] += expected_cell
        Isum <- Isum + observed[item_score, rest_score]
        Esum <- Esum + expected[item_score, rest_score]
      end for
    end for

    print the table body including endpoint seed cells
    print observed "n =" from Isum and expected "n =" from Esum
    # Source trace: skbias15.pas::Calculate_item_restscore_gamma1 writes
    # Isum/Esum after printCrosstab; endpoint cells are not added to those sums.
  end for
end algorithm
```

The current R package does not expose global-invariance output. The tables above
document traced source mechanics that remain relevant to item-fit maintenance,
not a current public workflow.

## Local invariance restricted MCA

Local invariance is not part of the current gRm package scope. The package
focuses on the validated GLLRM/Rasch fitting, SCREEN J, item parameters,
item fit, local dependence, DIF, global homogeneity, score effects, and
model-graph workflows.
