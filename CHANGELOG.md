# gRm Changelog

## gRm 0.06 - 2026-06-21

### Changed

- The `m2()` and `m3()` help pages now mark the diagnostics as experimental in
  their titles and details text.
- The package test helper now recognizes the standard `R CMD check`
  `00_pkg_src/gRm` source layout, so source-faithfulness tests can run from a
  built package tarball as well as from the development source tree.
- Added source-backed `m2()` and `m3()` fit diagnostics for fitted `gRm_fit`
  objects. The diagnostics implement the deterministic/asymptotic CM2/CM3
  portions of DIGRAM's `CM3` command, including source-order margin
  preparation, current-GLLRM expected margins, public print/summary output, and
  structured programmatic values for validation.
- M2/M3 item-exogenous margin preparation now follows the preserved DIGRAM
  `Prepare_CM3tests` `ItemBias(.i2,i1.)` convention exactly. This keeps active
  DIF rows as source zero-test rows when DIGRAM does so and suppresses the
  corresponding transposed source-artifact rows.
- The oracle validation module now parses and compares every numeric value in
  the two checked-in MULTIEXO CM3 reports against the R-side `m3()` diagnostic,
  without requiring the full validation matrix.
- `gllrm()` model construction now stores LD and DIF terms in DIGRAM source
  order regardless of formula order: local-dependence terms follow the source
  `IJinfo` upper-triangle item scan, and DIF terms follow the source `IXinfo`
  item-major, exogenous-variable order.
- The SPADI vignette now keeps the DIF example inline instead of introducing a
  separate `Check for DIF` top-level heading.
- The package `Description` field now uses clearer, shorter package wording.
- `print()` and `summary()` for `gRm_model` and `gRm_fit` now use the
  user-facing `gRm:` model labels, including expanded
  `Graphical Log-Linear Rasch Model` titles, instead of the former
  DIGRAM-specific specification and fit labels.
- `print(gRm_model)` now uses the same complete human-facing model
  specification output as `summary(gRm_model)`, while still returning the model
  object invisibly.
- `summary(gRm_model)` now has one public summary surface:
  `which = "model"` prints and returns the complete model specification,
  including the model summary plus LD and DIF term tables. The former
  independent `"ld"` and `"dif"` summary sections are no longer accepted.
- `summary(gRm_model)` now prints the model overview as descriptive labeled
  lines while retaining the stable internal table columns for programmatic
  access.
- `summary(gRm_model)` now prints LD and DIF term sections as compact
  human-readable term lists, using `None` for empty sections, while retaining
  the full LD and DIF term tables on the returned summary object.
- `summary(gRm_fit)` now prints fitted item parameters, threshold rows, and
  fitted LD/DIF term lists by default. `which = "parameters"` and
  `which = "thresholds"` select the parameter-related sections directly, while
  `print(gRm_fit)` owns compact fit status.
- `summary(gRm_fit)` now prints item-parameter and threshold table headers with
  compact user-facing labels while retaining stable numeric table columns on
  the returned summary object.
- Removed the public `item_parameters()` accessor to avoid duplicating
  `summary(gRm_fit)` output. Item parameters and thresholds are now exposed
  through `summary(fit)` and `summary(fit, which = "parameters" /
  "thresholds")`, while the internal source-facing item-parameter value
  machinery remains available to fitting, summaries, and validation.
- `print(gRm_analysis)` now shows a compact analysis status overview with the
  source label, row count, ID, item count, exogenous-variable count, and
  score-group count, leaving variable names, level counts, and the score-group
  distribution to `summary()`.
- `summary(gRm_global_homogeneity)` now uses the `gRm: Global homogeneity test`
  title and exposes cleaned public views: `which = "test"` for the global test
  statistic, `which = "score_groups"` for fitted score groups, and `which =
  "item_means"` for observed and expected item means by score group. The
  printed tables use user-facing column headings and compact numeric
  formatting, while residual and marker source-status metadata remain available
  as `attr(summary(x), "source_status")`. `print(gRm_global_homogeneity)` now
  shows a compact status overview with the test statistic and non-converged
  group-fit count.
- `summary(gRm_global_homogeneity)` now exposes the source-backed uniform
  local-dependence and uniform DIF sections for GLLRM fits with LD and/or DIF
  terms. The public tables use score-group-labelled observed/expected gamma
  columns, user-facing variable names where available, and the source chi-square
  test columns; `print(gRm_global_homogeneity)` reports non-zero uniform LD/DIF
  test counts.
- `print(gRm_screen)` now shows compact SCREEN J model-search status. The
  `summary(gRm_screen)` output now has one table-based surface with
  source-backed local-dependence, DIF, and score-effect decision rows, compact
  diagnostic numeric formatting, and a blank-headed `*` marker for rows
  selected by the SCREEN J source path. The former independent screen summary
  views for `"selected"`, `"all"`, `"score_effects"`, and `"bh"` are no longer
  public summary views; selected model terms and BH thresholds remain available
  programmatically on the returned summary object.
- Removed the lone base-pipe use from internal table helpers so package checks
  do not infer an undeclared R (>= 4.1.0) dependency.
- The `screen()` and summary documentation now explain why SCREEN J DIF summary
  rows can have blank `Chisq` or `Gamma` statistic-family cells: DIGRAM's
  source-shaped rule uses gamma when an exogenous variable is binary by category
  count or ordinal by source variable type, and chi-square otherwise.
- The SCREEN J summary marker note now distinguishes the source-backed 5%
  decision path from the global LD/DIF Benjamini-Hochberg threshold table, so
  score-effect markers are not described as if they were determined by the
  global LD/DIF BH cutoff.
- Public summary/direct-table p-value headers now mirror their displayed test
  statistic: local-dependence, DIF, and score-effect chi-square columns use
  `Pr(>Chisq)`, exact score-effect output uses `Exact Pr(>Chisq)`,
  fitted-model likelihood comparison uses `Pr(>Chisq)`, and global homogeneity
  uses `Pr(>CLR)`.
- Public print and summary surfaces now consistently start with a `gRm:` family
  header. This includes ARI tables, both `item_fit()` direct table views, the
  local-dependence diagnostic family, and SCREEN J print/summary output.
- Single-view summary methods no longer expose a public `which` argument.
  `summary(gRm_analysis)`, `summary(gRm_model)`, `summary(gRm_screen)`,
  `summary(gRm_local_dependence)`, and `summary(gRm_dif)` now have one public
  summary surface each, and stale calls that pass `which` to those methods
  error clearly. The `which` argument remains only on multi-view summaries:
  `summary(gRm_fit)` and `summary(gRm_global_homogeneity)`.
- Removed the unused `jobs` argument from `screen()`, `score_effects()`, and
  `global_homogeneity()`. Multicore `jobs` remains public only for
  `local_dependence()` and `dif()`, where it is backed by candidate-refit
  parallelization.
- `score_effects()` now returns its score-effect test table directly as a
  data-frame-like `gRm_score_effects` object, using the same direct-output
  pattern as `item_fit()`. Source-shaped backend values, selected rows, and
  Benjamini-Hochberg threshold metadata are available through `attr(result,
  "values")`, `attr(result, "selected")`, and `attr(result, "bh")`.
  Score-effects-specific `summary()` and `print()` methods are no longer part
  of the public output surface.
- The public `score_effects()` table now uses user-facing column headings and
  compact numeric formatting. Exact/simulation p-value columns are printed only
  for exact or repeated inference results where those values are present,
  avoiding all-`NA` exact columns for asymptotic output.
- Printed `score_effects()` tables now use the `gRm: Score-effect tests` title,
  matching the other human-facing diagnostic outputs.
- The generated summary-method documentation now shows the default
  `summary(gRm_fit)` call as a whole-object summary, with `which` used only
  for optional fitted-model section-specific output.
- R documentation examples now use ordinary fitting calls without explicitly
  overriding `max_step` or `max_delta`.
- `summary(gRm_analysis)` now groups its printed header into `Data`,
  `Variables`, and `Score groups` sections with blank lines between them.
- `summary(gRm_local_dependence)` now exposes a single public `"tests"` view
  for local-dependence diagnostic output. The printed table uses item names, compact
  `summary.lm()`-style numeric formatting, yes/no convergence status, delta,
  and a blank-headed `*` marker for rows at or below the
  Benjamini-Hochberg FDR threshold, while keeping selected rows and the BH
  cutoff metadata available programmatically.
- `summary(gRm_local_dependence)` now prints all local-dependence test rows
  regardless of the session `max.print` option.
- `summary(gRm_dif)` now exposes a single public `"tests"` view for DIF output.
  The printed table uses item and exogenous names, compact `summary.lm()`-style
  numeric formatting, yes/no convergence and output-stability statuses, delta,
  and a blank-headed `*` marker for rows at or below the Benjamini-Hochberg FDR
  threshold, while keeping selected rows, included DIF rows, and the BH cutoff
  metadata available programmatically. DIF summaries now also print all test
  rows regardless of the session `max.print` option.
- `print(gRm_local_dependence)` and `print(gRm_dif)` now show compact
  diagnostic status overviews with candidate counts, BH-selected counts,
  non-converged fit counts, and the applied BH threshold. DIF print output also
  reports included DIF terms and output-stable fits, leaving the complete test
  tables to `summary()`.
- `print(gRm_fit)` now shows a compact fitted-model status overview with
  convergence, iterations, delta, and log likelihood, leaving item parameters
  and fitted terms to `summary()`.
- `item_fit()` now returns the selected item fit diagnostic table directly,
  with `which = "tests"` as the default and `which = "items"` for extended
  per-item diagnostics. The unchanged source-facing item fit backend values
  remain available as `attr(result, "values")`, and the BH threshold metadata
  remains available as `attr(result, "bh")`.
- Item fit diagnostic tables now print directly from `item_fit()` using a
  shared direct-table printer. The compact test table now uses user-facing
  column headers, omits the internal item label and item fit FDR code columns
  from the public table, keeps those source-facing values in
  `attr(result, "values")$items`, and prints both available
  Benjamini-Hochberg thresholds below the table. Item fit specific `print()`
  and `summary()` methods are no longer part of the public output surface.
- `item_fit(..., which = "items")` now returns and prints the extended
  per-item outfit/infit diagnostic summary with user-facing column headers,
  including `Item`, `Outfit total`, and `Infit ratio`, while keeping the
  original source-shaped extended summary columns available through
  `attr(result, "values")$extended$summaries`.
- `item_fit()` direct table results no longer retain the full originating
  `gRm_analysis` and `gRm_fit` objects as attributes; source-facing values,
  BH metadata, view metadata, and the call remain available.
- `summary(gRm_model)` now prints the model overview subsection header as
  `Model`.
- `summary(gRm_fit)` now prints the fit overview subsection header as `Fit`.
- `summary(gRm_analysis)` printing now avoids repeating the data table after
  the analysis header while keeping the `"data"` table available on the
  returned summary object. The header now reports the captured data argument
  label instead of the generic project name fallback, including when project
  objects are converted back to analyses, and shows item/exogenous project
  level counts. Analysis summaries also expose and print a compact
  score-group distribution derived from the constructed project score data and
  analysis cutpoints using the source report's score/count/percent/cumulative
  shape. `summary(gRm_analysis)` now has one output surface: requesting
  `"data"` or `"score_groups"` validates the section name but returns the same
  complete analysis summary.
- Clarified the `gRm()` help page that `item_levels` and `exogenous_levels`
  are construction-time encoding controls whose resolved order is baked into
  encoded project data and category counts for later analyses.
- `summary(gRm_analysis)` printing now starts with a compact
  `gRm: Graphical Log-Linear Rasch Model analysis` header showing the data
  name, row count, item and exogenous roles, ID variable, and score groups
  before the selected table.
- Clarified the `read_digram_project()` roxygen/manual documentation that `id`
  is optional and defaults to the first CSV column when omitted.
- Renamed the GLLRM implementation files, native backend, tests, summary
  section names, and documentation to use plain GLLRM terminology. Fitted LD/DIF
  terms are now described as included terms, and the maintained fitting path is
  simply the GLLRM path.

## gRm 0.04 - 2026-06-15

### Changed

- DIGRAM `.cmd` files are no longer parsed or exposed by the R package. The
  package keeps exact/asymptotic/repeated behavior under explicit R API control
  and ignores command-file lines in `DIGRAM.imp` import metadata.
- Removed unused base global-homogeneity expected-variance helper paths from
  the R package. Base global-homogeneity residual, marker, and expected-variance
  cells remain unreported unless a future source trace supports them.
- Removed unused DIF conditional-test helpers so the R package has a
  single DIF reporting path: source-backed candidate GLLRM refits with
  item-screening gamma provenance.
- Removed unused LD conditional-test helpers so the R package has a
  single LD reporting path: source-backed candidate GLLRM refits with
  item-screening WPG gamma provenance.
- Removed unused native validation/debug probes from the compiled package while
  keeping the production native entry points and the documented SCREEN J
  source-random compatibility alias.
- Removed retired score-spec constructor scaffolding from the R package.
  Score grouping remains controlled by `score_cuts = "auto"` or integer-like
  cut vectors.
- Removed the pass-through role resolver from `read_digram_project()`. The
  reader still uses explicit caller-supplied item, exogenous, and ID roles.
- Removed an unused generic variable-name validator from project input
  construction; existing constructor-local checks remain the input contract.
- Removed a dead item-analysis summary helper cluster; current summary output
  continues to use `make_gRm_summary()` and focused public table builders.
- Removed dead local invalid-vector storage from source bundle construction
  while preserving the zero-valued invalid fields and manifest counters.
- Removed the unused zero label from generated DIGRAM `.imv` category labels;
  generated label names now match the supported one-based category subset.
- Removed an unused base-DIF result preallocation block; DIF test rows continue
  to be returned by the source-ordered candidate dispatch path.
- Removed dead metadata from the GLLRM expected-margin native payload and
  replaced an unused native row-count assignment with an explicit `item_gamma`
  dimension check.
- Removed an unreferenced native extended item-parameter gamma helper; the
  source-shaped item-parameter arithmetic path is unchanged.
- Consolidated native top-ICE cancellation formatting behind the shared
  source-traced helper; the Pascal expression order and signed-zero behavior are
  unchanged.
- Removed the uncalled R wrapper for SCREEN J source-random draws while keeping
  both native source-random aliases for validation compatibility.
- Updated a stale SCREEN J conditional-native probe comment to describe
  optimized-route parity coverage instead of a closed native gate.
- Updated `ALGORITHMS.md` to reflect the current source-faithful package
  algorithms and native boundaries after the simplification pass.

### Fixed

- SCREEN J asymptotic output now matches the MULTIEXO DIGRAM oracle for the
  source-shaped weighted partial gamma average and post-screen spurious DIF
  summaries. The fixes follow the Pascal source conventions for
  `Largest_possible_score`, `No tests` rows in Benjamini-Hochberg thresholds,
  and the printed average absolute partial gamma calculation.
- SCREEN J exact and repeated MULTIEXO validation now uses the parity-gated
  native conditional route with double-precision source p-values, and exact
  post-screen stepwise DIF keeps DIGRAM's no-test candidate retention behavior.
  The focused MULTIEXO SCREEN J validator now covers asymptotic, exact, and
  repeated oracle reports.
- SCREEN J native exact chi-square simulation now accumulates generated
  chi-square cells in Pascal `SkStat.RCCHI` row-major order, preserving
  source-tied `CHI >= CHITOT` comparisons. This fixes the BIRT_disinh exact and
  repeated SCREEN J spurious DIF p-value mismatch without disabling the native
  path.

## gRm 0.04 - 2026-06-14

### Added

- Saved DIGRAM import bundles now write item/background columns from the
  encoded one-based `project$raw_data`, keeping generated `DIGRAM.csv` and
  `DIGRAM.imv` files internally consistent while preserving the original ID
  column.
- Local-independence and no-DIF diagnostics now expose source-shaped
  item-screening gamma provenance columns (`wpg_gamma_source` and
  `gamma_source`) alongside the existing candidate-fit statistics.
- GLLRM fits and candidate diagnostic tables now carry
  explanatory stop-reason and checkpoint metadata. Summary output reports
  non-converged candidate counts and stop-reason counts when available.
- Public item fit values now expose `incomplete_records_used = FALSE` and
  `incomplete_records_status = "not_source_backed"` to make the current
  source-gated incomplete-record restriction machine-readable.

### Changed

- The duplicated source-shaped 30-pass table-standardization loop now lives in
  a shared helper with Pascal source comments, while preserving the existing
  wrapper contracts and fixed-pass row/column order.
- The conditional SCREEN J implementation was split into smaller source-traced
  helpers. The registered conditional native route is now parity-gated and
  remains on the R source-shaped fallback when the probe detects native/R
  drift.
- `rbind_fill()` now preserves existing integer, double, logical, and character
  column types when filling missing columns.
- Summary printing now shows the Benjamini-Hochberg header only for summaries
  whose tests are actually BH-backed.

### Fixed

- `anova.gRm_fit()` now leaves signed invalid-order likelihood-ratio
  diagnostics visible but sets their p-values to `NA` with a warning.
- Native boundary code now validates GLLRM expected-margin matrix inputs before
  duplicating them and explicitly protects scalar SCREEN J exact-test
  attributes.
- Generated native build products in `gRm/src/` are cleaned after source
  installs, preventing stale `.o`, `.so`, and `symbols.rds` files from
  polluting direct `R CMD check gRm` runs.

## gRm 0.03 - 2026-06-14

### Added

- Summary output now remarks on non-converged candidate fits whenever an
  included summary table has `converged == FALSE`. The detailed tables continue
  to expose the machine-readable `converged` and `delta` columns.

### Fixed

- Local-dependence diagnostics now match DIGRAM's CHECK LID reporting
  for non-converged IJ candidate fits by reporting the first post-50
  checkpoint used by the source output.
- Item-parameter and GLLRM input statistics now match DIGRAM's
  `Nresponses` convention for useful incomplete rows, including the
  LD/DIF missing-item case observed in the SocCog oracle reports.

## gRm 0.02 - 2026-06-13

### Changed

- GLLRM LD and DIF candidate checks now use an internal lightweight
  source-ordered candidate fit path instead of constructing a full public
  `gRm_fit` object for every tested candidate. The source-faithful
  fitting loop, term ordering, likelihood-ratio statistic, convergence fields,
  and degrees-of-freedom conventions are preserved, while unused item-parameter
  output construction is skipped for candidate diagnostics.
- GLLRM expected-margin calculation now has a source-traced native C++
  backend for the full fitted item, LD, and DIF margin pass. The fitting loop,
  update rules, gauge adjustment, likelihood calculation, candidate ordering,
  and public API remain in R; the original R implementation is retained as the
  reference path and focused tests compare the native margins and candidate
  refits against it. Targeted benchmarking showed about `15x` speed-up for the
  EmoReg LD candidate and about `35.5x` for the SocCog LD
  candidate used in the performance plan.

## gRm 0.02 - 2026-06-12 16:29:53 CEST

### Added

- Added `ari()` for fitted gRm models. The function returns the DIGRAM ARI
  item-by-total-score table as an R data frame with native-precision numeric
  values instead of writing `Ari_dot.csv` or `Ari_comma.csv`.
- Implemented ARI observed category probabilities, expected category
  probabilities, means, variances, and standardized residuals from the fitted
  model using the source score window `1..highest_possible_score-1`.
- Added ARI validation support for all available `Ari_dot.csv` oracle files
  under `validation/BIRT`. The ARI-only validator now compares all 13 files and
  all 291,434 semantic numeric values with zero mismatches.
- Added `plot()` support for `gRm_ari` objects. The plot method uses ggplot2
  to reproduce the DIGRAM/SAS ARI mean-curve transformation from the in-memory
  `ari()` table, including source-style score interval collapsing, weighted
  observed and expected means, and normal critical-value confidence bands.

### Changed

- ARI plots now use ordinary ggplot2 facet-strip item titles and
  `theme_minimal()` instead of drawing item names inside each panel.

## gRm 0.01 - 2026-06-12 15:22:52 CEST

### Fixed

- LD/DIF diagnostic refits no longer rebuild model terms by pasting item
  and exogenous names into formulas. Candidate terms now preserve the
  canonical parsed LD/DIF term tables directly, so non-syntactic names and names
  containing `:` remain valid while the source ordering of candidate terms is
  unchanged.
- Item fit output no longer depends on R's partial `$` matching to resolve
  `label` to `label_code`. All item fit label output now reads the explicit
  DIGRAM-style `label_code` metadata column, avoiding warnings under
  `warnPartialMatchDollar = TRUE` and preventing future metadata columns from
  silently changing labels.
- The internal GLLRM likelihood comment now correctly states that the
  fitted state stores the DIGRAM/source-style negative conditional log
  likelihood. A focused test guards that `logLik()` returns the standard R sign
  by negating the stored source value.
- The GLLRM component-enumeration guard now reports itself as an R
  implementation guard rather than a source limit. This keeps source-backed
  behavior distinct from package-side resource protection.
- Summary output for item parameters and global homogeneity no longer exposes
  duplicate aliases that returned identical tables. Item-parameter summaries
  now use coefficient-oriented names, while global homogeneity keeps the
  diagnostic `tests`, `groups`, and `items` views.
- The internal summary-table normalization helper was renamed away from the
  retired `tidy()` API terminology. This is an internal cleanup only and does
  not change summary table contents.
- Global-homogeneity degrees-of-freedom flooring was not changed. The
  code now carries an explicit guardrail comment explaining that the convention
  is mathematically unusual and must not be altered without source or validator
  evidence.

## gRm 0.01 - 2026-06-10 12:10:23 CEST

### Fixed

- Improved internal source-trace documentation for calculation-heavy
  GLLRM, DIF, LD, and screen J helpers. The added comments identify the
  corresponding Pascal procedures and distinguish source-faithful calculations
  from R-only implementation caches. No numerical behavior was changed.
- Internal R source files were reorganized to better match their
  responsibilities. The source bundle, DIGRAM import I/O, project construction,
  summary table helpers, source score-group helpers, and model-graph plotting
  code now live in files named after their current roles. Public API names,
  numerical algorithms, and source-faithful output tables were not changed.
- The stale internal local-invariance implementation was removed after a
  reference scan confirmed it was outside the current package and validator
  call paths. Local invariance remains outside the current gRm package scope.
- Exact and repeated SCREEN J were optimized without changing the source
  algorithm. The native exact kernel now reuses per-slice work buffers,
  precomputes deterministic slice metadata, and computes generated-table
  chi/gamma statistics with less repeated allocation. Source-faithful native
  routing is also enabled for compatible partial item-pair exact tests after a
  cached native/R parity probe. Generic conditional bias tests remain on the R
  source-shaped path because the native shortcut can drift on a borderline
  exact chi-square comparison. Focused tests preserve exact p-values, repeated
  stopping, draw counts, final source seeds, and selected model terms. The
  targeted SCREEN J exact and repeated validators now pass for the BIRT oracle
  projects. Full validation was intentionally deferred because this change is
  limited to compatible exact SCREEN J simulation paths.
- Native SCREEN J chi-square arithmetic now follows the Pascal `RCCHI`
  operation order, computing each cell contribution as residual times
  residual-over-expected. This avoids borderline floating comparisons from
  drifting away from the source-shaped R path. Repeated chi-only exact
  inference now stays on the R source-shaped path for sequential stopping,
  while the fused SCREEN J exact chi/gamma path remains native.
- Fitted-model comparison now has standard R likelihood methods. `logLik()`
  on a `gRm_fit` converts the source-stored DIGRAM negative log likelihood to
  R's usual log-likelihood sign, and `anova()` compares two or more fitted
  models with a likelihood-ratio table using `Model Df`, `-2 logLik`, `Df`,
  `Chisq`, and `Pr(>Chisq)`. The likelihood-ratio statistic remains
  source-faithful, while the p-value uses `stats::pchisq()` so the output
  follows gRm's public p-value heading convention.
- Bundle construction grew `complete_item_scores` and source bundle
  `model_rows` inside loops. These paths now preallocate or collect rows before
  binding, preserving the source-shaped score-window and model-file semantics
  while avoiding repeated object growth.
- Superseded non-context LD/DIF helper implementations remained in the
  package beside the current cached context paths, duplicating source-shaped
  margin, expected-count, and likelihood logic. The unused non-context helpers
  were removed so the candidate diagnostic implementation has a single maintained
  runtime path for each calculation.
- Marginal gamma self-pair cells were exposed as numeric zero in the R-facing
  gamma matrix. The Pascal report renders these diagonal cells as not
  applicable, so the R numeric result now uses `NA_real_` for the diagonal
  while preserving source-faithful off-diagonal gamma calculations.
- Global-homogeneity uniform-DIF preparation rebuilt the same source
  score-group lookup, complete-row set, and score/exogenous grouping for each
  DIF term. The DIF path now has an all-spec score-group-table pass
  analogous to the existing LD path, so shared preparation is performed once
  while the source-faithful observed/expected table construction and summary
  statistics remain unchanged.

## gRm 0.01 - 2026-06-09 18:02:57 CEST

### Fixed

- Base local-dependence diagnostics used declared item category dimensions for
  candidate LD degrees of freedom. The DIGRAM source derives the IJ degrees of
  freedom from the observed nonzero item margins for the candidate pair. Base
  LD tests now use the same observed-margin helper already used by GLLRM
  LD checks, so sparse item-pair tables get source-faithful df and p-values.
- SCREEN J exact/repeated inference passed the DIGRAM repeated-Monte-Carlo
  command state to the partial-gamma branch but not to the chi-square-only
  exact branch. The Pascal source sets `SEQUENTIAL`, `SEQ_P0`, `SEQ_B`, and
  `seq_limit` before `SKbias3.XYZ_bias_ANALYSE` evaluates simulated chi
  results. The chi exact helper, native bridge, and exogenous SCREEN J call
  site now carry the same source command state.
- Base Rasch `item_fit()` recomputed a fresh default fit from the project
  instead of diagnosing the `gRm_fit` object supplied by the caller. DIGRAM
  item fit diagnostics run against the current fitted state. The base item fit
  path now uses the supplied fit object's bundle and fitted Rasch state
  directly, so non-default fitting controls and non-converged supplied fits are
  reflected in the diagnostic values.
- Public Rasch and GLLRM fitting could return a fitted object with
  `converged = TRUE` even when the source-shaped bundle had no valid complete
  response patterns inside the DIGRAM score window. The source report path
  exits these sections when `Nvalid < 1`, so `fit()` now fails before
  estimation with a clear non-estimable-data error instead of producing zero
  likelihood and invalid parameter-count metadata.
- Automatic score-cut construction accepted invalid generated cuts, including
  missing, decreasing, and duplicate cut values from data with no complete
  item patterns, only lower-boundary scores, or too little score range. The
  constructor now validates automatic cuts before storing them and requires
  them to define at least two usable source score groups, matching the score
  grouping convention used by the source global-homogeneity path. Item
  selection now reports a clear error when no complete item response pattern is
  available for deriving score groups.
- The retired graphical DIGRAM `SCREEN` implementation still shipped beside
  the validated SCREEN J API and compiled helper kernels at runtime with
  `Rcpp::sourceCpp()`. Graphical SCREEN is outside the current package scope,
  so the old implementation, its generated man pages, its stale tests, and the
  unused `Rcpp` import were removed. The public `screen()` API continues to run
  the source-shaped SCREEN J implementation.
- Internal global-invariance/global-DIF report code remained in the package
  even though this functionality is outside the currently validated gRm scope
  and has no oracle target. The implementation file, generated man
  pages, and package-facing global-homogeneity wording that implied global
  invariance support were removed. Global homogeneity remains implemented.
- Exogenous exact score-effect simulations still called the historical
  `screen_rc_chi_square()` helper after the retired graphical SCREEN file was
  removed. The helper name is now retained in the SCREEN J support code as a
  thin wrapper around the source-shaped Pearson chi-square calculation, so
  exact and repeated exogenous-selection tests can run without restoring the
  retired graphical SCREEN implementation.
- Item fit incomplete-response contribution code exists internally but is not
  enabled in the public item fit path because the package has no source-backed
  way to determine when DIGRAM's runtime `NincompleteRecs` list has been
  populated. The package and `item_fit()` documentation now mark incomplete
  item fit records as an implementation restriction instead of implying that
  missing item-response rows are automatically incorporated.
- Global-homogeneity item-level residual and marker cells remain intentionally
  `NA` because the recovered Pascal source does not source-back the hidden
  runtime residual variance materialization. The public
  `global_homogeneity()` documentation and package overview now state this
  restriction explicitly and point users to the `"not_source_backed"` status in
  summaries.
- Legacy DIGRAM `.imv` import documented support only for contiguous
  one-based category codes, but the parser accepted non-contiguous, duplicated,
  zero, negative, fractional, or non-numeric declarations by reducing them to
  `raw_max`. The `.imv` reader now enforces the documented source-faithful
  import subset before constructing variable metadata.
- GLLRM LD/DIF formula parsing used formatted `terms()` labels and split them
  on `:`, which rejected valid backticked non-syntactic variable names such as
  `item one`. The parser now walks the formula language object directly,
  supports backticked item and exogeneous names, and keeps rejecting syntax
  outside `:` interactions joined by `+`.
- Public integer-like validators for `jobs`, score cuts, and explicit score
  group cuts could coerce oversized numeric values with `as.integer()` before
  checking range, producing base R missing-value errors. These paths now share
  a pre-coercion validator that checks type, length, missingness, finiteness,
  integer-likeness, and R integer range before conversion.
- Internal implementation helpers still generated installed help pages even
  though they were not part of the user-facing gRm API. Helper documentation
  comments were moved out of roxygen, stale internal `.Rd` files were removed,
  and a focused guard test now keeps the installed help surface limited to the
  public workflow functions and package/method overview pages. The same
  documentation pass corrected S3 generic usage entries and a malformed
  non-syntactic-name example caught by package documentation checks.
- The namespace still registered stale S3 methods for retired class names
  `gRm_item_analysis`, `gRm_gllrm_spec`, and `gRm_gllrm_fit`. Newly-created
  objects already used the simplified class surface, so the old method aliases
  and compatibility branches were removed and a focused API test now checks
  both object class vectors and generated namespace registrations.
- Package-maintainer documentation still described removed files and helper
  APIs, including retired check/provenance orchestration, the old graphical
  SCREEN implementation, global-invariance internals, and obsolete vignette
  paths. `ALGORITHMS.md` and `docs/r-package/GRM_R_IMPLEMENTATION_FILES.md` now describe
  the current source-faithful package boundary: explicit public workflow
  functions, `summary()`/compact `print()` output, current implementation
  files, and retired functionality as out of scope.
- Generated native artifacts in `gRm/src/` and a Markdown implementation note
  under `gRm/R/` polluted source-package checks and production-copy workflows.
  The ignored native build products were removed from the working tree, the
  implementation note was moved out of the package source tree, build-ignore
  rules now exclude generated native products from source bundles, and a focused
  API guard now fails if Markdown files reappear in the R source directory.
- Repeated internal association-gamma helpers used separate implementations for
  local dependence, DIF, SCREEN J, marginal gamma, item fit gamma, and
  GLLRM standardized tables. These helpers are now routed through one
  source-traced RC gamma core that mirrors the Pascal all-cell
  `AIJ`/`DIJ`/`P`/`Q` convention while preserving each caller's return shape,
  inference wrapper, and public numeric output.
- Package-check-only test harness failures were corrected without changing
  package behavior: tiny toy-analysis fixtures that cannot use automatic source
  score groups now pass explicit score cuts, and source-tree scaffold checks now
  skip direct source-file assertions in installed-package test contexts.

## gRm 0.01 - 2026-06-08 16:19:13 CEST

### Fixed

- Local-dependence and DIF diagnostics used persistent in-session caches with
  weak keys. Different projects with the same dimensions and raw-data sum could
  receive stale p-values, BH thresholds, and suggested terms from an earlier
  analysis. The caches were removed so these diagnostics now recompute from the
  current project/model state, matching the source behavior.
- Item fit diagnostics synthesized source-eligible incomplete response records
  from the current data and included them in outfit, infit, and item-restscore
  gamma calculations. DIGRAM only includes those records when the runtime state
  has populated `NincompleteRecs`, and the corresponding source report prints
  explicit incomplete-record messages. In the absence of that source-backed
  runtime evidence, item fit calculations now use an empty incomplete-record
  set so the compact and extended item fit values match the recovered source
  convention.
- Automatic score-cut selection used all complete score counts as the
  denominator, including score `0` and the maximum score. The Pascal source uses
  only non-extreme scores for this median-like cut. The denominator now excludes
  extreme scores, so `score_cuts = "auto"` follows the source convention.
- `global_homogeneity(score_cuts = ...)` accepted malformed explicit overrides
  by truncating non-integer cuts, clipping out-of-range cuts later in the
  calculation, and falling back to default cuts when too few cuts were supplied.
  Explicit score-cut overrides are now validated before the diagnostic runs:
  they must be non-missing integer-like values, strictly increasing, inside the
  possible score range, and define at least two usable source score groups.
- The older graphical `screen_values.R` path still used base R `pchisq()` and
  `pnorm()` tail probabilities. Other validated diagnostics use Pascal-shaped
  `PFCHI` and `PNORMAL` ports. Screen table p-values now use `source_pfchi()`
  and `source_tail_norm()` so this path follows the same source convention.
- `score_effects()` reported the exogenous-selection score distribution from
  complete item-score rows, which included rows with missing exogenous values.
  The Pascal source increments and prints the temporary `scoredistribution0`
  table only for complete item-plus-exogenous records, with its missing bucket
  remaining zero on this path. `score_effects()` now uses that source
  convention and the score-distribution helper still returns an explicit empty
  summary when no usable known scores exist.
- `build_item_parameters_bundle()` derived the fitted upper score window from
  rows with both complete item responses and complete exogenous values. If the
  highest complete item score occurred on a row with missing exogenous data, the
  fitted score window was too small and the row was not counted as a missing
  background case. The bundle now derives the score window from complete item
  rows before exogenous missingness is applied, and records the complete
  item-plus-background count separately in the manifest.
- `score_effects()` exposed `score_cap` as a public argument even though the
  DIGRAM source uses a fixed compile-time `MAXDIM` convention and collapses
  scores above `maxdim - 1`. The public argument was removed and the package now
  uses an internal source constant of `56`, so callers cannot request
  non-DIGRAM score-collapsing behavior.
- Exact inference controls (`nsim`, `seed`, `critlevel`, and `risk`) were
  normalized by taking the first element and calling `as.integer()`. Fractional
  values could therefore be silently truncated and vector-valued inputs could
  silently drop all but their first value. The exact command-state helper now
  validates scalar integer-like inputs before conversion, while preserving the
  DIGRAM source convention that `nsim = 0` requests the default of `1000`.
- Legacy DIGRAM import partially ignored source import structure. The reader
  looked for `.csv` and `.imv` files from the R-supplied prefix rather than the
  project path/name stored in `DIGRAM.imp`, and treated `.imv` recursive-level
  marker rows beginning with `<` as malformed variables. The import reader now
  follows the `.imp` project path and project name, uses that project name for
  the public analysis identity, and keeps `.imv` marker rows as import metadata
  instead of variable definitions.
- Several public `print()` and `summary()` S3 methods had duplicate
  implementations in different R files. The duplicates were not identical:
  some represented older print output and others represented the retired
  `summary(details = ..., detail = ...)` accessor interface. The stale copies
  were removed so the documented `api-summary.R` methods are the sole public
  `print()` and `summary()` implementations.
- The package-facing output contract was simplified to `summary()` and
  compact `print()` methods. The retired `tidy()`, `glance()`, `details()`,
  and `detail_names()` accessors are now treated as absent from both exports
  and the package namespace.
- Global-homogeneity item residual and marker cells are intentionally `NA`
  because the recovered Pascal source does not source-back the exact runtime
  residual materialization. That limitation was visible only in item-level
  flags. `summary(global_homogeneity_result, which = "summary")` now includes
  item residual and marker source-status fields so the public summary interface
  clearly reports these cells as `not_source_backed`.
- Fitting controls were not validated consistently before entering the Rasch
  and GLLRM fitting paths. Negative, zero, fractional, `NA`, infinite, or
  vector-valued `max_step`/`max_delta` inputs could either run, skip fitting, or
  fail later with base R errors depending on the path. A shared public validator
  now requires `max_step` to be a single positive integer-like value and
  `max_delta` to be a single positive finite number for `fit()`,
  `local_dependence()`, `dif()`, and `global_homogeneity()`.
