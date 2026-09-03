#' gRm: source-faithful R implementation of non-GUI DIGRAM slices
#'
#' `gRm` is the maintained R implementation track for selected non-GUI DIGRAM
#' computations. The implementation grows against source-shaped Pascal harnesses
#' in this repository and exposes native numeric R results.
#'
#' The package is currently experimental. The Pascal harness remains the
#' executable source-faithful reference, while `gRm` implements matching
#' computations independently in R and is validated by parity tests.
#'
#' @section Current ordinal-only restriction:
#' The current `gRm` package interface supports ordinal item and exogeneous
#' variables only. Internally, variables are represented with DIGRAM's ordinal
#' source type. Nominal and mixed variable-type behavior from the original
#' DIGRAM source is not implemented in this version. In particular,
#' multi-category nominal exogeneous variables are outside the source-faithful
#' scope of the current package and should not be interpreted as validated
#' nominal-variable analyses.
#'
#' @section Legacy DIGRAM import category-code restriction:
#' [read_digram_project()] supports the simple ordinal DIGRAM import subset
#' where category codes in `DIGRAM.imv` are contiguous one-based codes matching
#' the values in `DIGRAM.csv`. The reader follows `DIGRAM.imp` for the source
#' directory and project prefix, and treats `.imv` recursive-level marker rows
#' as structural metadata rather than variables. Historical DIGRAM projects with
#' zero-based, non-contiguous, or otherwise recoded `.imv` category values are
#' not source-faithfully implemented in this package version. The package does
#' not currently preserve separate source category codes and internal one-based
#' analysis categories.
#'
#' @section Incomplete item fit records:
#' Incomplete item fit records are a documented implementation restriction.
#' DIGRAM's Pascal source can add such records to item fit and item-restscore
#' gamma calculations only when the runtime `NincompleteRecs` state has been
#' populated. gRm does not currently have source-backed evidence for when that
#' runtime state exists, so incomplete item fit records are not synthesized from
#' missing item-response rows.
#'
#' @section Global homogeneity residuals:
#' Global homogeneity residuals and marker cells are a documented
#' implementation restriction. gRm computes the source-backed likelihood-ratio
#' test, p-value, score groups, and item means, but the item-level residual and
#' marker cells remain `NA` because the recovered source does not make the
#' runtime residual variance materialization source-backed. Public summaries
#' report those cells as not source-backed.
#'
#' @section Pascal Extended fixed-field capability:
#' DIGRAM's item-parameter report performs ICE cancellation and fixed-width
#' formatting in Pascal `Extended` arithmetic. Exact parity for those character
#' fields is platform-gated: fitted item-parameter values carry an
#' `extended_capability` list reporting the compiled long-double mantissa,
#' exponent range, storage size, x87 controls/logarithm support, and the final
#' `pascal_extended_fixed_field_supported` flag. The flag is true only for the
#' audited x87 80-bit path. On platforms such as arm64 where `long double` does
#' not provide that path, the package does not claim exact Pascal-Extended
#' fixed-field parity; native numeric fitting remains independent of this
#' presentation capability.
#'
#' @section Workflow scope:
#' Public workflows read DIGRAM projects, build item-analysis objects, specify
#' and fit Rasch/GLLRM models, screen candidate terms, run explicit diagnostic
#' functions, and inspect results through [summary()].
#' The installed package does not generate DIGRAM runtime output files or parse
#' historical DIGRAM output. Repository validation against those files belongs
#' under `validation/digram_oracle/`.
#'
#' @section Canonical workflow:
#' A `gRm_analysis` stores the encoded DIGRAM project and score-group setup.
#' From there, either specify a model manually with [gllrm()] or run [screen()]
#' to discover candidate local-dependence and DIF terms. [screen()] returns a
#' screening result, not a fitted model; pass that result to [gllrm()] to create
#' the model selected by screening. [fit()] is the boundary between model
#' specification and model diagnostics.
#'
#' ```
#' analysis <- gRm(data, items = items, exogenous = exogenous)
#'
#' manual_model <- gllrm(analysis, ld = ~ item1:item2, dif = ~ item3:group)
#' screened <- screen(analysis)
#' screened_model <- gllrm(screened)
#'
#' fitted <- fit(screened_model)
#'
#' summary(fitted)
#' item_fit(fitted)
#' local_dependence(fitted)
#' dif(fitted)
#' global_homogeneity(fitted)
#' ```
#'
#' [score_effects()] is analysis-level and can be called directly on a
#' `gRm_analysis` object. Model diagnostics require a fitted `gRm_fit` object.
#'
#' @useDynLib gRm, .registration = TRUE
#' @importFrom stats convolve
#' @importFrom utils type.convert
#' @examples
#' \dontrun{
#' analysis <- gRm(data, items = c("I1", "I2", "I3"), exogenous = "site")
#' fit0 <- fit(gllrm(analysis))
#' fit0
#' summary(fit0)
#' item_fit(fit0)
#' }
#'
"_PACKAGE"

#' Internal gRm object schemas and invariants
#'
#' Main `gRm` functions exchange S3 lists with compact public classes. This
#' contributor-facing topic records the fields that package code may rely on.
#' Fields may be added compatibly, but the identity, ordering, and nesting
#' invariants below must not be weakened silently.
#'
#' @section Project and bundle objects:
#' A `gRm_project` contains `paths`, ordered `variables`, `items`,
#' `backgrounds`, `raw_data`, `category_levels`, import metadata, and source
#' trace information. The variable tables retain source label codes, raw-data
#' positions, declared category counts, ordinal source types, names, and roles.
#' `raw_data` is the encoded, one-based source table in declared variable
#' order. A computational bundle contains `model`, `manifest`, and `data`.
#' Its item and background matrices are zero-based; its rows remain aligned
#' with project rows; and `data$status == 1L` is the conditional-likelihood
#' estimation sample. Diagnostic routines use their explicitly named record
#' policies rather than assuming that every diagnostic shares this mask.
#'
#' @section Analysis object:
#' A `gRm_analysis` owns the original `data`, its encoded `project`, ordered
#' item/exogenous names, ID declaration, and normalized `score_groups`.
#' `analysis_identity` is the canonical schema-versioned payload covering
#' roles, score cuts, category maps, row order and identity, encoded values,
#' and missing/validity state. `analysis_fingerprint` identifies that exact
#' payload. `likelihood_sample` stores its own schema, fingerprint, logical row
#' mask, indices, and count. Constructors for downstream objects must carry
#' these fields unchanged; fitting rejects an analysis mutated after creation.
#'
#' @section Model and SCREEN objects:
#' A `gRm_model` owns its `analysis`, canonical source-ordered `ld` and `dif`
#' tables, the same tables under `terms`, a `model_type` (`"rasch"` or
#' `"gllrm"`), optional originating SCREEN object, trace/status fields, and the
#' unchanged analysis fingerprint and likelihood sample. A `gRm_screen` owns
#' the same analysis identity fields plus the complete source-shaped SCREEN
#' `values`, normalized exact-command state, and final canonical terms.
#' Provisional evidence rows are retained in `values`; only finalized selected
#' rows may appear in `model_terms()` or a model created from SCREEN.
#'
#' @section Fit object:
#' A `gRm_fit` owns the exact input `analysis` and `model` (also available as
#' transition alias `spec`), the source-shaped computational `bundle`, the
#' low-level attempted `fit` state, public numeric `values` (also available as
#' transition alias `parameters`), a `gRm_convergence_state`, trace/status
#' fields, and the unchanged identity fields. The positive conditional
#' objective in `values` is DIGRAM's negative log likelihood; only
#' `logLik.gRm_fit()` changes its sign to R's log-likelihood convention.
#' Likelihood comparisons require identical canonical analyses and likelihood
#' samples and consecutive models whose complete LD/DIF term sets are nested.
#'
#' @section Public convergence state:
#' `gRm_fit$convergence` is a `gRm_convergence_state` with schema identifier,
#' final and pre-post-acceptance source flags, iteration count, source stopping
#' `report_delta`, recomputed `final_delta`, strict `tolerance`, `max_step`,
#' stable `stop_reason`, post-stop-acceptance flag, attempted-fit metadata, and
#' `report_value_source`. These deltas are maximum sufficient-count
#' discrepancies, not parameter movements. `delta` and `max_delta` are
#' transition aliases for `report_delta` and `tolerance`.
#'
#' @section Value and diagnostic result objects:
#' Numerical value helpers return source-shaped named lists or tables with
#' native-precision statistics and source-order rows. Fit values include item
#' labels, gamma parameters, thresholds, locations, ICE/MICE effects, item
#' statistics, observed score range, non-negative integer parameter count, and
#' the finite positive negative log likelihood. A diagnostic result class such
#' as `gRm_local_dependence`, `gRm_dif`, `gRm_global_homogeneity`, `gRm_cm2`, or
#' `gRm_cm3` owns `analysis`, the supplied `fit`, raw `values`, a stable `result`
#' label, method metadata, trace/status fields, and unchanged identity fields.
#' Public diagnostics compute from the supplied fit and do not silently refit.
#'
#' @section Summary objects:
#' Every package summary is a list inheriting first from
#' `summary.<input-class>`, then `summary.gRm`. It contains `title`, selected
#' `which` sections, a named `tables` list, printable headers and remarks,
#' display flags, and the input `object_class`; selected tables are also exposed
#' by their section names for programmatic access. Attributes may retain
#' statistic-family metadata such as BH thresholds or source status. Printing
#' formats this object but does not alter the native numeric tables.
#'
#' @section Validation objects and artifacts:
#' The lightweight `validation` field on analysis/model objects is status
#' metadata only; it is not oracle evidence. Repository validation lives under
#' `validation/` and persists a schema-versioned run manifest before numerical
#' work, then raw R result objects, extracted value tables, rendered reports,
#' normalized oracle and R tables, complete comparisons, summaries, timings,
#' failure subsets, and an event log. Comparisons consume persisted values, and
#' failed runs retain every artifact completed before failure. Production
#' package computations never read these artifacts or invoke Pascal.
#'
#' @name gRm-object-shapes
#' @keywords internal
NULL
