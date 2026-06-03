# R Implementation Files

This directory contains the active native R implementation for the installed
`gRm` package.

The package boundary is numeric and R-native:

- public constructors and workflow functions build item-analysis objects,
  specify GLLRM models, fit models, screen candidate terms, and run checks;
- public output methods expose computed values through `summary()`, `tidy()`,
  `glance()`, `detail_names()`, and `details()`;
- source-trace, status, warning, and unmodeled-field helpers describe what is
  implemented and which source-backed values remain intentionally unavailable.

DIGRAM runtime-output parsing, output-file comparison, and validation target
selection do not belong in this directory. Repository validation code lives
under `validation/digram_oracle/` and uses only the public numeric package API
when it needs `gRm` results.

Main file groups:

- `api-constructors.R`, `api-model-spec.R`, `api-fit.R`, `api-screen.R`,
  `api-check.R`, `api-output.R`, and `api-provenance.R`: public item-analysis,
  model, diagnostic, accessor, and metadata APIs;
- `data_entry_points.R`: internal data-loading machinery used by `gRm()`
  and `read_digram_project()`;
- `read_bundle.R`: build/read `GLLRMdata.txt`, `model.tsv`, and
  `manifest.tsv`;
- `rasch_base_fit.R`: base Rasch IPF/CML estimation;
- `*_values.R` files: source-shaped numeric value computation;
- support files for model terms, scoring, source traces, warnings, status, and
  unmodeled fields.

General interface documentation:

- `../../docs/R_DIGRAM_DUAL_OUTPUT_API_DESIGN.md` records the current numeric
  package boundary;
- `../../docs/R_TIDY_GLANCE_PROPOSAL.md` defines the accessor contract;
- `../vignettes/gRm-api-example-emo.Rmd` demonstrates the implemented API on
  the example analysis without package-level runtime-output
  validation.

Do not add production functions here before adding and running the matching
failing `testthat` tests.

Documentation rules:

- add roxygen documentation before every new public function;
- use `@export` only for user-facing workflow functions;
- use `@keywords internal` for implementation helpers and reference helpers;
- document list return shapes explicitly for public functions;
- regenerate package documentation with
  `roxygen2::roxygenise('gRm')` from the repository root.
