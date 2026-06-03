# gRm Test Policy

Package tests must be self-contained or use small fixtures committed under
`gRm/tests/testthat`.

The historical top-level `examples/legacy_*` fixtures are obsolete and must not be
used by package tests. DIGRAM oracle validation belongs under
`validation/digram_oracle/`.

Tests must not use historical output files under `data/`.
Production R code must not read DIGRAM output text.
