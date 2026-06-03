# External Data

This package scaffold does not duplicate the example runtime data.

During development, tests should reference the repository copies under:

- `path/to/DIGRAM/`
- `pascal_harness/item_parameters_report/`

Small package-local fixtures can be added here when they are needed for stable
package tests that should not depend on the full repository layout.
