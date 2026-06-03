# IncompleteItemFits Gap

The R ItemFits implementation is not complete until the DIGRAM
`IncompleteItemFits` branch has been implemented and validated.

## Why This Is Still Missing

The supplied example `ItemFits-extended.txt` runtime output does not contain the
source line indicating that persons with incomplete item responses were included.
Therefore the example ItemFits example cannot validate the branch in
`skbias12a.pas::IncompleteItemfits` / `skbias15.pas::Calculate_residuals_and_item_fits`.

For the current example report, the source-faithful behavior is to avoid
synthesizing incomplete item-fit records merely because rows with missing item
responses exist in the bundle.

## Source References

- `source/PAS_scd/DGRirtD.pas` around lines 2656-2809 decides whether
  incomplete records are carried into the runtime fit-statistics path and prints
  the source message that incomplete response patterns were included.
- `source/PAS_skunits/skbias15.pas` around lines 3483-3576 calls
  `IncompleteItemfits` from the extended outfit/infit branch and prints
  `<n> persons with incomplete responses were included during calculation of
  Outfits and Infits`.
- `source/PAS_skunits/skbias12a.pas::IncompleteItemfits` around lines
  9324-9473 contains the contribution formulas that still need an independent R
  port and runtime validation.

## What Is Needed

To complete the R ItemFits implementation, we need a DIGRAM runtime example
where `NincompleteRecs > 0` and the original report explicitly includes the
incomplete item-response branch. That example should include:

- the input data and DIGRAM project files;
- compact and extended ItemFits output from original DIGRAM;
- evidence in the output that incomplete item responses were included.

Once available, implement the branch in independent R code from the Pascal
source formulas and test it against both the Pascal harness and the original
DIGRAM output. Do not infer or hard-code the branch from example missingness alone.
