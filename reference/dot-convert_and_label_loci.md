# Internal processing of loci sets and labels.

This allows to process equally lists of GRanges objects or BED files.
Note that GRanges objects need to be lists, as c(GRanges) concatenates
elements in single GRanges object

## Usage

``` r
.convert_and_label_loci(loci_sets, labels)
```

## Arguments

- loci_sets:

  One or more BED files or GRanges. These need to be either all BED
  files or all GRanges objects.

- labels:

  One or more labels. Lengths need to match loci_sets or be NULL, in
  which case loci_sets are assumed to be BED files.

## Value

a named list with ranges and labels.
