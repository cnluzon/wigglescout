# GRanges left join operation

Perform a left_bind operation on a GRanges list, appending scores to
mcols.

## Usage

``` r
.granges_left_join(grlist, labels, granges = NULL)
```

## Arguments

- grlist:

  A list of GRanges objects that have all the same fields.

- labels:

  Vector of names for the score columns.

- granges:

  Optional granges with name field on it

## Value

A Sorted GRanges object with all the columns.

## Details

This assumes that you are trying to join things that match (bins
generated from the same parameters, BED intersections from the same
BED.)
