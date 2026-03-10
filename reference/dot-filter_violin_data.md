# Filter violin plot data with quantile threshold

Removes top quantile loci from the distribution. Values are removed by
mean rather than maximum.

## Usage

``` r
.filter_violin_data(gr, remove_top, columns)
```

## Arguments

- gr:

  GRanges object

- remove_top:

  (1-remove_top) quantile will be removed.

- columns:

  mcols in GRanges to use for selection.

## Value

A named list with ranges = GRanges object and stats a named list with
some calculated values.
