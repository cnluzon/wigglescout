# Remove top loci in a GRanges column by mean of specified values

Remove top loci in a GRanges column by mean of specified values

## Usage

``` r
.remove_top_by_mean(granges, quantile, columns)
```

## Arguments

- granges:

  GRanges object. Must have numerical mcols, and names of those needs to
  match with columns.

- quantile:

  Top quantile. Must be a number between zero and 1.

- columns:

  Which columns to use. If more than one, quantile will be selected
  according to means.

## Value

A named list, fields ranges for the resulting GRanges, plus calculated
values: quantile_value, filtered, na_values.
