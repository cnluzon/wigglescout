# Convert a GRanges object with coverage values to estimated counts

This method uses estimate_read_counts across a GRanges object and
converts the values to read counts.

## Usage

``` r
gr_coverage_to_read_counts(gr, mcol_names, fraglen)
```

## Arguments

- gr:

  GRanges object

- mcol_names:

  Columns to calculate on

- fraglen:

  Fragment length to use for the estimate

## Value

A GRanges object with count values.
