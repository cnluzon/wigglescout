# Summarize a intersect_bw_and_granges matrix

Compute averages and standard error values for a matrix returned by
intersect_bw_and_granges.

## Usage

``` r
.summarize_matrix(matrix, label)
```

## Arguments

- matrix:

  A matrix returned by .intersect_bw_and_granges.

- label:

  Label for the sample.

## Value

A data frame with summarized values, stderr and medians, plus label.
