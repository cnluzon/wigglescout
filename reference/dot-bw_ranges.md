# Score a GRanges object against a BigWig file

Build a scored GRanges object from a GRanges object and a single BigWig
file. The aggregating function can be min, max, sd, mean.

## Usage

``` r
.bw_ranges(
  bwfile,
  granges,
  per_locus_stat = "mean",
  selection = NULL,
  default_na = NA_real_,
  scaling = "none"
)
```

## Arguments

- bwfile:

  Path to a single BigWig file to be summarized.

- granges:

  GRanges object file to be summarized.

- per_locus_stat:

  Aggregating function (per locus). Mean by default. Choices: min, max,
  sd, mean.

- selection:

  A GRanges object to restrict binning to a certain set of intervals.
  This is useful for debugging and improving performance of locus
  specific analyses.

- default_na:

  Default value for missing values

- scaling:

  If none, no operation is performed (default). If relative, values are
  divided by global mean (1x genome coverage).

## Value

GRanges with column score.
