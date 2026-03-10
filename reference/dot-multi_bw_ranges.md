# Intersect a list of bigWig files with a GRanges object

Build a binned-scored GRanges object from a list of bigWig files. The
aggregating function per locus can be min, max, sd, mean.

## Usage

``` r
.multi_bw_ranges(
  bwfiles,
  labels,
  granges,
  per_locus_stat = "mean",
  selection = NULL,
  remove_top = 0,
  default_na = NA_real_,
  scaling = "none"
)
```

## Arguments

- bwfiles:

  Path or array of paths to the bigWig files to be summarized.

- labels:

  List of names to give to the mcols of the returned GRanges object. If
  NULL, file names are used.

- granges:

  GRanges object to intersect

- per_locus_stat:

  Aggregating function (per locus). Mean by default. Choices: min, max,
  sd, mean.

- selection:

  A GRanges object to restrict binning to a certain set of intervals.
  This is useful for debugging and improving performance of locus
  specific analyses.

- remove_top:

  Return range 0-(1-remove_top). By default returns the whole
  distribution (remove_top == 0).

- default_na:

  Default value for missing values

- scaling:

  If none, no operation is performed (default). If relative, values are
  divided by global mean (1x genome coverage).

## Value

A sorted GRanges object.
