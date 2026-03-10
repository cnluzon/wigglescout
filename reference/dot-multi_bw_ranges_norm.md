# Intersect a list of bigWig files with a GRanges object (with background)

Intersect a list of bigWig files with a GRanges object and normalize
values with a set of background bigWig files.

## Usage

``` r
.multi_bw_ranges_norm(
  bwfilelist,
  bg_bwfilelist,
  labels,
  granges,
  per_locus_stat = "mean",
  selection = NULL,
  norm_func = identity,
  remove_top = 0,
  default_na = NA_real_,
  scaling = "none"
)
```

## Arguments

- bwfilelist:

  BigWig file list to be summarized.

- bg_bwfilelist:

  Background BigWig files.

- labels:

  List of names to give to the mcols of the returned GRanges object. If
  NULL, file names are used.

- granges:

  GRanges object to intersect

- per_locus_stat:

  Aggregating function (per locus). Mean by default. Choices: min, max,
  sd, mean.

- selection:

  A GRanges object to restrict analysis to.

- norm_func:

  Function to apply after bw/bg_bw.

- remove_top:

  Return range 0-(1-remove_top). By default returns the whole
  distribution (remove_top == 0).

- default_na:

  Default value for missing values

- scaling:

  If none, no operation is performed (default). If relative, values are
  divided by global mean (1x genome coverage).

## Value

a sorted GRanges object
