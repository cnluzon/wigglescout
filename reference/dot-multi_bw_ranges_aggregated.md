# Intersect a list of bigWig files with a GRanges object and aggregate by name

Intersect a list of bigWig files with a GRanges object and aggregate by
name

## Usage

``` r
.multi_bw_ranges_aggregated(
  bwfiles,
  labels,
  granges,
  per_locus_stat,
  aggregate_by,
  remove_top,
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

  GRanges object to summarize. Should have a valid name field.

- per_locus_stat:

  Aggregating function (per locus). Mean by default. Choices: min, max,
  sd, mean.

- aggregate_by:

  Statistic to aggregate per group. If NULL, values are not aggregated.
  This is the behavior by default.

- remove_top:

  Return range 0-(1-remove_top). By default returns the whole
  distribution (remove_top == 0).

- default_na:

  Default value for missing values

- scaling:

  If none, no operation is performed (default). If relative, values are
  divided by global mean (1x genome coverage).

## Value

An aggregated dataframe
