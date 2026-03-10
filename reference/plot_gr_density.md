# Density plot of a precalculated GRanges object

Plots a density plot from two mcols of a GRanges object files that were
already calculated using bw_loci or bw_bins

## Usage

``` r
plot_gr_density(
  gr,
  x,
  y,
  plot_binwidth = 0.05,
  remove_top = 0,
  verbose = TRUE,
  selection = NULL
)
```

## Arguments

- gr:

  Scored GRanges object

- x:

  Column in gr corresponding to the x axis

- y:

  Column in gr corresponding to the y axis.

- plot_binwidth:

  Resolution of the bins in the density histogram (different to genomic
  bin size)

- remove_top:

  Return range 0-(1-remove_top). By default returns the whole
  distribution (remove_top == 0).

- verbose:

  Put a caption with relevant parameters on the plot.

- selection:

  A GRanges object to restrict binning to a certain set of intervals.
  This is useful for debugging and improving performance of locus
  specific analyses.

## Value

A ggplot object.
