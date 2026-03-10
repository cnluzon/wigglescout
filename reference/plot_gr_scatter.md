# Scatterplot of a precalculated GRanges object

Plots a scatter plot from two mcols of a GRanges object files and an
optional set of BED files as highlighted annotations. Bins are
highlighted if there is at least minoverlap base pairs overlap with any
loci in BED file.

## Usage

``` r
plot_gr_scatter(
  gr,
  x,
  y,
  highlight = NULL,
  minoverlap = 0L,
  highlight_label = NULL,
  highlight_colors = NULL,
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

- highlight:

  List of bed files to use as highlight for subgroups.

- minoverlap:

  Minimum overlap required for a bin to be highlighted

- highlight_label:

  Labels for the highlight groups. If not provided, filenames are used.

- highlight_colors:

  Array of color values for the highlighting groups

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

## Details

If specifying minoverlap, you must take into account the bin_size
parameter and the size of the loci you are providing as BED file.

This function does not calculate background normalization or anything,
you can do that in the prior call to bw_bins or bw_loci
