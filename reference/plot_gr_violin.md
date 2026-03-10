# Violin plot of a precalculated GRanges object

Plots a violin plot of bin distribution of a set of bigWig files
optionally overlaid with annotated bins. Bins overlapping loci of the
provided BED file will be shown as a jitter plot on top of the violin
plot.

## Usage

``` r
plot_gr_violin(
  gr,
  columns,
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

- columns:

  Columns in gr to plot

- highlight:

  BED file to use as highlight for subgroups.

- minoverlap:

  Minimum overlap required for a bin to be highlighted.

- highlight_label:

  Label for the highlighted loci set

- highlight_colors:

  Array of color values for the highlighted groups.

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
