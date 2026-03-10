# Scatterplot of values in a GRanges object with mcols x and y

Plots a scatter plot from a GRanges object that already has x and y
columns.

## Usage

``` r
.scatterplot_body(
  gr,
  highlight = NULL,
  minoverlap = 0L,
  highlight_label = NULL,
  highlight_colors = NULL
)
```

## Arguments

- gr:

  GRanges

- highlight:

  List of GRanges to use as highlight for subgroups.

- minoverlap:

  Minimum overlap required for a locus to be highlighted

- highlight_label:

  Labels for the highlight groups. If not provided, column names are
  used.

- highlight_colors:

  Array of color values for the highlighting groups

## Value

A named list where plot is a ggplot object and calculated is a list of
calculated values (for verbose mode).
