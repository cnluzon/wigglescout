# Helper function for calculating guide lines and labels in heatmap

Helper function for calculating guide lines and labels in heatmap

## Usage

``` r
.heatmap_lines(
  loci,
  nbins,
  bin_size,
  upstream,
  downstream,
  mode,
  expand = TRUE
)
```

## Arguments

- loci:

  Number of loci (rows)

- nbins:

  Number of bins (columns)

- bin_size:

  Bin size. Length of bin in base pairs. The lower, the higher the
  resolution.

- upstream:

  Number of base pairs to include upstream of loci.

- downstream:

  Number of base pairs to include downstream of loci.

- mode:

  How to handle differences in lengths across loci:

  stretch: Anchor each locus on both sides.

  start: Anchor all loci on start.

  end: Anchor all loci on end.

  center: Center all loci.

- expand:

  Whether to apply padding to the ticks. Otherwise a strange box shows
  on heatmap.

## Value

A list of ggproto objects to be plotted.
