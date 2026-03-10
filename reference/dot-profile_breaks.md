# Compute where the break ticks go in heatmap and profile

Compute where the break ticks go in heatmap and profile

## Usage

``` r
.profile_breaks(nrows, upstream, downstream, bin_size, mode)
```

## Arguments

- nrows:

  Number of bins

- upstream:

  Number of base pairs to include upstream of loci.

- downstream:

  Number of base pairs to include downstream of loci.

- bin_size:

  Bin size. Length of bin in base pairs. The lower, the higher the
  resolution.

- mode:

  How to handle differences in lengths across loci:

  stretch: Anchor each locus on both sides.

  start: Anchor all loci on start.

  end: Anchor all loci on end.

  center: Center all loci.

## Value

Array of numeric
