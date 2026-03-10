# Intersect a BigWig file over loci on a GRanges object

Intersect a BigWig file over loci on a GRanges object, taking npoints
points per locus.

## Usage

``` r
.intersect_bw_and_granges(
  bw,
  granges,
  npoints,
  ignore_strand = FALSE,
  default_na = NA_real_,
  scaling = "none"
)
```

## Arguments

- bw:

  BigWigFile object (rtracklayer).

- granges:

  GRanges object.

- npoints:

  How many points to take (different to bin size!).

- ignore_strand:

  Ignore strand information in granges.

- default_na:

  Value to use as default for missing values

- scaling:

  none if no scaling is done, relative if scaled to global mean of 1

## Value

A value matrix of dimensions len(granges) x npoints.
