# Calculate matrix for stretch mode

Calculate matrix for stretch mode

## Usage

``` r
.calculate_stretch_matrix(
  bw,
  granges,
  bin_size = 100,
  upstream = 2500,
  downstream = 2500,
  middle = NULL,
  ignore_strand = FALSE,
  default_na = NA_real_,
  scaling = "none"
)
```

## Arguments

- bw:

  BigWigFile object

- granges:

  GRanges object

- bin_size:

  Bin size

- upstream:

  Number of basepairs upstream

- downstream:

  Number of basepairs downstream

- middle:

  Number of base pairs that the middle section has. If not provided,
  median is used.

- ignore_strand:

  Ignore strand (bool)

- default_na:

  Default value for missing values

- scaling:

  none if no scaling is performed (default), relative if 1x scaling is
  done (divide by global coverage).

## Value

Summary matrix
