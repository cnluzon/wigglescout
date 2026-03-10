# Calculate profile values of a bigWig file over GRanges object.

Calculates profile values of a bigWig file over the loci speficied in a
GRanges object. Data points are taken each bin_size base pairs.

## Usage

``` r
.calculate_bw_profile(
  bw,
  granges,
  bg_bw = NULL,
  label = NULL,
  mode = "stretch",
  bin_size = 100,
  upstream = 2500,
  downstream = 2500,
  middle = NULL,
  ignore_strand = FALSE,
  norm_func = identity,
  remove_top = 0,
  default_na = NA_real_,
  scaling = "none"
)
```

## Arguments

- bw:

  BigWig file to be summarized.

- granges:

  GRanges object

- bg_bw:

  BigWig file to be used as background.

- label:

  Name to give to the values

- mode:

  How to handle differences in lengths across loci:

  stretch: Anchor each locus on both sides.

  start: Anchor all loci on start.

  end: Anchor all loci on end.

  center: Center all loci.

- bin_size:

  Bin size. Length of bin in base pairs. The lower, the higher the
  resolution.

- upstream:

  Number of base pairs to include upstream of loci.

- downstream:

  Number of base pairs to include downstream of loci.

- middle:

  Number of base pairs that the middle section has (in stretch mode). If
  not provided, median length of all loci is used.

- ignore_strand:

  Whether to use strand information in BED file.

- norm_func:

  Normalization function

- remove_top:

  Return range 0-(1-remove_top). By default returns the whole
  distribution (remove_top == 0).

- default_na:

  Default value for missing values

- scaling:

  Whether to use the bigWig values as they are (none - default) or
  calculate relative enrichment (relative) by dividing values by global
  average.

## Value

A DataFrame with the aggregated scores

## Details

Loci are aligned depending on mode parameter:

- stretch. Aligns all starts and all ends, sort of stretching the loci.
  The median of these lenghts is taken as the pseudo-length in order to
  show a realistic plot when displayed.

- start. All loci are aligned by start.

- end. All loci are aligned by end.

- center. All loci are aligned by center.

Adapted from seqplots rtracklayer wrapper functions to compute coverage:
For more on seqplots:
https://www.bioconductor.org/packages/release/bioc/html/seqplots.html
