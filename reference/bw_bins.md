# Build a binned-scored GRanges object from a set of bigWig files

Build a binned-scored GRanges object from one or many bigWig files. The
aggregating function per bin can be min, max, sd, mean.

## Usage

``` r
bw_bins(
  bwfiles,
  bg_bwfiles = NULL,
  labels = NULL,
  per_locus_stat = "mean",
  bin_size = 10000,
  genome = "mm9",
  selection = NULL,
  norm_mode = "fc",
  remove_top = 0,
  default_na = NA_real_,
  scaling = "none",
  canonical = FALSE
)
```

## Arguments

- bwfiles:

  Path or array of paths to the bigWig files to be summarized.

- bg_bwfiles:

  Path or array of paths to the bigWig files to be used as background.

- labels:

  List of names to give to the mcols of the returned GRanges object. If
  NULL, file names are used.

- per_locus_stat:

  Aggregating function (per locus). Mean by default. Choices: min, max,
  sd, mean.

- bin_size:

  Bin size.

- genome:

  Genome. Available choices are mm9, hg38.

- selection:

  A GRanges object to restrict binning to a certain set of intervals.
  This is useful for debugging and improving performance of locus
  specific analyses.

- norm_mode:

  Function to apply to normalize bin values. Default fc: divides bw /
  bg. Alternative: log2fc: returns log2(bw/bg).

- remove_top:

  Return range 0-(1-remove_top). By default returns the whole
  distribution (remove_top == 0).

- default_na:

  Default value for missing values

- scaling:

  If none, no operation is performed (default). If relative, values are
  divided by global mean (1x genome coverage).

- canonical:

  Use only canonical chromosomes (default: FALSE)

## Value

A GRanges object with each bwfile as a metadata column named after
labels, if provided, or after filenames otherwise.

## Details

bwfiles and bg_bwfiles must have the same length. If you are using same
background for several files, then file paths must be repeated
accordingly.

norm_mode can be either "fc", where values will represent bigwig /
bg_bigwig, or "log2fc", values will represent log2(bigwig / bg_bigwig)
per bin.

## Examples

``` r
# Get the raw files
bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
bw2 <- system.file("extdata",
                   "sample_H3K9me3_ChIP.bw", package="wigglescout")
bw_inp <- system.file("extdata", "sample_Input.bw", package="wigglescout")

# Sample bigWig files only have valid values on this region
locus <- GenomicRanges::GRanges(
  seqnames = "chr15",
  IRanges::IRanges(102600000, 103100000)
)

# Run single bw. The larger the bin size, the faster this is.
bw_bins(bw, bin_size = 100000, labels = c("H33"), selection = locus)
#> Warning: Not all bigWig files or the genome tiling provided share the same sequence info. 
#>    Common to all (1): chr15  
#>    Missing in some of the files: chr1, chr2, chr3, chr4, chr5 ... 
#>    This can be due to different versions of the same reference genome, or to completely different organisms. Make sure these match!
#> GRanges object with 6 ranges and 1 metadata column:
#>       seqnames              ranges strand |       H33
#>          <Rle>           <IRanges>  <Rle> | <numeric>
#>   [1]    chr15 102500001-102600000      * |  1.343961
#>   [2]    chr15 102600001-102700000      * |  0.934239
#>   [3]    chr15 102700001-102800000      * |  1.956845
#>   [4]    chr15 102800001-102900000      * |  2.130207
#>   [5]    chr15 102900001-103000000      * |  2.010765
#>   [6]    chr15 103000001-103100000      * |  3.355660
#>   -------
#>   seqinfo: 35 sequences from an unspecified genome; no seqlengths

# Use of some parameters
bw_bins(bw, bg_bwfiles = bw_inp, bin_size = 100000,
       labels = c("H33"), norm_mode = "log2fc", selection = locus)
#> Warning: Not all bigWig files or the genome tiling provided share the same sequence info. 
#>    Common to all (1): chr15  
#>    Missing in some of the files: chr1, chr2, chr3, chr4, chr5 ... 
#>    This can be due to different versions of the same reference genome, or to completely different organisms. Make sure these match!
#> GRanges object with 6 ranges and 1 metadata column:
#>       seqnames              ranges strand |       H33
#>          <Rle>           <IRanges>  <Rle> | <numeric>
#>   [1]    chr15 102500001-102600000      * |  0.384133
#>   [2]    chr15 102600001-102700000      * | -0.160072
#>   [3]    chr15 102700001-102800000      * |  0.912455
#>   [4]    chr15 102800001-102900000      * |  1.014908
#>   [5]    chr15 102900001-103000000      * |  0.781969
#>   [6]    chr15 103000001-103100000      * |  1.590060
#>   -------
#>   seqinfo: 35 sequences from an unspecified genome; no seqlengths

# Multiple bigWig
bw_bins(c(bw, bw2),
        bin_size = 100000,
        bg_bwfiles = c(bw_inp, bw_inp),
        norm_mode = "log2fc",
        selection = locus)
#> Warning: Not all bigWig files or the genome tiling provided share the same sequence info. 
#>    Common to all (1): chr15  
#>    Missing in some of the files: chr1, chr2, chr3, chr4, chr5 ... 
#>    This can be due to different versions of the same reference genome, or to completely different organisms. Make sure these match!
#> GRanges object with 6 ranges and 2 metadata columns:
#>       seqnames              ranges strand | sample_H33_ChIP sample_H3K9me3_ChIP
#>          <Rle>           <IRanges>  <Rle> |       <numeric>           <numeric>
#>   [1]    chr15 102500001-102600000      * |        0.384133            -1.50269
#>   [2]    chr15 102600001-102700000      * |       -0.160072            -1.32163
#>   [3]    chr15 102700001-102800000      * |        0.912455            -2.07487
#>   [4]    chr15 102800001-102900000      * |        1.014908            -2.66062
#>   [5]    chr15 102900001-103000000      * |        0.781969            -2.39259
#>   [6]    chr15 103000001-103100000      * |        1.590060            -3.99871
#>   -------
#>   seqinfo: 35 sequences from an unspecified genome; no seqlengths
```
