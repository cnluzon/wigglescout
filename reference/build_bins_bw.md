# Build a tiles GRanges out of a bigWig file on a given bin size.

It does not need the reference genome as that information can be parsed
from the bigwig file

## Usage

``` r
build_bins_bw(bwfile, bin_size)
```

## Arguments

- bwfile:

  File to tile

- bin_size:

  Size of the tile in basepairs

## Value

A GRanges object with bins of bin_size bins matching the genome the
bigWig file was mapped to.

## Examples

``` r

bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
build_bins_bw(bw, 50000)
#> GRanges object with 2070 ranges and 0 metadata columns:
#>          seqnames              ranges strand
#>             <Rle>           <IRanges>  <Rle>
#>      [1]    chr15             1-50000      *
#>      [2]    chr15        50001-100000      *
#>      [3]    chr15       100001-150000      *
#>      [4]    chr15       150001-200000      *
#>      [5]    chr15       200001-250000      *
#>      ...      ...                 ...    ...
#>   [2066]    chr15 103250001-103300000      *
#>   [2067]    chr15 103300001-103350000      *
#>   [2068]    chr15 103350001-103400000      *
#>   [2069]    chr15 103400001-103450000      *
#>   [2070]    chr15 103450001-103494974      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
