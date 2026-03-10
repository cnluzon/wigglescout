# Build a tiles GRanges out of a bigWig file on a given bin size.

It does not need the reference genome as that information can be parsed
from the bigwig file

## Usage

``` r
tile_seqinfo(seqinfo, bin_size)
```

## Arguments

- seqinfo:

  A Seqinfo object

- bin_size:

  Size of the tile in basepairs

## Value

A GRanges object with bins of bin_size bins matching the genome the
bigWig file was mapped to.

## Examples

``` r
tile_seqinfo(GenomeInfoDb::Seqinfo(genome = "mm9"), 50000)
#> GRanges object with 54531 ranges and 0 metadata columns:
#>               seqnames          ranges strand
#>                  <Rle>       <IRanges>  <Rle>
#>       [1]         chr1         1-50000      *
#>       [2]         chr1    50001-100000      *
#>       [3]         chr1   100001-150000      *
#>       [4]         chr1   150001-200000      *
#>       [5]         chr1   200001-250000      *
#>       ...          ...             ...    ...
#>   [54527] chrUn_random 5700001-5750000      *
#>   [54528] chrUn_random 5750001-5800000      *
#>   [54529] chrUn_random 5800001-5850000      *
#>   [54530] chrUn_random 5850001-5900000      *
#>   [54531] chrUn_random 5900001-5900358      *
#>   -------
#>   seqinfo: 35 sequences from an unspecified genome
```
