# Filter out unplaced contigs and uncharacterized chromosomes from a GRanges object.

This utility will keep only the sequences that belong to canonical
chromosomes returned by GenomeInfoDb::standardChromosomes function.

## Usage

``` r
keep_canonical(gr)
```

## Arguments

- gr:

  A filtered GRanges object

## Value

A filtered GRanges object

## Examples

``` r
bins_gr <- build_bins(bin_size = 50000, genome = "hg38")
keep_canonical(bins_gr)
#> GRanges object with 61776 ranges and 0 metadata columns:
#>           seqnames            ranges strand
#>              <Rle>         <IRanges>  <Rle>
#>       [1]     chr1           1-50000      *
#>       [2]     chr1      50001-100000      *
#>       [3]     chr1     100001-150000      *
#>       [4]     chr1     150001-200000      *
#>       [5]     chr1     200001-250000      *
#>       ...      ...               ...    ...
#>   [61772]     chrY 57050001-57100000      *
#>   [61773]     chrY 57100001-57150000      *
#>   [61774]     chrY 57150001-57200000      *
#>   [61775]     chrY 57200001-57227415      *
#>   [61776]     chrM           1-16569      *
#>   -------
#>   seqinfo: 25 sequences from an unspecified genome
```
