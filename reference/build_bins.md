# Build a unscored bins GRanges object.

Build a GRanges of bins of a given size, for a specific genome.
Supported genomes rely on
[Seqinfo](https://rdrr.io/pkg/Seqinfo/man/Seqinfo-class.html). This
requires internet access to work.

## Usage

``` r
build_bins(bin_size = 10000, genome = "mm9", canonical = FALSE)
```

## Arguments

- bin_size:

  Bin size.

- genome:

  Genome. Supported: any registered NCBI or UCSC genome see
  registered_UCSC_genomes(), registered_NCBI_assemblies() from
  GenomeInfoDb

- canonical:

  Use only canonical chromosomes (default: FALSE)

## Value

A GRanges object with a tiled genome

## Examples

``` r
build_bins(bin_size = 50000, genome = "mm9")
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
build_bins(bin_size = 50000, genome = "hg38")
#> GRanges object with 66390 ranges and 0 metadata columns:
#>                      seqnames        ranges strand
#>                         <Rle>     <IRanges>  <Rle>
#>       [1]                chr1       1-50000      *
#>       [2]                chr1  50001-100000      *
#>       [3]                chr1 100001-150000      *
#>       [4]                chr1 150001-200000      *
#>       [5]                chr1 200001-250000      *
#>       ...                 ...           ...    ...
#>   [66386] chrX_MU273397v1_alt 100001-150000      *
#>   [66387] chrX_MU273397v1_alt 150001-200000      *
#>   [66388] chrX_MU273397v1_alt 200001-250000      *
#>   [66389] chrX_MU273397v1_alt 250001-300000      *
#>   [66390] chrX_MU273397v1_alt 300001-330493      *
#>   -------
#>   seqinfo: 711 sequences from an unspecified genome
build_bins(bin_size = 50000, genome = "mm10", canonical = TRUE)
#> GRanges object with 54521 ranges and 0 metadata columns:
#>           seqnames            ranges strand
#>              <Rle>         <IRanges>  <Rle>
#>       [1]     chr1           1-50000      *
#>       [2]     chr1      50001-100000      *
#>       [3]     chr1     100001-150000      *
#>       [4]     chr1     150001-200000      *
#>       [5]     chr1     200001-250000      *
#>       ...      ...               ...    ...
#>   [54517]     chrY 91550001-91600000      *
#>   [54518]     chrY 91600001-91650000      *
#>   [54519]     chrY 91650001-91700000      *
#>   [54520]     chrY 91700001-91744698      *
#>   [54521]     chrM           1-16299      *
#>   -------
#>   seqinfo: 22 sequences from an unspecified genome
```
