# Summarizing bigWig files

## Calculating functions

This vignette uses the bundled data to showcase these functions. Please
refer to the quick start guide for more information about the bundled
data.

bigWig data can be summarized genome-wide using calculation functions
that leverage `rtracklayer` for handling bigWig files:
`bw_global_coverage`, `bw_chr_coverage`, `bw_bins`, `bw_loci`,
`bw_heatmap` or `bw_profile`.

### Genome-wide analysis

You can get a single global coverage value from a bigWig file with
`bw_global_coverage`, or compute it per chromosome or contig by calling
`bw_chr_coverage`. The latter returns a dataframe with matched
`seqnames`, `width` and a mean coverage per element.

``` r

# This returns a single value
bw_global_coverage(h33_chip)
#> [1] 2.030646

# This returns a data frame with one entry per contig
bw_chr_coverage(h33_chip)
#>   seqnames start       end     width strand    score
#> 1    chr15     1 103494974 103494974      * 2.030646
```

`bw_bins` function returns a `GRanges` object where each entry is a
locus of length `bin_size` and its score is the mean coverage.

**Note:** Small bin sizes (\< 5000bp) take longer to run.

``` r

# Several bwfiles can be provided at once
bw_bins(c(h33_chip, input_chip),
        bin_size = 2000,
        genome = "mm9",
        selection = locus)
#> GRanges object with 251 ranges and 2 metadata columns:
#>         seqnames              ranges strand | sample_H33_ChIP sample_Input
#>            <Rle>           <IRanges>  <Rle> |       <numeric>    <numeric>
#>     [1]    chr15 102598001-102600000      * |        0.548169     0.999076
#>     [2]    chr15 102600001-102602000      * |        2.382783     0.987250
#>     [3]    chr15 102602001-102604000      * |        1.872876     0.793113
#>     [4]    chr15 102604001-102606000      * |        0.640021     0.720190
#>     [5]    chr15 102606001-102608000      * |        1.035948     1.134485
#>     ...      ...                 ...    ... .             ...          ...
#>   [247]    chr15 103090001-103092000      * |         1.33303     0.961300
#>   [248]    chr15 103092001-103094000      * |         1.04270     0.977078
#>   [249]    chr15 103094001-103096000      * |         1.47173     1.212554
#>   [250]    chr15 103096001-103098000      * |         1.44952     1.117081
#>   [251]    chr15 103098001-103100000      * |         1.34690     0.749071
#>   -------
#>   seqinfo: 35 sequences from an unspecified genome; no seqlengths
```

It is also possible to provide a list of bigWig files to be used as
background to normalize the bin values to. For instance, in the previous
case, one could want to use the input values to normalize the H3.3 ChIP
data:

``` r

bw_bins(h33_chip, bg_bwfiles = input_chip,
        bin_size = 2000,
        genome = "mm9",
        selection = locus)
#> GRanges object with 251 ranges and 1 metadata column:
#>         seqnames              ranges strand | sample_H33_ChIP
#>            <Rle>           <IRanges>  <Rle> |       <numeric>
#>     [1]    chr15 102598001-102600000      * |        0.548676
#>     [2]    chr15 102600001-102602000      * |        2.413555
#>     [3]    chr15 102602001-102604000      * |        2.361423
#>     [4]    chr15 102604001-102606000      * |        0.888684
#>     [5]    chr15 102606001-102608000      * |        0.913143
#>     ...      ...                 ...    ... .             ...
#>   [247]    chr15 103090001-103092000      * |         1.38670
#>   [248]    chr15 103092001-103094000      * |         1.06716
#>   [249]    chr15 103094001-103096000      * |         1.21375
#>   [250]    chr15 103096001-103098000      * |         1.29760
#>   [251]    chr15 103098001-103100000      * |         1.79810
#>   -------
#>   seqinfo: 35 sequences from an unspecified genome; no seqlengths
```

If you want to get the log2 ratio sample / input, select
`norm_mode = "log2fc"`:

``` r

bw_bins(h33_chip, bg_bwfiles = input_chip,
        bin_size = 2000,
        genome = "mm9",
        selection = locus,
        norm_mode = "log2fc")
#> GRanges object with 251 ranges and 1 metadata column:
#>         seqnames              ranges strand | sample_H33_ChIP
#>            <Rle>           <IRanges>  <Rle> |       <numeric>
#>     [1]    chr15 102598001-102600000      * |       -0.865973
#>     [2]    chr15 102600001-102602000      * |        1.271160
#>     [3]    chr15 102602001-102604000      * |        1.239656
#>     [4]    chr15 102604001-102606000      * |       -0.170258
#>     [5]    chr15 102606001-102608000      * |       -0.131087
#>     ...      ...                 ...    ... .             ...
#>   [247]    chr15 103090001-103092000      * |       0.4716506
#>   [248]    chr15 103092001-103094000      * |       0.0937752
#>   [249]    chr15 103094001-103096000      * |       0.2794686
#>   [250]    chr15 103096001-103098000      * |       0.3758455
#>   [251]    chr15 103098001-103100000      * |       0.8464701
#>   -------
#>   seqinfo: 35 sequences from an unspecified genome; no seqlengths
```

`norm_mode` is set to `"fc"` by default, which means it shows fold
change sample / input, if input (background) is provided.

#### Supported genomes

`bw_bins` function is based on `tileGenome` function from
`GenomicRanges`, which requires information about chromosome length.
This is obtained relying on function `Seqinfo` from `GenomeInfoDb`
package. So, in theory, any genome ID that would return a `Seqinfo`
object by doing `Seqinfo(genome="your_genome")` should work.

For more custom applications, there are also auxiliary functions
`build_bins` to obtain a `GRanges` object that can then be used in
`bw_loci`, and also a further `build_bins_bw` which takes a bigWig file
and tiles it based on the references found in it. This is useful when
you are not certain exactly which reference version was used in a bigWig
or when you do not want to rely on internet connection for tiling.

### Locus-based analysis

These functions are useful when you have some genomic annotations you
want to look at, rather than general genome-wide trends.

#### Individual loci

To calculate values on a *locus* specific manner you can use `bw_loci`
function. For instance, we can look at the genes at this region:

``` r

bw_loci(c(h33_chip, input_chip), loci = genes)
#> GRanges object with 26 ranges and 3 metadata columns:
#>        seqnames              ranges strand | sample_H33_ChIP sample_Input
#>           <Rle>           <IRanges>  <Rle> |       <numeric>    <numeric>
#>    [1]    chr15 102751562-102759245      + |         2.91694     0.989028
#>    [2]    chr15 102767284-102768948      + |         3.61492     1.082491
#>    [3]    chr15 102774493-102778377      - |         3.97043     1.135737
#>    [4]    chr15 102784957-102787132      + |         2.17507     1.122432
#>    [5]    chr15 102797227-102802328      + |         1.45977     1.072421
#>    ...      ...                 ...    ... .             ...          ...
#>   [22]    chr15 103078643-103082118      - |         5.51339      1.22775
#>   [23]    chr15 103078643-103083241      - |         5.49223      1.25069
#>   [24]    chr15 103078643-103083583      - |         5.45826      1.25206
#>   [25]    chr15 103078643-103085870      - |         4.87355      1.23020
#>   [26]    chr15 103078643-103088834      - |         4.42227      1.20272
#>               name
#>        <character>
#>    [1]      Hoxc13
#>    [2]      Hoxc12
#>    [3]      Hotair
#>    [4]      Hoxc11
#>    [5]      Hoxc10
#>    ...         ...
#>   [22]        Nfe2
#>   [23]        Nfe2
#>   [24]        Nfe2
#>   [25]        Nfe2
#>   [26]        Nfe2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

Many of the options available for `bw_bins` are shared with `bw_loci`,
so you can also use background bigWig files:

``` r

bw_loci(h33_chip, bg_bwfiles = input_chip, loci = genes)
#> GRanges object with 26 ranges and 2 metadata columns:
#>        seqnames              ranges strand | sample_H33_ChIP        name
#>           <Rle>           <IRanges>  <Rle> |       <numeric> <character>
#>    [1]    chr15 102751562-102759245      + |         2.94930      Hoxc13
#>    [2]    chr15 102767284-102768948      + |         3.33944      Hoxc12
#>    [3]    chr15 102774493-102778377      - |         3.49591      Hotair
#>    [4]    chr15 102784957-102787132      + |         1.93782      Hoxc11
#>    [5]    chr15 102797227-102802328      + |         1.36119      Hoxc10
#>    ...      ...                 ...    ... .             ...         ...
#>   [22]    chr15 103078643-103082118      - |         4.49064        Nfe2
#>   [23]    chr15 103078643-103083241      - |         4.39137        Nfe2
#>   [24]    chr15 103078643-103083583      - |         4.35942        Nfe2
#>   [25]    chr15 103078643-103085870      - |         3.96159        Nfe2
#>   [26]    chr15 103078643-103088834      - |         3.67690        Nfe2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

A name column is provided if the given `BED` file includes *loci* IDs.

You can also provide `GRanges` objects instead of `BED` files:

``` r

genes_granges <- rtracklayer::import(genes, format = "BED")

bw_loci(h33_chip, bg_bwfiles = input_chip, loci = genes_granges)
#> GRanges object with 26 ranges and 2 metadata columns:
#>        seqnames              ranges strand | sample_H33_ChIP        name
#>           <Rle>           <IRanges>  <Rle> |       <numeric> <character>
#>    [1]    chr15 102751562-102759245      + |         2.94930      Hoxc13
#>    [2]    chr15 102767284-102768948      + |         3.33944      Hoxc12
#>    [3]    chr15 102774493-102778377      - |         3.49591      Hotair
#>    [4]    chr15 102784957-102787132      + |         1.93782      Hoxc11
#>    [5]    chr15 102797227-102802328      + |         1.36119      Hoxc10
#>    ...      ...                 ...    ... .             ...         ...
#>   [22]    chr15 103078643-103082118      - |         4.49064        Nfe2
#>   [23]    chr15 103078643-103083241      - |         4.39137        Nfe2
#>   [24]    chr15 103078643-103083583      - |         4.35942        Nfe2
#>   [25]    chr15 103078643-103085870      - |         3.96159        Nfe2
#>   [26]    chr15 103078643-103088834      - |         3.67690        Nfe2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

#### Aggregated loci

It is also possible to aggregate values by some kind of category. This
is useful for `BED` files where ID is not an individual entity but
rather a category, for instance, ChromHMM annotations:

``` r

chrom_ranges <- import(chromhmm)

chrom_ranges
#> GRanges object with 164 ranges and 2 metadata columns:
#>         seqnames              ranges strand |               name     score
#>            <Rle>           <IRanges>  <Rle> |        <character> <numeric>
#>     [1]    chr15 102602601-102606600      * |     12_Heterochrom         0
#>     [2]    chr15 102606601-102607600      * |       11_Repressed         0
#>     [3]    chr15 102615401-102615800      * |       15_Insulator         0
#>     [4]    chr15 102615801-102616800      * |       11_Repressed         0
#>     [5]    chr15 102616801-102626600      * |     12_Heterochrom         0
#>     ...      ...                 ...    ... .                ...       ...
#>   [160]    chr15 103086201-103088200      * |         2_Weak_Txn         0
#>   [161]    chr15 102627601-102636000      * |     12_Heterochrom         0
#>   [162]    chr15 102891201-102891600      * | 10_Poised_Promoter         0
#>   [163]    chr15 103017601-103040000      * |   1_Txn_Elongation         0
#>   [164]    chr15 102758201-102761800      * |       11_Repressed         0
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

It is possible to get a single value per category using `aggregate_by`
parameter.

``` r

# If aggregate_by is not provided, you get a non-aggregated GRanges
bw_loci(h33_chip,
        bg_bwfiles = input_chip,
        loci = chrom_ranges,
        aggregate_by = "mean") 
#>                    sample_H33_ChIP
#> 1_Txn_Elongation          4.089008
#> 2_Weak_Txn                2.438608
#> 4_Poised_Enhancer         3.298239
#> 5_Active_Promoter         3.807038
#> 6_Strong_Enhancer         4.691568
#> 7_Active_Promoter         5.561953
#> 8_Strong_Enhancer         6.381945
#> 9_Strong_Enhancer         3.320948
#> 9_Txn_Transition          3.458106
#> 10_Poised_Promoter        2.863107
#> 11_Repressed              1.961251
#> 12_Heterochrom            1.038954
#> 14_Heterochrom            3.597231
#> 15_Insulator              2.400263
```

**Note:** In this case you will not get a `GRanges` object but a
`data.frame`.

For these calculations you **always** get some kind of **average** per
locus. However, different ways of treating the average values are
available by changing the `aggregate_by` parameter:

- `true_mean`: This calculates a mean adding up all the values included
  per *locus*. It is most meaningful when loci represent categories
  rather than biological entities.
- `mean`: This calculates a mean-of-means instead, where each *locus*
  mean value is calculated and the mean of this distribution of means is
  returned as a result.
- `median`: This calculates a median-of-means. It is useful to perform
  an analysis that is more robust to outliers, but it does not work well
  with sparse data.

#### Profile averages

You may be interested in the **shape** the signal takes across a set of
*loci*. This is calculated by `bw_profile`. This function takes each
*locus* in the provided bed file and slices it into a fixed number of
bins (this is determined by `bin_size` parameter). Each bin is
aggregated together across loci.

``` r

profile <- bw_profile(h33_chip, loci = genes, bin_size = 100)

head(profile)
#>       mean   sderror   median index          sample
#> 1 2.651766 0.4211416 1.868056     1 sample_H33_ChIP
#> 2 2.353699 0.3952149 1.744718     2 sample_H33_ChIP
#> 3 2.372053 0.4093088 1.675340     3 sample_H33_ChIP
#> 4 2.516681 0.2644566 2.112161     4 sample_H33_ChIP
#> 5 2.288359 0.2839340 1.683047     5 sample_H33_ChIP
#> 6 2.415918 0.3109912 1.837221     6 sample_H33_ChIP
```

**Note:**. `bw_profile` returns data as long format, since it returns
more than one value per data point. You get: mean, standard error,
median value per point. Each point is represented by index. This index
determines its bin.

You can also provide `upstream` and `downstream` base pairs values to
include in this calculation.

``` r

profile <- bw_profile(h33_chip, loci = genes, bin_size = 100,
                      upstream = 500, downstream = 500)

head(profile)
#>       mean   sderror   median index          sample
#> 1 3.428869 0.5597654 2.405090     1 sample_H33_ChIP
#> 2 3.949016 0.6592606 2.949832     2 sample_H33_ChIP
#> 3 3.247899 0.4771775 2.590096     3 sample_H33_ChIP
#> 4 3.229546 0.4894977 2.058201     4 sample_H33_ChIP
#> 5 3.079227 0.4333806 2.497594     5 sample_H33_ChIP
#> 6 3.639487 0.3295900 3.067325     6 sample_H33_ChIP
```

Additionally, you can use `mode` parameter to define how these locus are
aligned before binning. Available modes are: `start`, `end`, `center`
(align the loci around these points) or `stretch` which anchors features
`start` to `end` and stretches them when length is different.

#### Per-locus profiles

This shape created by `bw_profile` function is still an average. So you
may be interested on a per-locus profile value. This is what you get
with `bw_heatmap`. `bw_heatmap` returns a `matrix` where each row is a
locus and each column is a bin. The rest of parameters are analogous to
`bw_profile`:

``` r

mat <- bw_heatmap(h33_chip, loci = genes, bin_size = 1000,
                  upstream = 5000, downstream = 5000)

# Note that bw_heatmap returns a list of matrices, so it supports arrays of 
# bwfiles as input as well.
head(mat[[1]])
#>          [,1]      [,2]      [,3]      [,4]     [,5]     [,6]     [,7]
#> [1,] 1.284604 1.3293102 1.7112652 2.8632013 3.007691 2.800395 1.766877
#> [2,] 4.952701 6.5925578 6.4507772 5.3514131 3.370813 2.863804 3.339249
#> [3,] 1.817719 1.8833939 2.0797163 1.5973820 3.622250 5.756086 3.083041
#> [4,] 1.666320 2.4688960 1.5384674 1.5898725 1.831000 2.566131 1.285138
#> [5,] 1.640442 2.2717496 1.9870765 1.9248061 2.705879 2.634167 1.318919
#> [6,] 1.259195 0.9254259 0.8669583 0.9295532 2.033445 2.275877 1.419363
#>           [,8]      [,9]    [,10]    [,11]    [,12]    [,13]    [,14]
#> [1,] 3.3833293 4.3805848 1.667400 1.561681 3.661995 4.761350 6.651459
#> [2,] 3.0416741 5.1382823 5.887521 9.196965 4.561860 3.057636 3.688586
#> [3,] 2.7693828 4.2945841 4.607247 3.440608 3.096550 7.539874 8.274333
#> [4,] 2.6849959 2.1768526 2.917202 1.981451 1.685489 5.213727 2.688018
#> [5,] 0.9578958 0.7213599 1.690166 2.632981 2.310449 1.297612 2.311900
#> [6,] 2.3248209 3.5973599 3.592973 1.099763 1.261885 4.124952 3.420550
```

### General shared options

Each `wigglescout` function has a different set of parameters. You can
check the documentation with `?` operator:
[`?bw_bins`](https://cnluzon.github.io/wigglescout/reference/bw_bins.md).

Here are a few that are commonly appearing in many of the functions:

- `labels`: This parameter gives name to the value column of the
  `data.frame` or `GRanges` object obtained. By default the names
  correspond to the file names.
- `remove_top`: Exclude a certain top quantile of the returning object.
  This value must be in \[0, 1\] and by default is 0 (no values
  removed). Underneath, threshold value is calculated by `quantile`
  function and used as a threshold, where the exact value is included
  (values \<= threshold are returned). If the calculated data structure
  has more than one column (i.e. multiple files are provided), quantiles
  are calculated over `rowMeans`. This prioritizes cases where outliers
  are such in more than one sample. In the case where `NA` values are
  found, these are also dropped if `remove_top > 0`.
- `per_locus_stat`: Aggregating function per locus. This is an interface
  with the underlying `rtracklayer` summary function. It defaults to
  `mean` and it is not recommended to change it.

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] wigglescout_0.21.2   rtracklayer_1.72.0   GenomicRanges_1.64.0
#> [4] Seqinfo_1.2.0        IRanges_2.46.0       S4Vectors_0.50.1    
#> [7] BiocGenerics_0.58.1  generics_0.1.4       ggplot2_4.0.3       
#> 
#> loaded via a namespace (and not attached):
#>  [1] SummarizedExperiment_1.42.0 beeswarm_0.4.0             
#>  [3] gtable_0.3.6                rjson_0.2.23               
#>  [5] xfun_0.58                   bslib_0.11.0               
#>  [7] Biobase_2.72.0              lattice_0.22-9             
#>  [9] vctrs_0.7.3                 tools_4.6.0                
#> [11] bitops_1.0-9                curl_7.1.0                 
#> [13] parallel_4.6.0              tibble_3.3.1               
#> [15] pkgconfig_2.0.3             Matrix_1.7-5               
#> [17] RColorBrewer_1.1-3          cigarillo_1.2.0            
#> [19] S7_0.2.2                    desc_1.4.3                 
#> [21] lifecycle_1.0.5             stringr_1.6.0              
#> [23] compiler_4.6.0              farver_2.1.2               
#> [25] Rsamtools_2.28.0            textshaping_1.0.5          
#> [27] Biostrings_2.80.1           codetools_0.2-20           
#> [29] vipor_0.4.7                 GenomeInfoDb_1.48.0        
#> [31] htmltools_0.5.9             sass_0.4.10                
#> [33] RCurl_1.98-1.19             yaml_2.3.12                
#> [35] tidyr_1.3.2                 pillar_1.11.1              
#> [37] pkgdown_2.2.0               crayon_1.5.3               
#> [39] jquerylib_0.1.4             BiocParallel_1.46.0        
#> [41] cachem_1.1.0                DelayedArray_0.38.2        
#> [43] abind_1.4-8                 tidyselect_1.2.1           
#> [45] digest_0.6.39               stringi_1.8.7              
#> [47] purrr_1.2.2                 dplyr_1.2.1                
#> [49] restfulr_0.0.16             fastmap_1.2.0              
#> [51] grid_4.6.0                  SparseArray_1.12.2         
#> [53] cli_3.6.6                   magrittr_2.0.5             
#> [55] S4Arrays_1.12.0             XML_3.99-0.23              
#> [57] withr_3.0.2                 UCSC.utils_1.8.0           
#> [59] scales_1.4.0                ggbeeswarm_0.7.3           
#> [61] rmarkdown_2.31              XVector_0.52.0             
#> [63] httr_1.4.8                  matrixStats_1.5.0          
#> [65] ragg_1.5.2                  evaluate_1.0.5             
#> [67] ggrastr_1.0.2               knitr_1.51                 
#> [69] BiocIO_1.22.0               rlang_1.2.0                
#> [71] glue_1.8.1                  jsonlite_2.0.0             
#> [73] R6_2.6.1                    MatrixGenerics_1.24.0      
#> [75] GenomicAlignments_1.48.0    systemfonts_1.3.2          
#> [77] fs_2.1.0
```
