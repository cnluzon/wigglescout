# Calculate profile values of a bigWig file over a BED file.

Calculates profile values of a set of tracks over the loci speficied in
a BED file. Data points are taken each bin_size base pairs.

## Usage

``` r
bw_profile(
  bwfiles,
  bg_bwfiles = NULL,
  loci = NULL,
  labels = NULL,
  mode = "stretch",
  bin_size = 100,
  upstream = 2500,
  downstream = 2500,
  middle = NULL,
  ignore_strand = FALSE,
  norm_mode = "fc",
  remove_top = 0,
  default_na = NA_real_,
  scaling = "none"
)
```

## Arguments

- bwfiles:

  Path or array of paths to the bigWig files to be summarized.

- bg_bwfiles:

  Path or array of paths to the bigWig files to be used as background.

- loci:

  BED file or GRanges to summarize

- labels:

  List of names to give to the mcols of the returned GRanges object. If
  NULL, file names are used.

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

- norm_mode:

  Function to apply to normalize bin values. Default fc: divides bw /
  bg. Alternative: log2fc: returns log2(bw/bg).

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

a data frame in long format

## Details

Loci are aligned depending on mode parameter:

- stretch. Aligns all starts and all ends, sort of stretching the loci.
  The median of these lenghts is taken as the pseudo-length in order to
  show a realistic plot when displayed.

- start. All loci are aligned by start.

- end. All loci are aligned by end.

- center. All loci are aligned by center.

## Examples

``` r
# Get the raw files
bed <- system.file("extdata", "sample_genes_mm9.bed", package="wigglescout")
bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
bw2 <- system.file("extdata",
                   "sample_H3K9me3_ChIP.bw", package="wigglescout")

# Profiles are returned in long format and include mean, median and stderror
bw_profile(bw, loci = bed, mode = "stretch")
#>        mean   sderror   median index          sample
#> 1  2.651766 0.4211416 1.868056     1 sample_H33_ChIP
#> 2  2.353699 0.3952149 1.744718     2 sample_H33_ChIP
#> 3  2.372053 0.4093088 1.675340     3 sample_H33_ChIP
#> 4  2.516681 0.2644566 2.112161     4 sample_H33_ChIP
#> 5  2.288359 0.2839340 1.683047     5 sample_H33_ChIP
#> 6  2.415918 0.3109912 1.837221     6 sample_H33_ChIP
#> 7  2.256240 0.3219986 1.428663     7 sample_H33_ChIP
#> 8  2.446569 0.2598446 2.232930     8 sample_H33_ChIP
#> 9  2.670303 0.2688390 2.307446     9 sample_H33_ChIP
#> 10 2.729953 0.3336880 2.461619    10 sample_H33_ChIP
#> 11 2.589913 0.2979482 2.258625    11 sample_H33_ChIP
#> 12 2.836956 0.3623161 2.669751    12 sample_H33_ChIP
#> 13 2.503650 0.2598396 2.448771    13 sample_H33_ChIP
#> 14 2.965801 0.3900738 2.854759    14 sample_H33_ChIP
#> 15 2.999204 0.4036154 2.507870    15 sample_H33_ChIP
#> 16 3.028204 0.4874518 2.471897    16 sample_H33_ChIP
#> 17 3.386287 0.6245647 2.235501    17 sample_H33_ChIP
#> 18 3.700873 0.7401184 2.556693    18 sample_H33_ChIP
#> 19 3.609102 0.5534036 3.032059    19 sample_H33_ChIP
#> 20 3.701056 0.5778502 3.003792    20 sample_H33_ChIP
#> 21 3.428869 0.5597654 2.405090    21 sample_H33_ChIP
#> 22 3.949016 0.6592606 2.949832    22 sample_H33_ChIP
#> 23 3.247899 0.4771775 2.590096    23 sample_H33_ChIP
#> 24 3.229546 0.4894977 2.058201    24 sample_H33_ChIP
#> 25 3.079227 0.4333806 2.497594    25 sample_H33_ChIP
#> 26 3.639487 0.3295900 3.067325    26 sample_H33_ChIP
#> 27 3.615557 0.3111358 3.394724    27 sample_H33_ChIP
#> 28 3.114177 0.2877722 2.646299    28 sample_H33_ChIP
#> 29 3.087110 0.2202857 2.971876    29 sample_H33_ChIP
#> 30 2.889450 0.2188670 2.972835    30 sample_H33_ChIP
#> 31 2.875842 0.2616709 2.387392    31 sample_H33_ChIP
#> 32 2.929630 0.3522394 2.125642    32 sample_H33_ChIP
#> 33 3.015563 0.3073536 2.458524    33 sample_H33_ChIP
#> 34 3.262732 0.3419128 2.545657    34 sample_H33_ChIP
#> 35 3.247373 0.3754521 2.880982    35 sample_H33_ChIP
#> 36 2.866761 0.3246729 2.409001    36 sample_H33_ChIP
#> 37 3.005472 0.3119875 2.753080    37 sample_H33_ChIP
#> 38 2.837359 0.3052706 2.548127    38 sample_H33_ChIP
#> 39 2.783164 0.2945724 2.463294    39 sample_H33_ChIP
#> 40 2.764363 0.2670326 2.307549    40 sample_H33_ChIP
#> 41 2.304252 0.2265529 2.403486    41 sample_H33_ChIP
#> 42 2.184851 0.2097188 1.877814    42 sample_H33_ChIP
#> 43 2.464733 0.2289928 2.594486    43 sample_H33_ChIP
#> 44 2.691855 0.2693417 2.539616    44 sample_H33_ChIP
#> 45 2.486217 0.2913635 2.108200    45 sample_H33_ChIP
#> 46 2.935590 0.3138470 2.815516    46 sample_H33_ChIP
#> 47 2.968804 0.3694925 2.486643    47 sample_H33_ChIP
#> 48 3.270211 0.4223941 2.603432    48 sample_H33_ChIP
#> 49 3.315354 0.3763270 2.815344    49 sample_H33_ChIP
#> 50 3.184501 0.3060692 2.718847    50 sample_H33_ChIP
#> 51 3.199918 0.3176035 3.066064    51 sample_H33_ChIP
#> 52 3.483032 0.3648622 3.841766    52 sample_H33_ChIP
#> 53 3.384507 0.3127053 3.642152    53 sample_H33_ChIP
#> 54 3.817368 0.3883696 3.355289    54 sample_H33_ChIP
#> 55 3.434996 0.3796849 3.155996    55 sample_H33_ChIP
#> 56 3.562528 0.3944338 3.124369    56 sample_H33_ChIP
#> 57 3.833116 0.4650076 3.173232    57 sample_H33_ChIP
#> 58 3.619649 0.4051178 3.533681    58 sample_H33_ChIP
#> 59 3.932313 0.4788887 3.552670    59 sample_H33_ChIP
#> 60 4.024829 0.4890062 3.387589    60 sample_H33_ChIP
#> 61 3.781616 0.4874027 2.603114    61 sample_H33_ChIP
#> 62 3.744094 0.4913049 3.453766    62 sample_H33_ChIP
#> 63 3.630780 0.4449057 3.087544    63 sample_H33_ChIP
#> 64 4.185833 0.5204512 3.224522    64 sample_H33_ChIP
#> 65 3.549521 0.4470859 2.678337    65 sample_H33_ChIP
#> 66 3.624105 0.4487285 3.117284    66 sample_H33_ChIP
#> 67 4.186247 0.4988139 3.857999    67 sample_H33_ChIP
#> 68 4.398364 0.5437542 3.731710    68 sample_H33_ChIP
#> 69 4.186880 0.5684283 4.162147    69 sample_H33_ChIP
#> 70 4.541604 0.6031076 4.152863    70 sample_H33_ChIP
#> 71 4.882680 0.6583484 3.330124    71 sample_H33_ChIP
#> 72 4.434294 0.6184024 2.805938    72 sample_H33_ChIP
#> 73 4.465311 0.5204311 3.211924    73 sample_H33_ChIP
#> 74 4.586631 0.5179597 3.784932    74 sample_H33_ChIP
#> 75 5.412189 0.6808004 4.478709    75 sample_H33_ChIP
#> 76 5.130822 0.5404398 5.056854    76 sample_H33_ChIP
#> 77 4.772371 0.5835512 3.497144    77 sample_H33_ChIP
#> 78 4.788891 0.5987570 4.242311    78 sample_H33_ChIP
#> 79 4.642610 0.5283055 4.316828    79 sample_H33_ChIP
#> 80 4.774207 0.5970659 4.409330    80 sample_H33_ChIP
#> 81 5.120912 0.7053851 3.702708    81 sample_H33_ChIP
#> 82 4.547538 0.5918970 3.625621    82 sample_H33_ChIP
#> 83 4.387126 0.6446382 3.304428    83 sample_H33_ChIP
#> 84 4.410434 0.5801318 3.330123    84 sample_H33_ChIP
#> 85 4.598745 0.5873002 3.309566    85 sample_H33_ChIP
#> 86 4.879927 0.6950630 2.482175    86 sample_H33_ChIP
#> 87 4.236072 0.5543600 2.795659    87 sample_H33_ChIP
#> 88 4.839732 0.6169402 3.612772    88 sample_H33_ChIP
#> 89 4.661699 0.6013412 2.931846    89 sample_H33_ChIP
#> 90 4.404561 0.5651088 3.810628    90 sample_H33_ChIP
#> 91 4.019680 0.4757994 3.006362    91 sample_H33_ChIP
#> 92 4.426952 0.5882985 4.100985    92 sample_H33_ChIP
#> 93 4.569928 0.5732725 3.224773    93 sample_H33_ChIP
#> 94 4.322334 0.4856785 3.736111    94 sample_H33_ChIP
#> 95 4.874971 0.6007702 3.954522    95 sample_H33_ChIP

bw_profile(bw, loci = bed, mode = "start",
           upstream = 1000, downstream = 1500)
#>        mean   sderror   median index          sample
#> 1  3.028204 0.4874518 2.471897     1 sample_H33_ChIP
#> 2  3.386287 0.6245647 2.235501     2 sample_H33_ChIP
#> 3  3.700873 0.7401184 2.556693     3 sample_H33_ChIP
#> 4  3.609102 0.5534036 3.032059     4 sample_H33_ChIP
#> 5  3.701056 0.5778502 3.003792     5 sample_H33_ChIP
#> 6  3.428869 0.5597654 2.405090     6 sample_H33_ChIP
#> 7  3.949016 0.6592606 2.949832     7 sample_H33_ChIP
#> 8  3.247899 0.4771775 2.590096     8 sample_H33_ChIP
#> 9  3.229546 0.4894977 2.058201     9 sample_H33_ChIP
#> 10 3.079227 0.4333806 2.497594    10 sample_H33_ChIP
#> 11 3.460437 0.3494176 2.996084    11 sample_H33_ChIP
#> 12 3.469246 0.2405814 2.898442    12 sample_H33_ChIP
#> 13 3.381148 0.3136780 3.032057    13 sample_H33_ChIP
#> 14 3.660127 0.3422368 3.309567    14 sample_H33_ChIP
#> 15 3.733174 0.2724274 3.535686    15 sample_H33_ChIP
#> 16 3.662145 0.3113339 3.353250    16 sample_H33_ChIP
#> 17 3.398767 0.3298550 2.978097    17 sample_H33_ChIP
#> 18 3.589647 0.3335419 3.399501    18 sample_H33_ChIP
#> 19 3.318746 0.3828422 3.006362    19 sample_H33_ChIP
#> 20 3.207521 0.3832718 2.733991    20 sample_H33_ChIP
#> 21 3.188432 0.3926917 2.620930    21 sample_H33_ChIP
#> 22 3.133921 0.3670644 2.649196    22 sample_H33_ChIP
#> 23 2.840810 0.3200884 2.477036    23 sample_H33_ChIP
#> 24 2.856044 0.2712727 2.895871    24 sample_H33_ChIP
#> 25 2.928174 0.2922429 2.644057    25 sample_H33_ChIP
```
