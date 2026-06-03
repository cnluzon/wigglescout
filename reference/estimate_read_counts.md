# Estimate the number of reads in a locus

This is a helper tool to back-calculate a number that represents the raw
counts in a locus. It is an estimate so one should consider several
factors: a) fraglen parameter needs to be accurate; b) fraglen is a
constant, so if the original distribution of fragment length was very
heterogeneous, it will over/under estimate locus where fragments are
shorter or longer than the average. c) If there was a scaling done in
the original bigWig, this number will not have a 1 to 1 correspondence
with the original number of reads.

## Usage

``` r
estimate_read_counts(mean_cov, width, fraglen)
```

## Arguments

- mean_cov:

  Mean coverage

- width:

  Width of the locus

- fraglen:

  Fragment length

## Value

An integer representing the estimated number of reads (mean_cov\*width)
/ fragment_length
