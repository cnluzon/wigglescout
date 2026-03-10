# Calculate overlap of one GRanges with a main GRanges object

Calculate overlap of one GRanges with a main GRanges object

## Usage

``` r
.multi_ranges_overlap(main_ranges, other_ranges, labels, minoverlap)
```

## Arguments

- main_ranges:

  Main GRanges to which overlap is calculated.

- other_ranges:

  A list of GRanges objects

- labels:

  Labels to each of the other_ranges.

- minoverlap:

  Minimum overlap to consider an overlap.

## Value

A data.frame in tall format with the values of the overlapping loci.
Loci returned belong to other_ranges, NOT main_ranges. A group field is
added as factor.
