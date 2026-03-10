# Aggregate scores of a GRanges object on a field

Aggregates scores of a GRanges object on a specific field. Used for
summary functions.

## Usage

``` r
.aggregate_scores(scored_granges, group_col, aggregate_by)
```

## Arguments

- scored_granges:

  A GRanges object with numerical metadata columns

- group_col:

  A column among the mcols that can be seen as a factor.

- aggregate_by:

  Function used to aggregate: mean, median, true_mean.

      true_mean: Mean coverage taking all elements in a class as one large bin.

      mean: mean-of-distribution approach. The mean of the aggregated value per
        locus is reported.

      median: median-of-distribution. The median of the aggregated value per
        locus is reported.

  This does not pass the checks but does not break the build process
  either so I left it as is, so the note thrown by build is noted
  everytime.

## Value

A data frame with aggregated scores.
