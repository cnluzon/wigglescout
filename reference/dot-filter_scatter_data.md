# Filter scatterplot data with quantile threshold on both axes

Filter scatterplot data with quantile threshold on both axes

## Usage

``` r
.filter_scatter_data(x, y, remove_top, col_x = "score", col_y = "score")
```

## Arguments

- x:

  GRanges for x axis

- y:

  GRanges for y axis

- remove_top:

  Return range 0-(1-remove_top). By default returns the whole
  distribution (remove_top == 0).

- col_x:

  Name of column in x to filter

- col_y:

  Name of column in y to filter

## Value

Named list, where gr is a GRanges object and stats a named list of
values: points: number of points in the figure, NA.x = Number of NA
values on the x axis. filtered.x: number of points that were higher than
the threshold on the x axis. quant.x: Quantile used as a threhsold on
the x axis. Same for the y values.
