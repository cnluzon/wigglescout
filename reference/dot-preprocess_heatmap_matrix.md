# Preprocess heatmap values for plotting

This includes capping on zmin / zmax values. If not provided,
percentiles 0.01 and 0.99 are used, to match plotting of other similar
tools. If provided, a specific order of rows is used. Otherwise, rows
are sorted by rowMeans. If the number of rows in the matrix is larger
than a certain max_rows_allowed value, number of rows is reduced to
max_rows_allowed. Rather than subsampling each row then will represent
an average of the underlying rows.

## Usage

``` r
.preprocess_heatmap_matrix(values, zmin, zmax, order_by, max_rows_allowed)
```

## Arguments

- values:

  Matrix to preprocess

- zmin:

  Minimum value to show in color.

- zmax:

  Maximum value to show in color.

- order_by:

  Row order (or NULL).

- max_rows_allowed:

  Maximum heatmap resolution.

## Value

A named list with values and stats calculated.
