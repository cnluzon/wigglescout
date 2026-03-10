# Helper function for matrix heatmap plot

Helper function for matrix heatmap plot

## Usage

``` r
.heatmap_body(values, zmin, zmax, cmap)
```

## Arguments

- values:

  Matrix with values

- zmin:

  Minimum of the color scale. Majority of tools set this to 0.01
  percentile of the data distribution.

- zmax:

  Maximum of the color scale. Majority of tools set this to 0.99.

- cmap:

  Color map. Any RColorBrewer palette name is accepted here.

## Value

Named list plot and calculated values
