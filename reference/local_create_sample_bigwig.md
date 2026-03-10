# Create a self-destructing local temporary bigwig based on a BED file and chromsizes array.

Create a self-destructing local temporary bigwig based on a BED file and
chromsizes array.

## Usage

``` r
local_create_sample_bigwig(
  bed,
  chromsizes,
  f = tempfile(fileext = ".bw"),
  env = parent.frame()
)
```

## Arguments

- bed:

  BED filename

- chromsizes:

  Vector with length of chromosomes

- f:

  bigWig filename

- env:

  Running environment

## Value

Local filename
