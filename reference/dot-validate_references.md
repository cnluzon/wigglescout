# Check that a list of bigWig files have the same Seqinfo. Throws a warning if they don't.

Check that a list of bigWig files have the same Seqinfo. Throws a
warning if they don't.

## Usage

``` r
.validate_references(bwfiles, show_max = 25)
```

## Arguments

- bwfiles:

  List of bigWig files

- show_max:

  Max number of references to show in the warning

## Value

TRUE if they share all references, FALSE otherwise.
