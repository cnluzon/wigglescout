# Check that a list of bigWig files have the same Seqinfo as an external object (usually a pre-built tile set). Throws a warning if they don't.

Check that a list of bigWig files have the same Seqinfo as an external
object (usually a pre-built tile set). Throws a warning if they don't.

## Usage

``` r
.validate_references_with_external_seqinfo(bwfiles, seqinfo, show_max = 25)
```

## Arguments

- bwfiles:

  List of bigWig files

- seqinfo:

  Seqinfo object

- show_max:

  Max number of shared references to show in the warning

## Value

TRUE if they share all references, FALSE otherwise.
