# wigglescout 0.15.0

* Added bw_global_mean function. Takes a bigwig file and returns global mean
 coverage by taking the per-contig summary from rtracklayer, multiplying by
 contig length and dividing by global length.

# wigglescout 0.14.0

* Added a default_na parameter to all the functions. This makes use of 
 rtracklayer summary function defaultValue parameter. Bins or loci that have
 no values are replaced by the given value instead. This is particularly useful
 with very sparse bigWig files, where .aggregate_scores might kill majority
 of rows due to NA values.

# wigglescout 0.13.9

* .aggregate_scores is robust to presence of NAs. Only rows with no NA
 values are aggregated. This affect results from plot_bw_summary_heatmap and
 bw_loci when aggregate_by is not NULL.

# wigglescout 0.13.8

* .bw_ranges refactored to make use of .fetch_bigwig function.

# wigglescout 0.13.7

* bw_profile bug fix: remove_top now ignores NA values for calculating rowMeans
and quantile.
* Documentation fix: Removed unnecessary packages from README.

# wigglescout 0.13.6

* Refactor plot functions to separate outlier filtering from plotting and
reporting.
* reshape2 dependency melt is removed in favor of tidyr pivot_longer.
* Explicit imports of ggplot2 functions.
* Better handling of decimal positions in verbose plots.
* Remove soft-deprecated aes_string calls and replace with .data.
* Row aggregation is now more precise in plot_bw_heatmap.
* Missing tick 0 on heatmap on row aggregation fixed back. Replaced by 1 value.

# wigglescout 0.13.5

* Code style fixes:
    - All lines under 80 characters.
    - Replace sapply() -> vapply().
    - Replace calls to class() to compare with is(obj, "class").
    - Replace T/F values with TRUE/FALSE.
    - Calls like warning(), stop() do not contain paste() calls or
      error keywords.
    - Indent 4 spaces all.
    - Plotting functions now dynamically match inner parameters. Tests for 
      inner parameter calls are hence removed.

# wigglescout 0.13.4

* Added a `NEWS.md` file to track changes to the package.
* Fixed `DESCRIPTION` to comply with Bioconductor's `BiocCheck`:
    - Added `biocViews` terms.
    - Fixed inconsistencies in imports/suggests packages.
    - Expanded description of the package.
    - Correctly formated title of the tool.
* Made all non-exported functions start with dots.
* Added `sessionInfo()` to the end of Vignettes.
* Added examples to all exported functions.
