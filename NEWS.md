# wigglescout 0.15.4

* Fixed warnings from deprecations in ggplot - removed all appearances of 
aes_string and replaced those with the recommended .data$field notation instead.
Additionally, incorporated legend.position = inside to the profile plots.
* wigglescout now requires ggplot2 >= 3.5.0 due to these changes.
* Look of figures improved by better positioning of legend, more adequate base
fontsize and truncation of too long labels (>35 characters). Legend title in
summary heatmaps is also shortened (background renamed as bg).
* Bug fix #68 - colors are now correctly assigned to labels in plot_bw_profile.

# wigglescout 0.15.3

* Testing refactor to exclude as much global variables as possible. Now some of
them have been moved to fixture rds files that are loaded on demand, and also
some of the necessary files generated on the fly per test function call and
removed after with `withr`. This might generate a bit of overhead but improves
legibility and maintainability of the tests.

# wigglescout 0.15.2

* Fixed #96: Inconsistency when bw_loci was summarized, normalized to background
 and remove_top > 0. Now loci selection for exclusion is done at the same time
 to both background and main samples.
* Fixed deprecation warnings for tidyselect 1.1.0, 1.2.0, ggplot 3.4.0:
 aes_string replaced by aes as tidyselect now goes for strings instead of 
 .data[["string"]], as explained on: https://github.com/tidyverse/tidyverse.org/pull/600
 This impacts ggplot2, and there is no easy fix, since .data[[col]] is
 also deprecated.
 All calls for `dplyr::select` that use vectors have been replaced by the
 relevant `tidyselect::all_of`, `tidyselect::any_of` use.
* Extra deprecation warnings for ggplot 3.4.0: linewidth replaces size for
 line and rectangle functions.
* The versions of the packages where this is deprecated are added 
 on the DESCRIPTION file, since these are breaking changes for older
 versions.
* Added a warning when heatmap colorscale is min == max. This can happen due
 to extremely sparse heatmaps. For consistency, if this happens it is plotted
 using raw colorscale (instead quantiles (0.01, 0.99), use (min, max)).

# wigglescout 0.15.1

* Reimplemented granges_cbind function to make use of dplyr::left_join instead,
 improving robustness of the approach. Ranges are also deduplicated before
 merging. In the case where a granges locus has more than one name, the final
 granges object will have as many rows as names.
* .loci_to_granges now does not sort the sequences, since it is not necessary
 for having a proper merging and keeps consistency.
* .loci_to_granges now gets rid of any column that is not a name, so it does
 not interfere with scoring naming.

# wigglescout 0.15.0

* Added bw_global_mean function. Takes a bigwig file and returns global mean
 coverage by taking the per-contig summary from rtracklayer, multiplying by
 contig length and dividing the sum by global length.
* Added a scaling parameter that divides results by global mean calculated by
 bw_global_mean function, to simulate 1x genome coverage tracks. This affects
 all the calculating bw_ functions and the plotting functions: plot_bw_ but it
 is backwards compatible (default is not to do anything).
* [BUG fix]: Added a sorting in .multi_bw_granges function to make sure names
 are assigned correctly. This affected only instances where BED file provided
 was not sorted. Might have been introduced recently, since older calculations
 are still correct.

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
