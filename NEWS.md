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
