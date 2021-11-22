# wigglescout 0.13.5

* Code style fixes:
    - All lines under 80 characters.
    - Replace sapply() -> vapply().
    - Replace calls to class() to compare with is(obj, "class").
    - Replace T/F values with TRUE/FALSE.
    - Calls like warning(), stop() do not contain paste() calls or
      error keywords.

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
