context("Test functions for bigWig handling")
library(GenomicRanges)
library(rtracklayer)
library(future)

# Setup -------------------------------------------

bed_to_bw <- function(bed, bw, chromsizes) {
  ranges <- import(bed)
  seqlengths(ranges) <- chromsizes
  export(ranges, bw)
}

bw1 <- tempfile("bigwig", fileext = ".bw")
bw2 <- tempfile("bigwig", fileext = ".bw")
bw3_zeros <- tempfile("bigwig", fileext = ".bw")
bw4_nas <- tempfile("bigwig", fileext = ".bw")
bg1 <- tempfile("bigwig_bg1", fileext = ".bw")
bg2 <- tempfile("bigwig_bg2", fileext = ".bw")
bg3_zeros <- tempfile("bigwig_bg3", fileext = ".bw")
bw_special <- tempfile("bigwig-2ñ", fileext = ".bw")

bed_with_names <- system.file("testdata", "labeled.bed", package = "wigglescout")
bed_with_names_na <- system.file("testdata", "labeled_na.bed", package = "wigglescout")
unnamed_bed <- system.file("testdata", "not_labeled.bed", package = "wigglescout")

tiles <- tileGenome(c(chr1 = 200, chr2 = 200),
                    tilewidth = 20,
                    cut.last.tile.in.chrom = TRUE)

granges <- import(bed_with_names)

setup({
  chromsizes <- c(200, 200)
  # Read bed files and transform these to bigwig
  bed_to_bw(system.file("testdata", "bed1.bed", package = "wigglescout"), bw1, chromsizes)
  bed_to_bw(system.file("testdata", "bed2.bed", package = "wigglescout"), bw2, chromsizes)
  bed_to_bw(system.file("testdata", "bed3.bed", package = "wigglescout"), bw3_zeros, chromsizes)
  bed_to_bw(system.file("testdata", "bed_na.bed", package = "wigglescout"), bw4_nas, chromsizes)
  bed_to_bw(system.file("testdata", "bg1.bed", package = "wigglescout"), bg1, chromsizes)
  bed_to_bw(system.file("testdata", "bg2.bed", package = "wigglescout"), bg2, chromsizes)
  bed_to_bw(system.file("testdata", "bg3.bed", package = "wigglescout"), bg3_zeros, chromsizes)
  bed_to_bw(system.file("testdata", "bg2.bed", package = "wigglescout"), bw_special, chromsizes)
})

teardown({
  unlink(bw1)
  unlink(bw2)
  unlink(bw3_zeros)
  unlink(bg1)
  unlink(bg2)
  unlink(bg3_zeros)
  unlink(bw4_nas)
  unlink(bw_special)
})

test_that("Setup files exist", {
  expect_true(file_test("-f", bw1))
  expect_true(file_test("-f", bw2))
  expect_true(file_test("-f", bw3_zeros))
  expect_true(file_test("-f", bg1))
  expect_true(file_test("-f", bg2))
  expect_true(file_test("-f", bg3_zeros))
  expect_true(file_test("-f", bed_with_names))
  expect_true(file_test("-f", bw_special))
})

# Core functions -----------------------------------------------

## bw_ranges ---------------------------------------------------

test_that(".bw_ranges returns a GRanges object", {
  bins <- .bw_ranges(bw1, tiles, per_locus_stat = "mean")
  expect_is(bins, "GRanges")

})

test_that(".bw_ranges returns correct values", {
  bins <- .bw_ranges(bw1, tiles, per_locus_stat = "mean")

  expect_equal(bins[1]$score, 1)
  expect_equal(bins[2]$score, 2)
  expect_equal(bins[10]$score, 10)
})

test_that("bw_ranges returns correct values on subset", {
  subset <- GRanges(seqnames = c("chr1"),
                    ranges = IRanges(10, 40))

  bins <- .bw_ranges(bw1, tiles, per_locus_stat = "mean", selection = subset)

  expect_equal(bins[1]$score, 1)
  expect_equal(bins[2]$score, 2)
  expect_equal(length(bins), 2)
})

## multi_bw_ranges ------------------------------------------------

test_that(".multi_bw_ranges returns correct values", {
  values <- .multi_bw_ranges(c(bw1, bw2), c("bw1", "bw2"), tiles)

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[1]$bw2, 20)
  expect_equal(values[2]$bw1, 2)
  expect_equal(values[2]$bw2, 19)
})

test_that(".multi_bw_ranges with zeros returns correct values", {
  subset <- GRanges(seqnames = c("chr2"), ranges = IRanges(1, 40))
  values <- .multi_bw_ranges(c(bw1, bw3_zeros), c("bw1", "bw3_zeros"), tiles,
                             selection = subset)

  expect_equal(values[1]$bw1, 11)
  expect_equal(values[1]$bw3_zeros, 0)
  expect_equal(values[2]$bw1, 12)
  expect_equal(values[2]$bw3_zeros, 0)
})

test_that(".multi_bw_ranges several processors returns correct values", {
  future::plan(multisession, workers=2)
  values <- .multi_bw_ranges(c(bw1, bw2), c("bw1", "bw2"), tiles)

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[1]$bw2, 20)
  expect_equal(values[2]$bw1, 2)
  expect_equal(values[2]$bw2, 19)
  future::plan(sequential)
})

test_that(".multi_bw_ranges returns correct values for single bigWig", {
  values <- .multi_bw_ranges(bw1, "bw1", tiles)

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[2]$bw1, 2)
})

test_that(".multi_bw_ranges returns correct values on subset", {
  subset <- GRanges(seqnames = c("chr1"), ranges = IRanges(30, 50))
  values <- .multi_bw_ranges(c(bw1, bw2),
                             c("bw1", "bw2"),
                             tiles,
                             selection = subset
  )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[1]$bw2, 19)
  expect_equal(values[2]$bw1, 3)
  expect_equal(values[2]$bw2, 18)
})

test_that(".multi_bw_ranges removes value even if quantile very small, due to interpolation", {
  # Note the use of bw1 twice. bw1 + bw2 means are always 10.5, so it would not
  # remove any rows.
  values <- .multi_bw_ranges(c(bw1, bw1),
                             c("bw1", "bw2"),
                             tiles,
                             remove_top = 0.01
  )

  expect_is(values, "GRanges")
  expect_equal(max(values$bw1), 19)
  expect_equal(max(values$bw2), 19)
})

test_that(".multi_bw_ranges removes percentile", {
  values <- .multi_bw_ranges(c(bw1, bw1),
                             c("bw1", "bw2"),
                             tiles,
                             remove_top = 0.05
  )

  expect_is(values, "GRanges")
  expect_equal(max(values$bw1), 19)
  expect_equal(max(values$bw2), 19)
})

test_that(".multi_bw_ranges removes percentile single column", {
  values <- .multi_bw_ranges(c(bw1),
                             c("bw1"),
                             tiles,
                             remove_top = 0.1
  )

  expect_is(values, "GRanges")
  expect_equal(max(values$bw1), 18)
})

## multi_bw_ranges_norm -------------------------------------------

test_that(
  ".multi_bw_ranges_norm with bwfiles == background returns all 1 values", {
    values <- .multi_bw_ranges_norm(c(bw1, bw2),
                                   bg_bwfilelist = c(bw1, bw2),
                                   c("bw1", "bw2"),
                                   tiles,
                                   norm_func = identity
    )

    expect_is(values, "GRanges")
    expect_equal(values[1]$bw1, 1)
    expect_equal(values[1]$bw2, 1)
    expect_equal(values[2]$bw1, 1)
    expect_equal(values[2]$bw2, 1)
  })

test_that(
  ".multi_bw_ranges_norm with background == 0 returns Infinite values", {
    values <- .multi_bw_ranges_norm(c(bw1, bw2),
                                    bg_bwfilelist = c(bg3_zeros, bg3_zeros),
                                    c("bw1", "bw2"),
                                    tiles,
                                    norm_func = identity
    )

    expect_is(values, "GRanges")
    expect_equal(values[1]$bw1, Inf)
    expect_equal(values[1]$bw2, Inf)
    expect_equal(values[2]$bw1, Inf)
    expect_equal(values[2]$bw2, Inf)
  })

test_that(
  ".multi_bw_ranges_norm with 0/0 returns NaN values", {
    values <- .multi_bw_ranges_norm(
      c(bw3_zeros, bw2),
      bg_bwfilelist = c(bg3_zeros, bg3_zeros),
      c("bw3_zeros", "bw2"),
      tiles,
      selection = GRanges(seqnames = c("chr1"), ranges = IRanges(1, 20)),
      norm_func = identity
    )

    expect_is(values, "GRanges")
    expect_equal(values[1]$bw3_zeros, NaN)
    expect_equal(values[1]$bw2, Inf)
  })

test_that(".multi_bw_ranges_norm returns correct values", {
  values <- .multi_bw_ranges_norm(c(bw1, bw2),
                                 bg_bwfilelist = c(bg1, bg2),
                                 c("bw1", "bw2"),
                                 tiles,
                                 norm_func = identity
  )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[2]$bw1, 2)
  expect_equal(values[1]$bw2, 10)
  expect_equal(values[2]$bw2, 9.5)
})

test_that(
  ".multi_bw_ranges_norm fails on background length not matching bwlist", {
    expect_error({
      values <- .multi_bw_ranges_norm(c(bw1, bw2),
                                     bg_bwfilelist = c(bg1),
                                     c("bw1", "bw2"),
                                     tiles,
                                     norm_func = identity
      )
    },
    "Background and signal bwfile lists must have the same length.")
  })

test_that(".multi_bw_ranges_norm removes percentile", {
  values <- .multi_bw_ranges_norm(c(bw1, bw1),
                            c(bg1, bg1),
                            c("bw1", "bw2"),
                            tiles,
                            remove_top = 0.1
  )

  expect_is(values, "GRanges")
  expect_equal(max(values$bw1), 18)
  expect_equal(max(values$bw2), 18)
})

test_that(".multi_bw_ranges_norm removes percentile single column", {
  values <- .multi_bw_ranges_norm(c(bw1), c(bg1),
                            c("bw1"),
                            tiles,
                            remove_top = 0.1
  )

  expect_is(values, "GRanges")
  expect_equal(max(values$bw1), 18)
})

## calculate_matrix_norm -------------------------------------

test_that("calculate_matrix_norm removes percentile", {
  values <- .calculate_matrix_norm(bw1, import(bed_with_names), bin_size = 2,
                                   upstream = 6, downstream = 6,
                                   remove_top = 0.05)

  expect_is(values, "matrix")
  expect_equal(nrow(values), 4)
  expect_equal(max(values), 17)

})


# Exported functions -----------------------------------------

## bw_global_coverage

test_that("bw_global_coverage returns correct value", {
    value <- bw_global_coverage(bw1)
    expect_equal(value, 10.5)
})

## bw_loci ---------------------------------------------------

test_that("bw_loci returns correct per locus values", {
  values <- bw_loci(bw1, bed_with_names, labels = "bw1", per_locus_stat = "mean")
  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[2]$bw1, 4.5)
})

test_that("bw_loci accepts GRanges objects", {
  values <- bw_loci(bw1, granges, labels = "bw1", per_locus_stat = "mean")
  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[2]$bw1, 4.5)
})

test_that("bw_loci returns correct per locus values on multiple files", {
  values <- bw_loci(c(bw1, bw2),
                    bed_with_names,
                    labels = c("bw1", "bw2"),
                    per_locus_stat = "mean"
  )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[2]$bw1, 4.5)
  expect_equal(values[1]$bw2, 19)
  expect_equal(values[2]$bw2, 16.5)
})

test_that(
  "bw_loci returns correct per locus values on multiple files, with bg", {
    values <- bw_loci(c(bw1, bw2),
                      bg_bwfiles = c(bg1, bg2),
                      bed_with_names,
                      labels = c("bw1", "bw2"),
                      per_locus_stat = "mean"
    )

    expect_is(values, "GRanges")
    expect_equal(values[1]$bw2, 9.5)
    expect_equal(values[2]$bw2, 8.25)
  })

test_that("bw_loci handles default names with special characters", {
  values <- bw_loci(bw_special,
                    bed_with_names,
                    aggregate_by = "true_mean"
  )

  expect_is(values, "data.frame")
})

test_that("bw_loci crashes on wrong number of labels for multiple files", {
  expect_error({
    values <- bw_loci(c(bw1, bw2), bed_with_names,
                      labels = "bw1",
                      per_locus_stat = "mean"
    )
  },
  "BigWig file list and column names must have the same length.")

})

test_that("bw_loci returns correct mean-of-means aggregated values", {
  values <- bw_loci(bw1,
                    bed_with_names,
                    labels = "bw1",
                    per_locus_stat = "mean",
                    aggregate_by = "mean"
  )

  expect_is(values, "data.frame")
  expect_equal(values["typeA", "bw1"], 7)
  expect_equal(values["typeB", "bw1"], 13.3333333333)
})

test_that("bw_loci returns correct true_mean aggregated values", {
  values <- bw_loci(bw1, bed_with_names,
                    labels = "bw1",
                    per_locus_stat = "mean",
                    aggregate_by = "true_mean"
  )

  expect_is(values, "data.frame")
  expect_equal(values["typeA", "bw1"], 7)
  expect_equal(values["typeB", "bw1"], 11.125)
})

test_that("bw_loci on an empty list throws an error", {
  expect_error({
    values <- bw_loci(c(), bed_with_names,
                      per_locus_stat = "mean",
                      aggregate_by = "true_mean"
    )
  },
  "File list provided is empty.")

})

test_that("bw_loci errors on non existing files on bwlist", {
  expect_error({
    values <- bw_loci(c(bw1, "invalidname.bw"),
                      bed_with_names,
                      per_locus_stat = "mean",
                      aggregate_by = "true_mean"
    )
  },
  "Files not found: invalidname.bw")
})

test_that("bw_loci returns correct median-of-means aggregated values", {
  values <- bw_loci(bw1, bed_with_names,
                    labels = "bw1",
                    per_locus_stat = "mean",
                    aggregate_by = "median"
  )

  expect_is(values, "data.frame")
  expect_equal(values["typeA", "bw1"], 7)
  expect_equal(values["typeB", "bw1"], 16.5)
})

test_that("bw_loci throws error on not implemented aggregate_by", {
  expect_error(
    values <- bw_loci(bw1, bed_with_names,
                      labels = "bw1",
                      per_locus_stat = "mean",
                      aggregate_by = "max"
    )
  )
})

test_that("bw_loci runs with background and aggregate_by parameter", {
  values <- bw_loci(bw1, bed_with_names, bg_bwfiles = bw2,
                    labels = "bw1",
                    per_locus_stat = "mean",
                    aggregate_by = "mean"
  )

  expect_is(values, "data.frame")
})

test_that("bw_loci runs with background == 0 and aggregate_by parameter", {
  values <- bw_loci(bw1, bed_with_names, bg_bwfiles = bg3_zeros,
                    labels = "bw1",
                    per_locus_stat = "mean",
                    aggregate_by = "mean"
  )

  expect_is(values, "data.frame")
  expect_equal(values["typeA", "bw1"], Inf)
  expect_equal(values["typeB", "bw1"], Inf)
})

test_that("bw_loci runs with 0/0 and aggregated values", {
  values <- bw_loci(bw3_zeros, bed_with_names, bg_bwfiles = bg3_zeros,
                    labels = "bw3_zeros",
                    per_locus_stat = "mean",
                    aggregate_by = "mean"
  )

  expect_equal(values["typeA", "bw3_zeros"], NaN)
})

test_that("bw_loci excludes NA values from true_mean aggregation", {
  values <- suppressWarnings(bw_loci(bw4_nas, bed_with_names,
                    labels = "bw4_nas",
                    per_locus_stat = "mean",
                    aggregate_by = "true_mean"
  ))

  expect_equal(values["typeA", "bw4_nas"], 12)
})

test_that("bw_loci excludes NA values from mean aggregation", {
  values <- suppressWarnings(bw_loci(bw4_nas, bed_with_names,
                                     labels = "bw4_nas",
                                     per_locus_stat = "mean",
                                     aggregate_by = "mean"
  ))

  expect_equal(values["typeA", "bw4_nas"], 12)
})

test_that("bw_loci excludes NA values from single locus", {
    # Here, exclude means that NAs are not counted as zeros, but overlapping
    # length that has NA value is not counted as valid length.
    # In this test, a locus 21-50 that has NA on 21-40 and 3 on 41-50 yields
    # a mean of 3.
    values <- suppressWarnings(bw_loci(bw4_nas, bed_with_names_na,
                                       labels = "bw4_nas",
                                       per_locus_stat = "mean")
    )
    expect_equal(values[1]$bw4_nas, 3)
})


test_that("bw_loci with default_na = 0 does NOT include NA values in mean calculation for single locus", {
    # NOTE: This test is here as documentation of somewhat counterintuitive
    # behavior: We could expect rtracklayer summary function to actually include
    # values in mean calculation if a different defaultValue to NA is used, but it
    # is not true. It still counts only the part of the sequence that has values,
    # excluding missing values from the length of the locus, which means the result
    # for a single loci is the same whether default_na is 0 or NA. default_na is
    # only presented when the FULL locus is empty (and a warning is also printed)

    # Here, a locus 21-50 that has NA on 21-40 and 3 on 41-50 STILL
    # yields a mean of 3.
    values <- suppressWarnings(bw_loci(bw4_nas, bed_with_names_na,
                                       labels = "bw4_nas",
                                       per_locus_stat = "mean",
                                       default_na = 0)
    )
    expect_equal(values[1]$bw4_nas, 3)
})


test_that("bw_loci fails if aggregate_by in an unnamed bed file", {
  expect_error({
    values <- bw_loci(bw1, unnamed_bed, bg_bwfiles = bw2,
                      labels = "bw1",
                      per_locus_stat = "mean",
                      aggregate_by = "mean"
    )},
    "missing values in 'row.names' are not allowed"
  )

})

## bw_bins ---------------------------------------------------

test_that("bw_bins returns correct per locus values", {
  values <- bw_bins(bw1,
                    selection = import(bed_with_names),
                    labels = "bw1",
                    per_locus_stat = "mean"
  )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 5.5)
  expect_equal(values[2]$bw1, 15.5)
})

test_that("bw_bins returns 1 when bwfile == bg_bwfile", {
  values <- bw_bins(bw1,
                    bg_bwfiles = bw1,
                    selection = import(bed_with_names),
                    labels = "bw1",
                    per_locus_stat = "mean"
  )

  expect_is(values, 'GRanges')
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[2]$bw1, 1)
})

## bw_profile -------------------------

test_that("bw_profile on an empty list throws an error", {
  expect_error({
    values <- bw_profile(c(),
                         bed_with_names,
                         labels = NULL
    )
  },
  "File list provided is empty.")

})

test_that("bw_loci on non-existing bed file throws an error", {
  expect_error({
    values <- bw_loci(bw1,
                     "invalidname.bed",
                     per_locus_stat = "mean",
                     aggregate_by = "true_mean"
    )
  },
  "Files not found: invalidname.bed")
})

test_that("bw_profile on non-existing bed file throws an error", {
  expect_error({
    values <- bw_profile(bw1,
                         loci = "invalidname.bed",
                         labels = NULL
    )
  },
  "Files not found: invalidname.bed")
})

test_that("bw_profile errors on non existing files on bwlist", {
  expect_error({
    values <- bw_profile(c(bw1, "invalidname.bw"),
                         loci = bed_with_names,
                         labels = NULL
    )
  },
  "Files not found: invalidname.bw")
})

test_that("bw_profile runs quiet on valid parameters", {
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         loci = bed_with_names,
                         upstream = 1,
                         downstream = 1,
                         bin_size = 1
    )
  })

})

test_that("bw_profile runs on GRanges object", {
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         loci = granges,
                         upstream = 1,
                         downstream = 1,
                         bin_size = 1
    )
  })
})

test_that("bw_profile runs quiet on valid parameters, mode start", {
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         loci = bed_with_names,
                         upstream = 1,
                         downstream = 1,
                         bin_size = 1,
                         mode = "start"
    )
  })
})

test_that(
  "bw_profile runs quiet on valid parameters, mode start, with background", {
    expect_silent({
      values <- bw_profile(c(bw1, bw2),
                           bg_bwfiles = c(bg1, bg2),
                           loci = bed_with_names,
                           upstream = 1,
                           downstream = 1,
                           bin_size = 1,
                           mode = "start"
      )
    })
})

test_that("bw_profile normalized returns 1 when fg == bg", {
  values <- bw_profile(c(bw1, bw2),
                       bg_bwfiles = c(bw1, bw2),
                       loci = bed_with_names,
                       upstream = 1,
                       downstream = 1,
                       bin_size = 1,
                       mode = "start"
  )

  expect_is(values, "data.frame")
  expect_equal(values[1, "mean"], 1)
  expect_equal(values[2, "mean"], 1)
})

test_that("bw_profile fails if labels and bwfiles have different length", {
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                         loci = bed_with_names,
                         labels = c("bw1"),
                         upstream = 1,
                         downstream = 1,
                         bin_size = 1,
                         mode = "start"
    )
  },
  "labels and bwfiles must have the same length")
})

test_that("bw_profile runs quiet on valid parameters, mode end", {
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         loci = bed_with_names,
                         upstream = 1,
                         downstream = 1,
                         bin_size = 1,
                         mode = "end"
    )
  })
})

test_that("bw_profile runs quiet on valid parameters, middle parameter", {
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         loci = bed_with_names,
                         upstream = 1,
                         downstream = 1,
                         middle = 1,
                         bin_size = 1,
                         mode = "end"
    )
  })
})

test_that("bw_profile runs quiet on valid parameters with background", {
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         bg_bwfiles = c(bg1, bg2),
                         loci = bed_with_names,
                         upstream = 1,
                         downstream = 1,
                         bin_size = 1
    )
  })

})

test_that("bw_profile throws error on flanking region smaller than bin size", {
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                         loci = bed_with_names,
                         upstream = 1,
                         downstream = 1,
                         bin_size = 10
    )
  },
  "bin size must be smaller than flanking regions")

})

test_that("bw_profile throws error on negative bin size", {
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                         loci = bed_with_names,
                         upstream = 1,
                         downstream = 1,
                         bin_size = -10
    )
  },
  "bin size must be a positive value: -10")

})

test_that("bw_profile throws error on negative upstream value", {
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                         loci = bed_with_names,
                         upstream = -10,
                         downstream = 10,
                         bin_size = 10
    )
  },
  "upstream size must be a positive value: -10")
})

test_that("bw_profile throws error on negative downstream value", {
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                         loci = bed_with_names,
                         upstream = 10,
                         downstream = -10,
                         bin_size = 10
    )
  },
  "downstream size must be a positive value: -10")
})

## bw_heatmap -------------------------------------------

test_that("bw_heatmap returns correct values odd bin size", {
  values <- bw_heatmap(bw1,
                       bg_bwfiles = NULL,
                       loci = bed_with_names,
                       upstream = 5,
                       downstream = 5,
                       bin_size = 2,
                       mode = "start"
  )

  expect_is(values[[1]], "matrix")
  expect_equal(nrow(values[[1]]), 5)
  expect_equal(ncol(values[[1]]), 4)

  expect_equal(values[[1]][1, ], c(rep(1,2), rep(2,2)))
  expect_equal(values[[1]][2, ], c(rep(3,2), rep(4,2)))
  expect_equal(values[[1]][3, ], c(rep(11, 2), rep(12, 2)))
  expect_equal(values[[1]][4, ], rep(16, 4))
  expect_equal(values[[1]][5, ], c(rep(18,2), rep(19,2)))

})

test_that("bw_heatmap returns correct values odd bin size stretch", {
  values <- bw_heatmap(bw1,
                       bg_bwfiles = NULL,
                       loci = bed_with_names,
                       upstream = 5,
                       downstream = 5,
                       middle = 4,
                       bin_size = 2,
                       mode = "stretch"
  )

  expect_is(values[[1]], "matrix")
  expect_equal(nrow(values[[1]]), 5)
  expect_equal(ncol(values[[1]]), 6)

  expect_equal(values[[1]][1, ], c(1, 1, 2, 2, 3, 3))
  expect_equal(values[[1]][2, ], c(3, 3, 4, 5, 6, 6))
  expect_equal(values[[1]][3, ], c(11,11,12,12,13,13))
  expect_equal(values[[1]][4, ], c(16,16,16,17,17,17))
  expect_equal(values[[1]][5, ], c(18,18,19,19,20,20))

})

test_that("bw_heatmap returns correct values", {
  values <- bw_heatmap(bw1,
                       bg_bwfiles = NULL,
                       loci = bed_with_names,
                       upstream = 5,
                       downstream = 5,
                       bin_size = 1,
                       mode = "start"
  )

  expect_is(values[[1]], "matrix")
  expect_equal(nrow(values[[1]]), 5)
  expect_equal(ncol(values[[1]]), 10)

  expect_equal(values[[1]][1, ], c(rep(1,5), rep(2,5)))
  expect_equal(values[[1]][2, ], c(rep(3,5), rep(4,5)))
  expect_equal(values[[1]][3, ], c(rep(11, 5), rep(12, 5)))
  expect_equal(values[[1]][4, ], rep(16, 10))
  expect_equal(values[[1]][5, ], c(rep(18,5), rep(19,5)))

})

test_that("bw_heatmap returns correct values on GRanges object", {
  values <- bw_heatmap(bw1,
                       bg_bwfiles = NULL,
                       loci = granges,
                       upstream = 5,
                       downstream = 5,
                       bin_size = 1,
                       mode = "start"
  )

  expect_is(values[[1]], "matrix")
  expect_equal(nrow(values[[1]]), 5)
  expect_equal(ncol(values[[1]]), 10)

  expect_equal(values[[1]][1, ], c(rep(1,5), rep(2,5)))
  expect_equal(values[[1]][2, ], c(rep(3,5), rep(4,5)))
  expect_equal(values[[1]][3, ], c(rep(11, 5), rep(12, 5)))
  expect_equal(values[[1]][4, ], rep(16, 10))
  expect_equal(values[[1]][5, ], c(rep(18,5), rep(19,5)))

})

test_that("bw_heatmap with bg returns 1 when fg == bg", {
  values <- bw_heatmap(bw1,
                       bg_bwfiles = bw1,
                       loci = bed_with_names,
                       upstream = 5,
                       downstream = 5,
                       bin_size = 1,
                       mode = "start"
  )

  expect_is(values[[1]], "matrix")
  expect_equal(nrow(values[[1]]), 5)
  expect_equal(ncol(values[[1]]), 10)

  expect_equal(values[[1]][1, ], c(rep(1,10)))
  expect_equal(values[[1]][2, ], c(rep(1,10)))
  expect_equal(values[[1]][3, ], c(rep(1,10)))
  expect_equal(values[[1]][4, ], c(rep(1,10)))
  expect_equal(values[[1]][5, ], c(rep(1,10)))
})

## build_bins --------------------------------------------

test_that("build_bins crashes on unknown or not included genome", {
  expect_error({
    build_bins(bin_size = 10000, genome = "not_a_genome")
  })
})

test_that("build_bins runs for mm9", {
  values <- build_bins(bin_size = 50000, genome = "mm9")
  expect_is(values, "GRanges")
})

test_that("build_bins creates bins of correct size", {
  values <- build_bins(bin_size = 50000, genome = "mm9")
  expect_is(values, "GRanges")
  expect_equal(width(ranges(values[1])), 50000)
})

test_that("build_bins runs for hg38", {
  values <- build_bins(bin_size = 50000, genome = "hg38")
  expect_is(values, "GRanges")
})
