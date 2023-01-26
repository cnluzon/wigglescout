context("Test functions for bigWig handling")
library(rtracklayer)
library(future)

# Setup -------------------------------------------
chromsizes <- c(200, 200)

# Core functions -----------------------------------------------

## bw_ranges ---------------------------------------------------

test_that(".bw_ranges returns a GRanges object", {
  bw <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bins <- .bw_ranges(bw, make_test_tiles(), per_locus_stat = "mean")
  expect_is(bins, "GRanges")
})

test_that(".bw_ranges returns correct values", {
  bw <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bins <- .bw_ranges(bw, make_test_tiles(), per_locus_stat = "mean")

  expect_equal(bins[1]$score, 1)
  expect_equal(bins[2]$score, 2)
  expect_equal(bins[10]$score, 10)
})

test_that(".bw_ranges returns correct values when scaled", {
    bw <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    bins <- .bw_ranges(bw, make_test_tiles(), per_locus_stat = "mean", scaling = "relative")

    expect_equal(bins[1]$score, 0.0952381, tolerance = 0.00001)
    expect_equal(bins[2]$score, 0.1904762, tolerance = 0.00001)
    expect_equal(bins[10]$score, 0.952381, tolerance = 0.00001)
})

test_that("bw_ranges returns correct values on subset", {
  subset <- GRanges(seqnames = c("chr1"),
                    ranges = IRanges(10, 40))
  bw <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bins <- .bw_ranges(bw, make_test_tiles(), per_locus_stat = "mean", selection = subset)

  expect_equal(bins[1]$score, 1)
  expect_equal(bins[2]$score, 2)
  expect_equal(length(bins), 2)
})

## multi_bw_ranges ------------------------------------------------

test_that(".multi_bw_ranges returns correct values", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  values <- .multi_bw_ranges(c(bw1, bw2), c("bw1", "bw2"), make_test_tiles())

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[1]$bw2, 20)
  expect_equal(values[2]$bw1, 2)
  expect_equal(values[2]$bw2, 19)
})

test_that(".multi_bw_ranges with sorted ranges same as shuffled ranges and sort after", {
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)

    gr_shuffled <- rtracklayer::import(get_testfile("labeled_shuffled.bed"))
    gr <- rtracklayer::import(get_testfile("labeled.bed"))
    shuffled <- .multi_bw_ranges(c(bw1, bw2), c("bw1", "bw2"), gr_shuffled)
    values <- .multi_bw_ranges(c(bw1, bw2), c("bw1", "bw2"), gr)

    expect_is(values, "GRanges")
    re_sorted <- sortSeqlevels(shuffled)
    re_sorted <- sort(re_sorted)
    expect_equal(values, re_sorted)
})


test_that(".multi_bw_ranges with names returns correct values in same order as input", {
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
    gr_shuffled <- rtracklayer::import(get_testfile("labeled_shuffled.bed"))
    values <- .multi_bw_ranges(c(bw1, bw2), c("bw1", "bw2"), gr_shuffled)

    expect_is(values, "GRanges")

    expect_equal(values[1]$bw1, 4.5)
    expect_equal(start(ranges(values[1])), 61)
    expect_equal(end(ranges(values[1])), 100)
    expect_true(data.frame(values)[1, "seqnames"] == "chr1")
    expect_equal(values[1]$name, "typeB")

    expect_equal(values[2]$bw1, 2)
    # 1-based
    expect_equal(start(ranges(values[2])), 21)
    expect_equal(end(ranges(values[2])), 40)
    expect_true(data.frame(values)[2, "seqnames"] == "chr1")
    expect_equal(values[2]$name, "typeA")

})

test_that(".multi_bw_ranges with zeros returns correct values", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed3.bed"), chromsizes)
  subset <- GRanges(seqnames = c("chr2"), ranges = IRanges(1, 40))
  values <- .multi_bw_ranges(c(bw1, bw2), c("bw1", "bw_zeros"), make_test_tiles(),
                             selection = subset)

  expect_equal(values[1]$bw1, 11)
  expect_equal(values[1]$bw_zeros, 0)
  expect_equal(values[2]$bw1, 12)
  expect_equal(values[2]$bw_zeros, 0)
})

test_that(".multi_bw_ranges several processors returns correct values", {
  future::plan(multisession, workers=2)
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  values <- .multi_bw_ranges(c(bw1, bw2), c("bw1", "bw2"), make_test_tiles())

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[1]$bw2, 20)
  expect_equal(values[2]$bw1, 2)
  expect_equal(values[2]$bw2, 19)
  future::plan(sequential)
})

test_that(".multi_bw_ranges returns correct values for single bigWig", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- .multi_bw_ranges(bw1, "bw1", make_test_tiles())

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[2]$bw1, 2)
})

test_that(".multi_bw_ranges returns correct values on subset", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  subset <- GRanges(seqnames = c("chr1"), ranges = IRanges(30, 50))
  values <- .multi_bw_ranges(c(bw1, bw2),
                             c("bw1", "bw2"),
                             make_test_tiles(),
                             selection = subset
  )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[1]$bw2, 19)
  expect_equal(values[2]$bw1, 3)
  expect_equal(values[2]$bw2, 18)
})

test_that(".multi_bw_ranges removes value even if quantile very small, due to interpolation", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  # Note the use of bw1 twice. bw1 + bw2 means are always 10.5, so it would not
  # remove any rows.
  values <- .multi_bw_ranges(c(bw1, bw1),
                             c("bw1", "bw2"),
                             make_test_tiles(),
                             remove_top = 0.01
  )

  expect_is(values, "GRanges")
  expect_equal(max(values$bw1), 19)
  expect_equal(max(values$bw2), 19)
})

test_that(".multi_bw_ranges removes percentile", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- .multi_bw_ranges(c(bw1, bw1),
                             c("bw1", "bw2"),
                             make_test_tiles(),
                             remove_top = 0.05
  )

  expect_is(values, "GRanges")
  expect_equal(max(values$bw1), 19)
  expect_equal(max(values$bw2), 19)
})

test_that(".multi_bw_ranges removes percentile single column", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- .multi_bw_ranges(c(bw1),
                             c("bw1"),
                             make_test_tiles(),
                             remove_top = 0.1
  )

  expect_is(values, "GRanges")
  expect_equal(max(values$bw1), 18)
})

## multi_bw_ranges_norm -------------------------------------------

test_that(
  ".multi_bw_ranges_norm with bwfiles == background returns all 1 values", {
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
    values <- .multi_bw_ranges_norm(c(bw1, bw2),
                                   bg_bwfilelist = c(bw1, bw2),
                                   c("bw1", "bw2"),
                                   make_test_tiles(),
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
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
    bg_zeros <- local_create_sample_bigwig(get_testfile("bg3.bed"), chromsizes)

    values <- .multi_bw_ranges_norm(c(bw1, bw2),
                                    bg_bwfilelist = c(bg_zeros, bg_zeros),
                                    c("bw1", "bw2"),
                                    make_test_tiles(),
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
    bw1 <- local_create_sample_bigwig(get_testfile("bed3.bed"), chromsizes)
    bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
    bg_zeros <- local_create_sample_bigwig(get_testfile("bg3.bed"), chromsizes)
    values <- .multi_bw_ranges_norm(
      c(bw1, bw2),
      bg_bwfilelist = c(bg_zeros, bg_zeros),
      c("bw_zeros", "bw2"),
      make_test_tiles(),
      selection = GRanges(seqnames = c("chr1"), ranges = IRanges(1, 20)),
      norm_func = identity
    )

    expect_is(values, "GRanges")
    expect_equal(values[1]$bw_zeros, NaN)
    expect_equal(values[1]$bw2, Inf)
  })

test_that(".multi_bw_ranges_norm returns correct values", {

  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  bg1 <- local_create_sample_bigwig(get_testfile("bg1.bed"), chromsizes)
  bg2 <- local_create_sample_bigwig(get_testfile("bg2.bed"), chromsizes)

  values <- .multi_bw_ranges_norm(c(bw1, bw2),
                                 bg_bwfilelist = c(bg1, bg2),
                                 c("bw1", "bw2"),
                                 make_test_tiles(),
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
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  bg1 <- local_create_sample_bigwig(get_testfile("bg1.bed"), chromsizes)
  expect_error({
      values <- .multi_bw_ranges_norm(
          c(bw1, bw2),
          bg_bwfilelist = c(bg1),
          c("bw1", "bw2"),
          make_test_tiles(),
          norm_func = identity
      )
  },
  "Background and signal bwfile lists must have the same length.")
})

test_that(".multi_bw_ranges_norm removes percentile", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bg1 <- local_create_sample_bigwig(get_testfile("bg1.bed"), chromsizes)
  values <- .multi_bw_ranges_norm(c(bw1, bw1),
                            c(bg1, bg1),
                            c("bw1", "bw2"),
                            make_test_tiles(),
                            remove_top = 0.1
  )

  expect_is(values, "GRanges")
  expect_equal(max(values$bw1), 18)
  expect_equal(max(values$bw2), 18)
})

test_that(".multi_bw_ranges_norm removes percentile single column", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bg1 <- local_create_sample_bigwig(get_testfile("bg1.bed"), chromsizes)
  values <- .multi_bw_ranges_norm(c(bw1), c(bg1),
                            c("bw1"),
                            make_test_tiles(),
                            remove_top = 0.1
  )

  expect_is(values, "GRanges")
  expect_equal(max(values$bw1), 18)
})

## calculate_matrix_norm -------------------------------------

test_that("calculate_matrix_norm removes percentile", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- .calculate_matrix_norm(bw1,
                                   import(get_testfile("labeled.bed")),
                                   bin_size = 2,
                                   upstream = 6, downstream = 6,
                                   remove_top = 0.05)

  expect_is(values, "matrix")
  expect_equal(nrow(values), 4)
  expect_equal(max(values), 17)

})


# Exported functions -----------------------------------------

## bw_global_coverage ----------------------------------------

test_that("bw_global_coverage returns correct value", {
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    value <- bw_global_coverage(bw1)
    expect_equal(value, 10.5)
})

## bw_loci ---------------------------------------------------

test_that("bw_loci returns correct per locus values", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_loci(bw1, get_testfile("labeled.bed"), labels = "bw1", per_locus_stat = "mean")
  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[2]$bw1, 4.5)
})

test_that("bw_loci accepts GRanges objects", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  granges <- rtracklayer::import(get_testfile("labeled.bed"))
  values <- bw_loci(bw1, granges, labels = "bw1", per_locus_stat = "mean")
  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[2]$bw1, 4.5)
})

test_that("bw_loci does not crash when a scored bed file is provided as loci", {
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)

    # bed files for generating bigWig have a score, which should be ignored
    values <- bw_loci(bw1, get_testfile("bed1.bed"), labels = "bw1", per_locus_stat = "mean")
    expect_is(values, "GRanges")
})

test_that("bw_loci with multiple columns, name and duplicated loci returns correct dimensions", {
    gr <- rtracklayer::import(get_testfile("labeled.bed"))
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
    # repeat one value
    gr <- append(gr, gr[length(gr)])
    values <- bw_loci(c(bw1, bw2), gr, labels = c("bw1","bw2"), per_locus_stat = "mean")
    expect_equal(length(values), length(gr)-1)
    expect_equal(ncol(mcols(values)), 3)
})

test_that("bw_loci with duplicated locus returns as many rows as unique loci", {
    gr <- rtracklayer::import(get_testfile("bed1.bed"))
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    # repeat one value
    gr <- append(gr, gr[length(gr)])
    values <- bw_loci(bw1, gr, labels = "bw1", per_locus_stat = "mean")
    expect_equal(length(values), length(gr)-1)
})


test_that("bw_loci with named loci and multi named locus returns as many rows as different names", {
    gr <- rtracklayer::import(get_testfile("bed1.bed"))
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    # repeat one value
    gr <- append(gr, gr[length(gr)])
    gr$name <- paste0("name", seq(1, length(gr)))

    values <- bw_loci(bw1, gr, labels = "bw1", per_locus_stat = "mean")
    expect_is(values, "GRanges")
    expect_equal(length(gr), 21)
    expect_equal(length(gr), length(values))
    expect_equal(length(values), length(unique(gr$name)))
})


test_that("bw_loci returns correct per locus values on multiple files", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  values <- bw_loci(c(bw1, bw2),
                    get_testfile("labeled.bed"),
                    labels = c("bw1", "bw2"),
                    per_locus_stat = "mean"
  )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[2]$bw1, 4.5)
  expect_equal(values[1]$bw2, 19)
  expect_equal(values[2]$bw2, 16.5)
})


test_that("bw_loci returns correct per locus values on multiple files, with bg",
{
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  bg1 <- local_create_sample_bigwig(get_testfile("bg1.bed"), chromsizes)
  bg2 <- local_create_sample_bigwig(get_testfile("bg2.bed"), chromsizes)

  values <- bw_loci(
      c(bw1, bw2),
      bg_bwfiles = c(bg1, bg2),
      get_testfile("labeled.bed"),
      labels = c("bw1", "bw2"),
      per_locus_stat = "mean"
  )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw2, 9.5)
  expect_equal(values[2]$bw2, 8.25)
})

test_that("bw_loci handles default names with special characters", {
  bw_special <- tempfile("bigwig-2ñ", fileext = ".bw")
  bw <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes, f = bw_special)
  values <- bw_loci(bw,
                    get_testfile("labeled.bed"),
                    aggregate_by = "true_mean"
  )
  expect_is(values, "data.frame")
})

test_that("bw_loci crashes on wrong number of labels for multiple files", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  expect_error({
    values <- bw_loci(c(bw1, bw2), get_testfile("labeled.bed"),
                      labels = "bw1",
                      per_locus_stat = "mean"
    )
  },
  "BigWig file list and column names must have the same length.")
})

test_that("bw_loci returns correct mean-of-means aggregated values", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_loci(bw1,
                    get_testfile("labeled.bed"),
                    labels = "bw1",
                    per_locus_stat = "mean",
                    aggregate_by = "mean"
  )

  expect_is(values, "data.frame")
  expect_equal(values["typeA", "bw1"], 7)
  expect_equal(values["typeB", "bw1"], 13.3333333333)
})

test_that("bw_loci returns correct true_mean aggregated values", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_loci(bw1, get_testfile("labeled.bed"),
                    labels = "bw1",
                    per_locus_stat = "mean",
                    aggregate_by = "true_mean"
  )

  expect_is(values, "data.frame")
  expect_equal(values["typeA", "bw1"], 7)
  # I keep forgetting why this has to be 11.125, so I write a note here:
  # this is true mean, and the length of type b intervals affects this.
  # It is not due to zero vs one-based systems, as BED files are nice like that
  # they have all the good properties length = (end - start) and the intervals
  # are square on the intervals except the ones that explicitly overlap.
  expect_equal(values["typeB", "bw1"], 11.125)
})

test_that("bw_loci returns correct true_mean aggregated values shuffled", {
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    values <- bw_loci(bw1, get_testfile("labeled_shuffled.bed"),
                      labels = "bw1",
                      per_locus_stat = "mean",
                      aggregate_by = "true_mean"
    )
    expect_is(values, "data.frame")
    expect_equal(values["typeA", "bw1"], 7)
    expect_equal(values["typeB", "bw1"], 11.125)
})

test_that("bw_loci returns correct true_mean aggregated with remove_top != 0", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_loci(bw1, get_testfile("labeled.bed"),
                    labels = "bw1",
                    per_locus_stat = "mean",
                    aggregate_by = "true_mean",
                    remove_top = 0.20
  )

  expect_is(values, "data.frame")
  # Does not change bc top 20% out of 5 values is 19 from type B
  expect_equal(values["typeA", "bw1"], 7)
  expect_equal(values["typeB", "bw1"], 8.5)
})

test_that("bw_loci aggregated returns 1 when fg == bg", {
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    values <- bw_loci(bw1, get_testfile("labeled.bed"), bg_bwfiles = bw1,
                      labels = "bw1",
                      per_locus_stat = "mean",
                      aggregate_by = "true_mean",
                      norm_mode = "fc"
    )

    expect_is(values, "data.frame")
    expect_equal(values["typeA", "bw1"], 1)
    expect_equal(values["typeB", "bw1"], 1)
})

test_that("bw_loci aggregated returns 1 when fg == bg and remove_top != 0", {
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    values <- bw_loci(bw1, get_testfile("labeled.bed"), bg_bwfiles = bw1,
                      labels = "bw1",
                      per_locus_stat = "mean",
                      aggregate_by = "true_mean",
                      norm_mode = "fc",
                      remove_top = 0.40
    )

    expect_is(values, "data.frame")
    expect_equal(values["typeA", "bw1"], 1)
    expect_equal(values["typeB", "bw1"], 1)
})


test_that("bw_loci on an empty list throws an error", {
  expect_error({
    values <- bw_loci(c(), get_testfile("labeled.bed"),
                      per_locus_stat = "mean",
                      aggregate_by = "true_mean"
    )
  },
  "File list provided is empty.")

})

test_that("bw_loci errors on non existing files on bwlist", {
  expect_error({
    values <- bw_loci(c("invalidname.bw"),
                      get_testfile("labeled.bed"),
                      per_locus_stat = "mean",
                      aggregate_by = "true_mean"
    )
  },
  "Files not found: invalidname.bw")
})

test_that("bw_loci returns correct median-of-means aggregated values", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_loci(bw1, get_testfile("labeled.bed"),
                    labels = "bw1",
                    per_locus_stat = "mean",
                    aggregate_by = "median"
  )

  expect_is(values, "data.frame")
  expect_equal(values["typeA", "bw1"], 7)
  expect_equal(values["typeB", "bw1"], 16.5)
})

test_that("bw_loci throws error on not implemented aggregate_by", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  expect_error(
    values <- bw_loci(bw1, get_testfile("labeled.bed"),
                      labels = "bw1",
                      per_locus_stat = "mean",
                      aggregate_by = "max"
    )
  )
})

test_that("bw_loci runs with background and aggregate_by parameter", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  values <- bw_loci(bw1, get_testfile("labeled.bed"), bg_bwfiles = bw2,
                    labels = "bw1",
                    per_locus_stat = "mean",
                    aggregate_by = "mean"
  )
  expect_is(values, "data.frame")
})

test_that("bw_loci returns correct values with background == 0 and aggregate_by", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bg_zeros <- local_create_sample_bigwig(get_testfile("bg3.bed"), chromsizes)
  values <- bw_loci(bw1, get_testfile("labeled.bed"), bg_bwfiles = bg_zeros,
                    labels = "bw1",
                    per_locus_stat = "mean",
                    aggregate_by = "mean"
  )
  expect_is(values, "data.frame")
  expect_equal(values["typeA", "bw1"], Inf)
  expect_equal(values["typeB", "bw1"], Inf)
})

test_that("bw_loci returns correct values with 0/0 and aggregate_by", {
  bw_zeros <- local_create_sample_bigwig(get_testfile("bed3.bed"), chromsizes)
  bg_zeros <- local_create_sample_bigwig(get_testfile("bg3.bed"), chromsizes)
  values <- bw_loci(bw_zeros, get_testfile("labeled.bed"), bg_bwfiles = bg_zeros,
                    labels = "bw_zeros",
                    per_locus_stat = "mean",
                    aggregate_by = "mean"
  )

  expect_equal(values["typeA", "bw_zeros"], NaN)
})

test_that("bw_loci excludes NA values from true_mean aggregation", {
  bw_na <- local_create_sample_bigwig(get_testfile("bed_na.bed"), chromsizes)
  values <- suppressWarnings(bw_loci(bw_na, get_testfile("labeled.bed"),
                    labels = "bw_na",
                    per_locus_stat = "mean",
                    aggregate_by = "true_mean"
  ))
  expect_equal(values["typeA", "bw_na"], 12)
})

test_that("bw_loci excludes NA values from mean aggregation", {
  bw_na <- local_create_sample_bigwig(get_testfile("bed_na.bed"), chromsizes)
  values <- suppressWarnings(bw_loci(bw_na, get_testfile("labeled.bed"),
                                     labels = "bw_na",
                                     per_locus_stat = "mean",
                                     aggregate_by = "mean"
  ))
  expect_equal(values["typeA", "bw_na"], 12)
})

test_that("bw_loci excludes NA values from single locus", {
    # Here, exclude means that NAs are not counted as zeros, but overlapping
    # length that has NA value is not counted as valid length.
    # In this test, a locus 21-50 that has NA on 21-40 and 3 on 41-50 yields
    # a mean of 3.
    bw_na <- local_create_sample_bigwig(get_testfile("bed_na.bed"), chromsizes)
    values <- suppressWarnings(bw_loci(bw_na, get_testfile("labeled_na.bed"),
                                       labels = "bw_na",
                                       per_locus_stat = "mean")
    )
    expect_equal(values[1]$bw_na, 3)
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
  bw_na <- local_create_sample_bigwig(get_testfile("bed_na.bed"), chromsizes)
  values <- suppressWarnings(
      bw_loci(bw_na,
              get_testfile("labeled_na.bed"),
              labels = "bw_na",
              per_locus_stat = "mean",
              default_na = 0
      )
  )
  expect_equal(values[1]$bw_na, 3)
})

test_that("bw_loci fails if aggregate_by in an unnamed bed file", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  expect_error({
    values <- bw_loci(bw1, get_testfile("not_labeled.bed"), bg_bwfiles = bw2,
                      labels = "bw1",
                      per_locus_stat = "mean",
                      aggregate_by = "mean"
    )},
    "missing values in 'row.names' are not allowed"
  )
})

## bw_bins ---------------------------------------------------

test_that("bw_bins returns correct per locus values", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_bins(bw1,
                    selection = import(get_testfile("labeled.bed")),
                    labels = "bw1",
                    per_locus_stat = "mean"
  )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 5.5)
  expect_equal(values[2]$bw1, 15.5)
})

test_that("bw_bins returns 1 when bwfile == bg_bwfile", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_bins(bw1,
                    bg_bwfiles = bw1,
                    selection = import(get_testfile("labeled.bed")),
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
                         get_testfile("labeled.bed"),
                         labels = NULL
    )
  },
  "File list provided is empty.")
})

test_that("bw_loci on non-existing bed file throws an error", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
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
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  expect_error({
    values <- bw_profile(bw1,
                         loci = "invalidname.bed",
                         labels = NULL
    )
  },
  "Files not found: invalidname.bed")
})

test_that("bw_profile errors on non existing files on bwlist", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  expect_error({
    values <- bw_profile(c(bw1, "invalidname.bw"),
                         loci = get_testfile("labeled.bed"),
                         labels = NULL
    )
  },
  "Files not found: invalidname.bw")
})

test_that("bw_profile runs quiet on valid parameters", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         loci = get_testfile("labeled.bed"),
                         upstream = 1,
                         downstream = 1,
                         bin_size = 1
    )
  })
})

test_that("bw_profile runs on GRanges object", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)

  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         loci = import(get_testfile("labeled.bed")),
                         upstream = 1,
                         downstream = 1,
                         bin_size = 1
    )
  })
})

test_that("bw_profile runs quiet on valid parameters, mode start", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         loci = get_testfile("labeled.bed"),
                         upstream = 1,
                         downstream = 1,
                         bin_size = 1,
                         mode = "start"
    )
  })
})

test_that(
  "bw_profile runs quiet on valid parameters, mode start, with background", {
    bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
    bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
    bg1 <- local_create_sample_bigwig(get_testfile("bg1.bed"), chromsizes)
    bg2 <- local_create_sample_bigwig(get_testfile("bg2.bed"), chromsizes)
    expect_silent({
      values <- bw_profile(c(bw1, bw2),
                           bg_bwfiles = c(bg1, bg2),
                           loci = get_testfile("labeled.bed"),
                           upstream = 1,
                           downstream = 1,
                           bin_size = 1,
                           mode = "start"
      )
    })
})

test_that("bw_profile normalized returns 1 when fg == bg", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  values <- bw_profile(c(bw1, bw2),
                       bg_bwfiles = c(bw1, bw2),
                       loci = get_testfile("labeled.bed"),
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
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                         loci = get_testfile("labeled.bed"),
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
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         loci = get_testfile("labeled.bed"),
                         upstream = 1,
                         downstream = 1,
                         bin_size = 1,
                         mode = "end"
    )
  })
})

test_that("bw_profile runs quiet on valid parameters, middle parameter", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         loci = get_testfile("labeled.bed"),
                         upstream = 1,
                         downstream = 1,
                         middle = 1,
                         bin_size = 1,
                         mode = "end"
    )
  })
})

test_that("bw_profile runs quiet on valid parameters with background", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  bg1 <- local_create_sample_bigwig(get_testfile("bg1.bed"), chromsizes)
  bg2 <- local_create_sample_bigwig(get_testfile("bg2.bed"), chromsizes)
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                         bg_bwfiles = c(bg1, bg2),
                         loci = get_testfile("labeled.bed"),
                         upstream = 1,
                         downstream = 1,
                         bin_size = 1
    )
  })

})

test_that("bw_profile throws error on flanking region smaller than bin size", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                         loci = get_testfile("labeled.bed"),
                         upstream = 1,
                         downstream = 1,
                         bin_size = 10
    )
  },
  "bin size must be smaller than flanking regions")

})

test_that("bw_profile throws error on negative bin size", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                         loci = get_testfile("labeled.bed"),
                         upstream = 1,
                         downstream = 1,
                         bin_size = -10
    )
  },
  "bin size must be a positive value: -10")
})

test_that("bw_profile throws error on negative upstream value", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                         loci = get_testfile("labeled.bed"),
                         upstream = -10,
                         downstream = 10,
                         bin_size = 10
    )
  },
  "upstream size must be a positive value: -10")
})

test_that("bw_profile throws error on negative downstream value", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  bw2 <- local_create_sample_bigwig(get_testfile("bed2.bed"), chromsizes)
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                         loci = get_testfile("labeled.bed"),
                         upstream = 10,
                         downstream = -10,
                         bin_size = 10
    )
  },
  "downstream size must be a positive value: -10")
})

## bw_heatmap -------------------------------------------

test_that("bw_heatmap returns correct values odd bin size", {
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_heatmap(bw1,
                       bg_bwfiles = NULL,
                       loci = get_testfile("labeled.bed"),
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
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_heatmap(bw1,
                       bg_bwfiles = NULL,
                       loci = get_testfile("labeled.bed"),
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
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_heatmap(bw1,
                       bg_bwfiles = NULL,
                       loci = get_testfile("labeled.bed"),
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
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_heatmap(bw1,
                       bg_bwfiles = NULL,
                       loci = import(get_testfile("labeled.bed")),
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
  bw1 <- local_create_sample_bigwig(get_testfile("bed1.bed"), chromsizes)
  values <- bw_heatmap(bw1,
                       bg_bwfiles = bw1,
                       loci = get_testfile("labeled.bed"),
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
