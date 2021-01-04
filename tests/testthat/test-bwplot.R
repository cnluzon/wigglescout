context("Test functions for bigWig plots")
library(GenomicRanges)
library(testthat)
library(mockery)


get_file_path <- function(filename) {
  system.file("extdata", filename, package = "wigglescout")
}

bw1 <- get_file_path("sample_H33_ChIP.bw")
bw2 <- get_file_path("sample_H3K9me3_ChIP.bw")
bg_bw <- get_file_path("sample_Input.bw")
bed <- get_file_path("sample_genes_mm9.bed")
bed_summary <- get_file_path("sample_chromhmm.bed")

bw_limits <- GRanges(seqnames = c("chr15"),
                     ranges = IRanges(c(102723600, 102959000)))

reduced_bins <- bw_bins(bw1, selection = bw_limits, labels = "x")
reduced_bins_2 <- bw_bins(bw2, selection = bw_limits, labels = "y")
summary_values <- bw_bed(c(bw1, bw2), bed_summary, aggregate_by = "mean")
profile_values <- bw_profile(bw1,
                             bedfile = bed,
                             upstream = 1000,
                             downstream = 1000
)

heatmap_values <- bw_heatmap(bw1,
                             bedfile = bed,
                             upstream = 1000,
                             downstream = 1000)

test_that("Setup files exist", {
  expect_true(file_test("-f", bw1))
  expect_true(file_test("-f", bw2))
  expect_true(file_test("-f", bg_bw))
  expect_true(file_test("-f", bed))
})


test_that("plot_bw_bins_scatter with defaults returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_scatter(bw1, bw2)
    expect_is(p, "ggplot")
  })
})


test_that("plot_bw_loci_scatter with defaults returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_loci_scatter(bw1, bw2, bed)
    expect_is(p, "ggplot")
  })
})


test_that("plot_bw_bins_scatter with highlight set returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_scatter(bw1, bw2, highlight = bed)
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_loci_scatter with highlight set returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_loci_scatter(bw1, bw2, bed, highlight = bed)
    expect_is(p, "ggplot")
  })
})


test_that("plot_bw_bins_scatter with verbose set returns a plot with a caption", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_scatter(bw1, bw2, highlight = bed, verbose = TRUE)
    expect_is(p, "ggplot")
    expect_true("caption" %in% names(p$labels))
  })
})

test_that("plot_bw_loci_scatter with verbose set returns a plot with a caption", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_loci_scatter(bw1, bw2, loci = bed, verbose = TRUE)
    expect_is(p, "ggplot")
    expect_true("caption" %in% names(p$labels))
  })
})

test_that("plot_bw_loci_scatter with remove_top returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_loci_scatter(bw1, bw2, loci = bed, verbose = TRUE, remove_top=0.01)
    expect_is(p, "ggplot")
    expect_true("caption" %in% names(p$labels))
  })
})


test_that("plot_bw_bins_scatter with highlight colors set returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_scatter(bw1, bw2, highlight = bed,
            highlight_color = "green")
    expect_is(p, "ggplot")

  })
})

test_that("plot_bw_bins_scatter with GRanges and label returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_scatter(bw1, bw2, highlight = rtracklayer::import(bed),
                              highlight_label = "A_label")
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_scatter with several GRanges and labels returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_scatter(bw1, bw2,
           highlight = list(rtracklayer::import(bed), rtracklayer::import(bed)),
           highlight_label = c("A_label", "Another label"))
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_scatter crashes with unmatched labels/highlight", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    expect_error(p <- plot_bw_bins_scatter(bw1, bw2,
                        highlight = c(rtracklayer::import(bed)),
                        highlight_label = c("A_label", "Another label")),
                 "Highlight loci sets don't match the number of labels provided"
    )
  })
})


test_that("plot_bw_bins_scatter with GRanges and no label crashes", {
  m <- mock(reduced_bins, reduced_bins_2)
    expect_error(
    with_mock(bw_bins = m, {
      p <- plot_bw_bins_scatter(bw1, bw2, highlight = rtracklayer::import(bed))

    }), "GRanges used as highlight loci but no labels provided")
})


test_that("plot_bw_bins_scatter with bg files passes on parameters", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_scatter(bw1, bw2, bg_x = bg_bw, bg_y = bg_bw)
  })

  expect_call(m, 1,
              bw_bins(
                x,
                bg_bwfiles = bg_x,
                bin_size = bin_size,
                genome = genome,
                norm_func = norm_func_x,
                labels = "x"
              )
  )

  expect_args(m, 1,
              x = bw1,
              bg_bwfiles = bg_bw,
              bin_size = 10000,
              genome = "mm9",
              norm_func = identity,
              labels = "x"
  )

  expect_args(m, 2,
              x = bw2,
              bg_bwfiles = bg_bw,
              bin_size = 10000,
              genome = "mm9",
              norm_func = identity,
              labels = "y"
  )
})


