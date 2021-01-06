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

reduced_bg <-  bw_bins(bg_bw, selection = bw_limits, labels = "x_bg")
reduced_bg_2 <-  bw_bins(bg_bw, selection = bw_limits, labels = "x_bg2")

reduced_bins_all <- bw_bins(c(bw1, bw2), bg_bwfiles=c(bg_bw, bg_bw),
                            selection = bw_limits, labels = c("x", "y"))

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

test_that("plot_bw_loci_scatter no verbose returns a plot with no caption", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_loci_scatter(bw1, bw2, loci = bed, verbose = FALSE)
    expect_is(p, "ggplot")
    expect_false("caption" %in% names(p$labels))
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


test_that("plot_bw_bins_violin with defaults returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_violin(c(bw1, bw2))
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_violin with bg files returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_violin(c(bw1, bw2), bg_bwfiles=c(bw1, bw2))
    expect_is(p, "ggplot")
  })
})


test_that("plot_bw_bins_violin verbose returns a plot with a caption", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_violin(c(bw1, bw2), verbose = TRUE)
    expect_is(p, "ggplot")
    expect_true("caption" %in% names(p$labels))
  })
})

test_that("plot_bw_bins_violin not verbose returns a plot with no caption", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_violin(c(bw1, bw2), verbose = FALSE)
    expect_is(p, "ggplot")
    expect_false("caption" %in% names(p$labels))
  })
})

test_that("plot_bw_bins_violin with highlight returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_violin(c(bw1, bw2), verbose = FALSE, highlight = bed)
    expect_is(p, "ggplot")
    expect_false("caption" %in% names(p$labels))
  })
})

test_that("plot_bw_bins_violin with highlight GRanges returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_violin(c(bw1, bw2), verbose = FALSE, highlight = import(bed), highlight_label="Bedfile")
    expect_is(p, "ggplot")
    expect_false("caption" %in% names(p$labels))
  })
})

test_that("plot_bw_bins_violin with highlight and remove top returns a plot", {
  m <- mock(reduced_bins_all)
  with_mock(bw_bins = m, {
    p <- p <- plot_bw_bins_violin(c(bw1, bw2),
                                  bg_bwfiles = c(bg_bw, bg_bw),
                                  labels = c("A", "B"),
                                  highlight = bed,
                                  bin_size = 5000,
                                  norm_func = log2,
                                  genome = "hg38",
                                  remove_top = 0.01,
                                  verbose = FALSE)

    expect_is(p, "ggplot")
    expect_false("caption" %in% names(p$labels))
  })
})


test_that("plot_bw_bins_violin passes on parameters", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_violin(c(bw1, bw2),
                             bg_bwfiles = c(bg_bw, bg_bw),
                             labels = c("A", "B"),
                             highlight = bed,
                             bin_size = 5000,
                             norm_func = log2,
                             genome = "hg38",
                             remove_top = 0
    )
  })

  expect_call(m, 1,
              bw_bins(bwfiles,
                      bg_bwfiles = bg_bwfiles,
                      labels = labels,
                      bin_size = bin_size,
                      genome = genome,
                      per_locus_stat = per_locus_stat,
                      norm_func = norm_func,
                      # FIXME: Remove top is done outside this function
                      remove_top = 0
              )
  )

  expect_args(m, 1,
              bwfiles = c(bw1, bw2),
              bg_bwfiles = c(bg_bw, bg_bw),
              labels = c("A", "B"),
              bin_size = 5000,
              genome = "hg38",
              per_locus_stat = "mean",
              norm_func = log2,
              remove_top = 0
  )

})



test_that(
  "plot_bw_loci_summary_heatmap with defaults returns a ggplot object", {
    m <- mock(summary_values)
    with_mock(bw_bed = m, {
      p <- plot_bw_loci_summary_heatmap(c(bw1, bw2), loci = bed)
      expect_is(p, "ggplot")
    })
  })

test_that("plot_bw_loci_summary_heatmap with verbose set returns a plot with a caption", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_loci_summary_heatmap(c(bw1, bw2), loci = bed, verbose = TRUE)
    expect_is(p, "ggplot")
    expect_true("caption" %in% names(p$labels))
  })
})

test_that("plot_bw_loci_summary_heatmap with verbose unset returns a plot without a caption", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_loci_summary_heatmap(c(bw1, bw2), loci = bed, verbose = FALSE)
    expect_is(p, "ggplot")
    expect_false("caption" %in% names(p$labels))
  })
})



test_that(
  "plot_bw_loci_summary_heatmap passes parameters on", {
    m <- mock(summary_values)
    with_mock(bw_bed = m, {
      p <- plot_bw_loci_summary_heatmap(c(bw1, bw2),
                                       loci = bed,
                                       bg_bwfiles <- c(bg_bw, bg_bw),
                                       aggregate_by = "median",
                                       norm_func = log2,
                                       labels = c("bw1", "bw2")
      )
    })

    expect_call(m, 1,
                bw_bed(bwfiles,
                       loci,
                       bg_bwfiles = bg_bwfiles,
                       aggregate_by = aggregate_by,
                       norm_func = norm_func,
                       labels = labels,
                       remove_top = remove_top
                )
    )

    expect_args(m, 1,
                bwfiles = c(bw1, bw2),
                loci = bed,
                bg_bwfiles = c(bg_bw, bg_bw),
                aggregate_by = "median",
                norm_func = log2,
                labels = c("bw1", "bw2"),
                remove_top = 0
    )
  })

