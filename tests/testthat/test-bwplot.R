context("Test functions for bigWig plots")
library(GenomicRanges)
library(testthat)
library(mockery)
library(withr)

# Setup -------------------------------------------------------------------
chromsizes <- c(200, 200)

# Bins scatter tests ----------------------------------------------

test_that("plot_bw_bins_scatter with defaults returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins(), {
    p <- plot_bw_bins_scatter(bw1, bw2)
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_scatter with highlight set returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins(), {
    p <- plot_bw_bins_scatter(bw1, bw2, highlight = get_datafile("sample_genes_mm9.bed"))
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_scatter with verbose set returns a plot with a caption", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins(), {
    p <- plot_bw_bins_scatter(bw1, bw2, highlight = get_datafile("sample_genes_mm9.bed"), verbose = TRUE)
    expect_is(p, "ggplot")
    expect_false(is.null(p$labels$caption))
  })
})

test_that("plot_bw_bins_scatter with highlight colors set returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins(), {
    p <- plot_bw_bins_scatter(
        bw1, bw2,
        highlight = get_datafile("sample_genes_mm9.bed"),
        highlight_color = "green"
    )
    expect_is(p, "ggplot")

  })
})

test_that("plot_bw_bins_scatter with GRanges and label returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  gr_loci <- rtracklayer::import(get_datafile("sample_genes_mm9.bed"))
  with_mock(bw_bins = make_mock_bins(), {
    p <- plot_bw_bins_scatter(
        bw1, bw2, highlight = gr_loci,
        highlight_label = "A_label"
    )
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_scatter with several GRanges and labels returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  gr_loci <- rtracklayer::import(get_datafile("sample_genes_mm9.bed"))
  with_mock(bw_bins = make_mock_bins(), {
    p <- plot_bw_bins_scatter(bw1, bw2,
           highlight = list(gr_loci, gr_loci),
           highlight_label = c("A_label", "Another label"))
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_scatter crashes with unmatched labels/highlight", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  gr_loci <- rtracklayer::import(get_datafile("sample_genes_mm9.bed"))
  with_mock(bw_bins = make_mock_bins(), {
    expect_error(
        p <- plot_bw_bins_scatter(
            bw1, bw2,
            highlight = c(gr_loci),
            highlight_label = c("A_label", "Another label")),
        "Highlight loci sets don't match the number of labels provided"
    )
  })
})

test_that("plot_bw_bins_scatter with GRanges and no label crashes", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  gr_loci <- rtracklayer::import(get_datafile("sample_genes_mm9.bed"))
  expect_error(
    with_mock(bw_bins = make_mock_bins(), {
      p <- plot_bw_bins_scatter(bw1, bw2, highlight = gr_loci)
    }),
    "GRanges used as highlight loci but no labels provided")
})

test_that("plot_bw_bins_scatter with density prints a deprecation message", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins(), {
    expect_message(
      p <- plot_bw_bins_scatter(bw1, bw2, density = TRUE),
      paste("plot_bw_bins_scatter with density = TRUE is deprecated.",
            "Please use plot_bw_bins_density instead."
      )
    )
  })
})

# Bins density tests ---------------------------------------------

test_that("plot_bw_bins_density with defaults returns a plot with tile layer", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins(), {
    p <- plot_bw_bins_density(bw1, bw2)
    expect_is(p, "ggplot")
    expect_true("GeomTile" %in% sapply(p$layers, function(x) class(x$geom)[1]))
  })
})

test_that("plot_bw_bins_density with verbose set to false returns a plot with no caption", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  bed <- local_file("bed.bed")
  with_mock(bw_bins = make_mock_bins(), {
    p <- plot_bw_bins_density(bw1, bw2, verbose = FALSE)
    expect_is(p, "ggplot")
    expect_true(is.null(p$labels$caption))
    expect_true("GeomTile" %in% sapply(p$layers, function(x) class(x$geom)[1]))
  })
})

# Loci scatter tests ----------------------------------------------

test_that("plot_bw_loci_scatter with defaults returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  bed <- local_file("bed.bed")
  with_mock(bw_loci = make_mock_bins(), {
    p <- plot_bw_loci_scatter(bw1, bw2, bed)
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_loci_scatter with highlight set returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  bed <- local_file("bed.bed")
  with_mock(bw_loci = make_mock_bins(), {
    p <- plot_bw_loci_scatter(
      bw1, bw2,
      loci = bed,
      highlight = get_datafile("sample_genes_mm9.bed")
    )
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_loci_scatter with verbose set returns a plot with a caption", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  bed <- local_file("bed.bed")
  with_mock(bw_loci = make_mock_bins(), {
    p <- plot_bw_loci_scatter(bw1, bw2, loci = bed, verbose = TRUE)
    expect_is(p, "ggplot")
    expect_false(is.null(p$labels$caption))
  })
})

test_that("plot_bw_loci_scatter no verbose returns a plot with no caption", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  bed <- local_file("bed.bed")
  with_mock(bw_loci = make_mock_bins(), {
    p <- plot_bw_loci_scatter(bw1, bw2, loci = bed, verbose = FALSE)
    expect_is(p, "ggplot")
    expect_true(is.null(p$labels$caption))
  })
})

test_that("plot_bw_loci_scatter with remove_top returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  bed <- local_file("bed.bed")
  with_mock(bw_loci = make_mock_bins(), {
    p <- plot_bw_loci_scatter(bw1, bw2, loci = bed, verbose = TRUE, remove_top = 0.01)
    expect_is(p, "ggplot")
    expect_false(is.null(p$labels$caption))
  })
})

# Bins violin tests ----------------------------------------------

test_that("plot_bw_bins_violin with defaults returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins_all(), {
    p <- plot_bw_bins_violin(c(bw1, bw2), labels = c("x", "y"))
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_violin with bg files returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins_all(), {
    p <- plot_bw_bins_violin(c(bw1, bw2), bg_bwfiles = c(bw1, bw2), labels = c("x", "y"))
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_violin verbose returns a plot with a caption", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins_all(), {
    p <- plot_bw_bins_violin(c(bw1, bw2), verbose = TRUE, labels=c("x", "y"))
    expect_is(p, "ggplot")
    expect_false(is.null(p$labels$caption))
  })
})

test_that("plot_bw_bins_violin not verbose returns a plot with no caption", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins_all(), {
    p <- plot_bw_bins_violin(c(bw1, bw2), verbose = FALSE, labels=c("x", "y"))
    expect_is(p, "ggplot")
    expect_true(is.null(p$labels$caption))
  })
})

test_that("plot_bw_bins_violin with highlight returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins_all(), {
    p <- plot_bw_bins_violin(c(bw1, bw2), verbose = FALSE, highlight = get_datafile("sample_genes_mm9.bed"), labels=c("x", "y"))
    expect_is(p, "ggplot")
    expect_true(is.null(p$labels$caption))

  })
})

test_that("plot_bw_bins_violin with highlight GRanges returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  gr_loci <- rtracklayer::import(get_datafile("sample_genes_mm9.bed"))
  with_mock(bw_bins = make_mock_bins_all(), {
    p <- plot_bw_bins_violin(
        c(bw1, bw2),
        verbose = FALSE,
        highlight = gr_loci,
        highlight_label = "loci",
        labels = c("x", "y")
    )
    expect_is(p, "ggplot")
    expect_true(is.null(p$labels$caption))
  })
})

test_that("plot_bw_bins_violin with highlight and remove top returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  with_mock(bw_bins = make_mock_bins_all(), {
    p <- plot_bw_bins_violin(
        c(bw1, bw2),
        bg_bwfiles = c(bw1, bw1),
        labels = c("x", "y"),
        highlight = get_datafile("sample_genes_mm9.bed"),
        bin_size = 5000,
        norm_mode = "log2fc",
        genome = "hg38",
        remove_top = 0.01,
        verbose = FALSE
    )

    expect_is(p, "ggplot")
    expect_true(is.null(p$labels$caption))
  })
})

# Summary heatmap tests -------------------------------------------

test_that("plot_bw_loci_summary_heatmap with defaults returns a ggplot object", {
    bw1 <- local_file("bw1.bw")
    bw2 <- local_file("bw2.bw")
    bed <- local_file("loci.bed")
    with_mock(bw_loci = make_mock_summary(), {
      p <- plot_bw_loci_summary_heatmap(c(bw1, bw2), loci = bed)
      expect_is(p, "ggplot")
    })
  })

test_that("plot_bw_loci_summary_heatmap with labels returns a ggplot object", {
    bw1 <- local_file("bw1.bw")
    bw2 <- local_file("bw2.bw")
    bed <- local_file("loci.bed")
    with_mock(bw_loci = make_mock_summary(), {
      p <- plot_bw_loci_summary_heatmap(c(bw1, bw2), loci = bed, labels = c("A", "B"))
      expect_is(p, "ggplot")
    })
  })

test_that("plot_bw_loci_summary_heatmap with labels including invalid chars returns a ggplot object", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  bed <- local_file("loci.bed")
  with_mock(bw_loci = make_mock_summary(), {
    p <- plot_bw_loci_summary_heatmap(c(bw1, bw2), loci = bed, labels = c("A-1", "B-2"))
    expect_is(p, "ggplot")
  })
  })

test_that("plot_bw_loci_summary_heatmap with verbose set returns a plot with a caption", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  bed <- local_file("loci.bed")
  with_mock(bw_loci = make_mock_summary(), {
    p <- plot_bw_loci_summary_heatmap(c(bw1, bw2), loci = bed, verbose = TRUE)
    expect_is(p, "ggplot")
    expect_true("caption" %in% names(p$labels))
  })
})

test_that("plot_bw_loci_summary_heatmap with verbose unset returns a plot without a caption", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  bed <- local_file("loci.bed")
  with_mock(bw_loci = make_mock_summary(), {
    p <- plot_bw_loci_summary_heatmap(c(bw1, bw2), loci = bed, verbose = FALSE)
    expect_is(p, "ggplot")
    expect_true(is.null(p$labels$caption))
  })
})

# Profile tests ----------------------------------------------

test_that("plot_bw_profile with defaults returns a ggplot object", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    p <- plot_bw_profile(c(bw1, bw2), loci = bed)
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_profile with GRanges returns a ggplot object", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    p <- plot_bw_profile(c(bw1, bw2), loci = import(bed))
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_profile with GRanges list returns a ggplot object", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    p <- plot_bw_profile(bw1, loci = list(import(bed), import(bed)), labels = c("A", "B"))
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_profile with GRanges list and no labels throws warning", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    expect_warning(p <- plot_bw_profile(bw1, loci = list(import(bed), import(bed))),
                   "Unlabeled objects or repeated labels. Adding numeric indices.")

  })
})

test_that("plot_bw_profile with loci file list returns a ggplot object", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    p <- plot_bw_profile(bw1, loci = c(bed, bed), labels = c("A", "B"))
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_profile with loci file list and bwlist fails", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    expect_error(p <- plot_bw_profile(c(bw1, bw2), loci = c(bed, bed)),
                 "If multiple loci provided only a single bwfile is allowed")
  })
})

test_that("plot_bw_profile with defaults has a line layer", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    p <- plot_bw_profile(c(bw1, bw2), loci = bed)
    expect_is(p, "ggplot")
    expect_true("GeomLine" %in% sapply(p$layers, function(x) class(x$geom)[1]))
  })
})

test_that("plot_bw_profile with show_error has a ribbon layer", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    p <- plot_bw_profile(c(bw1, bw2), loci = bed, show_error=TRUE)
    expect_is(p, "ggplot")
    expect_true("GeomRibbon" %in% sapply(p$layers, function(x) class(x$geom)[1]))
  })
})

test_that("plot_bw_profile with show_error and background does not have ribbon layer", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    expect_warning({p <- plot_bw_profile(bw1, bg_bwfiles = bw2, loci = bed, show_error = TRUE)},
                   "Stderr estimate not available when normalizing by input")

    expect_is(p, "ggplot")
    expect_false("GeomRibbon" %in% sapply(p$layers, function(x) class(x$geom)[1]))
  })
})

test_that("plot_bw_profile verbose returns a plot with a caption", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    p <- plot_bw_profile(c(bw1, bw2), loci = bed, verbose = TRUE)
    expect_is(p, "ggplot")
    expect_false(is.null(p$labels$caption))
  })
})

test_that("plot_bw_profile mode start valid parameters returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    p <- plot_bw_profile(c(bw1, bw2), loci = bed, mode = "start")
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_profile mode end valid parameters returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    p <- plot_bw_profile(c(bw1, bw2), loci = bed, mode = "end")
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_profile mode center valid parameters returns a plot", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    p <- plot_bw_profile(c(bw1, bw2), loci = bed, mode = "center")
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_profile not verbose returns a plot without a caption", {
  bw1 <- local_file("bw1.bw")
  bw2 <- local_file("bw2.bw")
  # bw_profile also calls a function to count the number of loci, so it is not
  # all covered by mock
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_profile = make_mock_profile(), {
    p <- plot_bw_profile(c(bw1, bw2), loci = bed, verbose = FALSE)
    expect_is(p, "ggplot")
    expect_true(is.null(p$labels$caption))
  })
})

# Heatmap tests ----------------------------------------------

test_that("plot_bw_heatmap with defaults returns a ggplot object", {
  bw1 <- local_file("bw1.bw")
  bed <- local_file("bed.bed")
  with_mock(bw_heatmap = mock(get_heatmap_values(), cycle = TRUE), {
    p <- plot_bw_heatmap(bw1, loci = bed)
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_heatmap with GRanges returns a ggplot object", {
  bw1 <- local_file("bw1.bw")
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_heatmap = mock(get_heatmap_values(), cycle = TRUE), {
    p <- plot_bw_heatmap(bw1, loci = import(bed))
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_heatmap sorts decreasingly by mean", {
  bw1 <- local_file("bw1.bw")
  bed <- get_datafile("sample_genes_mm9.bed")
  # y axis works reversed here
  heatmap_values <- get_heatmap_values()
  other_order <- order(rowMeans(heatmap_values[[1]]), decreasing = F)
  with_mock(bw_heatmap = mock(heatmap_values, cycle = TRUE), {
    p <- plot_bw_heatmap(bw1, loci = bed)
    p2 <- plot_bw_heatmap(bw1, loci = bed, order_by = other_order)
    expect_is(p, "ggplot")
    expect_identical(p$data, p2$data)
  })
})

test_that("plot_bw_heatmap with order returns a ggplot object", {
  bw1 <- local_file("bw1.bw")
  bed <- get_datafile("sample_genes_mm9.bed")
  heatmap_values <- get_heatmap_values()
  other_order <- order(rowMeans(heatmap_values[[1]]), decreasing = T)
  with_mock(bw_heatmap = mock(heatmap_values, cycle = TRUE), {
    p <- plot_bw_heatmap(bw1, loci = bed, order_by = other_order)
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_heatmap with different order returns a different ggplot", {
  bw1 <- local_file("bw1.bw")
  bed <- get_datafile("sample_genes_mm9.bed")
  heatmap_values <- get_heatmap_values()
  other_order <- order(rowMeans(heatmap_values[[1]]), decreasing = T)
  with_mock(bw_heatmap = mock(heatmap_values, cycle = TRUE), {
    p <- plot_bw_heatmap(bw1, loci = bed, order_by = other_order)
    p2 <- plot_bw_heatmap(bw1, loci = bed)
    expect_is(p, "ggplot")
    expect_false(p$data[1, "value"] == p2$data[1, "value"])
  })
})

test_that("plot_bw_heatmap with verbose returns a ggplot object with a caption", {
  bw1 <- local_file("bw1.bw")
  bed <- get_datafile("sample_genes_mm9.bed")
  with_mock(bw_heatmap = mock(get_heatmap_values(), cycle = TRUE), {
    p <- plot_bw_heatmap(bw1, loci = bed)
    expect_is(p, "ggplot")
    expect_true("caption" %in% names(p$labels))
  })
})
