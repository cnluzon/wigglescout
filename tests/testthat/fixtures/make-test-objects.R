library(wigglescout)

make_reduced_bins_fixture <- function(bw, out) {
  bw_limits <- GenomicRanges::GRanges(
    seqnames = c("chr15"),
    ranges = IRanges::IRanges(c(102723600, 102959000))
  )

  reduced_bins <- bw_bins(
    system.file("extdata", bw, package = "wigglescout"),
    selection = bw_limits, labels = "score"
  )
  saveRDS(reduced_bins, testthat::test_path("fixtures", out))
}

make_summary_values_fixture <- function(out) {
  bw1 <- get_datafile("sample_H33_ChIP.bw")
  bw2 <- get_datafile("sample_H3K9me3_ChIP.bw")
  bed_summary <- get_datafile("sample_chromhmm.bed")
  values <- bw_loci(c(bw1, bw2), bed_summary, aggregate_by = "mean")
  saveRDS(values, testthat::test_path("fixtures", out))
}

make_profile_values_fixture <- function(out) {
  bw1 <- get_datafile("sample_H33_ChIP.bw")
  bed <- get_datafile("sample_genes_mm9.bed")
  values <- bw_profile(bw1, loci = bed, upstream = 1000, downstream = 1000)
  saveRDS(values, testthat::test_path("fixtures", out))
}

make_heatmap_values_fixture <- function(out) {
  bw1 <- get_datafile("sample_H33_ChIP.bw")
  bed <- rtracklayer::import(get_datafile("sample_genes_mm9.bed"))[1:6]

  heatmap_values <- bw_heatmap(bw1, loci = bed, upstream = 1000, downstream = 1000)
  saveRDS(heatmap_values, testthat::test_path("fixtures", out))
}

make_bins_multiple_fixture <- function(out) {
  bw1 <- get_datafile("sample_H33_ChIP.bw")
  bw2 <- get_datafile("sample_H3K9me3_ChIP.bw")

  bw_limits <- GenomicRanges::GRanges(
    seqnames = c("chr15"),
    ranges = IRanges::IRanges(c(102723600, 102959000))
  )

  reduced_bins_all <- bw_bins(
    c(bw1, bw2),
    selection = bw_limits,
    labels = c("x", "y")
  )
  saveRDS(reduced_bins_all, testthat::test_path("fixtures", out))
}

make_reduced_bins_fixture("sample_H33_ChIP.bw", "mini_h33_bins.rds")
make_reduced_bins_fixture("sample_H3K9me3_ChIP.bw", "mini_h3k9_bins.rds")
make_reduced_bins_fixture("sample_Input.bw", "mini_input_bins.rds")
make_summary_values_fixture("mini_summary.rds")
make_profile_values_fixture("mini_profiles.rds")
make_heatmap_values_fixture("mini_heatmap.rds")
make_bins_multiple_fixture("mini_bins_all.rds")
