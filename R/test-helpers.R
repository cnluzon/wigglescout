#' Parse a scored BED file and exports a bigWig file with given chromsizes
#'
#' @param bed BED filename
#' @param bw bigWig filename
#' @param chromsizes Vector with length of chromosomes
#' @importFrom GenomeInfoDb `seqlengths<-`
#' @importFrom rtracklayer import export
#'
#' @return NULL
bed_to_bw <- function(bed, bw, chromsizes) {
    ranges <- import(bed)
    seqlengths(ranges) <- chromsizes
    export(ranges, bw)
}

#' Create a self-destructing local temporary bigwig based on a BED file and
#' chromsizes array.
#'
#' @param bed BED filename
#' @param f bigWig filename
#' @param env Running environment
#' @param chromsizes Vector with length of chromosomes
#'
#' @return Local filename
local_create_sample_bigwig <- function(bed, chromsizes, f = tempfile(fileext=".bw"), env = parent.frame()) {
    bed_to_bw(bed, f, chromsizes)

    withr::defer(
        unlink(f),
        envir = env
    )
    f
}

#' Return the full path of a test file in wigglescout package
#'
#' @param fname Name of the file (basename)
#'
#' @return Full path of the file
get_testfile <- function(fname) {
    system.file("testdata", fname, package = "wigglescout")
}

#' Return the full path of a test file in wigglescout package
#'
#' @param fname Name of the file (basename)
#'
#' @return Full path of the file
get_datafile <- function(fname) {
    system.file("extdata", fname, package = "wigglescout")
}

#' Make example test tiles
#'
#' @return genomic ranges object
#' @importFrom GenomicRanges tileGenome
make_test_tiles <- function() {
    tileGenome(
        c(chr1 = 200, chr2 = 200),
        tilewidth = 20,
        cut.last.tile.in.chrom = TRUE
    )
}

#' Make a mock object that returns bins GRanges for test plotting
#'
#' @return GRanges object
make_mock_bins <- function() {
    bins_obj_1 <- readRDS(testthat::test_path("fixtures", "mini_h33_bins.rds"))
    bins_obj_2 <- readRDS(testthat::test_path("fixtures", "mini_h3k9_bins.rds"))
    mockery::mock(bins_obj_1, bins_obj_2)
}

#' Make a mock object that returns bins GRanges for test plotting
#'
#' @return GRanges object
#' @importFrom GenomicRanges GRanges
make_mock_bins_all <- function() {
  values <- readRDS(testthat::test_path("fixtures", "mini_bins_all.rds"))
  mockery::mock(values)
}

#' Make a mock object from a fixture that has a profile result from bw_profile
#' for test plotting
#'
make_mock_profile <- function() {
  values <- readRDS(testthat::test_path("fixtures", "mini_profiles.rds"))
  mockery::mock(values, cycle = TRUE)
}

#' Make a mock object from a fixture that has a summary result from bw_heatmap
#' for test plotting
#'
make_mock_summary <- function() {
  values <- readRDS(testthat::test_path("fixtures", "mini_summary.rds"))
  mockery::mock(values)
}

#' Get heatmap fixture from bw_heatmap for test plotting
#'
get_heatmap_values <- function() {
  readRDS(testthat::test_path("fixtures", "mini_heatmap.rds"))
}
