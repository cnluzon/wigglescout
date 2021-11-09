# Copyright (C) 2020 Carmen Navarro Luzón <carmen.navarro@scilifelab.se>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Main functions ---------------------------------------------------

#' Build a bed-scored GRanges object from a bigWig file list and a BED file.
#'
#' Build a scored GRanges object from a bigWig file list and a BED file.
#' The aggregating function (per locus) can be min, max, sd, mean.
#'
#' bwfiles and bg_bwfiles must have the same length. If you are using same
#' background for several files, then file paths must be repeated accordingly.
#'
#' Values can be normalized using background bigWig files (usually input
#' files). By default, the value obtained will be bigwig / bg_bigwig per bin,
#' per bigWig.
#'
#' If norm_func is specified, this can be changed to any given function, for
#' instance, if norm_func = log2, values will represent log2(bigwig / bg_bigwig)
#' per bin.
#'
#' @param loci GRanges or BED file to summarize the BigWig file at.
#' @param aggregate_by Statistic to aggregate per group. If NULL, values are
#'    not aggregated. This is the behavior by default.
#' @export
#' @inheritParams bw_bins
#' @importFrom rtracklayer import BigWigFile
#' @importFrom GenomeInfoDb sortSeqlevels
bw_loci <- function(bwfiles,
                    loci,
                    bg_bwfiles = NULL,
                    labels = NULL,
                    per_locus_stat = "mean",
                    aggregate_by = NULL,
                    norm_mode = "fc",
                    remove_top = 0) {

  validate_filelist(bwfiles)
  validate_locus_parameter(loci)
  norm_func <- .process_norm_mode(norm_mode)

  if (is.null(labels)) {
    labels <- make_label_from_object(bwfiles)
  } else {
    # Ensures later on we only try to access valid labels
    labels <- make.names(labels)
  }

  bed <- loci_to_granges(loci)

  result <- NULL
  if (is.null(aggregate_by)) {
    if (is.null(bg_bwfiles)) {
      result <- .multi_bw_ranges(bwfiles, labels,
        granges = bed,
        per_locus_stat = per_locus_stat,
        remove_top = remove_top
      )
    }
    else {
      # Only want to normalize per-locus if not aggregating
      result <- .multi_bw_ranges_norm(
        bwfiles,
        bg_bwfilelist = bg_bwfiles,
        labels = labels,
        granges = bed,
        per_locus_stat = per_locus_stat,
        norm_func = norm_func,
        remove_top = remove_top,
      )
    }
  } else {
    result <- .multi_bw_ranges_aggregated(bwfiles,
      labels = labels,
      granges = bed,
      per_locus_stat = per_locus_stat,
      aggregate_by = aggregate_by,
      remove_top = remove_top
    )

    if (!is.null(bg_bwfiles)) {
      bg <- .multi_bw_ranges_aggregated(bg_bwfiles,
        labels = labels,
        granges = bed,
        per_locus_stat = per_locus_stat,
        aggregate_by = aggregate_by,
        remove_top = remove_top
      )

      rows <- rownames(result)
      result <- data.frame(norm_func(result[rows, labels] / bg[rows, labels]))
      rownames(result) <- rows
      colnames(result) <- labels
    }
  }

  result
}


#' Build a binned-scored GRanges object from a bigWig file
#'
#' Build a binned-scored GRanges object from a bigWig file. The aggregating
#' function per bin can be min, max, sd, mean.
#'
#' bwfiles and bg_bwfiles must have the same length. If you are using same
#' background for several files, then file paths must be repeated accordingly.
#'
#' Values can be normalized using background bigWig files (usually input
#' files). By default, the value obtained will be bigwig / bg_bigwig per bin,
#' per bigWig.
#'
#' If norm_func is specified, this can be changed to any given function, for
#' instance, if norm_func = log2, values will represent log2(bigwig / bg_bigwig)
#' per bin.
#'
#' @param bwfiles Path or array of paths to the bigWig files to be summarized.
#' @param bg_bwfiles Path or array of paths to the bigWig files to be used as
#'   background.
#' @param labels List of names to give to the mcols of the returned GRanges
#'     object. If NULL, file names are used.
#' @param per_locus_stat Aggregating function (per locus). Mean by default.
#'     Choices: min, max, sd, mean.
#' @param bin_size Bin size.
#' @param genome Genome. Available choices are mm9, hg38.
#' @param selection A GRanges object to restrict binning to a certain set of
#'     intervals. This is useful for debugging and improving performance of
#'     locus specific analyses.
#' @param norm_mode Function to apply to normalize bin values. Default fc:
#'     divides bw / bg. Alternative: log2fc: returns log2(bw/bg).
#' @param remove_top Return range 0-(1-remove_top). By default returns the
#'     whole distribution (remove_top == 0).
#' @return A GRanges object with each bwfile as a metadata column named
#'     after labels, if provided, or after filenames otherwise.
#' @export
bw_bins <- function(bwfiles,
                    bg_bwfiles = NULL,
                    labels = NULL,
                    per_locus_stat = "mean",
                    bin_size = 10000,
                    genome = "mm9",
                    selection = NULL,
                    norm_mode = "fc",
                    remove_top = 0) {

  validate_filelist(bwfiles)
  norm_func <- .process_norm_mode(norm_mode)

  if (is.null(labels)) {
    labels <- make_label_from_object(bwfiles)
  }

  tiles <- build_bins(bin_size = bin_size, genome = genome)

  if (is.null(bg_bwfiles)) {
    result <- .multi_bw_ranges(bwfiles, labels, tiles,
      per_locus_stat = per_locus_stat,
      selection = selection,
      remove_top = remove_top
    )
  } else {
    result <- .multi_bw_ranges_norm(bwfiles, bg_bwfiles, labels, tiles,
      per_locus_stat = per_locus_stat,
      selection = selection,
      norm_func = norm_func,
      remove_top = remove_top
    )
  }

  result
}

#' Calculate heatmap matrix of a bigWig file over a BED file
#'
#' @inheritParams bw_profile
#' @importFrom furrr future_map future_map2
#' @importFrom rtracklayer import
#' @importFrom purrr partial
#' @export
bw_heatmap <- function(bwfiles,
                       bg_bwfiles = NULL,
                       loci = NULL,
                       labels = NULL,
                       mode = "stretch",
                       bin_size = 100,
                       upstream = 2500,
                       downstream = 2500,
                       middle = NULL,
                       ignore_strand = FALSE,
                       norm_mode = "fc") {

  validate_filelist(bwfiles)
  validate_locus_parameter(loci)
  granges <- loci_to_granges(loci)
  norm_func <- .process_norm_mode(norm_mode)

  validate_profile_parameters(bin_size, upstream, downstream)

  if (is.null(labels)) {
    labels <- basename(bwfiles)
  }

  if (length(bwfiles) != length(labels)) {
    stop("labels and bwfiles must have the same length")
  }

  calculate_matrix_norm_fixed <- partial(.calculate_matrix_norm,
    granges = granges,
    mode = mode,
    bin_size = bin_size,
    upstream = upstream,
    downstream = downstream,
    middle = middle,
    ignore_strand = ignore_strand,
    norm_func = norm_func,
    remove_top = 0
  )

  if (is.null(bg_bwfiles)) {
    values_list <- future_map(bwfiles, calculate_matrix_norm_fixed, bg_bw = NULL)
  } else {
    values_list <- future_map2(bwfiles, bg_bwfiles, calculate_matrix_norm_fixed)
  }

  values_list
}


#' Calculate profile values of a bigWig file over a BED file.
#'
#' Calculates profile values of a set of tracks over the loci speficied in a
#' BED file. Data points are taken each bin_size base pairs.
#'
#' Loci are aligned depending on mode parameter:
#'
#' - stretch. Aligns all starts and all ends, sort of stretching the loci.
#' The median of these lenghts is taken as the pseudo-length in order to show
#' a realistic plot when displayed.
#'
#' - start. All loci are aligned by start.
#'
#' - end. All loci are aligned by end.
#'
#' - center. All loci are aligned by center.
#'
#' @param loci BED file or GRanges to summarize
#' @param mode How to handle differences in lengths across loci:
#'
#'   stretch: Anchor each locus on both sides.
#'
#'   start: Anchor all loci on start.
#'
#'   end: Anchor all loci on end.
#'
#'   center: Center all loci.
#'
#' @param bin_size Bin size. Length of bin in base pairs. The lower, the higher
#'   the resolution.
#' @param upstream Number of base pairs to include upstream of loci.
#' @param downstream Number of base pairs to include downstream of loci.
#' @param middle Number of base pairs that the middle section has (in stretch
#'  mode). If not provided, median length of all loci is used.
#' @param ignore_strand Whether to use strand information in BED file.
#' @inheritParams bw_bins
#' @return a data frame in long format
#' @export
bw_profile <- function(bwfiles,
                       bg_bwfiles = NULL,
                       loci = NULL,
                       labels = NULL,
                       mode = "stretch",
                       bin_size = 100,
                       upstream = 2500,
                       downstream = 2500,
                       middle = NULL,
                       ignore_strand = FALSE,
                       norm_mode = "fc",
                       remove_top = 0) {

  validate_filelist(bwfiles)
  validate_locus_parameter(loci)
  granges <- loci_to_granges(loci)
  norm_func <- .process_norm_mode(norm_mode)

  validate_profile_parameters(bin_size, upstream, downstream)

  if (is.null(labels)) {
    labels <- make_label_from_object(bwfiles)
  }

  if (length(bwfiles) != length(labels)) {
    stop("labels and bwfiles must have the same length")
  }

  calculate_bw_profile_fixed <- purrr::partial(.calculate_bw_profile,
    granges = granges,
    mode = mode,
    bin_size = bin_size,
    upstream = upstream,
    downstream = downstream,
    middle = middle,
    ignore_strand = ignore_strand,
    norm_func = norm_func,
    remove_top = remove_top
  )

  if (is.null(bg_bwfiles)) {
    values_list <- furrr::future_map2(bwfiles, labels,
      calculate_bw_profile_fixed,
      bg_bw = NULL
    )
  } else {
    values_list <- furrr::future_pmap(
      list(bwfiles, bg_bwfiles, labels),
      calculate_bw_profile_fixed
    )
  }

  values <- do.call(rbind, values_list)

  values
}


#' Build a unscored bins GRanges object.
#'
#' Build a GRanges of bins of a given size, for a specific genome. Supported
#' genomes (required for the package): mm9, mm10, hg38.
#'
#' @param bin_size Bin size.
#' @param genome Genome. Supported: mm9, mm10, hg38, hg38_latest.
#' @importFrom GenomicRanges tileGenome
#' @importFrom GenomeInfoDb Seqinfo seqlengths
#' @return A GRanges object
#' @export
build_bins <- function(bin_size = 10000, genome = "mm9") {
  seqinfo <- seqlengths(Seqinfo(genome = genome))
  tileGenome(seqinfo, tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)
}

# Helpers ---------------------------------------------------


#' Score a GRanges object against a BigWig file
#'
#' Build a scored GRanges object from a GRanges object and a single BigWig file.
#' The aggregating function can be min, max, sd, mean.
#'
#' @param bwfile Path to a single BigWig file to be summarized.
#' @param granges GRanges object file to be summarized.
#' @importFrom rtracklayer BigWigFile
#' @importFrom IRanges subsetByOverlaps
#' @importFrom methods getMethod
#' @importFrom utils download.file
#' @inheritParams bw_bins
#' @return GRanges with column score.
.bw_ranges <- function(bwfile,
                       granges,
                       per_locus_stat = "mean",
                       selection = NULL) {

  valid_bwfile <- bwfile
  if (RCurl::url.exists(bwfile)) {
    valid_bwfile <- tempfile()
    download.file(bwfile, valid_bwfile)
  }

  bw <- BigWigFile(valid_bwfile)
  explicit_summary <- getMethod("summary", "BigWigFile")

  if (!is.null(selection)) {
    granges <- subsetByOverlaps(granges, selection)
  }

  unlist(explicit_summary(bw, granges, type = per_locus_stat))
}

#' Intersect a list of bigWig files with a GRanges object
#'
#' Build a binned-scored GRanges object from a list of bigWig files. The
#' aggregating function per locus can be min, max, sd, mean.
#'
#' @param granges GRanges object to intersect
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom stats quantile
#' @inheritParams bw_bins
#' @return A sorted GRanges object.
.multi_bw_ranges <- function(bwfiles,
                             labels,
                             granges,
                             per_locus_stat = "mean",
                             selection = NULL,
                             remove_top = 0) {

  if (length(bwfiles) != length(labels)) {
    stop("BigWig file list and column names must have the same length.")
  }

  summaries <- furrr::future_map(bwfiles, .bw_ranges,
    granges = granges,
    per_locus_stat = per_locus_stat,
    selection = selection
  )

  # granges_cbind sorts each element so it's safer to merge and no need to
  # sort after
  result <- granges_cbind(summaries, labels)

  # Include names if granges has them
  if ("name" %in% names(mcols(granges))) {
    result$name <- granges$name
  }

  result <- remove_top_by_mean(
    result, remove_top,
    !names(mcols(result)) %in% c("name")
  )

  result$ranges
}

#' Intersect a list of bigWig files with a GRanges object and aggregate by name
#'
#' @inheritParams bw_loci
#' @param granges GRanges object to summarize. Should have a valid name field.
#' @return An aggregated dataframe
.multi_bw_ranges_aggregated <- function(bwfiles,
                                        labels,
                                        granges,
                                        per_locus_stat,
                                        aggregate_by,
                                        remove_top) {

  result <- .multi_bw_ranges(bwfiles, labels,
    granges = granges,
    per_locus_stat = per_locus_stat,
    remove_top = remove_top
  )

  df <- .aggregate_scores(
    result,
    group_col = "name",
    aggregate_by = aggregate_by
  )

  natural_sort_by_field(df, "name")
}


#' Intersect a list of bigWig files with a GRanges object (with background)
#'
#' Intersect a list of bigWig files with a GRanges object and normalize values
#' with a set of background bigWig files.
#'
#' @param bwfilelist BigWig file list to be summarized.
#' @param bg_bwfilelist Background BigWig files.
#' @param granges GRanges object to intersect
#' @param selection A GRanges object to restrict analysis to.
#' @param norm_func Function to apply after bw/bg_bw.
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom stats quantile
#' @inheritParams bw_bins
#' @return a sorted GRanges object
.multi_bw_ranges_norm <- function(bwfilelist,
                                  bg_bwfilelist,
                                  labels,
                                  granges,
                                  per_locus_stat = "mean",
                                  selection = NULL,
                                  norm_func = identity,
                                  remove_top = 0) {
  if (length(bwfilelist) != length(bg_bwfilelist)) {
    stop("Background and signal bwfile lists must have the same length.")
  }

  result <- .multi_bw_ranges(bwfilelist, labels, granges,
    per_locus_stat = per_locus_stat,
    selection = selection
  )

  bg <- .multi_bw_ranges(bg_bwfilelist, labels, granges,
    per_locus_stat = per_locus_stat,
    selection = selection
  )

  result_df <- data.frame(result)
  bg_df <- data.frame(bg)
  result_df[, labels] <- norm_func(result_df[, labels] / bg_df[, labels])

  result <- makeGRangesFromDataFrame(result_df, keep.extra.columns = TRUE)
  result <- remove_top_by_mean(result, remove_top, labels)

  result$ranges
}

# FIXME: importFrom tidyselect where does not work according to:
# https://github.com/r-lib/tidyselect/issues/201
# This removes the note in R CMD CHECK
utils::globalVariables("where")

#' Aggregate scores of a GRanges object on a field
#'
#' Aggregates scores of a GRanges object on a specific field. Used for summary
#' functions.
#'
#' @param scored_granges A GRanges object with numerical metadata columns
#' @param group_col A column among the mcols that can be seen as a factor.
#' @param aggregate_by Function used to aggregate: mean, median, true_mean.
#'
#'     true_mean: Mean coverage taking all elements in a class as one large bin.
#'
#'     mean: mean-of-distribution approach. The mean of the aggregated value per
#'       locus is reported.
#'
#'     median: median-of-distribution. The median of the aggregated value per
#'       locus is reported.
#'
#' This does not pass the checks but does not break the build process either
#' so I left it as is, so the note thrown by build is noted everytime.
#'
#' @importFrom dplyr group_by_at summarise across select `%>%`
#' @importFrom rtracklayer mcols
#' @return A data frame with aggregated scores.
.aggregate_scores <- function(scored_granges, group_col, aggregate_by) {
  validate_group_col(scored_granges, group_col)

  score_cols <- names(mcols(scored_granges))
  score_cols <- score_cols[!score_cols %in% c(group_col)]

  df <- data.frame(scored_granges) %>%
    select(c(score_cols, group_col, width))

  validate_categories(df[, group_col])

  if (aggregate_by == "true_mean") {
    sum_vals <- df[, score_cols, drop = F] * df$width
    colnames(sum_vals) <- score_cols
    sum_vals[, group_col] <- df[, group_col]
    sum_vals$width <- df$width

    # Summarize SUM only
    sum_vals <- sum_vals %>%
      group_by_at(group_col) %>%
      summarise(across(where(is.numeric), sum))

    # Divide sum(scores) by sum(length) and keep only scores
    df <- sum_vals[, score_cols, drop = F] / sum_vals$width
    df[, group_col] <- sum_vals[, group_col]
  } else if (aggregate_by %in% c("mean", "median")) {
    f <- get(aggregate_by)
    df <- df %>%
      group_by_at(group_col) %>%
      summarise(across(where(is.numeric), f))
  } else {
    stop(paste("Function not implemented as aggregate_by:", aggregate_by))
  }

  score_cols <- score_cols[!score_cols %in% c("width")]
  df <- df[, c(score_cols, group_col), drop = FALSE]
  data.frame(df)
}


#' Calculate profile values of a bigWig file over GRanges object.
#'
#' Calculates profile values of a bigWig file over the loci speficied in a
#' GRanges object. Data points are taken each bin_size base pairs.
#'
#' Loci are aligned depending on mode parameter:
#'
#' - stretch. Aligns all starts and all ends, sort of stretching the loci.
#' The median of these lenghts is taken as the pseudo-length in order to show
#' a realistic plot when displayed.
#'
#' - start. All loci are aligned by start.
#'
#' - end. All loci are aligned by end.
#'
#' - center. All loci are aligned by center.
#'
#' Adapted from seqplots rtracklayer wrapper functions to compute coverage:
#' For more on seqplots:
#' https://www.bioconductor.org/packages/release/bioc/html/seqplots.html
#'
#' @param bw BigWig file to be summarized.
#' @param granges GRanges object
#' @param bg_bw BigWig file to be used as background.
#' @param label Name to give to the values
#' @param norm_func Normalization function
#' @importFrom rtracklayer BigWigFile import
#' @importFrom utils download.file
#' @inheritParams bw_profile
#' @return A DataFrame with the aggregated scores
.calculate_bw_profile <- function(bw,
                                  granges,
                                  bg_bw = NULL,
                                  label = NULL,
                                  mode = "stretch",
                                  bin_size = 100,
                                  upstream = 2500,
                                  downstream = 2500,
                                  middle = NULL,
                                  ignore_strand = FALSE,
                                  norm_func = identity,
                                  remove_top = 0) {
  if (is.null(label)) {
    label <- basename(bw)
  }

  fg <- .calculate_matrix_norm(bw, granges,
    mode = mode,
    bin_size = bin_size,
    upstream = upstream,
    downstream = downstream,
    middle = middle,
    ignore_strand = ignore_strand,
    norm_func = identity,
    remove_top = remove_top
  )

  fg_sum <- .summarize_matrix(fg, label)
  full <- fg_sum

  if (! is.null(bg_bw)) {
    bg <- .calculate_matrix_norm(
      bg_bw,
      granges,
      mode = mode,
      bin_size = bin_size,
      upstream = upstream,
      downstream = downstream,
      middle = middle,
      ignore_strand = ignore_strand,
      norm_func = identity,
      remove_top = remove_top
    )

    bg_sum <- .summarize_matrix(bg, "bg")

    # FIXME: median and sderror? Does sderror even make sense here?
    full <- data.frame(index = fg_sum$index, sample = fg_sum$sample,
                       mean = norm_func(fg_sum$mean / bg_sum$mean),
                       sderror = norm_func(fg_sum$mean / bg_sum$mean),
                       median = norm_func(fg_sum$median / bg_sum$median))
  }

  full
}

#' Calculate a normalized heatmap matrix for a bigWig file over a BED file
#'
#' @inheritParams .calculate_bw_profile
#' @param norm_func Normalization function
#' @importFrom stats quantile
.calculate_matrix_norm <- function(bw,
                                   granges,
                                   bg_bw = NULL,
                                   mode = "stretch",
                                   bin_size = 100,
                                   upstream = 2500,
                                   downstream = 2500,
                                   middle = NULL,
                                   ignore_strand = FALSE,
                                   norm_func = identity,
                                   remove_top = 0) {
  if (mode == "stretch") {
    full <- .calculate_stretch_matrix(bw, granges,
      bin_size = bin_size,
      upstream = upstream,
      downstream = downstream,
      middle = middle,
      ignore_strand = ignore_strand
    )

    if (!is.null(bg_bw)) {
      bg <- .calculate_stretch_matrix(bg_bw, granges,
        bin_size = bin_size,
        upstream = upstream,
        downstream = downstream,
        ignore_strand = ignore_strand
      )

      full <- norm_func(full / bg)
    }
  } else {
    start_pos <- GenomicRanges::resize(granges, 1, fix = mode)
    granges <- GenomicRanges::promoters(start_pos, upstream, downstream)

    # To properly center one needs to floor separately upstream and downstream.
    # This way the tick will always be in between bins.
    npoints <- floor(upstream / bin_size) + floor(downstream / bin_size)

    full <- .intersect_bw_and_granges(
      bw,
      granges,
      npoints = npoints,
      ignore_strand = FALSE
    )

    if (!is.null(bg_bw)) {
      bg <- .intersect_bw_and_granges(
        bg_bw,
        granges,
        npoints = npoints,
        ignore_strand = FALSE
      )

      full <- norm_func(full / bg)
    }
  }

  if (remove_top > 0) {
    top_quantile <- quantile(rowMeans(full), probs = c(1 - remove_top))
    full <- full[rowMeans(full) <= top_quantile, ]
  }

  full
}


#' Calculate matrix for stretch mode
#'
#' @param bw BigWigFile object
#' @param granges GRanges object
#' @param bin_size Bin size
#' @param upstream Number of basepairs upstream
#' @param downstream Number of basepairs downstream
#' @param middle Number of base pairs that the middle section has. If not
#'   provided, median is used.
#' @param ignore_strand Ignore strand (bool)
#'
#' @return Summary matrix
.calculate_stretch_matrix <- function(bw,
                                      granges,
                                      bin_size = 100,
                                      upstream = 2500,
                                      downstream = 2500,
                                      middle = NULL,
                                      ignore_strand = FALSE) {
  left_npoints <- floor(upstream / bin_size)
  right_npoints <- floor(downstream / bin_size)

  if (is.null(middle)) {
    # Stretch to the median value of the GR object
    middle <- floor(median(GenomicRanges::width(granges)))
  }

  middle_npoints <- floor(middle / bin_size)

  left <- .intersect_bw_and_granges(bw,
    GenomicRanges::flank(granges, upstream, start = TRUE),
    npoints = left_npoints,
    ignore_strand = ignore_strand
  )

  right <- .intersect_bw_and_granges(bw,
    GenomicRanges::flank(granges, downstream, start = FALSE),
    npoints = right_npoints,
    ignore_strand = ignore_strand
  )

  middle <- .intersect_bw_and_granges(bw,
    granges,
    npoints = middle_npoints,
    ignore_strand = ignore_strand
  )

  cbind(left, middle, right)
}


#' Intersect a BigWig file over loci on a GRanges object
#'
#' Intersect a BigWig file over loci on a GRanges object, taking npoints points
#' per locus.
#'
#' @param bw BigWigFile object (rtracklayer).
#' @param granges GRanges object.
#' @param npoints How many points to take (different to bin size!).
#' @param ignore_strand Ignore strand information in granges.
#' @return A value matrix of dimensions len(granges) x npoints.
.intersect_bw_and_granges <- function(bw,
                                      granges,
                                      npoints,
                                      ignore_strand = FALSE) {
  bwfile <- fetch_bigwig(bw)

  values <- rtracklayer::summary(bwfile,
    which = granges,
    as = "matrix",
    size = npoints
  )

  # Reverse minus strand rows
  if (!ignore_strand) {
    indices <- seq(ncol(values), 1)
    values[as.character(GenomicRanges::strand(granges)) == "-", ] <-
      values[as.character(GenomicRanges::strand(granges)) == "-", indices]
  }

  values
}


#' Get a function out of a norm mode string.
#'
#' @param mode String representing mode. Valid values are fc and log2fc
#'
#' @return A function
.process_norm_mode <- function(mode) {
  switch(mode,
    "fc" = identity,
    "log2fc" = log2
  )
}


#' Summarize a intersect_bw_and_granges matrix
#'
#' Compute averages and standard error values for a matrix returned by
#' intersect_bw_and_granges.
#'
#' @param matrix A matrix returned by .intersect_bw_and_granges.
#' @param label Label for the sample.
#' @importFrom stats median sd
#' @return A data frame with summarized values, stderr and medians, plus label.
.summarize_matrix <- function(matrix, label) {
  # Ignore Inf and NaN in the computation of means/SD
  matrix[is.infinite(matrix)] <- NA

  omitted_vals <- sum(is.na(matrix))
  if (omitted_vals > 100) {
    mean_per_locus <- omitted_vals / nrow(matrix)
    warning(paste(
      "Profile plot:",
      omitted_vals, "generated (",
      mean_per_locus, "per locus)"
    ))
  }

  df <- data.frame(
    mean = colMeans(matrix, na.rm = TRUE),
    sderror = apply(
      matrix, 2,
      function(n) {
        sd(n, na.rm = TRUE) / sqrt(sum(!is.na(n)))
      }
    ),
    median = apply(matrix, 2, median, na.rm = TRUE)
  )

  df$index <- as.integer(rownames(df))
  df$sample <- label
  df
}
