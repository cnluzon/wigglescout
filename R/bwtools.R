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

#' Calculate bigWig global coverage.
#' @param bwfile bigWig file
#' @param default_na Value to replace missing values.
#' @importFrom dplyr mutate select `%>%`
#' @importFrom rtracklayer BigWigFile
#' @export
bw_global_coverage <- function(bwfile, default_na = NA_real_) {
    bw <- .fetch_bigwig(bwfile)
    if (!is.null(bw)) {
        explicit_summary <- getMethod("summary", "BigWigFile")
        df <- data.frame(
            unlist(explicit_summary(bw, type="mean", default_na = default_na))
        )
        result <- sum(df %>%
                          mutate(weighted=.data[["score"]]*.data[["width"]]) %>%
                          select("weighted")
                      ) / sum(df$width)
        result
    }
}


#' Score a bigWig file list and a BED file or GRanges object.
#'
#' Build a scored GRanges object from a bigWig file list and a BED file or
#' GRanges object. The aggregating function (per locus) can be min, max,
#' sd, mean.
#'
#' Values can be normalized using background bigWig files (usually input
#' files). By default, the value obtained will be bigwig / bg_bigwig per locus,
#' per bigWig.
#'
#' bwfiles and bg_bwfiles must have the same length. If you are using same
#' background for several files, then file paths must be repeated accordingly.
#'
#' norm_mode can be either "fc", where values will represent bigwig / bg_bigwig,
#' or "log2fc", values will represent log2(bigwig / bg_bigwig) per locus.
#'
#' @param loci GRanges or BED file to summarize the BigWig file at.
#' @param aggregate_by Statistic to aggregate per group. If NULL, values are
#'    not aggregated. This is the behavior by default.
#' @examples
#' # Get the raw files
#' bed <- system.file("extdata", "sample_genes_mm9.bed", package="wigglescout")
#' bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
#' bw2 <- system.file("extdata",
#'                    "sample_H3K9me3_ChIP.bw", package="wigglescout")
#' bw_inp <- system.file("extdata", "sample_Input.bw", package="wigglescout")
#'
#' # Run single bw with single bed
#' bw_loci(bw, bed)
#'
#' # Use of some parameters
#' bw_loci(bw, bed, labels = c("H33"), remove_top = 0.01)
#'
#' # Log2 fold change
#' bw_loci(bw, loci = bed, bg_bwfiles = bw_inp, norm_mode = "log2fc")
#'
#' # Multiple bigWig
#' bw_loci(c(bw, bw2),
#'         loci = bed,
#'         bg_bwfiles = c(bw_inp, bw_inp),
#'         norm_mode = "log2fc")
#' @export
#' @inheritParams bw_bins
#' @return A GRanges object with each bwfile as a metadata column named
#'     after labels, if provided, or after filenames otherwise.
#' @importFrom rtracklayer import BigWigFile
#' @importFrom GenomeInfoDb sortSeqlevels
bw_loci <- function(bwfiles,
                    loci,
                    bg_bwfiles = NULL,
                    labels = NULL,
                    per_locus_stat = "mean",
                    aggregate_by = NULL,
                    norm_mode = "fc",
                    remove_top = 0,
                    default_na = NA_real_,
                    scaling = "none") {

    .validate_filelist(bwfiles)
    .validate_locus_parameter(loci)
    norm_func <- .process_norm_mode(norm_mode)

    if (is.null(labels)) {
        labels <- .make_label_from_object(bwfiles)
    } else {
        # Ensures later on we only try to access valid labels
        labels <- make.names(labels)
    }

    bed <- .loci_to_granges(loci)

    result <- NULL
    if (is.null(aggregate_by)) {
        if (is.null(bg_bwfiles)) {
            result <- .multi_bw_ranges(
                bwfiles, labels,
                granges = bed,
                per_locus_stat = per_locus_stat,
                remove_top = remove_top,
                default_na = default_na,
                scaling = scaling
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
                default_na = default_na,
                scaling = scaling
            )
        }
    } else {
        labels_bg <- NULL
        if (!is.null(bg_bwfiles)) {
            labels_bg <- paste0("bg_", labels)
        }

        result <- .multi_bw_ranges_aggregated(
            c(bwfiles, bg_bwfiles),
            granges = bed,
            labels = c(labels, labels_bg),
            per_locus_stat = per_locus_stat,
            aggregate_by = aggregate_by,
            remove_top = remove_top,
            default_na = default_na,
            scaling = scaling
        )

        if (!is.null(bg_bwfiles)) {
            rows <- rownames(result)
            result <- data.frame(
                norm_func(result[rows, labels] / result[rows, labels_bg])
            )
            rownames(result) <- rows
            colnames(result) <- labels
        }
    }

    result
}


#' Build a binned-scored GRanges object from a set of bigWig files
#'
#' Build a binned-scored GRanges object from one or many bigWig files.
#' The aggregating function per bin can be min, max, sd, mean.
#'
#' bwfiles and bg_bwfiles must have the same length. If you are using same
#' background for several files, then file paths must be repeated accordingly.
#'
#' norm_mode can be either "fc", where values will represent bigwig / bg_bigwig,
#' or "log2fc", values will represent log2(bigwig / bg_bigwig) per bin.
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
#' @param scaling If none, no operation is performed (default). If relative,
#'    values are divided by global mean (1x genome coverage).
#' @param default_na Default value for missing values
#' @return A GRanges object with each bwfile as a metadata column named
#'     after labels, if provided, or after filenames otherwise.
#' @export
#' @examples
#' # Get the raw files
#' bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
#' bw2 <- system.file("extdata",
#'                    "sample_H3K9me3_ChIP.bw", package="wigglescout")
#' bw_inp <- system.file("extdata", "sample_Input.bw", package="wigglescout")
#'
#' # Sample bigWig files only have valid values on this region
#' locus <- GenomicRanges::GRanges(
#'   seqnames = "chr15",
#'   IRanges::IRanges(102600000, 103100000)
#' )
#'
#' # Run single bw. The larger the bin size, the faster this is.
#' bw_bins(bw, bin_size = 100000, labels = c("H33"), selection = locus)
#'
#' # Use of some parameters
#' bw_bins(bw, bg_bwfiles = bw_inp, bin_size = 100000,
#'        labels = c("H33"), norm_mode = "log2fc", selection = locus)
#'
#' # Multiple bigWig
#' bw_bins(c(bw, bw2),
#'         bin_size = 100000,
#'         bg_bwfiles = c(bw_inp, bw_inp),
#'         norm_mode = "log2fc",
#'         selection = locus)
bw_bins <- function(bwfiles,
                    bg_bwfiles = NULL,
                    labels = NULL,
                    per_locus_stat = "mean",
                    bin_size = 10000,
                    genome = "mm9",
                    selection = NULL,
                    norm_mode = "fc",
                    remove_top = 0,
                    default_na = NA_real_,
                    scaling = "none") {

    .validate_filelist(bwfiles)
    norm_func <- .process_norm_mode(norm_mode)

    if (is.null(labels)) {
        labels <- .make_label_from_object(bwfiles)
    }

    tiles <- build_bins(bin_size = bin_size, genome = genome)

    if (is.null(bg_bwfiles)) {
        result <- .multi_bw_ranges(
            bwfiles,
            labels,
            tiles,
            per_locus_stat = per_locus_stat,
            selection = selection,
            remove_top = remove_top,
            default_na = default_na,
            scaling = scaling
        )
    } else {
        result <- .multi_bw_ranges_norm(
            bwfiles,
            bg_bwfiles,
            labels,
            tiles,
            per_locus_stat = per_locus_stat,
            selection = selection,
            norm_func = norm_func,
            remove_top = remove_top,
            default_na = default_na,
            scaling = scaling
        )
    }
    result
}

#' Calculate heatmap matrix of a bigWig file over a GRanges or BED file
#'
#' Calculates a matrix of values that corresponds to the usual heatmaps where
#' each row is a locus and columns are data points taken each bin_size base
#' pairs.
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
#' @inheritParams bw_profile
#' @importFrom furrr future_map future_map2
#' @importFrom rtracklayer import
#' @importFrom purrr partial
#' @return A list of matrices where each element correspond to each bigWig file.
#' @export
#' @examples
#' # Get the raw files
#' bed <- system.file("extdata", "sample_genes_mm9.bed", package="wigglescout")
#' bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
#' bw2 <- system.file("extdata",
#'                    "sample_H3K9me3_ChIP.bw", package="wigglescout")
#'
#' # Heatmaps with a single bigWig
#' h <- bw_heatmap(bw, loci = bed, mode = "stretch")
#'
#' # h is a list
#' h[[1]]
#'
#' # Heatmaps with multiple bigWig
#' h <- bw_heatmap(c(bw, bw2), loci = bed, mode = "stretch")
#'
#' # h is a list
#' h[[2]]
#'
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
                        norm_mode = "fc",
                        default_na = NA_real_,
                        scaling = "none") {

    .validate_filelist(bwfiles)
    .validate_locus_parameter(loci)
    granges <- .loci_to_granges(loci)
    norm_func <- .process_norm_mode(norm_mode)

    .validate_profile_parameters(bin_size, upstream, downstream)

    if (is.null(labels)) {
        labels <- basename(bwfiles)
    }

    if (length(bwfiles) != length(labels)) {
        stop("labels and bwfiles must have the same length")
    }

    calculate_matrix_norm_fixed <- partial(
        .calculate_matrix_norm,
        granges = granges,
        mode = mode,
        bin_size = bin_size,
        upstream = upstream,
        downstream = downstream,
        middle = middle,
        ignore_strand = ignore_strand,
        norm_func = norm_func,
        remove_top = 0,
        default_na = default_na,
        scaling = scaling
    )

    if (is.null(bg_bwfiles)) {
        values_list <- future_map(
            bwfiles,
            calculate_matrix_norm_fixed,
            bg_bw = NULL
        )
    } else {
        values_list <- future_map2(
            bwfiles,
            bg_bwfiles,
            calculate_matrix_norm_fixed
        )
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
#' @param scaling Whether to use the bigWig values as they are (none - default)
#'    or calculate relative enrichment (relative) by dividing values by global
#'    average.
#' @inheritParams bw_bins
#' @return a data frame in long format
#' @importFrom purrr partial
#' @importFrom furrr future_pmap future_map2
#' @export
#' @examples
#' # Get the raw files
#' bed <- system.file("extdata", "sample_genes_mm9.bed", package="wigglescout")
#' bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
#' bw2 <- system.file("extdata",
#'                    "sample_H3K9me3_ChIP.bw", package="wigglescout")
#'
#' # Profiles are returned in long format and include mean, median and stderror
#' bw_profile(bw, loci = bed, mode = "stretch")
#'
#' bw_profile(bw, loci = bed, mode = "start",
#'            upstream = 1000, downstream = 1500)
#'
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
                        remove_top = 0,
                        default_na = NA_real_,
                        scaling = "none") {

    .validate_filelist(bwfiles)
    .validate_locus_parameter(loci)
    granges <- .loci_to_granges(loci)
    norm_func <- .process_norm_mode(norm_mode)
    .validate_profile_parameters(bin_size, upstream, downstream)

    if (is.null(labels)) {
        labels <- .make_label_from_object(bwfiles)
    }

    if (length(bwfiles) != length(labels)) {
        stop("labels and bwfiles must have the same length")
    }

    calculate_bw_profile_fixed <- partial(.calculate_bw_profile,
        granges = granges,
        mode = mode,
        bin_size = bin_size,
        upstream = upstream,
        downstream = downstream,
        middle = middle,
        ignore_strand = ignore_strand,
        norm_func = norm_func,
        remove_top = remove_top,
        default_na = default_na,
        scaling = scaling
    )

    if (is.null(bg_bwfiles)) {
        values_list <- future_map2(bwfiles, labels,
            calculate_bw_profile_fixed,
            bg_bw = NULL
        )
    } else {
        values_list <- future_pmap(list(bwfiles, bg_bwfiles, labels),
            calculate_bw_profile_fixed
        )
    }
    do.call(rbind, values_list)
}


#' Build a unscored bins GRanges object.
#'
#' Build a GRanges of bins of a given size, for a specific genome. Supported
#' genomes rely on \link[GenomeInfoDb]{Seqinfo}. This requires internet access
#' to work.
#'
#' @param bin_size Bin size.
#' @param genome Genome. Supported: mm9, mm10, hg38, hg38_latest.
#' @importFrom GenomicRanges tileGenome
#' @importFrom GenomeInfoDb Seqinfo seqlengths
#' @return A GRanges object with a tiled genome
#' @export
#' @examples
#'
#' build_bins(bin_size = 50000, genome = "mm9")
#' build_bins(bin_size = 50000, genome = "hg38")
#' build_bins(bin_size = 50000, genome = "mm10")
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
#' @param default_na Default value for missing values
#' @param scaling If none, no operation is performed (default). If relative,
#'    values are divided by global mean (1x genome coverage).
#' @importFrom rtracklayer BigWigFile
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges mcols
#' @importFrom methods getMethod
#' @importFrom utils download.file
#' @importFrom RCurl url.exists
#' @inheritParams bw_bins
#' @return GRanges with column score.
.bw_ranges <- function(bwfile, granges,
                        per_locus_stat = "mean",
                        selection = NULL,
                        default_na = NA_real_,
                        scaling = "none") {
    bw <- .fetch_bigwig(bwfile)
    if (!is.null(bw)) {
        explicit_summary <- getMethod("summary", "BigWigFile")

        if (!is.null(selection)) {
            granges <- subsetByOverlaps(granges, selection)
        }
        result <- unlist(explicit_summary(bw, granges, type = per_locus_stat, defaultValue = default_na))

        if (scaling == "relative") {
            global_mean <- bw_global_coverage(bwfile)
            values <- data.frame(mcols(result)) / global_mean

            GenomicRanges::mcols(result)[, colnames(values)] <- values
        }
        result
    }
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
#' @importFrom furrr future_map
#' @inheritParams bw_bins
#' @return A sorted GRanges object.
.multi_bw_ranges <- function(bwfiles, labels, granges,
                                per_locus_stat = "mean",
                                selection = NULL,
                                remove_top = 0,
                                default_na = NA_real_,
                                scaling = "none") {

    if (length(bwfiles) != length(labels)) {
        stop("BigWig file list and column names must have the same length.")
    }

    summaries <- future_map(
        bwfiles,
        .bw_ranges,
        granges = granges,
        per_locus_stat = per_locus_stat,
        selection = selection,
        default_na = default_na,
        scaling = scaling
    )

    result <- .granges_left_join(summaries, labels, granges)

    result <- .remove_top_by_mean(
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
.multi_bw_ranges_aggregated <- function(bwfiles, labels, granges,
                                        per_locus_stat,
                                        aggregate_by,
                                        remove_top,
                                        default_na = NA_real_,
                                        scaling = "none") {

    result <- .multi_bw_ranges(
        bwfiles,
        labels,
        granges = granges,
        per_locus_stat = per_locus_stat,
        remove_top = remove_top,
        default_na = default_na,
        scaling = scaling
    )

    df <- .aggregate_scores(
        result,
        group_col = "name",
        aggregate_by = aggregate_by
    )

    .natural_sort_by_field(df, "name")
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
                                    remove_top = 0,
                                    default_na = NA_real_,
                                    scaling = "none") {
    if (length(bwfilelist) != length(bg_bwfilelist)) {
        stop("Background and signal bwfile lists must have the same length.")
    }

    result <- .multi_bw_ranges(
        bwfilelist,
        labels,
        granges,
        per_locus_stat = per_locus_stat,
        selection = selection,
        default_na = default_na,
        scaling = scaling
    )

    bg <- .multi_bw_ranges(
        bg_bwfilelist,
        labels,
        granges,
        per_locus_stat = per_locus_stat,
        selection = selection,
        default_na = default_na,
        scaling = scaling
    )

    result_df <- data.frame(result)
    bg_df <- data.frame(bg)
    result_df[, labels] <- norm_func(result_df[, labels] / bg_df[, labels])

    result <- makeGRangesFromDataFrame(result_df, keep.extra.columns = TRUE)
    result <- .remove_top_by_mean(result, remove_top, labels)

    result$ranges
}

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
#' @importFrom dplyr group_by summarise across select if_all `%>%` everything
#' @importFrom tidyselect all_of where
#' @importFrom rtracklayer mcols
#' @return A data frame with aggregated scores.
.aggregate_scores <- function(scored_granges, group_col, aggregate_by) {
    .validate_group_col(scored_granges, group_col)

    score_cols <- names(mcols(scored_granges))
    score_cols <- score_cols[!score_cols %in% c(group_col)]

    df <- data.frame(scored_granges) %>%
        select(tidyselect::all_of(c(score_cols, group_col, "width")))

    .validate_categories(df[, group_col])

    if (aggregate_by == "true_mean") {
        sum_vals <- df[, score_cols, drop = FALSE] * df$width
        colnames(sum_vals) <- score_cols
        sum_vals[, group_col] <- df[, group_col]
        sum_vals$width <- df$width

        # Summarize SUM only
        # In this case one would expect that it's either all NAs or none,
        # since NA's should come from non-existent values in a given experiment.
        # We expect bw tracks to be filled with zeros.
        sum_vals <- sum_vals %>%
            dplyr::group_by(.data[[group_col]]) %>%
            dplyr::filter(if_all(everything(), ~ !is.na(.))) %>%
            dplyr::summarise(across(where(is.numeric), sum))


        # Divide sum(scores) by sum(length) and keep only scores
        df <- sum_vals[, score_cols, drop = FALSE] / sum_vals$width
        df[, group_col] <- sum_vals[, group_col]
    } else if (aggregate_by %in% c("mean", "median")) {
        f <- get(aggregate_by)
        df <- df %>%
            dplyr::group_by(.data[[group_col]]) %>%
            dplyr::filter(if_all(everything(), ~ !is.na(.))) %>%
            dplyr::summarise(across(where(is.numeric), f))
    } else {
        stop("Function not implemented as aggregate_by: ", aggregate_by)
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
#' @param default_na Default value for missing values
#' @importFrom rtracklayer BigWigFile import
#' @importFrom utils download.file
#' @inheritParams bw_profile
#' @return A DataFrame with the aggregated scores
.calculate_bw_profile <- function(bw, granges,
                                    bg_bw = NULL,
                                    label = NULL,
                                    mode = "stretch",
                                    bin_size = 100,
                                    upstream = 2500,
                                    downstream = 2500,
                                    middle = NULL,
                                    ignore_strand = FALSE,
                                    norm_func = identity,
                                    remove_top = 0,
                                    default_na = NA_real_,
                                    scaling = "none") {
    if (is.null(label)) {
        label <- basename(bw)
    }

    fg <- .calculate_matrix_norm(
        bw,
        granges,
        mode = mode,
        bin_size = bin_size,
        upstream = upstream,
        downstream = downstream,
        middle = middle,
        ignore_strand = ignore_strand,
        norm_func = identity,
        remove_top = remove_top,
        default_na = default_na,
        scaling = scaling
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
            remove_top = remove_top,
            default_na = default_na,
            scaling = scaling
        )

        bg_sum <- .summarize_matrix(bg, "bg")

        # FIXME: median and sderror? Does sderror even make sense here?
        full <- data.frame(
            index = fg_sum$index, sample = fg_sum$sample,
            mean = norm_func(fg_sum$mean / bg_sum$mean),
            sderror = norm_func(fg_sum$mean / bg_sum$mean),
            median = norm_func(fg_sum$median / bg_sum$median)
        )
    }

    full
}

#' Calculate a normalized heatmap matrix for a bigWig file over a BED file
#'
#' @inheritParams .calculate_bw_profile
#' @param norm_func Normalization function
#' @importFrom stats quantile
#' @return Summary matrix
.calculate_matrix_norm <- function(bw, granges,
                                    bg_bw = NULL,
                                    mode = "stretch",
                                    bin_size = 100,
                                    upstream = 2500,
                                    downstream = 2500,
                                    middle = NULL,
                                    ignore_strand = FALSE,
                                    norm_func = identity,
                                    remove_top = 0,
                                    default_na = NA_real_,
                                    scaling = "none") {
    if (mode == "stretch") {
        full <- .calculate_stretch_matrix(
            bw,
            granges,
            bin_size = bin_size,
            upstream = upstream,
            downstream = downstream,
            middle = middle,
            ignore_strand = ignore_strand,
            default_na = default_na,
            scaling = scaling
        )

        if (!is.null(bg_bw)) {
            bg <- .calculate_stretch_matrix(
                bg_bw,
                granges,
                bin_size = bin_size,
                upstream = upstream,
                downstream = downstream,
                ignore_strand = ignore_strand,
                default_na = default_na,
                scaling = scaling
            )

            full <- norm_func(full / bg)
        }
    } else {
        start_pos <- GenomicRanges::resize(granges, 1, fix = mode)
        granges <- GenomicRanges::promoters(start_pos, upstream, downstream)

        # To properly center one needs to floor separately upstream and
        # downstream.
        # This way the tick will always be in between bins.
        npoints <- floor(upstream / bin_size) + floor(downstream / bin_size)

        full <- .intersect_bw_and_granges(
            bw,
            granges,
            npoints = npoints,
            ignore_strand = FALSE,
            default_na = default_na,
            scaling = scaling
        )

        if (!is.null(bg_bw)) {
            bg <- .intersect_bw_and_granges(
                bg_bw,
                granges,
                npoints = npoints,
                ignore_strand = FALSE,
                default_na = default_na,
                scaling = scaling
            )

            full <- norm_func(full / bg)
        }
    }

    if (remove_top > 0) {
        top_quantile <- quantile(rowMeans(full, na.rm = TRUE),
                                 probs = c(1 - remove_top), na.rm = TRUE)
        full <- full[rowMeans(full, na.rm = TRUE) <= top_quantile, ]
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
#' @param default_na Default value for missing values
#' @param middle Number of base pairs that the middle section has. If not
#'   provided, median is used.
#' @param ignore_strand Ignore strand (bool)
#' @param scaling none if no scaling is performed (default), relative if 1x
#'   scaling is done (divide by global coverage).
#' @importFrom GenomicRanges flank width
#'
#' @return Summary matrix
.calculate_stretch_matrix <- function(bw, granges,
                                        bin_size = 100,
                                        upstream = 2500,
                                        downstream = 2500,
                                        middle = NULL,
                                        ignore_strand = FALSE,
                                        default_na = NA_real_,
                                        scaling = "none") {

    left_npoints <- floor(upstream / bin_size)
    right_npoints <- floor(downstream / bin_size)

    if (is.null(middle)) {
        # Stretch to the median value of the GR object
        middle <- floor(median(width(granges)))
    }

    middle_npoints <- floor(middle / bin_size)

    left <- .intersect_bw_and_granges(
        bw,
        flank(granges, upstream, start = TRUE),
        npoints = left_npoints,
        ignore_strand = ignore_strand,
        default_na = default_na,
        scaling = scaling
    )

    right <- .intersect_bw_and_granges(
        bw,
       flank(granges, downstream, start = FALSE),
       npoints = right_npoints,
       ignore_strand = ignore_strand,
       default_na = default_na,
       scaling = scaling
    )

    middle <- .intersect_bw_and_granges(
        bw,
        granges,
        npoints = middle_npoints,
        ignore_strand = ignore_strand,
        default_na = default_na,
        scaling = scaling
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
#' @param default_na Value to use as default for missing values
#' @param scaling none if no scaling is done, relative if scaled to global
#'   mean of 1
#' @importFrom rtracklayer summary
#' @importFrom GenomicRanges strand
#' @return A value matrix of dimensions len(granges) x npoints.
.intersect_bw_and_granges <- function(bw, granges, npoints,
                                        ignore_strand = FALSE,
                                        default_na = NA_real_,
                                        scaling = "none") {
    bwfile <- .fetch_bigwig(bw)
    values <- summary(bwfile, which = granges, as = "matrix", size = npoints, defaultValue = default_na)

    # Reverse minus strand rows
    if (!ignore_strand) {
        indices <- seq(ncol(values), 1)
        values[as.character(strand(granges)) == "-", ] <-
            values[as.character(strand(granges)) == "-", indices]
    }

    if (scaling == "relative") {
        global = bw_global_coverage(bw)
        values = values / global
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
        warning(
            "Profile plot: ",
            omitted_vals, " generated ( ",
            mean_per_locus, " per locus)"
        )
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
