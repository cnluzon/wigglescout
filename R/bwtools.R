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
#' @param bedfile BED file to intersect with the BigWig file.
#' @param aggregate_by Statistic to aggregate per group. If NULL, values are
#'    not aggregated. This is the behavior by default.
#' @export
#' @inheritParams bw_bins
#' @importFrom rtracklayer import BigWigFile
#' @importFrom GenomeInfoDb sortSeqlevels
bw_bed <- function(bwfiles,
                   bedfile,
                   bg_bwfiles = NULL,
                   labels = NULL,
                   per_locus_stat = "mean",
                   aggregate_by = NULL,
                   norm_mode = "fc",
                   remove_top = 0) {

  validate_filelist(bwfiles)
  validate_locus_parameter(bedfile)
  norm_func <- process_norm_mode(norm_mode)

  if (is.null(labels)) {
    labels <- make_label_from_filename(bwfiles)
  }

  bed <- loci_to_granges(bedfile)

  result <- NULL
  if (is.null(aggregate_by)) {
    if (is.null(bg_bwfiles)) {
      result <- multi_bw_ranges(bwfiles, labels,
                  granges = bed,
                  per_locus_stat = per_locus_stat,
                  remove_top = remove_top
      )
    }
    else {
      # Only want to normalize per-locus if not aggregating
      result <- multi_bw_ranges_norm(
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
    result <- multi_bw_ranges_aggregated(bwfiles,
                labels = labels,
                granges = bed,
                per_locus_stat = per_locus_stat,
                aggregate_by = aggregate_by,
                remove_top = remove_top
              )

    if (!is.null(bg_bwfiles)) {
      bg <- multi_bw_ranges_aggregated(bg_bwfiles,
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

process_norm_mode <- function(mode) {
  switch(mode,
         "fc" = identity,
         "log2fc" = log2
         )
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
  norm_func <- process_norm_mode(norm_mode)

  if (is.null(labels)) {
    labels <- make_label_from_filename(bwfiles)
  }

  tiles <- build_bins(bin_size = bin_size, genome = genome)

  if (is.null(bg_bwfiles)) {
    result <- multi_bw_ranges(bwfiles, labels, tiles,
                per_locus_stat = per_locus_stat,
                selection = selection,
                remove_top = remove_top
              )
  } else {
    result <- multi_bw_ranges_norm(bwfiles, bg_bwfiles, labels, tiles,
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
                       bedfile = NULL,
                       labels = NULL,
                       mode = "stretch",
                       bin_size = 100,
                       upstream = 2500,
                       downstream = 2500,
                       middle = NULL,
                       ignore_strand = FALSE,
                       norm_mode = "fc") {

  validate_filelist(bwfiles)
  validate_locus_parameter(bedfile)
  granges <- loci_to_granges(bedfile)
  norm_func <- process_norm_mode(norm_mode)

  validate_profile_parameters(bin_size, upstream, downstream)

  if (is.null(labels)) {
    labels <- basename(bwfiles)
  }

  if (length(bwfiles) != length(labels)) {
    stop("labels and bwfiles must have the same length")
  }

  calculate_matrix_norm_fixed <- partial(calculate_matrix_norm,
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
#' @param bedfile BED file to summarize
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
                       bedfile = NULL,
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
  validate_locus_parameter(bedfile)
  granges <- loci_to_granges(bedfile)
  norm_func <- process_norm_mode(norm_mode)

  validate_profile_parameters(bin_size, upstream, downstream)

  if (is.null(labels)) {
    labels <- make_label_from_filename(bwfiles)
  }

  if (length(bwfiles) != length(labels)) {
    stop("labels and bwfiles must have the same length")
  }

  calculate_bw_profile_fixed <- purrr::partial(calculate_bw_profile,
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
    values_list <- furrr::future_pmap(list(bwfiles, bg_bwfiles, labels),
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
  seqinfo <- seqlengths(Seqinfo(genome=genome))
  tileGenome(seqinfo, tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)
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
multi_bw_ranges <- function(bwfiles,
                            labels,
                            granges,
                            per_locus_stat = "mean",
                            selection = NULL,
                            remove_top = 0) {


  if (length(bwfiles) != length(labels)) {
    stop("BigWig file list and column names must have the same length.")
  }

  summaries <- furrr::future_map(bwfiles, bw_ranges,
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

  result <- remove_top_by_mean(result, remove_top,
                               !names(mcols(result)) %in% c("name"))

  result$ranges
}

#' Intersect a list of bigWig files with a GRanges object and aggregate by name
#'
#' @inheritParams bw_bed
#' @param granges GRanges object to summarize. Should have a valid name field.
#' @return An aggregated dataframe
multi_bw_ranges_aggregated <- function(bwfiles,
                                       labels,
                                       granges,
                                       per_locus_stat,
                                       aggregate_by,
                                       remove_top) {

  result <- multi_bw_ranges(bwfiles, labels,
              granges = granges,
              per_locus_stat = per_locus_stat,
              remove_top = remove_top
            )

  df <- aggregate_scores(
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
multi_bw_ranges_norm <- function(bwfilelist,
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

  result <- multi_bw_ranges(bwfilelist, labels, granges,
              per_locus_stat = per_locus_stat,
              selection = selection
            )

  bg <- multi_bw_ranges(bg_bwfilelist, labels, granges,
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
bw_ranges <- function(bwfile,
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

  if (! is.null(selection)) {
    granges <- subsetByOverlaps(granges, selection)
  }

  unlist(explicit_summary(bw, granges, type = per_locus_stat))
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
#' FIXME: importFrom tidyselect where does not work according to:
#' https://github.com/r-lib/tidyselect/issues/201
#' This does not pass the checks but does not break the build process either
#' so I left it as is, so the note thrown by build is noted everytime.
#'
#' @importFrom dplyr group_by_at summarise across `%>%`
#' @importFrom rtracklayer mcols
#' @return A data frame with aggregated scores.
aggregate_scores <- function(scored_granges, group_col, aggregate_by) {
  validate_group_col(scored_granges, group_col)

  # GRanges objects are 1-based and inclusive [start, end]
  end_value <- GenomicRanges::end(scored_granges)
  start_value <- GenomicRanges::start(scored_granges)
  scored_granges$length <- end_value - start_value + 1

  df <- data.frame(mcols(scored_granges))
  validate_categories(df[, group_col])
  # Make sure special characters are taken into account
  score_cols <- colnames(mcols(scored_granges))
  score_cols <- make.names(score_cols)
  score_cols <- score_cols[!score_cols %in% c(group_col)]

  if (aggregate_by == "true_mean") {
    sum_vals <- df[, score_cols] * df$length
    colnames(sum_vals) <- score_cols
    sum_vals[, group_col] <- df[, group_col]
    sum_vals$length <- df$length

    # Summarize SUM only
    sum_vals <- sum_vals %>%
      group_by_at(group_col) %>%
      summarise(across(where(is.numeric), sum))

    # Divide sum(scores) by sum(length) and keep only scores
    df <- sum_vals[, score_cols] / sum_vals$length
    df[, group_col] <- sum_vals[, group_col]

  } else if (aggregate_by %in% c("mean", "median")) {
    f <- get(aggregate_by)
    df <- df %>%
      group_by_at(group_col) %>%
      summarise(across(where(is.numeric), f))

  } else {
    stop(paste("Function not implemented as aggregate_by:", aggregate_by))
  }

  score_cols <- score_cols[!score_cols %in% c("length")]
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
calculate_bw_profile <- function(bw,
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

  full <- calculate_matrix_norm(bw, granges,
            bg_bw = bg_bw,
            mode = mode,
            bin_size = bin_size,
            upstream = upstream,
            downstream = downstream,
            middle = middle,
            ignore_strand = ignore_strand,
            norm_func = identity,
            remove_top = remove_top
          )

  summarize_matrix(full, label)
}

#' Calculate a normalized heatmap matrix for a bigWig file over a BED file
#'
#' @inheritParams calculate_bw_profile
#' @param norm_func Normalization function
#' @importFrom stats quantile
calculate_matrix_norm <- function(bw,
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
    full <- calculate_stretch_matrix(bw, granges,
              bin_size = bin_size,
              upstream = upstream,
              downstream = downstream,
              middle = middle,
              ignore_strand = ignore_strand
            )

    if (!is.null(bg_bw)) {
      bg <- calculate_stretch_matrix(bg_bw, granges,
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

    full <- intersect_bw_and_granges(
      bw,
      granges,
      npoints = npoints,
      ignore_strand = FALSE
    )

    if (!is.null(bg_bw)) {
      bg <- intersect_bw_and_granges(
              bg_bw,
              granges,
              npoints = npoints,
              ignore_strand = FALSE
            )

      full <- norm_func(full / bg)
    }
  }

  if (remove_top > 0) {
      top_quantile <- quantile(rowMeans(full), probs = c(1-remove_top))
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
calculate_stretch_matrix <- function(bw,
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

  left <- intersect_bw_and_granges(bw,
            GenomicRanges::flank(granges, upstream, start = TRUE),
            npoints = left_npoints,
            ignore_strand = ignore_strand
          )

  right <- intersect_bw_and_granges(bw,
            GenomicRanges::flank(granges, downstream, start = FALSE),
            npoints = right_npoints,
            ignore_strand = ignore_strand
          )

  middle <- intersect_bw_and_granges(bw,
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
intersect_bw_and_granges <- function(bw,
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


#' Remove top loci in a GRanges column by mean of specified values
#'
#' @param granges GRanges object. Must have numerical mcols, and names of those
#'   needs to match with columns.
#' @param quantile Top quantile. Must be a number between zero and 1.
#' @param columns Which columns to use. If more than one, quantile will be
#'   selected according to means.
#'
#' @return A named list, fields ranges for the resulting GRanges, plus calculated
#'   values: quantile_value, filtered, na_values.
#' @export
remove_top_by_mean <- function(granges, quantile, columns) {
  # !names(valid_columns) %in% c("name")
  n_na <- 0
  n_filtered <- 0
  top_quantile <- NULL

  if (quantile > 0) {
    if (ncol(mcols(granges)) > 1) {
      valid_columns <- data.frame(mcols(granges))
      valid_columns <- valid_columns[, columns]
      means <- rowMeans(valid_columns)
      top_quantile <- quantile(means, probs = c(1-quantile), na.rm=TRUE)
      granges$means <- means

      n_na <- sum(is.na(granges$means))
      granges <- granges[!is.na(granges$means), ]

      n_filtered <- length(granges[granges$means > top_quantile, ])
      granges <- granges[granges$means <= top_quantile, ]
      granges$means <- NULL
    }
    else {
      top_quantile <- quantile(mcols(granges)[, 1],
                               probs = c(1-quantile), na.rm = TRUE)

      n_na <- sum(is.na(mcols(granges)[, 1]))
      granges <- granges[!is.na(mcols(granges)[, 1]), ]

      n_filtered <- length(granges[mcols(granges)[, 1] <= top_quantile, ])
      granges <- granges[mcols(granges)[, 1] <= top_quantile, ]
    }
  }
  list(ranges=granges,
       calculated=list(na=n_na, filtered=n_filtered,
                       quantile=unname(top_quantile))
      )
}

#' Summarize a intersect_bw_and_granges matrix
#'
#' Compute averages and standard error values for a matrix returned by
#' intersect_bw_and_granges.
#'
#' @param matrix A matrix returned by intersect_bw_and_granges.
#' @param label Label for the sample.
#' @importFrom stats median sd
#' @return A data frame with summarized values, stderr and medians, plus label.
summarize_matrix <- function(matrix, label) {
  # Ignore Inf and NaN in the computation of means/SD
  matrix[is.infinite(matrix)] <- NA

  omitted_vals <- sum(is.na(matrix))
  if (omitted_vals > 100) {
    mean_per_locus <- omitted_vals / nrow(matrix)
    warning(paste("Profile plot:",
                  omitted_vals, "generated (",
                  mean_per_locus, "per locus)"
            )
    )
  }

  df <- data.frame(
          mean = colMeans(matrix, na.rm = TRUE),
          sderror = apply(matrix, 2,
            function(n) {
              sd(n, na.rm = TRUE) / sqrt(sum(!is.na(n)))
            }),
          median = apply(matrix, 2, median, na.rm = TRUE)
        )

  df$index <- as.integer(rownames(df))
  df$sample <- label
  df
}


#' Make a BigWigFile object out of a path or an URL
#'
#' @param bw BigWig file path or URL
#'
#' @return BigWigFile object
#' @export
fetch_bigwig <- function(bw) {
  if (!is.null(bw)) {
    valid_bwfile <- bw
    if (RCurl::url.exists(bw)) {
      valid_bwfile <- tempfile()
      download.file(bw, valid_bwfile)
    }
    BigWigFile(path = valid_bwfile)
  }
}

#' GRanges cbind-like operation
#'
#' Perform a cbind operation on a GRanges list, appending scores to mcols.
#' It will sort the GRanges elements in order to ensure the match is proper.
#'
#' This assumes that you are trying to cbind things that match (bins generated
#' from the same parameters, BED intersections from the same BED.)
#'
#' @param grlist A list of GRanges objects that have all the same fields.
#' @param labels Vector of names for the score columns.
#' @importFrom GenomeInfoDb sortSeqlevels
granges_cbind <- function(grlist, labels) {
  fixed_fields <- c("seqnames", "start", "end", "width", "strand")

  grlist[[1]] <- sortSeqlevels(grlist[[1]])
  grlist[[1]] <- sort(grlist[[1]])

  result <- data.frame(grlist[[1]])[, fixed_fields]
  for (i in seq(1, length(grlist))) {
    grlist[[i]] <- sortSeqlevels(grlist[[i]])
    grlist[[i]] <- sort(grlist[[i]])

    result[, labels[[i]]] <- grlist[[i]]$score
  }

  result <- makeGRangesFromDataFrame(result, keep.extra.columns = TRUE)
  result
}


#' Processes a loci and returns a GRanges object, sorted and with its seqlevels
#' also sorted.
#'
#' @param loci Either a BED file or a GRanges object
#'
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb sortSeqlevels
loci_to_granges <- function(loci) {
  bed <- loci
  if (class(loci) == "character") {
    bed <- import(loci, format = "BED")
  }

  bed <- sortSeqlevels(bed)
  bed <- sort(bed, ignore.strand = FALSE)
  bed
}


#' Get a valid label from a filename
#'
#' @param filename File to convert to label
#'
#' @return A valid label name
make_label_from_filename <- function(filename) {
  filename_clean <- basename(tools::file_path_sans_ext(filename))
  make.names(filename_clean)
}


#' Take a column of a data frame to use as row names, drop such column and
#' reorder the data frame so row names follow natural order.
#'
#' @param df A data frame.
#' @param col The column (usually name).
#' @importFrom stringr str_sort
#' @return A sorted df
natural_sort_by_field <- function(df, col) {
  rownames(df) <- df[, col]
  order <- str_sort(df[, col], numeric = TRUE)
  df[, col] <- NULL
  df[order, , drop = FALSE]
}


#' Validate a category array
#'
#' Checks whether the number of categories is reasonable for an aggregation.
#' Throws a warning if it finds more than 50 different values.
#'
#' @param cat_values An array of values
validate_categories <- function(cat_values) {
  max_categories <- 50
  # Test number of values in group_col
  ncat <- length(levels(as.factor(cat_values)))
  if (ncat > max_categories) {
    warning(paste(
      "Number of values in group column field very large:", ncat,
      "(does BED file have unique IDs instead of categories?)")
    )
  }
}


#' Validate an array of paths
#'
#' Check that a list of files is valid: not empty and contents exist.
#'
#' @param filelist An array of files
#' @importFrom RCurl url.exists
#' @return NULL
validate_filelist <- function(filelist) {
  if (length(filelist) == 0) {
    stop("File list provided is empty.")
  }

  existence_flag <- file.exists(filelist) | RCurl::url.exists(filelist)
  if (! all(existence_flag)) {
    msg <- paste("Files not found:", filelist[!existence_flag])
    stop(msg)
  }
}

#' Validate that a locus parameter is valid. Checks for paths and also whether
#' the parameter is otherwise a GRanges object, which is also valid.
#'
#' @param locus_param Parameter to validate
#'
validate_locus_parameter <- function(locus_param) {
  if (class(locus_param) == "character") {
    validate_filelist(locus_param)
  }
  else {
    if (class(locus_param) != "GRanges") {
      msg <- paste0("Unexpected type: ", class(locus_param))
      stop(msg)
    }
  }
}


#' Validate that group col exists in granges
#'
#' @param granges GRanges object to check
#' @param group_col Group column name. Usually, name.
validate_group_col <- function(granges, group_col) {
  if (!group_col %in% names(mcols(granges))) {
    stop(paste("Invalid group column not present in granges", group_col))
  }
}


#' Validate profile and heatmap relevant parameters
#'
#' @param bin_size Bin size. Must be a positive number.
#' @param upstream Upstream bp. Must be positive and larger than bin size.
#' @param downstream Downstream bp. Must be positive and larger than bin size.
#'
validate_profile_parameters <- function(bin_size, upstream, downstream) {
  if (bin_size <= 0) {
    stop(paste("bin size must be a positive value:", bin_size))
  }

  if (upstream <= 0) {
    stop(paste("upstream size must be a positive value:", upstream))
  }

  if (downstream <= 0) {
    stop(paste("downstream size must be a positive value:", downstream))
  }

  if (bin_size > upstream || bin_size > downstream) {
    stop("bin size must be smaller than flanking regions")
  }
}
