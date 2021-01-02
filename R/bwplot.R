#' Bin-based scatterplot of a pair of bigWig files
#'
#' Plots a scatter plot from two given bigWig files and an optional set of BED
#' files as highlighted annotations. Bins are highlighted if there is at least
#' minoverlap base pairs overlap with any loci in BED file.
#'
#' If specifying minoverlap, you must take into account the bin_size parameter
#' and the size of the loci you are providing as BED file.
#'
#' Values in x and y axis can be normalized using background bigWig files
#' (usually input files). By default, the value shown will be x / bg_x per bin.
#' If norm_func_x or norm_func_y are provided, this can be changed to any given
#' function, for instance, if norm_func_x = log2, values on the x axis will
#' represent log2(x / bg_x) for each bin.
#'
#' Values that are invalid (NaN, Inf, -Inf) in doing such normalization will
#' be ignored and shown as warnings, as this is ggplot default behavior.
#'
#' @param x BigWig file corresponding to the x axis.
#' @param y BigWig file corresponding to the y axis.
#' @param bg_x BigWig file to be used as x axis background (us. input).
#' @param bg_y BigWig file to be used as y axis background (us. input).
#' @param highlight List of bed files to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted
#' @param highlight_label Labels for the highlight groups.
#'  If not provided, filenames are used.
#' @param norm_func_x Function to use after x / x_bg.
#' @param norm_func_y Function to use after y / y_bg.
#' @param highlight_colors Array of color values for the highlighting groups
#' @import ggplot2
#' @inheritParams bw_bins
#' @return A ggplot object.
#' @export
plot_bw_bins_scatter <- function(x, y,
                                 bg_x = NULL, bg_y = NULL,
                                 bin_size = 10000,
                                 genome = "mm9",
                                 highlight = NULL,
                                 minoverlap = 0L,
                                 highlight_label = NULL,
                                 norm_func_x = identity,
                                 norm_func_y = identity,
                                 highlight_colors = NULL,
                                 remove_top = 0,
                                 verbose = TRUE) {

  bins_values_x <- bw_bins(x,
                           bg_bwfiles = bg_x,
                           bin_size = bin_size,
                           genome = genome,
                           norm_func = norm_func_x,
                           labels = "x"
  )

  bins_values_y <- bw_bins(y,
                           bg_bwfiles = bg_y,
                           bin_size = bin_size,
                           genome = genome,
                           norm_func = norm_func_y,
                           labels = "y"
  )

  if (!is.null(highlight)) {
    gr_list <- lapply(highlight, rtracklayer::import, format = "BED")
    if (is.null(highlight_label)) {
      highlight_label <- basename(highlight)
    }
  }

  x_label <- paste(make_label_from_filename(x), "-",
                   make_norm_label(substitute(norm_func_x), bg_x))

  y_label <- paste(make_label_from_filename(y), "-",
                   make_norm_label(substitute(norm_func_y), bg_y))

  plot_results <- plot_ranges_scatter(bins_values_x, bins_values_y,
         highlight = gr_list,
         minoverlap = minoverlap,
         highlight_label = highlight_label,
         highlight_colors = highlight_colors,
         remove_top = remove_top
      )

  verbose_tag <- ""
  if (verbose) {
    # Show parameters and relevant values
    relevant_params <- list(genome=genome,
                            bin_size=bin_size,
                            minoverlap=minoverlap,
                            remove_top=remove_top)

    crop_values <- list(points=length(bins_values_x),
                        removed=plot_results$removed,
                        NAs=plot_results$na_points)

    verbose_tag <- make_caption(relevant_params, crop_values)
  }

  plot_results$plot + ggtitle(paste("Genome-wide bin coverage (", bin_size, "bp)", sep = "")) +
    xlab(x_label) + ylab(y_label) + default_theme() + labs(caption=verbose_tag)
}

#' Bin-based scatterplot of a pair of bigWig files
#'
#' Plots a scatter plot from two given bigWig files and an optional set of BED
#' files as highlighted annotations. Bins are highlighted if there is at least
#' minoverlap base pairs overlap with any loci in BED file.
#'
#' If specifying minoverlap, you must take into account the bin_size parameter
#' and the size of the loci you are providing as BED file.
#'
#' Values in x and y axis can be normalized using background bigWig files
#' (usually input files). By default, the value shown will be x / bg_x per bin.
#' If norm_func_x or norm_func_y are provided, this can be changed to any given
#' function, for instance, if norm_func_x = log2, values on the x axis will
#' represent log2(x / bg_x) for each bin.
#'
#' Values that are invalid (NaN, Inf, -Inf) in doing such normalization will
#' be ignored and shown as warnings, as this is ggplot default behavior.
#'
#' @param x BigWig file corresponding to the x axis.
#' @param y BigWig file corresponding to the y axis.
#' @param bg_x BigWig file to be used as x axis background (us. input).
#' @param bg_y BigWig file to be used as y axis background (us. input).
#' @param highlight List of bed files to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted
#' @param highlight_label Labels for the highlight groups.
#'  If not provided, filenames are used.
#' @param norm_func_x Function to use after x / x_bg.
#' @param norm_func_y Function to use after y / y_bg.
#' @param highlight_colors Array of color values for the highlighting groups
#' @import ggplot2
#' @inheritParams bw_bins
#' @return A ggplot object.
#' @export
plot_bw_loci_scatter <- function(x, y,
                                 loci,
                                 bg_x = NULL, bg_y = NULL,
                                 highlight = NULL,
                                 minoverlap = 0L,
                                 highlight_label = NULL,
                                 norm_func_x = identity,
                                 norm_func_y = identity,
                                 highlight_colors = NULL,
                                 remove_top = 0,
                                 verbose = TRUE) {

  values_x <- bw_bed(x, bg_bwfiles = bg_x, bedfile = loci,
                     norm_func = norm_func_x,
                     labels = "x"
              )

  values_y <- bw_bed(y, bg_bwfiles = bg_y, bedfile = loci,
                     norm_func = norm_func_y,
                     labels = "y"
              )

  gr_list <- NULL
  if (!is.null(highlight)) {
    gr_list <- lapply(highlight, rtracklayer::import, format = "BED")
    if (is.null(highlight_label)) {
      highlight_label <- basename(highlight)
    }
  }

  x_label <- paste(make_label_from_filename(x), "-",
                   make_norm_label(substitute(norm_func_x), bg_x))

  y_label <- paste(make_label_from_filename(y), "-",
                   make_norm_label(substitute(norm_func_y), bg_y))

  plot_results <- plot_ranges_scatter(values_x, values_y,
                                      highlight = gr_list,
                                      minoverlap = minoverlap,
                                      highlight_label = highlight_label,
                                      highlight_colors = highlight_colors,
                                      remove_top = remove_top
  )

  verbose_tag <- ""
  # Show parameters and relevant values
  loci_name <- "GRanges object"
  if (class(loci) == "character") {
    loci_name <- basename(loci)
  }

  if (verbose) {
    relevant_params <- list(loci=loci_name,
                            minoverlap=minoverlap,
                            remove_top=remove_top)

    crop_values <- list(points=length(values_x),
                        removed=plot_results$removed,
                        NAs=plot_results$na_points)

    verbose_tag <- make_caption(relevant_params, crop_values)
  }

  plot_results$plot + ggtitle(paste("Per-locus coverage (", loci_name, ")", sep = "")) +
    xlab(x_label) + ylab(y_label) + default_theme() + labs(caption=verbose_tag)
}



#' Scatterplot of values in GRanges objects. Loci must match.
#'
#' Plots a scatter plot from two given GRanges objects and an optional set of
#' GRanges as highlighted annotations.
#'
#' Values in x and y axis can be normalized using background.
#'
#' Values that are invalid (NaN, Inf, -Inf) in doing such normalization will
#' be ignored and shown as warnings, as this is ggplot default behavior.
#'
#' @param x GRanges x axis
#' @param y GRanges y axis
#' @param highlight List of GRanges to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a locus to be highlighted
#' @param highlight_label Labels for the highlight groups.
#'  If not provided, column names are used.
#' @param highlight_colors Array of color values for the highlighting groups
#' @import ggplot2
#' @return A ggplot object.
plot_ranges_scatter <- function(x, y,
                                highlight = NULL,
                                minoverlap = 0L,
                                highlight_label = NULL,
                                highlight_colors = NULL,
                                remove_top = 0) {


  label_df <- function(df, name) {
    data.frame(df, group = name)
  }

  values_x <- data.frame(x)
  values_y <- data.frame(y)

  df <- merge(values_x, values_y, by=c("seqnames", "start", "end", "width", "strand"))

  n_removed <- 0
  n_filtered <- 0

  if (remove_top > 0) {
    column_names <- c("x", "y")
    df$means <- rowMeans(df[, column_names])
    df_clean <- df[!is.na(df$means), ]
    n_removed <- nrow(df) - nrow(df_clean)

    top_quantile = quantile(df_clean$means, probs=c(1-remove_top))
    df <- df_clean[df_clean$means <= top_quantile, ]

    n_filtered <- nrow(df_clean) - nrow(df)
  }

  ranges_values <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)

  extra_plot <- NULL
  extra_colors <- NULL

  if (!is.null(highlight)) {
    subset_func <- purrr::partial(
      IRanges::subsetByOverlaps,
      x = ranges_values,
      minoverlap = minoverlap
    )

    ranges_subset <- lapply(highlight, subset_func)
    subset_df <- lapply(ranges_subset, data.frame)

    df_values_labeled <- purrr::map2(subset_df, highlight_label, label_df)
    highlight_values <- do.call(rbind, df_values_labeled)

    # Order of factors need to match to assign properly colors to points
    highlight_values$group <- factor(highlight_values$group,
                                     levels = highlight_label)

    extra_plot <- geom_point(
      data = highlight_values,
      aes_string(x = "x", y = "y", color = "group"),
      alpha = 0.8
    )

    if (!is.null(highlight_colors)) {
      extra_colors <- scale_color_manual(values=highlight_colors)
    }
  }

  x_label <- "x"
  y_label <- "y"

  p <- ggplot(df, aes_string(x = "x", y = "y")) +
    geom_point(color = "#cccccc", alpha = 0.8) +
    xlab(x_label) +
    ylab(y_label) +
    extra_plot +
    extra_colors

  list(plot=p, na_points=n_removed, removed=n_filtered)
}


#' Generate a human-readable normalization function string
#'
#' @param f String representing normalization function.
#' @param bg Background file.
#'
#' @return A string describing normalization.
make_norm_label <- function(f, bg) {
  label <- "RPGC"
  if (!is.null(bg)) {
    if (f != "identity") {
      label <- paste(f, "(", label, " / background)", sep = "")
    } else {
      label <- paste(label, " / background", sep = "")
    }
  }
  label
}


#' Set default theme as classic with larger font size
default_theme <- function() {
  ggplot2::theme_classic(base_size = 18)
}

#' Make a string to put as capsion in verbose mode. Includes system date.
#'
#' @param params Named list with relevant parameters and their values
#' @param outcome Named values with relevant outcomes and their values
#' @return A caption string
make_caption <- function(params, outcome) {
  verbose_params <- paste(names(params),
                          params, sep = ":", collapse = ", ")

  verbose_crop <- paste(names(outcome),
                        outcome, sep = ":", collapse = ", ")

  date <- format(Sys.time(), "%a %b %d %X %Y")
  paste(verbose_params, verbose_crop, date, sep="\n")
}

