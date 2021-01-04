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
#' @param norm_func_x Function to use after x / x_bg.
#' @param norm_func_y Function to use after y / y_bg.
#' @param highlight List of bed files to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted
#' @param highlight_label Labels for the highlight groups.
#'  If not provided, filenames are used.
#' @param highlight_colors Array of color values for the highlighting groups
#' @param verbose Put a caption with relevant parameters on the plot.
#' @import ggplot2
#' @inheritParams bw_bins
#' @return A ggplot object.
#' @export
plot_bw_bins_scatter <- function(x, y,
                                 bg_x = NULL, bg_y = NULL,
                                 norm_func_x = identity,
                                 norm_func_y = identity,
                                 bin_size = 10000,
                                 genome = "mm9",
                                 highlight = NULL,
                                 minoverlap = 0L,
                                 highlight_label = NULL,
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

  highlight_data <- process_highlight_loci(highlight, highlight_label)

  x_label <- paste(make_label_from_filename(x), "-",
                   make_norm_label(substitute(norm_func_x), bg_x))

  y_label <- paste(make_label_from_filename(y), "-",
                   make_norm_label(substitute(norm_func_y), bg_y))

  plot_results <- plot_ranges_scatter(bins_values_x, bins_values_y,
         highlight = highlight_data$ranges,
         minoverlap = minoverlap,
         highlight_label = highlight_data$labels,
         highlight_colors = highlight_colors,
         remove_top = remove_top
      )

  plot <- plot_results$plot + ggtitle(paste("Genome-wide bin coverage (", bin_size, "bp)", sep = "")) +
    xlab(x_label) + ylab(y_label) + default_theme()

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
    plot <- plot + labs(caption=verbose_tag)
  }

  plot
}

#' Locus-based scatterplot of a pair of bigWig files
#'
#' Plots a scatter plot from two given bigWig files, a locus object and a
#' set of BED files as highlighted annotations.
#' Loci are highlighted if there is at least minoverlap base pairs overlap
#' with any loci in BED file.
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
#' @param loci Bed file or GRanges object to be plotted.
#' @param bg_x BigWig file to be used as x axis background (us. input).
#' @param bg_y BigWig file to be used as y axis background (us. input).
#' @param norm_func_x Function to use after x / x_bg.
#' @param norm_func_y Function to use after y / y_bg.
#' @param highlight List of bed files to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted
#' @param highlight_label Labels for the highlight groups.
#'  If not provided, filenames are used.
#' @param highlight_colors Array of color values for the highlighting groups
#' @param remove_top Return range 0-(1-remove_top). By default returns the
#'     whole distribution (remove_top == 0).
#' @param verbose Verbose plot. Returns a plot with all relevant parameters in
#'   a caption.
#' @import ggplot2
#' @return A ggplot object.
#' @export
plot_bw_loci_scatter <- function(x, y,
                                 loci,
                                 bg_x = NULL, bg_y = NULL,
                                 norm_func_x = identity,
                                 norm_func_y = identity,
                                 highlight = NULL,
                                 minoverlap = 0L,
                                 highlight_label = NULL,
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

  highlight_data <- process_highlight_loci(highlight, highlight_label)

  x_label <- paste(make_label_from_filename(x), "-",
                   make_norm_label(substitute(norm_func_x), bg_x))

  y_label <- paste(make_label_from_filename(y), "-",
                   make_norm_label(substitute(norm_func_y), bg_y))

  plot_results <- plot_ranges_scatter(values_x, values_y,
                                      highlight = highlight_data$ranges,
                                      minoverlap = minoverlap,
                                      highlight_label = highlight_data$labels,
                                      highlight_colors = highlight_colors,
                                      remove_top = remove_top
  )

  # Show parameters and relevant values
  loci_name <- "GRanges object"
  if (class(loci) == "character") {
    loci_name <- basename(loci)
  }

  plot <- plot_results$plot + ggtitle(paste("Per-locus coverage (", loci_name, ")", sep = "")) +
    xlab(x_label) + ylab(y_label) + default_theme()

  if (verbose) {
    relevant_params <- list(loci=loci_name,
                            minoverlap=minoverlap,
                            remove_top=remove_top)

    crop_values <- list(points=length(values_x),
                        removed=plot_results$removed,
                        NAs=plot_results$na_points)

    verbose_tag <- make_caption(relevant_params, crop_values)
    plot <- plot + labs(caption=verbose_tag)
  }

  plot
}


#' Summary heatmap of a categorized BED or GRanges object
#'
#' Make a summary heatmap where each cell contains an aggregated value of a
#' bigWig file from bwfiles and a category of a BED file or GRanges (loci). The
#' provided loci must have a name field that is valid (i.e. can be grouped,
#' representing some type of category).
#'
#' @param labels Labels to use for in the plot for the bw files.
#' @param loci BED file or GRanges object.
#' @param verbose Put a caption with relevant parameters on the plot.
#' @inheritParams bw_bed
#' @return A ggplot object
#' @export
plot_bw_loci_summary_heatmap <- function(bwfiles,
                                         loci,
                                         bg_bwfiles = NULL,
                                         labels = NULL,
                                         aggregate_by = "true_mean",
                                         norm_func = identity,
                                         remove_top = 0,
                                         verbose = TRUE) {

  summary_values <- bw_bed(bwfiles, loci,
                           bg_bwfiles = bg_bwfiles,
                           aggregate_by = aggregate_by,
                           norm_func = norm_func,
                           labels = labels,
                           remove_top = remove_top
  )


  if (sum(summary_values) == 0) {
    warning("All zero-values matrix. Using same background as bw input?")
  }



  legend_label = make_norm_label(substitute(norm_func), bg_bwfiles)
  colorscale <- scale_fill_gradient(name=legend_label, low = "white", high="#B22222")
  if (!is.null(bg_bwfiles)) {
    colorscale <- scale_fill_gradient2(name=legend_label, low = "#2e6694", mid="white", high="#B22222")
  }

  summary_values$type <- rownames(summary_values)
  vals_melted <- reshape2::melt(summary_values, id.vars="type")

  title <- paste("Coverage per region (", aggregate_by, ")")
  title <- paste(title, "-", make_norm_label(substitute(norm_func), bg_bwfiles))

  plot <- ggplot(vals_melted, aes_string("type", "variable", fill="value")) +
    geom_tile() + geom_text(aes(label=round(value, 2)), size=3.5) +
    coord_fixed() + colorscale +
    scale_y_discrete(position="right") + theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          legend.position=c(1,1.2),
          legend.direction="horizontal",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + xlab("") + ylab("") +
    ggtitle(title)

  if (verbose) {
    # Show parameters and relevant values
    relevant_params <- list(aggregate_by=aggregate_by,
                            remove_top=remove_top)

    verbose_tag <- make_caption(relevant_params, list())
    plot <- plot + labs(caption=verbose_tag)
  }

  plot
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
#' @param remove_top Return range 0-(1-remove_top). By default returns the
#'     whole distribution (remove_top == 0).
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
    geom_point(color = "#bbbbbb", alpha = 0.8) +
    xlab(x_label) +
    ylab(y_label) +
    extra_plot +
    extra_colors

  list(plot=p, na_points=n_removed, removed=n_filtered)
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


#' Internal processing of loci sets and labels.
#'
#' This allows to process equally lists of GRanges objects or BED files. Note
#' that GRanges objects need to be lists, as c(GRanges) concatenates elements
#' in single GRanges object
#'
#' @param loci_sets One or more BED files or GRanges. These need to be either
#'   all BED files or all GRanges objects.
#' @param labels One or more labels. Lengths need to match loci_sets or be NULL,
#'   in which case loci_sets are assumed to be BED files.

#' @return a named list with ranges and labels.
process_highlight_loci <- function(loci_sets, labels) {
  gr_list <- loci_sets
  lab_list <- labels
  if (!is.null(gr_list)) {
    if (class(gr_list) == "character") {
      # list is used instead of c() because GRanges c() concatenates ranges
      gr_list <- lapply(loci_sets, rtracklayer::import, format = "BED")
      if (is.null(labels)) {
        lab_list <- basename(loci_sets)
      }
    }
    else {
      if (is.null(labels)) {
        stop("GRanges used as highlight loci but no labels provided")
      }
      if (class(gr_list) != "list") {
        # If single GRanges was passed it needs to be converted to list
        gr_list <- list(gr_list)
      }
    }
  }
  if (length(gr_list) != length(lab_list)) {
    stop("Highlight loci sets don't match the number of labels provided")
  }

  list(ranges=gr_list, labels=lab_list)
}
