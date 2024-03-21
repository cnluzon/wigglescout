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

# Main plot functions -----------------------------------------------------

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
#' @param norm_mode_x Normalization mode for x axis.
#' @param norm_mode_y Normalization mode for y axis.
#' @param highlight List of bed files to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted
#' @param highlight_label Labels for the highlight groups.
#'  If not provided, filenames are used.
#' @param highlight_colors Array of color values for the highlighting groups
#' @param verbose Put a caption with relevant parameters on the plot.
#' @param density Plot density tiles for global distribution instead of points.
#' @import ggplot2
#' @importFrom purrr partial
#' @inheritParams bw_bins
#' @return A ggplot object.
#' @examples
#' # Get the raw files
#' bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
#' bw2 <- system.file("extdata",
#'                    "sample_H3K9me3_ChIP.bw", package="wigglescout")
#' bed <- system.file("extdata", "sample_genes_mm9.bed", package="wigglescout")
#'
#' # Sample bigWig files only have valid values on this region
#' locus <- GenomicRanges::GRanges(
#'   seqnames = "chr15",
#'   IRanges::IRanges(102600000, 103100000)
#' )
#'
#' plot_bw_bins_scatter(bw, bw2, bin_size = 50000, selection = locus)
#' plot_bw_bins_scatter(bw, bw2, bin_size = 50000, selection = locus,
#'     remove_top = 0.01)
#' plot_bw_bins_scatter(bw, bw2, bin_size = 50000, selection = locus,
#'     highlight = bed, highlight_color = "black")
#'
#' @export
plot_bw_bins_scatter <- function(x, y,
                                bg_x = NULL, bg_y = NULL,
                                norm_mode_x = "fc", norm_mode_y = "fc",
                                bin_size = 10000,
                                genome = "mm9",
                                highlight = NULL,
                                minoverlap = 0L,
                                highlight_label = NULL,
                                highlight_colors = NULL,
                                remove_top = 0,
                                verbose = TRUE,
                                density = FALSE,
                                selection = NULL,
                                default_na = NA_real_,
                                scaling = "none") {

    partial_bw_bins <- partial(
        bw_bins, bin_size = bin_size, genome = genome,
        selection = selection, labels = "score", scaling = scaling,
        default_na = default_na
    )
    bins_x <- partial_bw_bins(
        bwfiles = x, bg_bwfiles = bg_x, norm_mode = norm_mode_x
    )
    bins_y <- partial_bw_bins(
        bwfiles = y, bg_bwfiles = bg_y, norm_mode = norm_mode_y
    )
    highlight_data <- .convert_and_label_loci(highlight, highlight_label)
    clean_gr <- .filter_scatter_data(bins_x, bins_y, remove_top)

    main_plot <- .scatterplot_body(clean_gr$ranges,
        highlight = highlight_data$ranges,
        highlight_label = highlight_data$labels,
        highlight_colors = highlight_colors,
        minoverlap = minoverlap,
        density = density
    )

    title <- paste("Genome-wide bin coverage (", bin_size, "bp)", sep = "")
    x_lab <- .make_norm_file_label(norm_mode_x, x, bg_x)
    y_lab <- .make_norm_file_label(norm_mode_y, y, bg_y)
    params <- mget(c("genome", "bin_size", "minoverlap", "remove_top", "scaling"))
    caption <- .make_caption(params, clean_gr$stats, verbose = verbose)

    labels <- labs(title = title, x = x_lab, y = y_lab, caption = caption)
    main_plot + labels + .theme_default()
}

#' Bin-based violin plot of a set of bigWig files
#'
#' Plots a violin plot of bin distribution of a set of bigWig files optionally
#' overlaid with annotated bins. Bins overlapping loci of the provided BED
#' file will be shown as a jitter plot on top of the violin plot.
#'
#' @param highlight BED file to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted.
#' @param highlight_label Label for the highlighted loci set
#' @param highlight_colors Array of color values for the highlighted groups.
#' @param verbose Put a caption with relevant parameters on the plot.
#' @inheritParams bw_bins
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @return A ggplot object.
#' @examples
#' # Get the raw files
#' bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
#' bw2 <- system.file("extdata", "sample_H3K9me3_ChIP.bw",
#'                    package="wigglescout")
#'
#' bed <- system.file("extdata", "sample_genes_mm9.bed", package="wigglescout")
#'
#' # Sample bigWig files only have valid values on this region
#' locus <- GenomicRanges::GRanges(
#'   seqnames = "chr15",
#'   IRanges::IRanges(102600000, 103100000)
#' )
#'
#' plot_bw_bins_violin(c(bw, bw2), bin_size = 50000, selection = locus)
#'
#' # You need to provide as many color values as bigWig files when highlighting
#' plot_bw_bins_violin(c(bw, bw2), bin_size = 50000, selection = locus,
#'   highlight = bed, highlight_label = "Genes",
#'   highlight_color = c("black", "black"))
#'
#' @export
plot_bw_bins_violin <- function(bwfiles,
                                bg_bwfiles = NULL,
                                labels = NULL,
                                bin_size = 10000,
                                per_locus_stat = "mean",
                                genome = "mm9",
                                highlight = NULL,
                                minoverlap = 0L,
                                norm_mode = "fc",
                                highlight_label = NULL,
                                highlight_colors = NULL,
                                remove_top = 0,
                                verbose = TRUE,
                                selection = NULL,
                                default_na = NA_real_,
                                scaling = "none") {

    # Get parameter values that are relevant to the underlying function
    par <- .get_wrapper_parameter_values(bw_bins, mget(names(formals())))
    values <- do.call(bw_bins, mget(par))

    columns <- labels
    if (is.null(labels))
        columns <- .make_label_from_object(bwfiles)

    clean_gr <- .filter_violin_data(values, remove_top, columns)

    main_plot <- .violin_body(clean_gr$ranges,
        highlight = highlight,
        minoverlap = minoverlap,
        highlight_label = highlight_label,
        highlight_colors = highlight_colors
    )

    title <- paste("Genome-wide bin distribution (", bin_size, "bp)", sep = "")
    y_label <- .make_norm_label(norm_mode, bg_bwfiles)
    params <- mget(c("genome", "bin_size", "minoverlap", "remove_top", "scaling"))
    caption <- .make_caption(params, clean_gr$stats, verbose = verbose)
    labels <- labs(title = title, x = "", y = y_label, caption = caption)

    rotate_x <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
    no_legend <- theme(legend.position = "none")

    main_plot + labels + .theme_default() + rotate_x + no_legend
}

#' Plot a heatmap of a given bigWig file over a set of loci
#'
#' @param bwfiles BigWig file to plot
#' @param bg_bwfile Background bw file. Use this with care. Depending on bin
#'   size and actual values, this may result in a very noisy plot.
#' @param zmin Minimum of the color scale. Majority of tools set
#'   this to 0.01 percentile of the data distribution.
#' @param zmax Maximum of the color scale. Majority of tools set
#'   this to 0.99.
#' @param cmap Color map. Any RColorBrewer palette name is accepted here.
#' @param max_rows_allowed Maximum number of loci that will be allowed in the
#'   plot. If the amount of loci exceeds this value, the plot will be binned
#'   on the y axis until it fits max_rows_allowed. This speeds up plotting of
#'   very large matrices, where higher resolution would not be perceivable by
#'   eye.
#' @param verbose Put a caption with relevant parameters on the plot.
#' @param order_by Specific order to display rows. By default rows are sorted
#'   decreasingly by mean. Order should be an array of integers of the same
#'   length as number of rows. These can be obtained as a result of order()
#'   function, if one would want to sort one heatmap by values on another one.
#' @importFrom dplyr group_by summarise
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @inheritParams plot_bw_profile
#' @return A ggplot object
#' @examples
#' # Get the raw files
#' bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
#' bed <- system.file("extdata", "sample_genes_mm9.bed", package="wigglescout")
#'
#'
#' plot_bw_heatmap(bw, loci = bed,
#'                 mode = "center", upstream = 1000, downstream = 1500)
#'
#' # If resolution is lower than number of rows, rows are aggregated
#' plot_bw_heatmap(bw, loci = bed,
#'                 mode = "center", upstream = 1000, downstream = 1500,
#'                 max_rows_allowed = 14)
#' @export
plot_bw_heatmap <- function(bwfiles, loci,
                            bg_bwfile = NULL,
                            mode = "stretch",
                            bin_size = 100,
                            upstream = 2500,
                            downstream = 2500,
                            middle = NULL,
                            ignore_strand = FALSE,
                            norm_mode = "fc",
                            cmap = "Reds",
                            zmin = NULL,
                            zmax = NULL,
                            max_rows_allowed = 10000,
                            order_by = NULL,
                            verbose = TRUE,
                            default_na = NA_real_,
                            scaling = "none") {
    # Get parameter values that are relevant to the underlying function
    par <- .get_wrapper_parameter_values(bw_heatmap, mget(names(formals())))

    values <- do.call(bw_heatmap, mget(par))
    df <- .preprocess_heatmap_matrix(values[[1]],
        zmin, zmax,
        order_by,
        max_rows_allowed
    )

    main_plot <- .heatmap_body(df$values, df$stats$zmin, df$stats$zmax, cmap)
    params <- mget(c("mode", "bin_size", "middle", "ignore_strand",
                    "max_rows_allowed", "scaling"))

    if (! is.null(order_by)) {
      params$order_by = "user_defined"
    } else {
      params$order_by = "default"
    }

    caption <- .make_caption(params, df$stats, verbose = verbose)

    nloci <- nrow(values[[1]])
    y_label <- paste(.make_label_from_object(loci), "-", nloci, "loci")
    x_title <- .make_label_from_object(bwfiles)
    thin_rect <- element_rect(color = "black", fill = NA, linewidth = 0.1)
    thin_lines <- theme(axis.line = element_blank(), panel.border = thin_rect)
    plot_lines <- .heatmap_lines(nloci, ncol(values[[1]]), bin_size,
                                 upstream, downstream, mode)
    y_scale_pos <- nloci
    if (nloci > max_rows_allowed) {
        y_scale_pos <- max_rows_allowed
    }
    y_scale <- scale_y_continuous(breaks = c(1, y_scale_pos),
                                  labels = c(nloci, "1"),
                                  expand = c(0, 0))
    labels <- labs(title = "Heatmap", x = x_title, y = y_label,
                   caption = caption,
                   fill = .make_norm_label(norm_mode, bg_bwfile))
    main_plot + .theme_default() + thin_lines + plot_lines + y_scale + labels
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
#' @param norm_mode_x Normalization mode for x axis.
#' @param norm_mode_y Normalization mode for y axis.
#' @param bg_x BigWig file to be used as x axis background (us. input).
#' @param bg_y BigWig file to be used as y axis background (us. input).
#' @param highlight List of bed files to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted
#' @param highlight_label Labels for the highlight groups.
#'  If not provided, filenames are used.
#' @param highlight_colors Array of color values for the highlighting groups
#' @param remove_top Return range 0-(1-remove_top). By default returns the
#'     whole distribution (remove_top == 0).
#' @param default_na Default value for missing values
#' @param verbose Verbose plot. Returns a plot with all relevant parameters in
#'   a caption.
#' @import ggplot2
#' @importFrom methods is
#' @return A ggplot object.
#' @examples
#' # Get the raw files
#' bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
#' bw2 <- system.file("extdata", "sample_H3K9me3_ChIP.bw",
#'                    package="wigglescout")
#' bed <- system.file("extdata", "sample_genes_mm9.bed", package="wigglescout")
#'
#' plot_bw_loci_scatter(bw, bw2, loci = bed)
#'
#' # It is also possible to plot from a GRanges object
#' gr <- rtracklayer::import(bed)
#' plot_bw_loci_scatter(bw, bw2, loci = gr)
#' @export
plot_bw_loci_scatter <- function(x,
                                y,
                                loci,
                                bg_x = NULL,
                                bg_y = NULL,
                                norm_mode_x = "fc",
                                norm_mode_y = "fc",
                                highlight = NULL,
                                minoverlap = 0L,
                                highlight_label = NULL,
                                highlight_colors = NULL,
                                remove_top = 0,
                                verbose = TRUE,
                                default_na = NA_real_) {

    values_x <- bw_loci(x, bg_bwfiles = bg_x,
        loci = loci,
        norm_mode = norm_mode_x,
        labels = "score",
        default_na = default_na
    )

    values_y <- bw_loci(y, bg_bwfiles = bg_y,
        loci = loci,
        norm_mode = norm_mode_y,
        labels = "score",
        default_na = default_na
    )

    highlight_data <- .convert_and_label_loci(highlight, highlight_label)

    clean_gr <- .filter_scatter_data(values_x, values_y, remove_top)
    main_plot <- .scatterplot_body(clean_gr$ranges,
                                   highlight = highlight_data$ranges,
                                   minoverlap = minoverlap,
                                   highlight_label = highlight_data$labels,
                                   highlight_colors = highlight_colors,
    )

    loci_name <- .make_label_from_object(loci)
    params <- mget(c("loci_name", "minoverlap", "remove_top"))
    caption <- .make_caption(params, clean_gr$stats, verbose = verbose)

    title <- paste("Per-locus coverage (", loci_name, ")", sep = "")
    x_label <- .make_norm_file_label(norm_mode_x, x, bg_x)
    y_label <- .make_norm_file_label(norm_mode_y, y, bg_y)
    labels <- labs(title = title, x = x_label, y = y_label, caption = caption)

    main_plot + .theme_default() + labels
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
#' @inheritParams bw_loci
#' @return A ggplot object
#' @examples
#' # Get the raw files
#' bw <- system.file("extdata", "sample_H33_ChIP.bw", package = "wigglescout")
#' bw2 <- system.file("extdata",
#'                    "sample_H3K9me3_ChIP.bw", package = "wigglescout")
#' bed <- system.file("extdata", "sample_chromhmm.bed", package = "wigglescout")
#'
#' plot_bw_loci_summary_heatmap(c(bw, bw2), loci = bed,
#'                              labels = c("H33", "H3K9m3"))
#'
#' plot_bw_loci_summary_heatmap(c(bw, bw2), loci = bed, remove_top = 0.001,
#'                              labels = c("H33", "H3K9m3"))
#' @export
plot_bw_loci_summary_heatmap <- function(bwfiles, loci,
                                        bg_bwfiles = NULL,
                                        labels = NULL,
                                        aggregate_by = "true_mean",
                                        norm_mode = "fc",
                                        remove_top = 0,
                                        verbose = TRUE,
                                        default_na = NA_real_,
                                        scaling = "none") {

    # Get parameter values that are relevant to the underlying function
    par <- .get_wrapper_parameter_values(bw_loci, mget(names(formals())))
    summary_values <- do.call(bw_loci, mget(par))
    colorscale <- .colorscale(norm_mode, bg_bwfiles)
    plot <- .summary_body(summary_values)
    title <- paste0("Coverage per region (", aggregate_by, ")")
    params <- mget(c("aggregate_by", "remove_top", "scaling"))
    caption <- .make_caption(params, list(), verbose = verbose)
    labels <- labs(title = title, caption = caption, x = "", y = "")

    plot + colorscale + labels
}


#' Profile plot of a set of bigWig files
#'
#' Plots a profile of a set of bigWig files over a set of loci in a BED file.
#'
#' @param show_error Show standard error.
#' @param colors Array of colors that will  be assigned to labels or files
#'    (in that order)
#' @param verbose Put a caption with relevant parameters on the plot.
#' @param scaling Whether to use the bigWig values as they are (none - default)
#'    or calculate relative enrichment (relative) by dividing values by global
#'    average.
#' @inheritParams bw_profile
#' @import ggplot2
#' @importFrom methods is
#' @importFrom purrr partial map2
#' @return A ggplot object.
#' @examples
#' # Get the raw files
#' bw <- system.file("extdata", "sample_H33_ChIP.bw", package="wigglescout")
#' bed <- system.file("extdata", "sample_genes_mm9.bed", package="wigglescout")
#'
#'
#' plot_bw_profile(bw, loci = bed,
#'                 mode = "stretch", upstream = 1000, downstream = 1000)
#'
#' # It is also possible to run same bigWig across several BED files or GRanges
#' plot_bw_profile(bw, loci = c(bed, bed), labels = c("A", "B"))
#' @export
plot_bw_profile <- function(bwfiles, loci,
                            bg_bwfiles = NULL,
                            mode = "stretch",
                            bin_size = 100,
                            upstream = 2500,
                            downstream = 2500,
                            middle = NULL,
                            ignore_strand = FALSE,
                            show_error = FALSE,
                            norm_mode = "fc",
                            labels = NULL,
                            colors = NULL,
                            remove_top = 0,
                            verbose = TRUE,
                            default_na = NA_real_,
                            scaling = "none") {
    .validate_input_numbers(bwfiles, loci)
    values <- NULL
    nloci <- NULL
    x_label <- ""
    # Get parameter values that are relevant to the underlying function
    par <- .get_wrapper_parameter_values(bw_profile, mget(names(formals())))
    if (.is_multiple_loci(loci)) {
        labels <- .label_multiple_loci(loci, labels)

        profile_function <- partial(bw_profile, bwfile = bwfiles,
            bg_bwfiles = bg_bwfiles,
            mode = mode,
            bin_size = bin_size,
            upstream = upstream,
            downstream = downstream,
            middle = middle,
            ignore_strand = ignore_strand,
            norm_mode = norm_mode,
            remove_top = remove_top,
            default_na = default_na,
            scaling = scaling
        )
        value_list <- map2(loci, labels, profile_function)
        values <- do.call(rbind, value_list)
        x_tit <- "Multiple loci groups"
    }
    else {
        values <- do.call(bw_profile, mget(par))
        nloci <- .loci_length(loci)
        x_tit <- paste(
          .make_label_from_object(loci), "-", nloci, "loci",
          sep = " "
        )
    }
    y_lab <- .make_norm_label(norm_mode, bg_bwfiles)
    params <- mget(c("bin_size", "middle", "mode", "ignore_strand",
                     "remove_top", "scaling"))
    caption <- .make_caption(params, list(), verbose = verbose)
    if (!is.null(bg_bwfiles) && show_error == TRUE) {
        warning("Stderr estimate not available when normalizing by input")
        show_error <- FALSE
    }
    fig_labels <- labs(title = "Profile", x = x_tit, y = y_lab, caption = caption)
    .profile_body(values, show_error, colors, labels) +
        .heatmap_lines(
          nloci,
          max(values$index),
          bin_size,
          upstream,
          downstream,
          mode,
          expand = FALSE
        ) +
        fig_labels
}

# Helper plot functions ---------------------------------------------------

#' Helper function for matrix heatmap plot
#'
#' @param values Matrix with values
#' @inheritParams plot_bw_heatmap
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate row_number `%>%` summarise group_by
#' @importFrom ggplot2 ggplot scale_fill_gradientn geom_raster
#'
#' @return Named list plot and calculated values
.heatmap_body <- function(values, zmin, zmax, cmap) {
    gcol <- colorRampPalette(brewer.pal(n = 8, name = cmap))

    colorscale <- scale_fill_gradientn(colours = gcol(100),
        limits = c(zmin, zmax), breaks = c(zmin, zmax),
        labels = format(c(zmin, zmax), digits = 2),
        na.value = "#cccccc"
    )

    ggplot(values, aes(x = .data$x, y = .data$y, fill = .data$value)) +
        geom_raster() + colorscale
}

#' Helper function plots a profile from a dataframe
#'
#' @param values Values data frame in long format
#' @param colors Alternative colors to plot the lines
#' @param labels Labels order so label-i has colors-i.
#' @param show_error Boolean whether tho show error.
#'
#' @return A ggplot object
.profile_body <- function(values, show_error, colors, labels) {

    if (is.null(labels)) {
      labels <- unique(values$sample)
    }

    labels_short <- sapply(labels, .trunc_str)

    values <- values %>%
      mutate(
        min_error = .data$mean - .data$sderror,
        max_error = .data$mean + .data$sderror,
        sample = factor(.data$sample, levels = labels)
      )

    # The first step reorders, this one renames if needed
    levels(values$sample) <- labels_short

    p <- ggplot(
        values,
        aes(x = .data$index, y = .data$mean, color = .data$sample,  fill = .data$sample)
    ) + geom_line(linewidth = 1) + .theme_default() +
        theme(
            legend.position = "inside",
            legend.position.inside = c(0.8, 0.9),
            legend.direction = "vertical",
            legend.title = element_blank(),
            legend.background = element_rect(fill = alpha("white", 0.3))
        )

    if (!is.null(colors)) {
        p <- p +
            scale_color_manual(values = colors) +
            scale_fill_manual(values = colors)
    }

    if (show_error) {
        p <- p + geom_ribbon(aes(
            x = .data$index,
            ymin = .data$min_error,
            ymax = .data$max_error
        ),
        color = NA, alpha = 0.3
        )
    }

    p
}


#' Helper function for plotting a summary matrix
#' @param values Summary matrix
#' @return A ggplot object.
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_sort
.summary_body <- function(values) {
    values <- round(values, 2)

    # y axis goes from bottom to top!
    sample_names <- rev(colnames(values))
    values$type <- rownames(values)

    # Natural sort
    ordered_levels <- str_sort(values$type, numeric = TRUE)
    values$type <- factor(values$type, levels=ordered_levels)

    vals_long <- pivot_longer(values, !"type",
        values_to = "value", names_to = "variable")
    vals_long$variable <- factor(vals_long$variable, levels=sample_names)

    # Make sure NaN values will be written
    vals_long$text_value <- sprintf("%0.2f", round(vals_long$value, digits = 2))

    plot <- ggplot(vals_long, aes(x = .data$type, y = .data$variable, fill = .data$value)) +
        geom_tile(color = "white", linewidth = 0.6) +
        geom_text(aes(label = .data$text_value), size = 4) +
        coord_fixed() +
        scale_y_discrete(position = "right") +
        theme_classic(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "left",
            legend.justification.left = "top",
            legend.direction = "vertical",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks = element_blank()
        )

    plot
}


#' Scatterplot of values in a GRanges object with mcols x and y
#'
#' Plots a scatter plot from a GRanges object that already has x and y columns.
#'
#' @param gr GRanges
#' @param highlight List of GRanges to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a locus to be highlighted
#' @param highlight_label Labels for the highlight groups.
#'  If not provided, column names are used.
#' @param highlight_colors Array of color values for the highlighting groups
#' @param density Plot density tiles for global distribution instead of points.
#' @importFrom ggplot2 ggplot geom_point aes
#'  geom_bin2d scale_fill_gradient scale_color_manual
#' @return A named list where plot is a ggplot object and calculated is a list
#'   of calculated values (for verbose mode).
.scatterplot_body <- function(gr, highlight = NULL,
                                minoverlap = 0L,
                                highlight_label = NULL,
                                highlight_colors = NULL,
                                density = FALSE) {
    df <- data.frame(gr)

    extra_plot <- NULL
    extra_colors <- NULL

    if (!is.null(highlight)) {
        highlight_values <- .multi_ranges_overlap(
            gr,
            highlight,
            highlight_label,
            minoverlap
        )

        extra_plot <- geom_point(
            data = highlight_values,
            aes(x = .data$x, y = .data$y, color = .data$group),
            alpha = 0.8
        )

        if (!is.null(highlight_colors)) {
            extra_colors <- scale_color_manual(values = highlight_colors)
        }
    }

    points <- geom_point(color = "#bbbbbb", alpha = 0.7)
    if (density) {
        points <- list(
            geom_bin2d(binwidth=0.05),
            scale_fill_gradient(low="#dddddd", high="#B22222")
        )
    }

    p <- ggplot(df, aes(x = .data$x, y = .data$y)) +
        points +
        extra_plot +
        extra_colors

    p
}

#' Internal function to plot ranges in violin plot with a highlighted GRanges.
#'
#' @param gr GRanges object with as many columns as samples to plot.
#' @inheritParams .scatterplot_body
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect any_of
#' @importFrom ggplot2 ggplot aes geom_violin theme geom_jitter
#'   scale_color_manual
#' @return A named list where plot is a ggplot object and calculated is a list
#'   of calculated values (for verbose mode).
.violin_body <- function(gr,
                    highlight = NULL,
                    minoverlap = 0L,
                    highlight_label = NULL,
                    highlight_colors = NULL) {

    bwnames <- names(mcols(gr))
    df <- data.frame(gr)
    bin_id <- c("seqnames", "start", "end")

    bins_long <- pivot_longer(df[, c(bin_id, bwnames)], !any_of(bin_id), names_to = "variable", values_to = "value")
    # melted_bins <- melt(df[, c(bin_id, bwnames)], id.vars = bin_id)

    extra_plot <- NULL
    extra_colors <- NULL

    n_highlighted_points <- 0

    if (!is.null(highlight)) {
        highlight_data <- .convert_and_label_loci(highlight, highlight_label)

        highlight_values <- .multi_ranges_overlap(
            gr,
            highlight_data$ranges,
            highlight_data$labels,
            minoverlap
        )

        n_highlighted_points <- nrow(highlight_values)
        bin_id <- c("seqnames", "start", "end", "width", "strand", "group")
        highlight_long <- pivot_longer(
            highlight_values, !any_of(bin_id),
            values_to = "value", names_to="variable"
        )

        extra_plot <- geom_jitter(
            data = highlight_long,
            aes(x = .data$variable, y = .data$value, color = .data$variable),
            alpha = 0.7
        )

        if (!is.null(highlight_colors)) {
            extra_colors <- scale_color_manual(values = highlight_colors)
        }
    }

    p <- ggplot(bins_long, aes(x = .data$variable, y = .data$value)) +
        geom_violin(fill = "#cccccc") +
        extra_plot +
        extra_colors
    p
}


#' Helper function for calculating guide lines and labels in heatmap
#'
#' @param loci Number of loci (rows)
#' @param nbins Number of bins (columns)
#' @param expand Whether to apply padding to the ticks.
#'   Otherwise a strange box shows on heatmap.
#' @inheritParams plot_bw_heatmap
#' @return A list of ggproto objects to be plotted.
.heatmap_lines <- function(loci, nbins, bin_size, upstream, downstream, mode,
                            expand = TRUE) {

    axis_breaks <-
        .profile_breaks(nbins, upstream, downstream, bin_size, mode)
    axis_labels <- .profile_labels(upstream, downstream, mode)

    lines <- axis_breaks[2]
    if (mode == "stretch") {
        lines <- axis_breaks[2:3]
    }

    x <- scale_x_continuous(breaks = axis_breaks, labels = axis_labels)

    if (expand == TRUE) {
        x <- scale_x_continuous(
            breaks = axis_breaks,
            labels = axis_labels,
            expand = c(0, 0)
        )
    }

    gline <- geom_vline(
        xintercept = lines,
        linetype = "dashed",
        color = "#111111",
        linewidth = 0.2
    )

    list(x, gline)
}


#' Calculate overlap of one GRanges with a main GRanges object
#'
#' @param main_ranges Main GRanges to which overlap is calculated.
#' @param other_ranges A list of GRanges objects
#' @param labels Labels to each of the other_ranges.
#' @param minoverlap Minimum overlap to consider an overlap.
#'
#' @importFrom purrr partial map2
#' @return A data.frame in tall format with the values of the overlapping loci.
#'    Loci returned belong to other_ranges, NOT main_ranges. A group field
#'    is added as factor.
.multi_ranges_overlap <- function(main_ranges, other_ranges, labels,
                                    minoverlap) {
    .label_df <- function(df, name) {
        data.frame(df, group = name)
    }

    subset_func <- partial(
        IRanges::subsetByOverlaps,
        x = main_ranges,
        minoverlap = minoverlap
    )

    ranges_subset <- lapply(other_ranges, subset_func)
    subset_df <- lapply(ranges_subset, data.frame)

    df_values_labeled <- map2(subset_df, labels, .label_df)
    highlight_values <- do.call(rbind, df_values_labeled)

    # Order of factors need to match to assign properly colors to points
    highlight_values$group <- factor(highlight_values$group, levels = labels)

    highlight_values
}


#' Summary heatmap colorscale.
#'
#' @param norm_mode Type of normalization.
#' @param bg_bwfiles Background files.
#'
#' @return a ggproto object
.colorscale <- function(norm_mode, bg_bwfiles) {
    legend_label <- .make_norm_label(norm_mode, bg_bwfiles)
    colorscale <- scale_fill_gradient(
        name = legend_label,
        low = "white",
        high = "#b22222",
        na.value = "#cccccc"
    )
    if (!is.null(bg_bwfiles)) {
        colorscale <- scale_fill_gradient2(
            name = legend_label,
            low = "#2e6694",
            mid = "white",
            high = "#b22222",
            na.value = "#cccccc"
        )
    }
    colorscale
}

#' Calculate color limits from a value matrix and provided parameters
#'
#' @param m Value matrix.
#' @param zmin Min value. Overrides percentile.
#' @param zmax Max value. Overrides percentile.
#'
#' @return A pair of c(min, max)
.color_limits <- function(m, zmin, zmax) {
    # colorscale limits percentiles: 0.01 - 0.99
    zlim <- quantile(unlist(m), c(0.01, 0.99), na.rm = TRUE)

    if (zlim[[1]] == zlim[[2]]) {
        warning("Narrow colorscale (quantile(0.01) == quantile(0.99)). Heatmap is
                probably full of zero values. Leaving it raw.")
        zlim <- c(min(unlist(m), na.rm = TRUE),
                  max(unlist(m), na.rm = TRUE))
    }

    if (!is.null(zmin)) {
        zlim[[1]] <- zmin
    }

    if (!is.null(zmax)) {
        zlim[[2]] <- zmax
    }
    zlim
}

#' Compute where the break ticks go in heatmap and profile
#'
#' @param nrows Number of bins
#' @inheritParams plot_bw_profile
#'
#' @return Array of numeric
.profile_breaks <- function(nrows, upstream, downstream, bin_size, mode) {
    upstream_nbins <- floor(upstream / bin_size)
    downstream_nbins <- floor(downstream / bin_size)

    # index value starts at 1
    axis_breaks <- c(1, upstream_nbins + 1, nrows + 1)

    # Put ticks on the edges
    axis_breaks <- axis_breaks - 0.5

    if (mode == "stretch") {
        axis_breaks <- c(
            1,
            upstream_nbins + 1,
            nrows - downstream_nbins + 1,
            nrows + 1
        )
        # center ticks on the middle of the bins
        axis_breaks <- axis_breaks - 0.5
    }
    axis_breaks
}

#' Make profile labels
#'
#' @param upstream Basepairs upstream
#' @param downstream Basepairs downstream
#' @param mode Plot mode (stretch, start, end, center)
#'
#' @return Array of strings
.profile_labels <- function(upstream, downstream, mode) {
    if (mode == "stretch") {
        c(
            paste("-", upstream / 1000, "kb", sep = ""),
            "start", "end",
            paste("+", downstream / 1000, "kb", sep = "")
        )
    } else {
        c(
            paste("-", upstream / 1000, "kb", sep = ""),
            mode,
            paste("+", downstream / 1000, "kb", sep = "")
        )
    }
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
#' @importFrom methods is
#' @return a named list with ranges and labels.
.convert_and_label_loci <- function(loci_sets, labels) {
    gr_list <- loci_sets
    lab_list <- labels
    if (!is.null(gr_list)) {
        if (is(gr_list, "character")) {
            # list is used because GRanges c() concatenates ranges
            gr_list <- lapply(loci_sets, rtracklayer::import, format = "BED")
            if (is.null(labels)) {
                lab_list <- basename(loci_sets)
            }
        }
        else {
            if (is.null(labels)) {
                stop("GRanges used as highlight loci but no labels provided")
            }
            if (!is(gr_list, "list")) {
                # If single GRanges was passed it needs to be converted to list
                gr_list <- list(gr_list)
            }
        }
    }
    if (length(gr_list) != length(lab_list)) {
        stop("Highlight loci sets don't match the number of labels provided")
    }

    list(ranges = gr_list, labels = lab_list)
}
