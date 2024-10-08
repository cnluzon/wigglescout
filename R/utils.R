#' Make a BigWigFile object out of a path or an URL
#'
#' @param bw BigWig file path or URL
#'
#' @return BigWigFile object
.fetch_bigwig <- function(bw) {
    if (!is.null(bw)) {
        valid_bwfile <- bw
        if (RCurl::url.exists(bw)) {
            valid_bwfile <- tempfile()
            download.file(bw, valid_bwfile)
        }
        BigWigFile(path = valid_bwfile)
    }
}


#' Filter scatterplot data with quantile threshold on both axes
#'
#' @param x GRanges for x axis
#' @param y GRanges for y axis
#' @param col_x Name of column in x to filter
#' @param col_y Name of column in y to filter
#' @param remove_top Return range 0-(1-remove_top). By default returns the
#'     whole distribution (remove_top == 0).
#'
#' @return Named list, where gr is a GRanges object and stats a named list of
#'    values: points: number of points in the figure, NA.x = Number of NA values
#'    on the x axis. filtered.x: number of points that were higher than the
#'    threshold on the x axis. quant.x: Quantile used as a threhsold on the x
#'    axis. Same for the y values.
.filter_scatter_data <- function(x, y, remove_top, col_x = "score", col_y = "score") {
    clean_x <- .remove_top_by_mean(x, remove_top, c(col_x))
    clean_y <- .remove_top_by_mean(y, remove_top, c(col_y))

    clean_gr <- .granges_left_join(
        list(clean_x$ranges[, col_x], clean_y$ranges[, col_y]),
        c("x", "y")
    )

    stats <- list(
        points = length(clean_gr),
        NA.x = clean_x$calculated$na,
        filtered.x = clean_x$calculated$filtered,
        quant.x = .round_ignore_null(clean_x$calculated$quantile),
        NA.y = clean_y$calculated$na,
        filtered.y = clean_y$calculated$filtered,
        quant.y = .round_ignore_null(clean_y$calculated$quantile)
    )

    list(ranges = clean_gr, stats = stats)
}


#' Filter violin plot data with quantile threshold
#'
#' Removes top quantile loci from the distribution. Values are removed by mean
#' rather than maximum.
#'
#' @param gr GRanges object
#' @param remove_top (1-remove_top) quantile will be removed.
#' @param columns mcols in GRanges to use for selection.
#'
#' @return A named list with ranges = GRanges object and stats a named list with
#'   some calculated values.
.filter_violin_data <- function(gr, remove_top, columns) {
    clean_gr <- .remove_top_by_mean(gr, remove_top, columns)

    stats <- list(
        points = length(gr),
        removed = clean_gr$calculated$filtered,
        NAs = clean_gr$calculated$na,
        quantile_cutoff = .round_ignore_null(clean_gr$calculated$quantile)
    )

    list(ranges = clean_gr$ranges, stats = stats)
}


#' Preprocess heatmap values for plotting
#'
#' This includes capping on zmin / zmax values. If not provided, percentiles
#' 0.01 and 0.99 are used, to match plotting of other similar tools.
#' If provided, a specific order of rows is used. Otherwise, rows are sorted
#' by rowMeans.
#' If the number of rows in the matrix is larger than a certain max_rows_allowed
#' value, number of rows is reduced to max_rows_allowed. Rather than subsampling
#' each row then will represent an average of the underlying rows.
#'
#' @param values Matrix to preprocess
#' @param zmin Minimum value to show in color.
#' @param zmax Maximum value to show in color.
#' @param order_by Row order (or NULL).
#' @param max_rows_allowed Maximum heatmap resolution.
#'
#' @return A named list with values and stats calculated.
.preprocess_heatmap_matrix <- function(values, zmin, zmax, order_by, max_rows_allowed) {
    # Order matrix by mean and transpose it (image works flipped)
    if (is.null(order_by)) {
        order_by <- order(rowMeans(values), decreasing = FALSE)
    }
    m <- t(values[order_by, ])
    nvalues <- nrow(m) * ncol(m)
    n_non_finite <- length(m[!is.finite(m)])
    m[!is.finite(m)] <- NA

    zlim <- .color_limits(m, zmin, zmax)
    zmin <- zlim[[1]]
    zmax <- zlim[[2]]

    n_bottom_capped <- length(m[m < zmin])
    n_top_capped <- length(m[m > zmax])
    m[m < zmin] <- zmin
    m[m > zmax] <- zmax

    df <- data.frame(m) %>%
        mutate(x = row_number()) %>%
        pivot_longer(!"x", names_to = "col", values_to = "value") %>%
        mutate(y = as.numeric(gsub("X", "", .data$col)))

    n_rows <- max(df$y)
    downsample_factor <- NA
    if (n_rows > max_rows_allowed) {
        # Downsample rows only and downsample only enough to fit max_rows.
        warning("Large matrix: ", n_rows, ". Downscaled to ", max_rows_allowed)
        downsample_factor <- round(n_rows / max_rows_allowed)

        df <- df %>% group_by(
                x = .data$x,
                y = ((.data$y-1) %/% downsample_factor) + 1
            ) %>% summarise(value = mean(.data$value))
    }

    stats <- list(ncells = nvalues, zmin = zmin, zmax = zmax,
                  top_capped_vals = n_top_capped,
                  bottom_capped_vals = n_bottom_capped,
                  non_finite = n_non_finite,
                  downsample_factor = downsample_factor)

    list(values = df, stats = stats)
}


#' Get matching parameters of a target function with current context
#'
#' This is a helper function for wrapping functions such as plotting. It
#' gets the formal arguments of the target function func, and gets their
#' values from the context given in current, that will have the same name.
#' It returns a named list argument - value that can be passed to func
#' using do.call()
#'
#' @param func Target function
#' @param current Current context as a named list
#'
#' @return Named list argument - value
.get_wrapper_parameter_values <- function(func, current){
    f <- names(formals(func))
    n <- names(current)
    par <- n[n %in% f]
    par
}

#' GRanges left join operation
#'
#' Perform a left_bind operation on a GRanges list, appending scores to mcols.
#'
#' This assumes that you are trying to join things that match (bins generated
#' from the same parameters, BED intersections from the same BED.)
#'
#' @param grlist A list of GRanges objects that have all the same fields.
#' @param labels Vector of names for the score columns.
#' @param granges Optional granges with name field on it
#' @return A Sorted GRanges object with all the columns.
#' @importFrom tidyselect all_of
.granges_left_join <- function(grlist, labels, granges = NULL) {
    fixed_fields <- c("seqnames", "start", "end", "width", "strand")

    data_frames <- lapply(grlist, data.frame)
    dedup_frames <- lapply(data_frames, unique)

    result <- dedup_frames %>% purrr::reduce(dplyr::left_join, by = fixed_fields)
    colnames(result) <- c(fixed_fields, labels)
    # Include names if granges has them
    if (! is.null(granges)) {
        if ("name" %in% names(mcols(granges))) {
            # If there were equal loci with different names, this would
            # duplicate the corresponding entry
            dedup_granges <- unique(data.frame(granges))
            result <- dplyr::left_join(result, dedup_granges,
                                       by = fixed_fields, multiple = "all") %>%
                dplyr::select(all_of(c(fixed_fields, labels, "name")))
        }
    }


    makeGRangesFromDataFrame(result, keep.extra.columns = TRUE)
}

#' Make a string out of a named list.
#'
#' @param named_list A named list
#' @return A string
.key_value_string <- function(named_list) {
    paste(names(named_list), lapply(named_list, .format_value), sep = ":", collapse = ", ")
}

.format_value <- function(v) {
    if(is.numeric(v)) {
        v <- sprintf("%.4f", v)
        v <- sub("\\.?0+$", "", v)
    }
    v
}

#' Make a string out of a named list. Split into lines if too wide.
#'
#' @param named_list A named list
#'
#' @return A string
.limited_size_caption_line <- function(named_list) {
    size_limit <- 3
    chunks <- split(named_list, ceiling(seq_along(named_list)/size_limit))
    paste(vapply(chunks, .key_value_string, character(1)), collapse="\n")
}


#' Check if a variable contains a list of things or array of potentially file names
#'
#' @param loci A list of GRanges / files / mixed (or not)
#'
#' @return TRUE if variable contains multiple loci, FALSE otherwise
.is_multiple_loci <- function(loci) {
  (is(loci, "list") && length(loci) > 1) || (is(loci, "character") && length(loci) > 1)
}


#' Label loci with a list of labels
#'
#' If some labels are duplicated, numbers are added to the end
#'
#' @param loci List of loci
#' @param labels List of labels
#'
#' @return A list of labels
.label_multiple_loci <- function(loci, labels) {
    if (is.null(labels)) {
        labels <- lapply(loci, .make_label_from_object)
        if (length(unique(labels)) < length(loci)) {
            warning("Unlabeled objects or repeated labels. ",
                    "Adding numeric indices.")
            labels <- paste(labels, seq_len(length(labels)), sep = "_")
        }
    }
    labels
}

#' Checks how many loci
#'
#' @param loci GRanges or BED
#'
#' @return An integer
.loci_length <- function(loci) {
    if (is.character(loci)) {
        length(rtracklayer::import(loci))
    } else {
        length(loci)
    }
}


#' Processes a loci and returns a GRanges object. Removes anything that is not
#' a name field
#'
#' @param loci Either a BED file or a GRanges object
#' @return A GRanges object
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges `mcols<-`
#' @importFrom methods is
.loci_to_granges <- function(loci) {
    bed <- loci
    if (is(loci, "character")){
        bed <- import(loci)
    }
    if ("name" %in% names(mcols(bed))) {
        bed <- bed[, "name"]
    }
    else {
        mcols(bed) <- NULL
    }
    bed
}


#' Make a string to put as caption in verbose mode. Includes system date.
#'
#' @param params Named list with relevant parameters and their values
#' @param outcome Named values with relevant outcomes and their values
#' @param verbose Logical. If TRUE all information is printed. If FALSE returns
#'   a blank caption.
#' @importFrom utils packageVersion
#' @return A caption string
.make_caption <- function(params, outcome, verbose) {
    if (verbose) {
        verbose_params <- .limited_size_caption_line(params)
        verbose_crop <- .limited_size_caption_line(outcome)

        date <- format(Sys.time(), "%a %b %d %X %Y")
        pkg_version <- paste("wigglescout v.", packageVersion("wigglescout"))
        date <- paste(date, pkg_version, sep = ' - ')

        paste(verbose_params, verbose_crop, date, sep = "\n\n")
    }
}


#' Get a valid label
#'
#' @param obj Something to convert to label
#' @param max_length Max length allowed for a label to have (35)
#'
#' @return A valid label name
.make_label_from_object <- function(obj, max_length = 35) {
    if (is.character(obj)) {
        filename_clean <- basename(tools::file_path_sans_ext(obj))
        sapply(make.names(filename_clean), .trunc_str, max_length = max_length)
    } else {
      sapply(make.names(class(obj)), .trunc_str, max_length = max_length)
    }
}


#' Generate a human-readable normalization function string
#'
#' @param f String representing normalization function.
#' @param bg Background file.
#'
#' @return A string describing normalization.
.make_norm_label <- function(f, bg) {
    label <- "RPGC"
    if (!is.null(bg)) {
        label <- switch(f,
            "fc" = paste(label, " / bg", sep = ""),
            "log2fc" = paste("log2(", label, " / bg)", sep = "")
        )
    }
    label
}

#' Generate a human-readable normalization function string including
#'
#' @param f String representing normalization function.
#' @param fg Foreground file
#' @param bg Background file.
#'
#' @return A string describing normalization.
.make_norm_file_label <- function(f, fg, bg) {
    paste(.make_label_from_object(fg), "-", .make_norm_label(f, bg))
}



#' Take a column of a data frame to use as row names, drop such column and
#' reorder the data frame so row names follow natural order.
#'
#' @param df A data frame.
#' @param col The column (usually name).
#' @importFrom stringr str_sort
#' @return A sorted df
.natural_sort_by_field <- function(df, col) {
    rownames(df) <- df[, col]
    order <- str_sort(df[, col], numeric = TRUE)
    df[, col] <- NULL
    df[order, , drop = FALSE]
}


#' Remove top loci in a GRanges column by mean of specified values
#'
#' @param granges GRanges object. Must have numerical mcols, and names of those
#'   needs to match with columns.
#' @param quantile Top quantile. Must be a number between zero and 1.
#' @param columns Which columns to use. If more than one, quantile will be
#'   selected according to means.
#'
#' @return A named list, fields ranges for the resulting GRanges, plus
#'   calculated values: quantile_value, filtered, na_values.
.remove_top_by_mean <- function(granges, quantile, columns) {
    n_filtered <- 0
    top_quantile <- NA

    valid_columns <- data.frame(mcols(granges))
    valid_columns <- valid_columns[, columns, drop = FALSE]
    n_na <- sum(is.na(valid_columns))

    if (quantile > 0) {
        if (ncol(mcols(granges)) > 1) {
            means <- rowMeans(valid_columns)
            top_quantile <- quantile(means, probs = c(1 - quantile),
                                     na.rm = TRUE)
            granges$means <- means
            granges <- granges[!is.na(granges$means), ]

            n_filtered <- length(granges[granges$means > top_quantile, ])
            granges <- granges[granges$means <= top_quantile, ]
            granges$means <- NULL
        }
        else {
            top_quantile <- quantile(mcols(granges)[, 1],
            probs = c(1 - quantile), na.rm = TRUE
            )
            granges <- granges[!is.na(mcols(granges)[, 1]), ]
            n_filtered <- length(granges[mcols(granges)[, 1] > top_quantile, ])
            granges <- granges[mcols(granges)[, 1] <= top_quantile, ]
        }
    }
    list(
        ranges = granges,
        calculated = list(
            na = n_na, filtered = n_filtered,
            quantile = unname(top_quantile)
        )
    )
}

#' Round a value. If it's NULL, returns NULL.
#'
#' @param v Value
#' @param digits How many digits round
#'
#' @return Rounded value, or NULL
.round_ignore_null <- function(v, digits=3) {
    rounded <- v
    if (!is.null(v)) {
        rounded <- round(v, digits)
    }
    rounded
}

#' Set default theme as classic with larger font size
#' @importFrom ggplot2 theme_classic theme %+replace%
#' @return ggproto object
.theme_default <- function() {
    theme_classic(base_size = 14) %+replace%
      theme(plot.caption = element_text(size = 9, hjust = 1))
}

#' Truncate a string to a max length
#' @param s String to truncate
#' @param max_length Length to truncate s to (35)
.trunc_str <- function(s, max_length = 35) {
  if(nchar(s) > max_length) {
    s <- substr(s, 1, max_length)
  }
  s
}

#' Validate a category array
#'
#' Checks whether the number of categories is reasonable for an aggregation.
#' Throws a warning if it finds more than 50 different values.
#'
#' @param cat_values An array of values
#' @return Nothing (NULL) on success.
.validate_categories <- function(cat_values) {
    max_categories <- 50
    # Test number of values in group_col
    ncat <- length(levels(as.factor(cat_values)))
    if (ncat > max_categories) {
        warning("Number of values in group column field very large: ", ncat,
            " (does BED file have unique IDs instead of categories?)"
        )
    }
}

#' Validate an array of paths
#'
#' Check that a list of files is valid: not empty and contents exist.
#'
#' @param filelist An array of files
#' @importFrom RCurl url.exists
#' @return Nothing (NULL) on success.
.validate_filelist <- function(filelist) {
    if (length(filelist) == 0) {
        stop("File list provided is empty.")
    }
    existence_flag <- file.exists(filelist) | RCurl::url.exists(filelist)
    if (!all(existence_flag)) {
        msg <- paste("Files not found:", filelist[!existence_flag])
        stop(msg)
    }
}

#' Validate that a locus parameter is valid. Checks for paths and also whether
#' the parameter is otherwise a GRanges object, which is also valid.
#'
#' @param locus_param Parameter to validate
#' @importFrom methods is
#' @return Nothing (NULL) on success.
.validate_locus_parameter <- function(locus_param) {
    if (is(locus_param, "character")) {
        .validate_filelist(locus_param)
    }
    else {
        if (!is(locus_param, "GRanges")) {
            msg <- paste0("Unexpected type: ", class(locus_param))
            stop(msg)
        }
    }
}

.validate_input_numbers <- function(bw, loci) {
    if ((is(loci, "list") && length(loci) > 1) ||
        (is(loci, "character") && length(loci) > 1)) {
        if (length(bw) > 1) {
            stop("If multiple loci provided only a single bwfile is allowed")
        }
    }
}


#' Validate that group col exists in granges
#'
#' @param granges GRanges object to check
#' @param group_col Group column name. Usually, name.
#' @return Nothing (NULL) on success.
.validate_group_col <- function(granges, group_col) {
    if (!group_col %in% names(mcols(granges))) {
        stop("Invalid group column not present in granges", group_col)
    }
}


#' Validate profile and heatmap relevant parameters
#'
#' @param bin_size Bin size. Must be a positive number.
#' @param upstream Upstream bp. Must be positive and larger than bin size.
#' @param downstream Downstream bp. Must be positive and larger than bin size.
#' @return Nothing (NULL) on success.
#'
.validate_profile_parameters <- function(bin_size, upstream, downstream) {
    if (bin_size <= 0) {
        stop("bin size must be a positive value: ", bin_size)
    }

    if (upstream <= 0) {
        stop("upstream size must be a positive value: ", upstream)
    }

    if (downstream <= 0) {
        stop("downstream size must be a positive value: ", downstream)
    }

    if (bin_size > upstream || bin_size > downstream) {
        stop("bin size must be smaller than flanking regions")
    }
}


#' Check that a list of bigWig files have the same Seqinfo.
#' Throws a warning if they don't.
#'
#' @param bwfiles List of bigWig files
#' @param show_max Max number of references to show in the warning
#'
#' @return TRUE if they share all references, FALSE otherwise.
.validate_references <- function(bwfiles, show_max = 25) {
  seqinfos <- lapply(lapply(bwfiles, BigWigFile), GenomeInfoDb::seqinfo)
  seqnames_all <- lapply(seqinfos, GenomicRanges::seqnames)
  .match_seqnames(seqnames_all, show_max, 5)
}

#' Check that a list of bigWig files have the same Seqinfo as an external object
#' (usually a pre-built tile set). Throws a warning if they don't.
#'
#' @param bwfiles List of bigWig files
#' @param seqinfo Seqinfo object
#' @param show_max Max number of shared references to show in the warning
#'
#' @return TRUE if they share all references, FALSE otherwise.
.validate_references_with_external_seqinfo <- function(bwfiles, seqinfo, show_max = 25) {
  bw_seqinfos <- lapply(lapply(bwfiles, BigWigFile), rtracklayer::seqinfo)
  seqnames_all <- lapply(c(bw_seqinfos, seqinfo), GenomicRanges::seqnames)
  .match_seqnames(seqnames_all, show_max, 5)
}

#' Match seqname lists and throws the relevant warning
#'
#' @param seqnamelist Names list to compare
#' @param show_max_shared Maximum shared references to print
#' @param show_max_missing Maximum missing references to print
#'
#' @return TRUE if they match, FALSE if they don't
.match_seqnames <- function(seqnamelist, show_max_shared, show_max_missing) {
  all_equal <- Reduce(
    function(x, y) { if (base::setequal(x, y)) x else FALSE
    }, seqnamelist
  )

  if (isFALSE(all_equal)) {
    common_refs <- Reduce(intersect, seqnamelist)
    n_missing <- purrr::partial(setdiff, y = common_refs)
    missing <- lapply(seqnamelist, n_missing)
    missing_any <- Reduce(union, missing)
    ellipsis <- ""
    ellipsis_missing <- ""
    if(length(common_refs) > show_max_shared) {
      ellipsis <- "..."
    }
    if(length(missing_any) > show_max_missing) {
      ellipsis_missing <- "..."
    }

    msg <- paste(
      "Not all bigWig files or the genome tiling provided share the same sequence info. \n   Common to all",
      paste0("(",length(common_refs),"):"),
      paste(common_refs[1:min(length(common_refs), show_max_shared)], collapse = ", "),
      ellipsis,
      "\n  ",
      "Missing in some of the files:",
      paste(missing_any[1:min(length(missing_any), show_max_missing)], collapse = ", "),
      ellipsis_missing,
      "\n  ",
      "This can be due to different versions of the same reference genome, or",
      "to completely different organisms. Make sure these match!"
    )

    warning(msg)
    return(FALSE)
  }

  return(TRUE)
}
