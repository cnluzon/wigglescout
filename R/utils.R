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
#' @return A Sorted GRanges object with all the columns.
#' @importFrom GenomeInfoDb sortSeqlevels
.granges_cbind <- function(grlist, labels) {
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

#' Make a string out of a named list.
#'
#' @param named_list A named list
#' @return A string
.key_value_string <- function(named_list) {
    paste(names(named_list), named_list, sep = ":", collapse = ", ")
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


#' Checks how many loci
#'
#' @param loci GRanges or BED
#'
#' @return An integer
.loci_length <- function(loci) {
    if (is.character(loci)) {
        length(rtracklayer::import(loci, format="BED"))
    } else {
        length(loci)
    }
}


#' Processes a loci and returns a GRanges object, sorted and with its seqlevels
#' also sorted.
#'
#' @param loci Either a BED file or a GRanges object
#' @return A GRanges object
#' @importFrom rtracklayer import
#' @importFrom methods is
#' @importFrom GenomeInfoDb sortSeqlevels
.loci_to_granges <- function(loci) {
    bed <- loci
    if (is(loci, "character")){
        bed <- import(loci, format = "BED")
    }

    bed <- sortSeqlevels(bed)
    bed <- sort(bed, ignore.strand = FALSE)
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
#'
#' @return A valid label name
.make_label_from_object <- function(obj) {
    if (is.character(obj)) {
        filename_clean <- basename(tools::file_path_sans_ext(obj))
        make.names(filename_clean)
    } else {
        make.names(class(obj))
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
            "fc" = paste(label, " / background", sep = ""),
            "log2fc" = paste("log2(", label, " / background)", sep = "")
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
    top_quantile <- NULL

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
#' @import ggplot2
#' @return ggproto object
.theme_default <- function() {
    theme_classic(base_size = 18) +
        theme(plot.caption = element_text(size = 11))
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
