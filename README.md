# Summary and scope

![r-cmd-check](https://github.com/cnluzon/wigglescout/workflows/r-cmd-check/badge.svg)

`wigglescout` is an R library that allows you to calculate summary values 
across bigWig files and BED files and visualize them in a genomics-relevant
manner. It is based on broadly used libraries such as `rtracklayer` and
`GenomicRanges`, among others for calculation, and mostly `ggplot2` for
visualization. You can look at the `DESCRIPTION` file to get more information
about all the libraries that make this one possible. 

There are also many other tools whose functionality overlaps a little or
much with `wigglescout`, but there was no single tool that 
included all that I needed. The aim of this library is therefore not to 
replace any of those tools, or to provide a silver-bullet solution to genomics
data analysis, but to provide a comprehensive, yet simple enough set of tools
focused on bigWig files that can be used entirely from the R environment
without switching back and forth across tools.

Other tools and libraries for akin purposes that you may be looking for include:
`deepTools`, `SeqPlots`, `bwtool`, `wiggletools`, and the list is endless!

`wigglescout` allows you to summarize and visualize the contents of 
bigWig files in two main ways:

- **Genome-wide**. Genome is partitioned on equally-sized bins and
their aggregated value is calculated. Useful to get a general idea of the
signal distribution without looking at specific places.
- **Across sets of *loci***. This can be either summarized categories, or
individual values, as in genome-wide analyses.

`wigglescout` functionality is built in two layers. Names of functions that
calculate values over bigWig files start with `bw_`. These return `GRanges`
objects when possible, `data.frame` objects otherwise (i.e. when values are
summarized over some category, genomic location is lost in this process).

On the other hand, functions that plot such values and that usually make
internal use of `bw_` functions, start with `plot_`.

## Installation

`wigglescout` is a package under development. You will need `remotes` to
install it (and `devtools` if you plan to work on it):

    install.packages(c('devtools', 'remotes'))
    
Additionally, there was an issue in the past with installing dependencies that
come from `BioConductor` repository. This seems to have been fixed now, but if
you run into problems, I recommend installing manually these dependencies
before running the installation:

    install.packages(('BiocManager'))
    BiocManager::install(c('GenomeInfoDbData',
        'GenomeInfoDb',
        'GenomicRanges',
        'rtracklayer',
        'BSgenome.Mmusculus.UCSC.mm9',
        'BSgenome.Hsapiens.UCSC.hg38'))

Then you can install directly from this GitHub repository:

    library(remotes)
    install_github('cnluzon/wigglescout', build_vignettes=TRUE)

## Getting started

The vignettes as they can give a comprehensive overview of what is
available in the package. You can check the vignettes with
`browseVignettes("wigglescout")`.

These are the groups of functions that are included in `wigglescout`:

- `bwtools`. Functionality to handle `bigWig` files. Importing, binning
    and intersecting, aggregating `BED` files. This relies heavily on 
    `rtracklayer` calculations.
    
- `bwplot`. Make plots out of `bwtools` calculations.
    
## Troubleshooting

**Q**: When running `install_github` I get the following error:

    Error: package or namespace load failed for ‘GenomeInfoDb’ in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]):
    there is no package called ‘GenomeInfoDbData’
    Error: package ‘GenomeInfoDb’ could not be loaded
    Execution halted
    
**A**: This seemed to be a problem that came from installing `Bioconductor`
dependencies. A workaround is installing the `BioConductor` packages manually: 

    if (!requireNamespace('BiocManager', quietly = TRUE))
        install.packages('BiocManager')

    BiocManager::install(c('GenomeInfoDbData',
        'GenomeInfoDb',
        'GenomicRanges',
        'rtracklayer',
        'BSgenome.Mmusculus.UCSC.mm9',
        'BSgenome.Hsapiens.UCSC.hg38'))
 
