# wigglescout

![r-cmd-check](https://github.com/cnluzon/wigglescout/workflows/r-cmd-check/badge.svg)

R package to explore and visualize genomics wig and bigWig data based on robust
and broadly used BioConductor R packages such as `rtracklayer`, `GenomicRanges`,
among many other great packages (see dependencies in the `DESCRIPTION` file :) ).

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
 
