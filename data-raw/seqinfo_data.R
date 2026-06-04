# This code gets some sequence information for Seqinfo once, and then keeps it
# as internal, making bw_bins less dependent on the internet. If the genome
# looked for is not among this, it will make use of GenomeInfoDb API
genome_names <- c(
  "mm9", "mm10", "mm39", "hg19", "hg38", "dm6",
  "GRCm38.p6", "GRCm39", "GRCh38", "GRCh38.p14",
  "sacCer3", "ce11"
)

seqinfo_data <- lapply(genome_names, function(x) { Seqinfo::Seqinfo(genome = x) } )
names(seqinfo_data) <- genome_names
usethis::use_data(seqinfo_data, overwrite = TRUE, internal = TRUE)
