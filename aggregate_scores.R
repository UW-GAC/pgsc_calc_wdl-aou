library(optparse)
library(dplyr)

option_list <- list(
  make_option("--raw_files", type = "character"),
  make_option("--adjusted_files", type = "character"),
  make_option("--overlap_files", type = "character"),
  make_option("--out_prefix", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

raw.files <- strsplit(opt$raw_files, ",")[[1]]
adjusted.files <- strsplit(opt$adjusted_files, ",")[[1]]
overlap.files <- strsplit(opt$overlap_files, ",")[[1]]

# Read and bind all raw files
raw.data <- lapply(raw.files, read.table, comment.char = "", header = TRUE)
idcol <- raw.data[[1]][1]  # subject id
all.raw.avg <- do.call('cbind', list(idcol, lapply(raw.data, function(df) { dplyr::select(df[, -c(1:3)], contains('AVG')) })))
all.raw.sum <- do.call('cbind', list(idcol, lapply(raw.data, function(df) { dplyr::select(df[, -c(1:3)], contains('SUM')) })))
## the all.raw.sum data.frame includes duplicate 'NAMED_ALLELE_DOSAGE_SUM' columns.
## - Remove all but first occurence
write.table(all.raw.avg, paste0(opt$out_prefix, "_raw_avg_scores.tsv"),
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(all.raw.sum, paste0(opt$out_prefix, "_raw_sum_scores.tsv"),
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

if(!(all(adjusted.files == ""))) {

    # Read and bind all adjusted files
    adjusted.data <- lapply(adjusted.files, read.table, header = TRUE)
    idcol <- adjusted.data[[1]][1]
    all.adjusted <- do.call('cbind', list(idcol, lapply(adjusted.data, function(df) { df[, -1] } )))
    write.table(all.adjusted, paste0(opt$out_prefix, "_adjusted_scores.tsv"),
                append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

}

# Read and bind all overlap files
overlap_data <- do.call(rbind, lapply(overlap.files, read.table, header = TRUE))
write.table(overlap_data, paste0(opt$out_prefix, "_overlap.tsv"), append = FALSE, quote = FALSE,
            sep = "\t", row.names = FALSE)
