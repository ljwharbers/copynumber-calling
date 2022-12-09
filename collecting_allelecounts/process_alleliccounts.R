#!/usr/bin/env Rscript

## Author: Luuk Harbers
## Date: 2022-10-22
## Script to process allelic counts

## Load/install packages
packages = c("data.table", "argparser")
invisible(sapply(packages, function(x) suppressMessages(require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))))

## Parse arguments
parser = arg_parser("Process allelic counts")
parser = add_argument(parser, "--input", short = "-i", help = "Path to allelic counts", nargs = 1)
parser = add_argument(parser, "--output", short = "-o", help = "Path to output file", nargs = 1)
argv = parse_args(parser)

# In case of testing
# argv = list()
# argv$input = "/mnt/AchTeraD/data/scCUTseq/prostate/P3_subs/baf/NZ217_GTTGTACACTT_allelecounts.tsv"
# argv$output = "/mnt/AchTeraD/data/scCUTseq/prostate/P3_subs/baf/NZ217_GTTGTACACTT_baf.tsv"

# Load in file
counts = fread(argv$input)
counts[, Good_depth := as.numeric(Good_depth)]
counts = counts[Good_depth > 0]

# Write gzipped file
fwrite(counts, argv$output, row.names = F, col.names = T, quote = F, sep = "\t")
