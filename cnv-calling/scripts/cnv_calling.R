#!/usr/bin/env Rscript

## Author: Luuk Harbers
## Date: 2021-05-12
## Script to perform cnv calling from a list of binned read counts

## Load/install packages
packages = c("data.table", "argparser", "DNAcopy", "ParDNAcopy", "pbapply", "gtools")
invisible(sapply(packages, function(x) suppressMessages(require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))))

## Parse arguments
parser = arg_parser("Filter bincounts with the use of a blacklist and correct counts for GC content")
parser = add_argument(parser, "--counts", help = "Path to counts", nargs = 1)
parser = add_argument(parser, "--bins", help = "Path to bin file", nargs = 1)
parser = add_argument(parser, "--blacklist", help = "Path to blacklist file", nargs = 1)
parser = add_argument(parser, "--gc", help = "Path to GC content file", nargs = 1)
parser = add_argument(parser, "--alpha", help = "alpha parameter for segmantation with DNACopy", nargs = 1, type = "numeric")
parser = add_argument(parser, "--prune", help = "undo.prune parameter for segmentation with DNAcopy", nargs = 1, type = "numeric")
parser = add_argument(parser, "--type", help = "Type of sequencing 'single' or 'bulk'", nargs = 1, type = "character")
parser = add_argument(parser, "--minploidy", help = "Min ploidy of sample (single-cell only)", 
                      nargs = 1, type = "integer", default = 1.5)
parser = add_argument(parser, "--maxploidy", help = "Max ploidy of sample (single-cell only)", 
                      nargs = 1, type = "integer", default = 6)
parser = add_argument(parser, "--threads", help = "Number of threads to use", nargs = 1, type = "integer")
parser = add_argument(parser, "--output", help = "Path to output file", nargs = 1)
argv = parse_args(parser)

# argv = list()
# argv$counts = "/mnt/AchTeraD/data/CUTseq_2/results/cnv/100000/all-counts.tsv.gz"
# argv$gc = "/mnt/AchTeraD/Documents/Projects/scCUTseq/copy-number-pipeline/files/hg19/GC_variable_100000_48_bwa"
# argv$blacklist = "/mnt/AchTeraD/Documents/Projects/scCUTseq/copy-number-pipeline/files/hg19/hg19-blacklist.v2_adjusted.bed"
# argv$bins = "/mnt/AchTeraD/Documents/Projects/scCUTseq/copy-number-pipeline/files/hg19/variable_100000_48_bwa.bed"
# argv$alpha = 0.0001
# argv$undo.prune = 0.05
# argv$type = "bulk"
# argv$minploidy = 1.5
# argv$maxploidy = 6
# argv$threads = 4
# argv$output = "/mnt/AchTeraD/data/CUTseq_2/results/cnv/100000/cnv.rds"

# Check input parameters
if(!file.exists(argv$counts)) {
  cat("Error: --type argument has to be either 'single' or 'bulk'")
  stop()
}
if(!file.exists(argv$gc)) {
  cat("Error: --gc file does not exist")
  stop()
}
if(!file.exists(argv$blacklist)) {
  cat("Error: --blacklist file does not exist")
  stop()
}
if(!file.exists(argv$bins)) {
  cat("Error: --bins file does not exist")
  stop()
}
if(!argv$type %in% c("single", "bulk")) {
  cat("Error: --type argument has to be either 'single' or 'bulk'")
  stop()
}

# Make list of all data that needs to be included in final output
out = list()
out$counts = fread(argv$counts)
out$gc = fread(argv$gc)
out$blacklist = fread(argv$blacklist)
out$bins = fread(argv$bins)
out$alpha = argv$alpha
out$undo.prune = argv$prune

# Set parameters
type = argv$type
samples = colnames(out$counts)
minploidy = argv$minploidy
maxploidy = argv$maxploidy
threads = argv$threads
output = argv$output


## Define functions
# LOWESS regression for GC normalization
lowess.gc = function(x, y) {
  low = lowess(x, log(y), f=0.05);
  z = approx(low$x, low$y, x)
  return(exp(log(y) - z$y))
}

# Set names
setnames(out$bins, c("chr", "start", "end"))
setnames(out$blacklist, c("chr", "start", "end", "reason"))

# Get overlaps with blacklist
bins = copy(out$bins) # Make copy, do not assign by reference
blacklist = copy(out$blacklist) # Make copy, do not assign by reference
setkey(bins, chr, start, end)
setkey(blacklist, chr, start, end)

# Get the overlaps and retain bins that do not overlap
overlaps = foverlaps(bins, blacklist)
# Reorder
overlaps = overlaps[mixedorder(overlaps$chr)]
overlaps = unique(overlaps, by = c("chr", "i.start", "i.end"))
indices = which(is.na(overlaps$start))

# Filter reads by indices
out$counts = out$counts[indices]
out$gc = out$gc[indices]
out$bins = out$bins[indices]

# Normalize by mean and GC correction
out$counts_gc = out$counts[, pblapply(.SD, function(sample) {
  lowess.gc(out$gc$V1, (sample+1 / mean(sample+1)))
  }, cl = threads)]

# Add 0.001 to 0 ratios and make log ratios
out$counts_gc[out$counts_gc == 0] = 1e-3
out$counts_lrr = log2(out$counts_gc)

# Make CNA object, do smoothing and segmentation
cna = CNA(out$counts_lrr, out$bins$chr, out$bins$start, data.type="logratio", sampleid=colnames(out$counts_lrr)) 
cna_smooth = smooth.CNA(cna)
cna_segment = parSegment(cna_smooth, alpha=out$alpha, min.width=5, undo.splits="prune", undo.prune=out$undo.prune, 
                         njobs = threads, distrib = "Rparallel")
out$segments_long = data.table(cna_segment$output)

# Add segments to output list
segments_rep = out$segments_long[, .(chr = rep(chrom, num.mark), start = rep(loc.start, num.mark), 
                                     seg.mean = rep(seg.mean, num.mark)), by = .(ID)]
segments_rep = segments_rep[mixedorder(segments_rep$chr)]
setorder(segments_rep, ID)
segments_rep[, bin_num := rep(1:nrow(out$counts_lrr), ncol(out$counts_lrr))]
out$segments = dcast(segments_rep, bin_num ~ ID, value.var = "seg.mean")[, -1]

# Reorder out$segments to match sample orders of previous analysis
setcolorder(out$segments, colnames(out$counts_lrr))

# Modify log ratio segments to contain median normalized read counts per bin
median_segments = pblapply(colnames(out$counts_gc), function(cell) {
  # Subset segments_long to only have correct sample
  dt = data.table(num = out$segments_long[ID == cell]$num.mark)
  # Calculate start and end bin for segments
  dt[, end := cumsum(num)]
  dt[, start := shift(end, fill = 0) + 1]
  
  medians = sapply(1:nrow(dt), function(x) {
    return(median(out$counts_gc[[cell]][dt[x, start]:dt[x, end]]))
  })
  return(data.table(rep(medians, dt$num)))
}, cl = threads)
out$segments_read = do.call(cbind, median_segments)
setnames(out$segments_read, colnames(out$counts_gc))

# Renormalize
out$segments_read = data.table(sweep(out$segments_read, 2, colMeans(out$segments_read), '/'))

# Perform integer copy number calling if it's single cell data
if(type == "single"){
  
  # Initialize CN inference variables for SoS method -- Code adapted from ginkgo pipeline
  n_cells = ncol(out$segments_read)
  ploidy = rbind(c(0,0), c(0,0))
  CNgrid = seq(minploidy, maxploidy, by = 0.05)
  n_ploidy = length(CNgrid)  # Number of ploidy tests during CN inference
  CNmult = matrix(0, n_ploidy, n_cells)
  CNerror = matrix(0, n_ploidy, n_cells)
  outerColsums = matrix(0, n_ploidy, n_cells)
  CN = numeric(n_cells)
  
  # Loop through cells and get integer copy number
  res = pbsapply(1:n_cells, function(cell) {
    outerRaw = as.matrix(out$segments_read[, cell, with = F]) %o% CNgrid
    outerRound = round(outerRaw)
    outerDiff = (outerRaw - outerRound) ^ 2
    outerColsums[, cell] <<- colSums(outerDiff, na.rm = FALSE, dims = 1) ## Assignment outside scope
    CNmult[, cell] <<- CNgrid[order(outerColsums[, cell])] ## Assignment outside scope
    CNerror[, cell] <<- round(sort(outerColsums[, cell]), digits=2) ## Assignment outside scope
    
    # Define CN multiplier and assign to CN vector
    CN[cell] <<- CNmult[1, cell] ## Assignment outside scope
    # return integer copy number profile
    return(round(out$segments_read[, cell, with = F] * CN[cell]))
  })
  
  # Save to out
  out$copynumber = data.table(do.call(cbind, res))
  names(CN) = colnames(out$segments_read)
  out$cn = CN
  
  # Write stats
  out$stats = data.table(
    sample = samples
  )
  out$stats[, total_reads := sapply(sample, function(x) sum(out$counts[[x]]))]
  out$stats[, mean_reads := sapply(sample, function(x) mean(out$counts[[x]]))]
  out$stats[, median_reads := sapply(sample, function(x) median(out$counts[[x]]))]
  #out$stats[, spikiness := sapply(sample, function(x){ sum(abs(diff(x))) / sum(x)})]
  out$stats[, MSRSI_non_integerness := NULL]
  out$stats[, MBRSI_dispersion_non_integerness := NULL]
  out$stats[, MBRSM_dispersion := NULL]
  out$stats[, coef_variation := NULL]
  out$stats[, mean_absolute_deviation := NULL]
  out$stats[, total_halfiness := NULL]
  out$stats[, scaled_halfiness := NULL]
  out$stats[, mean_state_mads := NULL]
  out$stats[, mean_state_vars := NULL]
  out$stats[, breakpoints := NULL]
  out$stats[, mean_cn := NULL]
  out$stats[, mode_cn := NULL]
  
} else {
  # Write stats for bulk
  out$stats = data.table(
    sample = samples
  )
  out$stats[, total_reads := sapply(sample, function(x) sum(out$counts[[x]]))]
  out$stats[, mean_reads := sapply(sample, function(x) mean(out$counts[[x]]))]
  out$stats[, median_reads := sapply(sample, function(x) median(out$counts[[x]]))]
  out$stats[, coef_variation := NULL]
  out$stats[, mean_absolute_deviation := NULL]
  out$stats[, breakpoints := NULL]
  
}

# Write output
saveRDS(out, output)
