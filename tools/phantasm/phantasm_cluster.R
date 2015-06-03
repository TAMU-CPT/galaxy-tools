#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_hist <- args[2]
output_newick <- args[3]

if(is.na(output_hist)){
	output_hist <- 'hist.png'
}
if(is.na(output_newick)){
	output_newick <- 'clusters.newick'
}


# Read in data
d_main <- read.csv(args[1], sep="\t")
# Remove leading column (names)
d_main_new <- d_main[-1]
# Convert to a dist object
genomes <- dist(as.matrix(as.vector(d_main_new)))
# Apply labels to the dist object
attr(genomes,"Labels") <- colnames(d_main_new)
# Cluster
hc <- hclust(genomes, method="ave")
# Plot
png(output_hist, width=1500, height=700)
plot(hc)
dev <- dev.off()
# Install/include biocLite
tryCatch(
	library("ctc"),
	error=function(e){
		print("Downloading BioConductor and installing CTC")
		source("http://bioconductor.org/biocLite.R")
		biocLite("ctc")
	}
)

# CTC has a function to write the hc object as a newick tree
write.table(hc2Newick(hc), file=output_newick)
