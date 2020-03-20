library(data.table)

args <- commandArgs(TRUE)
my_sample_type <- args[2]

plot_figures <- function(case, output_name) {

case$log2fc <- as.numeric(case$log2fc)

fc_threshold <- sort(abs(case$log2fc), decreasing=T)[100]
case_filtered <- case[which(abs(case$log2fc) >= fc_threshold), ]

case_sorted <- case_filtered[order(case_filtered$log2fc), ]

case_sorted[which(case_sorted$Gene %in% cosmic$Gene_Symbol), ]$Gene <- paste0( "***   ", case_sorted[which(case_sorted$Gene %in% cosmic$Gene_Symbol), ]$Gene)

my_col <- rep("gray", dim(case_sorted)[2])
my_col[which(case_sorted$log2fc > 1)] <- "orange"
my_col[which(case_sorted$log2fc < -1)] <- "blue"

my_xmax <- round(max(abs(case_sorted$log2fc), na.rm=T)) + 3
my_xmin <- round(min(case_sorted$log2fc, na.rm=T)) - 3

    write('A', stderr())
    write(my_xmax, stderr())
    write(my_xmin, stderr())
    
pdf(paste0("barplot_", output_name, "_all_top100.pdf"))

barplot(case_sorted$log2fc, names.arg=case_sorted$Gene, horiz=T, las=1, col=my_col, xlab="log2fc", xlim=c(my_xmin, my_xmax), cex.axis=0.1, cex.names=0.3, main=paste0(dim(case_filtered)[1], " proteins"))
    
abline(v=fc_threshold, lty=3)
abline(v=-fc_threshold, lty=3)
dev.off()

}


list_tsv <- system(command="ls *_preprocessed.tsv", intern=T)

cosmic <- fread("/usr/local/data/Cancer_Census_all_072018_COSMIC.csv")

if(my_sample_type=="Melanoma") {
  list_prot <- fread("/usr/local/data/Melanoma_marker_IDs.csv", header=T)
} else if(my_sample_type=="OVCA") {
  list_prot <- fread("/usr/local/data/OVCA_marker_pre-selection_sgo.csv", header=T)
} else {
  cat("WARNING: unsupported sample type! Now set to Melanoma by default")
  list_prot <- fread("/usr/local/data/Melanoma_marker_IDs.csv", header=T)
}

#list_prot <- fread("/usr/local/data/Melanoma_marker_Anja_IDs.csv", header=T)
#list_prot <- fread("/usr/local/data/OVCA_makrer_pre-selection_sgo.csv", header=T)

for(i in 1:length(list_tsv)) {

filename <- gsub("_preprocessed.tsv", "", list_tsv[i])

message("Processing ", i, ": ", filename, "...")

raw <- fread(list_tsv[i], dec=",")

raw$log2fc <- as.numeric(raw$log2fc)

# Force to char
                                        # raw$Gene <- as.character(raw$Gene)
raw$Gene <- sapply( strsplit(raw$Gene, "\\|"), "[[", 1)
    
case <- copy(raw)

fc_threshold <- sort(abs(case$log2fc), decreasing=T)[100]
case_filtered <- case[which(abs(case$log2fc) >= fc_threshold), ]

case_sorted <- case_filtered[order(case_filtered$log2fc), ]

case_sorted[which(case_sorted$Gene %in% cosmic$Gene_Symbol), ]$Gene <- paste0( "***   ", case_sorted[which(case_sorted$Gene %in% cosmic$Gene_Symbol), ]$Gene)

my_col <- rep("gray", dim(case_sorted)[2])
my_col[which(case_sorted$log2fc > 1)] <- "orange"
my_col[which(case_sorted$log2fc < -1)] <- "blue"

my_xmax <- round(max(abs(case_sorted$log2fc), na.rm=T)) + 3
my_xmin <- round(min(case_sorted$log2fc, na.rm=T)) - 3

    write('B', stderr())
    write(my_xmax, stderr())
    write(my_xmin, stderr())

pdf(paste0(filename, "_barplot_top100.pdf"))

barplot(case_sorted$log2fc, names.arg=case_sorted$Gene, horiz=T, las=1, col=my_col, xlab="log2fc", xlim=c(my_xmin, my_xmax), cex.names=0.3, main=paste0(dim(case_sorted)[1], " proteins"))

abline(v=fc_threshold, lty=3)
abline(v=-fc_threshold, lty=3)
dev.off()




markers <- raw[which(raw$Entry %in% list_prot$Uniprot), ]

markers$log2fc <- as.numeric(markers$log2fc)

markers_sorted <- markers[order(markers$log2fc), ]

my_col <- rep("gray", dim(markers_sorted)[1])
my_col[which(markers_sorted$log2fc > 1)] <- "orange"
my_col[which(markers_sorted$log2fc < -1)] <- "blue"

my_xmax <- round(max(abs(markers_sorted$log2fc), na.rm=T)) + 1
my_xmin <- round(min(markers_sorted$log2fc, na.rm=T)) - 1

    write('C', stderr())
    write(markers_sorted$log2fc, stderr())
    write(my_xmax, stderr())
    write(my_xmin, stderr())

    
pdf(paste0(filename, "_barplot_marker.pdf"))

barplot(markers_sorted$log2fc, names.arg=markers_sorted$Gene, horiz=T, las=1, col=my_col, xlab="log2FC", xlim=c(my_xmin, my_xmax), cex.names=0.3, main=paste0(dim(markers_sorted)[1], " marker proteins"))
abline(v=1, lty=3)
abline(v=-1, lty=3)

dev.off()

message("Done.")

}
