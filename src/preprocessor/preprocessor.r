library(data.table)

fmi_gene <- fread("/usr/local/data/fmi_gene_list.txt", header=F)
uniprot2gene <- fread("/usr/local/data/uniprotEntry2Gene")

#fmi_gene <- fread("/IMSB/nas22/shaow/project/phrt/TpViz/docker-tp-viz/data/fmi_gene_list.txt", header=F)
#uniprot2gene <- fread("/IMSB/nas22/shaow/project/phrt/TpViz/docker-tp-viz/data/uniprotEntry2Gene")

alternatives <- fmi_gene[grepl("_", fmi_gene$V1), ]
alternatives[, gene1 := sapply( strsplit(alternatives$V1, "_"), "[[", 1  ) ]
alternatives[, gene2 := sapply( strsplit(alternatives$V1, "_"), "[[", 2  ) ]

alternatives_long <- melt(alternatives, id.vars="V1")

list_tsv <- system(command="ls *Normalized_Protein_Report.tsv", intern=T)

print(list_tsv)

for(i in 1:length(list_tsv)) {

    filename <- gsub("_Normalized_Protein_Report.tsv", "", list_tsv[i])
    pept_filename <- gsub("_Normalized_Protein_Report.tsv", "_Normalized_Peptides_Report.tsv", list_tsv[i])

    message("Processing ", i, ": ", filename, "...")

    raw <- fread(list_tsv[i], dec=",")

    raw <- merge(uniprot2gene, raw, by.x="Entry", by.y="PG.ProteinAccessions")
    raw <- raw[-which(raw$Gene==""), ]

    raw <- raw[which(grepl("HUMAN", raw$PG.ProteinNames)), ]
    raw$PG.ProteinNames <- gsub("_HUMAN", "", raw$PG.ProteinNames)
                                        #raw <- raw[-which(grepl("Keratin", raw$PG.ProteinNames)), ]
                                        #names(raw)[which( names(raw) == "PG.ProteinNames" )] <- "Gene"

    raw[which(!is.na(alternatives_long$V1[match(raw$Gene, alternatives_long$value)])), ]$Gene <- alternatives_long$V1[match(raw$Gene, alternatives_long$value)][which(!is.na(alternatives_long$V1[match(raw$Gene, alternatives_long$value)]))]

    pept <- fread(pept_filename, dec=",")
    pept[, numPept := length(EG.PrecursorId), by=(PG.ProteinAccessions)]
    
    #after discussion, we decied to remove this filtering and keep all proteins for the downstream analysis.
    #list_prot <- unique(pept[numPept > 2, ]$PG.ProteinAccessions)
    list_prot <- unique(pept$PG.ProteinAccessions)


                                        #confident <- raw[which(raw$PG.ProteinAccessions %in% list_prot), ]
    confident <- raw[which(raw$Entry %in% list_prot), ]
    confident[, log2fc := 0]
    confident$log2fc <- log2(as.numeric(unlist(confident[, (dim(confident)[2]-1), with=F]))) - log2(as.numeric(unlist(confident[, (dim(confident)[2]-2), with=F])))

    confident <- confident[which(!is.na(confident$log2fc)), ]

    write.table(confident, file=paste0(filename, "_preprocessed.tsv"), sep="\t", col.name=T, row.name=F, quote=F)

    message("Done.")

}
