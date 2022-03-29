#fetching matrix data from JASPAR
library(JASPAR2018)
library(TFBSTools)
library(Biostrings)
library(stringr)
library(dplyr)
library(biomaRt)
library(seqinr)



#USER INPUTS
#first column must be the EnsemblID
#valid csv format example:
#"ENSEMBL","symbol"
#"ENSMUSG00000073700","Klhl21"
#"ENSMUSG00000083774","Gm7180"
#"ENSMUSG00000024805","Pcgf5"

##add path to gene lists (seperated in up and downregulated genes)

upregulatedgenes <- read.csv(upregpath, #PATH1
                header = TRUE, skip = 0, sep = ",")
downregulatedgenes <- read.csv(downregpath, #PATH2
                header = TRUE, skip = 0, sep = ",")

#USER INPUTS
#choose the organism mathcing your data using the right format
#choose the region of interest as number of bps before TSS

speciesbiomart <- "mmusculus_gene_ensembl" #check available species with listDatasets(useMart('ensembl'))
speciesJASPAR <- "Mus musculus" #The species source for the sequences, in Latin (e.g. Homo sapiens)
upstreamlength <- 500

colnames(upregulatedgenes)[1] = "geneid"
colnames(downregulatedgenes)[1] = "geneid"

##NOW RUN THE FULL SCRIPT

#initialize mart
mart <- useEnsembl( "ensembl", dataset = speciesbiomart)

#Download upstreamflank
upregulated_upstreamflank <- getBM(
  attributes = c("external_gene_name",
                 "ensembl_gene_id",
                 "ensembl_gene_id_version",
                 "description", "ensembl_transcript_id",
                 "ensembl_transcript_id_version",
                 "coding_gene_flank",
                 "start_position",
                 "end_position"
  ),
  filters = c("ensembl_gene_id",
              "upstream_flank"
  ),
  values = list(upregulatedgenes$geneid, upstreamlength),
  mart = mart ,checkFilters = FALSE, bmHeader=TRUE
)

downregulated_upstreamflank <- getBM(
  attributes = c("external_gene_name",
                 "ensembl_gene_id",
                 "ensembl_gene_id_version",
                 "description", "ensembl_transcript_id",
                 "ensembl_transcript_id_version",
                 "coding_gene_flank",
                 "start_position",
                 "end_position"
  ),
  filters = c("ensembl_gene_id",
              "upstream_flank"
  ),
  values = list(downregulatedgenes$geneid, upstreamlength),
  mart = mart ,checkFilters = FALSE, bmHeader=TRUE
)

#set options for database search
opts <- list()
opts[["species"]] <- speciesJASPAR

#data retrieval and conversion to PWM
Hs_matrix <- getMatrixSet(JASPAR2018, opts)
Hs_PWM <- toPWM(Hs_matrix, pseudocounts = 0.8)

# delete NA sequences and convert to DNAStringset
upregulated_upstreamflank <- upregulated_upstreamflank[!(upregulated_upstreamflank$`Flank-coding region (Gene)`=="Sequence unavailable"),]
upregulated_stringset <- DNAStringSet(upregulated_upstreamflank$`Flank-coding region (Gene)`,use.names=TRUE)
upregulated_stringset@ranges@NAMES <- with(upregulated_upstreamflank, paste(`Gene stable ID version`,`Gene name`, `Gene description`, `Gene start (bp)`, `Gene end (bp)`, sep = " | "))
downregulated_upstreamflank <- downregulated_upstreamflank[!(downregulated_upstreamflank$`Flank-coding region (Gene)`=="Sequence unavailable"),]
downregulated_stringset <- DNAStringSet(downregulated_upstreamflank$`Flank-coding region (Gene)`,use.names=TRUE)
downregulated_stringset@ranges@NAMES <- with(downregulated_upstreamflank, paste(`Gene stable ID version`, `Gene name`, `Gene description`, `Gene start (bp)`, `Gene end (bp)`, sep = " | "))
#scan sequences
siteset_up <- searchSeq(Hs_PWM, upregulated_stringset, strand = "*")
siteset_down <- searchSeq(Hs_PWM, downregulated_stringset, strand = "*")

#write GFF
GFF_upseq <- writeGFF3(siteset_up)
GFF_downseq <- writeGFF3(siteset_down)

#calculate pvalues
#pvals_up <- pvalues(EN1ko_siteset_up, type = "sampling")
#pvals_down <- pvalues(EN1ko_siteset_down, type = "sampling")

#split GFF
names_upseq <- row.names(GFF_upseq)
GFF_upseq_names <- cbind(names_upseq, GFF_upseq)
GFF_upseq_split <- str_split_fixed(GFF_upseq_names$names_upseq, n = 2,
                                   pattern = ".E")
ID <- as.data.frame(GFF_upseq_split)$V1
GFF_upseq_IDs <- cbind(ID, GFF_upseq)
GFF_upseq_genes <-as.data.frame(str_split_fixed(GFF_upseq$seqname,n = 5,
                                                pattern = " \\| "))
GFF_upseq_IDS_genes <- cbind(GFF_upseq_genes$V2, GFF_upseq_IDs)
colnames(GFF_upseq_IDS_genes)[1] <- "gene"

names_downseq <- row.names(GFF_downseq)
GFF_downseq_names <- cbind(names_downseq, GFF_downseq)
GFF_downseq_split <- str_split_fixed(GFF_downseq_names$names_downseq, n = 2,
                                     pattern = ".E")
ID <- as.data.frame(GFF_downseq_split)$V1
GFF_downseq_IDs <- cbind(ID, GFF_downseq)
GFF_downseq_genes <-as.data.frame(str_split_fixed(GFF_downseq$seqname, n = 5,
                                                  pattern = " \\| "))
GFF_downseq_IDS_genes <- cbind(GFF_downseq_genes$V2, GFF_downseq_IDs)
colnames(GFF_downseq_IDS_genes)[1] <- "gene"

#filter results
GFF_upseq_filtered <- GFF_upseq_IDS_genes %>%
  arrange(desc(score)) %>%
  filter(score> quantile(score, .99))
GFF_downseq_filtered <- GFF_downseq_IDS_genes %>%
  arrange(desc(score)) %>%
  filter(score> quantile(score, .99))

#generate match name column
GFF_upseq_filtered$matchname <- paste0(GFF_upseq_filtered$gene, "-",
                                       GFF_upseq_filtered$ID)
GFF_downseq_filtered$matchname <- paste0(GFF_downseq_filtered$gene, "-",
                                         GFF_downseq_filtered$ID)



#mycode by Dominik
#Upseq
GFF_upseq_filtered_dupeless <- GFF_upseq_filtered[!duplicated(GFF_upseq_filtered$matchname),]
GFF_upseq_filtered_dupeless <- GFF_upseq_filtered_dupeless[,c(2,11,12)]
GFF_upseq_filtered_dupeless_2 <- GFF_upseq_filtered_dupeless[!duplicated(GFF_upseq_filtered_dupeless$ID),]
strsplit <- sapply(GFF_upseq_filtered_dupeless_2$attributes, function(x) strsplit(x, ";")[[1]], USE.NAMES=FALSE)
GFF_upseq_filtered_dupeless_2$TF <- substring(strsplit[1,],4)
GFF_upseq_filtered_dupeless_2$class <- substring(strsplit[2,],7)

GFF_upseq_filtered_dupeless_2$attributes <- NULL
GFF_upseq_filtered_dupeless_2$matchname <- NULL
tabelleup <- as.data.frame(table(GFF_upseq_filtered_dupeless$ID))
GFF_upseq_filtered_dupeless_2$count <- tabelleup$Freq[match(unlist(GFF_upseq_filtered_dupeless_2$ID), tabelleup$Var1)]

GFF_upseq_filtered_dupeless_2$ratio <- GFF_upseq_filtered_dupeless_2$count / length(upregulated_upstreamflank$`Gene name`) * 100

#to which upregulated Genes does the factor bind

GFF_upseq_filtered_dupeless_genelist <- data.frame(GFF_upseq_filtered_dupeless$ID)
names(GFF_upseq_filtered_dupeless_genelist) <- c("ID")
strsplit2 <- sapply(GFF_upseq_filtered_dupeless$attributes, function(x) strsplit(x, ";")[[1]], USE.NAMES=FALSE)
GFF_upseq_filtered_dupeless_genelist$TF <- substring(strsplit2[1,],4)
GFF_upseq_filtered_dupeless_genelist$class <- substring(strsplit2[2,],7)
GFF_upseq_filtered_dupeless_genelist$gene <- GFF_upseq_filtered_dupeless$matchname

#Downseq
GFF_downseq_filtered_dupeless <- GFF_downseq_filtered[!duplicated(GFF_downseq_filtered$matchname),]
GFF_downseq_filtered_dupeless <- GFF_downseq_filtered_dupeless[,c(2,11,12)]
GFF_downseq_filtered_dupeless_2 <- GFF_downseq_filtered_dupeless[!duplicated(GFF_downseq_filtered_dupeless$ID),]
strsplit <- sapply(GFF_downseq_filtered_dupeless_2$attributes, function(x) strsplit(x, ";")[[1]], USE.NAMES=FALSE)
GFF_downseq_filtered_dupeless_2$TF <- substring(strsplit[1,],4)
GFF_downseq_filtered_dupeless_2$class <- substring(strsplit[2,],7)

GFF_downseq_filtered_dupeless_2$attributes <- NULL
GFF_downseq_filtered_dupeless_2$matchname <- NULL
tabelledown <- as.data.frame(table(GFF_downseq_filtered_dupeless$ID))
GFF_downseq_filtered_dupeless_2$count <- tabelledown$Freq[match(unlist(GFF_downseq_filtered_dupeless_2$ID), tabelledown$Var1)]

GFF_downseq_filtered_dupeless_2$ratio <- GFF_downseq_filtered_dupeless_2$count / length(downregulated_upstreamflank$`Gene name`) * 100

#to which downregulated Genes does the factor bind

GFF_downseq_filtered_dupeless_genelist <- data.frame(GFF_downseq_filtered_dupeless$ID)
names(GFF_downseq_filtered_dupeless_genelist) <- c("ID")
strsplit2 <- sapply(GFF_downseq_filtered_dupeless$attributes, function(x) strsplit(x, ";")[[1]], USE.NAMES=FALSE)
GFF_downseq_filtered_dupeless_genelist$TF <- substring(strsplit2[1,],4)
GFF_downseq_filtered_dupeless_genelist$class <- substring(strsplit2[2,],7)
GFF_downseq_filtered_dupeless_genelist$gene <- GFF_downseq_filtered_dupeless$matchname

#export results
write.table(GFF_upseq_filtered_dupeless_2, "GFF_upregulated.csv", sep = ",",
            row.names = FALSE)
write.table(GFF_upseq_filtered_dupeless_genelist, "GFF_upregulated.csv", sep = ",",
            row.names = FALSE)
write.table(GFF_downseq_filtered_dupeless_2, "GFF_downregulated.csv", sep = ",",
            row.names = FALSE)
write.table(GFF_downseq_filtered_dupeless_genelist, "GFF_downregulated.csv", sep = ",",
            row.names = FALSE)


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")





