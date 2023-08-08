#-----------------------------------------------------
# Title: Phyloseq Object
# Date: 11/08/2023
# Author: Arif Ayman
# Description: Creating the PhyloSeq object
#-----------------------------------------------------

# Set working directory 
setwd("C:/Users/arifa/OneDrive/Desktop/DissertationCode/ReadCounts")

# Read input files
files <- list.files()
index <- grep("_profile_count.txt", files)
files <- as.matrix(files[index])

# Combine unique info from input files
uni_info <- NULL
for (i in 1:length(files)) {
  data_each <- as.matrix(read.csv(files[i], header = TRUE, skip = 5, sep = "\t"))
  info_each <- data_each[, c("X.clade_name", "clade_taxid")]
  uni_info <- rbind(uni_info, info_each)
}
uni_info <- as.matrix(unique(uni_info))

# Combine abundance values from input files
result <- NULL
for (j in 1:length(files)) {
  data_each <- as.matrix(read.csv(files[j], header = TRUE, skip = 5, sep = "\t"))
  loc <- match(uni_info[,"X.clade_name"], data_each[,"X.clade_name"])
  
  value_each <- as.matrix(data_each[loc,"relative_abundance"])
  colnames(value_each) <- paste0("sample_", gsub("_profile_count.txt","",files[j]))
  result <- cbind(result, value_each)
}

# Combine unique info and abundance values and write to file
result <- cbind(uni_info, result)
write.table(result, file = "merged_estimated_relative_abundance_2022.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Clear the environment
rm(list = ls())

# Read merged data file
main <- read.table('merged_estimated_relative_abundance_2022.txt', header = TRUE, sep = "\t")

# Split X.clade_name into taxa ranks
t <- separate(main, X.clade_name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|")
t <- t[!is.na(t$Species), ]
t <- t[!duplicated(t$Species), ]

# Generate OTUs column
x <- 1:dim(t)[1]
OTUs <- paste("OTU", x)
t <- add_column(t, OTUs = OTUs, .before = "clade_taxid")

# Prepare the OTU matrix
OTU <- as.data.frame(t[,-which(names(t) %in% c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Type","clade_taxid"))])
row.names(OTU) <- OTU$OTUs
OTU <- as.data.frame(OTU[,-which(names(OTU) %in% c("OTUs","Species"))])
OTU[is.na(OTU)] <- 0

# Convert OTU matrix to CountMatrix
CountMatrix <- OTU %>% as.matrix()
CountMatrix <- t(CountMatrix)
mode(CountMatrix) <- 'integer'

# Prepare the TAXA matrix
TAX <- as.data.frame(t[,which(names(t) %in% c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species","OTUs"))])
row.names(TAX) <- TAX$OTUs
TAX <- as.data.frame(TAX[,-which(names(TAX) %in% c("OTUs"))])
TaxaMatrix <- TAX %>% as.matrix()

# Read Metadata from file
Metadata <- read_xlsx("BiobankDatabase.xlsx", sheet = 1, range = "A1:B19") %>% as.data.frame()
rownames(Metadata) <- Metadata$SampleID

# Subset metadata to match rownames of CountMatrix
Metadata <- Metadata[rownames(CountMatrix), ]


# Generate OTU table, tax table, and sample data
otuTABLE <- otu_table(CountMatrix, taxa_are_rows = TRUE)
taxTABLE <- tax_table(TaxaMatrix)
sampleDATA <- sample_data(Metadata)

# Convert otu_table to a matrix
otu_mat <- as(otuTABLE, "matrix")

# Transpose the matrix
otu_mat_t <- t(otu_mat)

# Convert the transposed matrix back to otu_table
otuTABLE <- otu_table(otu_mat_t, taxa_are_rows = TRUE)

# Create phyloseq object
ps <- phyloseq(otuTABLE, taxTABLE, sampleDATA)

# Clean taxa names
tax_table(ps)[, colnames(tax_table(ps))] <- gsub(tax_table(ps)[, colnames(tax_table(ps))], pattern = "[a-z]__", replacement = "")

# Extract OTU data, sample data, and taxa data
OTUdata <- abundances(ps)
SampleData <- meta(ps)
TAXAData <- as.data.frame(tax_table(ps)@.Data)

# Save phyloseq object and associated data
save(file = "phyloseq_object_relative_abundance_2022.Rdata", ps, OTUdata, SampleData, TAXAData)

# Read Metadata again
Metadata <- read_xlsx("BiobankDatabase.xlsx", sheet = 1, range = "A1:B19") %>% as.data.frame()
rownames(Metadata) <- Metadata$SampleID
Metadata <- Metadata[rownames(CountMatrix), ]

# Fill empty slots with "n/a"
Metadata[Metadata == ""] <- "n/a"
