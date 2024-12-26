# Load required packages
library(GEOquery)
library(fastqcr)
library(limma)
library(FactoMineR)
library(DESeq2)
# Download data from GEO
GSE107490 <- getGEO("GSE107490", destdir = "data")
GSE61853 <- getGEO("GSE61853", destdir = "data")
# Read expression data for GSE107490
GSE107490_expr <- read.table("/Users/sruthi/Downloads/GSE107490_all_count.txt", skip=78, header=TRUE, row.names=1, sep="\t")
GSE61853_expr <- read.table("/Users/sruthi/Downloads/GSE61853_non-normalized.txt", skip=78, header=TRUE, row.names=1, sep="\t")
# Load required library
geo_id1 <- "GSE107490"  # Adjust this to match your dataset
gse1 <- getGEO(geo_id1)
annotation1 <- pData(gse1[[1]])
head(annotation1)

geo_id2 <- "GSE61853"  # Adjust this to match your dataset
gse2 <- getGEO(geo_id2)
annotation2 <- pData(gse2[[1]])
head(annotation2)

# Let's say you want to remove columns 3 and 5
columns_to_remove <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28)
GSE61853_expr <- GSE61853_expr[, -columns_to_remove]
columns_to_remove1 <- c(4, 5, 6)
GSE107490_expr <- GSE107490_expr[, -columns_to_remove1]
new_column_names1 <- c("CON1", "CON2",	"CON3",	"RAEB1",	"RAEB2",	"RAEB3",	"RCDM1",	"RCDM2", "RCDM3")
colnames(GSE107490_expr)[c(1, 2, 3, 4, 5, 6, 7, 8, 9)] <- new_column_names1

# Let's say you want to rename columns 1 and 3
new_column_names <- c("Control1", "Control2", "Control3", "Control4", "Control5", "Control6", "Control7")
colnames(GSE61853_expr)[c(1, 2, 3, 4, 5, 6, 7)] <- new_column_names

EXP1 <- GSE61853_expr
EXP <- GSE107490_expr
# Define row names to remove
remove1 <- c("GSM2868396", "GSM2868397", "GSM2868398")

# Remove rows from annotation1
annotation1 <- annotation1[!rownames(annotation1) %in% remove1, ]


# Let's say you want to rename columns 1 and 3
new_column_names <- c("MDS1", "MDS1", "MDS1", "MDS2", "MDS2", "MDS2", "MDS2")
colnames(GSE61853_expr)[c(8, 9, 10, 11, 12, 13, 14)] <- new_column_names


# Define the experimental groups in your study
group1 <- factor(c("Control", "Control", "Control", "RAEB", "RAEB", "RAEB", "RCMD", "RCMD", "RCMD"))
# Create the design matrix
design1 <- model.matrix(~0 + group1)
colnames(design1) <- levels(group1)  # Assign column names based on group levels
design1

# Define the experimental groups in your study
group <- factor(c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "RAEB", "RAEB", "RAEB", "RCMD", "RCMD", "RCMD", "RCMD"))
# Create the design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)  # Assign column names based on group levels
design

# Normalize expression data for GSE107490
norm <- normalizeBetweenArrays(EXP, method = "quantile")
# Normalize expression data for GSE61853
norm1 <- normalizeBetweenArrays(EXP1, method = "quantile")

dim(norm)
dim(norm1)


# Fit linear models
fit1 <- lmFit(norm, design = design1)
fit2 <- lmFit(norm1, design = design)

# Define contrasts
contrast1 <- makeContrasts(RAEB_vs_Control = RAEB - Control, RCMD_vs_Control = RCMD - Control, levels = design1)
contrast2 <- makeContrasts(RAEB_vs_Control = RAEB - Control, levels = design)
# Fit contrasts
fit1_contrast <- contrasts.fit(fit1, contrast1)
fit2_contrast <- contrasts.fit(fit2, contrast2)
# Apply eBayes moderation
fit1_ebayes <- eBayes(fit1_contrast)
fit2_ebayes <- eBayes(fit2_contrast)
# Extract top DEGs
topDEGs1 <- topTable(fit1_ebayes, coef = "RAEB_vs_Control", number = Inf)  # for GSE107490
topDEGs2 <- topTable(fit2_ebayes, coef = "RAEB_vs_Control", number = Inf)  # for GSE61853

# Filter top DEGs based on p-value and logFC threshold
topDEGs_filtered1 <- topDEGs1[abs(topDEGs1$logFC) >= 2 & topDEGs1$P.Value < 0.05, ]
topDEGs_filtered2 <- topDEGs2[abs(topDEGs2$logFC) >= 2 & topDEGs2$P.Value < 0.05, ]
print(topDEGs_filtered1)
print(topDEGs_filtered2)


# Install and load the annotation package if not already installed
if (!requireNamespace("illuminaHumanv4.db", quietly = TRUE)) {
  install.packages("illuminaHumanv4.db")
}
library(illuminaHumanv4.db)

# Extract probe IDs from topDEGs_filtered2
probe_ids2 <- rownames(topDEGs_filtered2)

# Retrieve gene symbols for Illumina probe IDs
gene_symbols2 <- select(illuminaHumanv4.db, keys = probe_ids2, columns = "SYMBOL", keytype = "PROBEID")

# Print the gene symbols
print(gene_symbols2)


# Merge gene_symbols2 with topDEGs_filtered2
topDEGs_filtered2_with_symbols <- merge(topDEGs_filtered2, gene_symbols2, by.x = "row.names", by.y = "PROBEID", all.x = TRUE)

# Remove the row names column since it's duplicated
row.names(topDEGs_filtered2_with_symbols) <- NULL

# Print the merged data frame
print(topDEGs_filtered2_with_symbols)

# Remove rows with NA values in the gene_symbols column
topDEGs_filtered2_with_symbols <- topDEGs_filtered2_with_symbols[complete.cases(topDEGs_filtered2_with_symbols$gene_symbols), ]

# Print the updated data frame
print(topDEGs_filtered2_with_symbols)

# Sort the data frames based on absolute log2FC and p-value
topDEGs_filtered1 <- topDEGs_filtered1[order(-abs(topDEGs_filtered1$logFC), topDEGs_filtered1$P.Value), ]
topDEGs_filtered2_with_symbols <- topDEGs_filtered2_with_symbols[order(-abs(topDEGs_filtered2_with_symbols$logFC), topDEGs_filtered2_with_symbols$P.Value), ]

# Subset the top 50 DEGs based on the given thresholds
top50_DEGs_1 <- topDEGs_filtered1[abs(topDEGs_filtered1$logFC) >= 2 & topDEGs_filtered1$P.Value < 0.05, ][1:50, ]
top50_DEGs_2 <- topDEGs_filtered2_with_symbols[abs(topDEGs_filtered2_with_symbols$logFC) >= 2 & topDEGs_filtered2_with_symbols$P.Value < 0.05, ][1:50, ]

# Write the top 50 DEGs data frames to CSV files
write.csv(top50_DEGs_1, file = "top50_DEGs_GSE107490.csv", row.names = FALSE)
write.csv(top50_DEGs_2, file = "top50_DEGs_GSE61853.csv", row.names = FALSE)


# Write the top 50 DEGs data frames to CSV files
write.csv(topDEGs_filtered1, file = "top_DEGs_GSE107490.csv", row.names = FALSE)
write.csv(topDEGs_filtered2_with_symbols, file = "top_DEGs_GSE61853.csv", row.names = FALSE)
# Define thresholds
logFC_threshold <- 2  # Log fold change threshold
p_value_threshold <- 0.05  # p-value threshold

# For GSE107490
upregulated_GSE107490 <- topDEGs_filtered1[topDEGs_filtered1$logFC > logFC_threshold & topDEGs_filtered1$P.Value < p_value_threshold, ]
downregulated_GSE107490 <- topDEGs_filtered1[topDEGs_filtered1$logFC < -logFC_threshold & topDEGs_filtered1$P.Value < p_value_threshold, ]

# For GSE61853
upregulated_GSE61853 <- topDEGs_filtered2[topDEGs_filtered2$logFC > logFC_threshold & topDEGs_filtered2$P.Value < p_value_threshold, ]
downregulated_GSE61853 <- topDEGs_filtered2[topDEGs_filtered2$logFC < -logFC_threshold & topDEGs_filtered2$P.Value < p_value_threshold, ]

# Print the number of upregulated and downregulated genes
cat("Number of upregulated genes in GSE107490:", nrow(upregulated_GSE107490), "\n")
cat("Number of downregulated genes in GSE107490:", nrow(downregulated_GSE107490), "\n")
cat("Number of upregulated genes in GSE61853:", nrow(upregulated_GSE61853), "\n")
cat("Number of downregulated genes in GSE61853:", nrow(downregulated_GSE61853), "\n")
# Assuming 'topDEGs_filtered1' and 'topDEGs_filtered2' are already filtered for significant DEGs


# Extract top 25 upregulated and downregulated genes for GSE107490
top25_upregulated_GSE107490 <- head(upregulated_GSE107490[order(-upregulated_GSE107490$logFC, upregulated_GSE107490$P.Value), ], 25)
top25_downregulated_GSE107490 <- head(downregulated_GSE107490[order(downregulated_GSE107490$logFC, upregulated_GSE107490$P.Value), ], 25)

# Extract top 25 upregulated and downregulated genes for GSE61853
top25_upregulated_GSE61853 <- head(upregulated_GSE61853[order(-upregulated_GSE61853$logFC, upregulated_GSE61853$P.Value), ], 25)
top25_downregulated_GSE61853 <- head(downregulated_GSE61853[order(downregulated_GSE61853$logFC, upregulated_GSE61853$P.Value), ], 25)

# Save these top 25 genes to CSV files
write.csv(upregulated_GSE107490, "top25_upregulated_GSE107490.csv", row.names = FALSE)
write.csv(downregulated_GSE107490, "top25_downregulated_GSE107490.csv", row.names = FALSE)
write.csv(upregulated_GSE61853, "top25_upregulated_GSE61853.csv", row.names = FALSE)
write.csv(downregulated_GSE61853, "top25_downregulated_GSE61853.csv", row.names = FALSE)

