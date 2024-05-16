# https://www.biostars.org/p/343138/

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("pheatmap")


library(DESeq2)
library(apeglm)
library("ggplot2")
library("pheatmap")
library("genefilter")


# Path to the HTSeq output data.
directory <- "data/RNA_s/"

# Using files from specified directory that contain "ERR" in their file name.
sampleFiles <- grep("ERR", list.files(directory), value=TRUE)
sampleNames <- sub(".*_(ERR[0-9]+)_.*", "\\1", sampleFiles)
sampleCondition <- sub("_.*", "", sampleFiles)

sampleFiles           # Filename of sample
sampleCondition       # Where is sample coming from BH or Serum

sampleTable <- data.frame(sampleName = sampleNames, 
                          fileName   = sampleFiles,
                          condition  = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable

# Preparing DESeq data set for analyzing differences between the different conditions {BH, Serum} using data from specified files
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)
# Contains raw count data and experimental setup information
ddsHTSeq



###### Differential expression analysis ######

# The standard differential expression analysis steps are wrapped into a single function, DESeq.
?DESeq

# Performs differential expression analysis on the ddsHTSeq.
#   Estimation of size factors : adjusts for differences in the sequencing depth across samples
#     Differences in how much RNA is sequenced can affect the number of reads that are mapped to genes.
#     They need to be normalized to make fair comparisons across samples
dds <- DESeq(ddsHTSeq)


# Results table : has log2 fold changes, p values, adjusted p values
res <- results(dds)
res

summary(res)

# Intercept and condition_Serum_vs_BH
resultsNames(dds)


####### Shrinkage ######

# Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes
resLFC <- lfcShrink(dds, coef="condition_Serum_vs_BH", type="apeglm")
resLFC


####### Plot counts #######

# Examine the counts of reads for a single gene across the groups
# Normalizes counts by the estimated size factors
# Counts grouped by the variables in intgroup.
# Gene that had the smallest p-value from the results table. 
#   Can be selected by row name or by numeric index
plotCounts(dds, gene=which.min(resLFC$padj), intgroup="condition")


# Customizes plotting. 
#   returnData=TRUE :: the function should return data frame for plotting the ggplot
# signGene is a table that consists of read counts and conditions for the gene with the lowest padj value.
signGene <- plotCounts(dds, gene=which.min(resLFC$padj), intgroup="condition", returnData=TRUE)

ggplot(signGene, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# Columns of res
mcols(resLFC)$description



####### Logarithmic transformations #######

# Performs a regularized logarithm transformation of the count data.
#   Stabilizes variance across the range of mean values - helpful in the methods that assume constant variance.
#   blind=FALSE the transformation should take into account the differences across conditions.
#     If TRUE the transformation would be done assuming that all samples are form the same experimental condition.
rld <- rlog(dds, blind = FALSE)

# The matrix of the transformed data from rld
assay(rld)
# Shows top 3 rows from the matrix  
head(assay(rld), 3)

colData(rld)


####### Heat maps #######

### Heat map of all genes

allGenes <- order(rowVars(assay(rld)), decreasing = TRUE)
matAll  <- assay(rld)[ allGenes, ]

total_columns <- ncol(matAll)
half_column_index <- total_columns %/% 2 
H1 <- matAll[, 1:half_column_index]
H2 <- matAll[, (half_column_index + 1):total_columns]
H1 <- H1 - rowMeans(H1)
H2 <- H2 - rowMeans(H2)
matAll <- cbind(H1, H2)

anno <- as.data.frame(colData(rld)[,c("condition")])
colnames(anno) <- "condition"
rownames(anno) <- colnames(matAll)
pheatmap(matAll, cluster_rows=FALSE, cluster_cols=FALSE, annotation_col = anno, main = "Heat map of all genes")


### Heat map of significant genes

topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50)
matSig  <- assay(rld)[ topVarGenes, ]

total_columns <- ncol(matSig)
half_column_index <- total_columns %/% 2 
H1 <- matSig[, 1:half_column_index]
H2 <- matSig[, (half_column_index + 1):total_columns]
H1 <- H1 - rowMeans(H1)
H2 <- H2 - rowMeans(H2)
matSig <- cbind(H1, H2)

anno <- as.data.frame(colData(rld)[,c("condition")])
colnames(anno) <- "condition"
rownames(anno) <- colnames(matSig)

#show_rownames=FALSE
pheatmap(matSig, cluster_rows=FALSE, cluster_cols=FALSE, annotation_col = anno, main = "Heat map of significant genes")



### Heat map of significant genes based on filter from the paper

# genes that show significant effect
sig001 <- subset(res, padj < 0.001 & (abs(log2FoldChange) > 2 | abs(log2FoldChange) < 0.5))
sig001


sg_filter <- rownames(sig001)
sg_data <- assay(rld)[sg_filter, ]
sg_data

total_columns <- ncol(sg_data)
half_column_index <- total_columns %/% 2 
H1 <- sg_data[, 1:half_column_index]
H2 <- sg_data[, (half_column_index + 1):total_columns]
H1 <- H1 - rowMeans(H1)
H2 <- H2 - rowMeans(H2)
sg_data <- cbind(H1, H2)

anno <- as.data.frame(colData(rld)[,c("condition")])
colnames(anno) <- "condition"
rownames(anno) <- colnames(sg_data)

#show_rownames=FALSE
pheatmap(sg_data, cluster_rows=FALSE, cluster_cols=FALSE, annotation_col = anno, main = "Heat map of significant genes")




### Heat map of significant genes based on filter from the paper

# genes that show significant effect
sig001 <- subset(res, padj < 0.001 & (abs(log2FoldChange) > 2 | abs(log2FoldChange) < 0.5))
sig001

upReg <- rownames(sig001[order(-sig001$log2FoldChange), ])[1:15]
downReg <-rownames(sig001[order(sig001$log2FoldChange), ])[1:15]
sig001Reg <- c(upReg, downReg)

#show_rownames=FALSE
pheatmap(assay(rld)[sig001Reg, ], cluster_rows=FALSE, cluster_cols=FALSE, annotation_col = anno, main = "Heat map of significant genes")





# Load the annotation data
prokka_data <- read.table("data/spades_prokka.tsv", sep="\t", header=TRUE, quote="")

# Subset CDS data (as that is what we counted)
prokka <- prokka_data[prokka_data$ftype == "CDS",]

# Filtering top 30 most significant up and down regulated genes.
filtered_prokka <- prokka[prokka$locus_tag %in% sig001Reg, ]
filtered_prokka

regulated_prokka <- filtered_prokka %>%
  mutate(Regulated = case_when(
    locus_tag %in% upReg ~ "Up regulated",
    locus_tag %in% downReg ~ "Down regulated"
  ))
regulated_prokka


results_dir <- "data/results/"
results_path <- file.path(results_dir, "12_Differential_expression_sign30genes.csv")
write.csv(regulated_prokka, results_path, row.names = FALSE)

