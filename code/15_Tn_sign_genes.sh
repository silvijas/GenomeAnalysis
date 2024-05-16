##############################################################################################
# Transposon sequencing (Tn-Seq) data is used to determine genes essential genes for growth  #
#   in human serum. Fitness is a measure of how essential is a gene for a growth in a given  #
#   environment (in this case human serum).                                                  #
##############################################################################################

# awk '{$5=$5-int(($5-$4+1)*0.1); print $0}' spades_prokka_noSeq.gff > spades_prokka_noSeq_90perc.gff

# Load necessary library
library(dplyr)
library(tidyr)
library(purrr)

# Input files are the HTSeq output data.
directory <- "data/Tn_s_90_egg/"
file_list <- grep("ERR", list.files(directory), value=TRUE)

# Data from the input file is merged based on the gene id.
merged_data <- lapply(file_list, function(file) {
  dat <- read.table(file, header = FALSE)
  
  sample_id <- sub(".*_(ERR[0-9]+)_.*", "\\1", file)
  
  names(dat)[2] <- sample_id
  names(dat)[1] <- "Genes"
  
  dat
}) %>% reduce(full_join, by = "Genes")

merged_data


# Genes are set to be rownames
rownames(merged_data) <- merged_data$Gene
mat <- merged_data[, -1]

# Removing rows where all elements are zeros
mat_nz <- mat[rowSums(mat) != 0, ]
mat_nz

# Calculating sum of read counts per sample
reads_per_sample <- colSums(mat_nz)


# Load the annotation data
prokka_data <- read.table("data/spades_prokka.tsv", sep="\t", header=TRUE, quote="")

# Subset CDS data (as that is what we counted)
prokka <- prokka_data[prokka_data$ftype == "CDS",]

# Select the 'locus_tag' and 'length_bp' columns
prokka <- prokka[, c("locus_tag", "length_bp")]

# Extracting lengths of a genes for which we have some read counts
gene_lengths <- prokka$length_bp[match(rownames(mat_nz), prokka$locus_tag)]


# Normalization with gene length (short genes will have less reads mapped to it)
norm_gl <- sweep(mat_nz, 1, gene_lengths, FUN="/") 

# Normalization with number of reads per sample (we have not sampled equally amount of times)
norm_rps <- sweep(norm_gl, 2, reads_per_sample, FUN="/")

# RPKM score is (read counts * 10**6) / (number of reads per sample * gene length)
scores <- norm_rps * 10**6
scores


# Calculating p values and fitness
p_values <- numeric(nrow(scores))
fitness  <- numeric(nrow(scores))

for (i in 1:nrow(scores)) {
  BHI_v <- c(scores[i, 1], scores[i, 2], scores[i, 3])
  HS_v <- c(scores[i, 4], scores[i, 5], scores[i, 6])
  test_result <- t.test(BHI_v, HS_v, paired=FALSE)
  p_values[i] <- test_result$p.value
  fitness[i]  <- mean(BHI_v) - mean(HS_v)
}

p_values
fitness

# Adjusting p-values for multiple testing
adjusted_p_values <- p.adjust(p_values, method = "fdr")
adjusted_p_values

results <- cbind(scores, p_values, adjusted_p_values, fitness)
results_df <- data.frame(results)

# Significant results
#results_sign <- results_df[results_df$adjusted_p_values < 0.05,]
results_sign <- results_df[results_df$p_values < 0.05,]
results_sign


# Significant genes
sign_genes <- rownames(results_sign["p_values"])
sign_genes

# Extracting information from annotations for significant genes
sign_genes_ann <- prokka_data[prokka_data$locus_tag %in% sign_genes, ]
sign_genes_ann



results_sign_df <- as.data.frame(results_sign)
results_sign_df$gene <- rownames(results_sign_df)
results_sign_df <- results_sign_df[, c("gene", "p_values", "fitness")] 

final <- merge(prokka_data, results_sign_df, by.x = "locus_tag", by.y = "gene", all.y = TRUE)
final <- final[order(final$p_value), ]
final

results_dir <- "data/results/"
results_path <- file.path(results_dir, "15_TnSeq_essential_genes_90_egg.csv")
write.csv(final, results_path, row.names = FALSE)



##############################################################################################
# If a gene has negative fitness, this generally means that the disruption of this gene has  #
#    led to a decrease in growth or survival rate compared to the control.                   #
#    The gene may be important for optimal growth or survival in the tested condition.       #
##############################################################################################






