##############################################################################
############ Differential gene expression analysis using Deseq2 ##############
##############################################################################

library(tximport)
library(readr)
library(sva)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

# Arguments
args <- commandArgs(trailingOnly = TRUE)
sample_table <- args[1]
quant_dir <- args[2]
design_formula <- args[3]
design_formula <- gsub(" ", "", design_formula) # Remove the whitespaces to avoid mistakes
design_formula_full <- as.formula(design_formula)

# Function to automatically get the design formula only with the covariates
remove_first_var <- function(design_formula) {
    terms <- strsplit(as.character(design_formula), "\\+")[[1]]
    terms <- gsub("~", "", terms)
    if (length(terms) > 1) {
        reformulated <- paste("~", paste(terms[-1], collapse = " + "))
    } else {
        reformulated <- paste("~", 1)
    }
    return(as.formula(reformulated))
}
design_formula_covariates <- remove_first_var(design_formula)


################################# 1. Create DESeq2 dataset ################################
# Load the sample table (metadata)
sample_table <- read.csv("/home/emcastillo/data/RNAseq/differential_gene_expression/RNA_pheno_data_male.txt", header = TRUE, sep = ",")
head(sample_table)

# Select the variables to take into account, format them (numeric or factor) and scale them (the numeric ones)
sample_table$control <- factor(sample_table$control, levels = c(1, 0), labels = c("Control", "Asthma"))
sample_table$Sex <- factor(sample_table$Sex)
sample_table$original_recruitment_yr.e <- factor(sample_table$original_recruitment_yr.e)
sample_table$child.ethnicity <- factor(
    sample_table$child.ethnicity,
    levels = c("Puerto Rican", "Other Latino"),
    labels = c("Puerto_Rican", "Other_Latino")
)
sample_table$age <- scale(sample_table$age)

# Define files and sample names
quant_dir <- "/home/emcastillo/data/RNAseq/preprocesado/RSEM_transcript_quant"
files <- list.files(quant_dir, pattern = "\\.genes\\.results$", full.names = TRUE, recursive = TRUE)
names(files) <- sub(".genes.results$", "", basename(files))

# Define a function to filter genes from a single file. The function change the value of 'effective_length'
# when is zero or less to 1e-6, as these values are invalid for downstream analysis. This may introduce noise
# in the results. Another alternative is to remove those genes from the samples and then just select the genes
# that are present in all the samples (because tximport needs all the samples to have the same genes)
filter_and_save_files <- function(file) {
    data <- read.table(file, header = TRUE) # Read the file
    data$effective_length[data$effective_length <= 0] <- 1e-6

    # Create a temporary file for the filtered data
    temp_file <- tempfile(fileext = ".txt")
    write.table(data, file = temp_file, quote = FALSE, row.names = FALSE, sep = "\t")

    return(temp_file) # Return the path to the temporary file
}

# Apply the function to all files
filtered_file_paths <- lapply(files, filter_and_save_files)

# Convert the list of filtered file paths to a named vector for tximport
names(filtered_file_paths) <- names(files)
filtered_file_paths <- unlist(filtered_file_paths)

# Use tximport with the paths to the filtered files
txi <- tximport(filtered_file_paths, type = "rsem", txIn = FALSE, txOut = FALSE)

# Convert to DESeq2 dataset
dds <- DESeqDataSetFromTximport(
    txi,
    colData = sample_table,
    design = design_formula_full
)


############################## 2. Filter low-expressed genes ###################################
# Criteria: Keep only rows that have a count of at least 10 for a minimal number of samples
# Find the smollest group size
groups <- table(sample_table$control)
smollest_group_size <- min(groups)

# Pre-filter genes
keep <- rowSums(counts(dds) >= 10) >= smollest_group_size

# Subset the dataset to include only the retained genes
dds <- dds[keep, ]


############################## 3. SVA analysis and covariates selection #############################
# Extract the normalized counts matrix
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized = TRUE)

# Filter low-expressed genes
idx <- rowMeans(dat) > 1
dat <- dat[idx, ]

# Define the full model
mod <- model.matrix(design_formula_full, colData(dds))

# Define the null model
mod0 <- model.matrix(design_formula_covariates, colData(dds))

# Perform the SVA
set.seed(123)
svseq <- svaseq(dat, mod, mod0)

###### Check which SVs are relevant
# First we check for correlation between the SVs and the main varibles. We need to binarize the variable control
control_binary <- as.numeric(colData(dds)$control) - 1 # Asegurarse que quede 0 y 1

# Calculate correlation (point-biserial correlation) of each SVs with control
correlations <- apply(svseq$sv, 2, function(sv) cor(sv, control_binary, method = "pearson"))

# Identify low-correlated SVs, as this means they are not takig away biological variation of the main variable (control)
selected_SVs <- which(abs(correlations) < 0.3)

# Keep only the low-correlated SVs
adjusted_svs <- svseq$sv[, selected_SVs]

# Number of SVs retained
ncol(adjusted_svs)

# Now let's check what % of variance explain each SV, to see how many SV are necessary to explain the majority of the variance
# Total variance (total sum of variances across all genes)
sv_matrix <- svseq$sv
variances <- apply(sv_matrix, 2, var) # Variance of each SV

# Percentage of variability explained by each SV
percent_variance <- variances / sum(variances) * 100

# Cumulative percentage of explained variance
cumulative_variance <- cumsum(percent_variance)

# Plot to decide wich number of SVs to choose
png("/home/emcastillo/data/RNAseq/differential_gene_expression/plots/plot_variance_explained.png")
plot(1:24, cumulative_variance,
    type = "b", pch = 16, xlab = "Número de SVs",
    ylab = "Varianza acumulada (%)", main = "Porcentaje acumulado de varianza explicada"
)
abline(h = 90, col = "red", lty = 2) # Línea para el 90%
abline(h = 95, col = "blue", lty = 2) # Línea para el 95%
dev.off()

##### Add the relevant SVs to the DeseqDataSet and to the design
# In this case, it seems like all SVs could be relevant
# Load the raw dds (without normalization)
load("/home/emcastillo/data/RNAseq/differential_gene_expression/dds_object.RData")

# Copy the original dds object
ddssva <- dds

# Extract SVs and create column names
svs <- svseq$sv
colnames(svs) <- paste0("SV", seq_len(ncol(svs))) # Name the columns as SV1, SV2, ...

# Add the SVs to colData
colData(ddssva) <- cbind(colData(ddssva), as.data.frame(svs))

# Update the design formula
terms <- strsplit(as.character(design_formula), "\\+")[[1]]
terms <- gsub("~", "", terms)
design_update <- paste("~", paste(c(terms[-1], colnames(svs), terms[1]), collapse = " + "))

design(ddssva) <- as.formula(design_update)

# Proceed with DESeq2 analysis
ddssva <- DESeq(ddssva)

# Extract the results
res <- results(ddssva, alpha = 0.05)
res_sig <- res[which(res$padj < 0.05), ]

# Vulcano plot
png("/home/emcastillo/data/RNAseq/differential_gene_expression/plots/Vulcano_plot.png")
EnhancedVolcano(res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "padj",
    title = "Volcano plot - DESeq2 results",
    pCutoff = 0.05, # Adjusted p-value threshold
    FCcutoff = 0.15, # Log2 fold change threshold
    pointSize = 2.0,
    labSize = 3.0,
    colAlpha = 0.75, # Transparency
    legendPosition = "bottom",
    legendLabSize = 10,
    legendIconSize = 3
)
dev.off()

save(ddssva, file = "/home/emcastillo/data/RNAseq/differential_gene_expression/ddssva_object_DESeq_done.RData")
