#Airway file

# Load the required libraries
library(Rsubread)
library(tidyverse)
air_read_counts <- read_rds('/home/daur2/rnaseq/rnaseq_airway/counts/read_counts.rds')

# Obtain the count matrix
air_count_matrix <- air_read_counts$counts
# Import the sample data
air_metadata <- read_csv("/home/daur2/rnaseq/rnaseq_airway/airway_sampledata.csv")
# Convert the metadata to dataframe object
air_metadata <- as.data.frame(air_metadata)

# Add rownames to the metadata dataframe
rownames(air_metadata) <- paste0(air_metadata$Run, ".bam")

# Check if column names of count table are the same as row names of metadata object
colnames(air_count_matrix) == rownames(air_metadata)

# Create a column specifying the test condition
air_metadata <- air_metadata %>% mutate(treatment = str_replace(dex, "trt", "treated"))
air_metadata$treatment <- air_metadata$treatment %>% factor(levels = c("untreated", "treated"))
library(DESeq2)

air_dds <- DESeqDataSetFromMatrix(
  countData = air_count_matrix,
  colData = air_metadata, 
  design = ~ treatment)

airway_dge <- DESeq(air_dds)

airway_dge_results <- results(airway_dge, lfcThreshold = 1, alpha = 0.05)

summary(airway_dge_results)

