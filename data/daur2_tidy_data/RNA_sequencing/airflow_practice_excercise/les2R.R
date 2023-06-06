
library(Rsubread)
library(tidyverse)
BiocManager::install("DESeq2")



bam_dir <- "./rnaseq_airway/bam/"
counts_dir <- "./rnaseq_airway/counts/"
bam_files <- list.files(bam_dir, pattern = ".*\\.bam$", full.names = TRUE)

read_counts <- featureCounts(
  
  files = bam_files,
  annot.inbuilt = "hg38",
  useMetaFeatures = TRUE,
  strandSpecific = 0,
  isPairedEnd = TRUE, 
  countReadPairs = TRUE, 
  nthreads = 10
)
View(read_counts)
read_counts

?featureCounts


read_count_new <-read_rds('/home/daur2/rnaseq/rnaseq_airway/counts/read_counts.rds')

str(read_count_new)
stat_table <-read_count_new$stat
stat_table
view(stat_table)
 new_stat_table <- stat_table %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(bamfile=colnames(stat_table))

 colnames(new_stat_table) <- c("Assigned", "Unassigned_Unmapped" ,"Unassigned_Read_Type" ,"Unassigned_Singleton" ,"Unassigned_MappingQuality" , "Unassigned_Chimera", "Unassigned_FragmentLength" ,
     "Unassigned_Duplicate" , "Unassigned_MultiMapping" ,"Unassigned_Secondary" , "Unassigned_NonSplit" , "Unassigned_NoFeatures" , "Unassigned_Overlapping_Length" , "Unassigned_Ambiguity" , "bamfiles")
view(new_stat_table)
new_stat_table <- new_stat_table[-1,]
parse_number(new_stat_table$Assigned , new_stat_table$Unassigned_Unmapped)
as_tibble(new_stat_table)
new_stat_table$Unassigned_Ambiguity <- parse_number(new_stat_table$Unassigned_Ambiguity)

stat_table_plot <- new_stat_table %>% mutate(percentage_assigend = round(Assigned/ (Assigned+Unassigned_NoFeatures+Unassigned_Ambiguity+Unassigned_Unmapped) *100 , digits = 2))

view(stat_table_plot)



stat_table_plot %>% ggplot(aes(x = bamfiles, y = percentage_assigend)) +
  geom_col() +
  ggtitle("Proportion of assigned fragments for each sample") +
  xlab("RNA-seq sample") +
  ylab("Percentage of assigned fragments") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_cartesian (ylim =c(0,100))


count_matrix <- read_counts$counts

metadata <- read_csv('/home/daur2/rnaseq/rnaseq_airway/airway_sampledata.csv')
metadata
metadata <- as.data.frame(metadata)
rownames(metadata) <- paste0(metadata$Run, ".bam")
head(metadata)

colnames(count_matrix) == rownames(metadata)
colnames(count_matrix) == rownames(metadata)
colnames(count_matrix)

metadata <- metadata %>% mutate(treatment = str_replace(dex, "trt", "treated"))
metadata$treatment <- metadata$treatment %>% factor(levels = c("untreated", "treated"))
levels(metadata$treatment)
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata, 
  design = ~ treatment)
dds

quantile(count_matrix)

count_matrix
metadata
#bijna de helft van de genen komen bijna niet tot expressie

quantile(count_matrix[, "SRR1039521.bam" ])
head(count_matrix)

dds_normalized <- rlog(dds)

dds_normalized


pca <- dds_normalized %>%  assay() %>% t() %>% prcomp()
pca

pca_summary <- summary(pca)$importance
pca_summary

pca_plotting <- cbind(metadata, pca$x)
view(pca_plotting)

PC1_var <- round(pca_summary["Proportion of Variance", "PC1"]*100, digits = 1)
PC2_var <- round(pca_summary["Proportion of Variance", "PC2"]*100, digits = 1)

ggplot(pca_plotting) + 
  geom_point(aes(x=PC1, y=PC2, color = treatment, shape = cell_line), size = 5) +
  ggtitle("PCA for airway study") +
  xlab(paste0("PC1 (", PC1_var, "%)")) +
  ylab(paste0("PC2 (", PC2_var, "%)")) +
  theme_bw()

PC3_var <- round(pca_summary["Proportion of Variance", "PC3"]*100, digits = 1)
PC4_var <- round(pca_summary["Proportion of Variance", "PC4"]*100, digits = 1)


ggplot(pca_plotting) + 
  geom_point(aes(x=PC1, y=PC3, color = treatment, shape = cell_line), size = 5) +
  ggtitle("PCA for airway study") +
  xlab(paste0("PC1 (", PC1_var, "%)")) +
  ylab(paste0("PC2 (", PC3_var, "%)")) +
  theme_bw()


dds_normalized_matrix <- assay(dds_normalized)
dds_normalized_matrix

airway_cor <- cor(dds_normalized_matrix)    
airway_cor

library(pheatmap)

pheatmap(airway_cor,annotation = metadata["treatment"] )

pheatmap(airway_cor,annotation = metadata["treatment"] , cluster_rows = FALSE, cluster_cols = FALSE)






















