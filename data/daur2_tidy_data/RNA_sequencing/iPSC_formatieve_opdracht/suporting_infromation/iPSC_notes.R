
library(Rsubread)


bam_dir <- "/home/daur2/rnaseq/rnaseq_ipsc/bam/"

# Create object with output dir for count tables
counts_dir <- "/home/daur2/rnaseq/rnaseq_ipsc/counts/"

# Create vector with names of bam files
bam_files <- list.files(bam_dir, pattern = ".*\\.bam$", full.names = TRUE)

# Count the reads per gene using the in-built NCBI RefSeq annotations
ipsc_counts <- featureCounts(
  
  files = bam_files,
  annot.inbuilt = "hg38",
  useMetaFeatures = TRUE,
  strandSpecific = 1,
  isPairedEnd = TRUE, 
  countReadPairs = TRUE, 
  nthreads = 10,
  reportReadsPath = 
)

library(tidyverse)

?featureCounts

ipsc_counts <- read_rds("/home/daur2/rnaseq/rnaseq_ipsc/counts/read_counts.rds")
ispc_counts_stats <- ispc_counts$stat

rownames(ispc_counts_stats) <- ispc_counts_stats$Status
ispc_counts_stats$Status <- NULL

ispc_counts_stats_tibble <- ispc_counts_stats %>% 
  t %>% 
  as_tibble() %>% 
  mutate(bamfile=colnames(ispc_counts_stats)) %>%
  mutate(Total=colSums(ispc_counts_stats)) %>%
  mutate(perc_assigned = Assigned/Total*100)

view(ispc_counts_stats_tibble)

ispc_counts_stats_tibble %>% ggplot(aes(x = bamfile , y = perc_assigned)) + 
  geom_col() +
  labs(title = "percentage assigned reads" , x = "bamfile" , y = "percentage assigned reads") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_cartesian(ylim = c(1,100))
  

BiocManager::install("DESeq2")


#countmatrix
ispc_count_matrix <- ispc_counts$counts

#dataframe with details (metadata)
ispc_metadata <- read_csv('/home/daur2/rnaseq/rnaseq_ipsc/ipsc_sampledata.csv')
ispc_metadata <- as.data.frame(ispc_metadata)
rownames(ispc_metadata) <- paste0(ispc_metadata$Run, ".bam")
# check of de namen van de countmatrix over een komen met die van de dataframe
colnames(ispc_count_matrix) == rownames(ispc_metadata)

#experimental design
view(ispc_metadata)
library(DESeq2)
                                                          
# Create a column specifying the test condition
ispc_metadata <- ispc_metadata %>% mutate(cell_lijn = str_replace(Cell_type, "Skin derived fibroblast", "fibroblast"))
ispc_metadata$cell_lijn <- ispc_metadata$cell_lijn %>% factor(levels = c("fibroblast", "iPSC"))


#create 
ispc_dds <- DESeqDataSetFromMatrix(
  countData = ispc_count_matrix,
  colData = ispc_metadata, 
  design = ~ cell_lijn 
)

ispc_dds


##############

ispc_dds_normalized <- rlog(ispc_dds)


ispc_pca <- ispc_dds_normalized %>% assay() %>% t() %>% prcomp()

ispc_pca_summary <- summary(ispc_pca)$importance

ispc_pca_plotting <- cbind(ispc_metadata, ispc_pca$x)

# Obtain the percentages of variation covered by PC1 and PC2
PC1_var <- round(ispc_pca_summary["Proportion of Variance", "PC1"]*100, digits = 1)
PC2_var <- round(ispc_pca_summary["Proportion of Variance", "PC2"]*100, digits = 1)

# Plot PC1 vs PC2
ggplot(ispc_pca_plotting) + 
  geom_point(aes(x=PC1, y=PC2, color = Cell_type, shape =source_name ), size = 5) +
  ggtitle("PCA for airway study") +
  xlab(paste0("PC1 (", PC1_var, "%)")) +
  ylab(paste0("PC2 (", PC2_var, "%)")) +
  theme_bw()

#####
ispc_dds_norm <- assay(ispc_dds_normalized)

ispc_cor <- cor(ispc_dds_norm)

library(pheatmap)
pheatmap(ispc_cor, annotation = ispc_metadata["source_name"])

#######

ispc_dge <- DESeq(ispc_dds)
results_ispc_dge <- results(ispc_dge, lfcThreshold = 1, alpha = 0.05)
results_ispc_dge
summary(results_ispc_dge)


####NOT USED

head(results_ispc_dge)
drf
padj <0.05
lfcSE >1

which(results_ispc_dge$padj <0.05) %>% which(results_ispc_dge$log2FoldChange >1)

results_ispc_dge[ which(results_ispc_dge$padj <0.05) | which(results_ispc_dge$log2FoldChange  >1),    ]
sign_1 <- results_ispc_dge[ which(results_ispc_dge$padj <0.05),]

sign_1 <- sign_1[ which(sign_1$lfcSE >1),    ]
head(sign_1)

highP_highLFC_names <- rownames(sign_1)

highP_highLFC_names

#####
library(tidyverse)

ispc_dge_plotting <- data.frame(results_ispc_dge) %>% filter(!is.na(padj))

ispc_dge_plotting <- ispc_dge_plotting %>% 
  mutate(valid = if_else(padj < 0.05 & log2FoldChange  >1 , "validP", "not_valid"))


head(ispc_dge_plotting)
ispc_dge_plotting %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = valid)) +
  geom_point() + 
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") + 
  theme_bw() + scale_colour_manual(values = c("grey", "darkblue"), name = "Significance")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept =  1, linetype = "dashed")
########
high_LFC_gene_names <- results_ispc_dge[ which(results_ispc_dge$padj <0.05),]

high_LFC_rownames <- rownames(high_LFC_gene_names[order(high_LFC_gene_names$log2FoldChange )[1:15],])

count_values_highLFC <- assay(ispc_dds)[high_LFC_rownames,]
colnames(count_values_highLFC)<- colData(ispc_dds)$source_name

pheatmap(count_values_highLFC, show_rownames = TRUE)


####### les 4A


head(count_values_highLFC)

BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

new_heath_map_data <- count_values_highLFC %>% data.frame()

new_heath_map_data<- mutate(new_heath_map_data, entrized = rownames(new_heath_map_data))

new_heath_map_data$symbol <- mapIds(org.Hs.eg.db,
                                    keys = new_heath_map_data$entrized,
                                    column = "SYMBOL",
                                    keytype = "ENTREZID",
                                    multiVals = "first")

new_heath_map_data_test <- new_heath_map_data
rownames(new_heath_map_data_test) <- new_heath_map_data$symbol
new_heath_map_data_test

pheatmap(new_heath_map_data_test[1:8] )


#####function 4B
head(results_ispc_dge)
#if upregulated = true do >
#if upregulated = false do <
library(GOstats)
library(GO.db)


gotermAnalysis <- function(DEseq_results, upregulated, LFC , Pvalue){
  all_names <- DEseq_results %>% data.frame() %>% rownames()
 
         if(upregulated == TRUE) { DEseq_new_filter<- DEseq_results %>% data.frame() %>% filter(log2FoldChange > LFC, padj < Pvalue) %>%rownames()
           } else if (upregulated == FALSE){ DEseq_new_filter<- DEseq_results %>% data.frame() %>%  filter(log2FoldChange < -LFC, padj < Pvalue) %>%rownames()   
    } else {return("NO")}
  
  test_object <- methods::new("GOHyperGParams",
                     geneIds = DEseq_new_filter,
                     universeGeneIds = all_names, 
                     annotation = "org.Hs.eg.db", 
                     ontology = "BP", 
                     pvalueCutoff = 1,
                     testDirection = "over")
  goterm_analysis <- hyperGTest(test_object)
  summary(goterm_analysis)
  
  
}
gotermAnalysis(results_ispc_dge, upregulated = FALSE , 1 , 0.05)


library("org.Hs.eg.db")


BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
library(tidyverse)
library(GO.db)

######
goterm_up_analysis_results <- gotermAnalysis(results_ispc_dge, upregulated = TRUE , 1 , 0.05)

goterm_up_analysis_results$padj <- p.adjust(goterm_up_analysis_results$Pvalue, method = "BH")

goterm_up_analysis_results <- goterm_up_analysis_results %>% filter(Count > 5) %>% filter(Count < 500)


goterm_up_analysis_top30 <- goterm_up_analysis_results[order(goterm_up_analysis_results$padj)[1:30],]
goterm_up_analysis_top30$Term <- factor(goterm_up_analysis_top30$Term, 
                                     levels = goterm_up_analysis_top30$Term[
                                       order(goterm_up_analysis_top30$padj, decreasing = TRUE)])

goterm_up_analysis_top30 %>% ggplot(aes(x = Term, y = -log10(padj))) +
  geom_point() +
  coord_flip() +
  ylab(expression(-log[10](adjusted~italic(P)~value))) + 
  xlab("GO terms") +
  ggtitle(" GO terms for upregulated genes") +
  theme_bw()


#######

goterm_down_analysis_results <- gotermAnalysis(results_ispc_dge, upregulated = FALSE , 1 , 0.05)

goterm_down_analysis_results$padj <- p.adjust(goterm_down_analysis_results$Pvalue, method = "BH")

goterm_down_analysis_results <- goterm_down_analysis_results %>% filter(Count > 5) %>% filter(Count < 500)


goterm_down_analysis_top30 <- goterm_down_analysis_results[order(goterm_down_analysis_results$padj)[1:30],]
goterm_down_analysis_top30$Term <- factor(goterm_down_analysis_top30$Term, 
                                        levels = goterm_down_analysis_top30$Term[
                                          order(goterm_down_analysis_top30$padj, decreasing = TRUE)])

goterm_down_analysis_top30 %>% ggplot(aes(x = Term, y = -log10(padj))) +
  geom_point() +
  coord_flip() +
  ylab(expression(-log[10](adjusted~italic(P)~value))) + 
  xlab("GO terms") +
  ggtitle(" GO terms for downregulated genes") +
  theme_bw()

##########

?feauturecounts
?featureCounts()
high_LFC_gene_names[order(abs(high_LFC_gene_names$log2FoldChange ))[1:15],]

