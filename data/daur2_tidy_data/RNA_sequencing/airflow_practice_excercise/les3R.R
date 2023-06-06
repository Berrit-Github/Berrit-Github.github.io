#les3
library(DESeq2)
airway_dge <- DESeq(dds)

airway_dge_results <- results(airway_dge)
airway_dge_results


summary(airway_dge_results)

2^0.1
2^1

#What does a LFC of 1 mean? How much did the gene expression then actually change?
# a 2x increase (2^1 =2)

new_airway_dge_results <- results(airway_dge, lfcThreshold = 1 , alpha = 0.05)
summary(new_airway_dge_results)

#How many genes are now considered to be upregulated? And how many genes are downregulated?
#LFC > 1.00 (up)    : 68, 0.32%
#LFC < -1.00 (down) : 36, 0.17%


sign_genes <- airway_dge_results[which(airway_dge_results$padj < 0.05),]


topGene <- sign_genes[which.max(sign_genes$log2FoldChange),]
topGene_name <- rownames(topGene)
topGene_name

geneCounts <- plotCounts(dds, gene = topGene_name, 
                         intgroup = c("treatment"), 
                         returnData = TRUE)

library(tidyverse)
ggplot(geneCounts, aes(x = treatment, y = count)) +
  scale_y_log10() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), 
             size = 3, colour = "darkgreen") +
  xlab("Dexamethasone treatment") +
  ylab("Fragment count") + 
  ggtitle(topGene_name) +
  theme_bw()




appel_bottemGene <- sign_genes[which.min(sign_genes$log2FoldChange),]
appel_bottemGene_name <- rownames(appel_bottemGene)
appel_bottemGene_name

appel_bottemGeneCounts <- plotCounts(dds, gene = appel_bottemGene_name, 
                         intgroup = c("treatment"), 
                         returnData = TRUE)

ggplot(appel_bottemGeneCounts, aes(x = treatment, y = count)) +
  scale_y_log10() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), 
             size = 3, colour = "darkgreen") +
  xlab("Dexamethasone treatment") +
  ylab("Fragment count") + 
  ggtitle(appel_bottemGene_name) +
  theme_bw()


data.frame(airway_dge_results) %>% filter(!is.na(padj))
airway_dge_plotting <- data.frame(airway_dge_results) %>% filter(!is.na(padj))
head(airway_dge_plotting)
airway_dge_plotting <- airway_dge_plotting %>% 
  mutate(signif = if_else(padj < 0.05, "padj < 0.05", "Not significant"))

airway_dge_plotting %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
  geom_point() + 
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") + 
  theme_bw() +
  scale_colour_manual(values = c("grey", "darkgreen"), name = "Significance")+
  annotate("text", x = topGene$log2FoldChange, y = -log10(topGene$padj)*0.8, 
           label = topGene_name, colour = "blue")

#verchil tussen de verschillende condities weten
#wat de counts betekenen 
#data normaalverdeeld tussen alle samples?
#




