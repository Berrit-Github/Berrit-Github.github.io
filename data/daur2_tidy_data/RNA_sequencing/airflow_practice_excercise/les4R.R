
BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
library(tidyverse)

columns(org.Hs.eg.db)
help("SYMBOL")

air_top10_genes <- airway_dge_results[order(airway_dge_results$padj)[1:10],] %>% data.frame()

air_top10_genes <- air_top10_genes %>% mutate(entrezid = rownames(air_top10_genes))


air_top10_genes$symbol <- mapIds(org.Hs.eg.db,
                             keys = air_top10_genes$entrezid,
                             column = "SYMBOL",
                             keytype = "ENTREZID",
                             multiVals = "first")
head(air_top10_genes)

?mapIds
# een data selectie methode met vele opties
#The function looks up for each gene identifier (the keys option can be used to specify the gene identifiers of interest, 
#and the keytype option to specify the type of gene identifier) the annotation of interest (which is specified by the column option). 
#In case of multiple matching annotations (in this case, one Entrez identifier matching more than one gene symbols),
#the multiVals = "first" option specifies that only the first match should be reported.

#nee omdat we de een deel vand e cellen met een glucocorticoide hebben gestimuleerd. dit zorgt er dus voor dat deze cellen meer tot expressie brengen

air_top10_genes$CHR <- mapIds(org.Hs.eg.db,
                                 keys = air_top10_genes$entrezid,
                                 column = "CHR",
                                 keytype = "ENTREZID",
                                 multiVals = "first")


air_top_upregulated <- air_top10_genes[which.max(air_top10_genes$log2FoldChange),"entrezid"]

air_top_upregulated_GOterms <- AnnotationDbi::select(org.Hs.eg.db, 
                                                     keys = air_top_upregulated, 
                                                     column = c("GO", "ONTOLOGY"),
                                                     multiVals = "list")

air_top_upregulated_GOterms <- air_top_upregulated_GOterms %>% filter(ONTOLOGY == "BP")
unique(air_top_upregulated_GOterms$GO)
library(GO.db)

air_GOterms_descriptions <- AnnotationDbi::select(GO.db, keys = unique(air_top_upregulated_GOterms$GO), 
                               columns = "DEFINITION", keytype = "GOID")
head(air_GOterms_descriptions$DEFINITION)
head(air_GOterms_descriptions)


DUSP1_gene <-filter(air_top10_genes, air_top10_genes$symbol =="DUSP1" )

DUSP1_gene_GOterms <- AnnotationDbi::select(org.Hs.eg.db, 
                                                     keys = DUSP1_gene$entrezid, 
                                                     column = c("GO", "ONTOLOGY"),
                                                     multiVals = "list")
DUSP1_gene_GOterms
DUSP1_gene_GOterms <- DUSP1_gene_GOterms %>% filter(ONTOLOGY == "BP")
DUSP1_gene_GOterms_descriptions <- AnnotationDbi::select(GO.db, keys = unique(DUSP1_gene_GOterms$GO), 
                                                  columns = "DEFINITION", keytype = "GOID")


DUSP1_gene_GOterms_descriptions <- DUSP1_gene_GOterms_descriptions[!is.na(DUSP1_gene_GOterms_descriptions$DEFINITION),]
head(DUSP1_gene_GOterms_descriptions$DEFINITION,10)


grep("cortico", DUSP1_gene_GOterms_descriptions$DEFINITION)

DUSP1_gene_GOterms_descriptions$DEFINITION[str_detect(DUSP1_gene_GOterms_descriptions, "cortico")]

library(GOstats)

airway_dge_results
typeof(air_upregulated_genes)
typeof(air_all_genes)


air_upregulated_genes <- airway_dge_results %>% data.frame() %>% 
  filter(log2FoldChange > 1, padj < 0.01) %>% rownames()

air_all_genes <- airway_dge_results %>% data.frame() %>% rownames()

air_test_object <- new("GOHyperGParams",
                   geneIds = air_upregulated_genes,
                   universeGeneIds = air_all_genes, 
                   annotation = "org.Hs.eg.db", 
                   ontology = "BP", 
                   pvalueCutoff = 1,
                   testDirection = "over")


new

goterm_analysis <- hyperGTest(air_test_object)
goterm_analysis


air_goterm_analysis_results <- summary(goterm_analysis)
air_goterm_analysis_results$padj <- p.adjust(air_goterm_analysis_results$Pvalue, method = "BH")

air_goterm_analysis_results <- air_goterm_analysis_results %>% filter(Count > 5) %>% filter(Count < 500)

air_goterm_analysis_top20 <- air_goterm_analysis_results[order(air_goterm_analysis_results$padj)[1:20],]

air_goterm_analysis_top20$Term <- factor(air_goterm_analysis_top20$Term, 
                                     levels = air_goterm_analysis_top20$Term[
                                       order(air_goterm_analysis_top20$padj, decreasing = TRUE)])

air_goterm_analysis_top20 %>% ggplot(aes(x = Term, y = -log10(padj))) +
  geom_point() +
  coord_flip() +
  ylab(expression(-log[10](adjusted~italic(P)~value))) + 
  xlab("GO terms") +
  ggtitle("Top 20 enriched GO terms\n for upregulated genes") +
  theme_bw()




air_downregulated_genes <- airway_dge_results %>% data.frame() %>% 
  filter(log2FoldChange < -1, padj < 0.01) %>% rownames()

air_all_genes <- airway_dge_results %>% data.frame() %>% rownames()

air_down_object <- methods::new("GOHyperGParams",
                       geneIds = air_downregulated_genes,
                       universeGeneIds = air_all_genes, 
                       annotation = "org.Hs.eg.db", 
                       ontology = "BP", 
                       pvalueCutoff = 1,
                       testDirection = "over")

goterm_analysis_down <- hyperGTest(air_down_object)
goterm_analysis_down


airdown_goterm_analysis_results <- summary(goterm_analysis_down)
airdown_goterm_analysis_results$padj <- p.adjust(airdown_goterm_analysis_results$Pvalue, method = "BH")

airdown_goterm_analysis_results <- airdown_goterm_analysis_results %>% filter(Count > 5) %>% filter(Count < 500)

airdown_goterm_analysis_top20 <- airdown_goterm_analysis_results[order(airdown_goterm_analysis_results$padj)[1:20],]

airdown_goterm_analysis_top20$Term <- factor(airdown_goterm_analysis_top20$Term, 
                                         levels = airdown_goterm_analysis_top20$Term[
                                           order(airdown_goterm_analysis_top20$padj, decreasing = TRUE)])

airdown_goterm_analysis_top20 %>% ggplot(aes(x = Term, y = -log10(padj))) +
  geom_point() +
  coord_flip() +
  ylab(expression(-log[10](adjusted~italic(P)~value))) + 
  xlab("GO terms") +
  ggtitle("Top 20 enriched GO terms\n for downregulated genes") +
  theme_bw()






























