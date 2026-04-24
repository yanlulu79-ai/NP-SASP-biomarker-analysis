
rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "03_Enrichment")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

library(clusterProfiler)
library(org.Hs.eg.db)  
library(ggplot2)
library(dplyr)
library(enrichplot)
library(stringr)
options(timeout = 300)

genes <- read.csv(file.path(ORIGINAL_DIR, "02_TargetGene/02.intersect_genes.csv"), header = TRUE)
gene_list <- genes$gene  

gene_list_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_results <- enrichGO(gene = gene_list_entrez$ENTREZID, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID", 
                       ont = "ALL",  
                       pAdjustMethod = "none", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 1)  

kegg_results <- enrichKEGG(gene = gene_list_entrez$ENTREZID, 
                           organism = "hsa",
                           pAdjustMethod = "none",
                           pvalueCutoff = 0.05)

go_ontology_counts <- table(go_results@result$ONTOLOGY)

print("GO 通路数量统计:")
print(go_ontology_counts)

kegg_pathway_counts <- nrow(kegg_results)

print("KEGG 显著通路数量统计:")
print(kegg_pathway_counts)

write.csv(go_results, file.path(output, "01.GO_results.csv"))
write.csv(kegg_results, file.path(output, "02.KEGG_results.csv"))

go_bubbleplot <- dotplot(go_results, showCategory = 5, orderBy = "GeneRatio", 
                         label_format = 100, color = "pvalue", split="ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 60))

ggsave(file.path(output, "03.GO_bubbleplot.png"), plot = go_bubbleplot, width = 10, height = 8)
ggsave(file.path(output, "04.GO_bubbleplot.pdf"), plot = go_bubbleplot, width = 10, height = 8)

kegg_bubbleplot <- dotplot(kegg_results, showCategory = 10, orderBy = "GeneRatio", 
                           label_format = 100, color = "pvalue")

ggsave(file.path(output, "05.KEGG_bubbleplot.png"), plot = kegg_bubbleplot, width = 10, height = 8)
ggsave(file.path(output, "06.KEGG_bubbleplot.pdf"), plot = kegg_bubbleplot, width = 10, height = 8)

cat("GO和KEGG富集分析完成，所有结果已保存至", output, "\n")

