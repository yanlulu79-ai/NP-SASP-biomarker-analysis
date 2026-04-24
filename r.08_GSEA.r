
rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "08_GSEA")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

set.seed(123)

library(limma)
library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)
library(enrichplot)
library(data.table)
library(ggplot2)
library(GseaVis)
library(GSVA)
library(GSEABase)
library(ComplexHeatmap)

expFile <- "00_rawdata/01.train_expr.csv"  
rt <- read.csv(expFile, header = TRUE, row.names = 1)
mat <- as.data.frame(rt)

gene_symbols_in_mat <- rownames(mat)

geneListFile <- "06_Expression_test/00.test_Wilcoxon_results.csv"
geneList <- read.csv(geneListFile, header = T, sep = ",")
genes_of_interest <- geneList$gene  

batch_cor <- function(gene, mat) {
  y <- as.numeric(mat[gene,])  
  rownames <- rownames(mat)    
  
  
  correlation_data <- do.call(rbind, lapply(rownames, function(x) {
    dd  <- cor.test(as.numeric(mat[x,]), y, type = "spearman", use = "complete.obs")  
    data.frame(gene = gene, mRNAs = x, cor = dd$estimate, p.value = dd$p.value)  
  }))
  
  return(correlation_data)
}

for (gene in genes_of_interest) {
  
  output_file <- paste0("08_GSEA/00.correlation_results_", gene, ".csv")
  
  if (file.exists(output_file)) {
    cat("文件已存在，跳过计算:", output_file, "\n")
    next  
  } else {
    cat("文件不存在，开始计算:", gene, "\n")
    
    system.time(correlation_results <- batch_cor(gene, mat))  
    
    print(head(correlation_results))  
    
    write.csv(correlation_results, file = output_file, row.names = FALSE)
    cat("计算完成，结果已保存:", output_file, "\n")
  }
}

gmtFile <- "08_GSEA/00.c2.cp.v2024.1.Hs.symbols.gmt"
gmt <- read.gmt(gmtFile)

dd<-read.csv("08_GSEA/00.correlation_results_CEACAM1.csv")

gene_df <- data.frame(cor=dd$cor,SYMBOL = dd$mRNAs)

gene_df <- gene_df[order(abs(gene_df$cor), decreasing = TRUE), ]

geneList <- gene_df$cor
names(geneList) = gene_df$SYMBOL
geneList = sort(geneList, decreasing = TRUE)

head(geneList)

KEGG<-GSEA(geneList,TERM2GENE = gmt,pvalueCutoff = 0.05)

kkTab <- as.data.frame(KEGG)

num_pathways <- nrow(kkTab)
cat("总共富集到", num_pathways, "条通路\n")

write.csv(kkTab, file = "08_GSEA/01.CEACAM1_GSEA_KEGG_results.csv", row.names = FALSE)

p<-gseaplot2(
  x = KEGG, 
  geneSetID = kkTab$ID[1:5],  
  title = "",                  
  color = c("#76BA99", "#EB4747", "#996699", "#5C88DA", "#FFCD00", 
            "#A65628", "#F781BF", "#999999", "#E41A1C", "#377EB8"),  
  base_size = 10,              
  rel_heights = c(2, 0.5, 1),  
  subplots = 1:3,             
  pvalue_table = FALSE,       
  ES_geom = "line"          
)

ggsave("08_GSEA/03.CEACAM1_top5.png", plot = p, width = 14, height = 8, bg = "white", dpi = 300)
ggsave("08_GSEA/03.CEACAM1_top5.pdf", plot = p, width = 14, height = 8, bg = "white")

dd<-read.csv("08_GSEA/00.correlation_results_RRAS2.csv")

gene_df <- data.frame(cor=dd$cor,SYMBOL = dd$mRNAs)
gene_df <- gene_df[order(abs(gene_df$cor), decreasing = TRUE), ]

geneList <- gene_df$cor
names(geneList) = gene_df$SYMBOL
geneList = sort(geneList, decreasing = TRUE)

head(geneList)

KEGG<-GSEA(geneList,TERM2GENE = gmt,pvalueCutoff = 0.05)

kkTab <- as.data.frame(KEGG)

num_pathways <- nrow(kkTab)
cat("总共富集到", num_pathways, "条通路\n")

write.csv(kkTab, file = "08_GSEA/02.RRAS2_GSEA_KEGG_results.csv", row.names = FALSE)

p<-gseaplot2(
  x = KEGG, 
  geneSetID = kkTab$ID[1:5],  
  title = "",                  
  color = c("#76BA99", "#EB4747", "#996699", "#5C88DA", "#FFCD00", 
            "#A65628", "#F781BF", "#999999", "#E41A1C", "#377EB8"),  
  base_size = 10,              
  rel_heights = c(2, 0.5, 1),  
  subplots = 1:3,             
  pvalue_table = FALSE,       
  ES_geom = "line"          
)

ggsave("08_GSEA/03.RRAS2_top5.png", plot = p, width = 14, height = 8, bg = "white", dpi = 300)
ggsave("08_GSEA/03.RRAS2_top5.pdf", plot = p, width = 14, height = 8, bg = "white")
