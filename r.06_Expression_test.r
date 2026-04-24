
rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "06_Expression_test")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(reshape2)
library(ggpubr)
library(corrplot)
library(data.table)
library(pROC)
library(ggsci)
library(Cairo)
library(glmnet)

test_matrix_raw <- read.csv("00_rawdata/02.test_expr.csv", row.names = 1)

test_group <- fread("00_rawdata/02.test_group.csv")

common_genes <- read.csv("05_Hub_corr/10.common_genes.csv")

common_genes <- common_genes$gene
str(common_genes)

test_matrix <- test_matrix_raw[rownames(test_matrix_raw) %in% common_genes, ]

test_matrix_long <- test_matrix %>%
  
  mutate(gene = rownames(test_matrix)) %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  left_join(test_group, by = "sample")

perform_wilcoxon_test <- function(data) {
  results <- data %>%
    group_by(gene) %>%
    summarise(
      p_value = wilcox.test(expression ~ group)$p.value,
      
      logFC = mean(expression[group == "NP"]) - mean(expression[group == "control"]),
      regulation = ifelse(logFC > 0, "Up", "Down")
    ) %>%
    filter(p_value < 0.05) %>%  
    arrange(p_value)
  return(results)
}

test_results <- perform_wilcoxon_test(test_matrix_long)

deg_results <- read.csv("01_DEG/00.train_expr_DEGs.csv", row.names = 1)

deg_filtered <- deg_results %>%
  filter(abs(logFC) > 0.5, P.Value < 0.05) %>%
  mutate(regulation = ifelse(logFC > 0, "Up", "Down"))

deg_filtered$gene <- rownames(deg_filtered)

merged_results <- inner_join(test_results, deg_filtered, by = "gene", suffix = c("_wilcoxon", "_deg"))

consistent_results <- merged_results %>%
  filter(regulation_wilcoxon == regulation_deg)
write.csv(consistent_results, file.path(output, "00.test_Wilcoxon_results.csv"), row.names = FALSE)

tain_matrix_raw <- read.csv("00_rawdata/01.train_expr.csv", row.names = 1)

tain_group <- fread("00_rawdata/01.train_group.csv")

tain_matrix <- tain_matrix_raw[rownames(tain_matrix_raw) %in% common_genes, ]

tain_matrix_long <- tain_matrix %>%
  
  mutate(gene = rownames(tain_matrix)) %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  left_join(tain_group, by = "sample")

str(tain_matrix_long)
test_matrix_long <- test_matrix_long %>% arrange(gene)
tain_matrix_long <- tain_matrix_long %>% arrange(gene)

p = ggboxplot(test_matrix_long, x = "gene", y = "expression", fill = "group", 
              xlab = "", ylab = "Gene expression", legend.title = "Group",
              palette = c("#0088FF", "#FF5555"), width = 0.75) 
p = p + rotate_x_text(45)

p1 = p + stat_compare_means(aes(group = group),
                            method = "wilcox.test",
                            symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                            label = "p.signif")

pdf(file = "06_Expression_test/03.test_boxplot.pdf", width = 6, height = 5)
print(p1)
dev.off()

png(file = "06_Expression_test/02.test_boxplot.png", width = 6*300, height = 5*300,res=300)
print(p1)
dev.off()

gene_p_values <- deg_filtered[, c("gene", "P.Value")]

gene_p_values <- gene_p_values[rownames(gene_p_values) %in% common_genes, ]

p = ggboxplot(tain_matrix_long, x = "gene", y = "expression", fill = "group", 
              xlab = "", ylab = "Gene expression", legend.title = "Group",
              palette = c("#0088FF", "#FF5555"), width = 0.75) 

p = p + rotate_x_text(45) 

p1 = p + geom_text(data = gene_p_values, aes(x = gene, y = max(tain_matrix_long$expression), 
                                             label = ifelse(P.Value < 0.001, "***", 
                                                            ifelse(P.Value < 0.01, "**", 
                                                                   ifelse(P.Value < 0.05, "*", "")))), 
                   vjust = -0.5, size = 4)

pdf(file = "06_Expression_test/04.tain_boxplot.pdf", width = 6, height = 5)
print(p1)
dev.off()

png(file = "06_Expression_test/05.tain_boxplot.png", width = 6*300, height = 5*300,res=300)
print(p1)
dev.off()

library(glmnet)
library(pROC)
library(ggsci)
library(Cairo)

exp=read.csv("00_rawdata/01.train_expr.csv", row.names = 1)     
gene=read.csv("06_Expression_test/00.test_Wilcoxon_results.csv",row.names = 1)          

rt<-as.data.frame(exp)

geneRT<-row.names(gene)

groupFile="00_rawdata/01.train_group.csv"  
groupInfo=read.csv(groupFile, header=T, stringsAsFactors=F)
y=ifelse(groupInfo$group=="control", 0, 1)

y = y[match(colnames(rt), groupInfo$sample)]
table(y)
geneRT<-as.data.frame(geneRT)

bioCol=pal_aaas("default")(nrow(geneRT))
str(bioCol)

aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
  k=k+1
  
  roc1=roc(y, as.numeric(rt[x,]))     
  if(k==1){
    pdf(file="06_Expression_test/01.tain_ROC.genes.pdf", width=5.5, height=5)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE, lwd=3)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}

legend("bottomright", aucText, lwd=3, bty="n", cex=1, col=bioCol[1:(ncol(rt)-1)])
dev.off()

aucData <- data.frame(Gene = sapply(aucText, function(x) strsplit(x, ", AUC=")[[1]][1]),
                      AUC = sapply(aucText, function(x) as.numeric(sub(".*AUC=", "", x))))

str(aucData)

write.csv(aucData, "06_Expression_test/04.tain_AUC_values.csv", row.names = FALSE)

aucText <- c()
k <- 0

for (x in as.vector(geneRT[, 1])) {
  k <- k + 1
  
  roc1 <- roc(y, as.numeric(rt[x, ]))  
  if (k == 1) {
    
    CairoPNG(file = "06_Expression_test/02.train_ROC.genes.png", width = 5.5 * 300, height = 5 * 300, dpi = 300)
    plot(roc1, print.auc = F, col = bioCol[k], legacy.axes = T, main = "", lwd = 3)
    aucText <- c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  } else {
    plot(roc1, print.auc = F, col = bioCol[k], legacy.axes = T, main = "", add = TRUE, lwd = 3)
    aucText <- c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  }
}

legend("bottomright", aucText, lwd = 3, bty = "n", cex = 1, col = bioCol[1:(ncol(rt) - 1)])

dev.off()

test_exp = read.csv("00_rawdata/02.test_expr.csv", row.names = 1)  
test_groupFile = "00_rawdata/02.test_group.csv"  

test_rt<-as.data.frame(test_exp)

test_groupInfo = read.csv(test_groupFile, header=T, stringsAsFactors=F)
test_y = ifelse(test_groupInfo$group == "control", 0, 1)

test_rt = test_rt[as.vector(geneRT[, 1]), ]
test_rt = t(test_rt)  
test_rt<-data.frame(test_rt)
selected_genes = geneRT[, 1]  
str(selected_genes)

bioCol = pal_aaas("default")(length(selected_genes))  
str(bioCol)

aucText = c()
k = 0

for (x in selected_genes) {
  k = k + 1
  
  roc1 = roc(test_y, as.numeric(test_rt[, x]))  
  
  if (k == 1) {
    pdf(file = "06_Expression_test/03.test_ROC_gene.pdf", width = 5.5, height = 5)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE, lwd=3)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}

legend("bottomright", aucText, lwd=3, bty="n", cex=1, col=bioCol[1:length(selected_genes)])
dev.off()

aucData <- data.frame(Gene = sapply(aucText, function(x) strsplit(x, ", AUC=")[[1]][1]),
                      AUC = sapply(aucText, function(x) as.numeric(sub(".*AUC=", "", x))))

str(aucData)

write.csv(aucData, "06_Expression_test/04.test_AUC_values.csv", row.names = FALSE)

aucText = c()
k = 0

for (x in selected_genes) {
  k = k + 1
  
  roc1 = roc(test_y, as.numeric(test_rt[, x]))  
  
  if (k == 1) {
    png(file = "06_Expression_test/03.test_ROC_gene.png", width = 5.5*300, height = 5*300,res=300)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE, lwd=3)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}

legend("bottomright", aucText, lwd=3, bty="n", cex=1, col=bioCol[1:length(selected_genes)])
dev.off()

test_data <- read.csv("06_Expression_test/04.test_AUC_values.csv")
train_data <- read.csv("06_Expression_test/04.tain_AUC_values.csv")

test_high_auc <- test_data[test_data$AUC > 0.7, ]
train_high_auc <- train_data[train_data$AUC > 0.7, ]

common_high_auc <- merge(test_high_auc, train_high_auc, by = "Gene")

print(common_high_auc)

write.csv(common_high_auc, "06_Expression_test/05.common_high_auc.csv", row.names = FALSE)
