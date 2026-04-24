
rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "01_DEG")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggrepel)

group_file <- "00_rawdata/01.train_group.csv"
group_info <- read.csv(group_file)
expr_data <-read.csv('00_rawdata/01.train_expr.csv', row.names = 1)

group_info <- group_info[match(colnames(expr_data), group_info$sample), ]

if (nrow(group_info) != ncol(expr_data)) {
  cat(paste("样本数不匹配，跳过数据集\n"))
  next
}

group <- factor(group_info$group)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)  

print(paste("Design matrix columns for:", colnames(design)))  

if (nrow(design) != ncol(expr_data)) {
  stop(paste("设计矩阵的行数与表达数据的列数不匹配，跳过数据集：", dataset_name))
}

fit <- lmFit(expr_data, design)

contrast_matrix <- makeContrasts(NP_vs_control = NP - control, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

DEGs <- topTable(fit2, coef = "NP_vs_control", number = Inf, adjust.method = "BH")

DEGs_filtered <- DEGs[abs(DEGs$logFC) > 0.5 & DEGs$P.Value < 0.05, ]
DEGs_filtered$regulation <- ifelse(DEGs_filtered$logFC > 0, "Up", "Down")

file_name <- sprintf("00.train_expr_DEGs.csv")

write.csv(DEGs_filtered, file.path(output, file_name))

cat(paste("差异分析完成，", file_name, "已保存。\n"))

DEGs_df <- data.frame(DEGs)

DEGs_df$Significance <- ifelse(
  DEGs_df$P.Value < 0.05 & DEGs_df$logFC > 0.5,
  "Up",
  ifelse(
    DEGs_df$P.Value < 0.05 & DEGs_df$logFC < -0.5,
    "Down",
    "Not"
  )
)

print(table(DEGs_df$Significance)) 

DEGs_df$Label <- ifelse(
  DEGs_df$P.Value < 0.05 & abs(DEGs_df$logFC) > 0.5,
  rownames(DEGs_df),
  NA
)

upregulated_genes <- DEGs_df[DEGs_df$logFC > 0.5, ]
upregulated_genes <- upregulated_genes[upregulated_genes$P.Value < 0.05,]
downregulated_genes <- DEGs_df[DEGs_df$logFC < -0.5, ]
downregulated_genes <- downregulated_genes[downregulated_genes$P.Value < 0.05,]

top5_upregulated <- upregulated_genes[order(upregulated_genes$P.Value), ][1:5, ]
top5_downregulated <- downregulated_genes[order(downregulated_genes$P.Value), ][1:5, ]

top10_genes <- rbind(top5_upregulated, top5_downregulated)
str(top10_genes)

DEGs_df$Label[rownames(DEGs_df) %in% rownames(top10_genes)] <- rownames(top10_genes)

p <- ggplot(DEGs_df, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = Significance)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  scale_color_manual(values = c("blue", "grey85", "darkred")) +
  ggtitle(paste("Volcano Plot:NP_DEGs")) +
  theme(plot.title = element_text(size = 14, hjust = 1, face = "bold"))

p <- p + geom_text_repel(
  data = top10_genes,
  aes(label = Label),
  box.padding = 0.3,
  point.padding = 0.4,
  min.segment.length = 0.01,
  size = 3,
  force = 4,
  na.rm = TRUE
) +
  theme_bw()

pdf_file <- paste0("01.NP_DEGs1_volcano_plot.pdf")
pdf(file.path(output, pdf_file), width = 6, height = 6)
print(p)
dev.off()

png_file <- paste0("01.NP_DEGs1_volcano_plot.png")
png(file.path(output, png_file), width = 6*300, height = 6*300, res = 300)
print(p)
dev.off()

cat(paste("火山图已保存为：", pdf_file, "和", png_file, "\n"))

top10_up <- DEGs_df %>% filter(Significance == "Up") %>% top_n(10, -P.Value)

top10_down <- DEGs_df %>% filter(Significance == "Down") %>% top_n(10, -P.Value)
top20 <- c(rownames(top10_up), rownames(top10_down))

de_expr <- expr_data[top20, ]
madt <- as.matrix(de_expr)
madt2 <- t(scale(t(madt)))
range(madt2)

col_fun <- colorRamp2(c(min(madt2), 0, max(madt2)), c("#4F94CD", "white", "#8B1A1A"))

ha1 = HeatmapAnnotation(
  group = group_info$group,  
  col = list(group = c("control" = "#4682B4", "NP" = "#E9967A"))
)

pdf_file <- paste0("02.NP_DEGs1_heatmap_plot.pdf")
pdf(file.path(output, pdf_file), width = 10, height = 8)
p2 <- densityHeatmap(madt2, 
                     col = colorRampPalette(c("#4F94CD", "white", "#8B1A1A"))(50), 
                     quantile_gp = gpar(fontsize = 9), 
                     title = paste("Heatmap Plot:NP_DEGs"), 
                     ylab = "Expression", 
                     top_annotation = ha1) %v% 
  Heatmap(madt2, 
          name = "expression", 
          col = col_fun, 
          show_row_names = TRUE, 
          show_column_names = FALSE, 
          cluster_rows = FALSE, 
          height = unit(15, "cm"))
print(p2)
dev.off()

png_file <- paste0("02.NP_DEGs1_heatmap_plot.png")
png(file.path(output, png_file), width = 10*300, height = 8*300, res = 300)
p2 <- densityHeatmap(madt2, 
                     col = colorRampPalette(c("#4F94CD", "white", "#8B1A1A"))(50), 
                     quantile_gp = gpar(fontsize = 9), 
                     title = paste("Heatmap Plot:NP_DEGs"), 
                     ylab = "Expression", 
                     top_annotation = ha1) %v% 
  Heatmap(madt2, 
          name = "expression", 
          col = col_fun, 
          show_row_names = TRUE, 
          show_column_names = FALSE, 
          cluster_rows = FALSE, 
          height = unit(15, "cm"))
print(p2)
dev.off()

cat(paste("热图已保存为：", pdf_file, "和", png_file, "\n"))
