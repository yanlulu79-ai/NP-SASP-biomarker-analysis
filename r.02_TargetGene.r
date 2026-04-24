
rm(list = ls()); gc()

ORIGINAL_DIR <- ""                           
output <- file.path(ORIGINAL_DIR, "02_TargetGene")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

library(ggvenn)

train_expr_DEGs <- read.csv("01_DEG/00.train_expr_DEGs.csv", row.names = 1)
train_genes <- rownames(train_expr_DEGs)

WH_genes <- read.csv("00_rawdata/SRGs.csv")
WH_genes <- WH_genes$gene

intersect_genes <- intersect(train_genes, WH_genes)

intersect_file <- file.path(output, "02.intersect_genes.csv")
write.csv(data.frame(gene = intersect_genes), intersect_file, row.names = FALSE)

cat("交集基因已保存至：", intersect_file, "\n")
cat("交集基因数量：", length(intersect_genes), "\n")
print(intersect_genes)

venn_data <- list(
  train_expr_DEGs = train_genes,
  SRGs_genes = WH_genes
)

venn_plot <- ggvenn(
  venn_data, 
  show_percentage = TRUE,  
  stroke_color = "white",  
  stroke_size = 0.5,       
  fill_color = c("#D8B7DD", "lightblue"),  
  set_name_color = c("black", "black"),  
  set_name_size = 6,       
  text_size = 4.5 ,  
  digits = 1   
)

venn_pdf <- file.path(output, "02.venn_diagram.pdf")
venn_png <- file.path(output, "03.venn_diagram.png")

pdf(file = venn_pdf, width = 6, height = 6)
grid::grid.draw(venn_plot)
dev.off()

png(file = venn_png, width = 6 * 300, height = 6 * 300, res = 300)
grid::grid.draw(venn_plot)
dev.off()

cat("韦恩图已保存至：", venn_png, "\n")
