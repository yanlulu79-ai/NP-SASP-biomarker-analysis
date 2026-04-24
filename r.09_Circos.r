
rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "09_Circos")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

library(AnnoProbe)
library(RCircos)

tain_matrix_raw <- read.csv("00_rawdata/01.train_expr.csv", row.names = 1)
consistent_results<-read.csv("06_Expression_test/00.test_Wilcoxon_results.csv", row.names = NULL)
tain_matrix <- tain_matrix_raw[rownames(tain_matrix_raw) %in% consistent_results$gene, ]

gene_symbols <- rownames(tain_matrix)
gene_pos <- annoGene(gene_symbols, "SYMBOL", "human")
expression_values <- rowMeans(tain_matrix)

biomarker_data <- data.frame(
  Chromosome = gene_pos$chr,
  chromStart = gene_pos$start,
  chromEnd = gene_pos$end,
  Gene = gene_pos$SYMBOL,
  Expression = expression_values
)

biomarker_data <- biomarker_data[complete.cases(biomarker_data), ]

data(UCSC.HG19.Human.CytoBandIdeogram)  
RCircos.Set.Core.Components(
  cyto.info = UCSC.HG19.Human.CytoBandIdeogram,
  chr.exclude = c("chrX", "chrY"),  
  tracks.inside = 4,    
  tracks.outside = 0
)

pdf(file.path(output, "00.Chromosome_Circos_Plot.pdf"), height=10, width=10)
RCircos.Set.Plot.Area()
par(mai = c(0.5, 0.5, 0.5, 0.5))
RCircos.Chromosome.Ideogram.Plot()

RCircos.Gene.Connector.Plot(
  genomic.data = biomarker_data[,1:4],  
  track.num = 1,
  side = "in",  
  genomic.columns = 3,  
  is.sorted = FALSE    
)

par(cex = 2)
RCircos.Gene.Name.Plot(
  gene.data = biomarker_data,
  name.col = 4,
  track.num = 2,
  side = "in",
  genomic.columns = 3,
  is.sorted = FALSE
)

dev.off()

png(file.path(output, "00.Chromosome_Circos_Plot.png"), height=10*300, width=10*300,res=300)
RCircos.Set.Plot.Area()
par(mai = c(0.5, 0.5, 0.5, 0.5))
RCircos.Chromosome.Ideogram.Plot()

RCircos.Gene.Connector.Plot(
  genomic.data = biomarker_data[,1:4],  
  track.num = 1,
  side = "in",  
  genomic.columns = 3,  
  is.sorted = FALSE    
)

par(cex = 2)
RCircos.Gene.Name.Plot(
  gene.data = biomarker_data,
  name.col = 4,
  track.num = 2,
  side = "in",
  genomic.columns = 3,
  is.sorted = FALSE
)

dev.off()
