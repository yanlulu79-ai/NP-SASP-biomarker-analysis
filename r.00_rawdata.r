
rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "00_rawdata")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

library(GEOquery)  
library(dplyr)     
library(Biobase)   
library(tidyr)     
library(stringr)
library(tibble)
library(data.table)
options(timeout=999999)

gse_data <- getGEO(filename = file.path(output, "GSE124272_series_matrix.txt.gz"),
                   GSEMatrix = TRUE,
                   AnnotGPL = FALSE,
                   getGPL = FALSE)
expr <- exprs(gse_data)

qx <- as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  expr[which(expr <= 0)] <- 0
  expr <- log2(expr + 1)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}

expr <- as.data.frame(expr)
expr <- na.omit(expr)
gse_phenodata <- pData(gse_data)

write.csv(gse_phenodata,file.path(output,"02.test_pd.csv"))

probe_annotation <- fread("00_rawdata/GPL21185-21174.txt", header = TRUE, sep = "\t", skip = "#")
colnames(probe_annotation)

probe_annotation$GeneSymbol <- trimws(sapply(strsplit(probe_annotation$GENE_SYMBOL, "///"), `[`, 1))

expr$GeneSymbol <- probe_annotation$GeneSymbol[match(rownames(expr), probe_annotation$ID)]

table(duplicated(expr$GeneSymbol))

remove_duplicate_genes <- function(expression_data) {
  expression_data <- expression_data[!is.na(expression_data$GeneSymbol), ]
  
  gene_avg_expr <- rowMeans(expression_data[, -which(names(expression_data) == "GeneSymbol")], na.rm = TRUE)
  
  sorted_indices <- order(gene_avg_expr, decreasing = TRUE)
  expression_data <- expression_data[sorted_indices, ]
  
  expression_data <- expression_data[!duplicated(expression_data$GeneSymbol), ]
  
  rownames(expression_data) <- expression_data$GeneSymbol
  
  expression_data <- expression_data[, -which(names(expression_data) == "GeneSymbol")]
  
  return(expression_data)
}

max_expr_data <- remove_duplicate_genes(expr)

gse_phenodata <- gse_phenodata %>%
  mutate(
    group = case_when(
      grepl("patient", title) ~ "NP",
      grepl("Healthy", title) ~ "control",
      TRUE ~ NA_character_  
    )
  ) %>%
  filter(!is.na(group))

table(gse_phenodata$group)

group_info <- data.frame(
  sample = rownames(gse_phenodata),
  group = gse_phenodata$group
)
group_info <- group_info %>%
  arrange(group)

max_expr_data <- max_expr_data[, colnames(max_expr_data) %in% group_info$sample]

max_expr_data <- max_expr_data[, group_info$sample]
write.csv(group_info, file.path(output, "02.test_group.csv"), row.names = FALSE)
write.csv(max_expr_data, file.path(output, "02.test_expr.csv"), row.names = TRUE)

gse_data <- getGEO(filename = file.path(output, "GSE150408_series_matrix.txt.gz"),
                   GSEMatrix = TRUE,
                   AnnotGPL = FALSE,
                   getGPL = FALSE)
expr <- exprs(gse_data)

qx <- as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  expr[which(expr <= 0)] <- 0
  expr <- log2(expr + 1)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}

expr <- as.data.frame(expr)
expr <- na.omit(expr)
gse_phenodata <- pData(gse_data)

write.csv(gse_phenodata,file.path(output,"01.train_pd.csv"))

probe_annotation <- fread("00_rawdata/GPL21185-21174.txt", header = TRUE, sep = "\t", skip = "#")
colnames(probe_annotation)

probe_annotation$GeneSymbol <- trimws(sapply(strsplit(probe_annotation$GENE_SYMBOL, "///"), `[`, 1))

expr$GeneSymbol <- probe_annotation$GeneSymbol[match(rownames(expr), probe_annotation$ID)]

table(duplicated(expr$GeneSymbol))

max_expr_data1 <- remove_duplicate_genes(expr)

gse_phenodata <- gse_phenodata[!grepl("treatment", gse_phenodata$title), ]
gse_phenodata$group <- ifelse(grepl("volunteer", gse_phenodata$title), "control", "NP")
table(gse_phenodata$group)

group_info <- data.frame(
  sample = rownames(gse_phenodata),
  group = gse_phenodata$group
)
group_info <- group_info %>%
  arrange(group)

max_expr_data1 <- max_expr_data1[, colnames(max_expr_data1) %in% group_info$sample]

max_expr_data1 <- max_expr_data1[, group_info$sample]

str(max_expr_data1)
str(max_expr_data)

write.csv(group_info, file.path(output, "01.train_group.csv"), row.names = FALSE)
write.csv(max_expr_data1, file.path(output, "01.train_expr.csv"), row.names = TRUE)
