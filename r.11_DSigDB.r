
rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "11_DSigDB")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
library(data.table)

setwd(ORIGINAL_DIR)

library(data.table)
library(ggplot2)
library(stringr)

data_file <- "11_DSigDB/00.DSigDB_All_detailed.txt"
data <- fread(data_file, header = TRUE, sep = "\t", select = c("Drug", "Gene", "Type"), stringsAsFactors = FALSE, fill = TRUE)
data <- unique(data, by = c("Drug", "Gene"))
str(data)

target_genes <- c("CEACAM1", "RRAS2")
filtered_data <- data[Gene %in% target_genes, ]

filtered_file <- file.path(output, "00.filtered_drugs.csv")
write.csv(filtered_data, filtered_file, row.names = FALSE, quote = FALSE)

network_file <- file.path(output, "01.gene_drug_network.tsv")
unique_network <- unique(filtered_data[, .(Source = Gene, Target = Drug)])
write.table(unique_network, network_file, sep = "\t", row.names = FALSE, quote = FALSE)
network_file <- file.path(output, "01.gene_drug_network.csv")
write.csv(unique_network, network_file, row.names = FALSE)

drugs_CEACAM1 <- unique(filtered_data[Gene == "CEACAM1", Drug])
drugs_RRAS2 <- unique(filtered_data[Gene == "RRAS2", Drug])
intersect_drugs <- intersect(drugs_CEACAM1, drugs_RRAS2)

attribute_file <- file.path(output, "02.gene_drug_attributes.tsv")

gene_attributes <- data.table(
  ID = unique(filtered_data$Gene),
  Type = "Gene"
)

drug_attributes <- data.table(
  ID = unique(filtered_data$Drug),
  Type = ifelse(unique(filtered_data$Drug) %in% intersect_drugs, "Intersect_Drug", "Drug")
)
all_attributes <- rbind(gene_attributes, drug_attributes)
write.table(all_attributes, attribute_file, sep = "\t", row.names = FALSE, quote = FALSE)

attribute_file <- file.path(output, "02.gene_drug_attributes.csv")
write.csv(all_attributes, attribute_file)

gene_drug_count <- unique_network[, .N, by = Source]
cat("每个基因对应的药物数量：\n")
print(gene_drug_count)

intersect_file <- file.path(output, "03.intersect_drugs.csv")
intersect_data <- data.table(intersect_drugs = intersect_drugs)
str(intersect_data)
print(intersect_data)
write.csv(intersect_data, intersect_file, row.names = FALSE, quote = FALSE)

cat("网络文件、属性文件和交集药物文件已生成并保存到", output, "\n")

