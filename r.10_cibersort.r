
rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "10_cibersort")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

library(IOBR)
library(limma)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(estimate)
library(ggh4x)

expFile <- "00_rawdata/01.train_expr.csv"
rt <- read.csv(expFile, header = TRUE, row.names = 1)
mat=as.matrix(rt)

groupFile <- "00_rawdata/01.train_group.csv"
groupData <- read.csv(groupFile, header = TRUE, sep = ",")
groupData$sample <- as.character(groupData$sample)  

rds_file <- "10_cibersort/afcibersort_results.rds"

if (file.exists(rds_file)) {
  
  afcibersort <- readRDS(rds_file)
  cat("Loaded afcibersort from RDS file.\n")
} else {
  
  afcibersort <- deconvo_tme(mat, method = "cibersort", perm = 1000)
  
  saveRDS(afcibersort, rds_file)
  cat("Saved afcibersort to RDS file.\n")
}

write.csv(afcibersort, "10_cibersort/01.CIBERSORT.csv", row.names = FALSE)  

mapping <- c(
  "B_cells_naive_CIBERSORT" = "Naive_B_cells_CIBERSORT",
  "B_cells_memory_CIBERSORT" = "Memory_B_cells_CIBERSORT",
  "Plasma_cells_CIBERSORT" = "Plasma_cells_CIBERSORT",
  "T_cells_CD8_CIBERSORT" = "CD8_T_cells_CIBERSORT",
  "T_cells_CD4_naive_CIBERSORT" = "Naive_CD4_T_cells_CIBERSORT",
  "T_cells_CD4_memory_resting_CIBERSORT" = "Resting_memory_CD4_T_cells_CIBERSORT",
  "T_cells_CD4_memory_activated_CIBERSORT" = "Activated_memory_CD4_T_cells_CIBERSORT",
  "T_cells_follicular_helper_CIBERSORT" = "T_follicular_helper_cells_CIBERSORT",
  "T_cells_regulatory_(Tregs)_CIBERSORT" = "Regulatory_T_cells_(Tregs)_CIBERSORT",
  "T_cells_gamma_delta_CIBERSORT" = "Gamma_delta_T_cells_CIBERSORT",
  "NK_cells_resting_CIBERSORT" = "Resting_NK_cells_CIBERSORT",
  "NK_cells_activated_CIBERSORT" = "Activated_NK_cells_CIBERSORT",
  "Monocytes_CIBERSORT" = "Monocytes_CIBERSORT",
  "Macrophages_M0_CIBERSORT" = "M0_macrophages_CIBERSORT",
  "Macrophages_M1_CIBERSORT" = "M1_macrophages_CIBERSORT",
  "Macrophages_M2_CIBERSORT" = "M2_macrophages_CIBERSORT",
  "Dendritic_cells_resting_CIBERSORT" = "Resting_dendritic_cells_CIBERSORT",
  "Dendritic_cells_activated_CIBERSORT" = "Activated_dendritic_cells_CIBERSORT",
  "Mast_cells_resting_CIBERSORT" = "Resting_mast_cells_CIBERSORT",
  "Mast_cells_activated_CIBERSORT" = "Activated_mast_cells_CIBERSORT",
  "Eosinophils_CIBERSORT" = "Eosinophils_CIBERSORT",
  "Neutrophils_CIBERSORT" = "Neutrophils_CIBERSORT"
)

new_colnames <- colnames(afcibersort)
new_colnames <- ifelse(new_colnames %in% names(mapping), mapping[new_colnames], new_colnames)
colnames(afcibersort) <- new_colnames

cibersort=as.data.frame(afcibersort)
rownames(cibersort)=cibersort[,1]
immune=cibersort[,-1]
immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data)=gsub("_CIBERSORT"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))

row.names(groupData) = groupData$sample

cluster = groupData[, "group", drop = FALSE]
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))
data$Type=factor(data$Type, levels=c(cluster[1,1],cluster[nrow(cluster),1]))
data$GSM=rownames(data)
data=melt(data,id.vars=c("Type","GSM"))
colnames(data)=c("Type","GSM","Celltype", "Freq")
Cellratio=data
colourCount = length(unique(Cellratio$Celltype))

colaa=colaa=col=c(pal_aaas()(10),pal_jco()(10),pal_nejm()(8),pal_d3()(10),pal_jama()(6))
ggplot(Cellratio) + 
  geom_bar(aes(x = GSM, y = Freq, fill = Celltype), 
           stat = "identity", width = 0.7, size = 0.5, colour = '#222222') + 
  facet_grid2(
    . ~ Type, 
    scales = "free",
    strip = strip_themed(
      background_x = elem_list_rect(fill = c("#0088FF", "#CDBE70")),
      text_x = elem_list_text(colour = "white")
    )
  ) +
  theme_classic() +
  labs(x = 'Cell cycle phase', y = 'Ratio') +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    legend.position = "right",
    strip.text.x = element_text(size = 20, colour = "white")
  ) +
  scale_fill_manual(values = colaa) +
  theme_bw() +
  xlab(NULL) +
  theme(axis.text.x = element_blank()) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))

ggsave("10_cibersort/02.immune.ration.pdf",width = 15,height = 8)        

ggsave("10_cibersort/03.immune.ration.png", width = 15, height = 8)

cibersort=as.data.frame(afcibersort)
rownames(cibersort)=cibersort[,1]
immune=cibersort[,-1]
immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data)=gsub("_CIBERSORT"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))
data=melt(data,id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
group=levels(factor(data$Type))
data$Type=factor(data$Type, c(cluster[1,1],cluster[nrow(cluster),1]))
bioCol=pal_jco()(6)
bioCol=bioCol[1:length(group)]

significant_cells <- data.frame(Immune = character(), p_value = numeric())

immune_cells <- unique(data$Immune)

for (immune_cell in immune_cells) {
  
  subset_data <- data %>% filter(Immune == immune_cell)
  
  test_result <- wilcox.test(Expression ~ Type, data = subset_data)
  
  p_value <- test_result$p.value
  
  significant_cells <- rbind(significant_cells, data.frame(Immune = immune_cell, p_value = p_value))
}

significant_cells_filtered <- significant_cells %>% filter(p_value < 0.05)

print(significant_cells_filtered)

write.csv(significant_cells_filtered, "10_cibersort/04.significant_immune_cells.csv")

boxplot=ggboxplot(data, x="Immune", y="Expression",  fill="Type",
                  xlab="",
                  ylab="CIBERSORT Fraction",
                  legend.title="Type", 
                  width=0.8,
                  palette=bioCol,add.params = list(size=0.1))
boxplot=boxplot+
  stat_compare_means(aes(group=Type), method = "wilcox",symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme_bw()+
  rotate_x_text(50)

pdf(file="10_cibersort/05.immune.diff.pdf", width=9, height=4.5)
print(boxplot)
dev.off()

png("10_cibersort/06.immune.diff.png", width = 9, height = 4.5, units = "in", res = 300)
print(boxplot)
dev.off()

exp = rt  
group = groupData$group  

group = as.factor(group)

geneRT <- read.csv("06_Expression_test/00.test_Wilcoxon_results.csv", row.names = 1)
exp = exp[rownames(exp) %in% rownames(geneRT), ]

data = t(exp) 
data<-as.data.frame(data)

afcibersort1 = afcibersort[, -((ncol(afcibersort)-2):ncol(afcibersort))]

colnames(afcibersort1)=gsub("_CIBERSORT"," ",colnames(afcibersort1))  
afcibersort1<-as.data.frame(afcibersort1)
rownames(afcibersort1) <- afcibersort1$ID

afcibersort1 <- afcibersort1[, -which(names(afcibersort1) == "ID")]

significant_cells <- trimws(significant_cells_filtered$Immune)
significant_cells <- gsub(" ", "_", significant_cells)

colnames(afcibersort1) <- trimws(colnames(afcibersort1))
significant_cells <- trimws(significant_cells)
matching_columns <- which(colnames(afcibersort1) %in% significant_cells)

afcibersort1 <- afcibersort1[, matching_columns]

outTab = data.frame()
significant_cells = list()  

for(cell in colnames(afcibersort1)){
  if(sd(afcibersort1[, cell]) == 0) { next }  
  cell_has_significant = FALSE  
  
  for(gene in colnames(data)){
    x = as.numeric(afcibersort1[, cell])
    y = as.numeric(data[, gene])
    corT = cor.test(x, y, method="spearman")
    cor = corT$estimate
    pvalue = corT$p.value
    
    if (pvalue < 1) {
      cell_has_significant = TRUE
    }
    text = ifelse(pvalue < 0.001, "***", 
                  ifelse(pvalue < 0.01, "**", 
                         ifelse(pvalue < 0.05, "*", "")))
    outTab = rbind(outTab, cbind(Gene = gene, Immune = cell, cor, text, pvalue))
  }
  
  if (cell_has_significant) {
    significant_cells[[cell]] = TRUE
  }
}

write.csv(outTab,file.path(output,"08.gene_cell_cor.csv"))

outTab$cor=as.numeric(outTab$cor)
pdf(file="10_cibersort/07.cor_gene_cell.pdf", width=7.5, height=4.5)
ggplot(outTab, aes(Immune, Gene)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +    
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),   
        axis.text.y = element_text(size = 8, face = "bold")) +       
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   
  scale_x_discrete(position = "bottom")      
dev.off()

png(file="10_cibersort/08.cor_gene_cell.png", width=7.5*300, height=4.5*300,res=300)
ggplot(outTab, aes(Immune, Gene)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +    
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),   
        axis.text.y = element_text(size = 8, face = "bold")) +       
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   
  scale_x_discrete(position = "bottom")      
dev.off()

