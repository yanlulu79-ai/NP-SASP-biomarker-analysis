
rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "05_Hub_corr")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

library(glmnet)
library(ggplot2)
library(e1071)
library(caret)
library(pheatmap)
library(ggvenn)
library(randomForest)
library(gbm)
library(gridExtra)
library(png)

set.seed(246)

groupFile <- "00_rawdata/01.train_group.csv"
groupData <- read.csv(groupFile, header = TRUE, sep = ",")
groupData$sample <- as.character(groupData$sample)  

expFile <- "00_rawdata/01.train_expr.csv"
expData <- read.csv(expFile, header = TRUE, row.names = 1)
expData <- expData[, groupData$sample]  

final_hub_genes_data <- read.csv("02_TargetGene/02.intersect_genes.csv", header = TRUE, stringsAsFactors = FALSE)

selectedGenes <- na.omit(final_hub_genes_data$gene)

print(selectedGenes)

selectedGenesLower <- tolower(selectedGenes)
geneNamesLower <- tolower(rownames(expData))
matchedGenes <- selectedGenes[which(selectedGenesLower %in% geneNamesLower)]

expData <- expData[matchedGenes, ]  

expData <- expData[!apply(expData, 1, function(x) all(is.na(x))), ]

expData <- expData[, groupData$sample]

x <- as.matrix(expData)

y <- groupData$group[match(colnames(x), groupData$sample)]  

y <- factor(y, levels = c("control", "NP"))

x <- t(x)

table(y)
str(x)

cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1, type.measure = 'deviance', nfolds = 5)
fit <- cvfit$glmnet.fit
cv_output_pdf <- file.path(output, "01.cvfit.pdf")
cv_output_png <- file.path(output, "02.cvfit.png")
lambda_min <- cvfit$lambda.min
log_lambda_min <- log(lambda_min)

pdf(cv_output_pdf, width = 6, height = 5.5)
par(cex.axis = 1.2, cex.lab = 1.2)  
plot(cvfit)
legend("topright", legend = paste0("log(lambda.min) = ", round(log_lambda_min, 4)), bty = "n")
dev.off()

png(cv_output_png, width = 6 * 300, height = 5.5 * 300, res = 300)
par(cex.axis = 1.2, cex.lab = 1.2)  
plot(cvfit)
legend("topright", legend = paste0("log(lambda.min) = ", round(log_lambda_min, 4)), bty = "n")
dev.off()

coef_output_pdf <- file.path(output, "03.coef_distribution.pdf")
coef_output_png <- file.path(output, "04.coef_distribution.png")

lambda_min <- cvfit$lambda.min
print(lambda_min)

coef_matrix <- as.matrix(coef(fit, s = cvfit$lambda.min))
non_zero_idx <- which(coef_matrix != 0)
lassoGenes <- rownames(coef_matrix)[non_zero_idx]
lassoGenes <- setdiff(lassoGenes, "(Intercept)")  
lassoCoef <- coef_matrix[non_zero_idx][-1]  

pdf(coef_output_pdf, width = 6, height = 5.5)
all_vars <- rownames(coef(fit))[-1]  
selected_idx <- match(lassoGenes, all_vars)  

default_colors <- rainbow(length(all_vars), alpha = 0.8)  
gene_colors <- default_colors[selected_idx]  

plot(fit, xvar = "lambda", col = gene_colors)  

abline(v = log(cvfit$lambda.min), lty = 2, col = "grey50", lwd = 1.8)

legend("topright",
       legend = lassoGenes,
       col = gene_colors,  
       lty = 1,            
       lwd = 1.5,          
       cex = 0.68,
       bty = "n",
       title = paste("Selected Features (n =", length(lassoGenes), ")"))

legend("bottomright", legend = paste0("log(lambda.min) = ", round(log_lambda_min, 4)), bty = "n")

dev.off()

png(coef_output_png, width = 6*300, height = 5.5*300,res=300)
all_vars <- rownames(coef(fit))[-1]  
selected_idx <- match(lassoGenes, all_vars)  
default_colors <- rainbow(length(all_vars), alpha = 0.8)  
gene_colors <- default_colors[selected_idx]  
plot(fit, xvar = "lambda", col = gene_colors)  
abline(v = log(cvfit$lambda.min), lty = 2, col = "grey50", lwd = 1.8)
legend("topright",
       legend = lassoGenes,
       col = gene_colors,  
       lty = 1,            
       lwd = 1.5,          
       cex = 0.68,
       bty = "n",
       title = paste("Selected Features (n =", length(lassoGenes), ")"))
legend("bottomright", legend = paste0("log(lambda.min) = ", round(log_lambda_min, 4)), bty = "n")
dev.off()

lassoGeneCoefFile <- file.path(output, "05.LASSO_genes_with_coef.csv")
lasso_results <- data.frame(Gene = lassoGenes, Coefficient = lassoCoef)
write.csv(lasso_results, lassoGeneCoefFile, row.names = FALSE)

cat("Selected genes and their coefficients: \n")
print(lasso_results)

lasso_genes <- lasso_results$Gene  

control <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
group <- y
num_features <- ncol(x) 
results <- rfe(x = x,  
               y = group, 
               sizes = c(1:num_features),  
               rfeControl = control,
               method = "svmRadial")
print(results)

data <- results$results
error_rate <- 1 - data$Accuracy

PlotErrors <- function(errors, errors2=NULL, no.info=0.5, 
                       ylim=range(c(errors, errors2), na.rm=T), 
                       xlab='Number of Features',  ylab='5 x CV Error') {
  AddLine <- function(x, col='#99CC00FF') {
    lines(which(!is.na(errors)), na.omit(x), col=col, lwd=3)
    points(which.min(x), min(x, na.rm=T), col='firebrick3')
    text(which.min(x), min(x, na.rm=T), 
         paste0("n=", which.min(x), " (", format(min(x, na.rm=T), dig=3), ")"), 
         pos=2, col='red', cex=1)
  }
  
  plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
  AddLine(errors)
  if(!is.null(errors2)) AddLine(errors2, 'gray30')
  abline(h=no.info, lty=2)
}

PlotAccuracy <- function(errors, errors2=NULL, no.info=0.5, 
                         ylim=range(c(errors, errors2), na.rm=T), 
                         xlab='Number of Features',  ylab='5 x CV Accuracy') {
  AddLine <- function(x, col='#99CC00FF') {
    lines(which(!is.na(errors)), na.omit(x), col=col, lwd=3)
    points(which.max(x), max(x, na.rm=T), col='firebrick3')
    text(which.max(x), max(x, na.rm=T), 
         paste0("n=", which.max(x), " (", format(max(x, na.rm=T), dig=3), ")"), 
         pos=2, col='red', cex=1)
  }
  
  plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
  AddLine(errors)
  if(!is.null(errors2)) AddLine(errors2, 'gray30')
  abline(h=no.info, lty=2)
}

data <- results$results
accuracy <- data$Accuracy
error_rate <- 1 - accuracy

pdf("05_Hub_corr/06.SVM_RFE_Error_Rate.pdf", width=7, height=6)
PlotErrors(error_rate, no.info=0.5, ylab="Cross-Validated Error Rate")
dev.off()

png("05_Hub_corr/06.SVM_RFE_Error_Rate.png", width=7, height=6, units="in", res=300)
PlotErrors(error_rate, no.info=0.5, ylab="Cross-Validated Error Rate")
dev.off()

pdf("05_Hub_corr/07.SVM_RFE_Accuracy.pdf", width=7, height=6)
PlotAccuracy(accuracy, no.info=0.5, ylab="Cross-Validated Accuracy")
dev.off()

png("05_Hub_corr/07.SVM_RFE_Accuracy.png", width=7, height=6, units="in", res=300)
PlotAccuracy(accuracy, no.info=0.5, ylab="Cross-Validated Accuracy")
dev.off()

img1 <- readPNG("05_Hub_corr/06.SVM_RFE_Error_Rate.png")
img2 <- readPNG("05_Hub_corr/07.SVM_RFE_Accuracy.png")

combined_plot <- grid.arrange(rasterGrob(img1), rasterGrob(img2), ncol = 2)

ggsave("05_Hub_corr/08.combined_plot.pdf", combined_plot, width = 12, height = 6)
ggsave("05_Hub_corr/08.combined_plot.png", combined_plot, width = 12, height = 6)

svmrfe_result <- predictors(results)
print(svmrfe_result)
length(svmrfe_result)

write.csv(svmrfe_result,file.path(output,"07.svm_gene.csv"), row.names = FALSE)

lasso_genes <- lasso_results$Gene  

rfe_genes <- svmrfe_result

geneList <- list(
  Lasso = lasso_genes,
  SVM_RFE = rfe_genes
)

venn_plot <- ggvenn(
  geneList, 
  show_percentage = T,  
  stroke_color = "white",  
  stroke_size = 0.5,       
  fill_color = c("#D8B7DD", "lightblue"),  
  set_name_color = c("black", "black"),  
  set_name_size = 6,       
  text_size = 4.5,
  digits = 1
)

venn_pdf <- file.path(output, "09.venn_diagram.pdf")
venn_png <- file.path(output, "09.venn_diagram.png")

pdf(file = venn_pdf, width = 6, height = 6)
grid::grid.draw(venn_plot)
dev.off()

png(file = venn_png, width = 6*300, height = 6*300, res = 300)
grid::grid.draw(venn_plot)
dev.off()

common_genes <- intersect(lasso_genes, rfe_genes)
cat("关键枢纽基因：", common_genes, "\n")

common_genes_df <- data.frame(gene = common_genes)

common_genes_csv <- file.path(output, "10.common_genes.csv")

write.csv(common_genes_df, file = common_genes_csv, row.names = FALSE)
