
rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "07_Nomogram")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
set.seed(123)

setwd(ORIGINAL_DIR)

library(data.table)
library(rms)
library(plotly)
library(regplot)
library(ggDCA)
library(pROC)
library(rmda)
library(magick)
library(ResourceSelection)

train_matrix_raw <- read.csv("00_rawdata/01.train_expr.csv", row.names = 1)
train_group <- fread("00_rawdata/01.train_group.csv")

common_genes <- read.csv("06_Expression_test/00.test_Wilcoxon_results.csv")
common_genes <- common_genes$gene

train_matrix <- train_matrix_raw[rownames(train_matrix_raw) %in% common_genes, ]

train_matrix <- as.data.frame(t(train_matrix))
train_matrix$sample <- rownames(train_matrix)  

merged_data <- merge(
  train_group,
  train_matrix,
  by = "sample"  
)

merged_data$group <- ifelse(merged_data$group == "NP", 1, 0)

dd <- datadist(merged_data)
options(datadist = "dd")

fit <- lrm(group ~ CEACAM1 + RRAS2, data = merged_data,x=TRUE,y=TRUE,penalty=0.1)
summary(fit)
str(merged_data)

summary(merged_data$RRAS2) 
summary(merged_data$CEACAM1)

regplot(
  fit,
  plots = c("density", "boxes"),  
  observation = TRUE,             
  title = "NP Risk Nomogram",
  points = TRUE,                  
  droplines = TRUE,               
  clickable = F            
)

img <- image_read("07_Nomogram/01.Nomogram.png")

image_write(img, path = "07_Nomogram/01.Nomogram.pdf", format = "pdf")

predicted_prob <- predict(fit, type = "fitted")

roc_obj <- roc(merged_data$group, predicted_prob)

auc_value <- auc(roc_obj)

png("07_Nomogram/02.roc_curve.png", width = 5.5*300, height = 5*300, res = 300)
plot(roc_obj, main = "ROC Curve for Nomogram", col = "#0072B2", lwd = 3,legacy.axes=T)
legend("bottomright", legend = paste("AUC =", round(auc_value, 2)), col = "#0072B2", lwd = 3)
dev.off()

pdf("07_Nomogram/02.roc_curve.pdf", width = 5.5, height = 5)
plot(roc_obj, main = "ROC Curve for Nomogram", col = "#0072B2", lwd = 3,legacy.axes=T)
legend("bottomright", legend = paste("AUC =", round(auc_value, 2)), col = "#0072B2", lwd = 3)
dev.off()

cal <- calibrate(fit, method = "boot", B = 1000)

y <- fit$y  
pred_prob <- predict(fit, type = "fitted")
hl_test <- hoslem.test(y, pred_prob, g = 10)

png("07_Nomogram/04.calibration_curve.png", width = 6 * 300, height = 5 * 300, res = 300)
plot(cal)  
text(
  x = 0.2, y = 0.9,  
  labels = paste("Hosmer-Lemeshow p =", format(hl_test$p.value, digits = 3)),
  cex = 0.9, col = "black", adj = 0  
)
dev.off()

pdf("07_Nomogram/04.calibration_curve.pdf", width = 6, height = 5)  
plot(cal)  
text(
  x = 0.2, y = 0.9,  
  labels = paste("Hosmer-Lemeshow p =", format(hl_test$p.value, digits = 3)),
  cex = 0.9, col = "black", adj = 0  
)
dev.off()

dca_data <- dca(fit, model.names = "Nomogram")
str(dca_data)
dca_plot<-ggplot(dca_data, aes(x = threshold, y = net_benefit)) +
  geom_line(aes(color = model), linewidth = 1,na.rm = TRUE) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(title = "Decision Curve Analysis", x = "Threshold", y = "Net Benefit")

fit_tkt <- lrm(group ~ CEACAM1, data = merged_data)  
dca_tkt <- dca(fit_tkt, model.names = "CEACAM1")

fit_cxcr6 <- lrm(group ~ RRAS2, data = merged_data)  
dca_cxcr6 <- dca(fit_cxcr6, model.names = "RRAS2")
str(dca_cxcr6)

dca_combined <- rbind(dca_data, dca_tkt, dca_cxcr6)

dca_plot<-ggplot(dca_combined, aes(x = threshold, y = net_benefit)) +
  geom_line(aes(color = model), linewidth = 1,na.rm = TRUE) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7","blue")) +
  labs(title = "Decision Curve Analysis", x = "Threshold", y = "Net Benefit")

ggsave(file.path(output, "03.dca.png"), dca_plot, width = 8, height = 6)
ggsave(file.path(output, "03.dca.pdf"), dca_plot, width = 8, height = 6)
