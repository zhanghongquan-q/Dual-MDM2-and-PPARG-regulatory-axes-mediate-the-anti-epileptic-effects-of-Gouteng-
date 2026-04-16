library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)
library(fs)
library(DALEXtra)
library(reshape2)
set.seed(2025)
expression_matrix <- read.csv("exp1.csv", row.names = 1)
core_genes <- read.csv("hubgene.csv")[,2] 
sample_classification <- read.csv("group1.csv", row.names = 1)
data <- expression_matrix[core_genes, , drop = FALSE]
data <- na.omit(data)
data <- t(data)
data <- as.data.frame(data)
data$Type <- sample_classification[rownames(data), "Title"]
colnames(data) <- make.names(colnames(data), unique = TRUE)  

inTrain <- createDataPartition(y = data$Type, p = 0.7, list = F)
train <- data[inTrain, ]
test <- data[-inTrain, ]

control <- trainControl(method = "repeatedcv", number = 5, savePredictions = TRUE)
mod_rf <- train(Type ~ ., data = train, method = 'rf',
                trControl = control)
mod_svm <- train(Type ~ ., data = train, method = 'svmRadial',
                 prob.model = TRUE, trControl = control)

mod_knn <- train(Type ~ ., data = train, method = 'knn',
                 trControl = control)

mod_nnet <- train(Type ~ ., data = train, method = 'nnet',
                  trControl = control)

mod_lasso <- train(Type ~ ., data = train, method = 'glmnet',
                   trControl = control)

mod_dt <- train(Type ~ ., data = train, method = 'rpart',
                trControl = control)

p_fun <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[,2]
}

yTest <- ifelse(test$Type == "Health", 0, 1)  
test_for_prediction <- test[, !colnames(test) %in% "Type"]


explainer_rf <- explain(mod_rf, label = "RF",
                        data = test_for_prediction, y = yTest,
                        predict_function = p_fun,
                        verbose = TRUE)
mp_rf <- model_performance(explainer_rf)
explainer_svm <- explain(mod_svm, label = "SVM",
                         data = test_for_prediction, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_svm <- model_performance(explainer_svm)

explainer_knn <- explain(mod_knn, label = "KNN",
                         data = test_for_prediction, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_knn <- model_performance(explainer_knn)

explainer_nnet <- explain(mod_nnet, label = "NNET",
                          data = test_for_prediction, y = yTest,
                          predict_function = p_fun,
                          verbose = FALSE)
mp_nnet <- model_performance(explainer_nnet)

explainer_lasso <- explain(mod_lasso, label = "LASSO",
                           data = test_for_prediction, y = yTest,
                           predict_function = p_fun,
                           verbose = FALSE)
mp_lasso <- model_performance(explainer_lasso)

explainer_dt <- explain(mod_dt, label = "DT",
                        data = test_for_prediction, y = yTest,
                        predict_function = p_fun,
                        verbose = FALSE)
mp_dt <- model_performance(explainer_dt)



pdf(file = "boxplot.pdf",width = 6,height = 6)
p2 <- plot(mp_rf,mp_svm,
           #mp_glm,
           mp_knn,
           mp_nnet,mp_lasso,
           mp_dt,
           geom = "boxplot")
print(p2)
dev.off()

pred1=predict(mod_rf,newdata=test,type="prob")
pred2=predict(mod_svm,newdata=test,type="prob")
pred3=predict(mod_knn,newdata=test,type="prob")
pred4=predict(mod_nnet,newdata=test,type="prob")
pred5=predict(mod_lasso,newdata=test,type="prob")
pred6=predict(mod_dt,newdata=test,type="prob")

roc1=roc(yTest,as.numeric(pred1[,2]))
roc2=roc(yTest,as.numeric(pred2[,2]))
roc3=roc(yTest,as.numeric(pred3[,2]))
roc4=roc(yTest,as.numeric(pred4[,2]))
roc5=roc(yTest,as.numeric(pred5[,2]))
roc6=roc(yTest,as.numeric(pred6[,2]))

pdf(file = "ROC.pdf",width = 5,height = 5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="chocolate", lwd=2)  # RF
plot(roc2, print.auc=F, legacy.axes=T, main="", col="aquamarine3", lwd=2, add=T)  # SVM
plot(roc3, print.auc=F, legacy.axes=T, main="", col="darkgoldenrod3", lwd=2, add=T)  # KNN
plot(roc4, print.auc=F, legacy.axes=T, main="", col="darkolivegreen3", lwd=2, add=T)  # NNET
plot(roc5, print.auc=F, legacy.axes=T, main="", col="dodgerblue4", lwd=2, add=T)  # LASSO
plot(roc6, print.auc=F, legacy.axes=T, main="", col="dodgerblue", lwd=2, add=T)  # DT


legend_labels <- c(
  paste0("RF (AUC = ", sprintf("%.03f", roc1$auc), ")"),
  paste0("SVM (AUC = ", sprintf("%.03f", roc2$auc), ")"),
  paste0("KNN (AUC = ", sprintf("%.03f", roc3$auc), ")"),
  paste0("NNET (AUC = ", sprintf("%.03f", roc4$auc), ")"),
  paste0("LASSO (AUC = ", sprintf("%.03f", roc5$auc), ")"),
  paste0("DT (AUC = ", sprintf("%.03f", roc6$auc), ")")
)
legend_cols <- c("chocolate", "aquamarine3", "darkgoldenrod3", "darkolivegreen3", "dodgerblue4", "dodgerblue")

legend(
  x = "bottomright",    
  legend = legend_labels,  
  col = legend_cols,    
  lwd = 2,             
  bty = "n"            
)
dev.off()



folds <- createFolds(data$Type, k = 5)
get_cv_prob <- function(model, data, folds) {
  prob <- numeric(nrow(data))
  for(i in 1:5){
    test_idx <- folds[[i]]
    prob[test_idx] <- predict(model, newdata = data[test_idx,], type = "prob")[,2]
  }
  return(prob)
}

prob_rf    <- get_cv_prob(mod_rf, data, folds)
prob_svm   <- get_cv_prob(mod_svm, data, folds)
prob_knn   <- get_cv_prob(mod_knn, data, folds)
prob_nnet  <- get_cv_prob(mod_nnet, data, folds)
prob_lasso <- get_cv_prob(mod_lasso, data, folds)
prob_dt    <- get_cv_prob(mod_dt, data, folds)


y_real <- as.factor(data$Type)
pdf("5_ROC.pdf", width=7, height=6)
roc_cv_rf    <- roc(y_real, prob_rf)
roc_cv_svm   <- roc(y_real, prob_svm)
roc_cv_knn   <- roc(y_real, prob_knn)
roc_cv_nnet  <- roc(y_real, prob_nnet)
roc_cv_lasso <- roc(y_real, prob_lasso)
roc_cv_dt    <- roc(y_real, prob_dt)

plot(roc_cv_rf,    col = "#E74C3C",      lwd=2, legacy.axes=T, main="5-Fold Cross-Validation ROC")
plot(roc_cv_svm,   col = "#3498DB",      lwd=2, add=T)
plot(roc_cv_knn,   col = "#F1C40F",      lwd=2, add=T)
plot(roc_cv_nnet,  col = "#2ECC71",      lwd=2, add=T)
plot(roc_cv_lasso, col = "#9B59B6",      lwd=2, add=T)
plot(roc_cv_dt,    col = "#E67E22",      lwd=2, add=T)

legend("bottomright",
       legend = c(paste0("RF (AUC=",round(roc_cv_rf$auc,3),")"),
                  paste0("SVM (AUC=",round(roc_cv_svm$auc,3),")"),
                  paste0("KNN (AUC=",round(roc_cv_knn$auc,3),")"),
                  paste0("NNET (AUC=",round(roc_cv_nnet$auc,3),")"),
                  paste0("LASSO (AUC=",round(roc_cv_lasso$auc,3),")"),
                  paste0("DT (AUC=",round(roc_cv_dt$auc,3),")")),
       col = c("#E74C3C","#3498DB","#F1C40F","#2ECC71","#9B59B6","#E67E22"),
       lwd=2, cex=0.8, bty="n")
dev.off()


auc_values <- c(round(roc_cv_rf$auc,3),round(roc_cv_svm$auc,3),round(roc_cv_knn$auc,3),
                round(roc_cv_nnet$auc,3),round(roc_cv_lasso$auc,3),round(roc_cv_dt$auc,3))
model_names <- c("RF","SVM","KNN","NNET","LASSO","DT")
auc_df <- data.frame(Model = model_names, AUC = auc_values)
pdf("5CV_AUC.pdf", width=8, height=5)
ggplot(auc_df, aes(x = reorder(Model, AUC), y = AUC, fill = Model)) +
  geom_col(alpha = 0.8, width = 0.6) + geom_text(aes(label = AUC), vjust = -0.3, size = 4) +
  scale_fill_brewer(palette = "Set2") + theme_bw() +
  labs(x = "Model", y = "5-Fold CV AUC", title = "5-Fold CV AUC Comparison") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
dev.off()

pdf("5CV_mix.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))
models_list <- list(RF = mod_rf, SVM = mod_svm, KNN = mod_knn, NNET = mod_nnet, LASSO = mod_lasso, DT = mod_dt)
for (name in names(models_list)) {
  cv_pred <- models_list[[name]]$pred$pred
  cv_obs <- models_list[[name]]$pred$obs
  cv_pred <- factor(cv_pred, levels = levels(cv_obs))
  cm <- confusionMatrix(cv_pred, cv_obs)
  fourfoldplot(cm$table, main = paste(name, "5-CV"), 
               col = c("#E8F4FD", "#2E86AB"))
}
dev.off()

############################ Clinical###############################################
library(ggpubr)   
library(corrplot) 
library(patchwork)

all_data <- data   
target_genes <- c("MDM2", "PPARG")
gene_data <- all_data[, c(target_genes, "Type")]
gene_data <- na.omit(gene_data)
gene_data$Group <- ifelse(gene_data$Type == "Health", "Health", "EP")
gene_data_long <- reshape2::melt(gene_data, 
                                 id.vars = c("Type", "Group"),
                                 variable.name = "Gene", 
                                 value.name = "Expression")
p_box <- ggboxplot(gene_data_long, 
                   x = "Group", 
                   y = "Expression",
                   fill = "Group",          
                   palette = c("#2E86AB", "#A23B72"),
                   facet.by = "Gene",       
                   ncol = 2,
                   xlab = "Group",
                   ylab = "Expression",
                   legend.title = "Type") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif",    
                     size = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12))
print(p_box)
dev.off()

gene_data$Phenotype <- ifelse(gene_data$Type == "Health", 0, 1)

p_scatter <- ggscatter(cor_data, 
                       x = "Phenotype", 
                       y = c("MDM2", "PPARG"),
                       combine = TRUE,
                       color = "#E63946",
                       add = "reg.line",  
                       conf.int = TRUE,   
                       xlab = "Clinical Phenotype (0=Health, 1=EP)",
                       ylab = "Gene Expression Level") +
  stat_cor(method = "spearman", label.x.npc = 0.2, size = 4) + 
  theme_bw()
print(p_scatter)
dev.off()



library(broom)       
library(ResourceSelection) 
library(ggplot2)
library(pROC)
library(ggpubr)
library(car)

lr_data <- all_data[, c("MDM2", "PPARG", "Type")]
lr_data <- na.omit(lr_data)
lr_data$y <- ifelse(lr_data$Type == "EP", 1, 0) 


lr_mdm2 <- glm(y ~ MDM2, data = lr_data, family = binomial())
hl_mdm2 <- hoslem.test(lr_data$y, fitted(lr_mdm2), g=10)
res_mdm2 <- tidy(lr_mdm2, conf.int = TRUE) 
res_mdm2$estimate <- exp(res_mdm2$estimate) 
res_mdm2$conf.low <- exp(res_mdm2$conf.low)
res_mdm2$conf.high <- exp(res_mdm2$conf.high)
res_mdm2$model <- "MDM2_Single"
res_mdm2$HosmerLemeshow_p <- hl_mdm2$p.value

lr_pparg <- glm(y ~ PPARG, data = lr_data, family = binomial())
hl_pparg <- hoslem.test(lr_data$y, fitted(lr_pparg), g=10)
res_pparg <- tidy(lr_pparg, conf.int = TRUE)
res_pparg$estimate <- exp(res_pparg$estimate)
res_pparg$conf.low <- exp(res_pparg$conf.low)
res_pparg$conf.high <- exp(res_pparg$conf.high)
res_pparg$model <- "PPARG_Single"
res_pparg$HosmerLemeshow_p <- hl_pparg$p.value


lr_combined <- glm(y ~ MDM2 + PPARG, data = lr_data, family = binomial())
hl_combined <- hoslem.test(lr_data$y, fitted(lr_combined), g=10)
res_combined <- tidy(lr_combined, conf.int = TRUE)
res_combined$estimate <- exp(res_combined$estimate)
res_combined$conf.low <- exp(res_combined$conf.low)
res_combined$conf.high <- exp(res_combined$conf.high)
res_combined$model <- "MDM2_PPARG_Combined"
res_combined$HosmerLemeshow_p <- hl_combined$p.value


lr_all_results <- rbind(res_mdm2, res_pparg, res_combined)
lr_all_results$estimate[lr_all_results$estimate < 1e-10] <- 1e-10
lr_all_results$conf.low[lr_all_results$conf.low < 1e-10] <- 1e-10
lr_all_results$conf.high[lr_all_results$conf.high < 1e-10] <- 1e-10
lr_all_results$Effect_Direction <- ifelse(
  lr_all_results$estimate < 1 & lr_all_results$p.value < 0.05,
  "Protective (Significant)",
  ifelse(
    lr_all_results$estimate < 1,
    "Protective (Not Significant)",
    "Risk"
  )
)


lr_data$pred_mdm2 <- predict(lr_mdm2, type = "response")
lr_data$pred_pparg <- predict(lr_pparg, type = "response")
lr_data$pred_comb <- predict(lr_combined, type = "response")

auc_mdm2 <- auc(roc(lr_data$y, lr_data$pred_mdm2, quiet=TRUE))
auc_pparg <- auc(roc(lr_data$y, lr_data$pred_pparg, quiet=TRUE))
auc_comb <- auc(roc(lr_data$y, lr_data$pred_comb, quiet=TRUE))

lr_perf <- data.frame(
  Model = c("MDM2 (Single)", "PPARG (Single)", "MDM2+PPARG (Combined)"),
  AUC = round(c(auc_mdm2, auc_pparg, auc_comb), 3),
  HosmerLemeshow_P = round(c(hl_mdm2$p.value, hl_pparg$p.value, hl_combined$p.value), 3),
  Effect = c("Significant Protective", "Significant Protective", "MDM2: Sig; PPARG: NS")
)


pdf("logfore.pdf", width=10, height=6)
plot_data <- lr_all_results[lr_all_results$term != "(Intercept)", ]
plot_data <- plot_data[plot_data$model %in% c("MDM2_Single", "PPARG_Single"), ]
plot_data$term <- factor(plot_data$term, levels = c("MDM2", "PPARG"))
ggplot(plot_data, 
       aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange(size=1.5, fatten=4, color="#E74C3C", linewidth=1.2) +
  geom_hline(yintercept = 1, linetype="dashed", color="red", linewidth=1.5) +
  scale_y_log10(limits = c(1e-10, 10), breaks = c(1e-10, 1e-5, 0.001, 0.1, 1, 10)) +
  coord_flip() +
  labs(x = "Protective Gene", 
       y = "Odds Ratio (OR, log10 scale, <1 = Protective)", 
       title = "MDM2 and PPARG are Independent Protective Genes for Epilepsy") +
  theme_bw(base_size=14) +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        panel.grid.minor = element_blank())
dev.off()

prob_mdm2 <- data.frame(
  Group = ifelse(lr_data$Type == "EP", "EP", "Health"),
  Model = "MDM2",
  Probability = lr_data$pred_mdm2
)
prob_pparg <- data.frame(
  Group = ifelse(lr_data$Type == "EP", "EP", "Health"),
  Model = "PPARG",
  Probability = lr_data$pred_pparg
)
prob_comb <- data.frame(
  Group = ifelse(lr_data$Type == "EP", "EP", "Health"),
  Model = "Combined",
  Probability = lr_data$pred_comb
)
lr_plot <- rbind(prob_mdm2, prob_pparg, prob_comb)

pdf("Logbox.pdf", width=10, height=6)
ggboxplot(lr_plot, 
          x="Group", 
          y="Probability", 
          fill="Group", 
          palette = c("Health"="#2ECC71", "EP"="#E74C3C"),
          facet.by="Model", 
          ncol=3) +
  stat_compare_means(method="wilcox.test", label="p.signif", size=5) +
  labs(x="Group", y="Predicted Probability of EP",
       title="Predicted EP Risk by Protective Genes") +
  theme_bw(base_size=12) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
dev.off()



pdf("MDM2.pdf", width=6, height=5)
ggplot(lr_data, aes(x=MDM2, fill=factor(ifelse(Type=="EP","EP","Health")), color=factor(ifelse(Type=="EP","EP","Health")))) +
  geom_density(alpha=0.7, size=1.2) +
  scale_fill_manual(values=c("Health"="#2ECC71","EP"="#E63946"), name="Group") +
  scale_color_manual(values=c("Health"="#2ECC71","EP"="#E63946"), name="Group") +
  labs(x="MDM2 Expression", y="Density", title="MDM2 (Protective Gene) Expression Density") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))
dev.off()

pdf("PPARG.pdf", width=6, height=5)
ggplot(lr_data, aes(x=PPARG, fill=factor(ifelse(Type=="EP","EP","Health")), color=factor(ifelse(Type=="EP","EP","Health")))) +
  geom_density(alpha=0.7, size=1.2) +
  scale_fill_manual(values=c("Health"="#3498DB","EP"="#E63946"), name="Group") +
  scale_color_manual(values=c("Health"="#3498DB","EP"="#E63946"), name="Group") +
  labs(x="PPARG Expression", y="Density", title="PPARG (Protective Gene) Expression Density") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))
dev.off()

