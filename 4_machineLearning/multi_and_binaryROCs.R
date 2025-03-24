library(caret)
library(MLmetrics)
library(pROC)
library(ggpubr)
library(tidyverse)
library(doParallel)
library(patchwork)
library(RColorBrewer)



# Loading data -----------------------------------------------------------------
physeq <- read_rds("../data/physeqMHO.rds")

# selects the same samples:
trainDescr <- readRDS("../results/trainDescr_multi.rds")
trainClass <- readRDS("../results/trainClass_multi.rds")
testDescr <- readRDS("../results/testDescr_multi.rds")
testClass <- readRDS("../results/testClass_multi.rds")

# these should be TRUE:
dim(trainDescr)[1] == length(trainClass)
dim(testDescr)[1] == length(testClass)
# OK!


## Load models -----------------------------------------------------------------
svmRadialMulti <- readRDS("../results/svmRadialMulti.rds")
svmRadialSigmaMulti <- readRDS("../results/svmRadialSigmaMulti.rds")

# Multiclass AUCs --------------------------------------------------------------
## Radial SVM ------------------------------------------------------------------
result.svmR <- predict(svmRadialMulti, testDescr, type="prob") 
roc.svmR <- multiclass.roc(testClass, result.svmR)
# Call:
#   multiclass.roc.default(response = testClass, predictor = result.svmR)
# 
# Data: multivariate predictor result.svmR with 4 levels of testClass: MHNO, MHO, MUNO, MUO.
# Multi-class area under the curve: 0.6164

## Radial + Sigma SVM ----------------------------------------------------------
result.svmRS <- predict(svmRadialSigmaMulti, testDescr, type="prob") 
roc.svmRS <- multiclass.roc(testClass, result.svmRS)
# Call:
#   multiclass.roc.default(response = testClass, predictor = result.svmRS)
# 
# Data: multivariate predictor result.svmRS with 4 levels of testClass: MHNO, MHO, MUNO, MUO.
# Multi-class area under the curve: 0.6164


# ROC Curves -------------------------------------------------------------------
# Make predictions
predictions <- predict(svmRadialMulti, testDescr, type = "prob")

# Initialize list to store ROC curves
roc_curves <- list()

# Calculate ROC curve for each class
for (class in colnames(predictions)) {
  binary_labels <- ifelse(testClass == class, 1, 0)
  roc_curve <- roc(binary_labels, predictions[[class]])
  roc_curves[[class]] <- roc_curve
}

## get AUCs:
for (i in 1:length(roc_curves)){
  roc_curves[[i]]$auc <- auc(roc_curves[[i]])
}

# extracts auc
data.auc.2 <- roc_curves %>% 
  map(~tibble(AUC = .x$auc)) %>%
  unlist() %>% 
  data.frame() %>%
  t()
data.auc.2 <- data.frame("AUC" = data.auc.2[1, ])

data.auc.2$label <- paste0(
  gsub("\\.AUC", " = ", rownames(data.auc.2)),
  round(data.auc.2$AUC, 2))

names(roc_curves) <- data.auc.2$label

# reorder so that legend is in decreasing auc order
data.auc.3 <- data.auc.2[order(data.auc.2[, 1], decreasing = TRUE), ]
roc_curves <- roc_curves[data.auc.3$label]


# Save ROC list for future:
saveRDS(roc_curves, "../results/multiROC_OvA.rds")

ovaROCs <- ggroc(c(roc_curves), legacy.axes = TRUE) + 
  labs(x = "False Positive Rate", 
       y = "True Positive Rate") +
  scale_color_manual(values = viridis::viridis(length(roc_curves))) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color = "black", linetype = "dotted") +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.key = element_rect(linewidth = .1, color = "white"),
        legend.key.size = unit(1, "lines"))  +
  labs(color = "One-vs-Rest AUCs") +
  coord_equal()

ovaROCs

# Shared figure: Binary & Multiclass Models ------------------------------------

# load Binary Model ROC list:
roclist.bi <- readRDS("../results/binaryROC.rds")

names(roclist.bi) <- gsub(", AUC", "", names(roclist.bi))

bi.ROCs <- ggroc(roclist.bi, legacy.axes = TRUE) + 
  labs(x = "False Positive Rate", 
       y = "True Positive Rate") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color = "black", linetype = "dotted") +
  scale_color_manual(values = viridis::viridis(length(roclist.bi))) +
  theme_bw() +
  coord_equal()  +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.key = element_rect(linewidth = .1, color = "white"),
        legend.key.size = unit(1, "lines"))
bi.ROCs

pALL <- bi.ROCs + ovaROCs

pALL2 <- pALL & theme(axis.title = element_text(size = 8),
                    legend.margin=margin(0,0,0,0),
                    legend.box.margin=margin(-5,-5,-5,-5))

ggsave("../figures/machineLearning.png", 
       plot = pALL2,
       width = 4, height = 3, 
       units = "in", dpi = 600)
