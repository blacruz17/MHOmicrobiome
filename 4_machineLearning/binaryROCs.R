library(caret)
library(MLmetrics)
library(pROC)
library(ggpubr)
library(tidyverse)
library(doParallel)
library(patchwork)
library(recipes)
# also required: ranger, randomForest, kernlab, LiblineaR

rm(list = ls())

# Loading data ----
physeq <- read_rds("../data/physeqMHO.rds")

# selects the same samples:
trainDescr <- readRDS("../results/trainDescr_binary.rds")
trainClass <- readRDS("../results/trainClass_binary.rds")
testDescr <- readRDS("../results/testDescr_binary.rds")
testClass <- readRDS("../results/testClass_binary.rds")

# these should be TRUE:
dim(trainDescr)[1] == length(trainClass)
dim(testDescr)[1] == length(testClass)
# OK!

## Load models ------

b.ranger <- readRDS("../results/basicBinary_ranger.rds") 
b.r.weights <- readRDS("../results/basicBinary_ranger_weighted.rds")
b.ranger.u <- readRDS("../results/upsampledBinary_ranger.rds")
b.ranger.d <- readRDS("../results/downsampledBinary_ranger.rds")
svmRadial <- readRDS("../results/svmRadial.rds")
svmRadialSigma <- readRDS("../results/svmRadialSigma.rds")
svmRadialWeights <- readRDS("../results/svmRadialWeights.rds")
svmRadialWeights.d <- readRDS("../results/svmRadialWeights_downsampl.rds")
svmRadialWeights.u <- readRDS("../results/svmRadialWeights_upsampl.rds")
svmLinearWeights <- readRDS("../results/svmLinearWeights.rds")
svmLinearWeights2 <- readRDS("../results/svmLinearWeights2.rds")
xgboost <- readRDS("../results/finalFitXGB.RDS")

# take a look at ranger models to get best version:
res.b <- lapply(basic, function(x) x$results)
res.b

# res.b: 1000 trees

lapply(b.ranger, function(x){
  res <- x$results
  fit <- x$bestTune
  return(res %>% 
           inner_join(fit, by = c("mtry", "splitrule", "min.node.size")))})
# b.ranger: 100 trees

lapply(b.r.weights, function(x){
  res <- x$results
  fit <- x$bestTune
  return(res %>% 
           inner_join(fit, by = c("mtry", "splitrule", "min.node.size")))})
# 100 trees 

lapply(b.ranger.d, function(x){
  res <- x$results
  fit <- x$bestTune
  return(res %>% 
           inner_join(fit, by = c("mtry", "splitrule", "min.node.size")))})
# 500 trees

lapply(b.ranger.u, function(x){
  res <- x$results
  fit <- x$bestTune
  return(res %>% 
           inner_join(fit, by = c("mtry", "splitrule", "min.node.size")))})
# 100 trees

best_models <- c(basic["1000"], b.ranger["100"], b.r.weights["100"],
                 b.ranger.d["500"], b.ranger.u["100"])
names(best_models) <- c("Original", "Ranger", "Ranger + Weights",
                        "Ranger + Downsampling", "Ranger + Upsampling")



# ROC Curve: Original model -------------
result.predicted.prob <- predict(best_models[["Original"]], 
                                 testDescr, type="prob") # Prediction
result.roc <- roc(testClass, result.predicted.prob$MU)
result.roc

(p.roc <- ggroc(result.roc,
                size = 2,
                legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(result.roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))

# ROC: Ranger -----
ranger.predicted.prob <- predict(best_models[["Ranger"]], 
                                 testDescr, type="prob") # Prediction
ranger.roc <- roc(testClass, ranger.predicted.prob$MU)
ranger.roc

(p.roc.ranger <- ggroc(ranger.roc,
                size = 2,
                legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(ranger.roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))

# ROC: Weighted -----
weighted.predicted.prob <- predict(best_models[["Ranger + Weights"]], 
                                 testDescr, type="prob") # Prediction
weighted.roc <- roc(testClass, weighted.predicted.prob$MU)
weighted.roc

(p.roc.weighted <- ggroc(weighted.roc,
                       size = 2,
                       legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(weighted.roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))

# Ranger: Up ---------------
ranger.u.predicted.prob <- predict(best_models[["Ranger + Upsampling"]], 
                                 testDescr, type="prob") # Prediction
ranger.u.roc <- roc(testClass, ranger.u.predicted.prob$MU)
ranger.u.roc

(p.roc.ranger.u <- ggroc(ranger.u.roc,
                       size = 2,
                       legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(ranger.roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))

# Ranger: Down ---------------
ranger.d.predicted.prob <- predict(best_models[["Ranger + Downsampling"]], 
                                   testDescr, type="prob") # Prediction
ranger.d.roc <- roc(testClass, ranger.d.predicted.prob$MU)
ranger.d.roc

(p.roc.ranger.d <- ggroc(ranger.d.roc,
                         size = 2,
                         legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(ranger.roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))

# Join together ------
allROCs <- ggroc(list(
            # "Weighted" = weighted.roc,
           "Ranger" = ranger.roc,
           # "Basic" = result.roc,
           "Ranger + Upsampling" = ranger.u.roc,
           "Ranger + Downsampling" = ranger.d.roc),
      size = 1,
      legacy.axes = TRUE) + 
  labs(x = "False Positive Rate", 
       y = "True Positive Rate") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color = "black", linetype = "dotted") +
  scale_color_manual(values = viridis::viridis(3)) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  coord_equal()  +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10))
allROCs

# Support Vector Machines ----------------------------
## Radial ---------------
svmRadial.prob <- predict(svmRadial, testDescr, type="prob") # Prediction
svmRadial.roc <- roc(testClass, svmRadial.prob$MU)
svmRadial.roc

(p.roc.svmRadial <- ggroc(svmRadial.roc,
                         size = 2,
                         legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(svmRadial.roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))

## Radial Sigma ---------------
svmRadialSigma.prob <- predict(svmRadialSigma, testDescr, type="prob") # Prediction
svmRadialSigma.roc <- roc(testClass, svmRadialSigma.prob$MU)
svmRadialSigma.roc

(p.roc.svmRadialSigma <- ggroc(svmRadialSigma.roc,
                          size = 2,
                          legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(svmRadialSigma.roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))


## Radial Weighted ---------------
svmRadialWeights.prob <- predict(svmRadialWeights, testDescr, type="prob") # Prediction
svmRadialWeights.roc <- roc(testClass, svmRadialWeights.prob$MU)
svmRadialWeights.roc

(p.roc.svmRadialWeights <- ggroc(svmRadialWeights.roc,
                               size = 2,
                               legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(svmRadialWeights.roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))

## Radial Weighted + Upsampling ---------------
svmRadialWeights.u.prob <- predict(svmRadialWeights.u, testDescr, type="prob") # Prediction
svmRadialWeights.u.roc <- roc(testClass, svmRadialWeights.u.prob$MU)
svmRadialWeights.u.roc

(p.roc.svmRadialWeights.u <- ggroc(svmRadialWeights.u.roc,
                                 size = 2,
                                 legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(svmRadialWeights.roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))

## Radial Weighted + Downsampling ---------------
svmRadialWeights.d.prob <- predict(svmRadialWeights.d, testDescr, type="prob") # Prediction
svmRadialWeights.d.roc <- roc(testClass, svmRadialWeights.d.prob$MU)
svmRadialWeights.d.roc

(p.roc.svmRadialWeights.d <- ggroc(svmRadialWeights.d.roc,
                                   size = 2,
                                   legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(svmRadialWeights.roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))

## Linear Weighted ---------------
svmLinearWeights.prob <- predict(svmLinearWeights, testDescr, type="prob") # Prediction
svmLinearWeights.roc <- roc(testClass, svmLinearWeights.prob$MU)
svmLinearWeights.roc

(p.roc.svmLinearWeights <- ggroc(svmLinearWeights.roc,
                                 size = 2,
                                 legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(svmLinearWeights.roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))

# Together --------------------------
allROCs <- ggroc(list(
  # "Weighted" = weighted.roc,
  "Radial" = svmRadial.roc,
  "Radial Sigma" = svmRadialSigma.roc,
  "Radial Weighted" = svmRadialWeights.roc,
  "Radial Weighted + Upsampling" = svmRadialWeights.u.roc,
  "Radial Weighted + Downsampling" = svmRadialWeights.d.roc,
  "Linear Weighted" = svmLinearWeights.roc),
  size = 1,
  legacy.axes = TRUE) + 
  labs(x = "False Positive Rate", 
       y = "True Positive Rate") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color = "black", linetype = "dotted") +
  scale_color_manual(values = viridis::viridis(6)) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  coord_equal()  +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10))
allROCs

# XGBoost -------------------------------
data_test <- testDescr
data_test$MetHealth <- testClass

xgb_prep <- recipe(MetHealth ~ ., data = data_test) %>%
  step_integer(all_nominal()) %>%
  prep(training = data_test, retain = TRUE) %>%
  juice()

X_t <- as.matrix(xgb_prep[setdiff(names(xgb_prep), "MetHealth")])
Y_t <- xgb_prep$MetHealth
Y_t[Y_t == 1] <- 0 # 0 = MH
Y_t[Y_t == 2] <- 1 # 1 = MU


predicted.prob.xgb <- predict(xgboost,
                              X_t, type="prob") # Prediction
roc.xgb <- roc(Y_t, predicted.prob.xgb)
roc.xgb

(p.roc.xgb <- ggroc(roc.xgb,
                size = 2,
                legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(roc.xgb$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))

# Figura final  -----------------------------
roclist <- list(
  "Random Forest" = ranger.roc,
  "Support Vector Machine" = svmRadialWeights.roc,
  "Gradient Boosting" = roc.xgb)

names(roclist) <- c(
  paste0(names(roclist)[1], ", AUC = ", round(ranger.roc$auc, 2)),
  paste0(names(roclist)[2], ", AUC = ", round(svmRadialWeights.roc$auc, 2)),
  paste0(names(roclist)[3], ", AUC = ", round(roc.xgb$auc, 2))
)


roclist.2 <- list(roclist[[2]], roclist[[1]], roclist[[3]])
names(roclist.2) <- c(names(roclist)[2], names(roclist)[1], names(roclist)[3])


# Save ROC list for future:
saveRDS(roclist.2, "../results/binaryROC.rds")
