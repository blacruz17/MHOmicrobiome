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

iniDir <- "~/ai4food"
iniDir <- "C:/Users/blancalacruz/OneDrive - FUNDACION IMDEA-ALIMENTACION/Escritorio/AI4Food"
setwd(paste0(iniDir, "/mho"))


# Loading data ----
physeq <- read_rds("./adj_physeq_MMUPHIN_mpa30_20240911.rds")

setwd(paste0(iniDir, "/mho/randomForest/mpav30_f1_mmuphin/cccFiles_jan25"))
# selects the corresponding samples:
trainDescr <- readRDS("trainDescr_basicBinary_ranger.rds")
trainDescr1 <- readRDS("trainDescr_basicBinary.rds")
trainClass <- readRDS("trainClass_basicBinary_ranger.rds")
testDescr <- readRDS("testDescr_basicBinary_ranger.rds")
testClass <- readRDS("testClass_basicBinary_ranger.rds")

# these should be TRUE:
dim(trainDescr)[1] == length(trainClass)
dim(testDescr)[1] == length(testClass)
# OK!

# chequeo:
all(trainDescr == trainDescr1)
# Genial, podemos comparar en el mismo script todos los modelos binarios.

## Load models ------

basic <- readRDS("basicBinary.rds") # F
b.ranger <- readRDS("basicBinary_ranger.rds") # Kappa
b.r.weights <- readRDS("basicBinary_ranger_weighted.rds")
b.ranger.u <- readRDS("upsampledBinary_ranger.rds")
b.ranger.d <- readRDS("downsampledBinary_ranger.rds")
svmRadial <- readRDS("svmRadial.rds")
svmRadialSigma <- readRDS("svmRadialSigma.rds")
svmRadialWeights <- readRDS("svmRadialWeights.rds")
svmRadialWeights.d <- readRDS("svmRadialWeights_downsampl.rds")
svmRadialWeights.u <- readRDS("svmRadialWeights_upsampl.rds")
svmLinearWeights <- readRDS("svmLinearWeights.rds")
svmLinearWeights2 <- readRDS("svmLinearWeights2.rds")
xgboost <- readRDS("fit_final_XGB_20250207.RDS")

# cotilleamos los modelos para quedarnos con lo mejor:
res.b <- lapply(basic, function(x) x$results)
res.b

# res.b: nos quedamos el de 1000 arboles.

lapply(b.ranger, function(x){
  res <- x$results
  fit <- x$bestTune
  return(res %>% 
           inner_join(fit, by = c("mtry", "splitrule", "min.node.size")))})
# b.ranger: nos quedamos con el de 100 arboles

lapply(b.r.weights, function(x){
  res <- x$results
  fit <- x$bestTune
  return(res %>% 
           inner_join(fit, by = c("mtry", "splitrule", "min.node.size")))})
# 100 arboles 

lapply(b.ranger.d, function(x){
  res <- x$results
  fit <- x$bestTune
  return(res %>% 
           inner_join(fit, by = c("mtry", "splitrule", "min.node.size")))})
# 500

lapply(b.ranger.u, function(x){
  res <- x$results
  fit <- x$bestTune
  return(res %>% 
           inner_join(fit, by = c("mtry", "splitrule", "min.node.size")))})
# 100

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

# write_rds(result.roc, "./ROCtesting.RDS")
# 
# ggsave("./originalmodelROC.png", plot = p.roc,
#        width = 15,
#        height = 14,
#        units = "cm",
#        dpi = 1200)


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

# ggsave("./upsampledmodelROC.png", plot = p.roc,
#        width = 15,
#        height = 14,
#        units = "cm",
#        dpi = 1200)

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

ggsave("./binaryModelsRF_ROC.png", 
       plot = allROCs,
       width = 10,
       height = 8,
       units = "cm",
       dpi = 1200)

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

ggsave("./binaryModels_SVM_ROC.png", 
       plot = allROCs,
       width = 10,
       height = 8,
       units = "cm",
       dpi = 1200)


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
saveRDS(roclist.2, "./binaryROC_20250212.RDS")

allROCs <- ggroc(roclist.2,
  # size = 1.5,
  legacy.axes = TRUE) + 
  labs(x = "False Positive Rate", 
       y = "True Positive Rate") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color = "black", linetype = "dotted",
               # linewidth = 1
               ) +
  scale_color_manual(values = viridis::viridis(length(roclist.2))) +
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
allROCs

ggsave("./binaryModels_RF_SVM_XGB_ROC.png", 
       plot = allROCs,
       width = 3, height = 3, units = "in", dpi = 1200)

# Variable importance plot -----
features_to_plot <- 20

svm.imp <- varImp(svmRadialWeights)$importance %>% 
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  arrange(MU) %>%
  slice_tail(n = features_to_plot)

# ranger.imp <- varImp(best_models[["Ranger"]])$importance %>% 
#   as.data.frame() %>%
#   tibble::rownames_to_column() %>%
#   arrange(MU) %>%
#   slice_tail(n = features_to_plot)
# 
# xgboost.imp <- varImp(xgboost)$importance %>% 
#   as.data.frame() %>%
#   tibble::rownames_to_column() %>%
#   arrange(MU) %>%
#   slice_tail(n = features_to_plot)
# 


## SVM -----------------------------------------------------------------------
svm.imp

p2 <- svm.imp %>%
  mutate(rowname = gsub("_", " ", rowname)) %>%
  mutate(rowname = fct_inorder(rowname)) %>%
  rename(feature = rowname) %>%
  ggplot(aes(x = feature, 
             y = MU))+
  geom_segment( aes(x=feature, xend=feature,
                    y=0, yend=MU), 
                color="#5B8EFD",
                linewidth = 1) +
  geom_point( color="#5B8EFD", size=3) +
  # scale_color_manual(values = '#21908CFF') +
  coord_flip() +
  labs(y = 'Feature Importance', x = 'Features',
       fill = 'Model') +
  theme_bw()  + 
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.text.y = element_text(face = "italic"),
        legend.position = 'none') +
  theme(axis.text.y = ggtext::element_markdown())

p2

df2 <- psmelt(physeq)

# boxplot:
p.b <- df2 %>%
  rename(id_voluntario = Sample,
         group = met_health,
         variable = Species,
         value = Abundance) %>%
  select(id_voluntario, group, variable, value) %>% 
  inner_join(svm.imp %>% rename(variable = rowname), by = "variable") %>% 
  arrange(MU) %>% 
  mutate(variable = gsub("_", " ", variable),
         variable = fct_inorder(variable),
         value = as.numeric(value)) %>% 
  ggplot(aes(x = variable, y = value,
             fill = group,
             color = group,
             shape = group)) +
  geom_boxplot(alpha=.8, outlier.shape = NA) +
  geom_point(size = .5, alpha = .3,
             position = position_jitterdodge()) +
  theme_bw() +
  coord_flip() +
  labs(y = 'Feature Value', x = 'Features',
       fill = 'Model') +
  scale_fill_manual(values = c( "#4C1D4BFF", "#F2704DFF")) +
  scale_color_manual(values = c( "#4C1D4BFF", "#F2704DFF")) +
  theme_bw()  + 
  theme(#plot.title = element_text(face = 'bold', size = 24),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.text.y = element_text(face = "italic"),
        axis.title.y = element_blank()) +
  theme(axis.text.y = ggtext::element_markdown())
p.b

(superplot <- p.b + 
    (p2 + theme(axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank())))

# ggsave("./SVM_VarImp.png", 
#        plot = superplot,
#        width = 28,
#        height = 18,
#        units = "cm",
#        dpi = 1200)
