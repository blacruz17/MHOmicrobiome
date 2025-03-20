library(caret)
library(tidyverse)
library(doParallel)
library(xgboost)
library(recipes)

rm(list = ls())

# Loading data ----
library(phyloseq)
physeq <- read_rds("./adj_physeq_MMUPHIN_mpa30_20240911.rds")
physeq <- microbiome::transform(physeq, "clr")
set.seed(505)

## relative abundances table:

data.RF <- data.frame(t(physeq@otu_table))
colnames(data.RF) <- data.frame(physeq@tax_table)$Species

# Labels: 
class.RF <- data.frame(physeq@sam_data) %>% 
  select(MetObesity) %>%
  rename(class = MetObesity)

class.RF <- data.frame("class" = class.RF[rownames(data.RF), ])
rownames(class.RF) <- rownames(data.RF)

## classes for binary predictor:
class.RF$class.binary <- gsub('O|NO', '', class.RF$class)
class.binary <- factor(class.RF$class.binary)

# Binary classifier -----------------------------------------------------------
# start parallel processing
# cl <- makePSOCKcluster(40)
# registerDoParallel(cl)

## Train/Test split ----
# creates split based on the class labels
inTrain <- createDataPartition(class.binary, 
                               p = 3/4, list = FALSE)
# selects the corresponding samples:
trainDescr <- readRDS("./trainDescr_basicBinary_ranger.rds")
testDescr <- readRDS("./testDescr_basicBinary_ranger.rds")
trainClass <- readRDS("./trainClass_basicBinary_ranger.rds")
testClass <- readRDS("./testClass_basicBinary_ranger.rds")

# these should be TRUE:
dim(trainDescr)[1] == length(trainClass)
dim(testDescr)[1] == length(testClass)
# OK!

## Model training --------------------------------------------------------------
# Running an XGBoost model with xgboost requires some additional data preparation.
# xgboost requires a matrix input for the features and the response to be a vector. 
# Consequently, to provide a matrix input of the features we need to encode our 
# categorical variables numerically (i.e. one-hot encoding, label encoding). 
# The following numerically label encodes all categorical features and converts 
# the training data frame to a matrix.

data_train <- trainDescr
data_train$MetHealth <- trainClass


xgb_prep <- recipe(MetHealth ~ ., data = data_train) %>%
  step_integer(all_nominal()) %>%
  prep(training = data_train, retain = TRUE) %>%
  juice()

X <- as.matrix(xgb_prep[setdiff(names(xgb_prep), "MetHealth")])
Y <- xgb_prep$MetHealth
Y[Y == 2] <- 0

# hyperparameter grid
hyper_grid <- readRDS("./hyper_grid_XGB_Reg_20250207.RDS")


# Once you’ve found the optimal hyperparameters, fit the final model with
# xgb.train or xgboost. Be sure to use the optimal number of trees found 
# during cross validation. In our example, adding regularization provides 
# no improvement so we exclude them in our final model.

# NOTA: por un error, en realidad la columna "rmse" se refiere al logloss.
# los mejores valores son:
bestval <- hyper_grid %>%
  filter(logloss > 0) %>%
  arrange(logloss) %>%
  slice_head(n = 1)


params <- list(
  eta = bestval$eta,
  max_depth = bestval$max_depth,
  min_child_weight = bestval$min_child_weight,
  subsample = bestval$subsample,
  colsample_bytree = bestval$colsample_bytree,
  gamma = bestval$gamma,
  lambda = bestval$lambda,
  alpha = bestval$alpha
  )

# train final model
set.seed(123)
xgb.fit.final <- xgboost(
  params = params,
  data = X,
  label = Y,
  nrounds = bestval$trees,
  objective = "binary:logistic",
  verbose = 1
)


saveRDS(xgb.fit.final, "./fit_final_XGB_20250304.RDS")


# stopCluster(cl)


# Test -----
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


predicted.prob <- predict(xgb.fit.final, 
                                   X_t, type="prob") # Prediction
roc <- roc(Y_t, predicted.prob)
roc

(p.roc <- ggroc(roc,
                         size = 2,
                         legacy.axes = TRUE, color = "#725DEF") + 
    labs(x = "False Positive Rate", 
         y = "True Positive Rate") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", linetype = "dotted",
                 linewidth = 2) +
    geom_text(aes( x=0.8, y=0.05,
                   label=paste0("AUC = ",
                                round(roc$auc, 2))),
              size = 6) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    coord_equal()  +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.title = element_text(size = 24),
          legend.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          legend.position = "none"))
