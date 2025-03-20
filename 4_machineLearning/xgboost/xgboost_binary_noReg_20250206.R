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
cl <- makePSOCKcluster(40)
registerDoParallel(cl)

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


# Next, we went through a series of grid searches and found the below model 
# hyperparameters (provided via the params argument) to perform quite well. 
# Our RMSE is slightly lower than the best regular and stochastic GBM models.

set.seed(123)
mh_xgb <- xgb.cv(
  data = X,
  label = Y,
  nrounds = 6000,
  objective = "binary:logistic",
  early_stopping_rounds = 50, 
  nfold = 10,
  params = list(
    eta = 0.1,
    max_depth = 3,
    min_child_weight = 3,
    subsample = 0.8,
    colsample_bytree = 1.0),
  verbose = 0
)  

# minimum test CV RMSE
min(mh_xgb$evaluation_log$test_logloss_mean)
## [1] 0.3537141


# Next, we assess if overfitting is limiting our model’s performance by performing 
# a grid search that examines various regularization parameters (gamma, lambda,
# and alpha). Our results indicate that the best performing models use lambda 
# equal to 1 and it doesn’t appear that alpha or gamma have any consistent 
# patterns. However, even when lambda equals 1, our CV RMSE has no improvement 
# over our previous XGBoost model.
# Due to the low learning rate (eta), this cartesian grid search takes a long time.
# We stopped the search after 2 hours and only 98 of the 245 models had completed. 

# hyperparameter grid
hyper_grid <- expand.grid(
  eta = c(0.3, 0.1, 0.05, 0.01, 0.005),
  max_depth = c(1, 3, 5, 6, 8),
  min_child_weight = c(3, 5, 7, 9, 11, 13, 15),
  subsample = c(0.5, 0.75, 1),
  colsample_bytree = c(0.5, 0.75, 1),
  logloss = 0,          # a place to dump logloss results
  trees = 0          # a place to dump required number of trees
)

# grid search
for(i in seq_len(nrow(hyper_grid))) {
  set.seed(123)
  m <- xgb.cv(
    data = X,
    label = Y,
    nrounds = 3000,
    objective = "binary:logistic",
    early_stopping_rounds = 50, 
    nfold = 10,
    verbose = 1,
    params = list( 
      eta = hyper_grid$eta[i], 
      max_depth = hyper_grid$max_depth[i],
      min_child_weight = hyper_grid$min_child_weight[i],
      subsample = hyper_grid$subsample[i],
      colsample_bytree = hyper_grid$colsample_bytree[i],
      nthread = 40
    ) 
  )
  hyper_grid$logloss[i] <- min(m$evaluation_log$test_logloss_mean)
  hyper_grid$trees[i] <- m$best_iteration
}
saveRDS(hyper_grid, "./hyper_grid_XGB_noReg_20250206.RDS")
# results
hyper_grid %>%
  filter(rmse > 0) %>%
  arrange(rmse) %>%
  glimpse()

# Once you’ve found the optimal hyperparameters, fit the final model with
# xgb.train or xgboost. Be sure to use the optimal number of trees found 
# during cross validation. In our example, adding regularization provides 
# no improvement so we exclude them in our final model.


stopCluster(cl)