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


# Rows: 1,575
# Columns: 7
# $ eta              <dbl> 0.050, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.100, 0.005, 0.0…
# $ max_depth        <dbl> 1, 1, 5, 6, 8, 3, 3, 6, 1, 3, 5, 6, 8, 1, 1, 3, 3, 5, 6, 8, 6, 3, …
# $ min_child_weight <dbl> 15, 15, 15, 15, 15, 3, 15, 9, 15, 11, 15, 15, 15, 15, 15, 11, 15, …
# $ subsample        <dbl> 0.75, 0.75, 0.75, 0.75, 0.75, 1.00, 0.75, 1.00, 0.75, 0.50, 0.75, …
# $ colsample_bytree <dbl> 0.50, 1.00, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.50, 0.75, 0.50, …
# $ logloss          <dbl> 0.3669600, 0.3681760, 0.3682180, 0.3682180, 0.3682180, 0.3683751, …
# $ trees            <dbl> 287, 1465, 707, 707, 707, 426, 708, 46, 2675, 850, 845, 845, 845, …

# Ahora podemos usar estos hiperparametros para probar la regularizacion!

# hyperparameter grid
hyper_grid <- expand.grid(
  eta = 0.050,
  max_depth = 1, 
  min_child_weight = 15,
  subsample = 0.75, 
  colsample_bytree = 0.50,
  gamma = c(0, 1, 10, 100, 1000),
  lambda = c(0, 1e-2, 0.1, 1, 100, 1000, 10000),
  alpha = c(0, 1e-2, 0.1, 1, 100, 1000, 10000),
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
      gamma = hyper_grid$gamma[i], 
      lambda = hyper_grid$lambda[i], 
      alpha = hyper_grid$alpha[i],
      nthread = 40
    ) 
  )
  hyper_grid$logloss[i] <- min(m$evaluation_log$test_logloss_mean)
  hyper_grid$trees[i] <- m$best_iteration
}
saveRDS(hyper_grid, "./hyper_grid_XGB_Reg_20250207.RDS")

# results
hyper_grid %>%
  arrange(logloss) %>%
  glimpse()


stopCluster(cl)
