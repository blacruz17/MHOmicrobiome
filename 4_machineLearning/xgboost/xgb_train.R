library(caret)
library(tidyverse)
library(doParallel)
library(xgboost)
library(recipes)

rm(list = ls())

# Loading data ----
library(phyloseq)
physeq <- read_rds("../data/physeqMHO.rds")
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
# selects the same samples:
trainDescr <- readRDS("../results/trainDescr_binary.rds")
trainClass <- readRDS("../results/trainClass_binary.rds")
testDescr <- readRDS("../results/testDescr_binary.rds")
testClass <- readRDS("../results/testClass_binary.rds")

# these should be TRUE:
dim(trainDescr)[1] == length(trainClass)
dim(testDescr)[1] == length(testClass)
# OK!

## Model training --------------------------------------------------------------
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
hyper_grid <- readRDS("./hyper_grid_XGB_Regularized.RDS")

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


saveRDS(xgb.fit.final, "./finalFitXGB.RDS")


# stopCluster(cl)
