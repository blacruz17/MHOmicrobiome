library(caret)
library(tidyverse)
library(doParallel)
library(kernlab)
rm(list = ls())

# Loading data -----------------------------------------------------------------
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

class.RF <- data.frame("class.multi" = class.RF[rownames(data.RF), ])
rownames(class.RF) <- rownames(data.RF)
class.multi <- factor(class.RF$class.multi)

# Model preparation ------------------------------------------------------------
# start parallel processing
cl <- makePSOCKcluster(4)
registerDoParallel(cl)

## Train/Test split ------------------------------------------------------------
# creates split based on the class labels
set.seed(505)
inTrain <- createDataPartition(class.multi, 
                               p = 3/4, list = FALSE)
# selects the corresponding samples:
trainDescr <- data.RF[inTrain, ]
saveRDS(trainDescr, "./trainDescr_multi.rds")

testDescr <- data.RF[-inTrain, ]
saveRDS(testDescr, "./testDescr_multi.rds")

trainClass <- class.multi[inTrain]
saveRDS(trainClass, "./trainClass_multi.rds")

testClass <- class.multi[-inTrain]
saveRDS(testClass, "./testClass_multi.rds")

# these should be TRUE:
dim(trainDescr)[1] == length(trainClass)
dim(testDescr)[1] == length(testClass)
# OK!

## Seeds -----------------------------------------------------------------------
seeds <- vector(mode = "list", length = 101)
for(i in 1:100) seeds[[i]] <- sample.int(1000, 196)
## For the last model:
seeds[[101]] <- sample.int(1000, 196)

## fitControl object -----------------------------------------------------------
fitControl <- trainControl(
  method = "repeatedcv",
  number = 10, 
  repeats = 10, 
  savePredictions = TRUE,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  verboseIter = TRUE,
  search = "grid",
  seeds = seeds
)



# Train SVMs -------------------------------------------------------------------
## Radial ----------------------------------------------------------------------
param.grid <- expand.grid(.sigma = 10^seq(-3, 3, length.out = 7),
                          .C = 10^seq(-3, 3, length.out = 7))

fit <- train(trainDescr, 
             trainClass,
             method = 'svmRadial',
             metric = "Kappa",
             tuneGrid = param.grid,
             trControl = fitControl)


saveRDS(fit, "./svmRadialMulti.rds")

## Radial + Sigma --------------------------------------------------------------
param.grid <- expand.grid(.sigma = 10^seq(-3, 3, length.out = 7),
                          .C = 10^seq(-3, 3, length.out = 7))
# same hyperparamete grid as in svmRadial

fit <- train(trainDescr, 
             trainClass,
             method = 'svmRadialSigma',
             metric = "Kappa",
             tuneGrid = param.grid,
             trControl = fitControl)

saveRDS(fit, "./svmRadialSigmaMulti.rds")



# Stop Cluster -----------------------------------------------------------------
stopCluster(cl)
