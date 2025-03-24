library(caret)
library(tidyverse)
library(doParallel)
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

# Model preparation ------------------------------------------------------------
# start parallel processing
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

## Train/Test split ----
trainDescr <- readRDS("../results/trainDescr_binary.rds")
trainClass <- readRDS("../results/trainClass_binary.rds")
testDescr <- readRDS("../results/testDescr_binary.rds")
testClass <- readRDS("../results/testClass_binary.rds")

# these should be TRUE:
dim(trainDescr)[1] == length(trainClass)
dim(testDescr)[1] == length(testClass)
# OK!

## Model training -----
# creates grid of mtry values:
mtry.grid <- expand.grid(.mtry = floor(seq(sqrt(ncol(trainDescr)), 
                                           ncol(trainDescr), 
                                           length.out = 5)),
                         .splitrule = "gini",
                         .min.node.size = seq(1,10))
# creates sequence of ntree values:
ntrees <- c(100, 500, 1000, 1500, 2000)

seeds <- vector(mode = "list", length = 101)
for(i in 1:100) seeds[[i]] <- sample.int(1000, 196)
## For the last model:
seeds[[101]] <- sample.int(1000, 196)

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


## Upsampling ------------------------------------------------------------------
fitControl$sampling <- "up"

upsampled <- list()

for (ntree in ntrees){
  print(ntree)
  fit <- train(trainDescr, trainClass,
               method = 'ranger',
               metric = "Kappa",
               tuneGrid = mtry.grid,
               trControl = fitControl,
               num.trees = ntree)
  
  key <- toString(ntree)
  upsampled[[key]] <- fit
}

saveRDS(upsampled, "../results/upsampledBinary_ranger.rds")

## Downsampling ----------------------------------------------------------------
fitControl$sampling <- "down"

downsampled <- list()

for (ntree in ntrees){
  print(ntree)
  fit <- train(trainDescr, trainClass,
               method = 'ranger',
               metric = "Kappa",
               tuneGrid = mtry.grid,
               trControl = fitControl,
               num.trees = ntree)
  
  key <- toString(ntree)
  downsampled[[key]] <- fit
}

saveRDS(downsampled, "../results/downsampledBinary_ranger.rds")


# stop parallelization ---------------------------------------------------------
stopCluster(cl)
