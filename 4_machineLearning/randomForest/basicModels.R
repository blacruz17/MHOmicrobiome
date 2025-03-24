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

# Binary classifier -----------------------------------------------------------
# start parallel processing
cl <- makePSOCKcluster(40)
registerDoParallel(cl)

## Train/Test split ----
# creates split based on the class labels
inTrain <- createDataPartition(class.binary, 
                               p = 3/4, list = FALSE)
# selects the corresponding samples:
trainDescr <- data.RF[inTrain, ]
saveRDS(trainDescr, "../results/trainDescr_binary.rds")

testDescr <- data.RF[-inTrain, ]
saveRDS(testDescr, "../results/testDescr_binary.rds")

trainClass <- class.binary[inTrain]
saveRDS(trainClass, "../results/trainClass_binary.rds")

testClass <- class.binary[-inTrain]
saveRDS(testClass, "../results/testClass_binary.rds")

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
for(i in 1:100) seeds[[i]] <- sample.int(1000, 100)
## For the last model:
seeds[[101]] <- sample.int(1000, 1)


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


## no subsampling
original <- list()

for (ntree in ntrees){
  print(ntree)
  fit <- train(trainDescr, trainClass,
               method = 'ranger',
               metric = "Kappa",
               tuneGrid = mtry.grid,
               trControl = fitControl,
               num.trees = ntree,
  	       importance = "impurity")
  
  key <- toString(ntree)
  original[[key]] <- fit
}

saveRDS(original, "../results/basicBinary_ranger.rds")


# with class weights:
class_counts <- table(class.binary)
total_samples <- sum(class_counts)
num_classes <- length(class_counts)
class_weights <- total_samples / (num_classes * class_counts)
names(class_weights) <- levels(class.binary)
print(class_weights)

weighted <- list()

for (ntree in ntrees){print(ntree)
  set.seed(505)
  fit <- train(trainDescr, trainClass,
               method = 'ranger',
               metric ="Kappa",
               tuneGrid = mtry.grid,
               trControl = fitControl,
               num.trees = ntree,
               classwt = class_weights)
  key <- toString(ntree)
  weighted[[key]] <- fit
}
saveRDS(weighted, "../results/basicBinary_ranger_weighted.rds")

# stop parallelization:
stopCluster(cl)
