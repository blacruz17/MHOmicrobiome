library(caret)
library(tidyverse)
library(doParallel)
library(kernlab)
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

# Model training ---------------------------------------------------------------
# start parallel processing
cl <- makePSOCKcluster(40)
registerDoParallel(cl)


# creates grid of mtry values:
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



## svmRadialWeights ------------------------------------------------------------
class_counts <- table(class.binary)
total_samples <- sum(class_counts)
num_classes <- length(class_counts)
class_weights <- total_samples / (num_classes * class_counts)
names(class_weights) <- levels(class.binary)
print(class_weights)

param.grid <- expand.grid(.sigma = 10^seq(-3, 3, length.out = 7),
                         .C = 10^seq(-3, 3, length.out = 7),
			 .Weight = class_weights)

fit <- train(trainDescr, 
             trainClass,
             method = 'svmRadialWeights',
             metric = "Kappa",
             tuneGrid = param.grid,
             trControl = fitControl)


saveRDS(fit, "./svmRadialWeights.rds")


# Stop Cluster ---------------------
stopCluster(cl)
