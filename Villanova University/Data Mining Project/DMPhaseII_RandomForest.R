#   MAT8480 Data Mining: Project Phase II
#   Lift Maximizers: Z.Scalyer.
#   DMPhaseII_RandomForest.R
#   Created 11/16/2019
#
#   Purpose - Parameter tuning and fitting of best Random Forest (RF) model.
#               .
#   input
#   %%%%%%%%
#   data_prep.RData
#     - data.prep: Rdata set with clean factor levels, new inputs, and recategorized inputs.
#     - inp.n: column index for non-input
#     - indx: not used
#     - split: index for 50% training data
#     - split.valid: index for 25% vaidation data
#     - split.test: index for 25% testing data
#   data_imp.RData
#     - data.imp: Rdata set with imputed data and na flag inputs.
#     - inp.n: column index for non-input
#     - indx: not used
#     - split: index for 50% training data
#     - split.valid: index for 25% vaidation data
#     - split.test: index for 25% testing data
#     - imp.parms.num: not used
#     - imp.parms.cat: not used
#
#   output
#   %%%%%%%%
#   best_RF.RData
#     - RF.Final: RF model object
#     - RFTuning: list of data.frames with validation and test fscores over
#         various tuning parameters.
#     - RF.best.VIM: data.frame for VIM of a RF fit used to subset inputs
#
#   Notes
#   %%%%%%%%
#   1. RF_quickfit.R is a minimal scrip file for updating simular fits
#        of RF.Final as updates where made to the LiftMaximizersPhaseI.R scrip.
#   2. The Gaussin clusters are no longer updated to RData files here and will
#        be conduted an independent scrip when methods for fiting scoring data
#        are provided.
#   Updates
#   %%%%%%%%
#
#%%%%%%%%%%%%#
#% Clean up %#
#%%%%%%%%%%%%#
rm(list = ls()) # clear the workspace
cat("\014")     # clear the console
graphics.off()  # clear all plots
#%%%%%%%%%%#
#% Set up %#
#%%%%%%%%%%#
getwd()                           # working directory set by R.Project: Project
par.default <- par(no.readonly=T) # save default graphic parameters
#%%%%%%%%%%%%#
#% Start up %#
#%%%%%%%%%%%%#
library(tidyverse)
library(caret)
library(randomForest)
## helper functions ----
inx <- function (data, inp.n) { # data: current dataframe; inp.n: position for non-inputs
  # numeric input indicator
  num <- sapply(data, is.numeric)
  num[inp.n] <- FALSE

  # nominal input indicator
  cat <- sapply(data, is.factor)
  cat[inp.n] <- FALSE

  # missing value indicator
  na <- sapply(data,function(x) any(is.na(x)))
  na[inp.n] <- FALSE

  data.frame(num, cat, na)
}

impute <- function(x,y) { # impute x with y
  x[is.na(x)]<-y
  x
}

FindMode <- function(x) { # find the first mode of x
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
} # Not used, will return <NA> when NA is the mode

z.mode <- function(x){ # find all modes of x
  return(names(table(x))[table(x)==max(table(x))])
}

# load data ----
#% Gaussin clusters ----
# Gaussian Mixture Model clusters made by Sebastian Tilson
# print(load("LiftMax_DB/data_imp_clusters.RData")) # load and print contents
# imp.clusters <- data.imp %>%                      # extract cluster input with EvalID for merge
#   dplyr::select(EvalID, cluster)
# rm(data.imp)                                      # remove data.frame

# load data.imp
print(load("LiftMax_DB/data_imp.RData")) # load and print contents
# data.imp <- data.imp %>%                 # update data.imp
#   inner_join(x=., y=imp.clusters,        # merge matching rows by EvalID
#              by = "EvalID")
rm(indx)                                 # remove for good house keeping practice

#% Update Data: imp ----
#% Routine: indx update %#
names(data.imp[inp.n]) # names of non-inputs
indx <- inx(data.imp, inp.n)
# names(data.imp[indx$num]) # numeric inputs
# names(data.imp[indx$cat]) # catagric inputs
# names(data.imp[indx$na]) # missing inputs

# save(inp.n, indx, split, split.valid, split.test,
#      imp.parms.num, imp.parms.cat, data.imp,  file="LiftMax_DB/data_imp.RData")

inp.n.imp = inp.n # save non-input indexes
rm(inp.n)

#% Update Data: down.imp ----
# print(load("LiftMax_DB/data_down_imp.RData")) # load RData and print contents
# data.down.imp <- data.down.imp %>%            # merge only matching EvalID
#   inner_join(y=imp.clusters,by = "EvalID")
#% Routine: indx update %#
# names(data.down.imp[inp.n]) # names of non-inputs
# indx <- inx(data.down.imp, inp.n)
# names(data.down.imp[indx$num]) # numeric inputs
# names(data.down.imp[indx$cat]) # catagric inputs
# names(data.down.imp[indx$na]) # missing inputs

# save(inp.n, indx, split.down, split.down.valid, split.down.test,
#      down.imp.parms.num, down.imp.parms.cat, data.down.imp,
#      file="LiftMax_DB/data_down_imp.RData")

# inp.n.down.imp <- inp.n # save non-input indexes
# indx.down.imp <- indx
# rm(inp.n,indx)

#% Update Data: data.prep ----
# getting inp.n, Splits, and data.prep
print(load("LiftMax_DB/data_prep.RData"))
rm(indx)
# data.prep <- data.prep  %>%                # merge only matching EvalID
#   inner_join(y=imp.clusters,by = "EvalID") # add cluster var from data_imp_clusters.RData
#% Routine: indx update %#
names(data.prep[inp.n]) # names of non-inputs
indx <- inx(data.prep, inp.n)
# names(data.prep[indx$num]) # numeric inputs
# names(data.prep[indx$cat]) # catagric inputs
# names(data.prep[indx$na]) # missing inputs

# save(inp.n, indx, split, split.valid, split.test, data.prep,  file="LiftMax_DB/data_prep.RData")

# Pass data: tree ----
data.tree <- data.prep

# update index for data.tree with inp.n in data_prep.RData
inp.n.dt <- which(colnames(data.tree) %in% colnames(data.prep[inp.n]))
names(data.tree[inp.n.dt]) # names of non-inputs for indexing

indx.dt <- inx(data.tree, inp.n.dt)
# names(data.tree[indx.dt$num])
# names(data.tree[indx.dt$cat])
# names(data.tree[indx.dt$na])

## Decision Tree ----
#%%%%%%%%%%%%%%%%%%#
#% Model building %#
#%%%%%%%%%%%%%%%%%%#
# not rejected varables
vars.dt <- -which(colnames(data.tree) %in% colnames(data.tree[inp.n.dt[-length(inp.n.dt)]]))

#% a simple tree with cp 0.001 %#
DT.001 <- rpart::rpart(formula = Book_12Mo ~ .,data = data.tree[split,vars.dt],
                       control = rpart::rpart.control(cp = 0.001))
summary(DT.001)
#
# #% scoring %#
# DT.001.predict <- predict(DT.001, data.tree[split.valid,vars.dt],
#                          type="class")
# DT.001.fscore <- caret::confusionMatrix(table(DT.001.predict,
#                                               data.tree[split.valid,]$Book_12Mo),
#                                         positive = "1")$byClass["F1"]
#
# #% list results %#
# list.fscores <- list(DT.001.fscore = DT.001.fscore)
# DT.VIMs <- list(c(DT.001[["variable.importance"]]))

#%%%%%%%%%%%%%%%%%#
#% Model pruning %#
#%%%%%%%%%%%%%%%%%#
# Model pruning on F-score USING alternative cutoff
cp.seq=DT.001$cptable[,1]
fscore<-numeric()
fscore[1]<-0  # Set root node F-score zero
for (i in 2:length(cp.seq)) {
  DT.001.prob = predict(rpart::prune(DT.001, cp=cp.seq[i]), data.tree[split.valid,vars.dt],type="prob")[,2]
  rocCurve.DT.001 <- pROC::roc(data.tree[split.valid,]$Book_12Mo, DT.001.prob, quiet=TRUE)
  DT.001Thresh <-  pROC::coords(rocCurve.DT.001, x = "best", best.method = "closest.topleft", transpose = FALSE)
  DT.001.class <- as.factor(ifelse(DT.001.prob >= DT.001Thresh$threshold, 1,0))
  fscore[i]<-caret::confusionMatrix(table(DT.001.class,data.tree[split.valid,]$Book_12Mo),
                                    positive = "1")$byClass["F1"]
}

pdf("ModelPruning_DT01.pdf")
plot(DT.001$cptable[,'nsplit']+1,fscore, type="o",
     xlab="Number of Leaves", ylab="F-score",
     main="Model pruning: DT.001")
points(x=DT.001$cptable[(max(fscore) == fscore),'nsplit']+1,
       y=max(fscore),pch=23,bg="blue")
dev.off()

#%%%%%%%%%%%%%%%#
#% Final model %#
#%%%%%%%%%%%%%%%#
# max(fscore) - fscore[2]  # 0.01619064
# DT.001.final <- rpart::prune(DT.001,cp=cp.seq[(max(fscore) == fscore)])  # max F-score=0.4647456
DT.001.final <- rpart::prune(DT.001,cp=cp.seq[2])

#% list results %#
list.fscores <- list(DT.001.fscore = fscore[2])
DT.VIMs <- list(DT.001.VIM = DT.001.final[["variable.importance"]])

#% plot tree %#
# pdf("best_DT01.pdf")
# plot(partykit::as.party(DT.001.2nd))
# dev.off()

#% Cost Matrix ----
prop <- 1/prop.table(table(data.tree[split,]$Book_12Mo))
costMatrix <- matrix(c(0,prop[2],prop[1],0), nrow=2)
costMatrix

DT.cm <- rpart::rpart(formula = Book_12Mo ~ .,data = data.tree[split,vars.dt],
                      parms=list(loss=costMatrix),
                      control=rpart::rpart.control(cp=0.001))
# summary(DT.cm)

#%%%%%%%%%%%%%%%%%#
#% Model pruning %#
#%%%%%%%%%%%%%%%%%#
# Model pruning on F-score USING alternative cutoff
cp.seq=DT.cm$cptable[,1]
fscore<-numeric()
fscore[1]<-0  # Set root node F-score zero
for (i in 2:length(cp.seq)) {
  DT.cm.prob = predict(rpart::prune(DT.cm, cp=cp.seq[i]), data.tree[split.valid,vars.dt],type="prob")[,2]
  rocCurve.DT.cm <- pROC::roc(data.tree[split.valid,]$Book_12Mo, DT.cm.prob, quiet=TRUE)
  DT.cmThresh <-  pROC::coords(rocCurve.DT.cm, x = "best", best.method = "closest.topleft", transpose = FALSE)
  DT.cm.class <- as.factor(ifelse(DT.cm.prob >= DT.cmThresh$threshold, 1,0))
  fscore[i]<-caret::confusionMatrix(table(DT.cm.class,data.tree[split.valid,]$Book_12Mo),
                             positive = "1")$byClass["F1"]
}

pdf("ModelPruning_DTcm.pdf")
plot(DT.cm$cptable[,'nsplit']+1,fscore, type="o",
     xlab="Number of Leaves", ylab="F-score",
     main="Model pruning: DT.cm")
points(x=DT.cm$cptable[(max(fscore) == fscore),'nsplit']+1,
       y=max(fscore),pch=23,bg="blue")
dev.off()

#%%%%%%%%%%%%%%%#
#% Final model %#
#%%%%%%%%%%%%%%%#
DT.cm.final <- rpart::prune(DT.cm,cp=cp.seq[(max(fscore) == fscore)])  # max F-score=0.4771969

#% list results %#
list.fscores = append(list.fscores, list(DT.cm.fscore = max(fscore)))
DT.VIMs <- append(DT.VIMs, list(DT.cm.VIM=DT.cm.final[["variable.importance"]]))

#% plot tree %#
# pdf("best_DTcm.pdf")
# plot(partykit::as.party(DT.cm.final))
# dev.off()

# Random Forest ----
#% Pass imp data
data.rf <- data.imp

#% setting index
inp.n.rf <- which(colnames(data.rf) %in% colnames(data.tree[inp.n.dt]))
names(data.rf[inp.n.rf]) # names of non-inputs for indexing

indx.rf <- inx(data.rf, inp.n.rf)
# names(data.rf[indx.rf$num])
# names(data.rf[indx.rf$cat])
# names(data.rf[indx.rf$na])

#% inputs for rf
# require(levels(input) < 33)
any(purrr::map(data.rf[,indx.rf$cat], ~length(levels(.))) > 32 ) # TRUE
# purrr::map(data.rf[,indx.rf$cat], ~length(levels(.))) # $TourDate_WeekYear

vars.rf <- -c(inp.n.rf[-length(inp.n.rf)],
             grep("^TourDate_WeekYear$", names(data.rf)))
#%%%%%%%%%%%%%%%%%%#
#% Model building %#
#%%%%%%%%%%%%%%%%%%#
#% down sampled %#
minor <- table(data.rf$Book_12Mo[split])[[2]] # stratified sample size minor = 2436 in training data
#% Parameter Tuning %#
# defalts:
#   mtry - m r.s. of perameters each split,
#          for classification floor(sqrt(ncol(data.rf[,vars1])-1))
#   ntree - n trees in forest, 500
#% ntree ----
n <- seq.int(500, to = 5000, by = 500)
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(n)){
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.rf],
                                   ntree=n[i],
                                   mtry=floor(sqrt(dim(data.rf[split, vars.rf])[2]-1)),
                                   strata= data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

n.best <- n[which.max(.5*(fscore.val+fscore.tes))] # 4000

pdf("ParameterTuning_ntree_mdefalt.pdf")
plot(n, fscore.val, pch=19 , col="blue", type="b",
     ylab="F-score",xlab="Number of Trees",
     main=paste0("randomForest: ntree \nmtry=",
                 floor(sqrt(dim(data.rf[split, vars.rf])[2]-1))))
points(n, fscore.tes, pch=19 , col="green", type="b")
points(x=n.best,
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFTuning <-list(rf.1 = data.frame(m=rep(floor(sqrt(dim(data.rf[split, vars.rf])[2]-1)),length(n)),
                                  n=n,
                                  fscore.val=fscore.val,
                                  fscore.tes=fscore.tes,
                                  fscore.mean=.5*(fscore.val+fscore.tes)))

#% mtry ----
m <- seq.int(2, to = floor(sqrt(dim(data.rf[split, vars.rf])[2]-1)+5)) # (2:14)
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(m)){
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.rf],
                                   ntree=n.best, mtry=m[i],
                                   strata= data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

m.best <- m[which.max(.5*(fscore.val+fscore.tes))] # 1500

pdf("ParameterTuning_mtry_mtreebest.pdf")
plot(m, fscore.val, pch=19 , col="blue", type="b",
     ylab="F-score",xlab="Number of Predictors considered at each split",
     main=paste0("randomForest: mtry \nntree=",n.best))
points(n, fscore.tes, pch=19 , col="green", type="b")
points(x=m.best,
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFTuning <- append(RFTuning,
                   list(rf.2 = data.frame(m=m,
                                          n=rep(n.best,length(m)),
                                          fscore.val=fscore.val,
                                          fscore.tes=fscore.tes,
                                          fscore.mean=.5*(fscore.val+fscore.tes))))

#% ntree for mbest ----
n <- seq.int(2000, to = 5000, by = 500)
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(n)){
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.rf],
                                   ntree=n[i],
                                   mtry=m.best,
                                   strata= data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

n.best <- n[which.max(.5*(fscore.val+fscore.tes))] # 3500

pdf("ParameterTuning_ntree_mbest.pdf")
plot(n, fscore.val, pch=19 , col="blue", type="b",
     ylab="F-score",xlab="Number of Trees",
     main=paste0("randomForest: ntree \nmtry=",
                 m.best))
points(n, fscore.tes, pch=19 , col="green", type="b")
points(x=n.best,
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFTuning <- append(RFTuning,
                   list(rf.3 = data.frame(m=rep(m.best,length(n)),
                                  n=n,
                                  fscore.val=fscore.val,
                                  fscore.tes=fscore.tes,
                                  fscore.mean=.5*(fscore.val+fscore.tes))))

#%%%%%%%%%%%%%%#
#% best model %#
#%%%%%%%%%%%%%%#
#% for VIMs %#
set.seed(77012)
RF.best <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.rf],
                                      ntree=3500, mtry=9, # n = 3500, m = 9
                                      strata= data.rf$Book_12Mo[split],
                                      sampsize=c(minor,minor),
                                      importance=T)

RF.best.VIM <- data.frame(Input = row.names(RF.best[["importance"]]),
                          MDA.Rank = rank(RF.best[["importance"]][,3]), # MeanDecreaseAccuracy
                          MDG.Rank = rank(RF.best[["importance"]][,4])) %>%  # MeanDecreaseGini
  mutate(Rank.Product = (MDA.Rank * MDG.Rank)) %>%
  arrange(desc(Rank.Product))

RF.best.vclass <- predict(RF.best, newdata=data.rf[split.valid,], type="class")
RF.best.vfscore <- caret::confusionMatrix(table(RF.best.vclass,data.rf[split.valid,]$Book_12Mo),
                                        positive = "1")$byClass["F1"]
RF.best.tclass <- predict(RF.best, newdata=data.rf[split.test,], type="class")
RF.best.tfscore <- caret::confusionMatrix(table(RF.best.tclass,data.rf[split.test,]$Book_12Mo),
                                        positive = "1")$byClass["F1"]
RF.best.vfscore # 0.5422285 (0.5480427 w/o VIM)
RF.best.tfscore # 0.5382698 (0.5383523 w/o VIM)


#% subset inputs on fixed m,n ----
p <- seq.int(dim(data.rf[split, vars.rf])[2]-1,
             to=dim(data.rf[split, vars.rf])[2]-11) # (117:107)
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(p)){
  vars.subset <- which(colnames(data.rf) %in% RF.best.VIM$Input[1:p[i]] | colnames(data.rf) == "Book_12Mo")
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.subset],
                                   ntree=n.best, mtry=m.best,
                                   strata= data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

pdf("ParameterTuning_subset.pdf")
plot(seq.int(0,length(p)-1), fscore.val,
     pch=19 , col="blue", type="b",
     ylim=c(min(fscore.tes), max(fscore.val)),
     ylab="F-score",xlab="Number of Inputs Excluded",
     main=paste0("randomForest: Subset Inputs by Importance \nmtry=",m.best,",ntree=",n.best))
points(seq.int(0,length(p)-1), fscore.tes, pch=19 , col="green", type="b")
points(x=which.max(.5*(fscore.val+fscore.tes))-1,
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFTuning <- append(RFTuning,
                   list(rf.4 = data.frame(m=rep(m.best,length(p)),
                                          n=rep(n.best,length(p)),
                                          p=p,
                                          fscore.val=fscore.val,
                                          fscore.tes=fscore.tes,
                                          fscore.mean=.5*(fscore.val+fscore.tes))))

p.lim <- p[7] # 111

vars.subset <- which(colnames(data.rf) %in% RF.best.VIM$Input[1:p.lim] | colnames(data.rf) == "Book_12Mo")

#% mtry on subset for n.best ----
m <- seq.int(4, to = 12) #
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(m)){
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.subset],
                                   ntree=n.best, mtry=m[i],
                                   strata=data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

pdf("ParameterTuning_mtry_subset.pdf")
plot(m, fscore.val, pch=19 , col="blue", type="b",
     ylim=c(min(fscore.tes), max(fscore.val)),
     ylab="F-score",xlab="Number of Predictors considered at each split",
     main=paste0("randomForest: mtry on subset \nntree=",n.best, "p.lim=",p.lim))
points(m, fscore.tes, pch=19 , col="green", type="b")
points(x=m[which.max(.5*(fscore.val+fscore.tes))],
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFTuning <- append(RFTuning,
                   list(rf.5 = data.frame(m=m,
                                          n=rep(n.best,length(m)),
                                          p=p.lim,
                                          fscore.val=fscore.val,
                                          fscore.tes=fscore.tes,
                                          fscore.mean=.5*(fscore.val+fscore.tes))))

#% ntree on !subset for m=7 ----
n <- seq.int(1000, to = 4000, by = 500)
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(n)){
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.rf], # should have been vars.split
                                   ntree=n[i],
                                   mtry=7,
                                   strata= data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

pdf("ParameterTuning_ntree_subset_m7.pdf")
plot(n, fscore.val, pch=19 , col="blue", type="b",
     ylab="F-score",xlab="Number of Trees",
     ylim=c(min(fscore.tes), max(fscore.val)),
     main="randomForest: ntree on subset \nmtry=7")
points(n, fscore.tes, pch=19 , col="green", type="b")
points(x=n[which.max(.5*(fscore.val+fscore.tes))],
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFTuning <- append(RFTuning,
                   list(rf.6 = data.frame(m=rep(7,length(n)),
                                          n=n,
                                          p=p[1], # edit from p.lim
                                          fscore.val=fscore.val,
                                          fscore.tes=fscore.tes,
                                          fscore.mean=.5*(fscore.val+fscore.tes))))
# 7 2500 111 for mtry 2:7

#% ntree on !subset for m=12 ----
n <- seq.int(1000, to = 5000, by = 500)
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(n)){
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.rf], # AGAIN vars.split was my intent
                                   ntree=n[i],
                                   mtry=12,
                                   strata= data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

pdf("ParameterTuning_ntree_subset_m12.pdf")
plot(n, fscore.val, pch=19 , col="blue", type="b",
     ylab="F-score",xlab="Number of Trees",
     ylim=c(min(fscore.tes), max(fscore.val)),
     main="randomForest: ntree on subset \nmtry=12")
points(n, fscore.tes, pch=19 , col="green", type="b")
points(x=n[which.max(.5*(fscore.val+fscore.tes))],
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFTuning <- append(RFTuning,
                   list(rf.7 = data.frame(m=rep(12,length(n)), # fucked up and used rf.7 twice
                                          n=n,
                                          p=p[1], # edit form p.lim
                                          fscore.val=fscore.val,
                                          fscore.tes=fscore.tes,
                                          fscore.mean=.5*(fscore.val+fscore.tes))))

#% ntree on subset for m=7 ----
n <- seq.int(1000, to = 4000, by = 500)
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(n)){
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.subset],
                                   ntree=n[i],
                                   mtry=7,
                                   strata= data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

pdf("ParameterTuning_ntree_plim_m7.pdf")
plot(n, fscore.val, pch=19 , col="blue", type="b",
     ylab="F-score",xlab="Number of Trees",
     ylim=c(min(fscore.tes), max(fscore.val)),
     main="randomForest: ntree on subset \nmtry=7")
points(n, fscore.tes, pch=19 , col="green", type="b")
points(x=n[which.max(.5*(fscore.val+fscore.tes))],
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFTuning <- append(RFTuning,
                   list(rf.8 = data.frame(m=rep(7,length(n)),
                                          n=n,
                                          p=p.lim,
                                          fscore.val=fscore.val,
                                          fscore.tes=fscore.tes,
                                          fscore.mean=.5*(fscore.val+fscore.tes))))
# we will see

#% ntree on plim for m=12 ----
n <- seq.int(1000, to = 5000, by = 500)
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(n)){
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.subset],
                                   ntree=n[i],
                                   mtry=12,
                                   strata= data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

pdf("ParameterTuning_ntree_plim_m12.pdf")
plot(n, fscore.val, pch=19 , col="blue", type="b",
     ylab="F-score",xlab="Number of Trees",
     ylim=c(min(fscore.tes), max(fscore.val)),
     main="randomForest: ntree on subset \nmtry=12")
points(n, fscore.tes, pch=19 , col="green", type="b")
points(x=n[which.max(.5*(fscore.val+fscore.tes))],
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFTuning <- append(RFTuning,
                   list(rf.9 = data.frame(m=rep(12,length(n)),
                                          n=n,
                                          p=p.lim,
                                          fscore.val=fscore.val,
                                          fscore.tes=fscore.tes,
                                          fscore.mean=.5*(fscore.val+fscore.tes))))
# fixing recrods
RFtuning <- RFTuning

RFtuning[[6]]$p <- rep(p[1],7)
RFtuning[[7]]<- NULL
# m,n,p = 7 1000 111

#% mtry on plim for n=1000 ----
m <- seq.int(4, to = 9)
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(m)){
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.subset],
                                   ntree=1000,
                                   mtry=m[i],
                                   strata= data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

pdf("ParameterTuning_mtry_plim_n1000.pdf")
plot(m, fscore.val, pch=19 , col="blue", type="b",
     ylab="F-score",xlab="Number of Predictors considered at each split",
     ylim=c(min(fscore.tes), max(fscore.val)),
     main="randomForest: mtry on subset \nntree=1000")
points(m, fscore.tes, pch=19 , col="green", type="b")
points(x=m[which.max(.5*(fscore.val+fscore.tes))],
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFtuning <- append(RFtuning,
                   list(rf.10 = data.frame(m=m,
                                          n=rep(1000,length(m)),
                                          p=p.lim,
                                          fscore.val=fscore.val,
                                          fscore.tes=fscore.tes,
                                          fscore.mean=.5*(fscore.val+fscore.tes))))
#% mtry on plim for n=1500 ----
m <- seq.int(4, to = 9)
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(m)){
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.subset],
                                   ntree=1500,
                                   mtry=m[i],
                                   strata= data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

pdf("ParameterTuning_mtry_plim_n1500.pdf")
plot(m, fscore.val, pch=19 , col="blue", type="b",
     ylab="F-score",xlab="Number of Predictors considered at each split",
     ylim=c(min(fscore.tes), max(fscore.val)),
     main="randomForest: mtry on subset \np=",p.lim," ntree=1500")
points(m, fscore.tes, pch=19 , col="green", type="b")
points(x=m[which.max(.5*(fscore.val+fscore.tes))],
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFtuning <- append(RFtuning,
                   list(rf.11 = data.frame(m=m,
                                           n=rep(1500,length(m)),
                                           p=p.lim,
                                           fscore.val=fscore.val,
                                           fscore.tes=fscore.tes,
                                           fscore.mean=.5*(fscore.val+fscore.tes))))

#% subset inputs on m,n=7,1000 ----
p <- seq.int(dim(data.rf[split, vars.rf])[2]-1,
             to=dim(data.rf[split, vars.rf])[2]-15) # (117:103)
fscore.val <- numeric()
fscore.tes <- numeric()
for(i in 1:length(p)){
  vars.subset <- which(colnames(data.rf) %in% RF.best.VIM$Input[1:p[i]] | colnames(data.rf) == "Book_12Mo")
  set.seed(77012)
  rf <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.subset],
                                   ntree=1000, mtry=7,
                                   strata= data.rf$Book_12Mo[split],
                                   sampsize=c(minor,minor))
  rf.vclass <- predict(rf, newdata=data.rf[split.valid,], type="class")
  fscore.val[i] <- caret::confusionMatrix(table(rf.vclass,data.rf[split.valid,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
  rf.tclass <- predict(rf, newdata=data.rf[split.test,], type="class")
  fscore.tes[i] <- caret::confusionMatrix(table(rf.tclass,data.rf[split.test,]$Book_12Mo),
                                          positive = "1")$byClass["F1"]
}

pdf("ParameterTuning_subset_n1000_m7.pdf")
plot(seq.int(0,length(p)-1), fscore.val,
     pch=19 , col="blue", type="b",
     ylim=c(min(fscore.tes), max(fscore.val)),
     ylab="F-score",xlab="Number of Inputs Excluded",
     main="randomForest: Subset Inputs by Importance \nmtry=7, ntree= 1000")
points(seq.int(0,length(p)-1), fscore.tes, pch=19 , col="green", type="b")
points(x=which.max(.5*(fscore.val+fscore.tes))-1,
       y=max(.5*(fscore.val+fscore.tes)),pch=23,bg="red")
dev.off()

RFTuning <- append(RFTuning,
                   list(rf.12 = data.frame(m=rep(7,length(p)),
                                           n=rep(1000,length(p)),
                                           p=p,
                                           fscore.val=fscore.val,
                                           fscore.tes=fscore.tes,
                                           fscore.mean=.5*(fscore.val+fscore.tes))))

vars.subset <- which(colnames(data.rf) %in% RF.best.VIM$Input[1:p.lim])
# names(data.rf[split, vars.subset[c(-57,-length(vars.subset))]])
# vars.subset <- vars.subset[c(-57,-length(vars.subset))]
#%%%%%%%%%%%%%%%#
#% Final model %#
#%%%%%%%%%%%%%%%#
#% Final model ----
vars.subset <- which(colnames(data.rf) %in% RF.best.VIM[1:111,1] |
                       colnames(data.rf) == "Book_12Mo")

set.seed(77012)
RF.Final <- randomForest::randomForest(Book_12Mo ~., data=data.rf[split, vars.rf],
                                       ntree=1000, mtry=7, # n = 1500, m = 6
                                       strata= data.rf$Book_12Mo[split],
                                       sampsize=c(minor,minor))
# rf.final<- RF.Final
# RF.Final <- caret::train(x=data.rf[split, vars.subset],y=data.rf[split,]$Book_12Mo, method = "rf",
#                          ntree=1000, mtry=7, # n = 1500, m = 6
#                          strata= data.rf$Book_12Mo[split],
#                          sampsize=c(minor,minor))


#% scoring %#
RF.Final.class.valid <- predict(RF.Final, newdata=data.rf[split.valid,], type="class")
RF.Final.fscore.valid <- caret::confusionMatrix(table(RF.Final.class.valid,data.rf[split.valid,]$Book_12Mo),
                                           positive = "1")$byClass["F1"]
RF.Final.fscore.valid # 0.5425799 {0.5333815 (0.5432452, 0.5421903,0.5489209)}

RF.Final.class.test <- predict(RF.Final, newdata=data.rf[split.test,], type="class")
RF.Final.fscore.test <- caret::confusionMatrix(table(RF.Final.class.test,data.rf[split.test,]$Book_12Mo),
                                                positive = "1")$byClass["F1"]
RF.Final.fscore.test # 0.5285611 {0.5302594 (0.5344333, 0.5398198, 0.5318841)}

# print(RF)
# plot(RF.Final)
# OOB error is OOB MISC (in black)
# gree is the change in class 1 error
#   do to few obs in training data with target 1 (minor = sum(RF[["confusion"]][2,1:2]))
#   there is many misclassified target 0 as target 1 (RF[["confusion"]][2,1])
#   why error rate is high (RF[["confusion"]][2,3])
#     not a problem/understandable
# red is the change in class error of 0

# SAVE Best ----
# save(RFTuning, RF.best.VIM, RF.Final, file="best_RF.RData")
save(RF.best.VIM, RF.Final, file="best_RF.RData")

# write.csv(RF.best.VIM, file="RF_VIM.csv")
