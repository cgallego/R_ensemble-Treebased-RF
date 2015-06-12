###################################################
### code Read and partition data 
###################################################
setwd("Z:/Cristina/MassNonmass/Section1 - ExperimentsUpToDate/experimentsRadiologypaper-revision/Tree-based-RF/ensemble-Treebased-RF")
library("RSQLite")

rpart_inputdata <- function(subdata) {
  sqlite <- dbDriver("SQLite")
  conn <- dbConnect(sqlite, "localData.db")
  
  # 2) all T1W features
  lesionsQuery <- dbGetQuery(conn, "SELECT *
                             FROM  lesion
                             INNER JOIN f_dynamic ON (lesion.lesion_id = f_dynamic.lesion_id)
                             INNER JOIN f_morphology ON (lesion.lesion_id = f_morphology.lesion_id)
                             INNER JOIN f_texture ON (lesion.lesion_id = f_texture.lesion_id)")
  
  # prune entries and extract feature subsets
  # corresponds to 5 entries lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
  lesionsfields = names(lesionsQuery[c(1,23,24,3,5,27:60,63:81,84:107)])
  lesioninfo = lesionsQuery[c(1,23,24,3,5)]
  dynfeatures = lesionsQuery[c(1,23,27:60)]
  morphofeatures = lesionsQuery[c(1,23,63:81)]
  texfeatures = lesionsQuery[c(1,23,84:107)]
  
  # combine all features
  allfeatures = cbind(dynfeatures, morphofeatures[3:ncol(morphofeatures)], texfeatures[3:ncol(morphofeatures)])  
  
  if(subdata=="mass"){
    # organized the data by subdata
    M<-subset(allfeatures, lesion_label=="massB" | lesion_label=="massM")
    ifelse( M$lesion_label == "massB", "NC", "C") -> M$lesion_label
    allfeatures = M
  }
  if(subdata=="nonmass"){
    # organized the data by subdata
    N<-subset(allfeatures, lesion_label=="nonmassB" | lesion_label=="nonmassM")
    ifelse( N$lesion_label == "nonmassB", "NC", "C") -> N$lesion_label
    allfeatures = N
  }
  if(subdata=="stage1"){
    # organized the data by subdata
    M<-subset(allfeatures, lesion_label=="massB" | lesion_label=="massM")
    ifelse( M$lesion_label == "massB", "mass", "mass") -> M$lesion_label
    N<-subset(allfeatures, lesion_label=="nonmassB" | lesion_label=="nonmassM")
    ifelse( N$lesion_label == "nonmassB", "nonmass", "nonmass") -> N$lesion_label
    allfeatures = data.frame(rbind(M,N)) 
  }
  if(subdata=="oneshot"){
    # organized the data by subdata
    M<-subset(allfeatures, lesion_label=="massB" | lesion_label=="massM")
    ifelse( M$lesion_label == "massB", "NC", "C") -> M$lesion_label
    N<-subset(allfeatures, lesion_label=="nonmassB" | lesion_label=="nonmassM")
    ifelse( N$lesion_label == "nonmassB", "NC", "C") -> N$lesion_label
    allfeatures = data.frame(rbind(M,N)) 
  }
  # procees data
  allfeatures$lesion_label <- as.factor(allfeatures$lesion_label)
  allfeatures$peakCr_inside <- as.factor(allfeatures$peakCr_inside)
  allfeatures$peakVr_inside <- as.factor(allfeatures$peakVr_inside)
  allfeatures$peakCr_countor <- as.factor(allfeatures$peakCr_countor)
  allfeatures$peakVr_countor <- as.factor(allfeatures$peakVr_countor)
  allfeatures$k_Max_Margin_Grad <- as.factor(allfeatures$k_Max_Margin_Grad)
  allfeatures$max_RGH_mean_k <- as.factor(allfeatures$max_RGH_mean_k)
  allfeatures$max_RGH_var_k <- as.factor(allfeatures$max_RGH_var_k)
  
  output <- allfeatures
  return(output)
}

###################################################
### code to create a cross-validation set up: 
### cvfoldk = number of cv folds typically 5 or 10
### out: particvfoldK = all cv-K ids
###################################################
library(MASS)
library(caret)

cvfold_partition <- function(dat, cvfoldK){
  ndat = nrow(dat)
  outcomesetDi  <- dat$lesion_label
  #For multiple k-fold cross-validation, completely independent folds are created.
  #when y is a factor in an attempt to balance the class distributions within the splits.
  #The names of the list objects will denote the fold membership using the pattern 
  #"Foldi.Repj" meaning the ith section (of k) of the jth cross-validation set (of times).
  partitionsetDi <- createFolds(y = outcomesetDi, ## the outcome data are needed
                                k = cvfoldK, ## The percentage of data in the training set
                                list = TRUE) ## The format of the results. 
  return(partitionsetDi)
}

###################################################
### code to sample kparti from a cross-validation set up: 
### kparti = k fold to exclude
### outs: cvTrainsetD, cvTestsetD
###################################################
kparti_sample <- function(dat, particvfoldK, cvfoldK, kparti){
  allparti = 1:cvfoldK
  allbutkparti = allparti[-kparti]
  cvfoldadd = c()
  for(i in 1:length(allbutkparti)){
    kadd = allbutkparti[i]
    cvfoldadd = c(cvfoldadd, particvfoldK[[kadd]])
  }
  # partition data
  cvTrainsetD <-  dat[ cvfoldadd ,]
  cvTestsetD <-   dat[-cvfoldadd ,]
  
  output <- list(cvTrainsetD=cvTrainsetD, cvTestsetD=cvTestsetD)
  return(output)
}

###################################################
### code Feature selection: 
### Boruta, cvfold, 
###################################################
library(Boruta)
require(data.table)
require(ggplot2)

# function to produce correlation coefficients on pair plots
panel.cor <- function(x, y, digits = 2, cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

subset_select <- function(setTrain){
  featsel_boruta <-Boruta(lesion_label ~., data=setTrain[,2:ncol(setTrain)], doTrace=2, ntree=1000)
  print(featsel_boruta)
  plot(featsel_boruta)
  
  relevant <- featsel_boruta$finalDecision[featsel_boruta$finalDecision == "Confirmed"]
  relevant_features = setTrain[c(names(relevant))]
  tentative <- featsel_boruta$finalDecision[featsel_boruta$finalDecision == "Tentative"]
  tentative_features = setTrain[c(names(tentative))]
  sel_features = cbind(setTrain[c(1,2)], relevant_features, tentative_features)
  
  super.sym <- trellis.par.get("superpose.symbol")
  ## pair plots for reatures
  setTrainrelevant = setTrain[c(names(relevant))]
  pairs(relevant_features, 
        upper.panel = panel.cor,
        pch = super.sym$pch[1:2],
        col = super.sym$col[1:2],
        text = list(levels(setTrainrelevant$lesion_label)),
        main = "Relevant features")
  
  return(sel_features)
}


###################################################
### code forest Train: 
### parameters, T= # of trees, D= tree depth, dat
###################################################
library(klaR)
library(rpart)
library(rpart.plot)

# bagged training was introduced as a way of reducing possible overfitting and 
# improving the generalization capabilities of random forests. 
# The idea is to train each tree in a forest on a different training subset, sampled at random from the same labeled database. 
rpart_looforestTrain <- function(T, D, dat) {
  # set control
  fitparm = rpart.control(maxdepth = D, minsplit = 5, minbucket = 4, cp = 0.00001,  xval = 5,
                          maxcompete = 0, maxsurrogate = 0, usesurrogate = 0, surrogatestyle = 0)
  
  # init forest
  forest = list()
  for (t in 1:T){
    #cat("Tree # ", t, "\n")  
    
    # build bagged trees from a bootstrap sample of trainSetD
    setD = dat[sample(1:nrow(dat), nrow(dat), replace=TRUE),]
    
    # find subsample of var
    # when training the ith tree we only make available a small random subset 
    subvar = sample(2:ncol(setD), sqrt(ncol(setD)-1), replace = FALSE)
    subfeat = colnames(setD)[subvar]
    
    # train tree
    treedata <- rpart(paste("lesion_label ~ ", paste(subfeat, collapse= "+")), 
                      method = "class", data = setD, control=fitparm)
    
    # display the probability per class of observations in the node (conditioned on the node, sum across a     node is 1) plus the percentage of observations in the node. 
    if (T==1){
      print(treedata)
      prp(treedata, type=2, digits=3, extra = 102, under=TRUE, nn=TRUE, col="black", 
          box.col=rainbow(2)[2], varlen=0, faclen=0, branch.type=0, gap=0, cex=.7,
          fallen.leaves=TRUE) # use fallen.leaves=TRUE, to plot at bottom  
    }  
    
    # append
    forest <- append(forest, list(tree = treedata))    
  }
  
  output <- list(forest=forest)
  return(output)
}

###################################################
### code forest Test: 
### parameters, T= # of trees, forest, TrainsetD, TestsetD
###################################################
library(pROC)
rpart_looforestTest <- function(T, TrainsetD, TestsetD, forest) {
  
  fclasspotrain=list()
  for (t in 1:T){
    # Calcultate posterior Probabilities on grid points
    temp <- predict(forest[t]$tree, newdata = TrainsetD) #
    fclasspotrain <- append(fclasspotrain, list(cpo = temp))
  }
  
  # run testing cases
  fclasspotest=list()
  for (t in 1:T){
    # Calcultate posterior Probabilities on grid points
    temp <- predict(forest[t]$tree, newdata = TestsetD) #
    fclasspotest <- append(fclasspotest, list(cpo = temp))
  }
  
  # performance on Train/Test set separately
  # extract ensamble class probabilities (when T > 1)
  trainpts = fclasspotrain[1]$cpo
  testpts = fclasspotest[1]$cpo
  # init ensample class posteriors
  enclasspotrain <- matrix(, nrow = nrow(as.data.frame(trainpts)), ncol = 2)
  enclasspotest <- matrix(, nrow = nrow(as.data.frame(testpts)), ncol = 2)
  enclasspotrain[,1] = fclasspotrain[1]$cpo[,1]
  enclasspotest[,1] = fclasspotest[1]$cpo[,1]
  enclasspotrain[,2] = fclasspotrain[1]$cpo[,2]
  enclasspotest[,2] = fclasspotest[1]$cpo[,2]
  if(T>=2){
    for (t in 2:T){
      #train
      enclasspotrain[,1] = enclasspotrain[,1]+fclasspotrain[t]$cpo[,1]
      enclasspotrain[,2] = enclasspotrain[,2]+fclasspotrain[t]$cpo[,2]
      #test
      enclasspotest[,1] = enclasspotest[,1]+fclasspotest[t]$cpo[,1]
      enclasspotest[,2] = enclasspotest[,2]+fclasspotest[t]$cpo[,2]
    }
  }
  # majority voting averaging
  enclasspotrain = (1/T)*enclasspotrain
  enclasspotest = (1/T)*enclasspotest
  
  # on training
  classes = levels(TrainsetD$lesion_label)
  trainprob = data.frame(C1=enclasspotrain[,1],
                         C2=enclasspotrain[,2],
                         pred=classes[apply(enclasspotrain, 1, which.max)], 
                         obs=TrainsetD$lesion_label)
  colnames(trainprob)[1:2] <- classes
  pred=as.factor(apply(enclasspotrain, 1, which.max))
  levels(pred) = levels(as.factor(unclass(TrainsetD$lesion_label)))
  perf_train = confusionMatrix(pred, as.factor(unclass(TrainsetD$lesion_label)))
  #print(perf_train)
  
  # on testing
  testprob = data.frame(C1=enclasspotest[,1],
                        C2=enclasspotest[,2],
                        pred=classes[apply(enclasspotest, 1, which.max)], 
                        obs=TestsetD$lesion_label)
  colnames(testprob)[1:2] <- classes
  pred=as.factor(apply(enclasspotest, 1, which.max))
  levels(pred) = levels(as.factor(unclass(TrainsetD$lesion_label)))
  pred[1]=as.factor(apply(enclasspotest, 1, which.max))
  
  groundT = as.factor(unclass(TestsetD$lesion_label))
  levels(groundT) = levels(as.factor(unclass(TrainsetD$lesion_label)))
  groundT[1] = as.factor(unclass(TestsetD$lesion_label))
  
  perf_test = confusionMatrix(pred, groundT)
  #print(perf_test)  
  
  output <- list(etrain = perf_train, etest=perf_test, trainprob=trainprob, testprob=testprob)
  return(output)
}



###################################################
### code for plotting perfm results: 
###statsAU
###################################################
create_ensemble <- function(dat, particvfoldK, cvK){
  #inint
  ensemblegrdperf=list()
  for(r in 1:cvK){
    ## pick one of cvfold for held-out test, train on the rest
    kparti_setdata = kparti_sample(dat, particvfoldK, cvK, r)
    
    # Boruta on $cvTrainsetD
    selfeatures_kfold = subset_select(kparti_setdata$cvTrainsetD)
    names(selfeatures_kfold)
    
    ###################################################
    # create grid of evaluation points
    gT = c(5,10,30,60,100,250,500,750) 
    gD = c(2,5,10) 
    grd <- expand.grid(x = gD, y = gT)
    
    ###################################################
    # for oneshot
    grdperf = data.frame(grd)
    grdperf$acuTrain =0
    grdperf$rocTrain =0
    grdperf$senTrain =0
    grdperf$speTrain =0
    
    grdperf$acuTest =0
    grdperf$rocTest =0
    grdperf$senTest =0
    grdperf$speTest =0
    
    for(k in 1:nrow(grd)){
      D=grd[k,1]
      T=grd[k,2]
      # Build in l
      cat("D: ", D, "T: ", T, "\n")
      TrainsetD <-  kparti_setdata$cvTrainsetD[c(names(selfeatures_kfold))]
      TestsetD <-  kparti_setdata$cvTestsetD[c(names(selfeatures_kfold))]
      fit <- rpart_looforestTrain(T, D, TrainsetD[c(2:ncol(TrainsetD))])
      # # predict
      perf <- rpart_looforestTest(T, TrainsetD[c(2:ncol(TrainsetD))], TestsetD[c(2:ncol(TestsetD))], fit$forest)
      # for train
      ROCF_train <- plot.roc(perf$trainprob$obs, perf$trainprob$C, col="#000086", main=paste0("ROC T=",T," D=",D," cv=",r))
      print(ROCF_train$auc)
      # collect data
      grdperf$acuTrain[k] = grdperf$acuTrain[k]+as.numeric(perf$etrain$overall[1])
      grdperf$rocTrain[k] = grdperf$rocTrain[k]+as.numeric(ROCF_train$auc)
      grdperf$senTrain[k] = grdperf$senTrain[k]+as.numeric(perf$etrain$byClass[1])
      grdperf$speTrain[k] = grdperf$speTrain[k]+as.numeric(perf$etrain$byClass[2])
      # for test
      par(new=TRUE)
      ROCF_test <- plot.roc(perf$testprob$obs, perf$testprob$C, col="#860000", main=paste0("ROC T=",T," D=",D," cv=",r))
      legend("bottomright", legend = c("train", "test"), col = c("#000086", "#860000"),lwd = 2)
      print(ROCF_test$auc)
      # collect data
      grdperf$acuTest[k] = grdperf$acuTest[k]+as.numeric(perf$etest$overall[1])
      grdperf$rocTest[k] = grdperf$rocTest[k]+as.numeric(ROCF_test$auc)
      grdperf$senTest[k] = grdperf$senTest[k]+as.numeric(perf$etest$byClass[1])
      grdperf$speTest[k] = grdperf$speTest[k]+as.numeric(perf$etest$byClass[2])
    }
    print(grdperf)
    ensemblegrdperf <- append(ensemblegrdperf, list(grdperf = grdperf))
  }
  return(ensemblegrdperf)
}    


surface_forestperfm <- function(grdperf){
  library(gridExtra) 
  library(base)
  library(lattice)
  
  
  graphlist<-list()
  count <- 1
  # design acuTrain
  z = grdperf$acuTrain
  gD=unique(grdperf$x)
  gT=unique(grdperf$y)
  dim(z) <- c(length(gD), length(gT))
  w1 <- wireframe(z, gD,gT,  box = FALSE,
                  xlab = "Depth of trees (D)", ylab = "Number of trees (T)",
                  main = "Influence of forest parameters on Accuracy Train",
                  drape = TRUE,
                  colorkey = TRUE,
                  light.source = c(10,0,10), 
                  col.regions = colorRampPalette(c("red", "blue"))(100),
                  screen = list(z = 30, x = -60))
  graphlist[[count]]<-w1
  count <- count+1
  
  # design rocTrain
  z = grdperf$rocTrain
  dim(z) <- c(length(gD), length(gT))
  w2 <- wireframe(z, gD,gT,  box = FALSE,
                  xlab = "Depth of trees (D)", ylab = "Number of trees (T)",
                  main = "Influence of forest parameters on ROC Train",
                  drape = TRUE,
                  colorkey = TRUE,
                  light.source = c(10,0,10), 
                  col.regions = colorRampPalette(c("red", "blue"))(100),
                  screen = list(z = 30, x = -60))
  graphlist[[count]]<-w2
  count <- count+1
  
  # design acuTest
  z = grdperf$acuTest
  dim(z) <- c(length(gD), length(gT))
  w3 <- wireframe(z, gD,gT,  box = FALSE,
                  xlab = "Depth of trees (D)", ylab = "Number of trees (T)",
                  main = "Influence of forest parameters on Accuracy Test",
                  drape = TRUE,
                  colorkey = TRUE,
                  light.source = c(10,0,10), 
                  col.regions = colorRampPalette(c("red", "blue"))(100),
                  screen = list(z = 30, x = -60))
  graphlist[[count]]<-w3
  count <- count+1
  
  # design rocTest
  z = grdperf$rocTest
  dim(z) <- c(length(gD), length(gT))
  w4 <- wireframe(z, gD,gT,  box = FALSE,
                  xlab = "Depth of trees (D)", ylab = "Number of trees (T)",
                  main = "Influence of forest parameters on ROC Test",
                  drape = TRUE,
                  colorkey = TRUE,
                  light.source = c(10,0,10), 
                  col.regions = colorRampPalette(c("red", "blue"))(100),
                  screen = list(z = 30, x = -60))
  graphlist[[count]]<-w4
  count <- count+1
  
  
  # finally plot in grid
  do.call("grid.arrange",c(graphlist,ncol=2))
}

###################################################
# read mass features
massdat = rpart_inputdata(subdata="mass")
# read nonmass features
nonmassdat = rpart_inputdata(subdata="nonmass")
# read stage1 features
stage1dat = rpart_inputdata(subdata="stage1")
# read oneshot features
oneshotdat = rpart_inputdata(subdata="oneshot")

###################################################
## create CV
# cvK = 5
# # run for mass
# particvfoldK = cvfold_partition(massdat, cvK)
# massensemblegrdperf = create_ensemble(massdat, particvfoldK, cvK)
# # plot
# accum_massgrdperf = massensemblegrdperf[1]$grdperf
# for(k in 1:3){
#   accum_massgrdperf = accum_massgrdperf+massensemblegrdperf[k]$grdperf
# }
# cvKmassgrdperf = accum_massgrdperf/3
# surface_forestperfm(cvKmassgrdperf)
# print(cvKmassgrdperf)

###################################################
# read nonmass features
nonmassdat = rpart_inputdata(subdata="nonmass")
## create CV
cvK = 5
# run for nonmass
particvfoldK = cvfold_partition(nonmassdat, cvK)
nonmassensemblegrdperf = create_ensemble(nonmassdat, particvfoldK, cvK)
# plot
accum_nonmassgrdperf = nonmassensemblegrdperf[1]$grdperf
for(k in 2:3){
  accum_nonmassgrdperf = accum_nonmassgrdperf+nonmassensemblegrdperf[k]$grdperf
}
cvKnonmassgrdperf = accum_nonmassgrdperf/3
surface_forestperfm(cvKnonmassgrdperf)
print(cvKnonmassgrdperf)

###################################################