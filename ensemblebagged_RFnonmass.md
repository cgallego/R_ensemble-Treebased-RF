
### code Read and partition data   

```r
setwd("Z:/Cristina/MassNonmass/Section1 - ExperimentsUpToDate/experimentsRadiologypaper-revision/Tree-based-RF/ensemble-Treebased-RF")
library("RSQLite")

rpart_inputdata <- function(subdata) {
    sqlite <- dbDriver("SQLite")
    conn <- dbConnect(sqlite, "stage1localData.db")
    
    # 2) all T1W features
    lesionsQuery <- dbGetQuery(conn, "SELECT *\n           FROM  stage1features\n           INNER JOIN lesion ON (stage1features.lesion_id = lesion.lesion_id)\n           INNER JOIN f_dynamic ON (stage1features.lesion_id = f_dynamic.lesion_id)\n           INNER JOIN f_morphology ON (stage1features.lesion_id = f_morphology.lesion_id)\n           INNER JOIN f_texture ON (stage1features.lesion_id = f_texture.lesion_id)")
    
    # prune entries and extract feature subsets corresponds to 5 entries
    # lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
    lesionfields = names(lesionsQuery)
    lesioninfo = lesionsQuery[c(1, 2, 150, 151)]
    stage1features = lesionsQuery[c(3:103, 124:127)]
    dynfeatures = lesionsQuery[c(154:187)]
    morphofeatures = lesionsQuery[c(190:208)]
    texfeatures = lesionsQuery[c(211:234)]
    
    # combine all features
    allfeatures = cbind(lesioninfo[1:272, 2:3], stage1features[1:272, 1:ncol(stage1features)], 
        dynfeatures[1:272, 1:ncol(dynfeatures)], morphofeatures[1:272, 1:ncol(morphofeatures)], 
        texfeatures[1:272, 1:ncol(texfeatures)])
    # combine all features
    
    
    if (subdata == "mass") {
        # organized the data by subdata
        M <- subset(allfeatures, lesion_label == "massB" | lesion_label == "massM")
        M$lesion_label <- ifelse(M$lesion_label == "massB", "NC", "C")
        allfeatures = M
    }
    if (subdata == "nonmass") {
        # organized the data by subdata
        N <- subset(allfeatures, lesion_label == "nonmassB" | lesion_label == 
            "nonmassM")
        N$lesion_label <- ifelse(N$lesion_label == "nonmassB", "NC", "C")
        allfeatures = N
    }
    if (subdata == "stage1") {
        # organized the data by subdata
        M <- subset(allfeatures, lesion_label == "massB" | lesion_label == "massM")
        M$lesion_label <- ifelse(M$lesion_label == "massB", "mass", "mass")
        N <- subset(allfeatures, lesion_label == "nonmassB" | lesion_label == 
            "nonmassM")
        N$lesion_label <- ifelse(N$lesion_label == "nonmassB", "nonmass", "nonmass")
        allfeatures = data.frame(rbind(M, N))
    }
    if (subdata == "oneshot") {
        # organized the data by subdata
        M <- subset(allfeatures, lesion_label == "massB" | lesion_label == "massM")
        M$lesion_label <- ifelse(M$lesion_label == "massB", "NC", "C")
        N <- subset(allfeatures, lesion_label == "nonmassB" | lesion_label == 
            "nonmassM")
        N$lesion_label <- ifelse(N$lesion_label == "nonmassB", "NC", "C")
        allfeatures = data.frame(rbind(M, N))
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
```


### code to create a cross-validation set up: 
### cvfoldk = number of cv folds typically 5 or 10
### out: particvfoldK = all cv-K ids

```r
library(MASS)
library(caret)
```

```
## Loading required package: cluster
## Loading required package: foreach
## Loading required package: lattice
## Loading required package: plyr
## Loading required package: reshape2
```

```r

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
```



### code to sample kparti from a cross-validation set up: 
### kparti = k fold to exclude
### outs: cvTrainsetD, cvTestsetD

```r
kparti_sample <- function(dat, particvfoldK, cvfoldK, kparti) {
    allparti = 1:cvfoldK
    allbutkparti = allparti[-kparti]
    cvfoldadd = c()
    for (i in 1:length(allbutkparti)) {
        kadd = allbutkparti[i]
        cvfoldadd = c(cvfoldadd, particvfoldK[[kadd]])
    }
    # partition data
    cvTrainsetD <- dat[cvfoldadd, ]
    cvTestsetD <- dat[-cvfoldadd, ]
    
    output <- list(cvTrainsetD = cvTrainsetD, cvTestsetD = cvTestsetD)
    return(output)
}
```



### code Feature selection: 
### Boruta, cvfold, 

```r
library(Boruta)
```

```
## Loading required package: randomForest
## randomForest 4.6-7
## Type rfNews() to see new features/changes/bug fixes.
```

```r
require(data.table)
```

```
## Loading required package: data.table
```

```r
require(ggplot2)
```

```
## Loading required package: ggplot2
```

```r

# function to produce correlation coefficients on pair plots
panel.cor <- function(x, y, digits = 2, cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
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
    if (p < 0.01) 
        txt2 <- paste("p= ", "<0.01", sep = "")
    text(0.5, 0.4, txt2)
}

subset_select <- function(setTrain) {
    featsel_boruta <- Boruta(lesion_label ~ ., data = setTrain[, 2:ncol(setTrain)], 
        doTrace = 2, ntree = 1000)
    print(featsel_boruta)
    plot(featsel_boruta)
    
    relevant <- featsel_boruta$finalDecision[featsel_boruta$finalDecision == 
        "Confirmed"]
    relevant_features = setTrain[c(names(relevant))]
    tentative <- featsel_boruta$finalDecision[featsel_boruta$finalDecision == 
        "Tentative"]
    tentative_features = setTrain[c(names(tentative))]
    sel_features = cbind(setTrain[c(1, 2)], relevant_features, tentative_features)
    
    super.sym <- trellis.par.get("superpose.symbol")
    ## pair plots for reatures
    setTrainrelevant = setTrain[c(names(relevant))]
    pairs(relevant_features, upper.panel = panel.cor, pch = super.sym$pch[1:2], 
        col = super.sym$col[1:2], text = list(levels(setTrainrelevant$lesion_label)), 
        main = "Relevant features")
    
    return(sel_features)
}
```



### code forest Train: 
### parameters, T= # of trees, D= tree depth, dat

```r
library(klaR)
library(rpart)
library(rpart.plot)

# bagged training was introduced as a way of reducing possible overfitting
# and improving the generalization capabilities of random forests.  The
# idea is to train each tree in a forest on a different training subset,
# sampled at random from the same labeled database.
rpart_looforestTrain <- function(T, D, dat) {
    # set control
    fitparm = rpart.control(maxdepth = D, minsplit = 5, minbucket = 4, cp = 1e-05, 
        xval = 5, maxcompete = 0, maxsurrogate = 0, usesurrogate = 0, surrogatestyle = 0)
    
    # init forest
    forest = list()
    for (t in 1:T) {
        # cat('Tree # ', t, '\n')
        
        # build bagged trees from a bootstrap sample of trainSetD
        setD = dat[sample(1:nrow(dat), nrow(dat), replace = TRUE), ]
        
        # find subsample of var when training the ith tree we only make available
        # a small random subset
        subvar = sample(2:ncol(setD), sqrt(ncol(setD) - 1), replace = FALSE)
        subfeat = colnames(setD)[subvar]
        
        # train tree
        treedata <- rpart(paste("lesion_label ~ ", paste(subfeat, collapse = "+")), 
            method = "class", data = setD, control = fitparm)
        
        # display the probability per class of observations in the node
        # (conditioned on the node, sum across a node is 1) plus the percentage of
        # observations in the node.
        if (T == 1) {
            print(treedata)
            prp(treedata, type = 2, digits = 3, extra = 102, under = TRUE, nn = TRUE, 
                col = "black", box.col = rainbow(2)[2], varlen = 0, faclen = 0, 
                branch.type = 0, gap = 0, cex = 0.7, fallen.leaves = TRUE)  # use fallen.leaves=TRUE, to plot at bottom  
        }
        
        # append
        forest <- append(forest, list(tree = treedata))
    }
    
    output <- list(forest = forest)
    return(output)
}
```



### code forest Test: 
### parameters, T= # of trees, forest, TrainsetD, TestsetD

```r

library(pROC)
```

```
## Type 'citation("pROC")' for a citation.
## 
## Attaching package: 'pROC'
## 
## The following object(s) are masked from 'package:stats':
## 
##     cov, smooth, var
```

```r
rpart_looforestTest <- function(T, TrainsetD, TestsetD, forest) {
    
    fclasspotrain = list()
    for (t in 1:T) {
        # Calcultate posterior Probabilities on grid points
        temp <- predict(forest[t]$tree, newdata = TrainsetD)  #
        fclasspotrain <- append(fclasspotrain, list(cpo = temp))
    }
    
    # run testing cases
    fclasspotest = list()
    for (t in 1:T) {
        # Calcultate posterior Probabilities on grid points
        temp <- predict(forest[t]$tree, newdata = TestsetD)  #
        fclasspotest <- append(fclasspotest, list(cpo = temp))
    }
    
    # performance on Train/Test set separately extract ensamble class
    # probabilities (when T > 1)
    trainpts = fclasspotrain[1]$cpo
    testpts = fclasspotest[1]$cpo
    # init ensample class posteriors
    enclasspotrain <- matrix(, nrow = nrow(as.data.frame(trainpts)), ncol = 2)
    enclasspotest <- matrix(, nrow = nrow(as.data.frame(testpts)), ncol = 2)
    enclasspotrain[, 1] = fclasspotrain[1]$cpo[, 1]
    enclasspotest[, 1] = fclasspotest[1]$cpo[, 1]
    enclasspotrain[, 2] = fclasspotrain[1]$cpo[, 2]
    enclasspotest[, 2] = fclasspotest[1]$cpo[, 2]
    if (T >= 2) {
        for (t in 2:T) {
            # train
            enclasspotrain[, 1] = enclasspotrain[, 1] + fclasspotrain[t]$cpo[, 
                1]
            enclasspotrain[, 2] = enclasspotrain[, 2] + fclasspotrain[t]$cpo[, 
                2]
            # test
            enclasspotest[, 1] = enclasspotest[, 1] + fclasspotest[t]$cpo[, 
                1]
            enclasspotest[, 2] = enclasspotest[, 2] + fclasspotest[t]$cpo[, 
                2]
        }
    }
    # majority voting averaging
    enclasspotrain = (1/T) * enclasspotrain
    enclasspotest = (1/T) * enclasspotest
    
    # on training
    classes = levels(TrainsetD$lesion_label)
    trainprob = data.frame(C1 = enclasspotrain[, 1], C2 = enclasspotrain[, 2], 
        pred = classes[apply(enclasspotrain, 1, which.max)], obs = TrainsetD$lesion_label)
    colnames(trainprob)[1:2] <- classes
    pred = as.factor(apply(enclasspotrain, 1, which.max))
    levels(pred) = levels(as.factor(unclass(TrainsetD$lesion_label)))
    perf_train = confusionMatrix(pred, as.factor(unclass(TrainsetD$lesion_label)))
    # print(perf_train)
    
    # on testing
    testprob = data.frame(C1 = enclasspotest[, 1], C2 = enclasspotest[, 2], 
        pred = classes[apply(enclasspotest, 1, which.max)], obs = TestsetD$lesion_label)
    colnames(testprob)[1:2] <- classes
    pred = as.factor(apply(enclasspotest, 1, which.max))
    levels(pred) = levels(as.factor(unclass(TrainsetD$lesion_label)))
    pred[1] = as.factor(apply(enclasspotest, 1, which.max))
    
    groundT = as.factor(unclass(TestsetD$lesion_label))
    levels(groundT) = levels(as.factor(unclass(TrainsetD$lesion_label)))
    groundT[1] = as.factor(unclass(TestsetD$lesion_label))
    
    perf_test = confusionMatrix(pred, groundT)
    # print(perf_test)
    
    output <- list(etrain = perf_train, etest = perf_test, trainprob = trainprob, 
        testprob = testprob)
    return(output)
}
```



### code for running and plotting perfm results: 
###statsAU

```r
create_ensemble <- function(dat, particvfoldK, cvK) {
    # inint
    ensemblegrdperf = list()
    maxM = list()
    for (r in 1:cvK) {
        ## pick one of cvfold for held-out test, train on the rest
        kparti_setdata = kparti_sample(dat, particvfoldK, cvK, r)
        
        # Boruta on $cvTrainsetD
        selfeatures_kfold = subset_select(kparti_setdata$cvTrainsetD)
        names(selfeatures_kfold)
        
        ################################################### create grid of
        ################################################### evaluation points
        gT = c(5, 10, 30, 60, 100, 250, 500, 750)
        gD = c(2, 5, 10, 20)
        grd <- expand.grid(x = gD, y = gT)
        
        ################################################### for oneshot
        grdperf = data.frame(grd)
        grdperf$acuTrain = 0
        grdperf$rocTrain = 0
        grdperf$senTrain = 0
        grdperf$speTrain = 0
        
        grdperf$acuTest = 0
        grdperf$rocTest = 0
        grdperf$senTest = 0
        grdperf$speTest = 0
        
        M = list()
        for (k in 1:nrow(grd)) {
            D = grd[k, 1]
            T = grd[k, 2]
            # Build in l
            cat("D: ", D, "T: ", T, "\n")
            TrainsetD <- kparti_setdata$cvTrainsetD[c(names(selfeatures_kfold))]
            TestsetD <- kparti_setdata$cvTestsetD[c(names(selfeatures_kfold))]
            fit <- rpart_looforestTrain(T, D, TrainsetD[c(2:ncol(TrainsetD))])
            # # predict
            perf <- rpart_looforestTest(T, TrainsetD[c(2:ncol(TrainsetD))], 
                TestsetD[c(2:ncol(TestsetD))], fit$forest)
            # for train
            ROCF_train <- roc(perf$trainprob$obs, perf$trainprob$C, col = "#000086", 
                main = paste0("ROC T=", T, " D=", D, " cv=", r))
            print(ROCF_train$auc)
            # collect data
            grdperf$acuTrain[k] = grdperf$acuTrain[k] + as.numeric(perf$etrain$overall[1])
            grdperf$rocTrain[k] = grdperf$rocTrain[k] + as.numeric(ROCF_train$auc)
            grdperf$senTrain[k] = grdperf$senTrain[k] + as.numeric(perf$etrain$byClass[1])
            grdperf$speTrain[k] = grdperf$speTrain[k] + as.numeric(perf$etrain$byClass[2])
            # for test par(new=TRUE)
            ROCF_test <- roc(perf$testprob$obs, perf$testprob$C, col = "#860000", 
                main = paste0("ROC T=", T, " D=", D, " cv=", r))
            # legend('bottomright', legend = c('train', 'test'), col = c('#000086',
            # '#860000'),lwd = 2)
            print(ROCF_test$auc)
            # collect data
            grdperf$acuTest[k] = grdperf$acuTest[k] + as.numeric(perf$etest$overall[1])
            grdperf$rocTest[k] = grdperf$rocTest[k] + as.numeric(ROCF_test$auc)
            grdperf$senTest[k] = grdperf$senTest[k] + as.numeric(perf$etest$byClass[1])
            grdperf$speTest[k] = grdperf$speTest[k] + as.numeric(perf$etest$byClass[2])
            
            # append perfm for ROC
            M = append(M, list(M = list(D = D, T = T, trainprob = perf$trainprob, 
                testprob = perf$testprob, forest = fit$forest)))
        }
        print(grdperf)
        surface_forestperfm(grdperf)
        index = which(grdperf$rocTest == max(grdperf$rocTest), arr.ind = TRUE)
        Dmax = grdperf$x[index]
        Tmax = grdperf$y[index]
        resamMax = M[index]$M$testprob
        # append
        maxM <- append(maxM, list(maxp = list(D = Dmax, T = Tmax, trainprob = M[index]$M$trainprob, 
            testprob = M[index]$M$testprob, forest = M[index]$M$forest)))
        ensemblegrdperf <- append(ensemblegrdperf, list(grdperf = grdperf))
    }
    output <- list(ensemblegrdperf = ensemblegrdperf, maxM = maxM)
    return(output)
}


surface_forestperfm <- function(grdperf) {
    library(gridExtra)
    library(base)
    library(lattice)
    
    
    graphlist <- list()
    count <- 1
    # design acuTrain
    z = grdperf$acuTrain
    gD = unique(grdperf$x)
    gT = unique(grdperf$y)
    dim(z) <- c(length(gD), length(gT))
    w1 <- wireframe(z, gD, gT, box = FALSE, xlab = "Depth of trees (D)", ylab = "Number of trees (T)", 
        main = "Influence of forest parameters on Accuracy Train", drape = TRUE, 
        colorkey = TRUE, light.source = c(10, 0, 10), col.regions = colorRampPalette(c("red", 
            "blue"))(100), screen = list(z = 30, x = -60))
    graphlist[[count]] <- w1
    count <- count + 1
    
    # design rocTrain
    z = grdperf$rocTrain
    dim(z) <- c(length(gD), length(gT))
    w2 <- wireframe(z, gD, gT, box = FALSE, xlab = "Depth of trees (D)", ylab = "Number of trees (T)", 
        main = "Influence of forest parameters on ROC Train", drape = TRUE, 
        colorkey = TRUE, light.source = c(10, 0, 10), col.regions = colorRampPalette(c("red", 
            "blue"))(100), screen = list(z = 30, x = -60))
    graphlist[[count]] <- w2
    count <- count + 1
    
    # design acuTest
    z = grdperf$acuTest
    dim(z) <- c(length(gD), length(gT))
    w3 <- wireframe(z, gD, gT, box = FALSE, xlab = "Depth of trees (D)", ylab = "Number of trees (T)", 
        main = "Influence of forest parameters on Accuracy Test", drape = TRUE, 
        colorkey = TRUE, light.source = c(10, 0, 10), col.regions = colorRampPalette(c("red", 
            "blue"))(100), screen = list(z = 30, x = -60))
    graphlist[[count]] <- w3
    count <- count + 1
    
    # design rocTest
    z = grdperf$rocTest
    dim(z) <- c(length(gD), length(gT))
    w4 <- wireframe(z, gD, gT, box = FALSE, xlab = "Depth of trees (D)", ylab = "Number of trees (T)", 
        main = "Influence of forest parameters on ROC Test", drape = TRUE, colorkey = TRUE, 
        light.source = c(10, 0, 10), col.regions = colorRampPalette(c("red", 
            "blue"))(100), screen = list(z = 30, x = -60))
    graphlist[[count]] <- w4
    count <- count + 1
    
    
    # finally plot in grid
    do.call("grid.arrange", c(graphlist, ncol = 2))
}
```


Run for nonmass lesions:
=====

```r

# read nonmass features
nonmassdat = rpart_inputdata(subdata = "nonmass")
## create CV
cvK = 10
# run for nonmass
particvfoldK = cvfold_partition(nonmassdat, cvK)
res = create_ensemble(nonmassdat, particvfoldK, cvK)
```

```
## Initial round 1: ..........
##  107  attributes rejected after this test:  V0 V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V12 V13 V14 V15 V16 V17 V18 V19 earlySE2 earlySE3 earlySE4 earlySE10 earlySE12 earlySE14 dce2SE3 dce2SE8 dce2SE10 dce2SE11 dce2SE12 dce2SE14 dce2SE18 dce3SE3 dce3SE7 dce3SE10 dce3SE11 dce3SE17 dce3SE18 lateSE4 lateSE5 lateSE10 lateSE11 lateSE12 lateSE13 lateSE15 lateSE16 lateSE19 degreeC closenessC betweennessC no_triangles no_con_comp A_inside alpha_inside beta_inside iAUC1_inside Tpeak_inside Kpeak_inside peakCr_inside washoutRate_inside maxVr_inside peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside Vr_post_1_inside Tpeak_countor Kpeak_countor maxCr_countor peakCr_countor UptakeRate_countor washoutRate_countor maxVr_countor peakVr_countor Vr_increasingRate_countor Vr_post_1_countor min_F_r_i skew_F_r_i kurt_F_r_i k_Max_Margin_Grad ivVariance circularity irregularity edge_sharp_mean max_RGH_mean max_RGH_mean_k max_RGH_var max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad texture_contrast_halfRad texture_contrast_threeQuaRad texture_homogeneity_zero texture_homogeneity_quarterRad texture_homogeneity_threeQuaRad texture_dissimilarity_zero texture_dissimilarity_quarterRad texture_dissimilarity_halfRad texture_dissimilarity_threeQuaRad texture_correlation_zero texture_correlation_quarterRad texture_correlation_threeQuaRad texture_ASM_zero texture_ASM_quarterRad texture_ASM_halfRad texture_energy_zero texture_energy_quarterRad texture_energy_halfRad 
## 
## Initial round 2: ..........
##  11  attributes rejected after this test:  V11 earlySE13 dce2SE4 lateSE0 lateSE7 lateSE18 iMax_Variance_uptake iiMin_change_Variance_uptake iiiMax_Margin_Gradient texture_homogeneity_halfRad texture_energy_threeQuaRad 
## 
## Initial round 3: ..........
##  15  attributes rejected after this test:  earlySE0 earlySE6 earlySE11 earlySE18 dce2SE0 dce2SE2 dce2SE6 dce2SE17 dce3SE9 dce3SE12 dce3SE16 SER_inside Vr_decreasingRate_countor texture_correlation_halfRad texture_ASM_threeQuaRad 
## 
## Final round: ..........
##  18  attributes rejected after this test:  earlySE9 earlySE19 dce2SE1 dce2SE7 dce2SE9 dce2SE13 dce2SE19 dce3SE4 dce3SE14 lateSE1 lateSE3 lateSE8 lateSE14 lateSE17 maxCr_inside beta_countor iAUC1_countor mean_F_r_i 
## ....
##  2  attributes confirmed after this test:  Slope_ini_countor edge_sharp_std 
## 
##  3  attributes rejected after this test:  earlySE7 dce3SE8 UptakeRate_inside 
## ....
##  1  attributes confirmed after this test:  earlySE17 
## ...
##  1  attributes rejected after this test:  dce3SE19 
## ...
##  1  attributes rejected after this test:  dce3SE1 
## ...
##  1  attributes rejected after this test:  dce3SE5 
## ...
##  1  attributes confirmed after this test:  max_F_r_i 
## 
##  2  attributes rejected after this test:  earlySE16 dce3SE15 
## ...
##  1  attributes confirmed after this test:  earlySE15 
## ..
##  1  attributes confirmed after this test:  earlySE5 
## ......
##  1  attributes confirmed after this test:  alpha_countor 
## ..
##  2  attributes confirmed after this test:  dce2SE16 dce3SE13 
## ...
##  1  attributes confirmed after this test:  dce2SE15 
## ..........
##  1  attributes rejected after this test:  earlySE8 
## .........................
##  1  attributes confirmed after this test:  dce2SE5 
## ............
##  1  attributes confirmed after this test:  dce3SE0 
## .......
## Boruta performed 130 randomForest runs in 43.33 secs.
##         12 attributes confirmed important: earlySE5 earlySE15
## earlySE17 dce2SE5 dce2SE15 dce2SE16 dce3SE0 dce3SE13 alpha_countor
## Slope_ini_countor max_F_r_i edge_sharp_std
##         160 attributes confirmed unimportant: V0 V1 V2 V3 V4 V5 V6
## V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE2
## earlySE3 earlySE4 earlySE6 earlySE7 earlySE8 earlySE9 earlySE10
## earlySE11 earlySE12 earlySE13 earlySE14 earlySE16 earlySE18
## earlySE19 dce2SE0 dce2SE1 dce2SE2 dce2SE3 dce2SE4 dce2SE6 dce2SE7
## dce2SE8 dce2SE9 dce2SE10 dce2SE11 dce2SE12 dce2SE13 dce2SE14
## dce2SE17 dce2SE18 dce2SE19 dce3SE1 dce3SE3 dce3SE4 dce3SE5 dce3SE7
## dce3SE8 dce3SE9 dce3SE10 dce3SE11 dce3SE12 dce3SE14 dce3SE15
## dce3SE16 dce3SE17 dce3SE18 dce3SE19 lateSE0 lateSE1 lateSE3
## lateSE4 lateSE5 lateSE7 lateSE8 lateSE10 lateSE11 lateSE12
## lateSE13 lateSE14 lateSE15 lateSE16 lateSE17 lateSE18 lateSE19
## degreeC closenessC betweennessC no_triangles no_con_comp A_inside
## alpha_inside beta_inside iAUC1_inside Tpeak_inside Kpeak_inside
## SER_inside maxCr_inside peakCr_inside UptakeRate_inside
## washoutRate_inside maxVr_inside peakVr_inside
## Vr_increasingRate_inside Vr_decreasingRate_inside Vr_post_1_inside
## beta_countor iAUC1_countor Tpeak_countor Kpeak_countor
## maxCr_countor peakCr_countor UptakeRate_countor
## washoutRate_countor maxVr_countor peakVr_countor
## Vr_increasingRate_countor Vr_decreasingRate_countor
## Vr_post_1_countor min_F_r_i mean_F_r_i skew_F_r_i kurt_F_r_i
## iMax_Variance_uptake iiMin_change_Variance_uptake
## iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance circularity
## irregularity edge_sharp_mean max_RGH_mean max_RGH_mean_k
## max_RGH_var max_RGH_var_k texture_contrast_zero
## texture_contrast_quarterRad texture_contrast_halfRad
## texture_contrast_threeQuaRad texture_homogeneity_zero
## texture_homogeneity_quarterRad texture_homogeneity_halfRad
## texture_homogeneity_threeQuaRad texture_dissimilarity_zero
## texture_dissimilarity_quarterRad texture_dissimilarity_halfRad
## texture_dissimilarity_threeQuaRad texture_correlation_zero
## texture_correlation_quarterRad texture_correlation_halfRad
## texture_correlation_threeQuaRad texture_ASM_zero
## texture_ASM_quarterRad texture_ASM_halfRad texture_ASM_threeQuaRad
## texture_energy_zero texture_energy_quarterRad
## texture_energy_halfRad texture_energy_threeQuaRad
##         10 tentative attributes left: earlySE1 dce3SE2 dce3SE6
## lateSE2 lateSE6 lateSE9 Slope_ini_inside A_countor SER_countor
## var_F_r_i
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-81.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-82.png) 

```
## D:  2 T:  5 
## Area under the curve: 0.945
## Area under the curve: 0.567
## D:  5 T:  5 
## Area under the curve: 0.984
## Area under the curve: 0.867
## D:  10 T:  5 
## Area under the curve: 0.956
## Area under the curve: 0.4
## D:  20 T:  5 
## Area under the curve: 0.907
## Area under the curve: 0.633
## D:  2 T:  10 
## Area under the curve: 0.91
## Area under the curve: 0.467
## D:  5 T:  10 
## Area under the curve: 0.961
## Area under the curve: 0.533
## D:  10 T:  10 
## Area under the curve: 0.983
## Area under the curve: 0.8
## D:  20 T:  10 
## Area under the curve: 0.954
## Area under the curve: 0.6
## D:  2 T:  30 
## Area under the curve: 0.926
## Area under the curve: 0.467
## D:  5 T:  30 
## Area under the curve: 0.972
## Area under the curve: 0.733
## D:  10 T:  30 
## Area under the curve: 0.978
## Area under the curve: 0.6
## D:  20 T:  30 
## Area under the curve: 0.973
## Area under the curve: 0.6
## D:  2 T:  60 
## Area under the curve: 0.94
## Area under the curve: 0.667
## D:  5 T:  60 
## Area under the curve: 0.976
## Area under the curve: 0.533
## D:  10 T:  60 
## Area under the curve: 0.972
## Area under the curve: 0.667
## D:  20 T:  60 
## Area under the curve: 0.973
## Area under the curve: 0.6
## D:  2 T:  100 
## Area under the curve: 0.944
## Area under the curve: 0.733
## D:  5 T:  100 
## Area under the curve: 0.981
## Area under the curve: 0.6
## D:  10 T:  100 
## Area under the curve: 0.984
## Area under the curve: 0.733
## D:  20 T:  100 
## Area under the curve: 0.969
## Area under the curve: 0.6
## D:  2 T:  250 
## Area under the curve: 0.936
## Area under the curve: 0.667
## D:  5 T:  250 
## Area under the curve: 0.972
## Area under the curve: 0.533
## D:  10 T:  250 
## Area under the curve: 0.979
## Area under the curve: 0.533
## D:  20 T:  250 
## Area under the curve: 0.981
## Area under the curve: 0.6
## D:  2 T:  500 
## Area under the curve: 0.939
## Area under the curve: 0.667
## D:  5 T:  500 
## Area under the curve: 0.98
## Area under the curve: 0.667
## D:  10 T:  500 
## Area under the curve: 0.985
## Area under the curve: 0.6
## D:  20 T:  500 
## Area under the curve: 0.984
## Area under the curve: 0.6
## D:  2 T:  750 
## Area under the curve: 0.941
## Area under the curve: 0.667
## D:  5 T:  750 
## Area under the curve: 0.981
## Area under the curve: 0.6
## D:  10 T:  750 
## Area under the curve: 0.981
## Area under the curve: 0.533
## D:  20 T:  750 
## Area under the curve: 0.982
## Area under the curve: 0.733
##     x   y acuTrain rocTrain senTrain speTrain acuTest rocTest senTest
## 1   2   5   0.8289   0.9449   0.9583   0.6071   0.500  0.5667     0.8
## 2   5   5   0.9079   0.9836   0.9167   0.8929   0.625  0.8667     1.0
## 3  10   5   0.8553   0.9557   0.8542   0.8571   0.500  0.4000     0.8
## 4  20   5   0.8947   0.9066   0.9375   0.8214   0.500  0.6333     0.8
## 5   2  10   0.8684   0.9103   0.9375   0.7500   0.500  0.4667     0.8
## 6   5  10   0.9079   0.9606   1.0000   0.7500   0.500  0.5333     0.8
## 7  10  10   0.8947   0.9829   0.9375   0.8214   0.500  0.8000     0.8
## 8  20  10   0.8684   0.9539   0.9583   0.7143   0.500  0.6000     0.8
## 9   2  30   0.8553   0.9263   0.9583   0.6786   0.500  0.4667     0.8
## 10  5  30   0.9211   0.9725   0.9792   0.8214   0.500  0.7333     0.8
## 11 10  30   0.8947   0.9777   0.9583   0.7857   0.500  0.6000     0.8
## 12 20  30   0.9342   0.9732   1.0000   0.8214   0.500  0.6000     0.8
## 13  2  60   0.9079   0.9397   1.0000   0.7500   0.500  0.6667     0.8
## 14  5  60   0.9079   0.9762   0.9792   0.7857   0.500  0.5333     0.8
## 15 10  60   0.9342   0.9725   1.0000   0.8214   0.500  0.6667     0.8
## 16 20  60   0.9211   0.9732   1.0000   0.7857   0.500  0.6000     0.8
## 17  2 100   0.8947   0.9442   0.9792   0.7500   0.500  0.7333     0.8
## 18  5 100   0.9211   0.9814   0.9792   0.8214   0.500  0.6000     0.8
## 19 10 100   0.9474   0.9844   1.0000   0.8571   0.500  0.7333     0.8
## 20 20 100   0.9079   0.9695   0.9792   0.7857   0.500  0.6000     0.8
## 21  2 250   0.8947   0.9360   0.9792   0.7500   0.500  0.6667     0.8
## 22  5 250   0.9342   0.9725   1.0000   0.8214   0.500  0.5333     0.8
## 23 10 250   0.9342   0.9792   1.0000   0.8214   0.500  0.5333     0.8
## 24 20 250   0.9342   0.9814   1.0000   0.8214   0.500  0.6000     0.8
## 25  2 500   0.8947   0.9390   0.9792   0.7500   0.500  0.6667     0.8
## 26  5 500   0.9342   0.9799   1.0000   0.8214   0.500  0.6667     0.8
## 27 10 500   0.9342   0.9851   1.0000   0.8214   0.500  0.6000     0.8
## 28 20 500   0.9342   0.9836   1.0000   0.8214   0.500  0.6000     0.8
## 29  2 750   0.8947   0.9412   0.9792   0.7500   0.500  0.6667     0.8
## 30  5 750   0.9342   0.9807   1.0000   0.8214   0.500  0.6000     0.8
## 31 10 750   0.9605   0.9807   1.0000   0.8929   0.500  0.5333     0.8
## 32 20 750   0.9342   0.9821   1.0000   0.8214   0.500  0.7333     0.8
##    speTest
## 1        0
## 2        0
## 3        0
## 4        0
## 5        0
## 6        0
## 7        0
## 8        0
## 9        0
## 10       0
## 11       0
## 12       0
## 13       0
## 14       0
## 15       0
## 16       0
## 17       0
## 18       0
## 19       0
## 20       0
## 21       0
## 22       0
## 23       0
## 24       0
## 25       0
## 26       0
## 27       0
## 28       0
## 29       0
## 30       0
## 31       0
## 32       0
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-83.png) 

```
## Initial round 1: ..........
##  122  attributes rejected after this test:  V0 V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE3 earlySE4 earlySE6 earlySE10 earlySE12 earlySE13 dce2SE1 dce2SE2 dce2SE3 dce2SE4 dce2SE5 dce2SE7 dce2SE9 dce2SE10 dce2SE12 dce2SE14 dce2SE15 dce2SE18 dce3SE3 dce3SE5 dce3SE7 dce3SE9 dce3SE10 dce3SE14 dce3SE16 dce3SE18 lateSE0 lateSE5 lateSE6 lateSE7 lateSE10 lateSE11 lateSE13 lateSE14 lateSE15 lateSE16 lateSE17 lateSE19 degreeC closenessC betweennessC no_triangles no_con_comp A_inside alpha_inside beta_inside iAUC1_inside Tpeak_inside Kpeak_inside SER_inside peakCr_inside washoutRate_inside maxVr_inside peakVr_inside Vr_increasingRate_inside Vr_post_1_inside A_countor iAUC1_countor Kpeak_countor peakCr_countor UptakeRate_countor washoutRate_countor maxVr_countor peakVr_countor Vr_increasingRate_countor Vr_decreasingRate_countor Vr_post_1_countor min_F_r_i kurt_F_r_i iiMin_change_Variance_uptake iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance circularity irregularity edge_sharp_mean max_RGH_mean max_RGH_mean_k max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad texture_contrast_halfRad texture_contrast_threeQuaRad texture_homogeneity_quarterRad texture_homogeneity_halfRad texture_homogeneity_threeQuaRad texture_dissimilarity_zero texture_dissimilarity_quarterRad texture_dissimilarity_halfRad texture_correlation_zero texture_correlation_quarterRad texture_correlation_halfRad texture_correlation_threeQuaRad texture_ASM_zero texture_ASM_quarterRad texture_ASM_halfRad texture_ASM_threeQuaRad texture_energy_zero texture_energy_quarterRad texture_energy_halfRad texture_energy_threeQuaRad 
## 
## Initial round 2: ..........
##  14  attributes rejected after this test:  earlySE14 dce2SE6 dce2SE8 dce2SE11 dce2SE17 dce3SE11 dce3SE15 lateSE4 lateSE18 Vr_decreasingRate_inside skew_F_r_i iMax_Variance_uptake texture_homogeneity_zero texture_dissimilarity_threeQuaRad 
## 
## Initial round 3: ..........
##  10  attributes rejected after this test:  earlySE2 earlySE7 earlySE9 dce3SE1 dce3SE12 lateSE1 lateSE3 lateSE12 Tpeak_countor maxCr_countor 
## 
## Final round: ..........
##  14  attributes rejected after this test:  earlySE11 earlySE16 earlySE18 dce2SE13 dce2SE16 dce2SE19 dce3SE8 dce3SE17 dce3SE19 lateSE2 lateSE8 SER_countor var_F_r_i max_RGH_var 
## ....
##  2  attributes rejected after this test:  dce2SE0 dce3SE2 
## ....
##  1  attributes confirmed after this test:  Slope_ini_countor 
## .........
##  1  attributes confirmed after this test:  max_F_r_i 
## ........
##  1  attributes rejected after this test:  beta_countor 
## ...
##  2  attributes confirmed after this test:  earlySE17 edge_sharp_std 
## ................
##  1  attributes confirmed after this test:  alpha_countor 
## .....
##  2  attributes confirmed after this test:  dce3SE4 maxCr_inside 
## .................
##  2  attributes confirmed after this test:  earlySE19 dce3SE13 
## ........................
##  1  attributes confirmed after this test:  earlySE8 
## 
## Boruta performed 130 randomForest runs in 45.08 secs.
##         10 attributes confirmed important: earlySE8 earlySE17
## earlySE19 dce3SE4 dce3SE13 maxCr_inside alpha_countor
## Slope_ini_countor max_F_r_i edge_sharp_std
##         163 attributes confirmed unimportant: V0 V1 V2 V3 V4 V5 V6
## V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE2
## earlySE3 earlySE4 earlySE6 earlySE7 earlySE9 earlySE10 earlySE11
## earlySE12 earlySE13 earlySE14 earlySE16 earlySE18 dce2SE0 dce2SE1
## dce2SE2 dce2SE3 dce2SE4 dce2SE5 dce2SE6 dce2SE7 dce2SE8 dce2SE9
## dce2SE10 dce2SE11 dce2SE12 dce2SE13 dce2SE14 dce2SE15 dce2SE16
## dce2SE17 dce2SE18 dce2SE19 dce3SE1 dce3SE2 dce3SE3 dce3SE5 dce3SE7
## dce3SE8 dce3SE9 dce3SE10 dce3SE11 dce3SE12 dce3SE14 dce3SE15
## dce3SE16 dce3SE17 dce3SE18 dce3SE19 lateSE0 lateSE1 lateSE2
## lateSE3 lateSE4 lateSE5 lateSE6 lateSE7 lateSE8 lateSE10 lateSE11
## lateSE12 lateSE13 lateSE14 lateSE15 lateSE16 lateSE17 lateSE18
## lateSE19 degreeC closenessC betweennessC no_triangles no_con_comp
## A_inside alpha_inside beta_inside iAUC1_inside Tpeak_inside
## Kpeak_inside SER_inside peakCr_inside washoutRate_inside
## maxVr_inside peakVr_inside Vr_increasingRate_inside
## Vr_decreasingRate_inside Vr_post_1_inside A_countor beta_countor
## iAUC1_countor Tpeak_countor Kpeak_countor SER_countor
## maxCr_countor peakCr_countor UptakeRate_countor
## washoutRate_countor maxVr_countor peakVr_countor
## Vr_increasingRate_countor Vr_decreasingRate_countor
## Vr_post_1_countor min_F_r_i var_F_r_i skew_F_r_i kurt_F_r_i
## iMax_Variance_uptake iiMin_change_Variance_uptake
## iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance circularity
## irregularity edge_sharp_mean max_RGH_mean max_RGH_mean_k
## max_RGH_var max_RGH_var_k texture_contrast_zero
## texture_contrast_quarterRad texture_contrast_halfRad
## texture_contrast_threeQuaRad texture_homogeneity_zero
## texture_homogeneity_quarterRad texture_homogeneity_halfRad
## texture_homogeneity_threeQuaRad texture_dissimilarity_zero
## texture_dissimilarity_quarterRad texture_dissimilarity_halfRad
## texture_dissimilarity_threeQuaRad texture_correlation_zero
## texture_correlation_quarterRad texture_correlation_halfRad
## texture_correlation_threeQuaRad texture_ASM_zero
## texture_ASM_quarterRad texture_ASM_halfRad texture_ASM_threeQuaRad
## texture_energy_zero texture_energy_quarterRad
## texture_energy_halfRad texture_energy_threeQuaRad
##         9 tentative attributes left: earlySE1 earlySE5 earlySE15
## dce3SE0 dce3SE6 lateSE9 Slope_ini_inside UptakeRate_inside
## mean_F_r_i
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-84.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-85.png) 

```
## D:  2 T:  5 
## Area under the curve: 0.844
## Area under the curve: 0.4
## D:  5 T:  5 
## Area under the curve: 0.939
## Area under the curve: 0.9
## D:  10 T:  5 
## Area under the curve: 0.927
## Area under the curve: 0.567
## D:  20 T:  5 
## Area under the curve: 0.973
## Area under the curve: 0.733
## D:  2 T:  10 
## Area under the curve: 0.91
## Area under the curve: 0.8
## D:  5 T:  10 
## Area under the curve: 0.967
## Area under the curve: 0.867
## D:  10 T:  10 
## Area under the curve: 0.965
## Area under the curve: 0.767
## D:  20 T:  10 
## Area under the curve: 0.979
## Area under the curve: 0.8
## D:  2 T:  30 
## Area under the curve: 0.91
## Area under the curve: 0.867
## D:  5 T:  30 
## Area under the curve: 0.983
## Area under the curve: 0.733
## D:  10 T:  30 
## Area under the curve: 0.984
## Area under the curve: 0.8
## D:  20 T:  30 
## Area under the curve: 0.985
## Area under the curve: 0.8
## D:  2 T:  60 
## Area under the curve: 0.913
## Area under the curve: 0.8
## D:  5 T:  60 
## Area under the curve: 0.986
## Area under the curve: 0.8
## D:  10 T:  60 
## Area under the curve: 0.982
## Area under the curve: 0.733
## D:  20 T:  60 
## Area under the curve: 0.988
## Area under the curve: 0.867
## D:  2 T:  100 
## Area under the curve: 0.931
## Area under the curve: 0.8
## D:  5 T:  100 
## Area under the curve: 0.985
## Area under the curve: 0.8
## D:  10 T:  100 
## Area under the curve: 0.987
## Area under the curve: 0.733
## D:  20 T:  100 
## Area under the curve: 0.992
## Area under the curve: 0.867
## D:  2 T:  250 
## Area under the curve: 0.921
## Area under the curve: 0.867
## D:  5 T:  250 
## Area under the curve: 0.99
## Area under the curve: 0.8
## D:  10 T:  250 
## Area under the curve: 0.992
## Area under the curve: 0.867
## D:  20 T:  250 
## Area under the curve: 0.987
## Area under the curve: 0.8
## D:  2 T:  500 
## Area under the curve: 0.928
## Area under the curve: 0.867
## D:  5 T:  500 
## Area under the curve: 0.99
## Area under the curve: 0.867
## D:  10 T:  500 
## Area under the curve: 0.99
## Area under the curve: 0.867
## D:  20 T:  500 
## Area under the curve: 0.991
## Area under the curve: 0.867
## D:  2 T:  750 
## Area under the curve: 0.92
## Area under the curve: 0.8
## D:  5 T:  750 
## Area under the curve: 0.988
## Area under the curve: 0.867
## D:  10 T:  750 
## Area under the curve: 0.99
## Area under the curve: 0.867
## D:  20 T:  750 
## Area under the curve: 0.99
## Area under the curve: 0.867
##     x   y acuTrain rocTrain senTrain speTrain acuTest rocTest senTest
## 1   2   5   0.8026   0.8445   0.9167   0.6071   0.750  0.4000       1
## 2   5   5   0.8553   0.9386   0.9375   0.7143   0.750  0.9000       1
## 3  10   5   0.8421   0.9271   0.8750   0.7857   0.750  0.5667       1
## 4  20   5   0.8947   0.9732   0.9583   0.7857   0.875  0.7333       1
## 5   2  10   0.8684   0.9103   0.9375   0.7500   0.750  0.8000       1
## 6   5  10   0.9079   0.9665   0.9792   0.7857   0.875  0.8667       1
## 7  10  10   0.9079   0.9650   0.9792   0.7857   0.750  0.7667       1
## 8  20  10   0.9342   0.9792   0.9792   0.8571   0.750  0.8000       1
## 9   2  30   0.8421   0.9100   0.9583   0.6429   0.750  0.8667       1
## 10  5  30   0.9474   0.9829   1.0000   0.8571   0.750  0.7333       1
## 11 10  30   0.9342   0.9836   1.0000   0.8214   0.750  0.8000       1
## 12 20  30   0.9474   0.9851   1.0000   0.8571   0.750  0.8000       1
## 13  2  60   0.8684   0.9129   0.9583   0.7143   0.750  0.8000       1
## 14  5  60   0.9474   0.9859   1.0000   0.8571   0.750  0.8000       1
## 15 10  60   0.9605   0.9821   1.0000   0.8929   0.750  0.7333       1
## 16 20  60   0.9474   0.9881   1.0000   0.8571   0.750  0.8667       1
## 17  2 100   0.8421   0.9308   0.9375   0.6786   0.750  0.8000       1
## 18  5 100   0.9474   0.9851   1.0000   0.8571   0.750  0.8000       1
## 19 10 100   0.9474   0.9866   1.0000   0.8571   0.750  0.7333       1
## 20 20 100   0.9342   0.9918   1.0000   0.8214   0.750  0.8667       1
## 21  2 250   0.8684   0.9211   0.9792   0.6786   0.750  0.8667       1
## 22  5 250   0.9737   0.9903   1.0000   0.9286   0.750  0.8000       1
## 23 10 250   0.9605   0.9918   1.0000   0.8929   0.750  0.8667       1
## 24 20 250   0.9342   0.9874   1.0000   0.8214   0.750  0.8000       1
## 25  2 500   0.8684   0.9278   0.9792   0.6786   0.750  0.8667       1
## 26  5 500   0.9474   0.9903   1.0000   0.8571   0.750  0.8667       1
## 27 10 500   0.9474   0.9903   1.0000   0.8571   0.750  0.8667       1
## 28 20 500   0.9605   0.9911   1.0000   0.8929   0.750  0.8667       1
## 29  2 750   0.8553   0.9204   0.9583   0.6786   0.750  0.8000       1
## 30  5 750   0.9474   0.9881   1.0000   0.8571   0.750  0.8667       1
## 31 10 750   0.9342   0.9903   1.0000   0.8214   0.750  0.8667       1
## 32 20 750   0.9474   0.9896   1.0000   0.8571   0.750  0.8667       1
##    speTest
## 1   0.3333
## 2   0.3333
## 3   0.3333
## 4   0.6667
## 5   0.3333
## 6   0.6667
## 7   0.3333
## 8   0.3333
## 9   0.3333
## 10  0.3333
## 11  0.3333
## 12  0.3333
## 13  0.3333
## 14  0.3333
## 15  0.3333
## 16  0.3333
## 17  0.3333
## 18  0.3333
## 19  0.3333
## 20  0.3333
## 21  0.3333
## 22  0.3333
## 23  0.3333
## 24  0.3333
## 25  0.3333
## 26  0.3333
## 27  0.3333
## 28  0.3333
## 29  0.3333
## 30  0.3333
## 31  0.3333
## 32  0.3333
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-86.png) 

```
## Initial round 1: ..........
##  130  attributes rejected after this test:  V0 V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE2 earlySE3 earlySE4 earlySE6 earlySE7 earlySE11 earlySE12 earlySE13 earlySE14 earlySE16 earlySE18 earlySE19 dce2SE1 dce2SE2 dce2SE3 dce2SE5 dce2SE10 dce2SE11 dce2SE12 dce2SE14 dce2SE18 dce2SE19 dce3SE1 dce3SE3 dce3SE5 dce3SE7 dce3SE10 dce3SE11 dce3SE12 dce3SE15 dce3SE17 dce3SE18 lateSE0 lateSE1 lateSE3 lateSE4 lateSE5 lateSE7 lateSE9 lateSE10 lateSE11 lateSE12 lateSE13 lateSE15 lateSE16 lateSE17 lateSE18 lateSE19 degreeC closenessC betweennessC no_triangles no_con_comp A_inside alpha_inside beta_inside iAUC1_inside Tpeak_inside Kpeak_inside peakCr_inside washoutRate_inside maxVr_inside peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside Vr_post_1_inside A_countor beta_countor Tpeak_countor Kpeak_countor maxCr_countor peakCr_countor washoutRate_countor maxVr_countor peakVr_countor Vr_increasingRate_countor Vr_decreasingRate_countor Vr_post_1_countor min_F_r_i var_F_r_i skew_F_r_i iMax_Variance_uptake iiMin_change_Variance_uptake iiiMax_Margin_Gradient k_Max_Margin_Grad circularity max_RGH_mean max_RGH_mean_k max_RGH_var max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad texture_contrast_threeQuaRad texture_homogeneity_zero texture_homogeneity_quarterRad texture_homogeneity_halfRad texture_dissimilarity_zero texture_dissimilarity_quarterRad texture_dissimilarity_threeQuaRad texture_correlation_zero texture_correlation_halfRad texture_correlation_threeQuaRad texture_ASM_zero texture_ASM_quarterRad texture_ASM_halfRad texture_ASM_threeQuaRad texture_energy_zero texture_energy_quarterRad texture_energy_threeQuaRad 
## 
## Initial round 2: ..........
##  11  attributes rejected after this test:  earlySE9 earlySE10 dce2SE9 dce2SE13 dce2SE17 dce3SE2 dce3SE14 dce3SE16 ivVariance irregularity texture_correlation_quarterRad 
## 
## Initial round 3: ..........
##  7  attributes rejected after this test:  dce2SE4 dce2SE6 dce3SE9 dce3SE19 kurt_F_r_i edge_sharp_mean texture_dissimilarity_halfRad 
## 
## Final round: ..........
##  1  attributes confirmed after this test:  edge_sharp_std 
## 
##  4  attributes rejected after this test:  lateSE2 lateSE14 SER_inside UptakeRate_countor 
## ....
##  1  attributes rejected after this test:  dce2SE7 
## ....
##  1  attributes rejected after this test:  texture_energy_halfRad 
## ...
##  3  attributes rejected after this test:  dce2SE15 iAUC1_countor SER_countor 
## ...
##  5  attributes rejected after this test:  earlySE1 dce2SE0 dce2SE8 lateSE6 texture_homogeneity_threeQuaRad 
## ...
##  2  attributes rejected after this test:  dce3SE6 texture_contrast_halfRad 
## ........
##  1  attributes confirmed after this test:  Slope_ini_countor 
## ...
##  1  attributes confirmed after this test:  Slope_ini_inside 
## ...
##  1  attributes confirmed after this test:  max_F_r_i 
## ....................
##  1  attributes confirmed after this test:  mean_F_r_i 
## ..........
##  1  attributes rejected after this test:  dce2SE16 
## .....
##  1  attributes rejected after this test:  earlySE5 
## ............
##  1  attributes rejected after this test:  lateSE8 
## ..........
##  1  attributes confirmed after this test:  earlySE8 
## ..
## Boruta performed 130 randomForest runs in 37.51 secs.
##         6 attributes confirmed important: earlySE8
## Slope_ini_inside Slope_ini_countor max_F_r_i mean_F_r_i
## edge_sharp_std
##         167 attributes confirmed unimportant: V0 V1 V2 V3 V4 V5 V6
## V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE1
## earlySE2 earlySE3 earlySE4 earlySE5 earlySE6 earlySE7 earlySE9
## earlySE10 earlySE11 earlySE12 earlySE13 earlySE14 earlySE16
## earlySE18 earlySE19 dce2SE0 dce2SE1 dce2SE2 dce2SE3 dce2SE4
## dce2SE5 dce2SE6 dce2SE7 dce2SE8 dce2SE9 dce2SE10 dce2SE11 dce2SE12
## dce2SE13 dce2SE14 dce2SE15 dce2SE16 dce2SE17 dce2SE18 dce2SE19
## dce3SE1 dce3SE2 dce3SE3 dce3SE5 dce3SE6 dce3SE7 dce3SE9 dce3SE10
## dce3SE11 dce3SE12 dce3SE14 dce3SE15 dce3SE16 dce3SE17 dce3SE18
## dce3SE19 lateSE0 lateSE1 lateSE2 lateSE3 lateSE4 lateSE5 lateSE6
## lateSE7 lateSE8 lateSE9 lateSE10 lateSE11 lateSE12 lateSE13
## lateSE14 lateSE15 lateSE16 lateSE17 lateSE18 lateSE19 degreeC
## closenessC betweennessC no_triangles no_con_comp A_inside
## alpha_inside beta_inside iAUC1_inside Tpeak_inside Kpeak_inside
## SER_inside peakCr_inside washoutRate_inside maxVr_inside
## peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside
## Vr_post_1_inside A_countor beta_countor iAUC1_countor
## Tpeak_countor Kpeak_countor SER_countor maxCr_countor
## peakCr_countor UptakeRate_countor washoutRate_countor
## maxVr_countor peakVr_countor Vr_increasingRate_countor
## Vr_decreasingRate_countor Vr_post_1_countor min_F_r_i var_F_r_i
## skew_F_r_i kurt_F_r_i iMax_Variance_uptake
## iiMin_change_Variance_uptake iiiMax_Margin_Gradient
## k_Max_Margin_Grad ivVariance circularity irregularity
## edge_sharp_mean max_RGH_mean max_RGH_mean_k max_RGH_var
## max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad
## texture_contrast_halfRad texture_contrast_threeQuaRad
## texture_homogeneity_zero texture_homogeneity_quarterRad
## texture_homogeneity_halfRad texture_homogeneity_threeQuaRad
## texture_dissimilarity_zero texture_dissimilarity_quarterRad
## texture_dissimilarity_halfRad texture_dissimilarity_threeQuaRad
## texture_correlation_zero texture_correlation_quarterRad
## texture_correlation_halfRad texture_correlation_threeQuaRad
## texture_ASM_zero texture_ASM_quarterRad texture_ASM_halfRad
## texture_ASM_threeQuaRad texture_energy_zero
## texture_energy_quarterRad texture_energy_halfRad
## texture_energy_threeQuaRad
##         9 tentative attributes left: earlySE15 earlySE17 dce3SE0
## dce3SE4 dce3SE8 dce3SE13 maxCr_inside UptakeRate_inside
## alpha_countor
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-87.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-88.png) 

```
## D:  2 T:  5 
## Area under the curve: 0.91
## Area under the curve: 0.933
## D:  5 T:  5 
## Area under the curve: 0.93
## Area under the curve: 0.933
## D:  10 T:  5 
## Area under the curve: 0.965
## Area under the curve: 0.867
## D:  20 T:  5 
## Area under the curve: 0.916
## Area under the curve: 0.867
## D:  2 T:  10 
## Area under the curve: 0.941
## Area under the curve: 0.867
## D:  5 T:  10 
## Area under the curve: 0.978
## Area under the curve: 0.933
## D:  10 T:  10 
## Area under the curve: 0.981
## Area under the curve: 0.867
## D:  20 T:  10 
## Area under the curve: 0.987
## Area under the curve: 1
## D:  2 T:  30 
## Area under the curve: 0.913
## Area under the curve: 0.933
## D:  5 T:  30 
## Area under the curve: 0.984
## Area under the curve: 1
## D:  10 T:  30 
## Area under the curve: 0.984
## Area under the curve: 0.933
## D:  20 T:  30 
## Area under the curve: 0.978
## Area under the curve: 0.933
## D:  2 T:  60 
## Area under the curve: 0.921
## Area under the curve: 0.933
## D:  5 T:  60 
## Area under the curve: 0.988
## Area under the curve: 0.867
## D:  10 T:  60 
## Area under the curve: 0.993
## Area under the curve: 0.867
## D:  20 T:  60 
## Area under the curve: 0.993
## Area under the curve: 1
## D:  2 T:  100 
## Area under the curve: 0.92
## Area under the curve: 0.867
## D:  5 T:  100 
## Area under the curve: 0.988
## Area under the curve: 0.933
## D:  10 T:  100 
## Area under the curve: 0.997
## Area under the curve: 0.933
## D:  20 T:  100 
## Area under the curve: 0.987
## Area under the curve: 1
## D:  2 T:  250 
## Area under the curve: 0.931
## Area under the curve: 0.933
## D:  5 T:  250 
## Area under the curve: 0.987
## Area under the curve: 1
## D:  10 T:  250 
## Area under the curve: 0.996
## Area under the curve: 0.867
## D:  20 T:  250 
## Area under the curve: 0.996
## Area under the curve: 0.867
## D:  2 T:  500 
## Area under the curve: 0.929
## Area under the curve: 0.867
## D:  5 T:  500 
## Area under the curve: 0.996
## Area under the curve: 0.933
## D:  10 T:  500 
## Area under the curve: 0.993
## Area under the curve: 0.933
## D:  20 T:  500 
## Area under the curve: 0.994
## Area under the curve: 0.867
## D:  2 T:  750 
## Area under the curve: 0.923
## Area under the curve: 0.867
## D:  5 T:  750 
## Area under the curve: 0.994
## Area under the curve: 1
## D:  10 T:  750 
## Area under the curve: 0.994
## Area under the curve: 0.933
## D:  20 T:  750 
## Area under the curve: 0.996
## Area under the curve: 0.867
##     x   y acuTrain rocTrain senTrain speTrain acuTest rocTest senTest
## 1   2   5   0.7895   0.9100   0.9167   0.5714   0.875  0.9333     1.0
## 2   5   5   0.8553   0.9301   0.9167   0.7500   0.750  0.9333     0.8
## 3  10   5   0.8684   0.9654   0.8958   0.8214   0.875  0.8667     0.8
## 4  20   5   0.8947   0.9163   0.9583   0.7857   0.875  0.8667     0.8
## 5   2  10   0.8684   0.9412   0.9167   0.7857   0.750  0.8667     0.8
## 6   5  10   0.9079   0.9777   0.9792   0.7857   0.875  0.9333     1.0
## 7  10  10   0.9474   0.9814   0.9792   0.8929   0.750  0.8667     0.8
## 8  20  10   0.9342   0.9874   0.9583   0.8929   0.875  1.0000     1.0
## 9   2  30   0.8289   0.9129   0.9583   0.6071   0.875  0.9333     1.0
## 10  5  30   0.9474   0.9844   1.0000   0.8571   0.750  1.0000     1.0
## 11 10  30   0.9474   0.9844   1.0000   0.8571   0.875  0.9333     0.8
## 12 20  30   0.8947   0.9784   1.0000   0.7143   0.875  0.9333     0.8
## 13  2  60   0.8026   0.9211   0.9792   0.5000   0.750  0.9333     0.8
## 14  5  60   0.9474   0.9881   1.0000   0.8571   0.875  0.8667     0.8
## 15 10  60   0.9605   0.9926   1.0000   0.8929   0.750  0.8667     0.8
## 16 20  60   0.9605   0.9933   1.0000   0.8929   1.000  1.0000     1.0
## 17  2 100   0.8684   0.9196   0.9792   0.6786   0.750  0.8667     0.8
## 18  5 100   0.9474   0.9881   1.0000   0.8571   0.750  0.9333     0.8
## 19 10 100   0.9474   0.9970   1.0000   0.8571   0.750  0.9333     0.8
## 20 20 100   0.9605   0.9866   1.0000   0.8929   0.875  1.0000     0.8
## 21  2 250   0.8421   0.9308   0.9792   0.6071   0.750  0.9333     0.8
## 22  5 250   0.9605   0.9874   1.0000   0.8929   1.000  1.0000     1.0
## 23 10 250   0.9605   0.9955   1.0000   0.8929   0.875  0.8667     0.8
## 24 20 250   0.9474   0.9963   1.0000   0.8571   0.875  0.8667     0.8
## 25  2 500   0.8289   0.9293   0.9792   0.5714   0.750  0.8667     0.8
## 26  5 500   0.9605   0.9955   1.0000   0.8929   0.875  0.9333     0.8
## 27 10 500   0.9605   0.9933   1.0000   0.8929   0.750  0.9333     0.8
## 28 20 500   0.9474   0.9940   1.0000   0.8571   0.875  0.8667     0.8
## 29  2 750   0.8421   0.9234   1.0000   0.5714   0.750  0.8667     0.8
## 30  5 750   0.9474   0.9940   1.0000   0.8571   0.875  1.0000     0.8
## 31 10 750   0.9605   0.9940   1.0000   0.8929   0.875  0.9333     0.8
## 32 20 750   0.9605   0.9963   1.0000   0.8929   0.875  0.8667     0.8
##    speTest
## 1   0.6667
## 2   0.6667
## 3   1.0000
## 4   1.0000
## 5   0.6667
## 6   0.6667
## 7   0.6667
## 8   0.6667
## 9   0.6667
## 10  0.3333
## 11  1.0000
## 12  1.0000
## 13  0.6667
## 14  1.0000
## 15  0.6667
## 16  1.0000
## 17  0.6667
## 18  0.6667
## 19  0.6667
## 20  1.0000
## 21  0.6667
## 22  1.0000
## 23  1.0000
## 24  1.0000
## 25  0.6667
## 26  1.0000
## 27  0.6667
## 28  1.0000
## 29  0.6667
## 30  1.0000
## 31  1.0000
## 32  1.0000
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-89.png) 

```
## Initial round 1: ..........
##  120  attributes rejected after this test:  V0 V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V12 V13 V14 V15 V16 V17 V18 V19 earlySE3 earlySE4 earlySE6 earlySE7 earlySE10 earlySE14 dce2SE7 dce2SE8 dce2SE10 dce2SE11 dce2SE12 dce2SE14 dce2SE15 dce2SE17 dce2SE18 dce3SE1 dce3SE3 dce3SE5 dce3SE7 dce3SE9 dce3SE10 dce3SE11 dce3SE12 dce3SE15 dce3SE18 lateSE0 lateSE1 lateSE3 lateSE4 lateSE5 lateSE6 lateSE7 lateSE10 lateSE11 lateSE12 lateSE13 lateSE14 lateSE15 lateSE16 lateSE17 lateSE19 degreeC closenessC betweennessC no_triangles no_con_comp A_inside beta_inside iAUC1_inside Tpeak_inside Kpeak_inside peakCr_inside maxVr_inside peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside A_countor beta_countor iAUC1_countor Tpeak_countor Kpeak_countor peakCr_countor washoutRate_countor maxVr_countor peakVr_countor Vr_increasingRate_countor Vr_post_1_countor min_F_r_i skew_F_r_i kurt_F_r_i iMax_Variance_uptake iiiMax_Margin_Gradient k_Max_Margin_Grad circularity irregularity edge_sharp_mean max_RGH_mean max_RGH_mean_k max_RGH_var max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad texture_contrast_halfRad texture_contrast_threeQuaRad texture_homogeneity_zero texture_homogeneity_quarterRad texture_homogeneity_halfRad texture_homogeneity_threeQuaRad texture_dissimilarity_zero texture_dissimilarity_quarterRad texture_dissimilarity_halfRad texture_dissimilarity_threeQuaRad texture_correlation_zero texture_correlation_quarterRad texture_correlation_halfRad texture_ASM_zero texture_ASM_quarterRad texture_ASM_halfRad texture_ASM_threeQuaRad texture_energy_halfRad texture_energy_threeQuaRad 
## 
## Initial round 2: ..........
##  22  attributes rejected after this test:  earlySE1 earlySE2 earlySE9 earlySE12 earlySE13 dce2SE1 dce2SE2 dce2SE3 dce2SE5 dce2SE6 dce3SE17 lateSE2 alpha_inside washoutRate_inside Vr_post_1_inside maxCr_countor UptakeRate_countor Vr_decreasingRate_countor ivVariance texture_correlation_threeQuaRad texture_energy_zero texture_energy_quarterRad 
## 
## Initial round 3: ..........
##  6  attributes rejected after this test:  earlySE0 earlySE18 dce2SE4 dce2SE9 dce3SE19 lateSE18 
## 
## Final round: ..........
##  3  attributes rejected after this test:  dce2SE19 SER_inside var_F_r_i 
## ....
##  1  attributes confirmed after this test:  edge_sharp_std 
## 
##  1  attributes rejected after this test:  earlySE11 
## ....
##  2  attributes rejected after this test:  dce2SE13 iiMin_change_Variance_uptake 
## ...
##  1  attributes confirmed after this test:  Slope_ini_inside 
## 
##  2  attributes rejected after this test:  earlySE5 earlySE8 
## ...
##  1  attributes confirmed after this test:  Slope_ini_countor 
## 
##  1  attributes rejected after this test:  dce3SE8 
## ...
##  3  attributes rejected after this test:  earlySE19 dce2SE0 lateSE8 
## ...
##  2  attributes rejected after this test:  earlySE15 dce3SE13 
## ...
##  2  attributes rejected after this test:  V11 dce3SE16 
## ..
##  2  attributes rejected after this test:  dce3SE6 dce3SE14 
## ........
##  1  attributes rejected after this test:  lateSE9 
## ......
##  1  attributes confirmed after this test:  UptakeRate_inside 
## ...............
##  1  attributes rejected after this test:  dce3SE2 
## ....................................
## Boruta performed 130 randomForest runs in 37.39 secs.
##         4 attributes confirmed important: Slope_ini_inside
## UptakeRate_inside Slope_ini_countor edge_sharp_std
##         168 attributes confirmed unimportant: V0 V1 V2 V3 V4 V5 V6
## V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE1
## earlySE2 earlySE3 earlySE4 earlySE5 earlySE6 earlySE7 earlySE8
## earlySE9 earlySE10 earlySE11 earlySE12 earlySE13 earlySE14
## earlySE15 earlySE18 earlySE19 dce2SE0 dce2SE1 dce2SE2 dce2SE3
## dce2SE4 dce2SE5 dce2SE6 dce2SE7 dce2SE8 dce2SE9 dce2SE10 dce2SE11
## dce2SE12 dce2SE13 dce2SE14 dce2SE15 dce2SE17 dce2SE18 dce2SE19
## dce3SE1 dce3SE2 dce3SE3 dce3SE5 dce3SE6 dce3SE7 dce3SE8 dce3SE9
## dce3SE10 dce3SE11 dce3SE12 dce3SE13 dce3SE14 dce3SE15 dce3SE16
## dce3SE17 dce3SE18 dce3SE19 lateSE0 lateSE1 lateSE2 lateSE3 lateSE4
## lateSE5 lateSE6 lateSE7 lateSE8 lateSE9 lateSE10 lateSE11 lateSE12
## lateSE13 lateSE14 lateSE15 lateSE16 lateSE17 lateSE18 lateSE19
## degreeC closenessC betweennessC no_triangles no_con_comp A_inside
## alpha_inside beta_inside iAUC1_inside Tpeak_inside Kpeak_inside
## SER_inside peakCr_inside washoutRate_inside maxVr_inside
## peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside
## Vr_post_1_inside A_countor beta_countor iAUC1_countor
## Tpeak_countor Kpeak_countor maxCr_countor peakCr_countor
## UptakeRate_countor washoutRate_countor maxVr_countor
## peakVr_countor Vr_increasingRate_countor Vr_decreasingRate_countor
## Vr_post_1_countor min_F_r_i var_F_r_i skew_F_r_i kurt_F_r_i
## iMax_Variance_uptake iiMin_change_Variance_uptake
## iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance circularity
## irregularity edge_sharp_mean max_RGH_mean max_RGH_mean_k
## max_RGH_var max_RGH_var_k texture_contrast_zero
## texture_contrast_quarterRad texture_contrast_halfRad
## texture_contrast_threeQuaRad texture_homogeneity_zero
## texture_homogeneity_quarterRad texture_homogeneity_halfRad
## texture_homogeneity_threeQuaRad texture_dissimilarity_zero
## texture_dissimilarity_quarterRad texture_dissimilarity_halfRad
## texture_dissimilarity_threeQuaRad texture_correlation_zero
## texture_correlation_quarterRad texture_correlation_halfRad
## texture_correlation_threeQuaRad texture_ASM_zero
## texture_ASM_quarterRad texture_ASM_halfRad texture_ASM_threeQuaRad
## texture_energy_zero texture_energy_quarterRad
## texture_energy_halfRad texture_energy_threeQuaRad
##         10 tentative attributes left: earlySE16 earlySE17 dce2SE16
## dce3SE0 dce3SE4 maxCr_inside alpha_countor SER_countor max_F_r_i
## mean_F_r_i
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-810.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-811.png) 

```
## D:  2 T:  5 
## Area under the curve: 0.901
## Area under the curve: 0.75
## D:  5 T:  5 
## Area under the curve: 0.881
## Area under the curve: 0.646
## D:  10 T:  5 
## Area under the curve: 0.983
## Area under the curve: 0.792
## D:  20 T:  5 
## Area under the curve: 0.953
## Area under the curve: 0.708
## D:  2 T:  10 
## Area under the curve: 0.928
## Area under the curve: 0.958
## D:  5 T:  10 
## Area under the curve: 0.976
## Area under the curve: 0.917
## D:  10 T:  10 
## Area under the curve: 0.958
## Area under the curve: 0.875
## D:  20 T:  10 
## Area under the curve: 0.991
## Area under the curve: 0.75
## D:  2 T:  30 
## Area under the curve: 0.911
## Area under the curve: 0.833
## D:  5 T:  30 
## Area under the curve: 0.976
## Area under the curve: 0.833
## D:  10 T:  30 
## Area under the curve: 0.982
## Area under the curve: 0.875
## D:  20 T:  30 
## Area under the curve: 0.991
## Area under the curve: 0.833
## D:  2 T:  60 
## Area under the curve: 0.922
## Area under the curve: 0.792
## D:  5 T:  60 
## Area under the curve: 0.987
## Area under the curve: 0.875
## D:  10 T:  60 
## Area under the curve: 0.991
## Area under the curve: 0.792
## D:  20 T:  60 
## Area under the curve: 0.989
## Area under the curve: 0.833
## D:  2 T:  100 
## Area under the curve: 0.918
## Area under the curve: 0.833
## D:  5 T:  100 
## Area under the curve: 0.989
## Area under the curve: 0.833
## D:  10 T:  100 
## Area under the curve: 0.987
## Area under the curve: 0.833
## D:  20 T:  100 
## Area under the curve: 0.987
## Area under the curve: 0.833
## D:  2 T:  250 
## Area under the curve: 0.919
## Area under the curve: 0.833
## D:  5 T:  250 
## Area under the curve: 0.989
## Area under the curve: 0.833
## D:  10 T:  250 
## Area under the curve: 0.99
## Area under the curve: 0.875
## D:  20 T:  250 
## Area under the curve: 0.989
## Area under the curve: 0.875
## D:  2 T:  500 
## Area under the curve: 0.921
## Area under the curve: 0.833
## D:  5 T:  500 
## Area under the curve: 0.988
## Area under the curve: 0.875
## D:  10 T:  500 
## Area under the curve: 0.991
## Area under the curve: 0.917
## D:  20 T:  500 
## Area under the curve: 0.987
## Area under the curve: 0.833
## D:  2 T:  750 
## Area under the curve: 0.928
## Area under the curve: 0.833
## D:  5 T:  750 
## Area under the curve: 0.989
## Area under the curve: 0.833
## D:  10 T:  750 
## Area under the curve: 0.987
## Area under the curve: 0.875
## D:  20 T:  750 
## Area under the curve: 0.99
## Area under the curve: 0.833
##     x   y acuTrain rocTrain senTrain speTrain acuTest rocTest senTest
## 1   2   5   0.8514   0.9015   0.8723   0.8148     0.5  0.7500  0.5000
## 2   5   5   0.8378   0.8814   0.9149   0.7037     0.7  0.6458  0.5000
## 3  10   5   0.8784   0.9827   1.0000   0.6667     0.8  0.7917  0.8333
## 4  20   5   0.8784   0.9527   0.9149   0.8148     0.7  0.7083  0.6667
## 5   2  10   0.8378   0.9283   0.9574   0.6296     0.9  0.9583  1.0000
## 6   5  10   0.9054   0.9764   0.9362   0.8519     0.8  0.9167  1.0000
## 7  10  10   0.8919   0.9582   0.9574   0.7778     0.8  0.8750  0.6667
## 8  20  10   0.9459   0.9905   0.9787   0.8889     0.8  0.7500  1.0000
## 9   2  30   0.8378   0.9106   0.8936   0.7407     0.7  0.8333  0.6667
## 10  5  30   0.9324   0.9764   0.9787   0.8519     0.7  0.8333  0.6667
## 11 10  30   0.9459   0.9819   0.9787   0.8889     0.9  0.8750  1.0000
## 12 20  30   0.9730   0.9913   1.0000   0.9259     0.8  0.8333  0.8333
## 13  2  60   0.8108   0.9220   0.8723   0.7037     0.9  0.7917  1.0000
## 14  5  60   0.9730   0.9866   1.0000   0.9259     0.9  0.8750  1.0000
## 15 10  60   0.9595   0.9913   0.9787   0.9259     0.7  0.7917  0.6667
## 16 20  60   0.9595   0.9890   1.0000   0.8889     0.7  0.8333  0.6667
## 17  2 100   0.8378   0.9180   0.8936   0.7407     0.9  0.8333  1.0000
## 18  5 100   0.9595   0.9890   1.0000   0.8889     0.9  0.8333  1.0000
## 19 10 100   0.9595   0.9874   0.9787   0.9259     0.7  0.8333  0.6667
## 20 20 100   0.9324   0.9874   0.9787   0.8519     0.7  0.8333  0.6667
## 21  2 250   0.8108   0.9188   0.8936   0.6667     0.8  0.8333  0.8333
## 22  5 250   0.9595   0.9890   1.0000   0.8889     0.9  0.8333  1.0000
## 23 10 250   0.9459   0.9898   1.0000   0.8519     0.7  0.8750  0.6667
## 24 20 250   0.9730   0.9890   1.0000   0.9259     0.8  0.8750  0.8333
## 25  2 500   0.8108   0.9212   0.8936   0.6667     0.8  0.8333  0.8333
## 26  5 500   0.9595   0.9882   1.0000   0.8889     0.8  0.8750  0.8333
## 27 10 500   0.9730   0.9913   1.0000   0.9259     0.7  0.9167  0.6667
## 28 20 500   0.9730   0.9874   1.0000   0.9259     0.8  0.8333  0.8333
## 29  2 750   0.8378   0.9275   0.9149   0.7037     0.9  0.8333  1.0000
## 30  5 750   0.9730   0.9890   1.0000   0.9259     0.7  0.8333  0.6667
## 31 10 750   0.9730   0.9866   1.0000   0.9259     0.9  0.8750  1.0000
## 32 20 750   0.9730   0.9898   1.0000   0.9259     0.8  0.8333  0.8333
##    speTest
## 1     0.50
## 2     1.00
## 3     0.75
## 4     0.75
## 5     0.75
## 6     0.50
## 7     1.00
## 8     0.50
## 9     0.75
## 10    0.75
## 11    0.75
## 12    0.75
## 13    0.75
## 14    0.75
## 15    0.75
## 16    0.75
## 17    0.75
## 18    0.75
## 19    0.75
## 20    0.75
## 21    0.75
## 22    0.75
## 23    0.75
## 24    0.75
## 25    0.75
## 26    0.75
## 27    0.75
## 28    0.75
## 29    0.75
## 30    0.75
## 31    0.75
## 32    0.75
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-812.png) 

```
## Initial round 1: ..........
##  137  attributes rejected after this test:  V0 V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE2 earlySE3 earlySE4 earlySE5 earlySE6 earlySE8 earlySE9 earlySE10 earlySE11 earlySE12 earlySE14 earlySE15 earlySE18 dce2SE0 dce2SE1 dce2SE2 dce2SE3 dce2SE4 dce2SE5 dce2SE6 dce2SE7 dce2SE8 dce2SE10 dce2SE11 dce2SE12 dce2SE14 dce2SE18 dce3SE3 dce3SE5 dce3SE10 dce3SE11 dce3SE12 dce3SE14 dce3SE15 dce3SE18 dce3SE19 lateSE0 lateSE1 lateSE4 lateSE5 lateSE6 lateSE7 lateSE11 lateSE12 lateSE15 lateSE17 lateSE18 degreeC closenessC betweennessC no_triangles no_con_comp A_inside alpha_inside beta_inside iAUC1_inside Tpeak_inside SER_inside peakCr_inside UptakeRate_inside washoutRate_inside maxVr_inside peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside Vr_post_1_inside A_countor beta_countor iAUC1_countor Tpeak_countor Kpeak_countor SER_countor peakCr_countor UptakeRate_countor washoutRate_countor maxVr_countor peakVr_countor Vr_increasingRate_countor Vr_post_1_countor skew_F_r_i kurt_F_r_i iMax_Variance_uptake iiMin_change_Variance_uptake iiiMax_Margin_Gradient k_Max_Margin_Grad circularity irregularity edge_sharp_mean max_RGH_mean max_RGH_mean_k max_RGH_var max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad texture_contrast_halfRad texture_contrast_threeQuaRad texture_homogeneity_zero texture_homogeneity_quarterRad texture_homogeneity_halfRad texture_homogeneity_threeQuaRad texture_dissimilarity_zero texture_dissimilarity_quarterRad texture_dissimilarity_halfRad texture_dissimilarity_threeQuaRad texture_correlation_zero texture_correlation_quarterRad texture_correlation_halfRad texture_correlation_threeQuaRad texture_ASM_zero texture_ASM_quarterRad texture_ASM_halfRad texture_ASM_threeQuaRad texture_energy_zero texture_energy_quarterRad texture_energy_halfRad texture_energy_threeQuaRad 
## 
## Initial round 2: ..........
##  5  attributes rejected after this test:  earlySE13 dce2SE19 dce3SE9 lateSE3 min_F_r_i 
## 
## Initial round 3: ..........
##  4  attributes rejected after this test:  dce2SE15 dce2SE17 lateSE10 lateSE16 
## 
## Final round: ..........
##  7  attributes rejected after this test:  earlySE16 dce3SE2 dce3SE8 dce3SE16 lateSE2 lateSE8 maxCr_countor 
## ....
##  1  attributes confirmed after this test:  edge_sharp_std 
## 
##  2  attributes rejected after this test:  lateSE14 ivVariance 
## ....
##  1  attributes rejected after this test:  dce3SE7 
## ...
##  8  attributes rejected after this test:  dce2SE9 dce3SE1 dce3SE17 lateSE13 lateSE19 Kpeak_inside Vr_decreasingRate_countor var_F_r_i 
## ...
##  1  attributes rejected after this test:  earlySE7 
## ...
##  1  attributes rejected after this test:  earlySE1 
## ...........
##  1  attributes rejected after this test:  lateSE9 
## ........
##  2  attributes confirmed after this test:  earlySE17 max_F_r_i 
## .....
##  1  attributes confirmed after this test:  dce3SE4 
## ...
##  1  attributes confirmed after this test:  dce3SE13 
## ............
##  1  attributes confirmed after this test:  mean_F_r_i 
## .....
##  1  attributes confirmed after this test:  maxCr_inside 
## .....
##  1  attributes confirmed after this test:  Slope_ini_countor 
## ........................
## Boruta performed 130 randomForest runs in 35.56 secs.
##         8 attributes confirmed important: earlySE17 dce3SE4
## dce3SE13 maxCr_inside Slope_ini_countor max_F_r_i mean_F_r_i
## edge_sharp_std
##         167 attributes confirmed unimportant: V0 V1 V2 V3 V4 V5 V6
## V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE1
## earlySE2 earlySE3 earlySE4 earlySE5 earlySE6 earlySE7 earlySE8
## earlySE9 earlySE10 earlySE11 earlySE12 earlySE13 earlySE14
## earlySE15 earlySE16 earlySE18 dce2SE0 dce2SE1 dce2SE2 dce2SE3
## dce2SE4 dce2SE5 dce2SE6 dce2SE7 dce2SE8 dce2SE9 dce2SE10 dce2SE11
## dce2SE12 dce2SE14 dce2SE15 dce2SE17 dce2SE18 dce2SE19 dce3SE1
## dce3SE2 dce3SE3 dce3SE5 dce3SE7 dce3SE8 dce3SE9 dce3SE10 dce3SE11
## dce3SE12 dce3SE14 dce3SE15 dce3SE16 dce3SE17 dce3SE18 dce3SE19
## lateSE0 lateSE1 lateSE2 lateSE3 lateSE4 lateSE5 lateSE6 lateSE7
## lateSE8 lateSE9 lateSE10 lateSE11 lateSE12 lateSE13 lateSE14
## lateSE15 lateSE16 lateSE17 lateSE18 lateSE19 degreeC closenessC
## betweennessC no_triangles no_con_comp A_inside alpha_inside
## beta_inside iAUC1_inside Tpeak_inside Kpeak_inside SER_inside
## peakCr_inside UptakeRate_inside washoutRate_inside maxVr_inside
## peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside
## Vr_post_1_inside A_countor beta_countor iAUC1_countor
## Tpeak_countor Kpeak_countor SER_countor maxCr_countor
## peakCr_countor UptakeRate_countor washoutRate_countor
## maxVr_countor peakVr_countor Vr_increasingRate_countor
## Vr_decreasingRate_countor Vr_post_1_countor min_F_r_i var_F_r_i
## skew_F_r_i kurt_F_r_i iMax_Variance_uptake
## iiMin_change_Variance_uptake iiiMax_Margin_Gradient
## k_Max_Margin_Grad ivVariance circularity irregularity
## edge_sharp_mean max_RGH_mean max_RGH_mean_k max_RGH_var
## max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad
## texture_contrast_halfRad texture_contrast_threeQuaRad
## texture_homogeneity_zero texture_homogeneity_quarterRad
## texture_homogeneity_halfRad texture_homogeneity_threeQuaRad
## texture_dissimilarity_zero texture_dissimilarity_quarterRad
## texture_dissimilarity_halfRad texture_dissimilarity_threeQuaRad
## texture_correlation_zero texture_correlation_quarterRad
## texture_correlation_halfRad texture_correlation_threeQuaRad
## texture_ASM_zero texture_ASM_quarterRad texture_ASM_halfRad
## texture_ASM_threeQuaRad texture_energy_zero
## texture_energy_quarterRad texture_energy_halfRad
## texture_energy_threeQuaRad
##         7 tentative attributes left: earlySE19 dce2SE13 dce2SE16
## dce3SE0 dce3SE6 Slope_ini_inside alpha_countor
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-813.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-814.png) 

```
## D:  2 T:  5 
## Area under the curve: 0.862
## Area under the curve: 0.767
## D:  5 T:  5 
## Area under the curve: 0.953
## Area under the curve: 0.7
## D:  10 T:  5 
## Area under the curve: 0.976
## Area under the curve: 0.867
## D:  20 T:  5 
## Area under the curve: 0.969
## Area under the curve: 0.8
## D:  2 T:  10 
## Area under the curve: 0.901
## Area under the curve: 0.867
## D:  5 T:  10 
## Area under the curve: 0.992
## Area under the curve: 0.8
## D:  10 T:  10 
## Area under the curve: 0.98
## Area under the curve: 0.733
## D:  20 T:  10 
## Area under the curve: 0.967
## Area under the curve: 0.867
## D:  2 T:  30 
## Area under the curve: 0.932
## Area under the curve: 0.8
## D:  5 T:  30 
## Area under the curve: 0.986
## Area under the curve: 0.8
## D:  10 T:  30 
## Area under the curve: 0.995
## Area under the curve: 0.867
## D:  20 T:  30 
## Area under the curve: 0.996
## Area under the curve: 0.867
## D:  2 T:  60 
## Area under the curve: 0.923
## Area under the curve: 0.8
## D:  5 T:  60 
## Area under the curve: 0.986
## Area under the curve: 0.867
## D:  10 T:  60 
## Area under the curve: 0.992
## Area under the curve: 0.933
## D:  20 T:  60 
## Area under the curve: 0.993
## Area under the curve: 0.8
## D:  2 T:  100 
## Area under the curve: 0.946
## Area under the curve: 0.933
## D:  5 T:  100 
## Area under the curve: 0.994
## Area under the curve: 0.933
## D:  10 T:  100 
## Area under the curve: 0.992
## Area under the curve: 0.867
## D:  20 T:  100 
## Area under the curve: 0.996
## Area under the curve: 0.933
## D:  2 T:  250 
## Area under the curve: 0.929
## Area under the curve: 0.867
## D:  5 T:  250 
## Area under the curve: 0.993
## Area under the curve: 0.933
## D:  10 T:  250 
## Area under the curve: 0.994
## Area under the curve: 0.933
## D:  20 T:  250 
## Area under the curve: 0.991
## Area under the curve: 0.867
## D:  2 T:  500 
## Area under the curve: 0.932
## Area under the curve: 0.867
## D:  5 T:  500 
## Area under the curve: 0.991
## Area under the curve: 0.867
## D:  10 T:  500 
## Area under the curve: 0.995
## Area under the curve: 0.867
## D:  20 T:  500 
## Area under the curve: 0.992
## Area under the curve: 0.933
## D:  2 T:  750 
## Area under the curve: 0.931
## Area under the curve: 0.867
## D:  5 T:  750 
## Area under the curve: 0.992
## Area under the curve: 0.867
## D:  10 T:  750 
## Area under the curve: 0.993
## Area under the curve: 0.933
## D:  20 T:  750 
## Area under the curve: 0.994
## Area under the curve: 0.867
##     x   y acuTrain rocTrain senTrain speTrain acuTest rocTest senTest
## 1   2   5   0.8289   0.8624   0.9167   0.6786   0.750  0.7667     0.8
## 2   5   5   0.8947   0.9531   0.9375   0.8214   0.750  0.7000     0.8
## 3  10   5   0.9079   0.9762   0.9583   0.8214   0.750  0.8667     0.8
## 4  20   5   0.9342   0.9695   0.9375   0.9286   0.625  0.8000     0.6
## 5   2  10   0.8816   0.9014   0.9583   0.7500   0.750  0.8667     0.8
## 6   5  10   0.9211   0.9918   0.9792   0.8214   0.875  0.8000     0.8
## 7  10  10   0.9474   0.9799   1.0000   0.8571   0.750  0.7333     0.8
## 8  20  10   0.9079   0.9665   0.9792   0.7857   0.750  0.8667     0.8
## 9   2  30   0.8553   0.9315   0.9583   0.6786   0.750  0.8000     0.8
## 10  5  30   0.9474   0.9859   1.0000   0.8571   0.750  0.8000     0.8
## 11 10  30   0.9342   0.9948   0.9792   0.8571   0.750  0.8667     0.8
## 12 20  30   0.9605   0.9955   1.0000   0.8929   0.750  0.8667     0.8
## 13  2  60   0.8816   0.9226   0.9375   0.7857   0.750  0.8000     0.8
## 14  5  60   0.9342   0.9859   1.0000   0.8214   0.750  0.8667     0.8
## 15 10  60   0.9474   0.9918   1.0000   0.8571   0.750  0.9333     0.8
## 16 20  60   0.9342   0.9926   1.0000   0.8214   0.750  0.8000     0.8
## 17  2 100   0.8947   0.9457   0.9792   0.7500   0.750  0.9333     0.8
## 18  5 100   0.9474   0.9940   1.0000   0.8571   0.750  0.9333     0.8
## 19 10 100   0.9342   0.9918   1.0000   0.8214   0.750  0.8667     0.8
## 20 20 100   0.9474   0.9955   1.0000   0.8571   0.750  0.9333     0.8
## 21  2 250   0.8816   0.9293   0.9792   0.7143   0.750  0.8667     0.8
## 22  5 250   0.9474   0.9933   1.0000   0.8571   0.750  0.9333     0.8
## 23 10 250   0.9474   0.9940   1.0000   0.8571   0.750  0.9333     0.8
## 24 20 250   0.9342   0.9911   1.0000   0.8214   0.750  0.8667     0.8
## 25  2 500   0.8816   0.9315   0.9792   0.7143   0.750  0.8667     0.8
## 26  5 500   0.9342   0.9911   1.0000   0.8214   0.750  0.8667     0.8
## 27 10 500   0.9342   0.9948   1.0000   0.8214   0.750  0.8667     0.8
## 28 20 500   0.9474   0.9918   1.0000   0.8571   0.750  0.9333     0.8
## 29  2 750   0.8816   0.9308   0.9792   0.7143   0.750  0.8667     0.8
## 30  5 750   0.9474   0.9918   1.0000   0.8571   0.750  0.8667     0.8
## 31 10 750   0.9474   0.9926   1.0000   0.8571   0.750  0.9333     0.8
## 32 20 750   0.9474   0.9940   1.0000   0.8571   0.750  0.8667     0.8
##    speTest
## 1   0.6667
## 2   0.6667
## 3   0.6667
## 4   0.6667
## 5   0.6667
## 6   1.0000
## 7   0.6667
## 8   0.6667
## 9   0.6667
## 10  0.6667
## 11  0.6667
## 12  0.6667
## 13  0.6667
## 14  0.6667
## 15  0.6667
## 16  0.6667
## 17  0.6667
## 18  0.6667
## 19  0.6667
## 20  0.6667
## 21  0.6667
## 22  0.6667
## 23  0.6667
## 24  0.6667
## 25  0.6667
## 26  0.6667
## 27  0.6667
## 28  0.6667
## 29  0.6667
## 30  0.6667
## 31  0.6667
## 32  0.6667
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-815.png) 

```
## Initial round 1: ..........
##  108  attributes rejected after this test:  V0 V1 V2 V3 V4 V6 V7 V8 V9 V10 V12 V13 V14 V16 V17 V18 V19 earlySE0 earlySE3 earlySE7 earlySE9 earlySE10 earlySE12 dce2SE0 dce2SE1 dce2SE3 dce2SE5 dce2SE7 dce2SE9 dce2SE11 dce2SE14 dce2SE17 dce2SE19 dce3SE3 dce3SE7 dce3SE9 dce3SE10 dce3SE11 dce3SE12 dce3SE15 lateSE0 lateSE5 lateSE7 lateSE10 lateSE11 lateSE12 lateSE13 lateSE14 lateSE15 lateSE16 degreeC closenessC betweennessC no_triangles no_con_comp A_inside iAUC1_inside Tpeak_inside Kpeak_inside SER_inside peakCr_inside UptakeRate_inside washoutRate_inside peakVr_inside Vr_decreasingRate_inside Vr_post_1_inside iAUC1_countor Tpeak_countor Kpeak_countor peakCr_countor washoutRate_countor maxVr_countor peakVr_countor Vr_increasingRate_countor Vr_post_1_countor min_F_r_i skew_F_r_i kurt_F_r_i iiiMax_Margin_Gradient k_Max_Margin_Grad circularity irregularity edge_sharp_mean max_RGH_mean_k max_RGH_var max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad texture_contrast_halfRad texture_contrast_threeQuaRad texture_homogeneity_zero texture_homogeneity_quarterRad texture_homogeneity_halfRad texture_homogeneity_threeQuaRad texture_dissimilarity_zero texture_dissimilarity_quarterRad texture_dissimilarity_halfRad texture_dissimilarity_threeQuaRad texture_correlation_zero texture_correlation_quarterRad texture_correlation_halfRad texture_correlation_threeQuaRad texture_ASM_zero texture_ASM_threeQuaRad texture_energy_zero texture_energy_quarterRad texture_energy_halfRad texture_energy_threeQuaRad 
## 
## Initial round 2: ..........
##  23  attributes rejected after this test:  V5 V15 earlySE4 earlySE13 earlySE14 earlySE18 dce2SE2 dce2SE6 dce2SE10 dce2SE12 dce2SE18 dce3SE5 dce3SE18 lateSE1 lateSE8 lateSE17 lateSE18 lateSE19 maxVr_inside Vr_increasingRate_inside UptakeRate_countor iiMin_change_Variance_uptake ivVariance 
## 
## Initial round 3: ..........
##  6  attributes rejected after this test:  earlySE11 dce2SE8 dce3SE17 Vr_decreasingRate_countor iMax_Variance_uptake texture_ASM_halfRad 
## 
## Final round: ..........
##  16  attributes rejected after this test:  V11 earlySE19 dce2SE4 dce3SE1 dce3SE2 dce3SE8 dce3SE14 dce3SE16 dce3SE19 lateSE3 alpha_inside beta_inside beta_countor var_F_r_i max_RGH_mean texture_ASM_quarterRad 
## ....
##  6  attributes rejected after this test:  earlySE2 earlySE6 earlySE8 dce2SE13 lateSE4 lateSE9 
## ....
##  1  attributes confirmed after this test:  Slope_ini_countor 
## 
##  2  attributes rejected after this test:  earlySE5 maxCr_countor 
## ...
##  2  attributes confirmed after this test:  Slope_ini_inside alpha_countor 
## ......
##  1  attributes confirmed after this test:  edge_sharp_std 
## ...................
##  1  attributes confirmed after this test:  max_F_r_i 
## ...............
##  1  attributes rejected after this test:  earlySE16 
## ........
##  1  attributes confirmed after this test:  dce3SE13 
## ..
##  1  attributes confirmed after this test:  earlySE15 
## .............................
## Boruta performed 130 randomForest runs in 40.68 secs.
##         7 attributes confirmed important: earlySE15 dce3SE13
## Slope_ini_inside alpha_countor Slope_ini_countor max_F_r_i
## edge_sharp_std
##         162 attributes confirmed unimportant: V0 V1 V2 V3 V4 V5 V6
## V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE2
## earlySE3 earlySE4 earlySE5 earlySE6 earlySE7 earlySE8 earlySE9
## earlySE10 earlySE11 earlySE12 earlySE13 earlySE14 earlySE16
## earlySE18 earlySE19 dce2SE0 dce2SE1 dce2SE2 dce2SE3 dce2SE4
## dce2SE5 dce2SE6 dce2SE7 dce2SE8 dce2SE9 dce2SE10 dce2SE11 dce2SE12
## dce2SE13 dce2SE14 dce2SE17 dce2SE18 dce2SE19 dce3SE1 dce3SE2
## dce3SE3 dce3SE5 dce3SE7 dce3SE8 dce3SE9 dce3SE10 dce3SE11 dce3SE12
## dce3SE14 dce3SE15 dce3SE16 dce3SE17 dce3SE18 dce3SE19 lateSE0
## lateSE1 lateSE3 lateSE4 lateSE5 lateSE7 lateSE8 lateSE9 lateSE10
## lateSE11 lateSE12 lateSE13 lateSE14 lateSE15 lateSE16 lateSE17
## lateSE18 lateSE19 degreeC closenessC betweennessC no_triangles
## no_con_comp A_inside alpha_inside beta_inside iAUC1_inside
## Tpeak_inside Kpeak_inside SER_inside peakCr_inside
## UptakeRate_inside washoutRate_inside maxVr_inside peakVr_inside
## Vr_increasingRate_inside Vr_decreasingRate_inside Vr_post_1_inside
## beta_countor iAUC1_countor Tpeak_countor Kpeak_countor
## maxCr_countor peakCr_countor UptakeRate_countor
## washoutRate_countor maxVr_countor peakVr_countor
## Vr_increasingRate_countor Vr_decreasingRate_countor
## Vr_post_1_countor min_F_r_i var_F_r_i skew_F_r_i kurt_F_r_i
## iMax_Variance_uptake iiMin_change_Variance_uptake
## iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance circularity
## irregularity edge_sharp_mean max_RGH_mean max_RGH_mean_k
## max_RGH_var max_RGH_var_k texture_contrast_zero
## texture_contrast_quarterRad texture_contrast_halfRad
## texture_contrast_threeQuaRad texture_homogeneity_zero
## texture_homogeneity_quarterRad texture_homogeneity_halfRad
## texture_homogeneity_threeQuaRad texture_dissimilarity_zero
## texture_dissimilarity_quarterRad texture_dissimilarity_halfRad
## texture_dissimilarity_threeQuaRad texture_correlation_zero
## texture_correlation_quarterRad texture_correlation_halfRad
## texture_correlation_threeQuaRad texture_ASM_zero
## texture_ASM_quarterRad texture_ASM_halfRad texture_ASM_threeQuaRad
## texture_energy_zero texture_energy_quarterRad
## texture_energy_halfRad texture_energy_threeQuaRad
##         13 tentative attributes left: earlySE1 earlySE17 dce2SE15
## dce2SE16 dce3SE0 dce3SE4 dce3SE6 lateSE2 lateSE6 maxCr_inside
## A_countor SER_countor mean_F_r_i
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-816.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-817.png) 

```
## D:  2 T:  5 
## Area under the curve: 0.875
## Area under the curve: 0.889
## D:  5 T:  5 
## Area under the curve: 0.991
## Area under the curve: 0.639
## D:  10 T:  5 
## Area under the curve: 0.949
## Area under the curve: 0.611
## D:  20 T:  5 
## Area under the curve: 0.975
## Area under the curve: 0.75
## D:  2 T:  10 
## Area under the curve: 0.936
## Area under the curve: 0.611
## D:  5 T:  10 
## Area under the curve: 0.975
## Area under the curve: 0.583
## D:  10 T:  10 
## Area under the curve: 0.986
## Area under the curve: 0.611
## D:  20 T:  10 
## Area under the curve: 0.975
## Area under the curve: 0.778
## D:  2 T:  30 
## Area under the curve: 0.926
## Area under the curve: 0.694
## D:  5 T:  30 
## Area under the curve: 0.99
## Area under the curve: 0.611
## D:  10 T:  30 
## Area under the curve: 0.984
## Area under the curve: 0.611
## D:  20 T:  30 
## Area under the curve: 0.998
## Area under the curve: 0.667
## D:  2 T:  60 
## Area under the curve: 0.94
## Area under the curve: 0.667
## D:  5 T:  60 
## Area under the curve: 0.994
## Area under the curve: 0.778
## D:  10 T:  60 
## Area under the curve: 0.992
## Area under the curve: 0.722
## D:  20 T:  60 
## Area under the curve: 0.997
## Area under the curve: 0.722
## D:  2 T:  100 
## Area under the curve: 0.929
## Area under the curve: 0.722
## D:  5 T:  100 
## Area under the curve: 0.993
## Area under the curve: 0.833
## D:  10 T:  100 
## Area under the curve: 0.996
## Area under the curve: 0.611
## D:  20 T:  100 
## Area under the curve: 0.998
## Area under the curve: 0.667
## D:  2 T:  250 
## Area under the curve: 0.935
## Area under the curve: 0.667
## D:  5 T:  250 
## Area under the curve: 0.995
## Area under the curve: 0.778
## D:  10 T:  250 
## Area under the curve: 0.994
## Area under the curve: 0.778
## D:  20 T:  250 
## Area under the curve: 0.995
## Area under the curve: 0.722
## D:  2 T:  500 
## Area under the curve: 0.932
## Area under the curve: 0.667
## D:  5 T:  500 
## Area under the curve: 0.995
## Area under the curve: 0.778
## D:  10 T:  500 
## Area under the curve: 0.995
## Area under the curve: 0.778
## D:  20 T:  500 
## Area under the curve: 0.993
## Area under the curve: 0.722
## D:  2 T:  750 
## Area under the curve: 0.935
## Area under the curve: 0.667
## D:  5 T:  750 
## Area under the curve: 0.995
## Area under the curve: 0.667
## D:  10 T:  750 
## Area under the curve: 0.994
## Area under the curve: 0.722
## D:  20 T:  750 
## Area under the curve: 0.994
## Area under the curve: 0.667
##     x   y acuTrain rocTrain senTrain speTrain acuTest rocTest senTest
## 1   2   5   0.7733   0.8746   0.8936   0.5714  0.7778  0.8889  1.0000
## 2   5   5   0.9333   0.9905   0.9149   0.9643  0.6667  0.6389  0.8333
## 3  10   5   0.8933   0.9487   0.8936   0.8929  0.7778  0.6111  1.0000
## 4  20   5   0.8800   0.9749   0.9574   0.7500  0.6667  0.7500  1.0000
## 5   2  10   0.8267   0.9358   0.9149   0.6786  0.7778  0.6111  1.0000
## 6   5  10   0.8933   0.9749   0.9149   0.8571  0.6667  0.5833  0.8333
## 7  10  10   0.9333   0.9856   0.9574   0.8929  0.7778  0.6111  0.8333
## 8  20  10   0.9200   0.9749   0.9787   0.8214  0.6667  0.7778  1.0000
## 9   2  30   0.8533   0.9255   0.9574   0.6786  0.7778  0.6944  1.0000
## 10  5  30   0.9600   0.9901   1.0000   0.8929  0.6667  0.6111  0.8333
## 11 10  30   0.9067   0.9840   0.9787   0.7857  0.7778  0.6111  1.0000
## 12 20  30   0.9867   0.9977   1.0000   0.9643  0.7778  0.6667  1.0000
## 13  2  60   0.8667   0.9400   0.9787   0.6786  0.6667  0.6667  1.0000
## 14  5  60   0.9333   0.9939   1.0000   0.8214  0.7778  0.7778  1.0000
## 15 10  60   0.9867   0.9924   1.0000   0.9643  0.6667  0.7222  1.0000
## 16 20  60   0.9467   0.9970   1.0000   0.8571  0.7778  0.7222  1.0000
## 17  2 100   0.8800   0.9293   0.9787   0.7143  0.6667  0.7222  1.0000
## 18  5 100   0.9467   0.9932   1.0000   0.8571  0.7778  0.8333  1.0000
## 19 10 100   0.9600   0.9962   1.0000   0.8929  0.7778  0.6111  1.0000
## 20 20 100   0.9733   0.9985   1.0000   0.9286  0.7778  0.6667  1.0000
## 21  2 250   0.8267   0.9347   0.9149   0.6786  0.7778  0.6667  1.0000
## 22  5 250   0.9600   0.9954   1.0000   0.8929  0.7778  0.7778  1.0000
## 23 10 250   0.9600   0.9939   1.0000   0.8929  0.7778  0.7778  1.0000
## 24 20 250   0.9600   0.9947   1.0000   0.8929  0.7778  0.7222  1.0000
## 25  2 500   0.8400   0.9316   0.9362   0.6786  0.7778  0.6667  1.0000
## 26  5 500   0.9600   0.9954   1.0000   0.8929  0.7778  0.7778  1.0000
## 27 10 500   0.9467   0.9947   1.0000   0.8571  0.7778  0.7778  1.0000
## 28 20 500   0.9467   0.9932   1.0000   0.8571  0.7778  0.7222  1.0000
## 29  2 750   0.8267   0.9347   0.9149   0.6786  0.7778  0.6667  1.0000
## 30  5 750   0.9467   0.9947   1.0000   0.8571  0.7778  0.6667  1.0000
## 31 10 750   0.9600   0.9939   1.0000   0.8929  0.7778  0.7222  1.0000
## 32 20 750   0.9600   0.9939   1.0000   0.8929  0.7778  0.6667  1.0000
##    speTest
## 1   0.3333
## 2   0.3333
## 3   0.3333
## 4   0.0000
## 5   0.3333
## 6   0.3333
## 7   0.6667
## 8   0.0000
## 9   0.3333
## 10  0.3333
## 11  0.3333
## 12  0.3333
## 13  0.0000
## 14  0.3333
## 15  0.0000
## 16  0.3333
## 17  0.0000
## 18  0.3333
## 19  0.3333
## 20  0.3333
## 21  0.3333
## 22  0.3333
## 23  0.3333
## 24  0.3333
## 25  0.3333
## 26  0.3333
## 27  0.3333
## 28  0.3333
## 29  0.3333
## 30  0.3333
## 31  0.3333
## 32  0.3333
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-818.png) 

```
## Initial round 1: ..........
##  125  attributes rejected after this test:  V0 V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE4 earlySE6 earlySE10 earlySE12 earlySE14 earlySE18 dce2SE1 dce2SE2 dce2SE3 dce2SE4 dce2SE6 dce2SE8 dce2SE9 dce2SE10 dce2SE11 dce2SE14 dce2SE17 dce2SE18 dce2SE19 dce3SE2 dce3SE3 dce3SE7 dce3SE9 dce3SE10 dce3SE11 dce3SE14 dce3SE15 dce3SE17 dce3SE18 lateSE0 lateSE1 lateSE5 lateSE8 lateSE10 lateSE11 lateSE12 lateSE14 lateSE15 lateSE16 lateSE17 lateSE18 lateSE19 degreeC closenessC betweennessC no_triangles no_con_comp A_inside beta_inside iAUC1_inside Tpeak_inside Kpeak_inside peakCr_inside washoutRate_inside maxVr_inside peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside Vr_post_1_inside iAUC1_countor Tpeak_countor Kpeak_countor SER_countor peakCr_countor UptakeRate_countor washoutRate_countor maxVr_countor peakVr_countor Vr_increasingRate_countor Vr_decreasingRate_countor Vr_post_1_countor min_F_r_i skew_F_r_i kurt_F_r_i iMax_Variance_uptake iiMin_change_Variance_uptake iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance irregularity max_RGH_mean_k max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad texture_contrast_halfRad texture_contrast_threeQuaRad texture_homogeneity_zero texture_homogeneity_quarterRad texture_homogeneity_halfRad texture_homogeneity_threeQuaRad texture_dissimilarity_zero texture_dissimilarity_quarterRad texture_dissimilarity_halfRad texture_dissimilarity_threeQuaRad texture_correlation_zero texture_correlation_quarterRad texture_correlation_halfRad texture_correlation_threeQuaRad texture_ASM_zero texture_ASM_halfRad texture_ASM_threeQuaRad texture_energy_quarterRad texture_energy_halfRad texture_energy_threeQuaRad 
## 
## Initial round 2: ..........
##  5  attributes rejected after this test:  lateSE7 lateSE13 alpha_inside circularity texture_energy_zero 
## 
## Initial round 3: ..........
##  5  attributes rejected after this test:  earlySE2 earlySE11 lateSE3 lateSE4 A_countor 
## 
## Final round: ..........
##  13  attributes rejected after this test:  earlySE13 earlySE16 dce2SE0 dce2SE5 dce2SE7 dce2SE12 dce2SE15 dce3SE12 dce3SE19 SER_inside beta_countor maxCr_countor edge_sharp_mean 
## ....
##  8  attributes rejected after this test:  earlySE7 earlySE17 earlySE19 dce3SE0 dce3SE5 dce3SE8 lateSE9 texture_ASM_quarterRad 
## ..........
##  1  attributes confirmed after this test:  edge_sharp_std 
## 
##  2  attributes rejected after this test:  earlySE8 lateSE2 
## ...
##  1  attributes rejected after this test:  lateSE6 
## ...
##  1  attributes rejected after this test:  dce2SE13 
## ................
##  1  attributes confirmed after this test:  max_F_r_i 
## ........
##  1  attributes confirmed after this test:  mean_F_r_i 
## .......
##  3  attributes confirmed after this test:  alpha_countor Slope_ini_countor max_RGH_var 
## ...
##  1  attributes confirmed after this test:  dce3SE13 
## ..........
##  1  attributes confirmed after this test:  dce2SE16 
## .......
##  1  attributes confirmed after this test:  UptakeRate_inside 
## .....
##  1  attributes confirmed after this test:  var_F_r_i 
## ..
##  2  attributes confirmed after this test:  earlySE1 max_RGH_mean 
## .....
##  1  attributes confirmed after this test:  Slope_ini_inside 
## .......
## Boruta performed 130 randomForest runs in 1.408 mins.
##         13 attributes confirmed important: earlySE1 dce2SE16
## dce3SE13 Slope_ini_inside UptakeRate_inside alpha_countor
## Slope_ini_countor max_F_r_i mean_F_r_i var_F_r_i edge_sharp_std
## max_RGH_mean max_RGH_var
##         160 attributes confirmed unimportant: V0 V1 V2 V3 V4 V5 V6
## V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE2
## earlySE4 earlySE6 earlySE7 earlySE8 earlySE10 earlySE11 earlySE12
## earlySE13 earlySE14 earlySE16 earlySE17 earlySE18 earlySE19
## dce2SE0 dce2SE1 dce2SE2 dce2SE3 dce2SE4 dce2SE5 dce2SE6 dce2SE7
## dce2SE8 dce2SE9 dce2SE10 dce2SE11 dce2SE12 dce2SE13 dce2SE14
## dce2SE15 dce2SE17 dce2SE18 dce2SE19 dce3SE0 dce3SE2 dce3SE3
## dce3SE5 dce3SE7 dce3SE8 dce3SE9 dce3SE10 dce3SE11 dce3SE12
## dce3SE14 dce3SE15 dce3SE17 dce3SE18 dce3SE19 lateSE0 lateSE1
## lateSE2 lateSE3 lateSE4 lateSE5 lateSE6 lateSE7 lateSE8 lateSE9
## lateSE10 lateSE11 lateSE12 lateSE13 lateSE14 lateSE15 lateSE16
## lateSE17 lateSE18 lateSE19 degreeC closenessC betweennessC
## no_triangles no_con_comp A_inside alpha_inside beta_inside
## iAUC1_inside Tpeak_inside Kpeak_inside SER_inside peakCr_inside
## washoutRate_inside maxVr_inside peakVr_inside
## Vr_increasingRate_inside Vr_decreasingRate_inside Vr_post_1_inside
## A_countor beta_countor iAUC1_countor Tpeak_countor Kpeak_countor
## SER_countor maxCr_countor peakCr_countor UptakeRate_countor
## washoutRate_countor maxVr_countor peakVr_countor
## Vr_increasingRate_countor Vr_decreasingRate_countor
## Vr_post_1_countor min_F_r_i skew_F_r_i kurt_F_r_i
## iMax_Variance_uptake iiMin_change_Variance_uptake
## iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance circularity
## irregularity edge_sharp_mean max_RGH_mean_k max_RGH_var_k
## texture_contrast_zero texture_contrast_quarterRad
## texture_contrast_halfRad texture_contrast_threeQuaRad
## texture_homogeneity_zero texture_homogeneity_quarterRad
## texture_homogeneity_halfRad texture_homogeneity_threeQuaRad
## texture_dissimilarity_zero texture_dissimilarity_quarterRad
## texture_dissimilarity_halfRad texture_dissimilarity_threeQuaRad
## texture_correlation_zero texture_correlation_quarterRad
## texture_correlation_halfRad texture_correlation_threeQuaRad
## texture_ASM_zero texture_ASM_quarterRad texture_ASM_halfRad
## texture_ASM_threeQuaRad texture_energy_zero
## texture_energy_quarterRad texture_energy_halfRad
## texture_energy_threeQuaRad
##         9 tentative attributes left: earlySE3 earlySE5 earlySE9
## earlySE15 dce3SE1 dce3SE4 dce3SE6 dce3SE16 maxCr_inside
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-819.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-820.png) 

```
## D:  2 T:  5 
## Area under the curve: 0.939
## Area under the curve: 0.6
## D:  5 T:  5 
## Area under the curve: 0.985
## Area under the curve: 0.667
## D:  10 T:  5 
## Area under the curve: 0.954
## Area under the curve: 0.933
## D:  20 T:  5 
## Area under the curve: 0.968
## Area under the curve: 0.6
## D:  2 T:  10 
## Area under the curve: 0.912
## Area under the curve: 0.8
## D:  5 T:  10 
## Area under the curve: 0.985
## Area under the curve: 0.667
## D:  10 T:  10 
## Area under the curve: 0.984
## Area under the curve: 0.6
## D:  20 T:  10 
## Area under the curve: 0.943
## Area under the curve: 0.533
## D:  2 T:  30 
## Area under the curve: 0.94
## Area under the curve: 0.667
## D:  5 T:  30 
## Area under the curve: 0.996
## Area under the curve: 0.667
## D:  10 T:  30 
## Area under the curve: 0.995
## Area under the curve: 0.6
## D:  20 T:  30 
## Area under the curve: 0.997
## Area under the curve: 0.467
## D:  2 T:  60 
## Area under the curve: 0.93
## Area under the curve: 0.733
## D:  5 T:  60 
## Area under the curve: 0.993
## Area under the curve: 0.533
## D:  10 T:  60 
## Area under the curve: 0.986
## Area under the curve: 0.667
## D:  20 T:  60 
## Area under the curve: 0.995
## Area under the curve: 0.6
## D:  2 T:  100 
## Area under the curve: 0.935
## Area under the curve: 0.6
## D:  5 T:  100 
## Area under the curve: 0.992
## Area under the curve: 0.667
## D:  10 T:  100 
## Area under the curve: 0.996
## Area under the curve: 0.667
## D:  20 T:  100 
## Area under the curve: 0.997
## Area under the curve: 0.8
## D:  2 T:  250 
## Area under the curve: 0.944
## Area under the curve: 0.667
## D:  5 T:  250 
## Area under the curve: 0.998
## Area under the curve: 0.6
## D:  10 T:  250 
## Area under the curve: 0.995
## Area under the curve: 0.467
## D:  20 T:  250 
## Area under the curve: 0.996
## Area under the curve: 0.733
## D:  2 T:  500 
## Area under the curve: 0.939
## Area under the curve: 0.667
## D:  5 T:  500 
## Area under the curve: 0.994
## Area under the curve: 0.667
## D:  10 T:  500 
## Area under the curve: 0.996
## Area under the curve: 0.667
## D:  20 T:  500 
## Area under the curve: 0.997
## Area under the curve: 0.667
## D:  2 T:  750 
## Area under the curve: 0.938
## Area under the curve: 0.733
## D:  5 T:  750 
## Area under the curve: 0.997
## Area under the curve: 0.667
## D:  10 T:  750 
## Area under the curve: 0.996
## Area under the curve: 0.733
## D:  20 T:  750 
## Area under the curve: 0.996
## Area under the curve: 0.667
##     x   y acuTrain rocTrain senTrain speTrain acuTest rocTest senTest
## 1   2   5   0.8289   0.9394   0.9792   0.5714   0.500  0.6000     0.6
## 2   5   5   0.9737   0.9851   0.9792   0.9643   0.625  0.6667     0.6
## 3  10   5   0.9342   0.9542   1.0000   0.8214   0.875  0.9333     1.0
## 4  20   5   0.9079   0.9684   1.0000   0.7500   0.375  0.6000     0.6
## 5   2  10   0.8421   0.9118   0.9375   0.6786   0.750  0.8000     0.8
## 6   5  10   0.9474   0.9851   0.9583   0.9286   0.750  0.6667     0.6
## 7  10  10   0.9211   0.9836   0.9792   0.8214   0.625  0.6000     0.8
## 8  20  10   0.8421   0.9427   0.9167   0.7143   0.500  0.5333     0.6
## 9   2  30   0.8553   0.9397   1.0000   0.6071   0.500  0.6667     0.8
## 10  5  30   0.9605   0.9955   1.0000   0.8929   0.625  0.6667     0.8
## 11 10  30   0.9737   0.9948   1.0000   0.9286   0.625  0.6000     0.8
## 12 20  30   0.9737   0.9970   1.0000   0.9286   0.500  0.4667     0.6
## 13  2  60   0.8553   0.9301   0.9792   0.6429   0.625  0.7333     0.8
## 14  5  60   0.9474   0.9933   1.0000   0.8571   0.500  0.5333     0.8
## 15 10  60   0.9474   0.9859   1.0000   0.8571   0.500  0.6667     0.6
## 16 20  60   0.9211   0.9948   1.0000   0.7857   0.625  0.6000     0.8
## 17  2 100   0.8421   0.9353   0.9792   0.6071   0.500  0.6000     0.8
## 18  5 100   0.9605   0.9918   1.0000   0.8929   0.750  0.6667     0.8
## 19 10 100   0.9737   0.9955   1.0000   0.9286   0.625  0.6667     0.8
## 20 20 100   0.9474   0.9970   1.0000   0.8571   0.625  0.8000     0.8
## 21  2 250   0.8684   0.9442   1.0000   0.6429   0.500  0.6667     0.8
## 22  5 250   0.9737   0.9978   1.0000   0.9286   0.625  0.6000     0.8
## 23 10 250   0.9474   0.9948   1.0000   0.8571   0.625  0.4667     0.8
## 24 20 250   0.9605   0.9963   1.0000   0.8929   0.625  0.7333     0.8
## 25  2 500   0.8684   0.9390   1.0000   0.6429   0.625  0.6667     0.8
## 26  5 500   0.9474   0.9940   1.0000   0.8571   0.625  0.6667     0.8
## 27 10 500   0.9737   0.9963   1.0000   0.9286   0.625  0.6667     0.8
## 28 20 500   0.9605   0.9970   1.0000   0.8929   0.625  0.6667     0.8
## 29  2 750   0.8684   0.9375   1.0000   0.6429   0.625  0.7333     0.8
## 30  5 750   0.9605   0.9970   1.0000   0.8929   0.625  0.6667     0.8
## 31 10 750   0.9342   0.9963   1.0000   0.8214   0.500  0.7333     0.8
## 32 20 750   0.9605   0.9955   1.0000   0.8929   0.625  0.6667     0.8
##    speTest
## 1   0.3333
## 2   0.6667
## 3   0.6667
## 4   0.0000
## 5   0.6667
## 6   1.0000
## 7   0.3333
## 8   0.3333
## 9   0.0000
## 10  0.3333
## 11  0.3333
## 12  0.3333
## 13  0.3333
## 14  0.0000
## 15  0.3333
## 16  0.3333
## 17  0.0000
## 18  0.6667
## 19  0.3333
## 20  0.3333
## 21  0.0000
## 22  0.3333
## 23  0.3333
## 24  0.3333
## 25  0.3333
## 26  0.3333
## 27  0.3333
## 28  0.3333
## 29  0.3333
## 30  0.3333
## 31  0.0000
## 32  0.3333
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-821.png) 

```
## Initial round 1: ..........
##  107  attributes rejected after this test:  V0 V1 V2 V3 V4 V6 V7 V8 V9 V10 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE4 earlySE9 earlySE10 earlySE14 dce2SE3 dce2SE7 dce2SE8 dce2SE9 dce2SE10 dce2SE11 dce2SE12 dce2SE15 dce2SE18 dce3SE3 dce3SE7 dce3SE9 dce3SE10 dce3SE11 dce3SE12 dce3SE14 dce3SE15 dce3SE18 lateSE0 lateSE4 lateSE5 lateSE7 lateSE9 lateSE10 lateSE11 lateSE12 lateSE14 degreeC closenessC betweennessC no_triangles no_con_comp A_inside alpha_inside beta_inside iAUC1_inside Tpeak_inside Kpeak_inside peakCr_inside washoutRate_inside peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside Vr_post_1_inside iAUC1_countor Kpeak_countor peakCr_countor UptakeRate_countor washoutRate_countor maxVr_countor peakVr_countor Vr_increasingRate_countor Vr_decreasingRate_countor Vr_post_1_countor skew_F_r_i kurt_F_r_i iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance circularity irregularity max_RGH_mean max_RGH_mean_k max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad texture_contrast_halfRad texture_contrast_threeQuaRad texture_homogeneity_zero texture_homogeneity_quarterRad texture_homogeneity_halfRad texture_homogeneity_threeQuaRad texture_dissimilarity_zero texture_dissimilarity_quarterRad texture_dissimilarity_halfRad texture_correlation_zero texture_correlation_quarterRad texture_correlation_threeQuaRad texture_ASM_zero texture_ASM_quarterRad texture_ASM_threeQuaRad texture_energy_zero texture_energy_halfRad texture_energy_threeQuaRad 
## 
## Initial round 2: ..........
##  22  attributes rejected after this test:  V5 V11 earlySE2 earlySE3 dce2SE6 dce3SE2 dce3SE5 dce3SE8 lateSE1 lateSE2 lateSE3 lateSE13 lateSE15 lateSE16 maxVr_inside beta_countor Tpeak_countor maxCr_countor min_F_r_i edge_sharp_mean texture_dissimilarity_threeQuaRad texture_ASM_halfRad 
## 
## Initial round 3: ..........
##  13  attributes rejected after this test:  earlySE6 earlySE11 earlySE12 dce2SE0 dce2SE1 dce2SE5 dce2SE14 lateSE18 lateSE19 A_countor var_F_r_i texture_correlation_halfRad texture_energy_quarterRad 
## 
## Final round: ..........
##  9  attributes rejected after this test:  earlySE7 earlySE13 dce2SE4 dce2SE19 dce3SE0 dce3SE6 lateSE8 lateSE17 SER_inside 
## ....
##  2  attributes confirmed after this test:  earlySE17 dce2SE16 
## 
##  4  attributes rejected after this test:  earlySE18 dce2SE2 dce2SE17 lateSE6 
## ..........
##  2  attributes confirmed after this test:  Slope_ini_countor max_F_r_i 
## 
##  1  attributes rejected after this test:  alpha_countor 
## ...
##  2  attributes confirmed after this test:  dce3SE13 Slope_ini_inside 
## ......
##  1  attributes rejected after this test:  SER_countor 
## ..
##  2  attributes rejected after this test:  iiMin_change_Variance_uptake max_RGH_var 
## ...
##  1  attributes rejected after this test:  maxCr_inside 
## .............................................
##  1  attributes confirmed after this test:  edge_sharp_std 
## .................
## Boruta performed 130 randomForest runs in 41.67 secs.
##         7 attributes confirmed important: earlySE17 dce2SE16
## dce3SE13 Slope_ini_inside Slope_ini_countor max_F_r_i
## edge_sharp_std
##         160 attributes confirmed unimportant: V0 V1 V2 V3 V4 V5 V6
## V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE2
## earlySE3 earlySE4 earlySE6 earlySE7 earlySE9 earlySE10 earlySE11
## earlySE12 earlySE13 earlySE14 earlySE18 dce2SE0 dce2SE1 dce2SE2
## dce2SE3 dce2SE4 dce2SE5 dce2SE6 dce2SE7 dce2SE8 dce2SE9 dce2SE10
## dce2SE11 dce2SE12 dce2SE14 dce2SE15 dce2SE17 dce2SE18 dce2SE19
## dce3SE0 dce3SE2 dce3SE3 dce3SE5 dce3SE6 dce3SE7 dce3SE8 dce3SE9
## dce3SE10 dce3SE11 dce3SE12 dce3SE14 dce3SE15 dce3SE18 lateSE0
## lateSE1 lateSE2 lateSE3 lateSE4 lateSE5 lateSE6 lateSE7 lateSE8
## lateSE9 lateSE10 lateSE11 lateSE12 lateSE13 lateSE14 lateSE15
## lateSE16 lateSE17 lateSE18 lateSE19 degreeC closenessC
## betweennessC no_triangles no_con_comp A_inside alpha_inside
## beta_inside iAUC1_inside Tpeak_inside Kpeak_inside SER_inside
## maxCr_inside peakCr_inside washoutRate_inside maxVr_inside
## peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside
## Vr_post_1_inside A_countor alpha_countor beta_countor
## iAUC1_countor Tpeak_countor Kpeak_countor SER_countor
## maxCr_countor peakCr_countor UptakeRate_countor
## washoutRate_countor maxVr_countor peakVr_countor
## Vr_increasingRate_countor Vr_decreasingRate_countor
## Vr_post_1_countor min_F_r_i var_F_r_i skew_F_r_i kurt_F_r_i
## iiMin_change_Variance_uptake iiiMax_Margin_Gradient
## k_Max_Margin_Grad ivVariance circularity irregularity
## edge_sharp_mean max_RGH_mean max_RGH_mean_k max_RGH_var
## max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad
## texture_contrast_halfRad texture_contrast_threeQuaRad
## texture_homogeneity_zero texture_homogeneity_quarterRad
## texture_homogeneity_halfRad texture_homogeneity_threeQuaRad
## texture_dissimilarity_zero texture_dissimilarity_quarterRad
## texture_dissimilarity_halfRad texture_dissimilarity_threeQuaRad
## texture_correlation_zero texture_correlation_quarterRad
## texture_correlation_halfRad texture_correlation_threeQuaRad
## texture_ASM_zero texture_ASM_quarterRad texture_ASM_halfRad
## texture_ASM_threeQuaRad texture_energy_zero
## texture_energy_quarterRad texture_energy_halfRad
## texture_energy_threeQuaRad
##         15 tentative attributes left: earlySE1 earlySE5 earlySE8
## earlySE15 earlySE16 earlySE19 dce2SE13 dce3SE1 dce3SE4 dce3SE16
## dce3SE17 dce3SE19 UptakeRate_inside mean_F_r_i
## iMax_Variance_uptake
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-822.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-823.png) 

```
## D:  2 T:  5 
## Area under the curve: 0.871
## Area under the curve: 0.556
## D:  5 T:  5 
## Area under the curve: 0.943
## Area under the curve: 0.889
## D:  10 T:  5 
## Area under the curve: 0.959
## Area under the curve: 0.639
## D:  20 T:  5 
## Area under the curve: 0.956
## Area under the curve: 0.833
## D:  2 T:  10 
## Area under the curve: 0.891
## Area under the curve: 0.556
## D:  5 T:  10 
## Area under the curve: 0.955
## Area under the curve: 0.611
## D:  10 T:  10 
## Area under the curve: 0.977
## Area under the curve: 0.583
## D:  20 T:  10 
## Area under the curve: 0.988
## Area under the curve: 0.611
## D:  2 T:  30 
## Area under the curve: 0.923
## Area under the curve: 0.722
## D:  5 T:  30 
## Area under the curve: 0.967
## Area under the curve: 0.611
## D:  10 T:  30 
## Area under the curve: 0.992
## Area under the curve: 0.778
## D:  20 T:  30 
## Area under the curve: 0.984
## Area under the curve: 0.667
## D:  2 T:  60 
## Area under the curve: 0.929
## Area under the curve: 0.611
## D:  5 T:  60 
## Area under the curve: 0.973
## Area under the curve: 0.667
## D:  10 T:  60 
## Area under the curve: 0.983
## Area under the curve: 0.444
## D:  20 T:  60 
## Area under the curve: 0.98
## Area under the curve: 0.556
## D:  2 T:  100 
## Area under the curve: 0.932
## Area under the curve: 0.667
## D:  5 T:  100 
## Area under the curve: 0.988
## Area under the curve: 0.5
## D:  10 T:  100 
## Area under the curve: 0.982
## Area under the curve: 0.611
## D:  20 T:  100 
## Area under the curve: 0.995
## Area under the curve: 0.556
## D:  2 T:  250 
## Area under the curve: 0.938
## Area under the curve: 0.444
## D:  5 T:  250 
## Area under the curve: 0.991
## Area under the curve: 0.444
## D:  10 T:  250 
## Area under the curve: 0.992
## Area under the curve: 0.444
## D:  20 T:  250 
## Area under the curve: 0.993
## Area under the curve: 0.556
## D:  2 T:  500 
## Area under the curve: 0.934
## Area under the curve: 0.611
## D:  5 T:  500 
## Area under the curve: 0.985
## Area under the curve: 0.444
## D:  10 T:  500 
## Area under the curve: 0.992
## Area under the curve: 0.556
## D:  20 T:  500 
## Area under the curve: 0.991
## Area under the curve: 0.5
## D:  2 T:  750 
## Area under the curve: 0.935
## Area under the curve: 0.5
## D:  5 T:  750 
## Area under the curve: 0.988
## Area under the curve: 0.444
## D:  10 T:  750 
## Area under the curve: 0.995
## Area under the curve: 0.611
## D:  20 T:  750 
## Area under the curve: 0.992
## Area under the curve: 0.5
##     x   y acuTrain rocTrain senTrain speTrain acuTest rocTest senTest
## 1   2   5   0.8533   0.8712   0.9149   0.7500  0.4444  0.5556  0.6667
## 2   5   5   0.8933   0.9434   0.9362   0.8214  0.5556  0.8889  0.8333
## 3  10   5   0.8800   0.9586   0.9362   0.7857  0.6667  0.6389  1.0000
## 4  20   5   0.8667   0.9559   0.8936   0.8214  0.7778  0.8333  0.8333
## 5   2  10   0.8533   0.8906   0.9787   0.6429  0.4444  0.5556  0.6667
## 6   5  10   0.8533   0.9552   0.9149   0.7500  0.4444  0.6111  0.6667
## 7  10  10   0.8933   0.9772   0.9787   0.7500  0.4444  0.5833  0.6667
## 8  20  10   0.9200   0.9878   0.9787   0.8214  0.4444  0.6111  0.6667
## 9   2  30   0.8667   0.9233   0.9362   0.7500  0.4444  0.7222  0.6667
## 10  5  30   0.9333   0.9673   1.0000   0.8214  0.4444  0.6111  0.6667
## 11 10  30   0.9067   0.9924   0.9787   0.7857  0.5556  0.7778  0.8333
## 12 20  30   0.9333   0.9840   1.0000   0.8214  0.5556  0.6667  0.6667
## 13  2  60   0.8667   0.9293   0.9362   0.7500  0.4444  0.6111  0.6667
## 14  5  60   0.9067   0.9734   1.0000   0.7500  0.4444  0.6667  0.6667
## 15 10  60   0.9200   0.9833   1.0000   0.7857  0.4444  0.4444  0.6667
## 16 20  60   0.9200   0.9802   0.9787   0.8214  0.4444  0.5556  0.6667
## 17  2 100   0.8800   0.9324   0.9574   0.7500  0.4444  0.6667  0.6667
## 18  5 100   0.9200   0.9878   1.0000   0.7857  0.4444  0.5000  0.6667
## 19 10 100   0.9333   0.9818   1.0000   0.8214  0.4444  0.6111  0.6667
## 20 20 100   0.9200   0.9954   1.0000   0.7857  0.4444  0.5556  0.6667
## 21  2 250   0.8800   0.9384   0.9787   0.7143  0.4444  0.4444  0.6667
## 22  5 250   0.9200   0.9909   1.0000   0.7857  0.4444  0.4444  0.6667
## 23 10 250   0.9200   0.9924   1.0000   0.7857  0.4444  0.4444  0.6667
## 24 20 250   0.9200   0.9932   1.0000   0.7857  0.4444  0.5556  0.6667
## 25  2 500   0.8800   0.9339   0.9574   0.7500  0.4444  0.6111  0.6667
## 26  5 500   0.9333   0.9848   1.0000   0.8214  0.4444  0.4444  0.6667
## 27 10 500   0.9200   0.9916   1.0000   0.7857  0.4444  0.5556  0.6667
## 28 20 500   0.9200   0.9909   1.0000   0.7857  0.4444  0.5000  0.6667
## 29  2 750   0.8667   0.9347   0.9574   0.7143  0.4444  0.5000  0.6667
## 30  5 750   0.9333   0.9878   1.0000   0.8214  0.4444  0.4444  0.6667
## 31 10 750   0.9333   0.9947   1.0000   0.8214  0.4444  0.6111  0.6667
## 32 20 750   0.9200   0.9916   1.0000   0.7857  0.4444  0.5000  0.6667
##    speTest
## 1   0.0000
## 2   0.0000
## 3   0.0000
## 4   0.6667
## 5   0.0000
## 6   0.0000
## 7   0.0000
## 8   0.0000
## 9   0.0000
## 10  0.0000
## 11  0.0000
## 12  0.3333
## 13  0.0000
## 14  0.0000
## 15  0.0000
## 16  0.0000
## 17  0.0000
## 18  0.0000
## 19  0.0000
## 20  0.0000
## 21  0.0000
## 22  0.0000
## 23  0.0000
## 24  0.0000
## 25  0.0000
## 26  0.0000
## 27  0.0000
## 28  0.0000
## 29  0.0000
## 30  0.0000
## 31  0.0000
## 32  0.0000
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-824.png) 

```
## Initial round 1: ..........
##  110  attributes rejected after this test:  V0 V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE2 earlySE3 earlySE4 earlySE10 earlySE12 earlySE14 dce2SE0 dce2SE2 dce2SE3 dce2SE4 dce2SE5 dce2SE8 dce2SE10 dce2SE11 dce2SE17 dce3SE3 dce3SE5 dce3SE7 dce3SE10 dce3SE11 dce3SE17 lateSE3 lateSE5 lateSE10 lateSE11 lateSE13 lateSE15 lateSE16 degreeC closenessC betweennessC no_triangles no_con_comp A_inside alpha_inside beta_inside iAUC1_inside Tpeak_inside SER_inside peakCr_inside washoutRate_inside maxVr_inside peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside Vr_post_1_inside beta_countor iAUC1_countor Tpeak_countor peakCr_countor washoutRate_countor maxVr_countor peakVr_countor Vr_increasingRate_countor Vr_decreasingRate_countor Vr_post_1_countor min_F_r_i skew_F_r_i kurt_F_r_i iMax_Variance_uptake iiMin_change_Variance_uptake iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance irregularity edge_sharp_mean max_RGH_mean_k max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad texture_contrast_halfRad texture_contrast_threeQuaRad texture_homogeneity_zero texture_homogeneity_quarterRad texture_homogeneity_halfRad texture_homogeneity_threeQuaRad texture_dissimilarity_zero texture_dissimilarity_quarterRad texture_dissimilarity_halfRad texture_dissimilarity_threeQuaRad texture_correlation_zero texture_correlation_quarterRad texture_correlation_halfRad texture_correlation_threeQuaRad texture_ASM_zero texture_ASM_halfRad texture_ASM_threeQuaRad texture_energy_zero texture_energy_halfRad texture_energy_threeQuaRad 
## 
## Initial round 2: ..........
##  19  attributes rejected after this test:  earlySE6 earlySE11 earlySE13 dce3SE12 dce3SE15 lateSE0 lateSE1 lateSE4 lateSE7 lateSE12 lateSE14 lateSE17 Kpeak_countor UptakeRate_countor mean_F_r_i circularity max_RGH_var texture_ASM_quarterRad texture_energy_quarterRad 
## 
## Initial round 3: ..........
##  5  attributes rejected after this test:  earlySE0 dce2SE12 dce2SE18 dce3SE18 lateSE8 
## 
## Final round: ..........
##  1  attributes confirmed after this test:  edge_sharp_std 
## 
##  19  attributes rejected after this test:  earlySE1 earlySE5 earlySE18 dce2SE1 dce2SE6 dce2SE13 dce2SE14 dce2SE19 dce3SE1 dce3SE2 dce3SE9 dce3SE14 dce3SE19 lateSE2 lateSE18 lateSE19 SER_countor var_F_r_i max_RGH_mean 
## ....
##  5  attributes rejected after this test:  earlySE7 earlySE19 dce2SE9 dce3SE8 Kpeak_inside 
## ..........
##  1  attributes confirmed after this test:  lateSE9 
## 
##  2  attributes rejected after this test:  dce2SE7 alpha_countor 
## ...
##  1  attributes rejected after this test:  dce3SE16 
## ........
##  1  attributes confirmed after this test:  Slope_ini_countor 
## ...........
##  1  attributes confirmed after this test:  dce3SE6 
## ...
##  1  attributes confirmed after this test:  earlySE17 
## ....................
##  1  attributes rejected after this test:  maxCr_countor 
## .................
##  1  attributes rejected after this test:  A_countor 
## ..............
## Boruta performed 130 randomForest runs in 1.399 mins.
##         5 attributes confirmed important: earlySE17 dce3SE6
## lateSE9 Slope_ini_countor edge_sharp_std
##         163 attributes confirmed unimportant: V0 V1 V2 V3 V4 V5 V6
## V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE1
## earlySE2 earlySE3 earlySE4 earlySE5 earlySE6 earlySE7 earlySE10
## earlySE11 earlySE12 earlySE13 earlySE14 earlySE18 earlySE19
## dce2SE0 dce2SE1 dce2SE2 dce2SE3 dce2SE4 dce2SE5 dce2SE6 dce2SE7
## dce2SE8 dce2SE9 dce2SE10 dce2SE11 dce2SE12 dce2SE13 dce2SE14
## dce2SE17 dce2SE18 dce2SE19 dce3SE1 dce3SE2 dce3SE3 dce3SE5 dce3SE7
## dce3SE8 dce3SE9 dce3SE10 dce3SE11 dce3SE12 dce3SE14 dce3SE15
## dce3SE16 dce3SE17 dce3SE18 dce3SE19 lateSE0 lateSE1 lateSE2
## lateSE3 lateSE4 lateSE5 lateSE7 lateSE8 lateSE10 lateSE11 lateSE12
## lateSE13 lateSE14 lateSE15 lateSE16 lateSE17 lateSE18 lateSE19
## degreeC closenessC betweennessC no_triangles no_con_comp A_inside
## alpha_inside beta_inside iAUC1_inside Tpeak_inside Kpeak_inside
## SER_inside peakCr_inside washoutRate_inside maxVr_inside
## peakVr_inside Vr_increasingRate_inside Vr_decreasingRate_inside
## Vr_post_1_inside A_countor alpha_countor beta_countor
## iAUC1_countor Tpeak_countor Kpeak_countor SER_countor
## maxCr_countor peakCr_countor UptakeRate_countor
## washoutRate_countor maxVr_countor peakVr_countor
## Vr_increasingRate_countor Vr_decreasingRate_countor
## Vr_post_1_countor min_F_r_i mean_F_r_i var_F_r_i skew_F_r_i
## kurt_F_r_i iMax_Variance_uptake iiMin_change_Variance_uptake
## iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance circularity
## irregularity edge_sharp_mean max_RGH_mean max_RGH_mean_k
## max_RGH_var max_RGH_var_k texture_contrast_zero
## texture_contrast_quarterRad texture_contrast_halfRad
## texture_contrast_threeQuaRad texture_homogeneity_zero
## texture_homogeneity_quarterRad texture_homogeneity_halfRad
## texture_homogeneity_threeQuaRad texture_dissimilarity_zero
## texture_dissimilarity_quarterRad texture_dissimilarity_halfRad
## texture_dissimilarity_threeQuaRad texture_correlation_zero
## texture_correlation_quarterRad texture_correlation_halfRad
## texture_correlation_threeQuaRad texture_ASM_zero
## texture_ASM_quarterRad texture_ASM_halfRad texture_ASM_threeQuaRad
## texture_energy_zero texture_energy_quarterRad
## texture_energy_halfRad texture_energy_threeQuaRad
##         14 tentative attributes left: earlySE8 earlySE9 earlySE15
## earlySE16 dce2SE15 dce2SE16 dce3SE0 dce3SE4 dce3SE13 lateSE6
## Slope_ini_inside maxCr_inside UptakeRate_inside max_F_r_i
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-825.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-826.png) 

```
## D:  2 T:  5 
## Area under the curve: 0.88
## Area under the curve: 0.8
## D:  5 T:  5 
## Area under the curve: 0.955
## Area under the curve: 0.467
## D:  10 T:  5 
## Area under the curve: 0.961
## Area under the curve: 0.8
## D:  20 T:  5 
## Area under the curve: 0.946
## Area under the curve: 0.567
## D:  2 T:  10 
## Area under the curve: 0.875
## Area under the curve: 0.7
## D:  5 T:  10 
## Area under the curve: 0.956
## Area under the curve: 0.6
## D:  10 T:  10 
## Area under the curve: 0.972
## Area under the curve: 0.6
## D:  20 T:  10 
## Area under the curve: 0.946
## Area under the curve: 0.667
## D:  2 T:  30 
## Area under the curve: 0.909
## Area under the curve: 0.733
## D:  5 T:  30 
## Area under the curve: 0.983
## Area under the curve: 0.733
## D:  10 T:  30 
## Area under the curve: 0.992
## Area under the curve: 0.8
## D:  20 T:  30 
## Area under the curve: 0.974
## Area under the curve: 0.733
## D:  2 T:  60 
## Area under the curve: 0.919
## Area under the curve: 0.667
## D:  5 T:  60 
## Area under the curve: 0.975
## Area under the curve: 0.733
## D:  10 T:  60 
## Area under the curve: 0.977
## Area under the curve: 0.667
## D:  20 T:  60 
## Area under the curve: 0.982
## Area under the curve: 0.733
## D:  2 T:  100 
## Area under the curve: 0.916
## Area under the curve: 0.667
## D:  5 T:  100 
## Area under the curve: 0.978
## Area under the curve: 0.8
## D:  10 T:  100 
## Area under the curve: 0.986
## Area under the curve: 0.667
## D:  20 T:  100 
## Area under the curve: 0.972
## Area under the curve: 0.8
## D:  2 T:  250 
## Area under the curve: 0.908
## Area under the curve: 0.667
## D:  5 T:  250 
## Area under the curve: 0.976
## Area under the curve: 0.733
## D:  10 T:  250 
## Area under the curve: 0.98
## Area under the curve: 0.733
## D:  20 T:  250 
## Area under the curve: 0.972
## Area under the curve: 0.667
## D:  2 T:  500 
## Area under the curve: 0.913
## Area under the curve: 0.667
## D:  5 T:  500 
## Area under the curve: 0.972
## Area under the curve: 0.733
## D:  10 T:  500 
## Area under the curve: 0.976
## Area under the curve: 0.667
## D:  20 T:  500 
## Area under the curve: 0.973
## Area under the curve: 0.667
## D:  2 T:  750 
## Area under the curve: 0.911
## Area under the curve: 0.667
## D:  5 T:  750 
## Area under the curve: 0.972
## Area under the curve: 0.667
## D:  10 T:  750 
## Area under the curve: 0.975
## Area under the curve: 0.667
## D:  20 T:  750 
## Area under the curve: 0.975
## Area under the curve: 0.733
##     x   y acuTrain rocTrain senTrain speTrain acuTest rocTest senTest
## 1   2   5   0.8158   0.8802   0.8333   0.7857   0.750  0.8000     0.8
## 2   5   5   0.8816   0.9550   0.8750   0.8929   0.625  0.4667     0.6
## 3  10   5   0.8816   0.9613   0.9583   0.7500   0.875  0.8000     1.0
## 4  20   5   0.8026   0.9457   0.8958   0.6429   0.500  0.5667     0.8
## 5   2  10   0.8158   0.8750   0.8958   0.6786   0.625  0.7000     0.8
## 6   5  10   0.8684   0.9561   0.9583   0.7143   0.750  0.6000     1.0
## 7  10  10   0.9079   0.9717   0.9583   0.8214   0.750  0.6000     1.0
## 8  20  10   0.8947   0.9464   0.9167   0.8571   0.500  0.6667     0.6
## 9   2  30   0.8026   0.9092   0.9375   0.5714   0.625  0.7333     0.8
## 10  5  30   0.9474   0.9829   1.0000   0.8571   0.500  0.7333     0.6
## 11 10  30   0.9474   0.9918   1.0000   0.8571   0.750  0.8000     0.8
## 12 20  30   0.9342   0.9740   0.9792   0.8571   0.625  0.7333     0.8
## 13  2  60   0.8289   0.9189   0.9167   0.6786   0.750  0.6667     1.0
## 14  5  60   0.9211   0.9754   1.0000   0.7857   0.500  0.7333     0.6
## 15 10  60   0.9211   0.9769   1.0000   0.7857   0.625  0.6667     0.8
## 16 20  60   0.9079   0.9821   1.0000   0.7500   0.750  0.7333     0.8
## 17  2 100   0.8421   0.9159   0.9167   0.7143   0.750  0.6667     1.0
## 18  5 100   0.9211   0.9777   1.0000   0.7857   0.625  0.8000     0.8
## 19 10 100   0.9474   0.9859   0.9792   0.8929   0.625  0.6667     0.6
## 20 20 100   0.9342   0.9717   0.9792   0.8571   0.750  0.8000     0.8
## 21  2 250   0.8289   0.9077   0.9167   0.6786   0.625  0.6667     0.8
## 22  5 250   0.9474   0.9762   1.0000   0.8571   0.625  0.7333     0.8
## 23 10 250   0.9342   0.9799   1.0000   0.8214   0.625  0.7333     0.8
## 24 20 250   0.9342   0.9725   1.0000   0.8214   0.625  0.6667     0.8
## 25  2 500   0.8421   0.9129   0.9375   0.6786   0.625  0.6667     0.8
## 26  5 500   0.9342   0.9725   0.9792   0.8571   0.625  0.7333     0.8
## 27 10 500   0.9211   0.9762   1.0000   0.7857   0.625  0.6667     0.8
## 28 20 500   0.9342   0.9732   1.0000   0.8214   0.625  0.6667     0.8
## 29  2 750   0.8553   0.9107   0.9583   0.6786   0.625  0.6667     0.8
## 30  5 750   0.9342   0.9717   1.0000   0.8214   0.625  0.6667     0.8
## 31 10 750   0.9474   0.9754   1.0000   0.8571   0.625  0.6667     0.8
## 32 20 750   0.9211   0.9754   0.9792   0.8214   0.625  0.7333     0.8
##    speTest
## 1   0.6667
## 2   0.6667
## 3   0.6667
## 4   0.0000
## 5   0.3333
## 6   0.3333
## 7   0.3333
## 8   0.3333
## 9   0.3333
## 10  0.3333
## 11  0.6667
## 12  0.3333
## 13  0.3333
## 14  0.3333
## 15  0.3333
## 16  0.6667
## 17  0.3333
## 18  0.3333
## 19  0.6667
## 20  0.6667
## 21  0.3333
## 22  0.3333
## 23  0.3333
## 24  0.3333
## 25  0.3333
## 26  0.3333
## 27  0.3333
## 28  0.3333
## 29  0.3333
## 30  0.3333
## 31  0.3333
## 32  0.3333
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-827.png) 

```
## Initial round 1: ..........
##  113  attributes rejected after this test:  V0 V2 V4 V5 V6 V7 V8 V9 V10 V12 V13 V14 V15 V16 V17 V18 V19 earlySE2 earlySE7 earlySE10 earlySE13 dce2SE1 dce2SE2 dce2SE3 dce2SE4 dce2SE5 dce2SE7 dce2SE8 dce2SE10 dce2SE11 dce2SE12 dce2SE14 dce2SE18 dce3SE1 dce3SE3 dce3SE5 dce3SE8 dce3SE10 dce3SE11 dce3SE14 dce3SE18 lateSE0 lateSE1 lateSE4 lateSE5 lateSE8 lateSE11 lateSE13 lateSE15 lateSE16 lateSE17 lateSE18 lateSE19 degreeC closenessC betweennessC no_triangles no_con_comp A_inside alpha_inside iAUC1_inside Tpeak_inside peakCr_inside washoutRate_inside maxVr_inside peakVr_inside Vr_decreasingRate_inside iAUC1_countor Tpeak_countor Kpeak_countor peakCr_countor washoutRate_countor maxVr_countor peakVr_countor Vr_increasingRate_countor Vr_post_1_countor min_F_r_i mean_F_r_i skew_F_r_i kurt_F_r_i iMax_Variance_uptake iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance circularity irregularity edge_sharp_mean max_RGH_mean_k max_RGH_var max_RGH_var_k texture_contrast_zero texture_contrast_quarterRad texture_contrast_halfRad texture_contrast_threeQuaRad texture_homogeneity_zero texture_homogeneity_quarterRad texture_homogeneity_halfRad texture_homogeneity_threeQuaRad texture_dissimilarity_zero texture_dissimilarity_quarterRad texture_dissimilarity_halfRad texture_dissimilarity_threeQuaRad texture_correlation_zero texture_correlation_quarterRad texture_correlation_halfRad texture_correlation_threeQuaRad texture_ASM_zero texture_ASM_quarterRad texture_ASM_halfRad texture_ASM_threeQuaRad texture_energy_zero texture_energy_halfRad texture_energy_threeQuaRad 
## 
## Initial round 2: ..........
##  15  attributes rejected after this test:  V1 V3 earlySE3 dce2SE6 dce3SE12 lateSE7 lateSE10 lateSE12 SER_inside Vr_increasingRate_inside Vr_post_1_inside UptakeRate_countor iiMin_change_Variance_uptake max_RGH_mean texture_energy_quarterRad 
## 
## Initial round 3: ..........
##  16  attributes rejected after this test:  V11 earlySE0 earlySE4 earlySE6 earlySE9 earlySE12 earlySE18 dce2SE13 dce2SE19 dce3SE15 dce3SE19 lateSE2 lateSE3 Kpeak_inside maxCr_countor var_F_r_i 
## 
## Final round: ..........
##  3  attributes rejected after this test:  earlySE11 earlySE14 lateSE14 
## ....
##  7  attributes rejected after this test:  earlySE1 earlySE5 earlySE8 dce2SE0 dce2SE17 lateSE6 Vr_decreasingRate_countor 
## ....
##  1  attributes rejected after this test:  dce3SE7 
## ...
##  2  attributes confirmed after this test:  Slope_ini_countor edge_sharp_std 
## 
##  2  attributes rejected after this test:  dce2SE9 dce2SE15 
## ...
##  2  attributes rejected after this test:  dce3SE9 SER_countor 
## ..............
##  1  attributes confirmed after this test:  lateSE9 
## .....
##  1  attributes rejected after this test:  dce3SE2 
## ......
##  3  attributes confirmed after this test:  earlySE17 dce2SE16 Slope_ini_inside 
## ..
##  2  attributes rejected after this test:  earlySE15 dce3SE13 
## ...
##  1  attributes rejected after this test:  beta_inside 
## .......................................
##  1  attributes rejected after this test:  earlySE19 
## ..
##  1  attributes rejected after this test:  dce3SE17 
## .....
## Boruta performed 130 randomForest runs in 41.29 secs.
##         6 attributes confirmed important: earlySE17 dce2SE16
## lateSE9 Slope_ini_inside Slope_ini_countor edge_sharp_std
##         165 attributes confirmed unimportant: V0 V1 V2 V3 V4 V5 V6
## V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 earlySE0 earlySE1
## earlySE2 earlySE3 earlySE4 earlySE5 earlySE6 earlySE7 earlySE8
## earlySE9 earlySE10 earlySE11 earlySE12 earlySE13 earlySE14
## earlySE15 earlySE18 earlySE19 dce2SE0 dce2SE1 dce2SE2 dce2SE3
## dce2SE4 dce2SE5 dce2SE6 dce2SE7 dce2SE8 dce2SE9 dce2SE10 dce2SE11
## dce2SE12 dce2SE13 dce2SE14 dce2SE15 dce2SE17 dce2SE18 dce2SE19
## dce3SE1 dce3SE2 dce3SE3 dce3SE5 dce3SE7 dce3SE8 dce3SE9 dce3SE10
## dce3SE11 dce3SE12 dce3SE13 dce3SE14 dce3SE15 dce3SE17 dce3SE18
## dce3SE19 lateSE0 lateSE1 lateSE2 lateSE3 lateSE4 lateSE5 lateSE6
## lateSE7 lateSE8 lateSE10 lateSE11 lateSE12 lateSE13 lateSE14
## lateSE15 lateSE16 lateSE17 lateSE18 lateSE19 degreeC closenessC
## betweennessC no_triangles no_con_comp A_inside alpha_inside
## beta_inside iAUC1_inside Tpeak_inside Kpeak_inside SER_inside
## peakCr_inside washoutRate_inside maxVr_inside peakVr_inside
## Vr_increasingRate_inside Vr_decreasingRate_inside Vr_post_1_inside
## iAUC1_countor Tpeak_countor Kpeak_countor SER_countor
## maxCr_countor peakCr_countor UptakeRate_countor
## washoutRate_countor maxVr_countor peakVr_countor
## Vr_increasingRate_countor Vr_decreasingRate_countor
## Vr_post_1_countor min_F_r_i mean_F_r_i var_F_r_i skew_F_r_i
## kurt_F_r_i iMax_Variance_uptake iiMin_change_Variance_uptake
## iiiMax_Margin_Gradient k_Max_Margin_Grad ivVariance circularity
## irregularity edge_sharp_mean max_RGH_mean max_RGH_mean_k
## max_RGH_var max_RGH_var_k texture_contrast_zero
## texture_contrast_quarterRad texture_contrast_halfRad
## texture_contrast_threeQuaRad texture_homogeneity_zero
## texture_homogeneity_quarterRad texture_homogeneity_halfRad
## texture_homogeneity_threeQuaRad texture_dissimilarity_zero
## texture_dissimilarity_quarterRad texture_dissimilarity_halfRad
## texture_dissimilarity_threeQuaRad texture_correlation_zero
## texture_correlation_quarterRad texture_correlation_halfRad
## texture_correlation_threeQuaRad texture_ASM_zero
## texture_ASM_quarterRad texture_ASM_halfRad texture_ASM_threeQuaRad
## texture_energy_zero texture_energy_quarterRad
## texture_energy_halfRad texture_energy_threeQuaRad
##         11 tentative attributes left: earlySE16 dce3SE0 dce3SE4
## dce3SE6 dce3SE16 maxCr_inside UptakeRate_inside A_countor
## alpha_countor beta_countor max_F_r_i
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-828.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-829.png) 

```
## D:  2 T:  5 
## Area under the curve: 0.867
## Area under the curve: 1
## D:  5 T:  5 
## Area under the curve: 0.933
## Area under the curve: 0.867
## D:  10 T:  5 
## Area under the curve: 0.995
## Area under the curve: 1
## D:  20 T:  5 
## Area under the curve: 0.951
## Area under the curve: 0.733
## D:  2 T:  10 
## Area under the curve: 0.893
## Area under the curve: 0.867
## D:  5 T:  10 
## Area under the curve: 0.978
## Area under the curve: 0.667
## D:  10 T:  10 
## Area under the curve: 0.962
## Area under the curve: 0.733
## D:  20 T:  10 
## Area under the curve: 0.959
## Area under the curve: 1
## D:  2 T:  30 
## Area under the curve: 0.93
## Area under the curve: 0.933
## D:  5 T:  30 
## Area under the curve: 0.984
## Area under the curve: 1
## D:  10 T:  30 
## Area under the curve: 0.982
## Area under the curve: 1
## D:  20 T:  30 
## Area under the curve: 0.987
## Area under the curve: 0.8
## D:  2 T:  60 
## Area under the curve: 0.932
## Area under the curve: 1
## D:  5 T:  60 
## Area under the curve: 0.984
## Area under the curve: 0.867
## D:  10 T:  60 
## Area under the curve: 0.981
## Area under the curve: 1
## D:  20 T:  60 
## Area under the curve: 0.991
## Area under the curve: 0.933
## D:  2 T:  100 
## Area under the curve: 0.935
## Area under the curve: 0.867
## D:  5 T:  100 
## Area under the curve: 0.983
## Area under the curve: 0.867
## D:  10 T:  100 
## Area under the curve: 0.986
## Area under the curve: 0.933
## D:  20 T:  100 
## Area under the curve: 0.987
## Area under the curve: 1
## D:  2 T:  250 
## Area under the curve: 0.917
## Area under the curve: 0.933
## D:  5 T:  250 
## Area under the curve: 0.987
## Area under the curve: 0.933
## D:  10 T:  250 
## Area under the curve: 0.982
## Area under the curve: 0.867
## D:  20 T:  250 
## Area under the curve: 0.99
## Area under the curve: 1
## D:  2 T:  500 
## Area under the curve: 0.929
## Area under the curve: 0.933
## D:  5 T:  500 
## Area under the curve: 0.987
## Area under the curve: 0.933
## D:  10 T:  500 
## Area under the curve: 0.988
## Area under the curve: 0.933
## D:  20 T:  500 
## Area under the curve: 0.987
## Area under the curve: 0.933
## D:  2 T:  750 
## Area under the curve: 0.931
## Area under the curve: 0.867
## D:  5 T:  750 
## Area under the curve: 0.986
## Area under the curve: 0.933
## D:  10 T:  750 
## Area under the curve: 0.987
## Area under the curve: 1
## D:  20 T:  750 
## Area under the curve: 0.984
## Area under the curve: 0.933
##     x   y acuTrain rocTrain senTrain speTrain acuTest rocTest senTest
## 1   2   5   0.8026   0.8672   0.9167   0.6071   0.750  1.0000     0.6
## 2   5   5   0.8684   0.9334   0.8958   0.8214   0.750  0.8667     0.6
## 3  10   5   0.9474   0.9948   0.9583   0.9286   0.750  1.0000     0.6
## 4  20   5   0.8684   0.9505   0.9167   0.7857   0.625  0.7333     0.6
## 5   2  10   0.8421   0.8929   0.9583   0.6429   0.750  0.8667     0.6
## 6   5  10   0.9605   0.9781   0.9792   0.9286   0.625  0.6667     0.4
## 7  10  10   0.9211   0.9621   0.9792   0.8214   0.625  0.7333     0.6
## 8  20  10   0.9079   0.9591   0.9583   0.8214   0.625  1.0000     0.4
## 9   2  30   0.8684   0.9301   0.9583   0.7143   0.750  0.9333     0.6
## 10  5  30   0.8947   0.9844   0.9792   0.7500   0.750  1.0000     0.6
## 11 10  30   0.9211   0.9821   0.9792   0.8214   0.625  1.0000     0.4
## 12 20  30   0.9079   0.9866   0.9583   0.8214   0.875  0.8000     0.8
## 13  2  60   0.8684   0.9315   0.9583   0.7143   0.750  1.0000     0.6
## 14  5  60   0.9211   0.9836   0.9583   0.8571   0.750  0.8667     0.6
## 15 10  60   0.9342   0.9814   0.9792   0.8571   0.750  1.0000     0.6
## 16 20  60   0.9474   0.9911   1.0000   0.8571   0.625  0.9333     0.4
## 17  2 100   0.8684   0.9353   0.9375   0.7500   0.750  0.8667     0.6
## 18  5 100   0.9605   0.9829   1.0000   0.8929   0.625  0.8667     0.4
## 19 10 100   0.9342   0.9859   0.9792   0.8571   0.625  0.9333     0.4
## 20 20 100   0.9211   0.9874   0.9792   0.8214   0.625  1.0000     0.4
## 21  2 250   0.8421   0.9174   0.9375   0.6786   0.750  0.9333     0.6
## 22  5 250   0.9474   0.9874   1.0000   0.8571   0.625  0.9333     0.4
## 23 10 250   0.9474   0.9821   0.9792   0.8929   0.750  0.8667     0.6
## 24 20 250   0.9605   0.9903   1.0000   0.8929   0.625  1.0000     0.4
## 25  2 500   0.8684   0.9293   0.9375   0.7500   0.625  0.9333     0.4
## 26  5 500   0.9342   0.9874   0.9792   0.8571   0.625  0.9333     0.4
## 27 10 500   0.9342   0.9881   0.9792   0.8571   0.625  0.9333     0.4
## 28 20 500   0.9474   0.9874   0.9792   0.8929   0.625  0.9333     0.4
## 29  2 750   0.8553   0.9308   0.9375   0.7143   0.625  0.8667     0.4
## 30  5 750   0.9342   0.9859   0.9792   0.8571   0.625  0.9333     0.4
## 31 10 750   0.9474   0.9874   0.9792   0.8929   0.625  1.0000     0.4
## 32 20 750   0.9342   0.9836   0.9792   0.8571   0.625  0.9333     0.4
##    speTest
## 1   1.0000
## 2   1.0000
## 3   1.0000
## 4   0.6667
## 5   1.0000
## 6   1.0000
## 7   0.6667
## 8   1.0000
## 9   1.0000
## 10  1.0000
## 11  1.0000
## 12  1.0000
## 13  1.0000
## 14  1.0000
## 15  1.0000
## 16  1.0000
## 17  1.0000
## 18  1.0000
## 19  1.0000
## 20  1.0000
## 21  1.0000
## 22  1.0000
## 23  1.0000
## 24  1.0000
## 25  1.0000
## 26  1.0000
## 27  1.0000
## 28  1.0000
## 29  1.0000
## 30  1.0000
## 31  1.0000
## 32  1.0000
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-830.png) 

```r
accum_nonmassgrdperf = res$ensemblegrdperf[1]$grdperf
for (k in 2:cvK) {
    accum_nonmassgrdperf = accum_nonmassgrdperf + res$ensemblegrdperf[k]$grdperf
}
cvKnonmassgrdperf = accum_nonmassgrdperf/cvK
print(cvKnonmassgrdperf)
```

```
##     x   y acuTrain rocTrain senTrain speTrain acuTest rocTest senTest
## 1   2   5   0.8175   0.8896   0.9118   0.6565  0.6597  0.7261  0.7767
## 2   5   5   0.8901   0.9494   0.9224   0.8347  0.6797  0.7574  0.7567
## 3  10   5   0.8889   0.9625   0.9330   0.8131  0.7619  0.7475  0.8833
## 4  20   5   0.8822   0.9514   0.9370   0.7886  0.6519  0.7225  0.7700
## 5   2  10   0.8505   0.9098   0.9393   0.6987  0.6997  0.7492  0.8267
## 6   5  10   0.9073   0.9722   0.9599   0.8173  0.7161  0.7178  0.8100
## 7  10  10   0.9166   0.9748   0.9706   0.8242  0.6772  0.7169  0.7967
## 8  20  10   0.9075   0.9688   0.9603   0.8175  0.6411  0.7606  0.7867
## 9   2  30   0.8466   0.9219   0.9516   0.6669  0.6672  0.7650  0.8133
## 10  5  30   0.9392   0.9822   0.9937   0.8459  0.6436  0.7722  0.7767
## 11 10  30   0.9312   0.9868   0.9853   0.8389  0.7108  0.7864  0.8233
## 12 20  30   0.9446   0.9863   0.9938   0.8604  0.7008  0.7367  0.8100
## 13  2  60   0.8557   0.9268   0.9516   0.6918  0.6886  0.7669  0.8467
## 14  5  60   0.9339   0.9842   0.9938   0.8319  0.6747  0.7519  0.8067
## 15 10  60   0.9471   0.9850   0.9958   0.8640  0.6436  0.7492  0.7733
## 16 20  60   0.9366   0.9881   0.9979   0.8317  0.6922  0.7644  0.7933
## 17  2 100   0.8650   0.9307   0.9538   0.7134  0.6761  0.7689  0.8467
## 18  5 100   0.9431   0.9871   0.9979   0.8496  0.6872  0.7767  0.8067
## 19 10 100   0.9484   0.9892   0.9937   0.8712  0.6547  0.7589  0.7533
## 20 20 100   0.9378   0.9881   0.9916   0.8459  0.6797  0.8056  0.7733
## 21  2 250   0.8544   0.9279   0.9558   0.6810  0.6647  0.7544  0.8100
## 22  5 250   0.9524   0.9880   1.0000   0.8710  0.6997  0.7589  0.8267
## 23 10 250   0.9458   0.9893   0.9979   0.8566  0.6797  0.7364  0.7933
## 24 20 250   0.9458   0.9892   1.0000   0.8533  0.6772  0.7686  0.7900
## 25  2 500   0.8583   0.9296   0.9579   0.6881  0.6647  0.7644  0.7900
## 26  5 500   0.9445   0.9879   0.9958   0.8567  0.6772  0.7764  0.7900
## 27 10 500   0.9445   0.9902   0.9979   0.8533  0.6547  0.7783  0.7733
## 28 20 500   0.9471   0.9890   0.9979   0.8604  0.6772  0.7589  0.7900
## 29  2 750   0.8584   0.9292   0.9600   0.6847  0.6747  0.7467  0.8067
## 30  5 750   0.9458   0.9881   0.9979   0.8569  0.6672  0.7544  0.7733
## 31 10 750   0.9498   0.9892   0.9979   0.8676  0.6747  0.7875  0.8067
## 32 20 750   0.9458   0.9892   0.9958   0.8604  0.6772  0.7667  0.7900
##    speTest
## 1   0.4500
## 2   0.5333
## 3   0.5417
## 4   0.4417
## 5   0.4750
## 6   0.5500
## 7   0.4667
## 8   0.3833
## 9   0.4083
## 10  0.4083
## 11  0.5083
## 12  0.5083
## 13  0.4083
## 14  0.4417
## 15  0.4083
## 16  0.5083
## 17  0.3750
## 18  0.4750
## 19  0.4750
## 20  0.5083
## 21  0.4083
## 22  0.4750
## 23  0.4750
## 24  0.4750
## 25  0.4417
## 26  0.4750
## 27  0.4417
## 28  0.4750
## 29  0.4417
## 30  0.4750
## 31  0.4417
## 32  0.4750
```

```r

# plot
surface_forestperfm(cvKnonmassgrdperf)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-831.png) 

```r

# plot ROC of resamples at max perf across cvFolds
resamROC_train = data.frame()
resamROC_test = data.frame()
for (k in 1:cvK) {
    resamROC_train = rbind(resamROC_train, res$maxM[k]$maxp$trainprob)
    resamROC_test = rbind(resamROC_test, res$maxM[k]$maxp$testprob)
}
# for resamROC
ROCF_train <- plot.roc(resamROC_train$obs, resamROC_train$C, col = "#000086", 
    lty = 1)
par(new = TRUE)
ROCF_test <- plot.roc(resamROC_test$obs, resamROC_test$C, col = "#860000", lty = 2, 
    main = "ROC for non-mass max cvFolds")
print(ROCF_train$auc)
```

```
## Area under the curve: 0.938
```

```r
print(ROCF_test$auc)
```

```
## Area under the curve: 0.722
```

```r
legend("bottomright", legend = c(paste0("train: AUC=", formatC(ROCF_train$auc, 
    digits = 2, format = "f")), paste0("cv.test: AUC=", formatC(ROCF_test$auc, 
    digits = 2, format = "f"))), col = c("#000086", "#860000"), lwd = 2, lty = c(1, 
    2))
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-832.png) 

```r
# save
save.image("Z:/Cristina/MassNonmass/Section1 - ExperimentsUpToDate/experimentsRadiologypaper-revision/Tree-based-RF/ensemble-Treebased-RF/results/cvKnonmassgrdperf.RData")
```


