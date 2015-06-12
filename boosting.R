library(mboost)
ctrl <- ctree_control(maxdepth = D)
imod <- mboost(lesion_label ~ btree(earlySE11, tree_controls = ctrl)
               +btree(dce2SE8, tree_controls = ctrl) 
               +btree(earlySE0, tree_controls = ctrl)
               +btree(dce3SE8, tree_controls = ctrl),
               data = setD, family = Binomial())[500]
layout(matrix(1:4, ncol = 2))
plot(imod)
predict(imod, newdata = TestsetD, type = "class")
mean(predict(imod, newdata = TestsetD, type = "class")==TestsetD$lesion_label)


predlda = predict(ldafit, newdata = TestsetD, type = "prob")
mean(apply(predlda, 1, which.max)==unclass(TestsetD$lesion_label))

### idea: to only include high strength treees
# run testing cases
fclasspotest=list()
for (t in 1:T){
  # Calcultate posterior Probabilities on grid points
  treepred <- predict(forest[t]$tree, newdata = TestsetD) #
  treepredC=apply(treepred, 1, which.max)
  faccu = mean(treepredC==as.factor(unclass(TestsetD$lesion_label)))
  print(faccu)
  if(faccu>=0.8){
    fclasspotest <- append(fclasspotest, list(cpo = temp))
  }



# set fit parameters
fitControl <- trainControl(method = "boost",
                           number = 10,
                           repeats = 1,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

gbmGrid <-  expand.grid(.maxdepth = c(10, 20),
                        .nu =0.1,
                        .iter  = (1:3)*50)
gbmFit3 <- train(lesion_label ~ ., data = setD[,c("lesion_label",subfeat)],
                 method = "ada",
                 trControl = fitControl,
                 verbose = TRUE,
                 tuneGrid = gbmGrid,
                 ## Specify which metric to optimize
                 metric = "ROC")



control <- rpart.control(cp = -1,minsplit = 0,xval = 0,maxdepth = 3)
adaFit3 <- ada(lesion_label ~ ., data = TrainsetD[2:ncol(TrainsetD)],
               iter=100,nu=1,type="discrete",control=control)
print(adaFit3)               

##add testing data set
adaFit3 = addtest(adaFit3,TestsetD[3:ncol(TestsetD)],TestsetD[,2])
###plot gdis
plot(adaFit3,TRUE,TRUE)
pred <- predict(adaFit3, TestsetD[3:ncol(TestsetD)])
groundT <- TestsetD$lesion_label
confusionMatrix(pred, groundT)
summary(adaFit3, n.iter = 98)