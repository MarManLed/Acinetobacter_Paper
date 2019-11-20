### Functions and packages ###

library(MALDIquant)
library(MALDIquantForeign)
library(dplyr)
library(purrr)
library(ggplot2)
library(DT)
library(rpart)
library(partykit)
library("pvclust")
library(randomForest)
library(xgboost)
library("cluster")
library("factoextra")
library("sda")
library(crossval)
library(mlr)
library(caret)
library("binda")
library(plotROC)
library(broom)
library(MASS)
library(statsr)
library(tidyr)

Avg.generation<-function(SO,transform.method,smoothmethod,halfwind,
                         labells,removemethod,iteration,calibrateMEto,
                         Snr,Tolerance, warpingmeto,avemeto,title){
    avgSpectra<-SO %>%
        trim() %>%
        transformIntensity( method=transform.method) %>%
        smoothIntensity(method=smoothmethod,
                        halfWindowSize=halfwind)  %>%
        removeBaseline( method=removemethod, iterations=iteration)%>%
        calibrateIntensity( method=calibrateMEto) %>%
        alignSpectra( halfWindowSize=halfwind,SNR=Snr,
                      tolerance=Tolerance,warpingMethod=warpingmeto)%>%
        averageMassSpectra( labels=labells, method=avemeto)
    
    save(avgSpectra, file=paste("avgSpectra", title,"rda", sep="."))
}

Clusters.Plots <- function (feature.matrix,
                            metodologia,
                            spot,
                            sensi,
                            other,
                            n.cluster,
                            title) {
    # Heat Map #
    rownames(feature.matrix) <- sensi
    dist.eucl <-dist(feature.matrix, method = metodologia)
    heat<-fviz_dist(dist.eucl)+
        ggtitle(paste("Heat Map",title,sep=" "))
    # K means clustering #
    rownames(feature.matrix) <- paste(spot,sensi,other,sep=".")
    Colmun.const<-which(apply(feature.matrix, 2, var)==0)
    feature.matrix.2<-feature.matrix[ , apply(feature.matrix, 2, var) != 0]
    set.seed(345)
    km.res <-kmeans(feature.matrix.2, n.cluster, nstart = 25)
    kmeansPlot<-fviz_cluster(km.res, 
                             data = feature.matrix.2, 
                             main = paste("K-means clustering", title, n.cluster, "cluster",sep=" "),
                             palette = "rdb",
                             ellipse.type = "euclid", 
                             star.plot = TRUE,
                             repel = TRUE, # 
                             ggtheme =theme_minimal())
    # Hierarchical k-means clustering #
    set.seed(345)
    km1.res <-hkmeans(feature.matrix, n.cluster,
                      hc.metric= metodologia,
                      hc.method="ward.D2",
                      iter.max=10)
    arbol<-fviz_dend(km1.res,
                     repel = TRUE,
                     type= "phylogenic",
                     cex = 0.7,
                     phylo_layout = "layout.auto")+
        ggtitle(paste("Cluster Dendrogram",title,n.cluster,"cluster",sep=" "))
    print(heat)
    print(kmeansPlot) 
    print(arbol)
    print(paste("Columna", Colmun.const, "con varianza cero descartada para k-means",sep=" "))
}

Binary.analisis<-function(feature.matrix,
                          Peaks,
                          spot,
                          sensi,
                          other,
                          title) {
    Ytrain1<-sensi
    categorias<-table(Ytrain1)
    total.de.NA<-sum(is.na(intensityMatrix(Peaks))) 
    promedio.de.NA<-mean(apply(is.na(intensityMatrix(Peaks)), 1, sum)) 
    
    mz.full <- as.double(colnames(feature.matrix))
    mz <- round( mz.full )
    picos.ducplicados<-any(duplicated(mz))
    
    thr <- optimizeThreshold(feature.matrix, Ytrain1)
    Xtrain.b <- dichotomize(feature.matrix, thr)
    colnames(Xtrain.b) <-mz
    rownames(Xtrain.b) <- paste(spot,other,sensi, sep=".")
    Chequeo.dicotomico<-is.binaryMatrix(Xtrain.b)
    
    Xtrain.b.naive <- ifelse(is.na(intensityMatrix(Peaks)), as.integer(0), as.integer(1))
    rownames(Xtrain.b.naive) <-paste(spot,other,sensi, sep=".")
    colnames(Xtrain.b.naive) <-mz
    Chequeo.binario<-is.binaryMatrix(Xtrain.b.naive)
    
    Feature.dicho3<-Xtrain.b
    Feature.bina3<-Xtrain.b.naive
    
    save(Feature.dicho3, Feature.bina3,
         file=paste("Matrices", title, "rda", sep="."))
    
    print(categorias)
    print(paste("total de NA",total.de.NA,sep=" "))
    print(paste("promedio de NA",promedio.de.NA,sep=" "))
    print(paste("Algun pico duplicado",picos.ducplicados,sep=" "))
    print(paste("Dicho es binario",Chequeo.dicotomico,sep=" "))
    print(paste("Bina es binario",Chequeo.binario,sep=" "))
    
    # Hierchical clustering with boots
    pv.d <- pvclust(t(Xtrain.b),
                    method.hclust="ward.D2",
                    method.dist="binary")
    pv.b <- pvclust(t(Xtrain.b.naive),
                    method.hclust="ward.D2",
                    method.dist="binary")
    
    plot(pv.b, main= "Binary data", print.num=FALSE)
    plot(pv.d, main= "Dichotomized data",print.num=FALSE)
    
}

Sampling.Random<-function(feature.Matrix,
                          title,
                          sensi,
                          seed,
                          train.p,
                          test.p){
    
    rownames(feature.Matrix)<-sensi
    set.seed(seed)
    ind <- sample(2, nrow(feature.Matrix), replace = TRUE, prob = c(train.p, test.p))
    Train <- feature.Matrix[ind==1, ] 
    Test <- feature.Matrix[ind==2, ]
    Ytrain11<-rownames(Train)
    Ytest1<-rownames(Test)
    Ytrain11<-as.factor(Ytrain11)
    Ytest1<-as.factor(Ytest1)
    print(table(Ytrain11))
    print(table(Ytest1))
    save(Train, Test, Ytrain11, Ytest1,
         file=paste("Train.test", title, "rda", sep="."))
}

Predicciones.binda<-function(feature.Matrix,
                             rowNames.featu,
                             TTrain,
                             YYtrain,
                             TTEst,
                             YYtest,
                             Negative,
                             title){
    # Binda ranking
    br = binda.ranking(feature.Matrix, rowNames.featu)
    Plot.binda<-plot(br, top=10, arrow.col="blue", zeroaxis.col="black", 
                     ylab="Peaks (m/z)",
                     main=paste("10 highly differentiated peaks", title, sep=" "))
    print(Plot.binda)
    print(br)
    # binda prediction hole spectrum
    Binda.out.hole <- binda(TTrain, YYtrain, verbose=FALSE)
    ynew.hole <- predict.binda(Binda.out.hole, TTEst, verbose=FALSE)$class
    cm.hole <- crossval::confusionMatrix(YYtest, ynew.hole, negative= Negative) 
    print("binda prediction hole spectrum")
    print(cm.hole)
    print(diagnosticErrors(cm.hole))
    c1<-tidy(diagnosticErrors(cm.hole))
    names(c1)<-c("Parametro","Valor.hole")
    
    # binda prediction 20 best peaks
    selPeaks20<-br[,"idx"][1:20]
    binda.out20 = binda(TTrain[, selPeaks20, drop=FALSE], YYtrain, verbose=FALSE)
    ynew20 = predict.binda(binda.out20, TTEst[, selPeaks20, drop=FALSE], verbose=FALSE)$class
    cm.20 = crossval::confusionMatrix(YYtest, ynew20, negative= Negative) 
    print("binda prediction 20 best peaks")
    print(cm.20)
    print(diagnosticErrors(cm.20))
    c2<-tidy(diagnosticErrors(cm.20))
    names(c2)<-c("Parametro","Valor.20")
    
    # binda prediction 15 best peaks
    selPeaks15<-br[,"idx"][1:15]
    binda.out15 = binda(TTrain[, selPeaks15, drop=FALSE], YYtrain, verbose=FALSE)
    ynew15 = predict.binda(binda.out15, TTEst[, selPeaks15, drop=FALSE], verbose=FALSE)$class
    cm.15 = crossval::confusionMatrix(YYtest, ynew15, negative= Negative) 
    print("binda prediction 15 best peaks")
    print(cm.15)
    print(diagnosticErrors(cm.15))
    c3<-tidy(diagnosticErrors(cm.15))
    names(c3)<-c("Parametro","Valor.15")
   
     # binda prediction 10 best peaks
    
    selPeaks10<-br[,"idx"][1:10]
    binda.out10 = binda(TTrain[, selPeaks10, drop=FALSE], YYtrain, verbose=FALSE)
    ynew10 = predict.binda(binda.out10, TTEst[, selPeaks10, drop=FALSE], verbose=FALSE)$class
    cm.10 = crossval::confusionMatrix(YYtest, ynew10, negative= Negative) 
    print("binda prediction 10 best peaks")
    print(cm.10)
    print(diagnosticErrors(cm.10))
    c4<-tidy(diagnosticErrors(cm.10))
    names(c4)<-c("Parametro","Valor.10")
    
    # binda prediction 5 best peaks
    
    selPeaks5<-br[,"idx"][1:5]
    binda.out5 = binda(TTrain[, selPeaks5, drop=FALSE], YYtrain, verbose=FALSE)
    ynew5 = predict.binda(binda.out5, TTEst[, selPeaks5, drop=FALSE], verbose=FALSE)$class
    cm.5 = crossval::confusionMatrix(YYtest, ynew5, negative= Negative) 
    print("binda prediction 5 best peaks")
    print(cm.5)
    print(diagnosticErrors(cm.5))
    
    c5<-tidy(diagnosticErrors(cm.5))
    names(c5)<-c("Parametro","Valor.5")
    Resumen.binda<-c1 %>%
        left_join(c2, by="Parametro")%>%
        left_join(c3, by="Parametro")%>%
        left_join(c4, by="Parametro")%>%
        left_join(c5, by="Parametro")
    
    save(Resumen.binda,ynew.hole,ynew15,ynew5,ynew10,ynew20,
         file=paste("Binda", title, "rda", sep="."))
    print(selPeaks10)
    
}

Crossvalidacion<-function(feature.Matrix,
                          rowNames.featu,
                          TTrain,
                          YYtrain,
                          TTEst,
                          YYtest,
                          Negative,
                          seed,
                          k,
                          b,
                          title)
{
    
    predfun1 <- function(Xtrain, Ytrain, Xtest, Ytest, selPeaks)
    {
        binda.out = binda(Xtrain[, selPeaks, drop=FALSE], Ytrain, verbose=FALSE)
        ynew = predict.binda(binda.out, Xtest[, selPeaks, drop=FALSE], verbose=FALSE)$class 
        cm = crossval::confusionMatrix(Ytest, ynew, negative=Negative) 
        
        return(cm)  
    }
    br = binda.ranking(feature.Matrix, rowNames.featu)
    
    # 20 peaks
    ourPeaks.20 = br[,"idx"][1:20]
    cvp.train.20 <- crossval::crossval(predfun1, TTrain, YYtrain, K=k, B=b, selPeaks=ourPeaks.20,verbose=FALSE)
    print("Cv binda prediction 20 best peaks on train")
    print(cvp.train.20$stat)
    print(diagnosticErrors(cvp.train.20$stat))
    
    c1<-tidy(diagnosticErrors(cvp.train.20$stat))
    names(c1)<-c("Parametro","Valor.20.train")
    
    cvp.test.20 <- crossval::crossval(predfun1, TTEst, YYtest, K=k, B=b, selPeaks=ourPeaks.20,verbose=FALSE)
    print("Cv binda prediction 20 best peaks on test")
    print(cvp.test.20$stat)
    print(diagnosticErrors(cvp.test.20$stat))
    
    c2<-tidy(diagnosticErrors(cvp.test.20$stat))
    names(c2)<-c("Parametro","Valor.20.test")
    Valor20Test<-c2$Valor.20.test
    
    # 15 peaks
    ourPeaks.15 = br[,"idx"][1:15]
    cvp.train.15 <- crossval::crossval(predfun1, TTrain, YYtrain, K=k, B=b, selPeaks=ourPeaks.15,verbose=FALSE)
    print("Cv binda prediction 15 best peaks on train")
    print(cvp.train.15$stat)
    print(diagnosticErrors(cvp.train.15$stat))
    
    c3<-tidy(diagnosticErrors(cvp.train.15$stat))
    names(c3)<-c("Parametro","Valor.15.train")
    Valor15Train<-c3$Valor.15.train
    
    cvp.test.15 <- crossval::crossval(predfun1, TTEst, YYtest, K=k, B=b, selPeaks=ourPeaks.15,verbose=FALSE)
    print("Cv binda prediction 15 best peaks on test")
    print(cvp.test.15$stat)
    print(diagnosticErrors(cvp.test.15$stat))
    
    c4<-tidy(diagnosticErrors(cvp.test.15$stat))
    names(c4)<-c("Parametro","Valor.15.test")
    Valor15Test<-c4$Valor.15.test
    
    # 10 peaks
    ourPeaks.10 = br[,"idx"][1:10]
    cvp.train.10 <- crossval::crossval(predfun1, TTrain, YYtrain, K=k, B=b, selPeaks=ourPeaks.10,verbose=FALSE)
    print("Cv binda prediction 10 best peaks on train")
    print(cvp.train.10$stat)
    print(diagnosticErrors(cvp.train.10$stat))
    
    c5<-tidy(diagnosticErrors(cvp.train.10$stat))
    names(c5)<-c("Parametro","Valor.10.train")
    Valor10Train<-c5$Valor.10.train
    
    cvp.test.10 <- crossval::crossval(predfun1, TTEst, YYtest, K=k, B=b, selPeaks=ourPeaks.10,verbose=FALSE)
    print("Cv binda prediction 10 best peaks on test")
    print(cvp.test.10$stat)
    print(diagnosticErrors(cvp.test.10$stat))
    
    c6<-tidy(diagnosticErrors(cvp.test.10$stat))
    names(c6)<-c("Parametro","Valor.10.test")
    Valor10Test<-c6$Valor.10.test
    
    # 5 peaks
    ourPeaks.5 = br[,"idx"][1:5]
    cvp.train.5<- crossval::crossval(predfun1, TTrain, YYtrain, K=k, B=b, selPeaks=ourPeaks.5,verbose=FALSE)
    print("Cv binda prediction 5 best peaks on train")
    print(cvp.train.5$stat)
    print(diagnosticErrors(cvp.train.5$stat))
    
    c7<-tidy(diagnosticErrors(cvp.train.5$stat))
    names(c7)<-c("Parametro","Valor.5.train")
    Valor5Train<-c7$Valor.5.train
    
    cvp.test.5 <- crossval::crossval(predfun1, TTEst, YYtest, K=k, B=b, selPeaks=ourPeaks.5,verbose=FALSE)
    print("Cv binda prediction 5 best peaks on test")
    print(cvp.test.5$stat)
    print(diagnosticErrors(cvp.test.5$stat))
    
    c8<-tidy(diagnosticErrors(cvp.test.5$stat))
    names(c8)<-c("Parametro","Valor.5.test")
    Valor5Test<-c8$Valor.5.test
    
    Resumen<-c1 %>%
        left_join(c2, by="Parametro")%>%
        left_join(c3, by="Parametro")%>%
        left_join(c4, by="Parametro")%>%
        left_join(c5, by="Parametro")%>%
        left_join(c6, by="Parametro")%>%
        left_join(c7, by="Parametro")%>%
        left_join(c8, by="Parametro")
    
    save(Resumen,
         file=paste("Cv", title, "rda", sep="."))
    
}


replace_missings <- function(x, replacement) {
    is_miss <- is.na(x)
    x[is_miss] <- replacement
    # Rewrite to use message()
    message(sum(is_miss), " missings replaced by the value ",  replacement)
    x
}

Other.ML.Arboles.RandomF<-function(feature.Matrix,
                                   rowName,
                                   seed,
                                   prob.train,
                                   prob.test,
                                   Negative,
                                   colnum,
                                   title){
    
    mz.full <- as.double(colnames(feature.Matrix))
    mz <- round( mz.full )
    colnames(feature.Matrix) <-mz
    
    feature.Matrix.Df<-as.data.frame(feature.Matrix)%>% 
        mutate(sensi=rowName)
    set.seed(seed)
    
    ind <- sample(2, nrow(feature.Matrix.Df), replace=TRUE, prob=c(prob.train,prob.test))
    Train <- feature.Matrix.Df[ind==1, ] 
    Test <- feature.Matrix.Df[ind==2, ] 
    
    # Arboles 
    tree.pros<- rpart(sensi~. , data=Train)
    print(tree.pros$cptable)
    arbol.plot<-plotcp(tree.pros)
    print(arbol.plot)
    
    cp<- min(tree.pros$cptable[1, ])
    prune.tree.biop<-prune(tree.pros, cp=cp)
    arbol.plot2<-plot(as.party(prune.tree.biop))
    print(arbol.plot2)
    
    party.pros.test<- predict(prune.tree.biop, newdata= Test, type="class")
    
    ca<-crossval::confusionMatrix(party.pros.test, Test$sensi, negative=Negative)
    caa<-diagnosticErrors(ca)
    print("Arboles de decisión")
    print(caa)
    
    # Random forest
    
    set.seed(71)
    rf.pros <- randomForest(Train[,-colnum],Train$sensi, importance=TRUE,
                            proximity=TRUE)
    
    min<-as.integer(which.min(rf.pros$err.rate[,1]))
    set.seed(71)
    rf.pros2 <- randomForest(Train[,-colnum],Train$sensi, ntree= min)
    rf.biop.test<- predict(rf.pros2, newdata=Test[,-colnum], type="response")
    
    cr<-crossval::confusionMatrix(rf.biop.test, Test$sensi, negative=Negative)
    crr<-diagnosticErrors(cr)
    print("Random Forest")
    print(crr)
    
    RandomPlot<-varImpPlot(rf.pros2, main="Random forest best peaks")
    
    random.Forest<-tidy(print(RandomPlot))%>%
        arrange(desc(MeanDecreaseGini))
    
    save(random.Forest,
         file=paste("random.Forest", title, "rda", sep="."))
    
    
}

Analisis.discriminante<-function(feature.Matrix,
                                 Ytain,
                                 Diagonal,
                                 Top,
                                 title) {
    
    ddar <- sda.ranking(Xtrain=feature.Matrix, L=Ytain, fdr=FALSE,
                        diagonal=Diagonal, verbose=FALSE)
    d1<-plot(ddar, top=Top, arrow.col="red", zeroaxis.col="black", ylab="Peaks (m/z)",
             main=paste("Differentially Expressed Peaks",title,sep=" "))
    top <- ddar[1:Top, "idx"]
    distanceMatrixTop <- dist(feature.Matrix[, top],
                              method="euclidean")
    hClustTop <- hclust(distanceMatrixTop, method="ward.D2")
    d2<-plot(hClustTop, hang=-1, main=paste(title,"top",Top,sep="."))
    
    
    print(d1)
    print(d2)
}

Predicciones.discriminantes<-function(feature.Matrix,Ytain,TraiN,
                                      YTraIn,TEsT,YTesT,Diagonal,
                                      TOp,seed,title,Negative){
    
    # Hole
    
    dda.fit <- sda(TraiN, YTraIn, diagonal=Diagonal, verbose=FALSE)
    ynew.hole.DA <- predict(dda.fit, TEsT, verbose=FALSE)$class
    cm1<-crossval::confusionMatrix(YTesT, ynew.hole.DA, negative=Negative)
    print(paste("Results hole spectrum",title,sep=" "))
    print(diagnosticErrors(cm1))
    
    # 15
    
    ddar <- sda.ranking(Xtrain=feature.Matrix, L=Ytain, fdr=FALSE,
                        diagonal=Diagonal, verbose=FALSE)
    selPeaks2n<-ddar[1:15, "idx"]
    dda.out2n = sda(TraiN[, selPeaks2n, drop=FALSE], YTraIn, diagonal=Diagonal,verbose=FALSE)
    ynew.15.DA = predict(dda.out2n, TEsT[, selPeaks2n, drop=FALSE], verbose=FALSE)$class
    cm4<-crossval::confusionMatrix(YTesT, ynew.15.DA, negative=Negative)
    print(paste("Results 15 best spectrum",title,sep=" "))
    print(diagnosticErrors(cm4))
    
    # 10
    
    ddar <- sda.ranking(Xtrain=feature.Matrix, L=Ytain, fdr=FALSE,
                        diagonal=Diagonal, verbose=FALSE)
    selPeaks2n<-ddar[1:10, "idx"]
    dda.out2n = sda(TraiN[, selPeaks2n, drop=FALSE], YTraIn, diagonal=Diagonal,verbose=FALSE)
    ynew.10.DA = predict(dda.out2n, TEsT[, selPeaks2n, drop=FALSE], verbose=FALSE)$class
    cm2<-crossval::confusionMatrix(YTesT, ynew.10.DA, negative=Negative)
    print(paste("Results 10 best spectrum",title,sep=" "))
    print(diagnosticErrors(cm2))
    
    # 5
    selPeaks3n<-ddar[1:5, "idx"]
    dda.out3n = sda(TraiN[, selPeaks3n, drop=FALSE], YTraIn,  diagonal=Diagonal, verbose=FALSE)
    ynew.5.DA = predict(dda.out3n, TEsT[, selPeaks3n, drop=FALSE], verbose=FALSE)$class
    cm3<-crossval::confusionMatrix(YTesT, ynew.5.DA, negative=Negative)
    
    print(paste("Results 5 best spectrum",title,sep=" "))
    print(diagnosticErrors(cm3))
    
    save(ynew.hole.DA,ynew.15.DA,ynew.10.DA,ynew.5.DA,
         file=paste("Predicciones",title,"rda",sep="."))
    
    print(ddar)
    
    predfun.dda <- function(Xtrain, Ytrain, Xtest, Ytest, negative) {
        dda.fit <- sda(Xtrain, Ytrain, diagonal=TRUE, verbose=FALSE)
        ynew <- predict(dda.fit, Xtest, verbose=FALSE)$class
        return(crossval::confusionMatrix(Ytest, ynew, negative=negative))
    }
    set.seed(seed)
    
    # top
    top <- ddar[1:TOp, "idx"]
    cv.out.dda.train <- crossval::crossval(predfun.dda, 
                                           X=TraiN[, top],
                                           Y=YTraIn,
                                           K=5, B=20,
                                           negative=Negative,
                                           verbose=FALSE)
    print(paste("Train results",TOp,title,sep=" "))
    print(diagnosticErrors(cv.out.dda.train$stat))
    
    cv.out.dda.test <- crossval::crossval(predfun.dda, 
                                          X=TEsT[, top],
                                          Y=YTesT,
                                          K=5, B=20,
                                          negative=Negative,
                                          verbose=FALSE)
    
    print(paste("Test results",TOp,title,sep=" "))
    print(diagnosticErrors(cv.out.dda.test$stat))
}

Extreme.g.boost<-function(feature.Matrix,rowName,colnum,
                          seed, prob.train, prob.test,
                          Negative,nRounds,colSample_bytree,
                          min_Child_weight,Eta,Gamma,
                          subSample,max_Depth){
    
    mz.full <- as.double(colnames(feature.Matrix))
    mz <- round( mz.full )
    colnames(feature.Matrix) <-mz
    
    feature.Matrix.d<-data.frame(feature.Matrix)
    feature.Matrix.d.l<-map(feature.Matrix.d, as.numeric)
    feature.Matrix.Df<-data.frame(feature.Matrix.d.l)
    
    feature.Matrix.Df <-feature.Matrix.Df %>% 
        mutate(sensi=rowName)
    
    set.seed(seed)
    
    ind <- sample(2, nrow(feature.Matrix.Df), replace=TRUE, prob=c(prob.train,prob.test))
    Train <- feature.Matrix.Df[ind==1, ] 
    Test <- feature.Matrix.Df[ind==2, ] 
    
    grid = expand.grid(
        nrounds = nRounds, # maximun number of iteration (number of trees in the final model)
        colsample_bytree = colSample_bytree, # number of features as ratio
        min_child_weight = min_Child_weight,
        eta = Eta, #0.3 is default, Learning rate, contribuiton of each tree to solution
        gamma = Gamma, # minimum loss reduction requeired to make another leaf partition
        subsample = subSample, # ratio of data observations 
        max_depth = max_Depth) # maximun depth of the individual trees
    
    cntrl = trainControl(
        method = "cv",
        number = 5,
        verboseIter = FALSE,
        returnData = FALSE,
        returnResamp = "final")
    
    set.seed(1)
    
    train.xgb = train(
        x = Train[, -colnum],
        y = factor(Train[, colnum]),
        trControl = cntrl,
        tuneGrid = grid,
        method = "xgbTree")
    
    train.xgb$results
    
    param <- list( objective = "binary:logistic", 
                   booster = "gbtree",
                   eval_metric = "error",
                   eta = Eta, 
                   max_depth = max_Depth, 
                   subsample = subSample,
                   colsample_bytree = colSample_bytree,
                   gamma = Gamma)
    
    xtrain <- as.matrix(Train[, -colnum])
    ytrain <- ifelse(Train$sensi == Negative, 0, 1)
    train.mat <- xgb.DMatrix(data = xtrain, label = ytrain)
    
    set.seed(1)
    xgb.fit <- xgb.train(params = param, data = train.mat, nrounds =nRounds)
    impMatrix <- xgb.importance(feature_names = dimnames(xtrain)[[2]], 
                                model = xgb.fit)   
    
    save(impMatrix,
         file=paste("Xg.boost", "rda", sep="."))
    
    xgb.plot.importance(impMatrix, main = "Gain by Feature")
    
    
    pred.train <- predict(xgb.fit, xtrain)
    threshold<-InformationValue::optimalCutoff(ytrain, pred.train)
    
    y.train <- ifelse(Train$sensi == Negative, 0, 1)    
    pred.train.b <- ifelse(pred.train > 0.5, 1, 0)    
    
    cm2<-crossval::confusionMatrix(y.train, pred.train.b,negative=0)
    print("Xgboost predictions train")
    print(diagnosticErrors(cm2))
    
    M_1.test_X <- as.matrix(Test[, -colnum])
    xgb.prima.test <- predict(xgb.fit, M_1.test_X)
    y.test <- ifelse(Test$sensi == Negative, 0, 1)    
    xgb.prima.test.b <- ifelse(xgb.prima.test > threshold, 1, 0)    
    
    cm1<-crossval::confusionMatrix(y.test, xgb.prima.test.b,negative=0)
    print("Xgboost predictions test")
    print(diagnosticErrors(cm1))
    print(impMatrix)
}

ROC.Curve.1<-function(yNEW,Ytest,Negative){
    TT<-yNEW %>%
        data.frame(Ytest) %>%
        mutate(Ytest = ifelse(Ytest==Negative,0,1),
               yNEW = ifelse(yNEW==Negative,0,1))
    
    basicplot<-ggplot(TT, aes(d = Ytest, m = yNEW)) + 
        geom_roc() + 
        geom_rocci(sig.level = 0.01) + 
        style_roc()
    direct_label(basicplot) 
    print(basicplot)
    
}

ROC.Curve.3<-function(yNEW1,d1.1,d2.1,d3.1,yNEW2, yNEW3,Ytest,title,Negative){
    
    TT<- data.frame(Ytest,yNEW2,yNEW3,yNEW1) %>%
        mutate(Ytest.1 = ifelse(Ytest==Negative,0,1),
               d1 = ifelse(yNEW1==Negative,0,1),
               d2 = ifelse(yNEW2==Negative,0,1),
               d3 = ifelse(yNEW3==Negative,0,1))
    
    TTtest <- melt_roc(TT, "Ytest.1", c("d1", "d2", "d3"))
    
    TTtest<-TTtest%>%
        mutate(
            name= case_when(
                name == "d1" ~  d1.1,
                name == "d2" ~  d2.1,
                name == "d3" ~  d3.1,
                TRUE ~ as.character(name)
            ))
    
    
    basicplot<-ggplot(TTtest, aes(d = D, m = M, color=name)) + 
        geom_roc() 
    basicplot + 
        geom_rocci(sig.level = 0.01) + 
        style_roc() +
        ggtitle(title) 
}

#### Reading bruker files #####

so4<-importBrukerFlex("~/ACINETO AGUA/Acineto.ID.agua", verbose=FALSE)

Names<-c(rep("3",5),rep("4",4),rep("11",4),rep("15",4),rep("16",4),
         rep("19",4),rep("20",4),rep("23",4),rep("24",4),rep("25",4),
         rep("26",4),rep("27",4), rep("28",4),rep("31",4),rep("38",4),
         rep("39",4),rep("40",4),rep("55",4),rep("57",4),rep("58",4))

Names.arr<-array(Names)
for (h in 1:length(so4)) {
    metaData(so4[[h]])$spot<-Names.arr[[h]]
}
table(lengths(so4))
any(sapply(so4, MALDIquant::isEmpty))
all(sapply(so4, isRegular))
so5<-so4
spot.factor <- factor(sapply(so4,function(x)metaData(x)$spot))
Avg.generation(SO=so4,
               transform.method="sqrt",
               labells=spot.factor,
               smoothmethod="SavitzkyGolay",
               halfwind=30,
               removemethod="SNIP",
               iteration=75,
               calibrateMEto="TIC",
               Snr=3,
               Tolerance=0.002,
               warpingmeto="quadratic",
               avemeto="sum",
               title="Acineto.agua.1")
load("avgSpectra.Acineto.agua.1.rda")
avgSpectra.Acineto.MBT.agua.1<-avgSpectra

so4<-importBrukerFlex("~/ACINETO AGUA/Acineto.ID.agua.1", verbose=FALSE)

Names<-c(rep("59",4),rep("60",4),rep("61",4),rep("62",4),rep("63",4),rep("70",5))
Names.arr<-array(Names)
for (h in 1:length(so4)) {
    metaData(so4[[h]])$spot<-Names.arr[[h]]
}
table(lengths(so4))
any(sapply(so4, MALDIquant::isEmpty))
all(sapply(so4, isRegular))

spot.factor <- factor(sapply(so4,function(x)metaData(x)$spot))
Avg.generation(SO=so4,
               transform.method="sqrt",
               labells=spot.factor,
               smoothmethod="SavitzkyGolay",
               halfwind=30,
               removemethod="SNIP",
               iteration=75,
               calibrateMEto="TIC",
               Snr=3,
               Tolerance=0.002,
               warpingmeto="quadratic",
               avemeto="sum",
               title="Acineto.agua.1.1")
load("avgSpectra.Acineto.agua.1.1.rda")

so4<-importBrukerFlex("~/ACINETO AGUA/Acineto.ID.agua.3", verbose=FALSE)
Names<-c(rep("9",4),rep("10",4),rep("11",4),rep("12",4),rep("14",4),
         rep("21",4),rep("22",4),rep("36",4),rep("37",4),rep("52",4),
         rep("56",4),rep("66",4), rep("68",4),rep("69",4))
Names.arr<-array(Names)
for (h in 1:length(so4)) {
    metaData(so4[[h]])$spot<-Names.arr[[h]]
}
table(lengths(so4))
any(sapply(so4, MALDIquant::isEmpty))
all(sapply(so4, isRegular))
so5<-so4
spot.factor <- factor(sapply(so4,function(x)metaData(x)$spot))
Avg.generation(SO=so4,
               transform.method="sqrt",
               labells=spot.factor,
               smoothmethod="SavitzkyGolay",
               halfwind=30,
               removemethod="SNIP",
               iteration=75,
               calibrateMEto="TIC",
               Snr=3,
               Tolerance=0.002,
               warpingmeto="quadratic",
               avemeto="sum",
               title="Acineto.agua.3")
load("avgSpectra.Acineto.agua.3.rda")
avgSpectra.Acineto.agua.MBT.3<-avgSpectra

so4<-importBrukerFlex("~/ACINETO AGUA/Acineto.agua.ID.4", verbose=FALSE)
Names<-c(rep("1",6),rep("7",6),rep("8",6),rep("12",6),rep("17",6),
         rep("34",6),rep("35",6),rep("43",6),rep("44",6),rep("45",6),
         rep("54",6),rep("64",6), rep("65",6))
Names.arr<-array(Names)
for (h in 1:length(so4)) {
    metaData(so4[[h]])$spot<-Names.arr[[h]]
}
table(lengths(so4))
any(sapply(so4, MALDIquant::isEmpty))
all(sapply(so4, isRegular))
so5<-so4
spot.factor <- factor(sapply(so4,function(x)metaData(x)$spot))
Avg.generation(SO=so4,
               transform.method="sqrt",
               labells=spot.factor,
               smoothmethod="SavitzkyGolay",
               halfwind=30,
               removemethod="SNIP",
               iteration=75,
               calibrateMEto="TIC",
               Snr=3,
               Tolerance=0.002,
               warpingmeto="quadratic",
               avemeto="sum",
               title="Acineto.agua.4")

load("avgSpectra.Acineto.agua.4.rda")
avgSpectra.Acineto.agua.MBT.4<-avgSpectra

AvgSpectra.Acineto.agua.Total.MBT<-c(avgSpectra.Acineto.MBT.agua.1,
                                     avgSpectra.Acineto.agua.MBT.1.1,
                                     avgSpectra.Acineto.agua.MBT.3,
                                     avgSpectra.Acineto.agua.MBT.4)

spot.factor.a <- factor(sapply(AvgSpectra.Acineto.agua.Total.MBT,function(x)metaData(x)$spot))
AvgSpectra.Acineto.agua.Total.MBT <-averageMassSpectra(AvgSpectra.Acineto.agua.Total.MBT, 
                                                       labels=spot.factor.a, method="sum")

spot.factor <- factor(sapply(AvgSpectra.Acineto.agua.Total.MBT,function(x)metaData(x)$spot))
spot.factor.C<-as.character(spot.factor)
spot.factor.I<-as.integer(spot.factor.C)

AvgSpectra.Acineto.agua.Total.MBT.df<-data.frame(spot.factor.I)
names(AvgSpectra.Acineto.agua.Total.MBT.df)<-"numero.MT"

#### Base.Acineto.tidy ####

Base.Acineto<-read.csv2("/Volumes/TOSHIBA EXT/Mac/trabajos MALDI/Resistencia.ATB/Acinetobacter/Análisis hernan/MBT/csv/MARTIN.csv", 
                        stringsAsFactors = FALSE)

Base.Acineto.tidy<-Base.Acineto
Base.Acineto.tidy$R.colistina<-factor(Base.Acineto.tidy$R.colistina)
Base.Acineto.tidy$ST<-factor(Base.Acineto.tidy$ST)
Base.Acineto.tidy$Mutacion<-factor(Base.Acineto.tidy$Mutacion)
Base.Acineto.tidy$COL<-factor(Base.Acineto.tidy$COL)
Base.Acineto.tidy <- Base.Acineto.tidy %>%
    mutate(ST = case_when(
        ST == "" ~ "ND",
        ST == "ind" ~ "ND",
        TRUE ~ as.character(ST)))
Base.Acineto.tidy$ST<-factor(Base.Acineto.tidy$ST)
Base.Acineto.tidy <- Base.Acineto.tidy %>%
    mutate(
        ORIGEN = case_when(
            ORIGEN == "BOLIVIA" ~ "Bolivia",
            ORIGEN == "BS AS" ~ "Bs.As",
            ORIGEN == "CATAMARCA" ~ "Catamarca.Arg",
            ORIGEN == "CHACO" ~ "Chaco.Arg",
            ORIGEN == "CHILE" ~ "Chile",
            ORIGEN == "ECUADOR" ~ "Ecuador",
            ORIGEN == "HDC" ~ "Htal.Clínicas.Arg",
            ORIGEN == "LAPAMPA" ~ "La Plata.Arg",
            
            ORIGEN == "MENDOZA" ~ "Mendoza.Arg",
            ORIGEN == "SJUAN" ~ "San Juan.Arg",
            ORIGEN == "TUCUMAN" ~ "Tucuman.Arg",
            ORIGEN == "URUGUAY" ~ "Uruguay",
            TRUE ~ as.character(ORIGEN)))
Base.Acineto.tidy$ORIGEN<-factor(Base.Acineto.tidy$ORIGEN)
Base.Acineto.tidy <- Base.Acineto.tidy %>%
    mutate(
        COL = case_when(
            COL == "S" ~ "0.06",
            TRUE ~ as.character(COL)))
Base.Acineto.tidy$COL<-factor(Base.Acineto.tidy$COL)
Base.Acineto.tidy$COL <- factor(Base.Acineto.tidy$COL, 
                                levels = c("0.06", "1", "4","8", "16", "32","64"))
Base.Acineto.tidy$R.colistina <- factor(Base.Acineto.tidy$R.colistina, 
                                        levels = c("S", "R"))

Base.Acineto.tidy <- Base.Acineto.tidy %>%
    mutate(
        Mutacion = case_when(
            Mutacion == "" ~ "ND",
            Mutacion == "sdatos" ~ "ND",
            TRUE ~ as.character(Mutacion)))
Base.Acineto.tidy$Mutacion<-factor(Base.Acineto.tidy$Mutacion,
                                   levels = c("LpxA", "LpxC", "LpxD","LpxA+LpxD",
                                              "LpxD+LpxC", "LpxA+pmrB", "LpxA+LpxC+LpxD+pmrB","pmrB", "ND"))

Base.Acineto.tidy <- Base.Acineto.tidy %>%
    mutate(
        Mutacion.2 = case_when(
            Mutacion == "" ~ "ND",
            Mutacion == "sdatos" ~ "ND",
            Mutacion == "LpxA+pmrB" ~ "Lpx+pmrB",
            Mutacion == "LpxA+LpxC+LpxD+pmrB" ~ "Lpx+pmrB",
            
            TRUE ~ as.character(Mutacion)))
Base.Acineto.tidy$Mutacion.2<-factor(Base.Acineto.tidy$Mutacion.2)
Base.Acineto.tidy$Mutacion.2<-factor(Base.Acineto.tidy$Mutacion.2,
                                     levels = c("LpxA", "LpxC", "LpxD","LpxA+LpxD",
                                                "LpxD+LpxC", "Lpx+pmrB", "pmrB", "ND"))

Base.Acineto.tidy <- Base.Acineto.tidy %>%
    mutate(
        Mutacion.red = case_when(
            Mutacion == "" ~ "ND",
            Mutacion == "sdatos" ~ "ND",
            Mutacion == "LpxA" ~ "Lpx",
            Mutacion == "LpxD" ~ "Lpx",
            Mutacion == "LpxC" ~ "Lpx",
            Mutacion == "LpxD+LpxC" ~ "Lpx",
            Mutacion == "LpxA+LpxD" ~ "Lpx",
            Mutacion == "LpxA+pmrB" ~ "Lpx+pmrB",
            Mutacion == "LpxA+LpxC+LpxD+pmrB" ~ "Lpx+pmrB",
            TRUE ~ as.character(Mutacion)))

Base.Acineto.tidy$Mutacion.red<-factor(Base.Acineto.tidy$Mutacion.red,
                                       levels = c("Lpx", "Lpx+pmrB", "pmrB","ND"))


#### Joining data with spectra data #####

Base.agua.MBT<-AvgSpectra.Acineto.agua.Total.MBT.df%>%
    left_join(Base.Acineto.tidy,by="numero.MT")

sensi.Ave.Factor.Avera<-factor(Base.agua.MBT$R.colistina)
ID.cepa<-Base.agua.MBT$numero.MT

peaks <- AvgSpectra.Acineto.agua.Total.MBT %>%
    detectPeaks( SNR=3, method="MAD", halfWindowSize=30) %>%
    binPeaks(tolerance=0.002) %>%
    filterPeaks(minFrequency=c(0.5, 0.5),
                 labels = sensi.Ave.Factor.Avera,
                 mergeWhitelists=TRUE)

featureMatrix <- intensityMatrix(peaks, AvgSpectra.Acineto.agua.Total.MBT)

mz.full <- as.double(colnames(featureMatrix))
mz <- round( mz.full )
aa<-picos.ducplicados<-any(duplicated(mz))
colnames(featureMatrix)<-mz
table(sensi.Ave.Factor.Avera)

## Clustering visualization

rownames(featureMatrix) <- paste(sensi.Ave.Factor.Avera,ID.cepa,sep=".")
Origen<-Base.agua.MBT$ORIGEN
mutacionn<-Base.agua.MBT$Mutacion.2

Clusters.Plots( feature.matrix= featureMatrix,
                metodologia= "euclidean",
                spot= spot.factor.C,
                sensi=sensi.Ave.Factor.Avera,
                other=paste(Origen,mutacionn,sep=","), 
                n.cluster= 4, 
                title="A. baumannii MBT-MALDI-TOF-MS Intensity data" )

##### Machine learning approachs #####

Sampling.Random(featureMatrix,
                title="Acineto.agua.MBT.Int",
                sensi=sensi.Ave.Factor.Avera,
                seed<-12345,
                train.p<-0.6,
                test.p<-0.4)
load("Train.test.Acineto.agua.MBT.Int.rda")

Test.Intens<-Test
Train.Intens<-Train
Ytest1.Int<-Ytest1
Ytrain11.Int<-Ytrain11

####### Intensity data ##########

### 1. Random forest

Other.ML.Arboles.RandomF(featureMatrix,
                         rowName=sensi.Ave.Factor.Avera,
                         seed=12345,
                         prob.train=0.6,
                         prob.test=0.4,
                         Negative="S",
                         colnum=66,
                         title="Acineto.intensity.MBT")
load("random.Forest.Acineto.intensity.MBT.rda")

Random.forest.peaks.10<-random.Forest[1:10,1]
picos.seleccionados<- as.character(Random.forest.peaks.10$.rownames)
Presentes<-picos.seleccionados %in% colnames(featureMatrix)
Picos<-picos.seleccionados[Presentes]
Seleccion<-featureMatrix[, Picos]

Other.ML.Arboles.RandomF(Seleccion,
                         rowName=sensi.Ave.Factor.Avera,
                         seed=12345,
                         prob.train=0.6,
                         prob.test=0.4,
                         Negative="S",
                         colnum=11,
                         title="Acineto.intensity.MBT")
load("random.Forest.Acineto.intensity.MBT.rda")

### 2. XG.boost

Extreme.g.boost(feature.Matrix=featureMatrix,
                rowName=sensi.Ave.Factor.Avera,
                colnum=66,
                seed=12345,
                prob.train=0.6,
                prob.test=0.4,
                Negative="S",
                nRounds=100,
                colSample_bytree=1,
                min_Child_weight=1,
                Eta=0.4,
                Gamma=0.4,
                subSample=0.5,
                max_Depth=2
)

XG.boost.best.10<-c("3449","2807","6048","8914","7999","8721","3074","2874","2585","8588")
Seleccion.xg<-featureMatrix[, XG.boost.best.10]

Extreme.g.boost(feature.Matrix=Seleccion.xg,
                rowName=sensi.Ave.Factor.Avera,
                colnum=11,
                seed=12345,
                prob.train=0.6,
                prob.test=0.4,
                Negative="S",
                nRounds=100,
                colSample_bytree=1,
                min_Child_weight=1,
                Eta=0.4,
                Gamma=0.4,
                subSample=0.5,
                max_Depth=2
)

### 3. L/D Discriminant analysis

rownames(featureMatrix) <- sensi.Ave.Factor.Avera
Analisis.discriminante(featureMatrix,
                       Ytain=sensi.Ave.Factor.Avera,
                       Diagonal=FALSE,
                       Top=10,
                       title="LDA.Acineto.MBT")
Analisis.discriminante(featureMatrix,
                       Ytain=sensi.Ave.Factor.Avera,
                       Diagonal=TRUE,
                       Top=10,
                       title="DDA.Acineto.LMW")
Predicciones.discriminantes(featureMatrix,
                            Ytain=sensi.Ave.Factor.Avera,
                            TraiN=Train,
                            YTraIn=Ytrain11,
                            TEsT=Test,
                            YTesT=Ytest1,
                            Diagonal=TRUE,
                            TOp=10,
                            seed=2324,
                            title= "Predicciones.DDA.Acineto.MBT",
                            Negative="S")
load(file=paste("Predicciones","Predicciones.DDA.Acineto.MBT","rda",sep="."))

ROC.Curve.3(Ytest=Ytest1,
            yNEW1=ynew.hole.DA,
            d1.1="hole",
            yNEW2=ynew.15.DA,
            d2.1="b.15",
            yNEW3=ynew.5.DA,
            d3.1="b.5",
            title="Intensity DDA analysis",
            Negative="S")

### Binary analysis ### 

Binary.analisis(feature.matrix= featureMatrix,
                Peaks=peaks,
                spot= spot.factor.C,
                sensi=sensi.Ave.Factor.Avera,
                other="",
                title="Acineto.MBT")
load("Matrices.Acineto.MBT.rda")

# Visualization #

Clusters.Plots( feature.matrix= Feature.bina3,
                metodologia= "binary",
                spot= spot.factor.C,
                sensi=sensi.Ave.Factor.Avera,
                other=Origen, 
                n.cluster= 7, 
                title="A. baumannii MBT-MALDI-TOF-MSbinary data" )
Clusters.Plots( feature.matrix= Feature.dicho3,
                metodologia= "binary",
                spot= spot.factor.C,
                sensi=sensi.Ave.Factor.Avera,
                other=paste(Origen,mutacionn,sep="."), 
                n.cluster= 5, 
                title="A. baumannii MBT-MALDI-TOF-MS dichotomized data" )

## Sampling random dichotomized ##

Sampling.Random(Feature.dicho3,
                title="Acineto.MBT.dicho",
                sensi=sensi.Ave.Factor.Avera,
                seed<-12345,
                train.p<-0.6,
                test.p<-0.4)
load("Train.test.Acineto.MBT.dicho.rda")

## BDA Predictions dichotomized ##

Predicciones.binda(Feature.dicho3,
                   sensi.Ave.Factor.Avera,
                   Train,
                   Ytrain11,
                   Test,
                   Ytest1,
                   Negative="S",
                   title= "Dicho.data")
load("Binda.Dicho.data.rda")

ROC.Curve.3(Ytest=Ytest1,
            yNEW1=ynew.hole,
            d1.1="hole",
            yNEW2=ynew15,
            d2.1="b.15",
            yNEW3=ynew5,
            d3.1="b.5",
            title="BDA dichotomized data",
            Negative="S")

## BDA Cross-validation dichotomized ##

Crossvalidacion(Feature.dicho3,
                sensi.Ave.Factor.Avera,
                Train,
                Ytrain11,
                Test,
                Ytest1,
                Negative="S",
                seed=12345,
                k=5,
                b=20,
                title="Acineto.Dicho")
load("Cv.Acineto.Dicho.rda")


## Random forest dichotomized ##

Other.ML.Arboles.RandomF(Feature.dicho3,
                         rowName=sensi.Ave.Factor.Avera,
                         seed=12345,
                         prob.train=0.6,
                         prob.test=0.4,
                         Negative="S",
                         colnum=66,
                         title="Acineto.dicho.MBT")

load("random.Forest.Acineto.dicho.MBT.rda")

Random.forest.peaks.10<-random.Forest[1:10,1]
picos.seleccionados<- as.character(Random.forest.peaks.10$.rownames)
Presentes<-picos.seleccionados %in% colnames(featureMatrix)
Picos<-picos.seleccionados[Presentes]
Seleccion<-Feature.dicho3[, Picos]

Other.ML.Arboles.RandomF(Seleccion,
                         rowName=sensi.Ave.Factor.Avera,
                         seed=12345,
                         prob.train=0.6,
                         prob.test=0.4,
                         Negative="S",
                         colnum=11,
                         title="Acineto.dicho.MBT")

## Xg.boost dichotomized ##

Extreme.g.boost(feature.Matrix=Feature.dicho3,
                rowName=sensi.Ave.Factor.Avera,
                colnum=66,
                seed=12345,
                prob.train=0.6,
                prob.test=0.4,
                Negative="S",
                nRounds=110,
                colSample_bytree=1,
                min_Child_weight=1,
                Eta=0.4,
                Gamma=0.4,
                subSample=0.5,
                max_Depth=2)

load("Xg.boost.rda")
Xg.boost.dich<-impMatrix[,1]
Xg.boost.dich<-separate(Xg.boost.dich, Feature, c("X","Picos"), sep="X") 

Xg.boost.dich<-Xg.boost.dich$Picos

Xg.boost.dich<-c("3449","3968","6048","11042","3189","6248","6642","8323","8489","5177")
Seleccion.xg<-Feature.dicho3[, Xg.boost.dich]

Extreme.g.boost(feature.Matrix=Seleccion.xg,
                rowName=sensi.Ave.Factor.Avera,
                colnum=11,
                seed=12345,
                prob.train=0.6,
                prob.test=0.4,
                Negative="S",
                nRounds=100,
                colSample_bytree=1,
                min_Child_weight=1,
                Eta=0.4,
                Gamma=0.4,
                subSample=0.5,
                max_Depth=2)

### Machine learning prospective validation ###

so4<-importBrukerFlex("~/ACINETO.1/cepas23_4_19", verbose=FALSE)

Names<-c(rep("162",9),rep("166",9),rep("167",6),rep("385",9),rep("h22",9))
Names.arr<-array(Names)
for (h in 1:length(so4)) {
    metaData(so4[[h]])$spot<-Names.arr[[h]]
}
table(lengths(so4))
any(sapply(so4, MALDIquant::isEmpty))
all(sapply(so4, isRegular))
so5<-so4
spot.factor <- factor(sapply(so4,function(x)metaData(x)$spot))
Avg.generation(SO=so4,
               transform.method="sqrt",
               labells=spot.factor,
               smoothmethod="SavitzkyGolay",
               halfwind=30,
               removemethod="SNIP",
               iteration=75,
               calibrateMEto="TIC",
               Snr=3,
               Tolerance=0.002,
               warpingmeto="quadratic",
               avemeto="sum",
               title="Nuevascepas")

load("avgSpectra.Nuevascepas.rda")

spot.factor.23_4_19 <- factor(sapply(avgSpectra,function(x)metaData(x)$spot))
peaks.23_4_19 <- avgSpectra %>%
    detectPeaks( SNR=3, method="MAD", halfWindowSize=30) %>%
    binPeaks(tolerance=0.002) %>%
    filterPeaks( minFrequency=c(0.5, 0.5),
                 mergeWhitelists=TRUE)
featureMatrix.23_4_19 <- intensityMatrix(peaks.23_4_19, avgSpectra)
mz.full <- as.double(colnames(featureMatrix.23_4_19))
mz <- round( mz.full )
picos.ducplicados<-any(duplicated(mz))
colnames(featureMatrix.23_4_19)<-mz
rownames(featureMatrix.23_4_19)<-spot.factor.23_4_19
Matrix.23_4_19<-featureMatrix.23_4_19

### BDA with best 10 peaks ###

Sampling.Random(Feature.dicho3,
                title="Acineto.MBT.dicho",
                sensi=sensi.Ave.Factor.Avera,
                seed<-12345,
                train.p<-0.6,
                test.p<-0.4)
load("Train.test.Acineto.MBT.dicho.rda")

br.d = binda.ranking(Feature.dicho3, sensi.Ave.Factor.Avera)
selPeaks10.d<-br.d[,"idx"][1:10]
picos<-data.frame(selPeaks20.d)
Picos.10.BDA<-rownames(picos)

# 1. Selection of model peaks that were present in newdata

picos.seleccionados.bd<-as.numeric(Picos.20.BDA)/10
picos.seleccionados.bd <- round(picos.seleccionados.bd)

Matrix.23_4_19.20.ceiling<-Matrix.23_4_19
mz.cep <- as.double(colnames(Matrix.23_4_19.20.ceiling))
mz.cep<-mz.cep/10
mz.cep <- round( mz.cep )
colnames(Matrix.23_4_19.20.ceiling)<-mz.cep

Matrix.23_4_19.20.ceiling.df<-as.data.frame(Matrix.23_4_19.20.ceiling)
Presentes.bd<-picos.seleccionados.bd %in% colnames(Matrix.23_4_19.20.ceiling.df)
Picos.bd.new<-as.character(picos.seleccionados.bd[Presentes.bd])
Seleccion.b.new<-Matrix.23_4_19.20.ceiling[, Picos.bd.new]
Seleccion.b.new.df<-as.data.frame(Seleccion.b.new)

means<-map(Seleccion.b.new.df,mean)
means.df<-as.data.frame(means)
MEAN<-as.numeric(means.df)
Seleccion.b.df.b <- dichotomize(Seleccion.b.new.df, MEAN)

# 2. Subsetting model peaks with present model peaks on newdata

mz.t<-colnames(Train)
mz.v<-mz.t
mz.v<-as.numeric(mz.v)/10
mz.v <- round(mz.v)
colnames(Train)<-mz.v

Train.df<-as.data.frame(Train)
Presentes.bd.train<-Picos.bd.new %in% colnames(Train.df)
Picos.bd.train<-as.character(Picos.bd.new[Presentes.bd.train])
Seleccion.b.train<-Train[, Picos.bd.train]

#. 3 Prediction
binda.out10.d = binda(Seleccion.b.train[, Picos.bd.new, drop=FALSE], Ytrain11, verbose=FALSE)
ynew10 = predict.binda(binda.out10.d, Seleccion.b.df.b, verbose=FALSE)$posterior
ynew10


## Random forest with 10 best peaks

Other.ML.Arboles.RandomF(Feature.dicho3,
                         rowName=sensi.Ave.Factor.Avera,
                         seed=12345,
                         prob.train=0.6,
                         prob.test=0.4,
                         Negative="S",
                         colnum=66,
                         title="Acineto.dicho.MBT")

load("random.Forest.Acineto.dicho.MBT.rda")

Random.forest.peaks.10<-random.Forest[1:25,1]
picos.seleccionados<- as.character(Random.forest.peaks.10$.rownames)

# 1. Selection of model peaks that were present in newdata

picos.seleccionados.rf<-as.numeric(picos.seleccionados)/10
picos.seleccionados.rf <- ceiling(picos.seleccionados.rf)

Matrix.23_4_19.20.ceiling<-Matrix.23_4_19
mz.cep <- as.double(colnames(Matrix.23_4_19.20.ceiling))
mz.cep<-mz.cep/10
mz.cep <- round( mz.cep )
colnames(Matrix.23_4_19.20.ceiling)<-mz.cep

Matrix.23_4_19.20.ceiling.df<-as.data.frame(Matrix.23_4_19.20.ceiling)
Presentes.rf<-picos.seleccionados.rf %in% colnames(Matrix.23_4_19.20.ceiling.df)
Picos.rf.new<-as.character(picos.seleccionados.rf[Presentes.rf])
Seleccion.r.new<-Matrix.23_4_19.20.ceiling[, Picos.rf.new]
Seleccion.r.new.df<-as.data.frame(Seleccion.r.new)

means<-map(Seleccion.r.new.df,mean)
means.df<-as.data.frame(means)
MEAN<-as.numeric(means.df)
Seleccion.rf.df.d <- dichotomize(Seleccion.r.new.df, MEAN)

# 2. Subsetting model peaks with present model peaks on newdata

Sampling.Random(Feature.dicho3,
                title="Acineto.MBT.dicho",
                sensi=sensi.Ave.Factor.Avera,
                seed<-12345,
                train.p<-0.6,
                test.p<-0.4)
load("Train.test.Acineto.MBT.dicho.rda")

mz.t<-colnames(Train)
mz.v<-mz.t
mz.v<-as.numeric(mz.v)/10
mz.v <- ceiling(mz.v)
colnames(Train)<-mz.v

Train.df<-as.data.frame(Train)
Presentes.rf.train<-Picos.rf.new %in% colnames(Train.df)
Picos.rf.train<-as.character(Picos.rf.new[Presentes.rf.train])
Seleccion.rf.train<-Train[, Picos.rf.train]

# 3. Prediction

set.seed(12345)
rf.pros <- randomForest(Seleccion.rf.train,Ytrain11, importance=TRUE,
                        proximity=TRUE)
rf.biop.test<- predict(rf.pros, newdata=Seleccion.rf.df.d, type="prob")
Random.f<-data.frame(rf.biop.test)
colnames(Random.f)<-c("Prob.R.rf","Prob.S.rf")

rf.biop.test.resp<- predict(rf.pros, newdata=Seleccion.rf.df.d, type="response")
Random.f.b<-data.frame(rf.biop.test.resp)
colnames(Random.f.b)<-c("Rf")

### Xg.boost dichotomized with 10 best peaks

Extreme.g.boost(feature.Matrix=Feature.dicho3,
                rowName=sensi.Ave.Factor.Avera,
                colnum=66,
                seed=12345,
                prob.train=0.6,
                prob.test=0.4,
                Negative="S",
                nRounds=110,
                colSample_bytree=1,
                min_Child_weight=1,
                Eta=0.4,
                Gamma=0.4,
                subSample=0.5,
                max_Depth=2
                
)

load("Xg.boost.rda")
Xg.boost.dich<-impMatrix[,1]
Xg.boost.dich<-separate(Xg.boost.dich, Feature, c("X","Picos"), sep="X") 
Xg.boost.dich<-Xg.boost.dich$Picos
Xg.boost.dich<-Xg.boost.dich[1:15]

# 1. Selection of model peaks that were present in newdata

picos.seleccionados.xg<-as.numeric(Xg.boost.dich)/10
picos.seleccionados.xg <- ceiling(picos.seleccionados.xg)

Matrix.23_4_19.20.ceiling<-Matrix.23_4_19
mz.cep <- as.double(colnames(Matrix.23_4_19.20.ceiling))
mz.cep<-mz.cep/10
mz.cep <- round( mz.cep )
colnames(Matrix.23_4_19.20.ceiling)<-mz.cep

Matrix.23_4_19.20.ceiling.df<-as.data.frame(Matrix.23_4_19.20.ceiling)
Presentes.xg<-picos.seleccionados.xg %in% colnames(Matrix.23_4_19.20.ceiling.df)
Picos.xg.new<-as.character(picos.seleccionados.xg[Presentes.xg])
Seleccion.xg.new<-Matrix.23_4_19.20.ceiling[, Picos.xg.new]
Seleccion.xg.new.df<-as.data.frame(Seleccion.xg.new)

means<-map(Seleccion.xg.new.df,mean)
means.df<-as.data.frame(means)
MEAN<-as.numeric(means.df)
Seleccion.xg.df.d <- dichotomize(Seleccion.xg.new.df, MEAN)

# 2. Subsetting model peaks with present model peaks on newdata

Sampling.Random(Feature.dicho3,
                title="Acineto.MBT.dicho",
                sensi=sensi.Ave.Factor.Avera,
                seed<-12345,
                train.p<-0.6,
                test.p<-0.4)
load("Train.test.Acineto.MBT.dicho.rda")

mz.t<-colnames(Train)
mz.v<-mz.t
mz.v<-as.numeric(mz.v)/10
mz.v <- ceiling(mz.v)
colnames(Train)<-mz.v

Train.df<-as.data.frame(Train)
Presentes.xg.train<-Picos.xg.new %in% colnames(Train.df)
Picos.xg.train<-as.character(Picos.xg.new[Presentes.xg.train])
Seleccion.xg.train<-Train[, Picos.xg.train]

# 3. Prediction

grid <- expand.grid(nrounds = 100, colsample_bytree = 0.8, min_child_weight = 1,eta = 0.3,
                    gamma = 0.4, subsample = 0.5, max_depth = 3) 
cntrl <- trainControl( method = "cv",number = 5,verboseIter = FALSE,returnData = FALSE,
                       returnResamp = "final")

set.seed(12345)
train.xgb = train(x = Seleccion.xg.train, y = Ytrain11, trControl = cntrl, tuneGrid = grid,
                  method = "xgbTree")
param <- list( objective = "binary:logistic", 
               booster = "gbtree",
               eval_metric = "error",
               eta = 0.3, 
               max_depth = 3, 
               subsample = 0.5,
               colsample_bytree = 0.8,
               gamma = 0.4)
x <- as.matrix(Seleccion.xg.train)
y <- ifelse(Ytrain11 == "S", 0, 1)
train.mat <- xgb.DMatrix(data = x, label = y)
set.seed(12345)
xgb.fit <- xgb.train(params = param, data = train.mat, nrounds =100)
pred <- predict(xgb.fit, x)
threshold<-InformationValue::optimalCutoff(y, pred)
threshold

M_1.test_X <- as.matrix(Seleccion.xg.df.d)
xgb.prima.test <- predict(xgb.fit, M_1.test_X, type="prob")
xgb.prima.test.b <- ifelse(xgb.prima.test > threshold, "R", "S")    

final.table<-data.frame(rownames(M_1.test_X),xgb.prima.test)
colnames(final.table)<-c("muestra","Prob.Xg")

final.table.b<-data.frame(rownames(M_1.test_X),xgb.prima.test.b)
colnames(final.table.b)<-c("muestra","Xg")



##### Chi-squared test #######

# 1. ST vs R.col.categor
xsq_Endo<-chisq.test(Base.Acineto.tidy$ST, 
                     Base.Acineto.tidy$R.colistina, 
                     simulate.p.value=TRUE)
desv<-xsq_Endo$stdres
desv<-as.data.frame(desv)
desv<-desv %>%
    mutate(ST= as.factor(desv$Base.Acineto.tidy.ST),
           Colistin= factor(desv$Base.Acineto.tidy.R.colistina,levels = c("R", "S")))
ggplot(data= desv, aes(x=ST, y= Freq,  fill=Colistin, color= Colistin)) +
    geom_point(size=4) +
    scale_fill_brewer() +
    theme(axis.text.x = element_text(family="Times",
                                     angle = 0, 
                                     hjust = 1,
                                     size=rel(1.1)),
          legend.position = "bottom",
          legend.title = element_text(colour="white"))+
    geom_hline(yintercept = 2, col="blue")+
    geom_hline(yintercept = -2, col="blue") +
    labs(x="ST", y="Pearson residuals-Chi squared test")  

# 2. ST vs Origen

xsq_Endo.1<-chisq.test(Base.Acineto.tidy$ST, 
                       Base.Acineto.tidy$ORIGEN, 
                       simulate.p.value=TRUE)
desv1<-xsq_Endo.1$stdres
desv1<-as.data.frame(desv1)
desv1<-desv1 %>%
    mutate(ST= as.factor(desv1$Base.Acineto.tidy.ST),
           Origen= factor(desv1$Base.Acineto.tidy.ORIGEN))
ggplot(data= desv1, aes(x=ST, y= Freq,  fill=Origen, color= Origen)) +
    geom_point(size=4) +
    scale_fill_brewer() +
    theme(axis.text.x = element_text(family="Times",
                                     angle = 0, 
                                     hjust = 1,
                                     size=rel(1.1)),
          legend.position = "bottom",
          legend.title = element_text(colour="white"))+
    geom_hline(yintercept = 2, col="blue")+
    geom_hline(yintercept = -2, col="blue") +
    labs(x="ST", y="Pearson residuals-Chi squared test")  

# 3. ST vs Mutación.red

xsq_Endo.1<-chisq.test(Base.Acineto.tidy$ST, 
                       Base.Acineto.tidy$Mutacion.red, 
                       simulate.p.value=TRUE)
desv1<-xsq_Endo.1$stdres
desv1<-as.data.frame(desv1)
desv1<-desv1 %>%
    mutate(ST= as.factor(desv1$Base.Acineto.tidy.ST),
           Mutation= factor(desv1$Base.Acineto.tidy.Mutacion.red))
ggplot(data= desv1, aes(x=ST, y= Freq,  fill=Mutation, color= Mutation)) +
    geom_point(size=4) +
    scale_fill_brewer() +
    theme(axis.text.x = element_text(family="Times",
                                     angle = 0, 
                                     hjust = 1,
                                     size=rel(1.1)),
          legend.position = "bottom",
          legend.title = element_text(colour="white"))+
    geom_hline(yintercept = 2, col="blue")+
    geom_hline(yintercept = -2, col="blue") +
    labs(x="ST", y="Pearson residuals-Chi squared test")  

# 4.1 R.col.CIM vs mutacion

Base.Acineto.tidy.R<-Base.Acineto.tidy %>%
    filter(R.colistina=="R")

xsq_Endo.1<-chisq.test(Base.Acineto.tidy.R$COL, 
                       Base.Acineto.tidy.R$Mutacion.2, 
                       simulate.p.value=TRUE)
desv1<-xsq_Endo.1$stdres
desv1<-as.data.frame(desv1)
desv1<-desv1 %>%
    mutate(Colistin.CIM= as.factor(desv1$Base.Acineto.tidy.R.COL),
           Mutation= factor(desv1$Base.Acineto.tidy.R.Mutacion.2))
ggplot(data= desv1, aes(x=Colistin.CIM, y= Freq,  fill=Mutation, color= Mutation)) +
    geom_point(size=4) +
    scale_fill_brewer() +
    theme(axis.text.x = element_text(family="Times",
                                     angle = 0, 
                                     hjust = 1,
                                     size=rel(1.1)),
          legend.position = "bottom",
          legend.title = element_text(colour="white"))+
    geom_hline(yintercept = 2, col="blue")+
    geom_hline(yintercept = -2, col="blue") +
    labs(x="CIM (µg/ml)", y="Pearson residuals-Chi squared test")  

# 4.2 R.col.CIM vs mutacion

xsq_Endo.1<-chisq.test(Base.Acineto.tidy.R$COL, 
                       Base.Acineto.tidy.R$Mutacion.red, 
                       simulate.p.value=TRUE)

desv1<-xsq_Endo.1$stdres
desv1<-as.data.frame(desv1)
desv1<-desv1 %>%
    mutate(Colistin.CIM= as.factor(desv1$Base.Acineto.tidy.R.COL),
           Mutation= factor(desv1$Var2))
ggplot(data= desv1, aes(x=Colistin.CIM, y= Freq,  fill=Mutation, color= Mutation)) +
    geom_point(size=4) +
    scale_fill_brewer() +
    theme(axis.text.x = element_text(family="Times",
                                     angle = 0, 
                                     hjust = 1,
                                     size=rel(1.1)),
          legend.position = "bottom",
          legend.title = element_text(colour="white"))+
    geom_hline(yintercept = 2, col="blue")+
    geom_hline(yintercept = -2, col="blue") +
    labs(x="CIM (µg/ml)", y="Pearson residuals-Chi squared test")  

###### Geolocalization ########

library(sp)
library(rgdal)
library(leaflet)
library(htmltools)

Geo.Acineto<-read.csv2("/Volumes/TOSHIBA EXT/Mac/trabajos MALDI/Resistencia.ATB/Acinetobacter/Análisis hernan/Analisis Final/Geo.Acineto.csv", stringsAsFactors = FALSE)
Geo.Acineto$lat<-as.numeric(Geo.Acineto$lat)
Geo.Acineto$lon<-as.numeric(Geo.Acineto$lon)

Geo.Acineto$REGION<-c("Bolivia","Argentina","Argentina","Argentina","Chile","Ecuador","Argentina",
                      "Argentina","Argentina","Argentina","Argentina","Uruguay")

Geo.Acineto$Provincia<-c("Bolivia","Buenos Aires","Catamarca","Chaco","Chile","Ecuador","Buenos Aires","Buenos Aires","Mendoza","San Juan","Tucumán","Uruguay")

Geo.Acineto.ST<-Geo.Acineto %>%
    group_by(Provincia,REGION) %>%
    summarise(Porcentaje=sum(Porcentaje),
              N.ST1=sum(ST.1),
              N.ST5=sum(ST.5),
              N.ST15=sum(ST.15),
              N.ST25=sum(ST.25),
              N.ST79=sum(ST.79),
              N.R=sum(X..R))%>%
    as.data.frame()

# Argentina
tmp <- tempdir()
url <- "http://biogeo.ucdavis.edu/data/diva/adm/ARG_adm.zip"
file <- basename(url)
download.file(url, file)
unzip(file, exdir = tmp)
argentina <- readOGR(dsn = tmp, layer = "ARG_adm1", use_iconv=TRUE, encoding='UTF-8')
Geo.Acineto.ST.Arg<-Geo.Acineto.ST %>%
    filter(REGION == "Argentina") 
argentina.2<-merge(argentina,Geo.Acineto.ST.Arg,by.x= "NAME_1",by.y= "Provincia")

# Chile
tmp <- tempdir()
url <- "http://biogeo.ucdavis.edu/data/diva/adm/CHL_adm.zip"
file <- basename(url)
download.file(url, file)
unzip(file, exdir = tmp)
Chile <- readOGR(dsn = tmp, layer = "CHL_adm1", use_iconv=TRUE, encoding='UTF-8')

Geo.Acineto.Chile<-Geo.Acineto %>%
    filter(REGION == "Chile") %>%
    group_by(Provincia) %>%
    summarise(cantidad=sum(Porcentaje))%>%
    as.data.frame()

Geo.Acineto.ST.Chile<-Geo.Acineto.ST %>%
    filter(REGION == "Chile") 

colnames(Geo.Acineto.Chile) <- c("NAME_1", "Cantidad")
Chile.2<-merge(Chile,Geo.Acineto.ST.Chile,by.x= "NAME_0",by.y= "REGION")

# Bolivia
tmp <- tempdir()
url <- "http://biogeo.ucdavis.edu/data/diva/adm/BOL_adm.zip"
file <- basename(url)
download.file(url, file)
unzip(file, exdir = tmp)
Bolivia <- readOGR(dsn = tmp, layer = "BOL_adm1", use_iconv=TRUE, encoding='UTF-8')

Geo.Acineto.Bol<-Geo.Acineto %>%
    filter(REGION == "Bolivia") %>%
    group_by(Provincia) %>%
    summarise(cantidad=sum(Porcentaje))%>%
    as.data.frame()

Geo.Acineto.ST.Bol<-Geo.Acineto.ST %>%
    filter(REGION == "Bolivia") 

colnames(Geo.Acineto.Bol) <- c("NAME_1", "Cantidad")
Bolivia.2<-merge(Bolivia,Geo.Acineto.ST.Bol,by.x= "NAME_0",by.y= "REGION")

# Uruguay

tmp <- tempdir()
url <- "http://biogeo.ucdavis.edu/data/diva/adm/URY_adm.zip"
file <- basename(url)
download.file(url, file)
unzip(file, exdir = tmp)
Uruguay <- readOGR(dsn = tmp, layer = "URY_adm1", use_iconv=TRUE, encoding='UTF-8')

Geo.Acineto.Uru<-Geo.Acineto %>%
    filter(REGION == "Uruguay") %>%
    group_by(Provincia) %>%
    summarise(cantidad=sum(Porcentaje))%>%
    as.data.frame()

Geo.Acineto.ST.Uru<-Geo.Acineto.ST %>%
    filter(REGION == "Uruguay") 

colnames(Geo.Acineto.Uru) <- c("NAME_1", "Cantidad")
Uruguay.2<-merge(Uruguay,Geo.Acineto.ST.Uru,by.x= "NAME_0",by.y= "REGION")

# Ecuador

tmp <- tempdir()
url <- "http://biogeo.ucdavis.edu/data/diva/adm/ECU_adm.zip"
file <- basename(url)
download.file(url, file)
unzip(file, exdir = tmp)
Ecuador <- readOGR(dsn = tmp, layer = "ECU_adm1", use_iconv=TRUE, encoding='UTF-8')

Geo.Acineto.Ecu<-Geo.Acineto %>%
    filter(REGION == "Ecuador") %>%
    group_by(Provincia) %>%
    summarise(cantidad=sum(Porcentaje))%>%
    as.data.frame()

Geo.Acineto.ST.Ecu<-Geo.Acineto.ST %>%
    filter(REGION == "Ecuador") 

colnames(Geo.Acineto.Ecu) <- c("NAME_1", "Cantidad")
Ecuador.2<-merge(Ecuador,Geo.Acineto.ST.Ecu,by.x= "NAME_0",by.y= "REGION")

## Map ## 

leaflet() %>% 
    addTiles() %>%
    addProviderTiles("CartoDB.Positron") %>%   
    addPolygons(data = argentina.2, 
                fillColor = ~pal(Porcentaje), 
                fillOpacity = 0.8, 
                color = "#BDBDC3", 
                weight = 1, 
                popup = paste0("<strong>Provincia: </strong>", 
                               argentina.2$NAME_1, 
                               "<br><strong>Cantidad de casos: </strong>", 
                               argentina.2$Porcentaje),
                group = "Argentina") %>% 
    addPolygons(data = Chile.2, fillColor = ~pal(Porcentaje), 
                fillOpacity = 0.8, 
                color = "#BDBDC3", 
                weight = 1, 
                popup = paste0("<strong>Estado: </strong>", 
                               Chile.2$NAME_0, 
                               "<br><strong>Cantidad de casos: </strong>", 
                               Chile.2$Porcentaje),
                group = "Chile") %>% 
    addPolygons(data = Bolivia.2, fillColor = ~pal(Porcentaje), 
                fillOpacity = 0.8, 
                color = "#BDBDC3", 
                weight = 1, 
                popup = paste0("<strong>Estado: </strong>", 
                               Bolivia.2$NAME_0, 
                               "<br><strong>Cantidad de casos: </strong>", 
                               Bolivia.2$Porcentaje),
                group = "Bolivia") %>% 
    addPolygons(data = Uruguay.2, fillColor = ~pal(Porcentaje), 
                fillOpacity = 0.8, 
                color = "#BDBDC3", 
                weight = 1, 
                popup = paste0("<strong>Estado: </strong>", 
                               Uruguay.2$NAME_0, 
                               "<br><strong>Cantidad de casos: </strong>", 
                               Uruguay.2$Porcentaje),
                group = "Uruguay") %>% 
    addPolygons(data = Ecuador.2, fillColor = ~pal(Porcentaje), 
                fillOpacity = 0.8, 
                color = "#BDBDC3", 
                weight = 1, 
                popup = paste0("<strong>Estado: </strong>", 
                               Ecuador.2$NAME_0, 
                               "<br><strong>Cantidad de casos: </strong>", 
                               Ecuador.2$Porcentaje),
                group = "Ecuador") 








