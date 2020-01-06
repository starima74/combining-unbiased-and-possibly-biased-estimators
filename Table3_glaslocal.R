#####################################################################################################
# This function reproduces results for global methods run locally (15%) presented in table 3,middle
#####################################################################################################


Table3_glaslocal <- function(){
  
  # distance  
  
  dist <- function(all, one, w=NA) {
    all <- all[,names(one)]
    ret = as.vector(rep(0,nrow(all)))
    Vs <- diag(var(all))  
    if(is.na(sum(w))) w = rep(1,ncol(all)) 
    for(d in 1:ncol(one))
      ret <- ret + as.numeric(w[d]) * ( as.vector(all[,d])  - as.numeric(one[d]))^2 / as.numeric(Vs[d]) 
    ret
  }
  
  f_full <- formula(paste("D6 ~ A1C0 + age  ", sep=""))
  
  
  taus <- c(0.5)
  
  a.train<- read.csv("a.train.csv")
  a.test <- read.csv("a.test.csv")
  a.train$newcurr_met <- as.numeric(a.train$newcurr_met)
  a.train$neworder_insulin <- as.numeric(a.train$neworder_insulin)
  
  a.test$newcurr_met <- as.numeric(a.test$newcurr_met)
  a.test$neworder_insulin <- as.numeric(a.test$neworder_insulin)
  
  
  BART <- NULL
  RR <- NULL
  NET <- NULL
  GBM <- NULL
  RF <- NULL
  GAM <- NULL 
  QR_nonlin <- NULL
  LM <- NULL  
  
  NT <- nrow(a.test)
  #To predict using a small subset of the test data set, uncomment 
  #NT <- 500
  
  for(id in 1:nrow(a.test[1:NT,])) { 
    t <- 1
    pred_vars <- c("A1C0","age")
    print(id)
    newperson <- a.test[id,pred_vars ]
    one <- a.test[id,c("age", "A1C0")]
    dists <- dist(a.train,one, w = c(1, 3))
    a.train$dists <- dists
    cutoff <- quantile(dists, probs = 0:20/20)[4] # 15% cutoff
    ind <- a.train$dists<cutoff
    
    #BART   
    
    dat <- a.train[ind,c("D6",pred_vars)]  
    x.train <- a.train[ind,pred_vars]
    x.test <- newperson
    y.train <-  a.train[ind,"D6"] 
    tau <- taus[t]
    
    a.trainnew <- a.train[ind,]
    
    fitBART <-  wbart( x.train , y.train,x.test)   
    BART[id] <- fitBART$yhat.test.mean
    
    a.trainnew$predBART <-  fitBART$yhat.train.mean  
    newperson$predBART <- fitBART$yhat.test.mean
    
    #Elastic net  
    fitNet <- train(D6~.,data=dat,method = "glmnet") 
    NET[id] <- predict(fitNet,x.test,submodels=NULL)
    a.trainnew$predNET <- predict(fitNet,x.train,submodels=NULL)
    newperson$predNET <- predict(fitNet,x.test,submodels=NULL)
    
    #Ridge Regression  
    fitRR <- train(D6~.,data=dat,method = "ridge") 
    RR[id] <- predict(fitRR, a.test[c(id,id),pred_vars],submodels=NULL)[1]
    a.trainnew$predRR <- predict(fitRR,a.train[ind,pred_vars],submodels=NULL)
    newperson$predRR <- predict(fitRR,a.test[c(id,id),pred_vars],submodels=NULL)[1]
    
    # GBM
    gbmGrid <- expand.grid(interaction.depth=c(1,5,9),n.trees=(1:50)*50,
                           shrinkage=0.1, n.minobsinnode=10) 
    fitControl <- trainControl(method="cv",number=10)  
    fitGBM <- train(D6~.,data=dat,method = "gbm",verbose=F,
                    trControl=fitControl,tuneGrid=gbmGrid) 
    GBM[id] <-  predict(fitGBM,x.test,submodels=NULL)
    
    a.trainnew$predGBM <- predict(fitGBM,x.train,submodels=NULL)
    newperson$predGBM <- predict(fitGBM,x.test,submodels=NULL)
    
    #RANDOM FOREST
    fitRF <- train(x.train,y.train,method = "rf",importance=T)
    RF[id] <- predict(fitRF,x.test,submodels=NULL)
    a.trainnew$predRF <- predict(fitRF,x.train,submodels=NULL)
    newperson$predRF <- predict(fitRF,x.test,submodels=NULL)
    
    #Generalized additive model
    fitGAM <- gam(f_full, data = dat)
    GAM[id] <-   predict(fitGAM,newdata=x.test) 
    a.trainnew$predGAM <- predict(fitGAM,x.train)
    newperson$predGAM <- predict(fitGAM,x.test)
    
    
    #QR_nonparametric predictors 
    fitQR_nonlin <- rq(D6~bs(A1C0,5)+bs(age,5) , tau=taus, data=dat)
    QR_nonlin[id] <- predict(fitQR_nonlin,newdata=x.test)
    a.trainnew$predQR_nonlin <- predict(fitQR_nonlin,x.train)
    newperson$predQR_nonlin <- predict(fitQR_nonlin,x.test)
    
    #LM with linear predictors
    y0hat <- predict(lm(f_full,data=dat), newdata = x.test, se=T) 
    LM[id] <- y0hat$fit + qnorm(taus) * y0hat$se.fit
    newperson$predLM <- y0hat$fit + qnorm(taus) * y0hat$se.fit
    
    y0hat_p <- predict(lm(f_full,data=dat), newdata = x.train, se=T)
    a.trainnew$predLM <- y0hat_p$fit + qnorm(taus) * y0hat_p$se.fit
    
  }
  
  ######################################################################
  # Estimate bias, L1 and L2 measures for the above estimators
  #########################################################################
  
  #BIAS  
  
  bias_BART <- mean( I(a.test$D6[1:NT]<BART)* (a.test$D6[1:NT]-BART) +I(a.test$D6[1:NT] >= BART)* (a.test$D6[1:NT]-BART))
  bias_NET<- mean( I(a.test$D6[1:NT]<NET)* (a.test$D6[1:NT]-NET) +I(a.test$D6[1:NT] >= NET)* (a.test$D6[1:NT]-NET))
  bias_RR <- mean( I(a.test$D6[1:NT]<RR)* (a.test$D6[1:NT]-RR) +I(a.test$D6[1:NT] >= RR)* (a.test$D6[1:NT]-RR))
  bias_GBM<- mean( I(a.test$D6[1:NT]<GBM)* (a.test$D6[1:NT]-GBM) +I(a.test$D6[1:NT] >= GBM)* (a.test$D6[1:NT]-GBM))
  bias_GAM <- mean( I(a.test$D6[1:NT]<GAM)* (a.test$D6[1:NT]-GAM) +I(a.test$D6[1:NT] >= GAM)* (a.test$D6[1:NT]-GAM))
  bias_QR_nonlin <- mean( I(a.test$D6[1:NT]<QR_nonlin)* (a.test$D6[1:NT]-QR_nonlin) +I(a.test$D6[1:NT] >= QR_nonlin)* (a.test$D6[1:NT]-QR_nonlin))
  bias_LM <- mean( I(a.test$D6[1:NT]<LM)* (a.test$D6[1:NT]-LM) +I(a.test$D6[1:NT] >= LM)* (a.test$D6[1:NT]-LM))
  bias_RF <- mean( I(a.test$D6[1:NT]<RF)* (a.test$D6[1:NT]-RF) +I(a.test$D6[1:NT] >= RF)* (a.test$D6[1:NT]-RF))
  
  #L1
  
  L1_BART <- mean( I(a.test$D6[1:NT]<BART)*abs(a.test$D6[1:NT]-BART) +I(a.test$D6[1:NT] >= BART)*abs(a.test$D6[1:NT]-BART))
  L1_NET<- mean( I(a.test$D6[1:NT]<NET)*abs(a.test$D6[1:NT]-NET) +I(a.test$D6[1:NT] >= NET)*abs(a.test$D6[1:NT]-NET))
  L1_RR <- mean( I(a.test$D6[1:NT]<RR)*abs(a.test$D6[1:NT]-RR) +I(a.test$D6[1:NT] >= RR)*abs(a.test$D6[1:NT]-RR))
  L1_GBM<- mean( I(a.test$D6[1:NT]<GBM)*abs(a.test$D6[1:NT]-GBM) +I(a.test$D6[1:NT] >= GBM)*abs(a.test$D6[1:NT]-GBM))
  L1_GAM <- mean( I(a.test$D6[1:NT]<GAM)*abs(a.test$D6[1:NT]-GAM) +I(a.test$D6[1:NT] >= GAM)*abs(a.test$D6[1:NT]-GAM))
  L1_QR_nonlin <- mean( I(a.test$D6[1:NT]<QR_nonlin)*abs(a.test$D6[1:NT]-QR_nonlin) +I(a.test$D6[1:NT] >= QR_nonlin)*abs(a.test$D6[1:NT]-QR_nonlin))
  L1_LM <- mean( I(a.test$D6[1:NT]<LM)*abs(a.test$D6[1:NT]-LM) +I(a.test$D6[1:NT] >= LM)*abs(a.test$D6[1:NT]-LM))
  L1_RF <- mean( I(a.test$D6[1:NT]<RF)*abs(a.test$D6[1:NT]-RF) +I(a.test$D6[1:NT] >= RF)*abs(a.test$D6[1:NT]-RF))
  
  
  #MSE:L2
  
  L2_BART <- mean( I(a.test$D6[1:NT]<BART)* (a.test$D6[1:NT]-BART)^2 + I(a.test$D6[1:NT] >= BART)* (a.test$D6[1:NT]-BART)^2)
  L2_NET<- mean( I(a.test$D6[1:NT]<NET)* (a.test$D6[1:NT]-NET)^2 +I(a.test$D6[1:NT] >= NET)* (a.test$D6[1:NT]-NET)^2)
  L2_RR <- mean( I(a.test$D6[1:NT]<RR)* (a.test$D6[1:NT]-RR)^2 +I(a.test$D6[1:NT] >= RR)* (a.test$D6[1:NT]-RR)^2)
  L2_GBM<- mean( I(a.test$D6[1:NT]<GBM)* (a.test$D6[1:NT]-GBM)^2 +I(a.test$D6[1:NT] >= GBM)* (a.test$D6[1:NT]-GBM)^2)
  L2_GAM <- mean( I(a.test$D6[1:NT]<GAM)* (a.test$D6[1:NT]-GAM)^2 +I(a.test$D6[1:NT] >= GAM)* (a.test$D6[1:NT]-GAM)^2)
  L2_QR_nonlin<- mean( I(a.test$D6[1:NT]<QR_nonlin)* (a.test$D6[1:NT]-QR_nonlin)^2 +I(a.test$D6[1:NT] >= QR_nonlin)* (a.test$D6[1:NT]-QR_nonlin)^2)
  L2_LM <- mean( I(a.test$D6[1:NT]<LM)* (a.test$D6[1:NT]-LM)^2 +I(a.test$D6[1:NT] >= LM)* (a.test$D6[1:NT]-LM)^2)
  L2_RF <- mean( I(a.test$D6[1:NT]<RF)* (a.test$D6[1:NT]-RF)^2 +I(a.test$D6[1:NT] >= RF)* (a.test$D6[1:NT]-RF)^2)
  
  L <- c( "L1","L2","BIAS" )
  TAU <- c( rep(0.5,3) )
  
  BARTS <- c( L1_BART,L2_BART ,bias_BART)
  NETS <- c( L1_NET,L2_NET ,bias_NET )
  RRS <- c( L1_RR,L2_RR ,bias_RR )
  RFS <- c( L1_RF,L2_RF ,bias_RF )
  GBMS <- c( L1_GBM,L2_GBM ,bias_GBM )
  GAMS <- c( L1_GAM,L2_GAM ,bias_GAM) 
  QR_nonlins <- c( L1_QR_nonlin,L2_QR_nonlin ,bias_QR_nonlin)
  LMs <- c( L1_LM,L2_LM ,bias_LM) 
  
  TABLE_MEDIANS_GLASLOCAL <- data.frame(Estimator = L,
                                        TAU       = TAU, 
                                        NET       = NETS,
                                        RR        = RRS,
                                        GBM       = GBMS,
                                        GAM       = GAMS,       
                                        QR_SP     = QR_nonlins,
                                        LM        = LMs, 
                                        RF        = RFS,
                                        BART      = BARTS)
  
  TABLE_MEDIANS_GLASLOCAL[,-c(1:2)] <- round(TABLE_MEDIANS_GLASLOCAL[,-c(1:2)],4) 
  write.csv(data.frame(t(TABLE_MEDIANS_GLASLOCAL)),"TABLE_MEDIANS_GLASLOCAL.csv")
  return(data.frame(t(TABLE_MEDIANS_GLASLOCAL)))
}

