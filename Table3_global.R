#######################################################################################
# This function reproduces results for global methods presented in table 3, top
#####################################################################################


Table3_global <- function(){
  f_full <- formula(paste("D6 ~ A1C0  + age  ", sep=""))
  
  taus <- c(0.5)
  #Data 
  
  a.train<- read.csv("a.train.csv")
  a.test <- read.csv("a.test.csv") 
  pred_vars <- c("A1C0","age") 
  
  #BART
  
  ndpost = 1000
  nskip = 100  
  a.train$newcurr_met <- as.numeric(a.train$newcurr_met)
  a.train$neworder_insulin <- as.numeric(a.train$neworder_insulin)
  
  a.test$newcurr_met <- as.numeric(a.test$newcurr_met)
  a.test$neworder_insulin <- as.numeric(a.test$neworder_insulin)
  
  dat <- a.train[,c("D6",pred_vars)] 
  x.train <- a.train[,pred_vars]
  x.test <- a.test[,pred_vars]
  y.train <-  dat[,"D6"] 
  
  bart = wbart( x.train=x.train, x.test=x.test,y.train=y.train,ndpost=ndpost, nskip=nskip)   
  BART =  bart$yhat.test.mean 
  a.train$predBART = bart$yhat.train.mean  
  
  #Elastic net 
  fitNet <- train(D6~.,data=dat,method = "glmnet") 
  NET <- predict(fitNet,x.test,submodels=NULL)
  a.train$predNET <- predict(fitNet,x.train,submodels=NULL)
  
  #Ridge Regression  
  fitRR <- train(D6~.,data=dat,method = "ridge") 
  RR <- predict(fitRR,x.test,submodels=NULL)
  a.train$predRR <- predict(fitRR,x.train,submodels=NULL)
  
  # GBM
  gbmGrid <- expand.grid(interaction.depth=c(1,5,9),n.trees=(1:50)*50,
                         shrinkage=0.1, n.minobsinnode=10)
  fitControl <- trainControl(method="cv",number=10) 
  fitGBM <- train(D6~.,data=dat,method = "gbm",verbose=F,
                  trControl=fitControl,tuneGrid=gbmGrid) 
  GBM <- predict(fitGBM,x.test,submodels=NULL)
  a.train$predGBM <- predict(fitGBM,x.train,submodels=NULL)
  
  #ADDITIONAL METHODS
  
  #RANDOM FOREST
  fitRF <- train(x.train,y.train,method = "rf",importance=T)
  RF <- predict(fitRF,x.test,submodels=NULL)
  a.train$predRF <- predict(fitRF,x.train,submodels=NULL)
  
  #Generalized additive model 
  fitgam <- gam(f_full, data = a.train)
  GAM  <-   predict(fitgam,newdata = a.test)
  a.train$predGAM <- predict(fitgam,newdata = a.train)
  
  #QR_linear 
  fitlin <- rq(f_full, tau=taus, data=a.train)
  QR_lin <- predict(fitlin,newdata=a.test)
  a.train$predQR_lin <- predict( fitlin,newdata=a.train)
  
  
  #QR_nonparametric predictors 
  fitnonlin <- rq(D6~bs(A1C0,15)+bs(age,15),tau=taus, data=a.train)
  QR_nonlin <- predict(fitnonlin,newdata=a.test)
  a.train$predQR_nonlin <- predict(fitnonlin,newdata=a.train)
  
  #LM with linear predictors
  
  y0hat <- predict(lm(f_full,data=a.train), newdata = a.test, se=T)
  LM <- y0hat$fit + qnorm(taus) * y0hat$se.fit
  
  y0hat_p <- predict(lm(f_full,data=a.train), newdata = a.train, se=T)
  a.train$predLM <- y0hat_p$fit + qnorm(taus) * y0hat_p$se.fit
  
  ## Superlearner
  set.seed(1)
  
  SL.library <- c("SL.gam", "SL.gam.3", "SL.gam.4",
                  "SL.gam.5", "SL.gbm.1","SL.gbm", "SL.glm",
                  "SL.glmnet", "SL.glmnet.0.25","SL.glmnet.alpha.0.5", "SL.glmnet.0.75",
                  "SL.polymars", "SL.randomForest","SL.ridge", "SL.svm",
                  "SL.bayesglm", "SL.step","SL.step.interaction","SL.bart")
  SL.library <- c("SL.gam", "SL.gbm", "SL.glm",
                  "SL.glmnet","SL.polymars", "SL.randomForest","SL.ridge", "SL.svm",
                  "SL.bayesglm", "SL.step","SL.step.interaction")
  
  sl = SuperLearner(Y=y.train, X=x.train, family = gaussian(),
                    SL.library = SL.library)
  SL <- predict(sl, newdata = x.test)$pred
  a.train$predSL <- predict(sl, newdata = x.train)$pred
   
######################################################################
# Estimate bias, L1 and L2 measures for the above estimators
#########################################################################

#Bias
  
  bias_BART <- mean( I(a.test$D6<BART)*(a.test$D6-BART) + I(a.test$D6 >= BART)*(a.test$D6-BART))
  bias_NET<- mean( I(a.test$D6<NET)*(a.test$D6-NET) + I(a.test$D6 >= NET)*(a.test$D6-NET))
  bias_RR <- mean( I(a.test$D6<RR)*(a.test$D6-RR) + I(a.test$D6 >= RR)*(a.test$D6-RR))
  bias_GBM<- mean( I(a.test$D6<GBM)*(a.test$D6-GBM) + I(a.test$D6 >= GBM)*(a.test$D6-GBM))
  bias_GAM <- mean( I(a.test$D6<GAM)*(a.test$D6-GAM) + I(a.test$D6 >= GAM)*(a.test$D6-GAM))
  bias_QR_lin <- mean( I(a.test$D6<QR_lin)*(a.test$D6-QR_lin) + I(a.test$D6 >= QR_lin)*(a.test$D6-QR_lin))
  bias_QR_nonlin <- mean( I(a.test$D6<QR_nonlin)*(a.test$D6-QR_nonlin) + I(a.test$D6 >= QR_nonlin)*(a.test$D6-QR_nonlin))
  bias_LM <- mean( I(a.test$D6<LM)*(a.test$D6-LM) + I(a.test$D6 >= LM)*(a.test$D6-LM))
  bias_RF <- mean( I(a.test$D6<RF)*(a.test$D6-RF) + I(a.test$D6 >= RF)*(a.test$D6-RF))
  bias_SL <- mean( I(a.test$D6<SL)*(a.test$D6-SL) + I(a.test$D6 >= SL)*(a.test$D6-SL))
  
#L1
  
  L1_BART <- mean( I(a.test$D6<BART)*abs(a.test$D6-BART) + I(a.test$D6 >= BART)*abs(a.test$D6-BART))
  L1_NET<- mean( I(a.test$D6<NET)*abs(a.test$D6-NET) + I(a.test$D6 >= NET)*abs(a.test$D6-NET))
  L1_RR <- mean( I(a.test$D6<RR)*abs(a.test$D6-RR) + I(a.test$D6 >= RR)*abs(a.test$D6-RR))
  L1_GBM<- mean( I(a.test$D6<GBM)*abs(a.test$D6-GBM) + I(a.test$D6 >= GBM)*abs(a.test$D6-GBM))
  L1_GAM <- mean( I(a.test$D6<GAM)*abs(a.test$D6-GAM) + I(a.test$D6 >= GAM)*abs(a.test$D6-GAM))
  L1_QR_lin <- mean( I(a.test$D6<QR_lin)*abs(a.test$D6-QR_lin) + I(a.test$D6 >= QR_lin)*abs(a.test$D6-QR_lin))
  L1_QR_nonlin <- mean( I(a.test$D6<QR_nonlin)*abs(a.test$D6-QR_nonlin) + I(a.test$D6 >= QR_nonlin)*abs(a.test$D6-QR_nonlin))
  L1_LM <- mean( I(a.test$D6<LM)*abs(a.test$D6-LM) + I(a.test$D6 >= LM)*abs(a.test$D6-LM))
  L1_RF <- mean( I(a.test$D6<RF)*abs(a.test$D6-RF) + I(a.test$D6 >= RF)*abs(a.test$D6-RF))
  L1_SL <- mean( I(a.test$D6<SL)*abs(a.test$D6-SL) + I(a.test$D6 >= SL)*abs(a.test$D6-SL))
  
  
#MSE:L2
  
  L2_BART <- mean( I(a.test$D6<BART)*(a.test$D6-BART)^2 + I(a.test$D6 >= BART)*(a.test$D6-BART)^2)
  L2_NET<- mean( I(a.test$D6<NET)*(a.test$D6-NET)^2 + I(a.test$D6 >= NET)*(a.test$D6-NET)^2)
  L2_RR <- mean( I(a.test$D6<RR)*(a.test$D6-RR)^2 + I(a.test$D6 >= RR)*(a.test$D6-RR)^2)
  L2_GBM<- mean( I(a.test$D6<GBM)*(a.test$D6-GBM)^2 + I(a.test$D6 >= GBM)*(a.test$D6-GBM)^2)
  L2_GAM <- mean( I(a.test$D6<GAM)*(a.test$D6-GAM)^2 + I(a.test$D6 >= GAM)*(a.test$D6-GAM)^2)
  L2_QR_lin <- mean( I(a.test$D6<QR_lin)*(a.test$D6-QR_lin)^2 + I(a.test$D6 >= QR_lin)*(a.test$D6-QR_lin)^2)
  L2_QR_nonlin<- mean( I(a.test$D6<QR_nonlin)*(a.test$D6-QR_nonlin)^2 + I(a.test$D6 >= QR_nonlin)*(a.test$D6-QR_nonlin)^2)
  L2_LM <- mean( I(a.test$D6<LM)*(a.test$D6-LM)^2 + I(a.test$D6 >= LM)*(a.test$D6-LM)^2)
  L2_RF <- mean( I(a.test$D6<RF)*(a.test$D6-RF)^2 + I(a.test$D6 >= RF)*(a.test$D6-RF)^2)
  L2_SL <- mean( I(a.test$D6<SL)*(a.test$D6-SL)^2 + I(a.test$D6 >= SL)*(a.test$D6-SL)^2)
  
  L <- c( "L1","L2","BIAS" )
  TAU <- c( rep(0.5,3) ) 
  BARTS <- c( L1_BART,L2_BART ,bias_BART)
  NETS <- c( L1_NET,L2_NET ,bias_NET )
  RRS <- c( L1_RR,L2_RR ,bias_RR )
  RFS <- c( L1_RF,L2_RF ,bias_RF )
  GBMS <- c( L1_GBM,L2_GBM ,bias_GBM ) 
  GAMS <- c( L1_GAM,L2_GAM ,bias_GAM)
  QR_lins <- c( L1_QR_lin,L2_QR_lin ,bias_QR_lin)
  QR_nonlins <- c( L1_QR_nonlin,L2_QR_nonlin ,bias_QR_nonlin)
  LMs <- c( L1_LM,L2_LM ,bias_LM)
  SLs <- c( L1_SL,L2_SL ,bias_SL)
  
  TABLE_MEDIANS_GLOBAL <- data.frame(Estimator = L,
                                     TAU       = TAU, 
                                     NET       = NETS,
                                     RR        = RRS,
                                     GBM       = GBMS,
                                     GAM       = GAMS,  
                                     QR        = QR_lins,
                                     QR_SP     = QR_nonlins,
                                     LM        = LMs, 
                                     RF        = RFS,
                                     BART      = BARTS,
                                     SL        = SLs) 
  
  TABLE_MEDIANS_GLOBAL[,-c(1:2)] <- round(TABLE_MEDIANS_GLOBAL[,-c(1:2)],4) 
  write.csv(data.frame(t(TABLE_MEDIANS_GLOBAL)),"TABLE_MEDIANS_GLOBAL.csv")
  return(data.frame(t(TABLE_MEDIANS_GLOBAL)))
}



