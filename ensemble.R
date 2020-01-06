
#Ensemble method

ensemble <- function(f0,taus,data,newdata){
  data$predhat <- NULL;  data$predtilde <- NULL
  ENS_pred <- matrix(NA,nrow=nrow(newdata),ncol=length(taus))
  ENS_predcov <- matrix(NA,nrow=nrow(newdata),ncol=length(taus))
  
  fitlm <- lm(f0,data=data)
  betatilde <- summary(fitlm)$coef[,1]
  sigmahat2 <- sigma(fitlm)^2
  X <- model.matrix(fitlm)
  XpXinv <- ginv(t(X)%*%X)   
  
  newdata$Y <- 1
  Xval <- model.matrix(fitlm,data=newdata)
  n <- nrow(data)
  p <- ncol(Xval)
  
  for(j in 1:length(taus)) {
    tau <- taus[j]  
    fitrq <-  rq(f0,tau=tau,data = data)
    data$predhat  <-  as.numeric(fitrq$fitted.values)
    newdata$predhat <- predict(fitrq,newdata=newdata)
    
    y0hat <-   X%*%as.matrix(betatilde)     
    critval <- qt(tau,n-p)     
    for(i in 1:nrow(X)){
      data$predtilde[i] <- y0hat[i] + critval*sqrt(sigmahat2*(1+t(X[i,])%*%XpXinv%*%as.matrix(X[i,])))
    }
    
    y0hat_new <-   Xval%*%as.matrix(betatilde)       
    for(n in 1:nrow(Xval)){
      newdata$predtilde[n] <- y0hat_new[n] + critval*sqrt(sigmahat2*(1+t(Xval[n,])%*%XpXinv%*%as.matrix(Xval[n,])))
    }
    
    ENS_pred[,j] <- predict(rq( Y~predhat+predtilde,tau=tau,data=data),newdata=newdata) 
    #ENS_predcov[,j] <- predict(rq(Y~ X.1 + X.2 +predhat + predtilde ,tau=tau,data=data),newdata=newdata)  
    
  }  
  return(list (ENS_pred = ENS_pred))
}  