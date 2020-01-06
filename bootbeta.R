 
bootbeta<- function(f0,taus,data,B)
{   
  fit <- lm(f0,data=data)
  Xp <- model.matrix(fit,data=data)
  n <- nrow(data)
  p <- ncol(Xp)
  BETAHAT<- array(NA,dim=c(B,length(taus),p))
  BETATILDE <- array(NA,dim=c(B,p))
  SIGMA2 <- NULL
  XpXinv <- array(NA,dim=c(B,p,p))
  bsobject <- list()
  bsobjectboot <- list()
  
  for(b in 1:B)
  { 
    newdatab <-  data[sample(nrow(data), dim(data)[1],replace=TRUE),]
    for(j in 1:length(taus)) {
      tau <- taus[j]  
      fitrq <-  rq(f0,tau=tau,data = newdatab)
      bsobject[[j]]<- fitrq
    } 
    bsobjectboot[[b]] <- bsobject
    fitlm <- lm(f0,data=newdatab)    
    BETATILDE[b,] <- summary(fitlm)$coef[,1]
    betatilde <- summary(fitlm)$coef[,1]
    SIGMA2[b] <- sigma(fitlm)^2   
    X <- model.matrix(fitlm)
    XpXinv[b,,] <- ginv(t(X)%*%X)
    
  }
  
  return(list( fitobjectboot= bsobjectboot,
               BETATILDE = BETATILDE,
               sigmahat2 = SIGMA2, 
               XpXinv = XpXinv,
               taus= taus,
               B=B, 
               f0=f0
  )
  )
}



