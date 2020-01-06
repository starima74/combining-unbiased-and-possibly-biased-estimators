###########################################################################
# This function takes in the booststrap samples from bootbeta
# and returns theta estimates and their MSEs
######################################################################
predvalK<- function(data,newdata,betas,bias=1)
{  
  f0 <- betas$f0  
  fit <- lm(betas$f0,data=data)
  newdata$Y <- 1
  Xval <- model.matrix(fit,data=newdata)
  B <- betas$B 
  taus <- betas$taus
  n <- nrow(data)
  p <- ncol(Xval)
  thetahat <- NULL
  thetatilde <- NULL;  betahat <- NULL;  betatilde <- NULL
  thetanewvar <- NULL;  thetanewmse <- NULL
  thetanewmse_c2 <- NULL;thetanewmse_c4 <- NULL;thetanewmse_c8 <- NULL
  THETAHAT <- array(NA,dim=c(B,length(taus)))
  THETATILDE<- array(NA,dim=c(B,length(taus)))
  YOHAT<- array(NA,dim=c(B,length(taus)))
  thetalist <- list()
  K11 <- array(NA,dim=c(nrow(newdata),length(taus)))
  K12  <- array(NA,dim=c(nrow(newdata),length(taus)))
  K22var <- array(NA,dim=c(nrow(newdata),length(taus)))
  K22mse  <- array(NA,dim=c(nrow(newdata),length(taus)))
  K22mse_c2  <- array(NA,dim=c(nrow(newdata),length(taus)))
  K22mse_c4  <- array(NA,dim=c(nrow(newdata),length(taus)))
  K22mse_c8  <- array(NA,dim=c(nrow(newdata),length(taus)))
  K22 <- array(NA,dim=c(nrow(newdata),length(taus)))
  newhat<- array(NA,dim=c(nrow(newdata),length(taus)))
  newtilde<- array(NA,dim=c(nrow(newdata),length(taus)))
  newvar<- array(NA,dim=c(nrow(newdata),length(taus)))
  newmse<- array(NA,dim=c(nrow(newdata),length(taus)))
  
  
  fitlm <- lm(f0,data=data)
  betatilde <- summary(fitlm)$coef[,1]
  sigmahat2 <- sigma(fitlm)^2
  X <- model.matrix(fitlm)
  XpXinv <- ginv(t(X)%*%X)  
  
  for(i in 1:nrow(Xval))
  {
    for(j in 1:length(taus)) {
      
      tau <- taus[j]  
      fitrq <-  rq(f0,tau=tau,data = data)
      thetahat[j] <-  predict(fitrq,newdata[i,])       
      y0hat <-   t(Xval[i,]) %*% betatilde      
      critval <- qt(tau,n-p)     
      
      thetatilde[j] <- y0hat + critval*sqrt(sigmahat2*(1+t(Xval[i,])%*%XpXinv%*%Xval[i,])) 
      
      for(b in 1:B)
      {
        THETAHAT[b,j] <-   predict(betas$fitobjectboot[[b]][[j]],newdata[i,]) 
        YOHAT[b]   <-   t(Xval[i,]) %*% betas$BETATILDE[b,] 
        THETATILDE[b,j] <- as.matrix( YOHAT[b]+  critval*sqrt(betas$sigmahat2[b]*(1+ t(Xval[i,])%*%betas$XpXinv[b,,]%*%Xval[i,])))
      }  
      
      delta <- THETAHAT[,j] - THETATILDE[,j]
      p  <- pchisq(mean(delta)^2/var(delta),df=1)
       
      c2 <- n^{-1/2}
      c4 <- n^{-1/4}
      c8 <- n^{-1/8}
      bias_c2 <- bias*c2
      bias_c4 <- bias*c4
      bias_c8 <- bias*c8
      
      K11[i,j] <- var(THETAHAT[,j])
      K12[i,j] <-  cov(THETAHAT[,j], THETAHAT[,j]-THETATILDE[,j])
      K22var[i,j] <- var(THETAHAT[,j] - THETATILDE[,j])  
      K22mse[i,j] <- var(THETAHAT[,j] - THETATILDE[,j]) + bias*(thetahat[j]-thetatilde[j])^2
      K22mse_c2[i,j] <- var(THETAHAT[,j] - THETATILDE[,j]) + bias_c2*(thetahat[j]-thetatilde[j])^2
      K22mse_c4[i,j] <- var(THETAHAT[,j] - THETATILDE[,j]) + bias_c4*(thetahat[j]-thetatilde[j])^2
      K22mse_c8[i,j] <- var(THETAHAT[,j] - THETATILDE[,j]) + bias_c8*(thetahat[j]-thetatilde[j])^2
      K22[i,j]<- var(THETATILDE[,j])
      
      K22varinv <- ginv(K22var[i,j])
      K22mseinv <- ginv(K22mse[i,j])
      K22mse_c2inv <- ginv(K22mse_c2[i,j])
      K22mse_c4inv <- ginv(K22mse_c4[i,j])
      K22mse_c8inv <- ginv(K22mse_c8[i,j])
      thetanewvar[j] <- thetahat[j] - K12[i,j]*K22varinv*(thetahat[j]-thetatilde[j])
      thetanewmse[j] <- thetahat[j] - K12[i,j]*K22mseinv*(thetahat[j]-thetatilde[j])
      thetanewmse_c2[j] <- thetahat[j] - K12[i,j]*K22mse_c2inv*(thetahat[j]-thetatilde[j])
      thetanewmse_c4[j] <- thetahat[j] - K12[i,j]*K22mse_c4inv*(thetahat[j]-thetatilde[j])
      thetanewmse_c8[j] <- thetahat[j] - K12[i,j]*K22mse_c8inv*(thetahat[j]-thetatilde[j])
      
      #MSEs
      newhat[i,j]<- K11[i,j]#+ (thetahat[j] - mean(THETAHAT[,j]))^2 
      newtilde[i,j] <- K22[i,j] + mean((THETAHAT[,j]-THETATILDE[,j]))^2-K22var[i,j]
      newvar[i,j] <- K11[i,j] - K12[i,j]*K22varinv*t(K12[i,j])  
      newmse[i,j] <- newhat[i,j] - K12[i,j]*K22mseinv*t(K12[i,j])
    }
    
    theta <- data.frame(thetahat =  thetahat,
                        thetatilde = thetatilde,
                        thetanewvar = thetanewvar,
                        thetanewmse = thetanewmse,
                        thetanewmse_c2 = thetanewmse_c2,
                        thetanewmse_c4 = thetanewmse_c4,
                        thetanewmse_c8 = thetanewmse_c8
    ) 
    thetalist[[i]]<- theta
    
  }
  
  
  return(list(thetas = thetalist, 
              newhat=newhat,
              newtilde=newtilde,
              newvar=newvar,
              newmse=newmse))
  
}