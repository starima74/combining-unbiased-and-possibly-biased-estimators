#############################################################################
# This function post-processes the simulations results to make summary tables 1 and 2
#############################################################################


post_proc <- function(X.1,
                      X.2,
                      eps){
  
  output_ <- data.frame(matrix(NA,nrow=nsims,ncol=108))
  
  ts<- array(NA,dim=c(nrow(newdata),length(taus),7))
  
  output<- as.matrix(output_)  
  
  for(s in 1:nsims){
    N <- 400
    B <- 200
    Betas <- c(-1,1)  
    Y <-  Betas[1]*X.1 + Betas[2]*X.2 + eps
    data <- data.frame(Y = Y, X.1 = X.1,X.2 = X.2)
    f0<- formula(Y ~ X.1 + X.2 )
    betas <- bootbeta(f0=f0, taus=taus, data=data, B=B)
    theta <- predvalK(data, newdata, betas, bias=1)
    
    for(i in 1:nrow(newdata))
    {
      ts[i,,] <- as.matrix(theta$thetas[[i]],nrow=length(taus),ncol=7)
    }
    #ensemble : QRA
    
    ens <- ensemble(f0,taus,data,newdata) 
    
    #Random forest : QRF
    
    xtrain <- data.frame( X.1 = X.1,X.2 = X.2)
    ytrain <- Y 
    qrf <- quantregForest(x=xtrain,y=ytrain)
    RF <- predict(qrf,newdata,what=taus)
    
    output[s,1:84] <-  c(ts) 
    output[s,85:96] <- c(ens$ENS_pred)
    output[s,97:108] <- c(RF)
    
  }
  return(output)
}
