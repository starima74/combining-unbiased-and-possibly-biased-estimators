#########################################################################
# This function implements local methods(15%) given in table 4
#########################################################################

Table4_local <- function(){
  
  ## distance 
  
  dist <- function(all, one, w=NA) {
    all <- all[,names(one)]
    ret = as.vector(rep(0,nrow(all)))
    Vs <- diag(var(all))
    # print(Vs)
    if(is.na(sum(w))) w = rep(1,ncol(all))
    #for(i in 1:nrow(all))
    for(d in 1:ncol(one))
      ret <- ret + as.numeric(w[d]) * ( as.vector(all[,d])  - as.numeric(one[d]))^2 / as.numeric(Vs[d]) 
    #sqrt(ret)
    ret
  }
  
  a.train <- read.csv("a.train.csv")
  a.test <- read.csv("a.test.csv")
  
  taus <- c(0.25,0.75)
  
  
  HAT <- matrix(NA,nrow=nrow(a.test),ncol = length(taus))
  TILDE <- matrix(NA,nrow=nrow(a.test),ncol = length(taus))
  HAT2 <- matrix(NA,nrow=nrow(a.test),ncol = length(taus))
  MMSE <- matrix(NA,nrow=nrow(a.test),ncol = length(taus))
  MMSE1 <- matrix(NA,nrow=nrow(a.test),ncol = length(taus))
  MMSE2 <- matrix(NA,nrow=nrow(a.test),ncol = length(taus))
  MMSE3 <- matrix(NA,nrow=nrow(a.test),ncol = length(taus))
  MMSE4 <- matrix(NA,nrow=nrow(a.test),ncol = length(taus))
  MMSE5 <- matrix(NA,nrow=nrow(a.test),ncol = length(taus))
  MMSE6 <- matrix(NA,nrow=nrow(a.test),ncol = length(taus))
  MVAR <- matrix(NA,nrow=nrow(a.test),ncol = length(taus))
   
  
  outcome <- "D6"
  f_full <- formula(paste(outcome," ~ A1C0 + newcurr_met + neworder_met + age + White", sep=""))
  f2 <- formula(paste(outcome," ~ A1C0 + age", sep=""))
  f3 <- formula(paste(outcome," ~ A1C0 + age + I(age^2)", sep=""))
  
  
  theta.f <- function(d) {
    d <- d[is.na(d[,outcome])==FALSE,]
    tmplm <- lm(f3,data=d)
    tmpeps <- min(svd((t(model.matrix(tmplm)) %*% (model.matrix(tmplm))))$d)
    if(tmpeps>0.0000001) ret <- predict(rq(f3, tau = tau, data = d), newdata = newperson)
    if(tmpeps<=0.0000001) ret <- NA
    ret
  }
  
  lm.quantile <- function(d) {
    d <- d[is.na(d[,outcome])==FALSE,]
    y0hat <- predict(lm(f2,data=d), 
                     newdata = newperson, se=T)
    y0hat$fit + qnorm(tau) * y0hat$se.fit
  }
  
  qr.quantile <- function(d) {
    d <- d[is.na(d[,outcome])==FALSE,]
    tmplm <- lm(f2,data=d)
    tmpeps <- min(svd((t(model.matrix(tmplm)) %*% (model.matrix(tmplm))))$d)
    if(tmpeps>0.0000001) ret <- predict(rq(f2, tau = tau, data = d), newdata = newperson)
    if(tmpeps<=0.0000001) ret <- NA
    ret
  }
  
  NT <- nrow(a.test)
  #To predict using a small subset of the test data set, uncomment 
  #NT <- 500
  
  for(id in 1:nrow(a.test[1:NT,])){ 
    for(t in 1:length(taus)){ 
    tau <- taus[t]
    if(id/50  == floor(id/50)) print( paste("id=",id,",tau=",0.5))
    newperson <- a.test[id, ]
    one <- a.test[id,c("age", "A1C0")]
    dists <- dist(a.train,one, w = c(1, 3))
    a.train$dists <- dists
    cutoff <- quantile(dists, probs = 0:20/20)[4]
    ind <- a.train$dists<cutoff
    
    Add.Info.Means <- list()
    Add.Info.Vars <- list()
    Add.Info.Functions <- list()
    Add.Info.Biases <- list()
    
    Add.Info.Means[[1]] <- lm.quantile(a.train[ind,])
    Add.Info.Vars[[1]] <- NA
    Add.Info.Functions[[1]] <- lm.quantile
    Add.Info.Biases[[1]] <- 1
    
    Add.Info.Means[[2]] <- qr.quantile(a.train[ind,])
    Add.Info.Vars[[2]] <- NA
    Add.Info.Functions[[2]] <- qr.quantile
    Add.Info.Biases[[2]] <- 1
    
    Add.Info <- data.frame(Means = rep(NA,2),
                           Vars = rep(NA,2),
                           Functions = rep(NA,2),
                           Biases = rep(NA,2))
    
    Add.Info$Means = Add.Info.Means
    Add.Info$Vars  = Add.Info.Vars
    Add.Info$Functions = Add.Info.Functions
    Add.Info$Biases = Add.Info.Biases
    n <- nrow(a.train[ind,])
    
    res <- MMSEnew(a.train[ind,], theta.f, 
                   Add.Info, nboots = 100, 
                   eig.cutoff = 1,
                   c1 = n^(-0.125),
                   c2 = n^(- 0.250),
                   c3 = n^(-0.5),
                   c4 = n^(-0.75),
                   c5 = n^(-1.0),
                   c6 = n^(-1.25))
    
    MMSE[id,t]  <- c(res$Theta.MMSE)
    MMSE1[id,t]  <- c(res$Theta.MMSE.c1)
    MMSE2[id,t]  <- c(res$Theta.MMSE.c2)
    MMSE3[id,t] <- c(res$Theta.MMSE.c3)
    MMSE4[id,t]  <- c(res$Theta.MMSE.c4)
    MMSE5[id,t]  <- c(res$Theta.MMSE.c5)
    MMSE6[id,t] <- c(res$Theta.MMSE.c6)
    MVAR[id,t] <- c(res$Theta.MVAR)
    HAT[id,t]   <- c(res$Theta.Hat)
    TILDE[id,t] <- Add.Info[1]$Means[[1]]
    HAT2[id,t] <- Add.Info[1]$Means[[2]]
    
  }
  }
  
  
  #Bias 
  
  bias_MMSE <- apply(I(a.test$D6[1:NT] < MMSE[1:NT,])* (a.test$D6[1:NT]-MMSE[1:NT,]) +
                      I(a.test$D6[1:NT] >= MMSE[1:NT,])*(a.test$D6[1:NT]-MMSE[1:NT,]),2,mean)
  bias_MMSE1 <- apply( I(a.test$D6[1:NT]<MMSE1[1:NT,])* (a.test$D6[1:NT]-MMSE1[1:NT,])+
                        I(a.test$D6[1:NT] >= MMSE1[1:NT,])*(a.test$D6[1:NT]-MMSE1[1:NT,]),2,mean)
  bias_MMSE2 <- apply( I(a.test$D6[1:NT]<MMSE2[1:NT,])*(a.test$D6[1:NT]-MMSE2[1:NT,]) +
                        I(a.test$D6[1:NT] >= MMSE2[1:NT,])*(a.test$D6[1:NT]-MMSE2[1:NT,]),2,mean)
  bias_MMSE3 <- apply( I(a.test$D6[1:NT]<MMSE3[1:NT,])*(a.test$D6[1:NT]-MMSE3[1:NT,]) +
                        I(a.test$D6[1:NT] >= MMSE3[1:NT,])*(a.test$D6[1:NT]-MMSE3[1:NT,]),2,mean)
  bias_MMSE4 <- apply( I(a.test$D6[1:NT]<MMSE4[1:NT,])*(a.test$D6[1:NT]-MMSE4[1:NT,]) +
                        I(a.test$D6[1:NT] >= MMSE4[1:NT,])*(a.test$D6[1:NT]-MMSE4[1:NT,]),2,mean)
  bias_MMSE5 <- apply( I(a.test$D6[1:NT]<MMSE5[1:NT,])*(a.test$D6[1:NT]-MMSE5[1:NT,]) +
                        I(a.test$D6[1:NT] >= MMSE5[1:NT,])*(a.test$D6[1:NT]-MMSE5[1:NT,]),2,mean)
  bias_MMSE6 <- apply( I(a.test$D6[1:NT]<MMSE6[1:NT,])*(a.test$D6[1:NT]-MMSE6[1:NT,]) +
                        I(a.test$D6[1:NT] >= MMSE6[1:NT,])* (a.test$D6[1:NT]-MMSE6[1:NT,]),2,mean)
  bias_MVAR <- apply( I(a.test$D6[1:NT]<MVAR[1:NT,])* (a.test$D6[1:NT]-MVAR[1:NT,]) +
                       I(a.test$D6[1:NT] >= MVAR[1:NT,])* (a.test$D6[1:NT]-MVAR[1:NT,]),2,mean)
  bias_HAT <- apply( I(a.test$D6[1:NT]<HAT[1:NT,])* (a.test$D6[1:NT]-HAT[1:NT,]) +
                      I(a.test$D6[1:NT] >= HAT[1:NT,])*(a.test$D6[1:NT]-HAT[1:NT,]),2,mean)
  bias_HAT2 <- apply( I(a.test$D6[1:NT]<HAT2[1:NT,])*(a.test$D6[1:NT]-HAT2[1:NT,]) +
                       I(a.test$D6[1:NT] >= HAT2[1:NT,])*(a.test$D6[1:NT]-HAT2[1:NT,]),2,mean)
  bias_TILDE <- apply( I(a.test$D6[1:NT]<TILDE[1:NT,])*(a.test$D6[1:NT]-TILDE[1:NT,])+
                        I(a.test$D6[1:NT] >= TILDE[1:NT,])* (a.test$D6[1:NT]-TILDE[1:NT,]),2,mean)
  
  
  
  #L1 
  
  L1_MMSE <- apply(I(a.test$D6[1:NT] < MMSE[1:NT,])* abs(a.test$D6[1:NT]-MMSE[1:NT,]) +
                    I(a.test$D6[1:NT] >= MMSE[1:NT,])*abs(a.test$D6[1:NT]-MMSE[1:NT,]),2,mean)
  L1_MMSE1 <- apply( I(a.test$D6[1:NT]<MMSE1[1:NT,])* abs(a.test$D6[1:NT]-MMSE1[1:NT,]) +
                      I(a.test$D6[1:NT] >= MMSE1[1:NT,])*abs(a.test$D6[1:NT]-MMSE1[1:NT,]),2,mean)
  L1_MMSE2 <- apply( I(a.test$D6[1:NT]<MMSE2[1:NT,])*abs(a.test$D6[1:NT]-MMSE2[1:NT,]) +
                      I(a.test$D6[1:NT] >= MMSE2[1:NT,])*abs(a.test$D6[1:NT]-MMSE2[1:NT,]),2,mean)
  L1_MMSE3 <- apply( I(a.test$D6[1:NT]<MMSE3[1:NT,])*abs(a.test$D6[1:NT]-MMSE3[1:NT,]) +
                      I(a.test$D6[1:NT] >= MMSE3[1:NT,])*abs(a.test$D6[1:NT]-MMSE3[1:NT,]),2,mean)
  L1_MMSE4 <- apply( I(a.test$D6[1:NT]<MMSE4[1:NT,])*abs(a.test$D6[1:NT]-MMSE4[1:NT,]) +
                      I(a.test$D6[1:NT] >= MMSE4[1:NT,])*abs(a.test$D6[1:NT]-MMSE4[1:NT,]),2,mean)
  L1_MMSE5 <- apply( I(a.test$D6[1:NT]<MMSE5[1:NT,])*abs(a.test$D6[1:NT]-MMSE5[1:NT,]) +
                      I(a.test$D6[1:NT] >= MMSE5[1:NT,])*abs(a.test$D6[1:NT]-MMSE5[1:NT,]),2,mean)
  L1_MMSE6 <- apply( I(a.test$D6[1:NT]<MMSE6[1:NT,])*abs(a.test$D6[1:NT]-MMSE6[1:NT,]) +
                      I(a.test$D6[1:NT] >= MMSE6[1:NT,])* (a.test$D6[1:NT]-MMSE6[1:NT,]),2,mean)
  L1_MVAR <- apply( I(a.test$D6[1:NT]<MVAR[1:NT,])* abs(a.test$D6[1:NT]-MVAR[1:NT,]) +
                     I(a.test$D6[1:NT] >= MVAR[1:NT,])* abs(a.test$D6[1:NT]-MVAR[1:NT,]),2,mean)
  L1_HAT <- apply( I(a.test$D6[1:NT]<HAT[1:NT,])* abs(a.test$D6[1:NT]-HAT[1:NT,]) +
                    I(a.test$D6[1:NT] >= HAT[1:NT,])*abs(a.test$D6[1:NT]-HAT[1:NT,]),2,mean)
  L1_HAT2 <- apply( I(a.test$D6[1:NT]<HAT2[1:NT,])*abs(a.test$D6[1:NT]-HAT2[1:NT,]) +
                     I(a.test$D6[1:NT] >= HAT2[1:NT,])*abs(a.test$D6[1:NT]-HAT2[1:NT,]),2,mean)
  L1_TILDE <- apply( I(a.test$D6[1:NT]<TILDE[1:NT,])*abs(a.test$D6[1:NT]-TILDE[1:NT,])+
                      I(a.test$D6[1:NT] >= TILDE[1:NT,])* abs(a.test$D6[1:NT]-TILDE[1:NT,]),2,mean)
  
  
  
  # L2  
  
  
  L2_MMSE <- apply(I(a.test$D6[1:NT] < MMSE[1:NT,])* (a.test$D6[1:NT]-MMSE[1:NT,])^2 +
                    I(a.test$D6[1:NT] >= MMSE[1:NT,])* (a.test$D6[1:NT]-MMSE[1:NT,])^2,2,mean)
  L2_MMSE1 <- apply( I(a.test$D6[1:NT]<MMSE1[1:NT,])* (a.test$D6[1:NT]-MMSE1[1:NT,])^2 +
                      I(a.test$D6[1:NT] >= MMSE1[1:NT,])* (a.test$D6[1:NT]-MMSE1[1:NT,])^2,2,mean)
  L2_MMSE2 <- apply( I(a.test$D6[1:NT]<MMSE2[1:NT,])* (a.test$D6[1:NT]-MMSE2[1:NT,])^2 +
                      I(a.test$D6[1:NT] >= MMSE2[1:NT,])* (a.test$D6[1:NT]-MMSE2[1:NT,])^2,2,mean)
  L2_MMSE3 <- apply( I(a.test$D6[1:NT]<MMSE3[1:NT,])* (a.test$D6[1:NT]-MMSE3[1:NT,])^2 +
                      I(a.test$D6[1:NT] >= MMSE3[1:NT,])* (a.test$D6[1:NT]-MMSE3[1:NT,])^2,2,mean)
  L2_MMSE4 <- apply( I(a.test$D6[1:NT]<MMSE4[1:NT,])* (a.test$D6[1:NT]-MMSE4[1:NT,])^2 +
                      I(a.test$D6[1:NT] >= MMSE4[1:NT,])* (a.test$D6[1:NT]-MMSE4[1:NT,])^2,2,mean)
  L2_MMSE5 <- apply( I(a.test$D6[1:NT]<MMSE5[1:NT,])* (a.test$D6[1:NT]-MMSE5[1:NT,])^2 +
                      I(a.test$D6[1:NT] >= MMSE5[1:NT,])* (a.test$D6[1:NT]-MMSE5[1:NT,])^2,2,mean)
  L2_MMSE6 <- apply( I(a.test$D6[1:NT]<MMSE6[1:NT,])* (a.test$D6[1:NT]-MMSE6[1:NT,])^2 +
                      I(a.test$D6[1:NT] >= MMSE6[1:NT,])* (a.test$D6[1:NT]-MMSE6[1:NT,])^2,2,mean)
  L2_MVAR <- apply( I(a.test$D6[1:NT]<MVAR[1:NT,])* (a.test$D6[1:NT]-MVAR[1:NT,])^2 +
                     I(a.test$D6[1:NT] >= MVAR[1:NT,])* (a.test$D6[1:NT]-MVAR[1:NT,])^2,2,mean)
  L2_HAT <- apply( I(a.test$D6[1:NT]<HAT[1:NT,])* (a.test$D6[1:NT]-HAT[1:NT,])^2 +
                    I(a.test$D6[1:NT] >= HAT[1:NT,])* (a.test$D6[1:NT]-HAT[1:NT,])^2,2,mean)
  L2_HAT2 <- apply( I(a.test$D6[1:NT]<HAT2[1:NT,])* (a.test$D6[1:NT]-HAT2[1:NT,])^2 +
                     I(a.test$D6[1:NT] >= HAT2[1:NT,])* (a.test$D6[1:NT]-HAT2[1:NT,])^2,2,mean)
  L2_TILDE <- apply( I(a.test$D6[1:NT]<TILDE[1:NT,])* (a.test$D6[1:NT]-TILDE[1:NT,])^2+
                      I(a.test$D6[1:NT] >= TILDE[1:NT,])* (a.test$D6[1:NT]-TILDE[1:NT,])^2,2,mean)
  
  
  L <- c( "L1","L2","BIAS" )
  TAU <- rep(taus,3)
    
  
  MMSES <- c( L1_MMSE,L2_MMSE ,bias_MMSE)
  MMSE1S <- c( L1_MMSE1,L2_MMSE1 ,bias_MMSE1) 
  MMSE2S <- c( L1_MMSE2,L2_MMSE2 ,bias_MMSE2)  
  MMSE3S <- c( L1_MMSE3,L2_MMSE3 ,bias_MMSE3)
  MMSE4S <- c( L1_MMSE4,L2_MMSE4 ,bias_MMSE4)
  MMSE5S <- c( L1_MMSE5,L2_MMSE5 ,bias_MMSE5)
  MMSE6S <- c( L1_MMSE6,L2_MMSE6 ,bias_MMSE6)
  MVARS <- c( L1_MVAR,L2_MVAR ,bias_MVAR)
  HATS <- c( L1_HAT,L2_HAT ,bias_HAT)
  HAT2S <- c( L1_HAT2,L2_HAT2 ,bias_HAT2)
  TILDES <- c( L1_TILDE,L2_TILDE ,bias_TILDE)
  
  ALLTAUS_LOCAL <- data.frame(Estimator = L,
                                    TAU       = TAU, 
                                    MMSE      = MMSES,
                                    MMSE1     = MMSE1S,
                                    MMSE2     = MMSE2S,
                                    MMSE3     = MMSE3S,
                                    MMSE4     = MMSE4S,
                                    MMSE5     = MMSE5S,
                                    MMSE6     = MMSE6S,
                                    HAT       = HATS,
                                    HAT2      = HAT2S,
                                    TILDE     = TILDES)
  
ALLTAUS_LOCAL <- ALLTAUS_LOCAL[order(TAU),]

  ALLTAUS_LOCAL[,-c(1:2)] <- round(ALLTAUS_LOCAL[,-c(1:2)],4) 
  write.csv(data.frame(t(ALLTAUS_LOCAL)),"ALLTAUS_LOCAL.csv")
  return(data.frame(t(ALLTAUS_LOCAL)))
}

