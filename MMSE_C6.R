

############ Function (MMSE internal) for local estimators #############
MMSEnew = function(dd, theta.f, Add.Inf,
                nboots = 500, eig.cutoff = 0.9, c1 = NA, 
                c2=NA, c3=NA, c4 = NA, c5 = NA, c6 =NA, 
                ...) {
  require(Matrix)
  p <- nrow(Add.Inf)
  n <- nrow(dd)
  thetahat <- theta.f(dd, ...)
  
  thetahat.length <- length(thetahat)
  betahat <- c()
  for(i in 1:p) {
    betahat.add <- Add.Inf$Functions[[i]](dd) - thetahat
    betahat <- c(betahat,betahat.add)
  }
  
  betahat.length <- length(betahat)
  bootres <- matrix(NA, nrow = nboots, 
                    ncol = thetahat.length + betahat.length)
 
  
  for(i in 1:nboots) {
    ind <- sample(1:n, replace = T)
    boot.hat <- theta.f(dd[ind,], ...)
    boot.hat1 <- theta.f(dd[ind,], ...)
    for(j in 1:p) 
      boot.hat <- c(boot.hat, 
                    Add.Inf$Functions[[j]](dd[ind,]) - 
                      boot.hat1)
      bootres[i,] <- boot.hat
  }
  
  K = var(bootres,use="pairwise.complete.obs")#/nrow(dd)
  ind1 <- 1:thetahat.length
  ind2 <- (thetahat.length+1):(betahat.length + 
                                 thetahat.length)
  K11 <- K[ind1,ind1]
  K12 <- K[ind1,ind2]
  K21 <- K[ind2,ind1]
  Biases <- colMeans(bootres,na.rm=TRUE)
  Biases <- Biases[ind2]
  Biases2 <- Biases %o% Biases
  Biases2[Add.Inf$Biases == 0, ] <- 0
  Biases2[, Add.Inf$Biases == 0] <- 0
  
  ### MMSE
  K22 <- K[ind2,ind2] + Biases2
  SVD_var <- svd(K22)
  ind.svd <- rep(1,length(SVD_var$d))
  if(length(SVD_var$d)>1)
    for(i in 2:length(SVD_var$d))
      ind.svd[i] <- sum(SVD_var$d[1:(i-1)])/sum(SVD_var$d) < eig.cutoff
  diag_elem <- rep(0,length(SVD_var$d))
  diag_elem[SVD_var$d>0 & ind.svd == 1] <- 
    1/SVD_var$d[SVD_var$d>0 & ind.svd == 1]
  if(length(SVD_var$d>0 & ind.svd == 1) == 1)  
    Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
  if(length(SVD_var$d>0 & ind.svd == 1) > 1)  
    Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
  thetahat.MMSE = matrix(thetahat,ncol=1) - K12 %*% Kinv %*% 
    matrix(betahat,ncol=1)
  thetahat.MMSE.Var = K11 - K12 %*% Kinv %*% K21
  K22inv_MMSE <- Kinv
  EVPROP_MMSE <- SVD_var$d/sum(SVD_var$d)

  ### MMSE(c1)
  thetahat.MMSE.c1 = NA
  thetahat.MMSE.Var.c1 = NA
  if(!is.na(c1)) {
  K22 <- K[ind2,ind2] + c1*Biases2
  SVD_var <- svd(K22)
  ind.svd <- rep(1,length(SVD_var$d))
  if(length(SVD_var$d)>1)
    for(i in 2:length(SVD_var$d))
      ind.svd[i] <- sum(SVD_var$d[1:(i-1)])/sum(SVD_var$d) < eig.cutoff
  diag_elem <- rep(0,length(SVD_var$d))
  diag_elem[SVD_var$d>0 & ind.svd == 1] <- 
    1/SVD_var$d[SVD_var$d>0 & ind.svd == 1]
  if(length(SVD_var$d>0 & ind.svd == 1) == 1)  
    Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
  if(length(SVD_var$d>0 & ind.svd == 1) > 1)  
    Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
  thetahat.MMSE.c1 = matrix(thetahat,ncol=1) - K12 %*% Kinv %*% 
    matrix(betahat,ncol=1)
  thetahat.MMSE.Var.c1 = K11 - K12 %*% Kinv %*% K21
  }
  
  ### MMSE(c2)
  thetahat.MMSE.c2 = NA
  thetahat.MMSE.Var.c2 = NA
  if(!is.na(c2)) {
    K22 <- K[ind2,ind2] + c2*Biases2
  SVD_var <- svd(K22)
  ind.svd <- rep(1,length(SVD_var$d))
  if(length(SVD_var$d)>1)
    for(i in 2:length(SVD_var$d))
      ind.svd[i] <- sum(SVD_var$d[1:(i-1)])/sum(SVD_var$d) < eig.cutoff
  diag_elem <- rep(0,length(SVD_var$d))
  diag_elem[SVD_var$d>0 & ind.svd == 1] <- 
    1/SVD_var$d[SVD_var$d>0 & ind.svd == 1]
  if(length(SVD_var$d>0 & ind.svd == 1) == 1)  
    Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
  if(length(SVD_var$d>0 & ind.svd == 1) > 1)  
    Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
  thetahat.MMSE.c2 = matrix(thetahat,ncol=1) - K12 %*% Kinv %*% 
    matrix(betahat,ncol=1)
  thetahat.MMSE.Var.c2 = K11 - K12 %*% Kinv %*% K21
  }
  
  ### MMSE(c3)
  thetahat.MMSE.c3 = NA
  thetahat.MMSE.Var.c3 = NA
  if(!is.na(c3)) {
    K22 <- K[ind2,ind2] + c3*Biases2
    SVD_var <- svd(K22)
    ind.svd <- rep(1,length(SVD_var$d))
    if(length(SVD_var$d)>1)
      for(i in 2:length(SVD_var$d))
        ind.svd[i] <- sum(SVD_var$d[1:(i-1)])/sum(SVD_var$d) < eig.cutoff
    diag_elem <- rep(0,length(SVD_var$d))
    diag_elem[SVD_var$d>0 & ind.svd == 1] <- 
      1/SVD_var$d[SVD_var$d>0 & ind.svd == 1]
    if(length(SVD_var$d>0 & ind.svd == 1) == 1)  
      Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
    if(length(SVD_var$d>0 & ind.svd == 1) > 1)  
      Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
    thetahat.MMSE.c3 = matrix(thetahat,ncol=1) - K12 %*% Kinv %*% 
      matrix(betahat,ncol=1)
    thetahat.MMSE.Var.c3 = K11 - K12 %*% Kinv %*% K21
  }

  ### MMSE(c4)
  thetahat.MMSE.c4 = NA
  thetahat.MMSE.Var.c4 = NA
  if(!is.na(c4)) {
    K22 <- K[ind2,ind2] + c4*Biases2
    SVD_var <- svd(K22)
    ind.svd <- rep(1,length(SVD_var$d))
    if(length(SVD_var$d)>1)
      for(i in 2:length(SVD_var$d))
        ind.svd[i] <- sum(SVD_var$d[1:(i-1)])/sum(SVD_var$d) < eig.cutoff
    diag_elem <- rep(0,length(SVD_var$d))
    diag_elem[SVD_var$d>0 & ind.svd == 1] <- 
      1/SVD_var$d[SVD_var$d>0 & ind.svd == 1]
    if(length(SVD_var$d>0 & ind.svd == 1) == 1)  
      Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
    if(length(SVD_var$d>0 & ind.svd == 1) > 1)  
      Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
    thetahat.MMSE.c4 = matrix(thetahat,ncol=1) - K12 %*% Kinv %*% 
      matrix(betahat,ncol=1)
    thetahat.MMSE.Var.c4 = K11 - K12 %*% Kinv %*% K21
  }

  ### MMSE(c5)
  thetahat.MMSE.c5 = NA
  thetahat.MMSE.Var.c5 = NA
  if(!is.na(c5)) {
    K22 <- K[ind2,ind2] + c5*Biases2
    SVD_var <- svd(K22)
    ind.svd <- rep(1,length(SVD_var$d))
    if(length(SVD_var$d)>1)
      for(i in 2:length(SVD_var$d))
        ind.svd[i] <- sum(SVD_var$d[1:(i-1)])/sum(SVD_var$d) < eig.cutoff
    diag_elem <- rep(0,length(SVD_var$d))
    diag_elem[SVD_var$d>0 & ind.svd == 1] <- 
      1/SVD_var$d[SVD_var$d>0 & ind.svd == 1]
    if(length(SVD_var$d>0 & ind.svd == 1) == 1)  
      Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
    if(length(SVD_var$d>0 & ind.svd == 1) > 1)  
      Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
    thetahat.MMSE.c5 = matrix(thetahat,ncol=1) - K12 %*% Kinv %*% 
      matrix(betahat,ncol=1)
    thetahat.MMSE.Var.c5 = K11 - K12 %*% Kinv %*% K21
  }

  ### MMSE(c6)
  thetahat.MMSE.c6 = NA
  thetahat.MMSE.Var.c6 = NA
  if(!is.na(c6)) {
    K22 <- K[ind2,ind2] + c6*Biases2
    SVD_var <- svd(K22)
    ind.svd <- rep(1,length(SVD_var$d))
    if(length(SVD_var$d)>1)
      for(i in 2:length(SVD_var$d))
        ind.svd[i] <- sum(SVD_var$d[1:(i-1)])/sum(SVD_var$d) < eig.cutoff
    diag_elem <- rep(0,length(SVD_var$d))
    diag_elem[SVD_var$d>0 & ind.svd == 1] <- 
      1/SVD_var$d[SVD_var$d>0 & ind.svd == 1]
    if(length(SVD_var$d>0 & ind.svd == 1) == 1)  
      Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
    if(length(SVD_var$d>0 & ind.svd == 1) > 1)  
      Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
    thetahat.MMSE.c6 = matrix(thetahat,ncol=1) - K12 %*% Kinv %*% 
      matrix(betahat,ncol=1)
    thetahat.MMSE.Var.c6 = K11 - K12 %*% Kinv %*% K21
  }
  ### MVAR
  K22 <- K[ind2,ind2]
  SVD_var <- svd(K22)
  ind.svd <- rep(1,length(SVD_var$d))
  if(length(SVD_var$d)>1)
    for(i in 2:length(SVD_var$d))
      ind.svd[i] <- sum(SVD_var$d[1:(i-1)])/sum(SVD_var$d) < eig.cutoff
  diag_elem <- rep(0,length(SVD_var$d))
  diag_elem[SVD_var$d>0 & ind.svd == 1] <- 
    1/SVD_var$d[SVD_var$d>0 & ind.svd == 1]
  if(length(SVD_var$d>0 & ind.svd == 1) == 1)  
    Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
  if(length(SVD_var$d>0 & ind.svd == 1) > 1)  
    Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
  thetahat.MVAR = matrix(thetahat,ncol=1) - K12 %*% Kinv %*% 
    matrix(betahat,ncol=1)
  thetahat.MVAR.Var = K11 - K12 %*% Kinv %*% K21
  K22inv_VAR <- Kinv
  EVPROP_MVAR <- SVD_var$d/sum(SVD_var$d)
  
  list(Theta.MMSE = thetahat.MMSE,
       Theta.MMSE.c1 = thetahat.MMSE.c1,
       Theta.MMSE.c2 = thetahat.MMSE.c2, 
       Theta.MMSE.c3 = thetahat.MMSE.c3, 
       Theta.MMSE.c4 = thetahat.MMSE.c4, 
       Theta.MMSE.c5 = thetahat.MMSE.c5, 
       Theta.MMSE.c6 = thetahat.MMSE.c6, 
       Theta.MMSE.Var = thetahat.MMSE.Var,
       Theta.MMSE.Var.c1 = thetahat.MMSE.Var.c1,
       Theta.MMSE.Var.c2 = thetahat.MMSE.Var.c2,
       Theta.MMSE.Var.c3 = thetahat.MMSE.Var.c3,
       Theta.MMSE.Var.c4 = thetahat.MMSE.Var.c4,
       Theta.MMSE.Var.c5 = thetahat.MMSE.Var.c5,
       Theta.MMSE.Var.c6 = thetahat.MMSE.Var.c6,
       Theta.MVAR = thetahat.MVAR, 
       Theta.MVAR.Var = thetahat.MVAR.Var,
       Theta.Hat = thetahat, 
       Theta.Hat.Var = K11,
       EVPROP_MVAR = EVPROP_MVAR,
       EVPROP_MMSE = EVPROP_MMSE
       )
}

