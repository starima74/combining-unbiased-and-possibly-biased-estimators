library(tools)
library(MASS)
library(quantreg)
library(plyr)
library(utils)
library(car)
library(MASS)
library(splines) 
library(quantregForest) 
### save the file in a working folder defined below in XXX
path <- "XXX"
setwd(path)
source(paste(path,"/ensemble.R",sep=""))
source(paste(path,"/bootbeta.R",sep=""))
source(paste(path,"/predvalK.R",sep=""))
source(paste(path,"/post_proc.R",sep=""))

################################################################################
# Run the simulation  to reproduce table 1 and 2 of MSEs for our estimators
################################################################################

seed<- 100
set.seed(seed)
nsims <- 500
N <- 400
B <- 200
f0 <- formula(Y ~ X.1 + X.2 ) 
Betas <- c(-1,1)
taus <- c(0.1,0.5,0.9)

newdata <- data.frame(X.1 = c(1,1,0,0),
                      X.2 = c(1,0,1,0))

 
  truetheta_corr <- matrix(NA, nrow=nrow(newdata),ncol = length(taus))
  truetheta_misp <- matrix(NA, nrow=nrow(newdata),ncol = length(taus))
  for(j in 1:length(taus))
  {
    # For correct model
    truetheta_corr[,j] <- Betas[1]*newdata$X.1 + Betas[2]*newdata$X.2 + qnorm(taus[j])
    truetheta_misp[,j] <- Betas[1]*newdata$X.1 + Betas[2]*newdata$X.2 + exp(qnorm(taus[j]))-1
  } 

################################################################
# Get output for misspecified and correct models, for binary and 
# Continuous covariates
################################################################
out_bin_corr <- post_proc(X.1 = rbinom(N,1,prob=0.5),
                          X.2 =  rbinom(N,1,prob=0.5),
                          eps = rnorm(N))
saveRDS(out_bin_corr,"out_bin_corr1.rds")

out_bin_misp <- post_proc(X.1 = rbinom(N,1,prob=0.5),
                     X.2 =  rbinom(N,1,prob=0.5),
                     eps =  exp(rnorm(N))-1)
saveRDS(out_bin_misp,"out_bin_misp1.rds")

out_cont_corr <- post_proc(X.1 = rnorm(N),
                      X.2 = rnorm(N),
                      eps = rnorm(N))
saveRDS(out_cont_corr,"out_cont_corr1.rds")

out_cont_misp <- post_proc(X.1 = rnorm(N),
                      X.2 = rnorm(N),
                      eps =  exp(rnorm(N))-1)
saveRDS(out_cont_misp,"out_cont_misp1.rds")

##################################################################################
#Reproduce tables 1 and 2
##################################################################################

########################################################
#Correct model, Binary covariates
###############################################################

thetaests <- matrix(colMeans(out_bin_corr),
                    nrow=12,ncol=9)

truemses <-t((t(out_bin_corr)-rep(truetheta_corr,9))^2)
truemse <- matrix(apply(truemses,2,mean),nrow=12,ncol=9)


TABLE1_cb <- data.frame(X1  = newdata$X.1,
                     X2  = newdata$X.2,
                     tau = c(rep(taus[1],4),rep(taus[2],4),rep(taus[3],4)),
                     QR  = round(truemse[,1],4),
                     MLE = round(truemse[,2],4),
                     QRA = round(truemse[,8],4),
                     QRF = round(truemse[,9],4)
)


########################################################
#Mispecified model, Binary covariates
###############################################################

thetaests <- matrix(colMeans(out_bin_misp),
                    nrow=12,ncol=9)

truemses <-t((t(out_bin_misp)-rep(truetheta_misp,9))^2)
truemse <- matrix(apply(truemses,2,mean),nrow=12,ncol=9)


TABLE1_mb <- data.frame(X1  = newdata$X.1,
                        X2  = newdata$X.2,
                        tau = c(rep(taus[1],4),rep(taus[2],4),rep(taus[3],4)),
                        QR  = round(truemse[,1],4),
                        MLE = round(truemse[,2],4),
                        QRA = round(truemse[,8],4),
                        QRF = round(truemse[,9],4)
)


########################################################
#Correct model, continuous covariates
###############################################################



thetaests <- matrix(colMeans(out_cont_corr),
                    nrow=12,ncol=9)

truemses <-t((t(out_cont_corr)-rep(truetheta_misp,9))^2)
truemse <- matrix(apply(truemses,2,mean),nrow=12,ncol=9)


TABLE1_cc <- data.frame(X1  = newdata$X.1,
                        X2  = newdata$X.2,
                        tau = c(rep(taus[1],4),rep(taus[2],4),rep(taus[3],4)),
                        QR  = round(truemse[,1],4),
                        MLE = round(truemse[,2],4),
                        QRA = round(truemse[,8],4),
                        QRF = round(truemse[,9],4)
)



########################################################
#Correct model, continuous covariates
###############################################################



thetaests <- matrix(colMeans(out_cont_misp),
                    nrow=12,ncol=9)

truemses <-t((t(out_cont_misp)-rep(truetheta_misp,9))^2)
truemse <- matrix(apply(truemses,2,mean),nrow=12,ncol=9)


TABLE1_mc <- data.frame(X1  = newdata$X.1,
                        X2  = newdata$X.2,
                        tau = c(rep(taus[1],4),rep(taus[2],4),rep(taus[3],4)),
                        QR  = round(truemse[,1],4),
                        MLE = round(truemse[,2],4),
                        QRA = round(truemse[,8],4),
                        QRF = round(truemse[,9],4)
)


TABLE <- rbind(cbind(TABLE1_cb,TABLE1_mb),
               cbind(TABLE1_cc,TABLE1_mc))
TABLE1 <- TABLE[!(TABLE$X1==0 & TABLE$X2==1),]

saveRDS(TABLE1,"TABLE1.rds")

