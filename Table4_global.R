#########################################################################
# This function implements global methods given in table 4
#########################################################################

Table4_global <- function(){
taus <- c(0.25,0.5,0.75)
a.train<- read.csv("a.train.csv")
a.test <- read.csv("a.test.csv")

f2 <- formula(paste("D6 ~ A1C0 + age", sep=""))

f3 <- formula(paste("D6 ~ A1C0 + age + I(age^2)", sep=""))

  

dat <- a.train[,c("D6","A1C0","age")] 
x.train <- a.train[,c("A1C0","age")]
x.test <- as.matrix(a.test[,c("A1C0","age")]) 
y.train <-  dat[,"D6"]
 
#QR_nonparametric predictors 
QR_nonlin <- predict(rq(D6~bs(A1C0,15)+bs(age,15), tau=taus, data=a.train),newdata=a.test)

#LM with linear predictors

y0hat <- predict(lm(f2,data=a.train), newdata = a.test, se=T)
LM <- matrix(NA,nrow=nrow(x.test),ncol=length(taus))

for(j in 1:length(taus)){
  LM[,j] <- y0hat$fit + qnorm(taus[j]) * y0hat$se.fit
}

#Random forest
 
wi<- 0.75


#LAD: L1 

 
L1_25_QR_nonlin <- 2*mean(wi*I(a.test$D6<QR_nonlin[,1])*abs(a.test$D6-QR_nonlin[,1]) + (1-wi)*I(a.test$D6 >= QR_nonlin[,1])*abs(a.test$D6-QR_nonlin[,1]))
L1_50_QR_nonlin <- 2*mean(0.5*I(a.test$D6<QR_nonlin[,2])*abs(a.test$D6-QR_nonlin[,2]) + 0.5*I(a.test$D6 >= QR_nonlin[,2])*abs(a.test$D6-QR_nonlin[,2]))
L1_75_QR_nonlin <- 2*mean((1-wi)*I(a.test$D6<QR_nonlin[,3])*abs(a.test$D6-QR_nonlin[,3]) +  wi*I(a.test$D6 >= QR_nonlin[,3])*abs(a.test$D6-QR_nonlin[,3]))

L1_25_LM <- 2*mean(wi*I(a.test$D6<LM[,1])*abs(a.test$D6-LM[,1]) + (1-wi)*I(a.test$D6 >= LM[,1])*abs(a.test$D6-LM[,1]))
L1_50_LM <- 2*mean(0.5*I(a.test$D6<LM[,2])*abs(a.test$D6-LM[,2]) + 0.5*I(a.test$D6 >= LM[,2])*abs(a.test$D6-LM[,2]))
L1_75_LM <- 2*mean((1-wi)*I(a.test$D6<LM[,3])*abs(a.test$D6-LM[,3]) +  wi*I(a.test$D6 >= LM[,3])*abs(a.test$D6-LM[,3]))
 
#MSE:L2

 L2_25_QR_nonlin <- 2*mean(wi*I(a.test$D6<QR_nonlin[,1])*(a.test$D6-QR_nonlin[,1])^2 + (1-wi)*I(a.test$D6 >= QR_nonlin[,1])*(a.test$D6-QR_nonlin[,1])^2)
L2_50_QR_nonlin <- 2*mean(0.5*I(a.test$D6<QR_nonlin[,2])*(a.test$D6-QR_nonlin[,2])^2 + 0.5*I(a.test$D6 >= QR_nonlin[,2])*(a.test$D6-QR_nonlin[,2])^2)
L2_75_QR_nonlin <- 2*mean((1-wi)*I(a.test$D6<QR_nonlin[,3])*(a.test$D6-QR_nonlin[,3])^2 +  wi*I(a.test$D6 >= QR_nonlin[,3])*(a.test$D6-QR_nonlin[,3])^2)

L2_25_LM <- 2*mean(wi*I(a.test$D6<LM[,1])*(a.test$D6-LM[,1])^2 + (1-wi)*I(a.test$D6 >= LM[,1])*(a.test$D6-LM[,1])^2)
L2_50_LM <- 2*mean(0.5*I(a.test$D6<LM[,2])*(a.test$D6-LM[,2])^2 + 0.5*I(a.test$D6 >= LM[,2])*(a.test$D6-LM[,2])^2)
L2_75_LM <- 2*mean((1-wi)*I(a.test$D6<LM[,3])*(a.test$D6-LM[,3])^2 +  wi*I(a.test$D6 >= LM[,3])*(a.test$D6-LM[,3])^2)
 
#BIAS

bias_25_QR_nonlin <- 2*mean(wi*I(a.test$D6<QR_nonlin[,1])* (a.test$D6-QR_nonlin[,1]) + (1-wi)*I(a.test$D6 >= QR_nonlin[,1])* (a.test$D6-QR_nonlin[,1]))
bias_50_QR_nonlin <- 2*mean(0.5*I(a.test$D6<QR_nonlin[,2])* (a.test$D6-QR_nonlin[,2]) + 0.5*I(a.test$D6 >= QR_nonlin[,2])* (a.test$D6-QR_nonlin[,2]))
bias_75_QR_nonlin <- 2*mean((1-wi)*I(a.test$D6<QR_nonlin[,3])* (a.test$D6-QR_nonlin[,3]) +  wi*I(a.test$D6 >= QR_nonlin[,3])* (a.test$D6-QR_nonlin[,3]))

bias_25_LM <- 2*mean(wi*I(a.test$D6<LM[,1])* (a.test$D6-LM[,1]) + (1-wi)*I(a.test$D6 >= LM[,1])* (a.test$D6-LM[,1]))
bias_50_LM <- 2*mean(0.5*I(a.test$D6<LM[,2])* (a.test$D6-LM[,2]) + 0.5*I(a.test$D6 >= LM[,2])* (a.test$D6-LM[,2]))
bias_75_LM <- 2*mean((1-wi)*I(a.test$D6<LM[,3])* (a.test$D6-LM[,3]) +  wi*I(a.test$D6 >= LM[,3])* (a.test$D6-LM[,3]))

 
############
L <- c("L1","L2","BIAS", "L1","L2","BIAS", "L1","L2","BIAS")
TAU <- c(rep(0.25,3),rep(0.5,3),rep(0.75,3))

 

QR_nonlinS <- c( L1_25_QR_nonlin,L2_25_QR_nonlin, bias_25_QR_nonlin,
                 L1_50_QR_nonlin,L2_50_QR_nonlin, bias_50_QR_nonlin,
                 L1_75_QR_nonlin,L2_75_QR_nonlin, bias_75_QR_nonlin) 

LMS <- c( L1_25_LM,L2_25_LM, bias_25_LM,
          L1_50_LM,L2_50_LM, bias_50_LM,
          L1_75_LM,L2_75_LM, bias_75_LM) 
 

ALLTAUS_GLOBAL <- data.frame(L = L,
                        TAU = TAU,   
                        QR_SP = QR_nonlinS,
                        LM = LMS) 

result <- cbind(ALLTAUS_GLOBAL[TAU!=0.5,1:2],round(ALLTAUS_GLOBAL[TAU!=0.5,-(1:2)],4))
 
  
write.csv(data.frame(t(result)),"ALLTAUS_GLOBAL.csv")
return(data.frame(t(result)))
}


