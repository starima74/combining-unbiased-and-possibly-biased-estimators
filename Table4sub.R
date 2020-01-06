
### save the file in a working folder defined below in XXX
path <- "XXX"
setwd(path)
source(paste(path,"/Table4_global.R",sep="")) 
source(paste(path,"/Table4_local.R",sep="")) 


source("MMSE_C6.R") 
library(quantreg)  
library(BART) 
library(mgcv)
library(car)
library(quantreg)
library(caret)
library(splines)
library(SuperLearner)
set.seed(100) 

local <- Table4_local()
global <- Table4_global() 

Table4 <- rbind(local,global)