
### save the file in a working folder defined below in XXX
path <- "XXX"
setwd(path)
source(paste(path,"/Table3_global.R",sep=""))
source(paste(path,"/Table3_glaslocal.R",sep=""))
source(paste(path,"/Table3_local.R",sep=""))


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

global <- Table3_global()
glaslocal <- Table3_glaslocal()
local <- Table3_local()

Table3 <- rbind(global,glaslocal,local)
Table3