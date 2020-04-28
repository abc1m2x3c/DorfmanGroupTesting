#################################################
##############this is an example#################
#################################################
rm(list=ls(all=TRUE))
library(maxLik)
library(rgl)
library(reshape)

source("src.R")

my.data=interative.input()


#The summary of the MLE of p, Se, and Sp and confidence intervals
#input:
#my.data = the input Dorfman testing data
#output:
#
DT.summary(my.data)


#3D plot of the joint confidence region
#require package rgl and reshape
#my.data = the input Dorfman testing data
#type    = the type of confidence region. Currently support 'Wald' (Wald confidence region)
#          and 'LR' (Likelihood confidence region).
DT.plot(my.data,type='Wald')

#The summary of charateristics of 
DT.chr(my.data)




