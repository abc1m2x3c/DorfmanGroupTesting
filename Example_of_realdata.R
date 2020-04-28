#################################################
##############this is an example#################
#################################################
rm(list=ls(all=TRUE))
library(maxLik)
library(rgl)
library(reshape)

source("https://raw.githubusercontent.com/abc1m2x3c/DorfmanGroupTesting/master/src.R")


#input the data
my.data=interative.input()


#The summary of the MLE of p, Se, and Sp and confidence intervals
#require package maxLik
#input:
#my.data = the input Dorfman testing data
#output:
#a matrix containing p, Se, and Sp's
#MLE, standard error, 95% Wald confidence interval, and 95% Profile confidence interval
DT.fit(my.data)


#3D plot of the joint confidence region
#require package rgl and reshape
#input
#my.data = the input Dorfman testing data
#type    = the type of confidence region. Currently support 'Wald' (Wald confidence region)
#          and 'LR' (Likelihood confidence region).
#output:
#a 3D plot that can be rotated interactively
#the volume of the joint area
DT.plot(my.data,type='Wald')

#The Wald confidence intervals of group testing charateristics
#input:
#my.data = the input Dorfman testing data
#output
#for each pool size, the Wald confidence interval of 
#E(T): expected number of tests per individual
#E(C): expected number of correct classifications per individual specimen
#PPV : The (pooling) positive predictive value
#NPV : The (pooling) negative predictive value
DT.chr(my.data)




