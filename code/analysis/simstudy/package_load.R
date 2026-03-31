rm(list = ls())

library(fields)
library(SpatialTools)
library(Rcpp)

#### Load simdata
load(file='./simdata.RData')
source('../../R/mcmc_cont_lambda.R', chdir = TRUE)
source('../../R/auxfunctions.R')
source('./max-stab/Bayes_GEV.R')
source('./max-stab/MCMC4MaxStable.R', chdir = TRUE)

options(warn=2)