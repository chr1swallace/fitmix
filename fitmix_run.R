#!/usr/bin/env Rscript
setwd("~/Projects/ase")
library(data.table)
x=fread("gtex_ebv_indpendent_snps.csv")
## rm(list = ls())

x <- x[!is.na(slope_se)]
x <- x[sample(1:nrow(x),100000)]
obsData=x$slope
obsVars=x$slope_se^2


library(randomFunctions)
args <- getArgs(defaults=list(id=99,
                              k=2,
                              nits=100,
                              thin=10),
                numeric=c("id","k","nits","thin"))
print(args)

## function arguments
fargs <- list(obsData=obsData,
              obsVars=obsVars,
              savefile=paste0("~/scratch/asemix-id",args$id,"-its",args$nits,"-k",args$k,".RData"),
              unique_id=args$id,
              nIts=args$nits,
              K=args$k, #floor(nrow(obsData)/2),
              thinningFreq=args$thin,
              quiet=TRUE)
              ## quiet=FALSE) 

## run
source("fitmix_cw.R") # functions to do the fit
ff=do.call("fitmix",fargs)

(load(fargs$savefile))

par(mfrow=c(1,2))
v <- do.call("rbind",savedClusterVars)
matplot(log(v),type="l")
m <- do.call("rbind",savedMixtureWeights)
matplot(m,type="l")

tail(m)
tail(v)

medtail <- function(M) {
    n <- nrow(M)
    apply(M[-c(1:n/2),],2,median)
}
cbind(m=medtail(m),v=medtail(v))

