#!/usr/bin/env Rscript
setwd("~/Projects/ase")
library(magrittr)
library(data.table)
## x=fread("gtex_ebv_indpendent_snps.csv")
## rm(list = ls())

## x <- x[!is.na(slope_se)]
## x <- x[sample(1:nrow(x),100000)]
## obsData=x$slope
## obsVars=x$slope_se^2

library(randomFunctions)
args <- getArgs(defaults=list(t1="blood",
                              id=99,
                              k=2,
                              nits=100,
                              thin=10),
                numeric=c("id","k","nits","thin"))
print(args)

files <- list(skin="Skin_Not_Sun_Exposed_Suprapubic_independent_common_snps.txt",
              ebv="gtex_ebv_independent_common_snps.txt",
              blood="Whole_Blood_independent_common_snps.txt")
if(!args$t1 %in% names(files))
  stop("type not known: ",args$t1)
filter <- function(x)
  x[!is.na(slope_se)]

x=fread(files[[args$t1]])  %>% filter()
x <- x[sample(1:nrow(x),100000)]
obsData=x$slope
obsVars=x$slope_se^2

## function arguments
fargs <- list(obsData=obsData,
              obsVars=obsVars,
              savefile=paste0("~/scratch/asemix-",args$t1,"-id",args$id,"-its",args$nits,"-k",args$k,".RData"),
              unique_id=args$id,
              nIts=args$nits,
              K=args$k, #floor(nrow(obsData)/2),
              thinningFreq=args$thin,
              quiet=TRUE)
              ## quiet=FALSE) 

if(interactive())
  attach(fargs)

print(fargs)
## run
source("fitmix_cw.R") # functions to do the fit
ff=do.call("fitmix",fargs)

if(!interactive())
  q("no")

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
    c(apply(M[-c(1:n/2),],2,median),
      apply(M[-c(1:n/2),],2,mean),
      apply(M[-c(1:n/2),],2,sd))
}
data.frame(stat=c("median","median","mean","mean","sd","sd"),
           m=medtail(m),
           v=medtail(v),
           logv=medtail(log(v)))

