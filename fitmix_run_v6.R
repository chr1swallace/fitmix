#!/usr/bin/env Rscript
setwd("~/Projects/ase")
library(data.table)
library(magrittr)
library(randomFunctions)
args <- getArgs(defaults=list(t1="blood",t2="ebv",
                              id=99,
                              k=4,
                              nits=100,
                              thin=10),
                numeric=c("id","k","nits","thin"))
print(args)

files <- list(skin="Skin_Not_Sun_Exposed_Suprapubic_independent_common_snps.txt",
              ebv="gtex_ebv_independent_common_snps.txt",
              blood="Whole_Blood_independent_common_snps.txt")
if(!args$t1 %in% names(files))
  stop("type not known: ",args$t1)
if(!args$t2 %in% names(files))
  stop("type not known: ",args$t2)

filter <- function(x)
  x[!is.na(slope_se)]

x=fread(files[[args$t1]])  %>% filter()
y=fread(files[[args$t2]])  %>% filter()
## rm(list = ls())

## check
m <- merge(x,y,by=c("gene_id","variant_id"),suffixes=c(".x",".y"))
## set.seed(42);
m <- m[sample(1:nrow(m),100000)]
obsData=m[,.(slope.x,slope.y)]  %>% as.matrix()
obsVars=m[,.(slope_se.x^2,slope_se.y^2)]  %>% as.matrix()

## function arguments
fargs <- list(obsData=obsData,
              obsVars=obsVars,
              savefile=paste0("~/scratch/asemix2dv6-",args$t1,"-",args$t2,"-id",args$id,"-its",args$nits,"-k",args$k,".RData"),
              unique_id=args$id,
              nIts=args$nits,
              K=args$k, #floor(nrow(obsData)/2),
              thinningFreq=args$thin,
              quiet=TRUE)
              ## quiet=FALSE) 

if(interactive()) # for testing
  attach(fargs)

              ## fargs$nIts <- 1000
## savefile=paste0("~/scratch/asemix2dv3-",args$t1,"-",args$t2,"-id",args$id,"-its",fargs$nIts,"-k",args$k,".RData"); unique_id=args$id; nIts=args$nits; K=args$k; thinningFreq=args$thin; quiet=TRUE;

## run
## source("fitmix_2d.R") # functions to do the fit with independent inv wishart priors
## source("fitmix_2d_alt.R") # functions to do the fit with 3 parameter variance model
## source("fitmix_2d_v3.R") # functions to do the fit with 3 parameter variance model, component wise updates
source("fitmix_2d_v6.R") # functions to do the fit with 4 parameter variance model, component wise updates
ff=do.call("fitmix",fargs)

if(!interactive())
  q("no")

(load(fargs$savefile))

## plot mixing weights
par(mfrow=c(1,3))
mixweights <- do.call("rbind",savedMixtureWeights) 
matplot((mixweights),type="l")
matplot((mixweights),type="l",ylim=c(0,0.1))
legend("topleft",col=1:4,lty=rep(1,4),legend=c("neither","1","2","12"))

params <- do.call("rbind",savedParams) 
matplot(params,type="l")
legend("topleft",col=1:4,lty=rep(1,4),legend=c("v1","v2","r1","r2"))

colMeans(mixweights)
colMeans(params)

## library(data.table)
## library(ggplot2)
## mixweights <- do.call("rbind",savedMixtureWeights)  %>% as.data.frame()  %>% as.data.table()
## mixweights$it <- 1:nrow(mixweights)
## m <- melt(mixweights)
## ggplot(m, aes(x=m,y=value,col=variable)) + geom_path()

library(ellipse)
library(scales)
plotvar <- function(v,add=FALSE,...) {
  m <- vec2mat(v)
  ## if(add)
    lines(ellipse(m),...)
  ## else
  ##   plot(ellipse(m),type="l",...)
}

savedClusterVars <- lapply(savedParams, calcClusterVars)
par(mfrow=c(3,2))
mx <- lapply(savedClusterVars, apply, 2, max)  %>%
  do.call("rbind",.)  %>%
  apply(., 2, max)
el <- ellipse(mx)
nj <- length(savedClusterVars)
for(i in 1:4) {
  plot(0,0,xlim=c(min(el[,1]),max(el[,1])),ylim=c(min(el[,2]),max(el[,2])))
  for(j in 1:length(savedClusterVars))
      plotvar(savedClusterVars[[j]][i,],col=alpha(i,j/nj),add=(j>1))
}
matplot(log(mixweights),type="l")
legend("topleft",col=1:4,lty=rep(1,4),legend=c("neither","1","2","12"))

params <- do.call("rbind",savedParams)
params[,3] <- params[,3]*params[,2]
matplot(params,type="l")
legend("topleft",col=1:3,lty=rep(1,3),legend=c("v0","v1","rho"))

plotmix <- function() {
}



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

