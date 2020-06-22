#!/usr/bin/env Rscript
setwd("~/Projects/ase")
library(data.table)
library(magrittr)
library(randomFunctions)
library(seaborn)
library(ggplot2); library(cowplot); theme_set(theme_cowplot())
list.files("~/scratch",pattern="asemix-")
## system("ls -lhtr ~/scratch/asemix2dv4*")
source("~/Projects/ase/fitmix_cw.R")
source("~/Projects/ase/fitmix_functions.R")

results <- loader("asemix-")

sumtail(results$weights) # 97% / 3 %
sumtail(results$params)

plotter(results)


library(ellipse)
library(scales)
plotvar <- function(v,add=FALSE,...) {
  m <- vec2mat(v)
  ## if(add)
    lines(ellipse(m),...)
  ## else
  ##   plot(ellipse(m),type="l",...)
}

par(mfrow=c(2,2))
#4,length(results$vars)))
mx <- lapply(results$vars, tail, 1) %>% # extract last element
  lapply(., "[[",1) %>% # make it not a list
  lapply(., "[", 4,1:3,drop=FALSE) %>% # take last row
  do.call("rbind",.)  %>%
  apply(., 2, max) ## rbind
el <- ellipse(mx)
nj <- length(vars)
for(i in 1:4) { # iterate group
  plot(0,0,xlim=c(min(el[,1]),max(el[,1])),ylim=c(min(el[,2]),max(el[,2])))
for(j in seq_along(results$vars)) { # iterate chain
  for(k in 1:length(results$vars[[j]]))
    plotvar(results$vars[[j]][[k]][i,],col=alpha(j,k/length(results$vars[[j]])))
}
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

