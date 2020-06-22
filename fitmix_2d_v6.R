## in v3 allow correlation in all sets, but higher correlation in group 4
## this correlation results from the same individuals being sampled for the two tissues

library("progress")
library("bivnorm")
f <- function(x,y) sum(y == x)
matlogdnorm <- function(dmat,mvec,vvec,vmat) {
    rmat <- function(x) matrix(x,nrow=nrow(dmat),ncol=ncol(dmat),byrow=TRUE)
    cmat <- function(x) matrix(x,nrow=nrow(dmat),ncol=ncol(dmat),byrow=FALSE)
    cmat(dnorm(dmat,
               mean = rmat(mvec),
               sd = sqrt( rmat(vvec) + vmat),
               log = T))
}
whichfirst <- function(x) {
    which(x)[1]
}

## Informative priors, from univariate analyses
##nonSymmetric Dirichlet parameter:
a <- 1/4 #c(0.97^2,rep(1-0.97^2,3)/3)
## sample clustervars 
## (97%, sd=0.03) and a broader distribution (3%, sd=0.35),
## log(0.03^2)=-7
## log(0.35^2)=-2.1
## m0 <- log(c(0.001,0.13)) # from univariate fits
## s0 <- c(0.12,0.09)
m0 <- log(c(0.15,0.15))
s0 <- c(2,2)

sampleParams <- function() {
  v <- c((exp(rnorm(2, mean = m0, sd = s0))),
        (runif(4,-0.5,0.5))) # cor between +/- 0.5
  v
}
calcClusterParams <- function(params) {
  matrix(c(params[1],params[1],params[3],
           params[2],params[1],params[4],
           params[1],params[2],params[5],
           params[2],params[2],params[6]),
         4,3,byrow=TRUE)
}

var.meta <- function(vlist,nlist) {
  vlist*(nlist-1)/(sum(nlist)-length(nlist))
}
calcCovar <- function(v1,v2,rho)
  sqrt(v1)*sqrt(v2)*rho

## select k which depend on param j
updates <- list(1:3, 2:4, 1, 2, 3, 4)

estParams <- function(obsData,clusterAllocations) {
  tmp <- sampleParams() # if no data will sample from prior
  v1 <- tmp[1]; v2 <- tmp[2]; r1 <- tmp[3]; r2 <- tmp[4]; r3 <- tmp[5]; r4 <- tmp[6]
  wh <- lapply(1:4, function(k) which(clusterAllocations==k))
  if(length(unlist(wh[1:3]))>1)
    v1 <-  var(c(obsData[unlist(wh[c(1,3)]),1],
                 obsData[unlist(wh[c(1,2)]),2]))
  if(length(unlist(wh[2:4]))>1)
    v2 <-  var(c(obsData[unlist(wh[c(2,4)]),1],
                 obsData[unlist(wh[c(3,4)]),2]))
  if(length((wh[[1]]))>2)
    r1 <- cor(obsData[wh[[1]],])[1,2]
  if(length((wh[[2]]))>2)
    r2 <- cor(obsData[wh[[2]],])[1,2]
  if(length((wh[[3]]))>2)
    r3 <- cor(obsData[wh[[3]],])[1,2]
  if(length(wh[[4]])>2)
    r4 <- cor(obsData[wh[[4]],])[1,2]
  c(v1,v2,r1,r2,r3,r4) # sigma0^1, sigma1^2, rho
}

## convenience function
## create a symmetric 2x2 matrix from the diagonal and off diagnonal entries
vec2mat <- function(A)
  matrix(c(A[1],A[3],A[3],A[2]),2,2)
pvec2sigma <- function(pvec) # pvec = var1, var2, cor
  vec2mat(c(pvec[1:2],pvec[3]*sqrt(prod(pvec[1:2])))) # sigma0^1, sigma1^2, sigma12^2

## convenience function
## create a vector containing the diagonal and off diagnonal entries from a symmetric 2x2 matrix 
mat2vec <- function(M)
  c(M[1,1],M[2,2],M[1,2])




## fits a K-gaussian 2-d mixture, means are fixed at 0, estimates
## proportions and cluster variances
fitmix <- function(obsData, obsVars, savefile,
                    unique_id=42, nIts=1000,K=4,thinningFreq=100,
                    quiet=FALSE,
                    init.m0=-8,init.s0=1) {
                                        #Max number of clusters: 
    set.seed(unique_id)
    d    <- 2 # number of dimensions
    n    <- nrow(obsData) # number of observations

    ## ##Parameters of normal prior on log(rho): centers on 0.0225, iqr=2x10-5 -- 8x10-5
    ## m0  = log(0.0225)  %>% rep(.,d+1) #rep(init.m0,d)
    ## s0  = 1 # broad variance - range 0.002 - 0.2

  clusterMeans <- matrix(0,K,d) 
  clusterParams <- numeric(6) 
  sumSqResiduals <- vector(mode="numeric",length=K)
  mixtureWeights     <- vector(mode = "numeric", length = K)
  clusterAllocations <- vector(mode = "numeric", length = n)
  nClusters          <- vector(mode = "numeric", length = nIts)
  
  thinner                <- thinningFreq
  nSaved                 <- round(nIts / thinner)
  counter                <- 1
  ## savedClusterAllocations <- matrix(nrow = nSaved, ncol = n)
  savedParams <-
    savedMixtureWeights <-
      vector(mode = "list", length = nSaved)
  
    ## kmeansInit             <- kmeans(obsData, centers = K)

  ## clusterAllocations <- sample(1:K,n,replace=TRUE)
  ##   nAllocsPerCluster      <- sapply(1:K, f, y = clusterAllocations)
  mixtureWeights         <-
     MCMCpack::rdirichlet(1, c(0.97^2, rep(1-0.97^2,3)/3)*n + a)
  ## colMeans(tmp)
  ## apply(tmp,2,quantile,c(0.01,0.1,0.25,0.5,0.75,0.9,0.99))
  
    varianceAcceptances <- 0
    ## meanAcceptances <- meanTrials <- matrix(0, nrow = K, ncol = d) 

  params <- sampleParams()
  
  if(quiet) {
    pb <- progress_bar$new(total = nIts/thinner)
    pb$tick(0)
  }

  for(its in 1:nIts) {

### Update the component allocations 
    clusterParams <- calcClusterParams(params)
    logLikeliVec <- matrix(0,n,K)
    for(k in 1:K) {
      logLikeliVec[,k] <- dbivnorm(obsData,
                                   mean=clusterMeans[k,],
                                   v1=clusterParams[k,1] + obsVars[,1],
                                   v2=clusterParams[k,2] + obsVars[,2],
                                   v12=calcCovar(clusterParams[k,1] + obsVars[,1],
                                                 clusterParams[k,2] + obsVars[,2],
                                                 clusterParams[k,3]),
                                   log=TRUE)
    }
    
    logAllocProbs          <- matrix(log(mixtureWeights),n,k,byrow=TRUE) + logLikeliVec
    allocProbs             <- exp(logAllocProbs - apply(logAllocProbs,1,coloc:::logsum))
    for(k in 2:K)
      allocProbs[,k] <- allocProbs[,k] + allocProbs[,k-1]
    clusterAllocations <- apply(runif(n) <= allocProbs, 1, whichfirst)

### Update the mixture weights

    nAllocsPerCluster <- sapply(1:K, f, y = clusterAllocations)
    mixtureWeights    <- MCMCpack::rdirichlet(1, nAllocsPerCluster + a)
    
### Update the cluster variances conditional on current allocations
    ## params[3:6] <- estParams(obsData, clusterAllocations)[3:6]
    params <- estParams(obsData, clusterAllocations)
    ## params[1:2] <- sampleParams()[1:2]
    proposedParams <- 0.95 * params + 0.05 * sampleParams() # move slowly to limit label switching
    ## update vars, cors
    for(j in 1:6) {
      if(j %in% c(3:6)) {
        logPriorCurrent <- logPriorProposed <- 0 # difference between two uniform is 0
      } else { # only changing one, so only care about prior difference at that one
        logPriorCurrent <- dnorm(params[j],mean=m0[j],sd=s0[j],log=TRUE) 
        logPriorProposed <- dnorm(proposedParams[j],mean=m0[j],sd=s0[j],log=TRUE)
      }
      proposedLogLikeli <- currentLogLikeli <- 0
      clusterParams <- calcClusterParams(params)
      proposedClusterParams <- calcClusterParams(ifelse(j==1:6,proposedParams,params))
      for(k in updates[[j]]) {
        allocInds <- which(clusterAllocations == k)
        if(length(allocInds)<=2) {
          ## here
          if(k>2) { # no information in other groups, sample from prior
            params[j] <- sampleParams()[j] # take correlation and mult by current sigma1^2
          }
          next
        }
        currentLogLikeli  <- currentLogLikeli + 
          sum(dbivnorm(obsData[allocInds,,drop=FALSE],
                       mean=clusterMeans[k,],
                       v1=clusterParams[k,1] + obsVars[allocInds,1],
                       v2=clusterParams[k,2] + obsVars[allocInds,2],
                       v12=calcCovar(clusterParams[k,1] + obsVars[allocInds,1],
                                     clusterParams[k,2] + obsVars[allocInds,2],
                                     clusterParams[k,3]),
                       log=TRUE))
        proposedLogLikeli  <- proposedLogLikeli +
          sum(dbivnorm(obsData[allocInds,,drop=FALSE],
                       mean=clusterMeans[k,],
                       v1=clusterParams[k,1] + obsVars[allocInds,1],
                       v2=clusterParams[k,2] + obsVars[allocInds,2],
                       v12=calcCovar(clusterParams[k,1] + obsVars[allocInds,1],
                                     clusterParams[k,2] + obsVars[allocInds,2],
                                     proposedClusterParams[k,3]),
                       log=TRUE))
      }
      logMHratio <- proposedLogLikeli + logPriorProposed -
        currentLogLikeli - logPriorCurrent 
      ## message(logMHratio)
      
      Accept <- log(runif(1)) < logMHratio
      if(Accept) {
        params[j] <- proposedParams[j]
        varianceAcceptances <- varianceAcceptances + 1 
      }
      ## if(j==2 & params[2] < params[1])
      ##   params[1:2] <- sort(params[1:2])
    }

    if(its%%thinner == 0){
      if(quiet)
        pb$tick()
      if(!quiet) {
        message("iteration ", its)
        message("params")
        print(params)
        message("mixweights")
        print(mixtureWeights)
      }
      ##   message(its, "\t", clusterVars[1],"\t",clusterVars[2],"\t",nAllocsPerCluster[1],"\t",nAllocsPerCluster[2])
      
      ## savedClusterMeans[[counter]]   <- clusterMeans
      savedParams[[counter]]    <- params
      savedMixtureWeights[[counter]] <- mixtureWeights
      ## savedClusterAllocations[counter,] <- clusterAllocations
      counter <- counter + 1
      if(its %% 50000 == 0) 
        print("Saving...")
        save(savedMixtureWeights,
             savedParams, file = savefile)
    }
  }
  if(its %% 50000 !=0)
    save(savedMixtureWeights,
         savedParams, file = savefile)
  return(savefile)
}

