library("progress")
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
    
## fits a K-gaussian 1-d mixture, means are fixed at 0, estimates
## proportions and cluster variances
fitmix <- function(obsData, obsVars, savefile,
                    unique_id=42, nIts=1000,K=2,thinningFreq=100,
                    quiet=FALSE,
                    init.m0=-8,init.s0=1) {
                                        #Max number of clusters: 
    set.seed(unique_id)
    d    <- 1
    n    <- length(obsData)

                                        #Symmetric Dirichlet parameter:
    a    <- 1/K

                                        #Parameters of normal prior on mu:
    ## mu0     = mean(obsData)
    ## sigma0  = var(obsData)


    ##Parameters of normal prior on log(rho): centers on 0.0225, iqr=2x10-5 -- 8x10-5
    m0  = log(0.0225) #rep(init.m0,d)
    s0  = 1 # broad variance - range 0.002 - 0.2

    clusterMeans       <- rep(0,K) #numeric(K)
    clusterVars <- c(0.1*exp(m0),rep(1,K-1))
    sumSqResiduals <- vector(mode="numeric",length=K)
    mixtureWeights     <- vector(mode = "numeric", length = K)
    clusterAllocations <- vector(mode = "numeric", length = n)
    nClusters          <- vector(mode = "numeric", length = nIts)

    thinner                <- thinningFreq
    nSaved                 <- round(nIts / thinner)
    counter                <- 1
    savedClusterAllocations <- matrix(nrow = nSaved, ncol = n)
    ## savedClusterMeans <-
        savedClusterVars <-
            savedMixtureWeights <-
                vector(mode = "list", length = nSaved)

    ## kmeansInit             <- kmeans(obsData, centers = K)
    clusterAllocations <- sample(1:K,n,replace=TRUE)
    nAllocsPerCluster      <- sapply(1:K, f, y = clusterAllocations)
    mixtureWeights         <- #c(0.9,0.1)
    MCMCpack::rdirichlet(1, nAllocsPerCluster + a)

    varianceAcceptances <- varianceTrials  <- vector(mode="numeric",length=d)
    meanAcceptances <- meanTrials <- matrix(0, nrow = K, ncol = d) 

    if(quiet) {
        pb <- progress_bar$new(total = nIts/thinner)
        pb$tick(0)
    }

    for(its in 1:nIts)
    {
###Update the component allocations 
        logLikeliVec <- matrix(0,n,K)
        for(k in 1:K) {
            logLikeliVec[,k]<- dnorm(obsData,
                                    mean = clusterMeans[k],
                                    sd = sqrt(clusterVars[k] + obsVars),
                                    log = T)
        }
        logAllocProbs          <- matrix(log(mixtureWeights),n,k,byrow=TRUE) + logLikeliVec
        allocProbs             <- exp(logAllocProbs - apply(logAllocProbs,1,coloc:::logsum))
        for(k in 2:K)
            allocProbs[,k] <- allocProbs[,k] + allocProbs[,k-1]
        clusterAllocations <- apply(runif(n) <= allocProbs, 1, whichfirst)

###Update the mixture weights
        nAllocsPerCluster <- sapply(1:K, f, y = clusterAllocations)
        mixtureWeights    <- MCMCpack::rdirichlet(1, nAllocsPerCluster + a)
        
###Update the cluster variances conditional on current allocations
        ## clusterVars[1] <- 0
        for(k in 1:K) {
            if(nAllocsPerCluster[k] == 0) {
                                        #Sample from prior 
                clusterVars[k] <- exp(rnorm(1, mean = m0, sd = s0))
            } else {
                allocInds <- which(clusterAllocations == k)
                ## for(j in 1:d){
                                        #Propose a new (log) variance:
                currentClusterVars     <- clusterVars[k]
                proposedClusterVars    <- exp(log(currentClusterVars) + rnorm(d, sd = 0.5))
                ## ## keep entry 1 < entry 2
                ## if(k<K && proposedClusterVars >= clusterVars[k+1]) 
                ##     proposedClusterVars <- clusterVars[k+1] * 0.9
                ## if(k>1 && proposedClusterVars <= clusterVars[k-1])
                ##     proposedClusterVars <- clusterVars[k-1] * 1.1

                logPriorCurrent     <-
                    dnorm( log(currentClusterVars), mean = m0, sd = s0, log = T)
                logPriorProposed    <-
                    dnorm(log(proposedClusterVars), mean = m0, sd = s0, log = T)
                currentLogLikeli  <-
                    sum(dnorm(obsData[allocInds],
                              mean=clusterMeans[k],
                              sd=sqrt(currentClusterVars + obsVars[allocInds]),
                              log=TRUE))
                proposedLogLikeli <-
                    sum(dnorm(obsData[allocInds],
                              clusterMeans[k],
                              sd=sqrt(proposedClusterVars + obsVars[allocInds]),
                              log=TRUE))
                logMHratio <- proposedLogLikeli+logPriorProposed - 
                                   currentLogLikeli- logPriorCurrent 
                    
                Accept <- log(runif(d)) < logMHratio
                if(Accept) {
                    clusterVars[k]      <- proposedClusterVars
                    varianceAcceptances[k] <- varianceAcceptances[k] + 1 
                }
            }
        }
        clusterVars <- sort(clusterVars)

## ###Update the cluster means:
##         for(k in 2:K) ## don't update null cluster 
##         {
##             if(nAllocsPerCluster[k] == 0)
##             {
##                                         #Sample from prior 
##                 clusterMeans[k] <- rnorm(d, mean = mu0, sd = sigma0)
##             } else{
##                 allocInds <- which(clusterAllocations == k)
##                 currentClusterMean     <- clusterMeans[k]
##                 currentClusterVars=exp(logClusterVars[k])
##                 ## proposedClusterMean    <- currentClusterMean  + rnorm(d, sd = 0.5)
##                 proposedClusterMean    <- currentClusterMean  + rnorm(d, sd = 1) # data limited by (-0.1,0.1), no point jumping about outside this
                
##                 logPriorCurrent     <- dnorm( currentClusterMean, mean = mu0, sd = sigma0, log = T)
##                 logPriorProposed    <- dnorm(proposedClusterMean, mean = mu0, sd = sigma0, log = T)
##                 currentLogLikeli  <-
##                     sum(dnorm(obsData[allocInds],
##                               currentClusterMean,
##                               sd=sqrt(currentClusterVars + obsVars[allocInds])))
##                 proposedLogLikeli <-
##                     sum(dnorm(obsData[allocInds],
##                               proposedClusterMean,
##                               sd=sqrt(currentClusterVars + obsVars[allocInds])))
##                 meanTrials[k] <- meanTrials[k] + 1
##                 MHratio <- exp(proposedLogLikeli+logPriorProposed - currentLogLikeli- logPriorCurrent )
##                 Accept <- runif(d) < MHratio
##                 if(Accept) {
##                     clusterMeans[k]    <- proposedClusterMean
##                     meanAcceptances[k] <- meanAcceptances[k] + 1
##                 }
##             }
##         }

        if(its%%thinner == 0){
            if(quiet)
                pb$tick()
            if(!quiet) 
                message(its, "\t", clusterVars[1],"\t",clusterVars[2],"\t",nAllocsPerCluster[1],"\t",nAllocsPerCluster[2])
            
            ## savedClusterMeans[[counter]]   <- clusterMeans
            savedClusterVars[[counter]]    <- clusterVars
            savedMixtureWeights[[counter]] <- mixtureWeights
            ## savedClusterAllocations[counter,] <- clusterAllocations
            counter <- counter + 1
            if(its %% 50000 == 0)
            {
                print("Saving...")
##                 save(varianceAcceptances,varianceTrials,
##          meanAcceptances,meanTrials,
## savedClusterAllocations, savedClusterMeans, savedMixtureWeights, savedLogClusterVars, file = savefile)
            }
        }
    }
    save(#varianceAcceptances,varianceTrials,
         #meanAcceptances,meanTrials,
        #savedClusterAllocations, #savedClusterMeans,
        savedMixtureWeights,
        savedClusterVars, file = savefile)
}

