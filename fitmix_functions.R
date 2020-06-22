
proc.weights <- function(obj,chain,thin=100) {
  ret <- do.call("rbind",obj)  %>%
    cbind(.,chain=chain)  %>%
    as.data.frame()  %>%
    as.data.table()
  ret[,it:=thin*(1:.N)]
  ret[!is.na(ret[[1]])]
}

loader <- function(pattern="asemix2dv2") {
  files <- list.files("~/scratch",pattern=pattern,full=TRUE)
  ss <- tstrsplit(files,"-")
  vars <- weights <- params <- comb <- vector("list",length(files))
  for(i in seq_along(files)) {
    message(files[i])
    obj <- load(files[i])
    if("savedClusterVars" %in% obj)
      savedParams <- savedClusterVars
    nulls <- sapply(savedParams,is.null)
    use <- which(!nulls)
    params[[i]] <- proc.weights(savedParams[use],i)
## vars[[i]] <- lapply(savedParams[use],calcClusterParams) ## %>%
##  lapply(., function(p) cbind(p[,1:2],p[,2]*p[,3]))
    weights[[i]] <- proc.weights(savedMixtureWeights[use],i)
    setnames(weights[[i]],c("V1","V2"),c("pi1","pi2"))
    comb[[i]] <- cbind(weights[[i]],params[[i]][,.(V1,V2)])[,tissue:=ss[[2]][i]]
  }
  return(list(x=comb,params=params,weights=weights,#vars=vars,
              files=basename(files)))
}

library(GGally)
library(seaborn)

gettails <- function(x,n) 
  lapply(x, function(y) tail(y,-n))  %>% rbindlist()

lowerfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_point(aes(col=pi1 < pi4), alpha=0.2)+
    geom_density2d() +
    geom_abline(col="grey") +
    scale_colour_seaborn()
}  
upperfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_density2d()
}
sortpairs <- function(results,what=c("pi1","pi4","v1","v2"),
                      do.sort=TRUE,do.pairs=TRUE,do.trace=FALSE) {
  nm <-deparse(substitute(results))
  weights <- gettails(results$weights, 1000)
  setnames(weights,setdiff(colnames(weights),c("it","chain")),c("pi1","pi2","pi3","pi4"))
  params <- gettails(results$params, 1000)
  params  <- params[,1:6]
  setnames(params,c("v1","v2","r1","r2","r3","r4"))
  x=cbind(weights,params)
## print(head(x))
  if(do.sort) {
  sw <- x$v2 < x$v1
  x[sw,c("pi4","pi1","v1","v2","r4","r1","r2","r3"):=list(pi1,pi4,v2,v1,r1,r4,r3,r2)]
  }
  if(do.pairs){
   print(ggpairs(x[,what,with=FALSE],
                 ## upper = list(continuous = wrap(upperfun)),
                 lower = list(continuous = wrap(lowerfun))) +
         ggtitle(nm)) 
}
  if(do.trace){
    m <- melt(x, c("it","chain"), what)
    print(ggplot(m,aes(x=it,y=value,col=factor(chain),group=chain)) +
         geom_path(alpha=0.2) +
         geom_smooth(se=FALSE) +
         scale_colour_seaborn("Chain") +
         facet_wrap(~variable) + #,scales="free_y") +
         ylim(0,1) +
         labs(x="Iteration") +
         background_grid(major="y") +
         theme(legend.position="none") +
         ggtitle(nm)) 
  }
  colMeans(x[,.(pi1,pi2,pi3,pi4,v1,v2,r1,r2,r3,r4)])
}


## pairsplot <- function(results) {
##    nm <-deparse(substitute(results))
##    params <- rbindlist(results$params)
##   weights <- rbindlist(results$weights)
##     ggpairs(weights[,1:4])
## }

plotter <- function(results) {
   nm <-deparse(substitute(results))

   ## traceplot mixing weights
   weights <- rbindlist(results$weights)
colMeans(weights)
weights %<>% melt(.,c("chain","it"))
weights[,variable:=sub("V","group",variable)]
p1 <- ggplot(weights,aes(x=it,y=value,col=factor(chain),group=chain)) +
  geom_path(alpha=0.2) +
  geom_smooth(se=FALSE) +
  scale_colour_seaborn("Chain") +
  facet_wrap(~variable,scales="free_y") +
  labs(x="Iteration") +
  background_grid(major="y") +
  theme(legend.position="none") +
  ggtitle("Mixing weights for each group")

## traceplot params
params <- rbindlist(results$params)
colMeans(params)
## params[,V3:=V3/V2]
params %<>% melt(.,c("chain","it"))
params[,variable:=sub("V","param ",variable)]
p2 <- ggplot(params,aes(x=it,y=value,col=factor(chain),group=chain)) +
  geom_path(alpha=0.2) +
  geom_smooth(se=FALSE) +
  scale_colour_seaborn("Chain") +
  facet_wrap(~variable,scales="free_y") +
  labs(x="Iteration") +
  background_grid(major="y") +
  theme(legend.position="none") +
  ggtitle("Variance and correlation parameters")

plot_grid(p1,p2,rel_heights=c(6,4)) + ggtitle(nm)
}


sumtail <- function(x) {
  tmp <- lapply(x, function(xx) tail(xx,nrow(xx)/2))  %>%
    rbindlist()  
  tmp <- rbind(apply(tmp,2,median),colMeans(tmp),apply(tmp,2,sd),apply(log(tmp),2,mean),apply(log(tmp),2,sd))  %>% as.data.frame()
  tmp$stat <- c("median","mean","sd","logmean","logsd")
  tmp
}

