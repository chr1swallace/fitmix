library(data.table)

library(rstan)
options(mc.cores = 3)
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
m <- m[sample(1:nrow(m),1000)] # temporarily small
obsData=m[,.(slope.x,slope.y)]  %>% as.matrix()
obsVars=m[,.(slope_se.x^2,slope_se.y^2)]  %>% as.matrix()


input_data=list(N=nrow(m),
                K=4,
                obsData=obsData,
                obsVars=obsVars)

rstan_options(auto_write = TRUE)
fit <- stan(file='stan_2d.stan', data=input_data,
            init=function(chain_id) {
              list(theta=c(0.7,0.1,0.1,0.1),
                   sigma=exp(rnorm(2,0,1)),
                   rho=runif(4)^5)
            },
            iter=100, warmup=50,sample_file="stan_2d_test",
            chains=1)

print(fit)
save(fit, file="stan_2d_test.RData")
