setwd("~/Projects/ase")
library(magrittr)
library(data.table)
library(randomFunctions)
args <- getArgs(defaults=list(t1="blood",
                              id=99,
                              k=2,
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

library(rstan)
options(mc.cores = 4)

input_data=list(N=nrow(x),
                K=args$k,
                obsData=x$slope,
                obsVars=x$slope_se^2)

rstan_options(auto_write = TRUE)
fit <- stan(file='stan_1d_v2.stan', data=input_data,
            seed=args$id,
            iter=10000, warmup=200,sample_file=paste0("/home/cew54/scratch/stanmix1d-",args$t1),
            chains=4,thin=args$thin)

print(fit)
save(fit, file=paste0("~/scratch/stanmix1d-",args$t1,".RData"))
     
if(!interactive())
  q("no")
