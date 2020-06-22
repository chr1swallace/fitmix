library(ggplot2)
library(data.table)
library(cowplot); theme_set(theme_cowplot())
library(seaborn)
library(magrittr)

################################################################################

## stan

################################################################################

files <- list.files("~/scratch",pattern="stanmix1d-.*.csv",full=TRUE)
files
reader <- function(f) {
  ## cat(f,"\t")
  x <- fread(cmd=paste('grep -v "^#"',f))[,.(theta.1,theta.2,sigma.1,sigma.2)]
  x[,f:=f][,c("junk","tissue","chain","junk2"):=tstrsplit(f,"[-_\\.]")]
  setnames(x,c("theta.1","theta.2","sigma.1","sigma.2"),
           c("pi1","pi2","V1","V2"))
  x[,it:=1:.N]
  x[it>50,]
}
stan <- lapply(files,reader)  %>% rbindlist()
stan[,.N,by=c("tissue","chain")]

################################################################################

## custom

################################################################################

source("~/Projects/ase/fitmix_functions.R")
x <- loader("asemix1dnosort-")$x  %>% rbindlist()
  
## print summary
stan[,lapply(.SD,mean),.SDcols=c("pi1","pi2","V1","V2"),by=c("tissue","chain")]
x[,lapply(.SD,mean),.SDcols=c("pi1","pi2","V1","V2"),by=c("tissue","chain")]

## traceplot
m <- melt(stan,c("tissue","chain","it"),c("pi1","pi2","V1","V2"))
p1 <- ggplot(m,aes(x=it,y=value,col=chain)) +
  geom_path() +
  facet_grid(variable~tissue,scales="free") +
  scale_colour_seaborn() +
  background_grid(major="y") +
  ggtitle("stan")

m <- melt(x,c("tissue","chain","it"),c("pi1","pi2","V1","V2"))
p2 <- ggplot(m,aes(x=it,y=value,col=factor(chain))) +
  geom_path() +
  facet_grid(variable~tissue,scales="free") +
  scale_colour_seaborn() +
  background_grid(major="y") +
  ggtitle("custom") 

plot_grid(p1,p2)
