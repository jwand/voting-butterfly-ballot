options(warn=1)
library("rgenoud")
source("../function.lqd8.R",echo=T,max.deparse.length = 1000)
library("KernSmooth");

source("../data.pbc.R",echo=T,max.deparse.length = 1000) ## MySQL retrieves data
vote <- pbc

indx <- as.logical(vote$useme.ss35)
str.dep         <- "pbuch"
outname         <- "ss35.buch"
initial.parms   <- c(-5.677547e+00,1.943731e+00,-5.453711e+01) ## from round3/v1 (pop 1500)
initial.domains <- NULL

#run specific variables
pop.size         <- 100;
wait.generations <- 1;
gradient.check   <- F;
max.generations  <- wait.generations*10;

str.pop   <- "ptotal"
str.indep <- c("pnelson","pdeckard")
statename <- "PBC"

## where we are going to store things...
ps.flag <- T
ps.dir  <- paste("postscript.lqd8.",str.dep,sep="");
project <- paste("genoud.pro.lqd8.",str.dep,".ps",pop.size,".wg",wait.generations,".",sep="");
project <- paste(project,system("date | sed 's/ //g'",T),sep="");
## ... and create the postscript directory if it does not exist.
if ( ! file.exists(ps.dir) )
  dir.create(ps.dir)

## loop over regions/states/samples
tmp.str.indep <- str.indep
tmp.out <- 
  lqd(mydata=vote[indx,],
      nparms = 3,
      str.pop,
      str.dep,
      tmp.str.indep,
      sample.name     = statename ,
      pop.size        = pop.size,
      wait.generations= wait.generations,
      max.generations = max.generations,
      gradient.check  = gradient.check,
      project   = project,
      outlier.threshold = 2.5,
      flag.print= T,
      ps.dir    = ps.dir ,
      ps.flag   = T,
      pstart   = initial.parms,
      mdomains = initial.domains 
      )

hist(tmp.out$student)
eval(parse(text=paste(outname,".",pop.size,".",wait.generations,"<- tmp.out",sep="")))
