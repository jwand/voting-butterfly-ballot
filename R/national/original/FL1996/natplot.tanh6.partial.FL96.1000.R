#
# LQD estimator: 
#
# Croux, Christophe, Peter J. Rousseeuw and Ola Hossjer.  1994.
# "Generalized S-Estimators."  JASA 89: 1271-1281.
#
# Rousseeuw, Peter J., and Christophe Croux.  1993.  "Alternatives to
# the Median Absolute Deviation," JASA 88: 1273-1283.
#

#let's only do 1 state!
#  FL
redo.states _ c("FL");

# pr3!!!

# uses 1999 race estimates data and the "parms.use.tanh5b.merge" file
# for starting values.

# This runs from an old natplot.tanhX.R run using the "parms.new" file
# (instead of parms.best).  

# This version adds the demographic principle components analysis to
# natplot.tanh.R

#two factor model

# edited ../lqd/natplot.tanh.R so that we read in the best parameter
# solutions from lqd and run genoud on them.  This leads to pure
# replication of the work.  See natplot.tanh2.R for the demographic
# add on.

#
# With comments from Walter.
# We add in the estimation of s0 (sd) with degress of freedom correction
#
# a2) This version saves results across runs (via eval,paste).
#

options(width=150);
options(warn=1)


library("rgenoud");
library("KernSmooth");
source("rewtanh.R");

#parms.use initialization
#parms.use _ read.table("parms.use.tanh8.merge",header=T,
#                       as.is=c(1,2));

#run specific variables
pop.size _ 1000;
wait.generations _ 100;
BFGS _ T;
gradient.check _ F;
max.generations  _ wait.generations*10;

nvars _ 4;
ncomponents _ 1;
runname _ paste("tanh6.partial.FL96.ps",pop.size,".wg",wait.generations,".mg",max.generations,".",
                system("hostname | sed 's/[.].*//g'",T),sep="");

ddir <- paste("postscript.",runname,sep="");
ppos <- T
#par(ask=F);

project _ paste("genoud.pro.natplot.",runname,".ps",pop.size,".wg",
                wait.generations,".",sep="");
project _ paste(project,system("date | sed 's/ //g'",T),sep="");

#create the postscript directory if it does not exist.
parms.use.name    _ paste("parms.use.",runname,sep="");
parms.use.command _ paste(parms.use.name," _ parms.use;",sep="");

county.results.name    _ paste("county.results.",runname,sep="");
county.results.command _ paste(county.results.name," _ county.results;",sep="");


if ( ! file.exists(ddir) )
  {
    dir.create(ddir)
  }

#load up new data
us _ read.csv(file="../data.rev.caci.only.07april01/us.mergedwithnewdata.02april.csv",header=T);
indx _ us$state!="CT";
us _ us[indx,];
indx _ us$state!="NH";
us _ us[indx,];

tmp <- as.factor(us$state)
states <- levels(tmp)
states <- states[ states != "AK"  & states != "CT" & states != "DC" & states != "DE"
                 & states != "HI" & states != "NH" & states != "MI" & states != "RI"]
nstates<- length(states)
#get rid of CT.  The other states have already been removed from us.
indx _ us$state!="CT";
ncounties _ length(us$name[indx]);
indx _ us$state!="NH";
ncounties _ length(us$name[indx]);

#1: state name
#2: sigma
state.results _ as.data.frame(us$gore[1:nstates]);
names(state.results) _ "state";
state.results$state _ states;
state.results$sigma _ data.frame(matrix(0,nrow=nstates,ncol=1));

county.results   _ as.data.frame(us$state);
names(county.results)   _ "state";
county.results$name     _ us$name;
county.results$student  _ matrix(nrow=ncounties,ncol=1);
county.results$standard _ matrix(nrow=ncounties,ncol=1);
county.results$pred     _ matrix(nrow=ncounties,ncol=1);
county.results$actual   _ matrix(nrow=ncounties,ncol=1);

states _  c("AL", "AR", "AZ", "CA", "CO", "FL", "GA", "IA", "ID", "IL",
            "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MN", "MO", "MS", "MT",
            "NC", "ND", "NE", "NJ", "NM", "NV", "NY", "OH", "OK", "OR",
            "PA", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV",
            "WY");

total.states  _ length(states);
parms.use _ as.data.frame(matrix(0,nrow=total.states,ncol=7));
names(parms.use) _ c("state","run","fit","x1","x2","x3","pr1");
parms.use$state _ states;
parms.use$run _ rep("none",total.states);
parms.use$fit _ rep(9999.9,total.states);
parms.use$x1  _ rep(0,total.states);
parms.use$x2  _ rep(0,total.states);
parms.use$x3  _ rep(0,total.states);
parms.use$pr1 _ rep(0,total.states);

print(parms.use);

parm.mat _ as.matrix(cbind(parms.use$x1, parms.use$x2,
                           parms.use$x3, parms.use$pr1));

#state counter
j _ 0;
#county counter
cc _ 1;
obs _ 0;
for (statename in states) 
  {
    j _ j+1;
    cc _ cc+obs;

    if (sum(statename == redo.states))
      {
        SKIP _ F;
        cat("STATE:",statename,"\n");
      }
    else
      {
        SKIP _ T;
        cat("STATE:",statename," skipping! \n");
      }

    indx _ us$state==statename;
    us.tmp _ us[indx,];
    obs _ sum(indx);

    nparms _ nvars;
    h _ ceiling((obs+nparms)/2);
    diflen _ (obs*(obs-1))/2;
    hidx _ (h*(h-1))/2;

    pop       _ us.tmp$vtot96;
    pperot96  _ us.tmp$pperot96;
    prep92    _ us.tmp$prep92;
    pperot92  _ us.tmp$pperot92;

############################################################################
# Principal components for each state...let's add one orthogonal component #
############################################################################

    source("prcomp2_1996.R");

############################################
# Reestimation Module                      #
############################################    

    fit.lqd1 _ function(foo)
      {
        return(.Call("fit_lqd1",as.integer(nparms), as.integer(ncomponents), as.integer(obs),
                     as.vector(foo), as.vector(prep92), as.vector(pperot92),
                     as.vector(pr1), as.vector(pr2), as.vector(pr3),
                     as.vector(pop), as.vector(pperot96)));
      }

    fit.lqd1.reload _ function()
      {
        if (is.loaded("fit.lqd1"))
          {
            dyn.unload("fit_lqd.so");
          }
        dyn.load("fit_lqd.so");
      }


#Do not use glm, use the best fits from "parms.use"

#    if (exists("NEVERDEFINED"))
#      {
        glm0 _ glm(pperot96 ~ prep92 + pperot92,
                   family=binomial(link="logit"), weights=vtot96,
                   data=us.tmp)
        cat("GLM0:\n");
        print(summary(glm0));
        
        glm1 _ glm(pperot96 ~ prep92 + pperot92 + pr1,
                   family=binomial(link="logit"), weights=vtot96,
                   data=us.tmp)
        cat("GLM1:\n");
        print(summary(glm1));
        
        glm2 _ glm(pperot96 ~ prep92 + pperot92 + pr1 + pr2,
                   family=binomial(link="logit"), weights=vtot96,
                   data=us.tmp)
        cat("GLM2:\n");
        print(summary(glm2));

        glm3 _ glm(pperot96 ~ prep92 + pperot92 + pr1 + pr2 + pr3,
                   family=binomial(link="logit"), weights=vtot96,
                   data=us.tmp)
        cat("GLM3:\n");
        print(summary(glm3));    
#      }

    tindx _ parms.use$state==statename;
    starting.values _
      starting.values _ as.vector(glm1$coeff);
    
    cat("Xstarting values:\n");
    print(starting.values);
    
    parms.use$run[tindx] _ "GLM";
    parms.use$fit[tindx] _ 9999.9;
    parms.use[tindx,4:(nvars+3)]  _ as.vector(glm1$coeff);

    if (!SKIP)
      {

    starting.values _
      starting.values _ c(-2.698299e+00, -1.419761e+00, 5.711554e+00, 7.155696e-02);
    
    cat("Xstarting values:\n");
    print(starting.values);
    
    parms.use$run[tindx] _ "TMP";
    parms.use$fit[tindx] _ 3.539749e+00;
    parms.use[tindx,4:(nvars+3)]  _ starting.values;
    
    min.domain.width _ 1;
    mdomains _ matrix(0,nparms,2);
    for (id in 1:nparms) {
      if (id < 4)
        {
          domain.scale _ 5;
        }
      else
        {
          domain.scale _ 10;
        }
      tmp _ domain.scale*max(abs(starting.values[id]),min.domain.width);
      mdomains[id,1] _
        starting.values[id] - tmp;
      mdomains[id,2] _
        starting.values[id] + tmp;
    }    

    fit.lqd1.reload();
    gen1 _ genoud(fit.lqd1,nvars=nvars,max=F,starting.values=starting.values,
                  Domains=mdomains,pop.size=pop.size,
                  wait.generations=wait.generations,
                  max.generations=max.generations,
                  BFGS=BFGS,
                  hard.generation.limit=T,
                  gradient.check=gradient.check,
                  project.path=project);

    min.domain.width _ 1;
    domain.scale _ 5;
    mdomains _ matrix(0,nparms,2);
    for (id in 1:nparms) {
      mdomains[id,1] _
        gen1$par[id] - domain.scale*max(abs(gen1$par[id]),min.domain.width);
      mdomains[id,2] _
        gen1$par[id] + domain.scale*max(abs(gen1$par[id]),min.domain.width);
    }    
    pop.boost _ 5;
    wg.boost _ 3;
    fit.lqd1.reload();
    gen2 _ genoud(fit.lqd1,nvars=nvars,max=F,starting.values=gen1$par,
                  Domains=mdomains,pop.size=pop.boost*pop.size,
                  wait.generations=wg.boost*wait.generations,
                  max.generations=wg.boost*wait.generations*10,
                  BFGS=T,
                  gradient.check=F,
                  project.path=project);

    min.domain.width _ 1;
    domain.scale _ 3;
    mdomains _ matrix(0,nparms,2);
    for (id in 1:nparms) {
      mdomains[id,1] _
        gen1$par[id] - domain.scale*max(abs(gen1$par[id]),min.domain.width);
      mdomains[id,2] _
        gen1$par[id] + domain.scale*max(abs(gen1$par[id]),min.domain.width);
    }

    pop.boost _ 5;
    wg.boost _ 5;
    fit.lqd1.reload();    
    gen3 _ genoud(fit.lqd1,nvars=nvars,max=F,starting.values=gen2$par,
                  Domains=mdomains,pop.size=pop.boost*pop.size,
                  wait.generations=wg.boost*wait.generations,
                  max.generations=wg.boost*wait.generations*10,
                  hard.generation.limit=T,
                  BFGS=BFGS,
                  gradient.check=gradient.check,
                  project.path=project);

    genoud.best _ gen3;

    tmp _ parse(text=paste("gen.",statename,".",pop.size,"_ gen3",sep="")); 
    eval(tmp);
    print(tmp);


#save our best results so far
    if (genoud.best$value < parms.use$fit[j] | parms.use$run=="GLM")
      {
        parm.mat[j,] _ genoud.best$par;
        parms.use$run[j] _ runname;
        parms.use$fit[j] _ genoud.best$value;
        parms.use$x1[j]  _ genoud.best$par[1];
        parms.use$x2[j]  _ genoud.best$par[2];
        parms.use$x3[j]  _ genoud.best$par[3];
        parms.use$pr1[j] _ genoud.best$par[4];
        
        cat(parms.use.command,"\n");
        eval(parse(text=parms.use.command));        
      }
  } # end of !SKIP

############################################
# ESTIMATION OF S0 (SD) WITH DF CORRECTION #
############################################

    residual.generator _ function (foo, x)
      {
        tprep92    _ x[,1];
        tpperot92  _ x[,2];
        tpperot96  _ x[,3];
        tpop       _ x[,4];
        tpr1       _ x[,5];
        
        mu _ foo[1] + foo[2]*tprep92 + foo[3]*tpperot92 +
          foo[4]*tpr1;
        
        y _ 1/(1+exp(-1*mu));

        raw.res _ ( tpop*(tpperot96-y) )/( sqrt(tpop*y*(1-y)) );
        res _ raw.res - median(raw.res)

        return(list(raw.res=raw.res,res=res));
      }

    s0.lqd1 _ function(rawres,nparms)
      {
        reslen _ length(rawres);
        diflen _ (reslen*(reslen-1))/2;
        h _ ceiling((reslen+nparms)/2)
        hidx _ (h*(h-1))/2;
        dif _ abs(outer(rawres,rawres,"-")[outer(1:reslen,1:reslen,">")]);
# dif _ vector(length=diflen,mode="numeric");
# for (i in 2:reslen) {
#   dif[((i-1)*(i-2))/2+(1:(i-1))] _ abs(rawres[i]-rawres[1:(i-1)]);
# }
# qrt1 _ sort(dif)[hidx] * 1/(sqrt(2)*qnorm(5/8));
        qrt1 _ sort(dif)[hidx] * 2.21914446599;

        return(qrt1)
      }

    robustified.leverage _ function (foo, x, w) {
      tprep92    _ x[,1];
      tpperot92  _ x[,2];
      tpperot96  _ x[,3];
      tpop       _ x[,4];
      tpr1       _ x[,5];
      
      mu _ foo[1] + foo[2]*tprep92 + foo[3]*tpperot92 +
        foo[4]*tpr1;;
      
      y _ 1/(1+exp(-1*mu));

      X _ cbind(1,x[,c(1,2,5)]);
      V2 _ diag(1/sqrt(tpop*y*(1-y)));
      V2X _ V2 %*% X;
# see, e.g., McCullagh and Nelder 1989, 397, for nonrobust studentization
# Hdiag _ diag(V2 %*% X %*% solve(t(X) %*% V2 %*% diag(w) %*% V2 %*% X) %*% t(X) %*% V2);
      Hdiag _ diag(V2X %*% solve(t(V2X) %*% diag(w) %*% V2X) %*% t(V2X));

# put the negative forecasting variance adjustment values in Hdiag[w==0]
      Hdiag[w==0] _ -Hdiag[w==0];

      return(Hdiag);
    }

    x _ matrix(data=c(prep92,pperot92,pperot96,us.tmp$vtot96,pr1),
               ncol=(nvars+1));
    residuals _ residual.generator(parm.mat[j,],x);
    print(residuals$raw.res)
    print(median(residuals$raw.res))
    print(residuals$res)
    
    s0  _ s0.lqd1(residuals$raw.res,nparms);

    print(paste("s0: ",s0,sep=""));

# tanh estimator
    print("tanh estimator (using LQD sigma^2)");
    rtanh <-
      rewtanh(parm.mat[j,], s0^2, residuals$res,
              x[,4], x[,3]*x[,4], cbind(1,prep92,pperot92,pr1));
    print(summary(rtanh$glmobj, dispersion=rtanh$disp));
    print("sandwich SEs:");
    print(sqrt(diag(rtanh$rcov)));

#final scale variable
    sigma2 _ rtanh$disp;
    sigma  _ sqrt(sigma2);

    print(paste("sigma: ",sigma,sep=""));

    state.results$sigma[j] _ sigma;

    Hdiag _ robustified.leverage(rtanh$glmobj$coef, x, ifelse(rtanh$w>0,1,0));
    print("Hdiag:");
    print(Hdiag);

#calculate county results
    fn.county.results _ function (foo, x, sigma, Hdiag)
      {
        mu _ foo[1] + foo[2]*x$prep92 + foo[3]*x$pperot92 +
          foo[4]*x$pr1;
        
        y _ 1/(1+exp(-1*mu));
        
        rawres _ ( pop*(x$pperot96-y) )/sqrt(pop*y*(1-y));
        standard _ rawres / sigma;
        student _ rawres / (sigma * sqrt(1-Hdiag));
        
        return(list(pred=y,student=student,standard=standard));
      }
    
    foo _ fn.county.results(rtanh$glmobj$coef, us.tmp, sigma, Hdiag);

    county.results$student[cc:(cc+obs-1)] _ foo$student;
    county.results$standard[cc:(cc+obs-1)] _ foo$standard;
    county.results$pred[cc:(cc+obs-1)]    _ foo$pred;
    county.results$actual[cc:(cc+obs-1)]  _ us.tmp$pperot96;

#update the master county.results matrix
    cat(county.results.command,"\n");
    eval(parse(text=county.results.command));

#print the county results if we are not skipoing
    if (!SKIP)
      {
        cat("BEGIN COUNTY RESULTS FOR ",statename,"\n");
        tindx2 _ county.results$state==statename;
        foo2 _ county.results[tindx2,];
        tindx2 _ order(county.results$student[tindx2]);
        foo2$order[tindx2] _ 1:obs;
        print(foo2[tindx2,c("state","name","student","pred","actual","order")]);        
        cat("END COUNTY RESULTS FOR ",statename,"\n");                       
      }

    if (ppos) {
      postscript(file=paste(ddir,"/",statename,".",parms.use$run[j],".ps",sep=""),
                 width=6,height=8,horizontal=F)
      tmain _ paste("Buchanan outliers in ",statename, "\n");
      par(mfrow=c(2,1));

      plot(foo$pred,foo$student,frame.plot=T,
           xlab = "Expected Proportion",
           ylab = "Studentized Residual")
      title(main=paste(tmain,"Vote Proportions"),
            sub="Proportion versus Studentized Residuals")
      
      qqnorm(foo$student,main="",
             ylab = "Studentized Residual Quantiles",
             sub="Normal Q-Q Plot",frame.plot=T);
      par(mfrow=c(1,1));
      dev.off()
    }
    
  } #end of state loop

##plot the national level analysis 
if (ppos) {
  postscript(file=paste(ddir,"/","national.student.bests.ps",sep=""),
             width=6,height=8,horizontal=F,onefile=F)
  x <- county.results$student;
  #use stdev because it is more robust!
  hh <- dpih(x,scale="stdev")
  bins <- seq(min(x)-0.1,max(x)+0.1+hh,by=hh)
  hist(x,breaks=bins,main="",xlab="Discrepancies from Expected Vote for Buchanan",
       ylab="Number of Reporting Units") #,xlim=c(-10,40))
  title(  main="Studentized");

  dev.off()

  postscript(file=paste(ddir,"/","national.standard.bests.ps",sep=""),
             width=6,height=8,horizontal=F,onefile=F)

  x <- county.results$standard;
  #use stdev because it is more robust!
  hh <- dpih(x,scale="stdev")
  bins <- seq(min(x)-0.1,max(x)+0.1+hh,by=hh)
  hist(x,breaks=bins,main="",xlab="Discrepancies from Expected Vote for Buchanan",
       ylab="Number of Reporting Units") #,xlim=c(-10,40))
  title(  main="Standardized");

  dev.off()
}

#copy the parms.use file to our the new data structure;
cat(parms.use.command,"\n");
#eval(parse(text=paste(parms.use.name," _ parms.use;",sep="")));
eval(parse(text=parms.use.command));
cat(county.results.command,"\n");
eval(parse(text=county.results.command));

#
# print the parms.use file to a file
#

#tfoo _ paste("tindx _ ",parms.use.name,"$state==redo.states;",sep="");
#eval(parse(text=tfoo));

tfoo _ paste("write.table(",parms.use.name,",file=\"",parms.use.name,"\",row.names=F)",sep="");
eval(parse(text=tfoo));

#print the reporting units (only for states which have been run!)
tindx _ county.results$state==redo.states;

county.results.use _ county.results[tindx,];

indx _ order(county.results.use$student);
county.results.use$order[indx] _ 1:ncounties;

cat("\nBEGIN PARMS FILES\n");
eval(parse(text=paste(print(parms.use.name,sep=""))));
cat("END PARMS FILE\n");

cat("\nBEGIN FINAL RESULTS\n");
print(county.results.use[indx,c("state","name","student","pred","actual","order")]);
cat("END FINAL RESULTS\n");
