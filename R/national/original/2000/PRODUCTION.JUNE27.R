options(warn=1)
options(width=120);

library("rgenoud");
library("KernSmooth");
source("rewtanh2.R");
#source("data3.R")

#parms.use initialization
parms.use _ read.table("parms.use.ak6.ps1000.wg100.mg1000.lapo",header=T,
                       as.is=c(1,2));

#run specific variables
USEALL  _ F;
PARTIAL _ T;
REVERSE _ F;
pop.size _ 0;
wait.generations _ 0;
gradient.check _ F;
max.generations  _ wait.generations*10;

#if we are a partial run, what states do we do?
redo.states _ c("");


nvars _ 4;
ncomponents _ 1;

runname _ paste("PRODUCTION.JUNE27.");
runname _ paste(runname,"ps",pop.size,".wg",wait.generations,".mg",
                max.generations,".",
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
us _ read.csv(file="../data.rev/us.with2000census.26june.csv",header=T);

tmp <- as.factor(us$state)
states <- levels(tmp)
states <- states[ states != "CT" & states != "DC" & states != "DE" & states != "HI" & states != "MI" &
                  states != "RI"];
nstates<- length(states);
ncounties _ dim(us)[1];

#1: state name
#2: sigma
state.results _ as.data.frame(us$gore[1:nstates]);
names(state.results) _ "state";
state.results$state _ states;
state.results$sigma _ data.frame(matrix(0,nrow=nstates,ncol=1));

if (!PARTIAL)
{
  redo.states _ states;
}

if (REVERSE)
{
  states.order <- states[nstates:1];
} else
{
  states.order <- states;
}

county.results   _ as.data.frame(us$state);
names(county.results)   _ "state";
county.results$name     _ us$name;
county.results$student  _ matrix(nrow=ncounties,ncol=1);
county.results$standard _ matrix(nrow=ncounties,ncol=1);
county.results$pred     _ matrix(nrow=ncounties,ncol=1);
county.results$actual   _ matrix(nrow=ncounties,ncol=1);

print(parms.use);

parm.mat _ as.matrix(cbind(parms.use$x1, parms.use$x2,
                           parms.use$x3, parms.use$pr1));

if (!PARTIAL)
  {
    cat("\nPARTIAL==F.  Doing all States\n");
  } else
{
  cat("\nPARTIAL==T.  Only doing the folling states: ",redo.states,"\n\n");
}

if (!REVERSE)
  {
    cat("\nREVERSE==F.  Proceeding alphabetically\n");
  } else
{
  cat("\nREVERSE==T.  Proceeding REVERSE-alphabetically\n");
}

if (!USEALL)
  {
    cat("\nUSEALL==F.  Only updating if fits are better than parms.use.\n");
  } else
{
  cat("\nUSEALL==T.  Using all obtained fits, even if they are worse than parms.use.\n");
}

cat("\n\n");

#state counter
j _ 0;
#county counter
cc _ 1;
obs _ 0;
for (statename in states.order) 
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

    pop       _ us.tmp$vtot00;
    pbuchanan _ us.tmp$pbuchanan;
    prep96    _ us.tmp$prep96;
    pperot96  _ us.tmp$pperot96;

############################################################################
# Principal components for each state...let's add one orthogonal component #
############################################################################

    source("prcomp4.R");

############################################
# Reestimation Module                      #
############################################    

    fit.lqd1 _ function(foo)
      {
        return(.Call("fit_lqd1",as.integer(nparms), as.integer(ncomponents), as.integer(obs),
                     as.vector(foo), as.vector(prep96), as.vector(pperot96),
                     as.vector(pr1), as.vector(pr2), as.vector(pr3),
                     as.vector(pop), as.vector(pbuchanan)));
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
        glm0 _ glm(pbuchanan ~ prep96 + pperot96,
                   family=binomial(link="logit"), weights=vtot00,
                   data=us.tmp)
        cat("GLM0:\n");
        print(summary(glm0));
        
        glm1 _ glm(pbuchanan ~ prep96 + pperot96 + pr1,
                   family=binomial(link="logit"), weights=vtot00,
                   data=us.tmp)
        cat("GLM1:\n");
        print(summary(glm1));
        
        glm2 _ glm(pbuchanan ~ prep96 + pperot96 + pr1 + pr2,
                   family=binomial(link="logit"), weights=vtot00,
                   data=us.tmp)
        cat("GLM2:\n");
        print(summary(glm2));

        glm3 _ glm(pbuchanan ~ prep96 + pperot96 + pr1 + pr2 + pr3,
                   family=binomial(link="logit"), weights=vtot00,
                   data=us.tmp)
        cat("GLM3:\n");
        print(summary(glm3));    
#      }

#    if(exists("NEVERDEFINED"))
#      {

    tindx _ parms.use$state==statename;
    
#let's use the glm starting values!
    if (parms.use$fit[tindx]!=88888)
      {
        starting.values _
          c(parms.use$x1[tindx],parms.use$x2[tindx],parms.use$x3[tindx],
            parms.use$pr1[tindx]);
      } else {
        starting.values _ glm1$coef;
      }
    
    cat("Xstarting values:\n");
    print(starting.values);
    cat("sv.end\n");

    if (!SKIP)
      {

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
                      gradient.check=gradient.check,
                      project.path=project);

        genoud.best _ gen3;

        tmp _ parse(text=paste("gen.",statename,".",pop.size,"_ gen3",sep="")); 
        eval(tmp);
        print(tmp);


#save our best results so far
        if (!USEALL)
          {
            if (genoud.best$value < parms.use$fit[j])
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
          }
        else
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
#  } #neverdefined

############################################
# ESTIMATION OF S0 (SD) WITH DF CORRECTION #
############################################

    residual.generator _ function (foo, x)
      {
        tprep96    _ x[,1];
        tpperot96  _ x[,2];
        tpbuchanan _ x[,3];
        tpop       _ x[,4];
        tpr1       _ x[,5];
        
        mu _ foo[1] + foo[2]*tprep96 + foo[3]*tpperot96 +
          foo[4]*tpr1;
        
        y _ 1/(1+exp(-1*mu));

        raw.res _ ( tpop*(tpbuchanan-y) )/( sqrt(tpop*y*(1-y)) );
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
      tprep96    _ x[,1];
      tpperot96  _ x[,2];
      tpbuchanan _ x[,3];
      tpop       _ x[,4];
      tpr1       _ x[,5];
      
        mu _ foo[1] + foo[2]*tprep96 + foo[3]*tpperot96 +
          foo[4]*tpr1;
      
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

    x _ matrix(data=c(prep96,pperot96,pbuchanan,us.tmp$vtot00,pr1),
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
              x[,4], x[,3]*x[,4], cbind(1,prep96,pperot96,pr1));
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
        mu _ foo[1] + foo[2]*x$prep96 + foo[3]*x$pperot96 +
          foo[4]*x$pr1;
        
        y _ 1/(1+exp(-1*mu));
        
        rawres _ ( pop*(x$pbuchanan-y) )/sqrt(pop*y*(1-y));
        standard _ rawres / sigma;
        student _ rawres / (sigma * sqrt(1-Hdiag));
        
        return(list(pred=y,student=student,standard=standard));
      }
    
    foo _ fn.county.results(rtanh$glmobj$coef, us.tmp, sigma, Hdiag);

    county.results$student[cc:(cc+obs-1)] _ foo$student;
    county.results$standard[cc:(cc+obs-1)] _ foo$standard;
    county.results$pred[cc:(cc+obs-1)]    _ foo$pred;
    county.results$actual[cc:(cc+obs-1)]  _ us.tmp$pbuchanan;

#update the master county.results matrix
    cat(county.results.command,"\n");
    eval(parse(text=county.results.command));    

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
#par(mfrow=c(1,1));
if (ppos) {
  postscript(file=paste(ddir,"/","national.student.bests.ps",sep=""),
             width=6,height=8,horizontal=F,onefile=F)
  x <- county.results$student;
  #use stdev because it is more robust!
  hh <- dpih(x,scale="stdev")
  bins <- seq(min(x)-0.1,max(x)+0.1+hh,by=hh)
  hist(x,breaks=bins,main="",
       xlab="Studentized Discrepancies from Expected Vote for Buchanan",
       ylab="Number of Counties",freq=T,xlim=c(-14,39));
#  title(  main="Studentized");

  text(34, 200, "Palm Beach, FL" ,cex=0.9)
  xx <- 190; yy <- 3.614137e+01; arrows( yy,xx,yy,4,length=0.1)
  
  text(28, 180,"Jasper, SC",cex=0.9)
  xx <- 170; yy <- 2.825638e+01; arrows( yy,xx,yy,4,length=0.1)

  text(-7.9, 200, "Orleans, LA" ,cex=0.9)
  xx <- 190; yy <- -8.890691e+00; arrows( yy,xx,yy,4,length=0.1)
  
  text(-5, 180,"Caddo\nLA",cex=0.9)
  xx <- 165; yy <- -6.205890e+00; arrows( yy,xx,yy,4,length=0.1)    

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

ttfoo _ paste("write.table(",parms.use.name,",file=\"",parms.use.name,"\",row.names=F)",sep="");
eval(parse(text=ttfoo));

#print the reporting units
indx _ order(county.results$student);
county.results$order[indx] _ 1:ncounties;

cat("\nBEGIN PARMS FILES\n");
eval(parse(text=paste(print(parms.use.name,sep=""))));
cat("END PARMS FILE\n");

cat("\nBEGIN FINAL RESULTS\n");
print(county.results[indx,c("state","name","student","pred","actual","order")]);
cat("END FINAL RESULTS\n");
