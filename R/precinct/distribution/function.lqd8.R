#
# LQD estimator: 
#
# Croux, Christophe, Peter J. Rousseeuw and Ola Hossjer.  1994.
# "Generalized S-Estimators."  JASA 89: 1271-1281.
#
# Rousseeuw, Peter J., and Christophe Croux.  1993.  "Alternatives to
# the Median Absolute Deviation," JASA 88: 1273-1283.
#
# With comments from Walter.
# We add in the estimation of s0 (sd) with degress of freedom correction
#
# a2) This version saves results across runs (via eval,paste).
#


lqd <- function( mydata,nparms=3,str.pop,str.dep,str.indep,sample.name,
                pop.size,wait.generations,max.generations=20,gradient.check,project,
                outlier.threshold=2.5,flag.print=T,ps.dir="postscript",ps.flag=T,
                pstart=NULL,mdomains=NULL) {

  ## constants
  obs    <- dim(mydata)[1]
  h      <- ceiling((obs+nparms)/2);
  diflen <- (obs*(obs-1))/2;
  hidx   <- (h*(h-1))/2;

  ## define the data here
  data.pop   <- eval(parse(text=paste("mydata$",str.pop,sep="")))
  data.dep   <- eval(parse(text=paste("mydata$",str.dep,sep="")))

  ## build matrix of indep vars
  data.indep <- matrix(NA,obs,length(str.indep)+1)
  data.indep[,1] <- 1.0
  for (i in 1:length(str.indep)) {
    statement <- paste("mydata$",str.indep[i],sep="")
    print(statement)
    tmp <- eval(parse(text=statement))
    ##print(tmp)
    ##print(data.indep[,i+1] )
    data.indep[,i+1] <- tmp
  }
  
  fit.lqd1 <- function(foo) {
    mu <- c(data.indep %*% foo)
    y  <- 1/(1+exp(-1*mu));
      
    rawres <- ( data.pop * (data.dep - y) ) / ( sqrt(data.pop * y * (1-y)) );
    dif <- abs(outer(rawres,rawres,"-")[outer(1:obs,1:obs,">")])
    ## dif <- vector(length=diflen,mode="numeric");
    ## for (i in 2:obs) {
    ##   dif[((i-1)*(i-2))/2+(1:(i-1))] <- abs(rawres[i]-rawres[1:(i-1)]);
    ## }
    qrt1 <- sort(dif)[hidx];

#        print(foo,digits=10)
#        print(qrt1,digits=10)

    
    return(qrt1)
  }

  ## estimate starting values using non-robust glm
  ## NOTE: remove intercept (-1) since we include a column of 1's in data.indep


  glm1 <- glm(data.dep ~ -1 + data.indep, family=binomial(link="logit"), weights=data.pop)
  if (flag.print)
    print(summary(glm1));

  if (is.null(mdomains) ) {
    tdomains <- 100;
    cmax <- max(abs(glm1$coeff));
    if (tdomains > cmax*10)
      {
        udomains <- tdomains;
      }
    else
      {
        udomains <- cmax*10;

      }
  }

  if (is.null(pstart) ) {
    gen1 <- genoud(fit.lqd1,
                   nvars=nparms,
                   max=F,
                   starting.values=glm1$coeff,
                   default.domains=udomains,
                   Domains=mdomains,
                   pop.size=pop.size,
                   wait.generations=wait.generations,
                   max.generations=20,
                   hard.generation.limit=T,
                   gradient.check=gradient.check,
                   project.path=project);
  
    min.domain.width <- 1;
    domain.scale     <- 5;
    mdomains         <- matrix(0,nparms,2);
    for (id in 1:nparms) {
      mdomains[id,1] <- gen1$par[id] - domain.scale*max(abs(gen1$par[id]),min.domain.width);
      mdomains[id,2] <- gen1$par[id] + domain.scale*max(abs(gen1$par[id]),min.domain.width);
    }
  
    pop.boost <- 5;
    wg.boost  <- 3;
  
    gen2 <- genoud(fit.lqd1,
                   nvars=nparms,
                   max=F,
                   starting.values= gen1$par,
                   Domains=mdomains,
                   pop.size=pop.boost*pop.size,
                   wait.generations=wg.boost*wait.generations,
                   max.generations=10*max.generations,
                   BFGS=F,
                   gradient.check=F,
                   project.path=project);
  
    min.domain.width <- 1;
    domain.scale <- 3;
    mdomains <- matrix(0,nparms,2);
    for (id in 1:nparms) {
      mdomains[id,1] <- gen2$par[id] - domain.scale*max(abs(gen2$par[id]),min.domain.width);
      mdomains[id,2] <- gen2$par[id] + domain.scale*max(abs(gen2$par[id]),min.domain.width);
    }
    newparm <- gen2$par
    pop.boost <- 5;
    wg.boost <- 5;

    pop.size         <- pop.boost*pop.size
    wait.generations <- wg.boost*wait.generations
    max.generations  <- 10*max.generations
  }

  else {
    newparm <- pstart
  }

  gen3 <- genoud(fit.lqd1,
                 nvars= nparms,
                 max  = F,
                 starting.values = newparm,
                 default.domains = udomains,
                 Domains         = mdomains,
                 pop.size        = pop.size,
                 wait.generations= wait.generations,
                 max.generations = max.generations,
                 gradient.check = gradient.check,
                 project.path   = project);


  finalpar <- gen3$par
  
  ## save all to <<- global scope
  tmp <- parse(text=paste("glm.",sample.name,"<<- glm1",sep=""));
  eval(tmp);
  if (flag.print) {
    print("Saving GLM and RGENOUD Results");
    print(tmp)
  }

  tmp <- parse(text=paste("gen.",sample.name,".",pop.size,"<<- gen3",sep="")); 
  eval(tmp);

  if (flag.print) {
    print(tmp);
  }


  #############################################
  ## ESTIMATION OF S0 (SD) WITH DF CORRECTION #
  #############################################

  residual.generator <- function (foo) {
    mu <- c(data.indep %*% foo)
    y  <- 1/(1+exp(-1*mu));

    rawres <- ( data.pop * (data.dep - y) ) / ( sqrt(data.pop * y * (1-y)) );
    res <- rawres - median(rawres)

    return(list(raw.res=rawres, res=res));
  }

  s0.lqd1 <- function(rawres,nparms)
    {
      reslen <- length(rawres);
      diflen <- (reslen*(reslen-1))/2;
      h <- ceiling((reslen+nparms)/2)
      hidx <- (h*(h-1))/2;
      dif <- abs(outer(rawres,rawres,"-")[outer(1:reslen,1:reslen,">")]);
      ## dif <- vector(length=diflen,mode="numeric");
      ## for (i in 2:reslen) {
      ##   dif[((i-1)*(i-2))/2+(1:(i-1))] <- abs(rawres[i]-rawres[1:(i-1)]);
      ## }
      ## qrt1 <- sort(dif)[hidx] * 1/(sqrt(2)*qnorm(5/8));
      qrt1 <- sort(dif)[hidx] * 2.21914446599;

      return(qrt1)
    }

  residuals <- residual.generator(gen3$par);
  s0  <- s0.lqd1(residuals$raw.res,nparms);

  if (flag.print) {
    print("residuals$raw.res")
    print(residuals$raw.res)
    cat("median:",median(residuals$raw.res),"\n")
    print("residuals$res")
    print(residuals$res)
    print(paste("s0: ",s0,sep=""));
  }
    

  ##weight variable
  w <- vector(length=obs,mode="numeric");
  indx <- abs(residuals$res/s0) < outlier.threshold;  
  w[indx] <- 1;

  ##final scale variable
  sigma2 <- sum(w * residuals$res^2) / (sum(w)-nparms) ;
  sigma  <- sqrt(sigma2);

  if (flag.print)
    cat("sigma: ",sigma,"\n");


  ##calculate county results
  fn.county.results <- function (foo, sigma)  {
    mu <- c(data.indep %*% foo)
    y  <- 1/(1+exp(-1*mu));

    rawres <- ( data.pop * (data.dep - y) ) / ( sqrt(data.pop * y * (1-y)) );
    student <- (rawres-median(rawres)) / sigma;
        
    return(list(pred=y,student=student));
  }
    
  foo <- fn.county.results(gen3$par,sigma);

  if (flag.print) {
    print("data: student, pred, dep")
    print(cbind(foo$student,foo$pred,data.dep))
  }

  
  if (ps.flag) {
    postscript(file=paste(ps.dir,"/",sample.name,".ps",pop.size,".wg",wait.generations,".ps",sep=""),
               width=6,height=8,horizontal=F)
      tmain <- paste(str.dep,"outliers in ",sample.name,"\n");
      par(mfrow=c(2,1));

      plot(foo$pred,foo$student,frame.plot=T,
           xlab = "Expected Proportion",
           ylab = "Standardized Residual")
      title(main=paste(tmain,"Vote Proportions"),
            sub="Proportion versus Studentized Residuals")
      
      qqnorm(foo$student,main="",
             ylab = "Standardized Residual Quantiles",
             sub="Normal Q-Q Plot",frame.plot=T);
      par(mfrow=c(1,1));
      dev.off()
  }

  ##  state.results$sigma[j] <- sigma;
  ##  county.results$student[cc:(cc+obs-1)] <- foo$student;
  ##  county.results$pred[cc:(cc+obs-1)]    <- foo$pred;
  ##  county.results$actual[cc:(cc+obs-1)]  <- deppercvote;

  return(list( w=w, s0=s0, sigma=sigma, student.res=foo$student, pred=foo$pred,
              dep=data.dep, indep=data.indep,
              stand.res=residuals$res, raw.res= residuals$raw.res, par=finalpar,
              domains=mdomains,pop=data.pop) )
  
} ## end of lqd()

