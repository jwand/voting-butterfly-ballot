# reweighted estimator, for overdispersed grouped binomial GLM

rewtanh _ function(bstart, sigma2start, resstart, nobsvec, y, RHSmat, itmax=100) {
 hessianfunc _
  function(xmat,bvec,Nyes,Nno,w) {
    nu _ c(xmat %*% bvec);
    phat _ 1/(1+exp(-nu));
    return( t(xmat) %*% (w * (-(Nyes + Nno) * phat*(1-phat)) * xmat) )
  }
 scorefunc _
  function(xmat,bvec,Nyes,Nno) {
    nu _ c(xmat %*% bvec);
    phat _ 1/(1+exp(-nu));
    return( -(Nyes * (1-phat) - Nno * phat) * xmat )
  }
 resfunc _
  function(xmat,bvec,Nyes,Nall) {
    nu _ c(xmat %*% bvec);
    phat _ 1/(1+exp(-nu));
    return( (Nyes - Nall*phat)/sqrt(Nall*phat*(1-phat)) )
  }
 psifunc _
  function(arg) {
  # Hampel, Rousseeuw and Ronchetti 1981.  constants are from Table 2, p. 645
  #                                                                          effic.
  # c _ 3.0;    k _ 5.0;    A _ 0.680593;    B _ 0.769313;    d _ 1.470089;  # 87%
    c _ 4.0;    k _ 5.0;    A _ 0.857044;    B _ 0.911135;    d _ 1.803134;  # 97%
    return(
      ifelse(abs(arg)<d, arg,
        ifelse(abs(arg)<c,
          sqrt(A*(k-1))*tanh(sqrt((k-1)*B^2/A)*(c-abs(arg))/2)*sign(arg), 0))
    )
  }
 nsfunc _
  function(opg,score) {
    iopg _ solve(opg)
    return( sqrt(apply(score,1,function(x,iA){ t(x) %*% iA %*% x },iA=iopg)) )
  }
 converged _
  function(bnew,bold) {
    return( !any((bnew-bold)/(bold + 1e-6) > 1e-6) )
  }
 robustified.leverage _
  function (X, bvec, Nall, w) {
    nu _ c(X %*% bvec);
    phat _ 1/(1+exp(-nu));
    V2 _ diag(1/sqrt(Nall*phat*(1-phat)));
    V2X _ V2 %*% X;
    # see, e.g., McCullagh and Nelder 1989, 397, for nonrobust studentization
  # Hdiag _
  # diag(V2 %*% X %*% solve(t(X) %*% V2 %*% diag(w) %*% V2 %*% X) %*% t(X) %*% V2);
    Hdiag _ diag(V2X %*% solve(t(V2X) %*% diag(w) %*% V2X) %*% t(V2X));
 
    # put the negative forecasting variance adjustment values in Hdiag[w==0]
    Hdiag[w==0] _ -Hdiag[w==0];
 
    return(Hdiag);
 }
 nregs _ dim(RHSmat)[2];
 bvec _ bstart;
 sigma2 _ sigma2start;
 res _ resstart;
 sres _ res/sqrt(sigma2);
 ipsi _ psifunc(abs(sres));
 w _ ifelse(sres==0,1,ipsi/abs(sres));
 yprop _ y/nobsvec;
 ny _ nobsvec-y;
 iters _ 0;
 for (i in 1:itmax) {
  iters _ iters + 1;
  wprev _ w;
  bprev _ bvec;
  s2prev _ sigma2;
  # grouped binomial
  iglm _ glm(yprop ~ -1 + RHSmat, subset = w > 0, start=bvec,
             family=binomial(link="logit"), weights=w*nobsvec);
  bvec _ iglm$coef;
#  sigma2 _ sum(iglm$weights * iglm$residuals^2)/(sum(w)-nregs)
  res _ resfunc(RHSmat,bvec,y,nobsvec);
  sres _ res/sqrt(sigma2);
  ipsi _ psifunc(abs(sres));
  w _ ifelse(sres==0,1,ipsi/abs(sres));
  if (converged(bvec,bprev) & converged(sigma2,s2prev)) break;
 }
 score _ scorefunc(RHSmat, bvec, y, ny) / sqrt(sigma2);
 Aopg _ t(w * score) %*% score ;
 infomat _ hessianfunc(RHSmat, bvec, y, ny, w) / sigma2 ;
 invB _ solve(infomat);
 rcovmat _ invB %*% Aopg %*% invB ;
 h _ robustified.leverage(RHSmat, bvec, nobsvec, w) ;
 return(
  list(glmobj=iglm, dispersion=sigma2, w=w, psi=ipsi, res=res, sres=sres, h=h,
       A=Aopg, B=infomat, rcovmat=rcovmat, iters=iters)
 )
}
