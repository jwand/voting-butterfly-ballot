/*

  RGENOUD (limited version)

  Jasjeet Singh Sekhon 
  Harvard University and Lamarck, Inc.
  http://jsekhon.fas.harvard.edu/
  jsekhon@fas.harvard.edu

  $Header: $

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

int bcacmp(double *a, double *b);

extern "C" 
{

  // <Rdefines.h> must appear inside the {extern "C"} declaration
  // because this header file, unlike <R.h>, does not have an {#ifdef
  // __cplusplus} statment included.
#include <Rdefines.h>

  SEXP fit_lqd1(SEXP I_nparms, SEXP I_components, SEXP I_nobs,
		SEXP I_X, SEXP prep96, SEXP pperot96, SEXP pr1, SEXP pr2, SEXP pr3,
		SEXP pop, SEXP pbuchanan)
  {
    SEXP ret;

    long nparms, ncomponents, nobs, i, j, k, diflen, h, hidx, length, ii;
    double *mu, *y, *rawres, *dif, qrt1, *X, tpop, ty;

    nparms = asInteger(I_nparms);
    ncomponents = asInteger(I_components);
    nobs = asInteger (I_nobs);

//    printf("nobs: %d\n", nobs);

    h = (int) (nobs+nparms+1)/2;
    diflen = (nobs*(nobs-1))/2;
    hidx = (h*(h-1))/2;

    mu = (double *) malloc(nobs*sizeof(double));
    y  = (double *) malloc(nobs*sizeof(double));
    rawres = (double *) malloc(nobs*sizeof(double));
    dif  = (double *) malloc(diflen*sizeof(double));
    X  = (double *) malloc(nparms*sizeof(double));

    for (i=0; i<nparms; i++)
      {
	X[i] = REAL(I_X)[i];

//	printf("X[%d]: \t %lf\n",i+1, X[i]);
      }

    if (ncomponents==1)
      {
	for (i=0; i<nobs; i++)
	  {
	    mu[i] = 
	      X[0] + X[1]*REAL(prep96)[i] + X[2]*REAL(pperot96)[i] + X[3]*REAL(pr1)[i];

//	    printf("mu[%d]: \t %lf\n",i+1, mu[i]);

	    //  y[i] = 1/(1+exp(-1*mu[i]));
	    
	    y[i] = mu[i]>0.0 ? 1/(1+exp(-1*mu[i])) : exp(mu[i])/(1+exp(mu[i]));

//	    printf("y[%d]: \t %lf\n",i+1, y[i]);

	    tpop = REAL(pop)[i];
	    ty   = y[i];

	    rawres[i] = tpop*(REAL(pbuchanan)[i]-ty) / (sqrt(tpop*ty*(1-ty)));

	  }
	
	for (i=2; i<=nobs; i++)
	  {
	    ii = ((i-1)*(i-2))/2;

	    for(j=1; j<=(i-1); j++)
	      {
		k = ii+j;

		dif[k-1] = fabs(rawres[i-1]-rawres[j-1]);

	      }
	  }

	/* R code check
	for (i=2; i<=nobs; i++)
	  {
	    for(j=1; j<=(i-1); j++)
	      {
		k = ((i-1)*(i-2))/2+j;

		dif[k-1] = fabs(rawres[i-1]-rawres[j-1]);

	      }
	  }
	*/
      } // end of ncomponents==2
    else if (ncomponents==2)
      {
	for (i=0; i<nobs; i++)
	  {
	    mu[i] = 
	      X[0] + X[1]*REAL(prep96)[i] + X[2]*REAL(pperot96)[i] + X[3]*REAL(pr1)[i] + X[4]*REAL(pr2)[i];

//	    printf("mu[%d]: \t %lf\n",i+1, mu[i]);

	    //  y[i] = 1/(1+exp(-1*mu[i]));
	    
	    y[i] = mu[i]>0.0 ? 1/(1+exp(-1*mu[i])) : exp(mu[i])/(1+exp(mu[i]));

//	    printf("y[%d]: \t %lf\n",i+1, y[i]);

	    tpop = REAL(pop)[i];
	    ty   = y[i];

	    rawres[i] = tpop*(REAL(pbuchanan)[i]-ty) / (sqrt(tpop*ty*(1-ty)));

	  }
	
	for (i=2; i<=nobs; i++)
	  {
	    ii = ((i-1)*(i-2))/2;

	    for(j=1; j<=(i-1); j++)
	      {
		k = ii+j;

		dif[k-1] = fabs(rawres[i-1]-rawres[j-1]);

	      }
	  }

	/* R code check
	for (i=2; i<=nobs; i++)
	  {
	    for(j=1; j<=(i-1); j++)
	      {
		k = ((i-1)*(i-2))/2+j;

		dif[k-1] = fabs(rawres[i-1]-rawres[j-1]);

	      }
	  }
	*/
      } // end of ncomponents==2
    else if (ncomponents==3)
      {
	for (i=0; i<nobs; i++)
	  {
	    mu[i] = 
	      X[0] + X[1]*REAL(prep96)[i] + X[2]*REAL(pperot96)[i] + X[3]*REAL(pr1)[i] + X[4]*REAL(pr2)[i] + X[5]*REAL(pr3)[i];

//	    printf("mu[%d]: \t %lf\n",i+1, mu[i]);

	    //  y[i] = 1/(1+exp(-1*mu[i]));
	    
	    y[i] = mu[i]>0.0 ? 1/(1+exp(-1*mu[i])) : exp(mu[i])/(1+exp(mu[i]));

//	    printf("y[%d]: \t %lf\n",i+1, y[i]);

	    tpop = REAL(pop)[i];
	    ty   = y[i];

	    rawres[i] = tpop*(REAL(pbuchanan)[i]-ty) / (sqrt(tpop*ty*(1-ty)));

	  }
	
	for (i=2; i<=nobs; i++)
	  {
	    ii = ((i-1)*(i-2))/2;

	    for(j=1; j<=(i-1); j++)
	      {
		k = ii+j;

		dif[k-1] = fabs(rawres[i-1]-rawres[j-1]);

	      }
	  }

	/* R code check
	for (i=2; i<=nobs; i++)
	  {
	    for(j=1; j<=(i-1); j++)
	      {
		k = ((i-1)*(i-2))/2+j;

		dif[k-1] = fabs(rawres[i-1]-rawres[j-1]);

	      }
	  }
	*/
      } // end of ncomponents==3
    else
      {
	printf("ERROR! illegal value of NCOMPONENTS\n");
	qrt1 = 912345;
      }
    

    qsort(dif, diflen, sizeof(double), (int (*)(const void *, const void *)) bcacmp);

    qrt1 = dif[hidx-1];

    if (!R_finite(qrt1))
      {
	printf("XXX\n");
	qrt1 = 999991234;
      }

    //free up ram
    free(X);
    free(dif);
    free(rawres);
    free(y);
    free(mu);

    length=1;
    PROTECT(ret=allocVector(REALSXP,1));
    REAL(ret)[0]=qrt1;
    UNPROTECT(1);
    return(ret);
  } // end of fit_lqd1
} // end of extern "C"


int bcacmp(double *a, double *b) 
{
  int i = 0;

  if (*a > *b) i = 1;
  else if (*a < *b) i = -1;
  return i;
}
