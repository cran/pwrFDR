#include <R.h>
#include <Rmath.h>

typedef struct{
  int index;
  double X;
  double pval;
  int HA;
} XPHindxd;

/* double rbinom(double nin, double pp) */
/* s2.power=v.S/(m.1[k]^2*nsim) */

int cmprXPH(const void *, const void *);

void pwrFDRsimCS(int *pnsim, double *pFDR, int *pNg, double *pr1, int *pn, double *ptheta,
		 double *prho, int *pnWC, int *pM, int *pJ, int *pS, double *pX_i, int *pM_i)
{
    int nsim=*pnsim, Ng=*pNg, nWC=*pnWC, m1, ii, j, J_, cS_, n, done, h, k, nC;
    double r1, xN, xnsim, xM1, X_j, FDR=*pFDR, *fdrcrit, theta=*ptheta, rtnth, U, xdf, Z;
    double abs_rho=fabs(*prho), sgn_rho=sign(*prho), tau, sig;
    XPHindxd *pXPH;

    fdrcrit = (double *)Calloc(Ng, double);
    pXPH = (XPHindxd *)Calloc(Ng, XPHindxd);

    GetRNGstate();
    xN = (double)Ng;
    xnsim = (double)nsim;
    r1 = *pr1;
    n = *pn;
    nC = Ng/nWC;
    rtnth = pow(n, 0.5)*theta;
    xdf = (double) 2*n - 2;
    tau = pow(abs_rho, 0.5);
    sig = pow(1.0-abs_rho, 0.5);
    // qnorm5(double p, double mu, double sigma, int lower_tail, int log_p)
    // pnorm5(double x, double mu, double sigma, int lower_tail, int log_p)

    for(h=0;h<nC;h++) 
    for(k=0;k<nWC;k++) *(fdrcrit + nWC*h + k) = FDR*((double)(nWC*h + k+1))/xN;
    for(ii=0;ii<nsim;ii++)
    {
      *(pM+ii) = m1 = (int) rbinom(xN, r1);
      xM1 = (double) m1;
      if(ii==0) *pM_i = m1;
      for(h=0;h<nC;h++)
      {
	/* simulate the Z for this cluster */
	U = unif_rand();
	Z = qnorm5(U, 0.0, tau, 0, 0);
        for(k=0;k<nWC;k++)
        {
	  j = nWC*h + k;
	  (pXPH+j)->index = (j+1);
	  U = unif_rand();
	  X_j = qnorm5(U, 0.0, sig, 0, 0) + Z;
	  if(j < m1) X_j = X_j + rtnth;
          if(ii==0) *(pX_i + j) = X_j;

	  (pXPH+j)->X = X_j;
	  (pXPH+j)->pval = 2.0*pnorm5(fabs(X_j), 0.0, 1.0, 0, 0);
	  (pXPH+j)->HA = ((j+1)<=m1);
	}
      }
      
      qsort(pXPH, Ng, sizeof(XPHindxd), &cmprXPH);
      
      J_=0;
      cS_=0;
      done=0;
      j=0;
      while(!done && j < Ng)
      {
	if((pXPH+Ng-j-1)->pval < (*(fdrcrit + Ng-j-1)))
	{
	  J_ = Ng - j;
	  done = 1;
	}
	j++;
      }
      for(j=0;j<J_;j++) cS_+= (pXPH+j)->HA;
      *(pJ + ii) = J_;
      *(pS + ii) = cS_;
    }

    PutRNGstate();

    Free(pXPH);
    Free(fdrcrit);
}
