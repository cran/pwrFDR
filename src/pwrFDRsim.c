#include <R.h>
#include <Rmath.h>

typedef struct{
  int opt;
  double p1;
  double p2;
  double p3;
} distpars;

typedef struct{
  int index;
  double X;
  double pval;
  int HA;
} XPHindxd;

double pdist(double x, distpars *par)
{
  double xncp, sig, xdf1, xdf2, ans;
  if(par->opt == 0)                                         
  {
    xncp=par->p1;
    sig=par->p2;
    ans=pnorm5(x, xncp, sig, 1, 0);
  }                                                         
  if(par->opt == 1)
  {                                                         
    xncp=par->p1;
    xdf1=par->p2;
    if(fabs(xncp)>1e-6) ans=pnt(x,xdf1,xncp,1,0);
    if(fabs(xncp)<=1e-6) ans=pt(x,xdf1,1,0);
  }
  if(par->opt == 2)
  {                                                         
    xncp=par->p1;
    xdf1=par->p2;
    xdf2=par->p3;
    if(fabs(xncp)>1e-6) ans=pnf(x,xdf1,xdf2,xncp,1,0);
    if(fabs(xncp)<=1e-6) ans=pf(x,xdf1,xdf2,1,0);
  }
  return(ans);
}

double qdist(double x, distpars *par)
{
  double xncp, sig, xdf1, xdf2, ans;
  if(par->opt == 0)                                         
  {
    xncp=par->p1;
    sig=par->p2;
    ans=qnorm5(x, xncp, sig, 0, 0);
  }                                                         
  if(par->opt == 1)
  {                                                         
    xncp=par->p1;
    xdf1=par->p2;
    if(fabs(xncp)>1e-6) ans=qnt(x,xdf1,xncp,1,0);
    if(fabs(xncp)<=1e-6) ans=qt(x,xdf1,1,0);
  }
  if(par->opt == 2)
  {                                                         
    xncp=par->p1;
    xdf1=par->p2;
    xdf2=par->p3;
    if(fabs(xncp)>1e-6) ans=qnf(x,xdf1,xdf2,xncp,1,0);
    if(fabs(xncp)<=1e-6) ans=qf(x,xdf1,xdf2,1,0);
  }
  return(ans);
}


// qnt(double p, double df, double ncp, int lower_tail, int log_p)
// pnt(double t, double df, double ncp, int lower_tail, int log_p)
// pt(double x, double n, int lower_tail, int log_p)
// qnorm5(double p, double mu, double sigma, int lower_tail, int log_p)
// pnorm5(double x, double mu, double sigma, int lower_tail, int log_p)

int cmprXPH(const void *, const void *);

void pwrFDRsim(int *pnsim, double *pFDR, int *pNg, double *pr1, int *pn, double *ptheta,
	       int *pdistopt, double *pgroups, int *pverb, int *pM, int *pJ, int *pS, 
	       double *pX_i, int *pM_i)
{
    int nsim=*pnsim, Ng=*pNg, m1, ii, j, J_, cS_, n=*pn, done, idistopt=*pdistopt;
    double r1=*pr1, xN, xnsim, xn, xgroups=*pgroups, TWO, xM1, X_j, FDR=*pFDR, *fdrcrit, theta=*ptheta, rtnth, U, xncp;
    XPHindxd *pXPH;
    distpars *par0, *par1;

    par0 = (distpars *) Calloc(1, distpars);
    par1 = (distpars *) Calloc(1, distpars);
    fdrcrit = (double *)Calloc(Ng, double);
    pXPH = (XPHindxd *)Calloc(Ng, XPHindxd);

    GetRNGstate();
    xN = (double)Ng;
    xnsim = (double)nsim;
    xn = (double)n;
    
    par0->opt  = par1->opt  = idistopt;
    /* Rprintf("nsim: %d, FDR: %g, Ng: %d, r1: %g, n: %d, theta: %g, distopt: %d, groups: %g\n",
               *pnsim, *pFDR, *pNg, *pr1, *pn, *ptheta, *pdistopt, *pgroups); */
    if(idistopt==0)
    {
      xncp = pow(xn/2.0, 0.5)*theta;
      par0->p1 = par1->p1 = 0.0;
      par0->p2 = par1->p2  = 1.0;
    }
    if(idistopt==1)
    {
      xncp = pow(xn/2.0, 0.5)*theta;
      par0->p1 = par1->p1 = 0.0;
      par0->p2 = par1->p2 = 2.0*xn - 2.0;
    }
    if(idistopt==2)
    {
      xncp = (double) xn/2.0*theta*theta;
      par0->p1 = par1->p1 = 0.0;
      par0->p2 = par1->p2 = xgroups-1.0;
      par0->p3 = par1->p3 = xgroups*(xn - 1.0);
    }
    /* Rprintf("xncp: %g, par0->opt: %d, par0->p1: %g, par0->p2: %g, par1->opt: %d, par1->p1: %g, par1->p2: %g\n",
               xncp, par0->opt, par0->p1, par0->p2,  par1->opt, par1->p1, par1->p2); */

    par0->p1 = par1->p1 = 0.0;
    par1->p1 = xncp;
    U = unif_rand();
    X_j = qdist(U, par1);
    /* Rprintf("xncp: %g, par0->opt: %d, par0->p1: %g, par0->p2: %g, par1->opt: %d, par1->p1: %g, par1->p2: %g\n",
	    xncp, par0->opt, par0->p1, par0->p2,  par1->opt, par1->p1, par1->p2);
       Rprintf("U: %g, X_j: %g\n", U, X_j); */
    
    for(j=0;j<Ng;j++) *(fdrcrit + j) = FDR*((double)(j+1))/xN;

    TWO = 2.0;
    if(idistopt==2) TWO = 1.0;
    
    for(ii=0;ii<nsim;ii++)
    {
      *(pM+ii) = m1 = (int) rbinom(xN, r1);
      xM1 = (double) m1;
      if(ii==0) *pM_i = m1;
      for(j=0;j<Ng;j++)
      {
        par0->p1 = par1->p1 = 0.0;
        if(j < m1) par1->p1 = xncp;

        U = unif_rand();
	X_j = qdist(U, par1);
	
        if(ii==0) *(pX_i + j) = X_j;

        (pXPH+j)->index = (j+1);
        (pXPH+j)->X = X_j;
        (pXPH+j)->pval = TWO*(1.0 - pdist(fabs(X_j), par0));
        (pXPH+j)->HA = ((j+1)<=m1);
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

int cmprXPH(const void *x, const void *y)
{
  XPHindxd *xx, *yy;
  xx = (XPHindxd *) x;
  yy = (XPHindxd *) y;

  return(1*((xx->pval) > (yy->pval)) - 1*((xx->pval) < (yy->pval)));
}
