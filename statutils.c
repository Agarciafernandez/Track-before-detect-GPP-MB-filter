#include "mex.h"
#include "matrix.h"
#include "statutils.h"
#include "matrixutils.h"
#include <math.h>
#include <errno.h>
#include <float.h>

double lntrectprob2(double a, double b, double v)
{
  double lnP, z1, z2;
  double alpha, beta, lnP1, lnP2;
  
  z1 = (2*(a>0)-1)*(v-0.6667+0.1/v)*sqrt(log(1+a*a/v)/(v-0.8333));   
  z2 = (2*(b>0)-1)*(v-0.6667+0.1/v)*sqrt(log(1+b*b/v)/(v-0.8333));   
  //lnP = log(normcdf(z2,0,1)-normcdf(z1,0,1));
  
  alpha = 0.344;
  beta = 5.334;
  lnP1 = -0.5*z1*z1-log(sqrt(2*PI)*((1-alpha)*fabs(z1)+alpha*sqrt(z1*z1+beta)));
  lnP2 = -0.5*z2*z2-log(sqrt(2*PI)*((1-alpha)*fabs(z2)+alpha*sqrt(z2*z2+beta)));
  
  if (z1<0)
  {
    if (z2<0) lnP = lnP2+log(1-exp(lnP1-lnP2));
    else lnP = log(1-exp(lnP1)-exp(lnP2));
  }
  else lnP = lnP1+log(1-exp(lnP2-lnP1));
  
  return lnP;
}

double lntrectprob(double a, double b, double v, int nreal)
{
  double muq, sigmaq, *w, wmax, addterm, lnP;
  double lnsigmaq, invsigmaqsqrd, xsamp, lngv;
  int i, cnt;
  
  w = (double *) mxCalloc(nreal,sizeof(double));
  
  muq = 0.5*(b+a);
  sigmaq = (b-a)/4;
  invsigmaqsqrd = 0.5/(sigmaq*sigmaq);
  lnsigmaq = -0.5*log(v*invsigmaqsqrd);
  lngv = lngamma(0.5*(v+1))-lngamma(0.5*v);
  
  cnt = 0;
  wmax = -DBL_MAX;
  for (i=0; i<nreal; i++)
  {
    xsamp = normrnd(muq,sigmaq);
    if (xsamp>a && xsamp<b)
    {
      w[cnt] = lnsigmaq+lngv-0.5*(v+1)*log(1+xsamp*xsamp/v)+invsigmaqsqrd*(xsamp-muq)*(xsamp-muq);
      if (w[cnt]>wmax) wmax = w[cnt];
      cnt++;
    }
  }
  
  addterm = 0;
  for (i=0; i<cnt; i++) addterm += exp(w[i]-wmax);
  lnP = wmax+log(addterm/((double) nreal));
  
  mxFree(w);
  
  return lnP;
}

double lnnormrectprob(double a, double b, int nreal)
{
  double muq, sigmaq, *w, wmax, addterm, lnP;
  double lnsigmaq, invsigmaqsqrd, xsamp;
  int i, cnt;
  
  w = (double *) mxCalloc(nreal,sizeof(double));
  
  muq = 0.5*(b+a);
  sigmaq = (b-a)/2;
  invsigmaqsqrd = 0.5/(sigmaq*sigmaq);
  lnsigmaq = log(sigmaq);
  cnt = 0;
  wmax = -DBL_MAX;
  for (i=0; i<nreal; i++)
  {
    xsamp = muq+sigmaq*genonerandn();
    if (xsamp>a && xsamp<b)
    {
      w[cnt] = lnsigmaq-0.5*xsamp*xsamp+invsigmaqsqrd*(xsamp-muq)*(xsamp-muq);
      if (w[cnt]>wmax) wmax = w[cnt];
      cnt++;
    }
  }
  
  addterm = 0;
  for (i=0; i<cnt; i++) addterm += exp(w[i]-wmax);
  lnP = wmax+log(addterm/((double) nreal));
  
  mxFree(w);
  
  return lnP;
}

double lnbinormprob(double *lbnds, double *ubnds, double *G, int nsamp)
{
   // Compute bivariate normal rectangle probs
  
  double y1, d1, c1, d2, c2, u, np1, *w, wmax, addterm;
  double muq, sigmaq, *lnv, *v, lnvmax, vsum;
  int i;
  double lnP;
  
  w = (double *) mxCalloc(nsamp,sizeof(double));
  v = (double *) mxCalloc(nsamp,sizeof(double));
  lnv = (double *) mxCalloc(nsamp,sizeof(double));

  d1 = ubnds[1]/G[0];
  c1 = lbnds[1]/G[0];
  muq = 0.5*(d1+c1);
  sigmaq = (d1-c1)/2;
  np1 = lnnormrectprob(c1,d1,100);
  lnvmax = -DBL_MAX;
  for (i=0; i<nsamp; i++)
  {
    y1 = muq+sigmaq*genonerandn();//normrnd(muq,sigmaq);
    if (y1>c1 & y1<d1) lnv[i] = log(sigmaq)-0.5*y1*y1+0.5*(y1-muq)*(y1-muq)/(sigmaq*sigmaq);
    else  lnv[i] = -DBL_MAX;
    if (lnv[i]>lnvmax) lnvmax = lnv[i];
  }
  vsum = 0;
  for (i=0; i<nsamp; i++)
  {
    v[i] = exp(lnv[i]-lnvmax);
    vsum += v[i];
  }
  wmax = -DBL_MAX;
  for (i=0; i<nsamp; i++)
  {
    lnv[i] = log(v[i]/vsum);
    d2 = (ubnds[0]-G[1]*y1)/G[3];
    c2 = (lbnds[0]-G[1]*y1)/G[3];
    w[i] = lnv[i]+lnnormrectprob(c2,d2,100);
    if (w[i]>wmax) wmax = w[i];
  }
  addterm = 0;
  for (i=0; i<nsamp; i++) addterm += exp(w[i]-wmax);
  lnP = np1+wmax+log(addterm);

  mxFree(w);
  mxFree(v);
  mxFree(lnv);
  
  return lnP;
  
}

double binormprob(double *lbnds, double *ubnds, double *Sigmasqrt, int nsamp)
{
  double y1, d1, c1, d2, c2, u, ncd1, ncc1;
  int i;
  double P;

  d1 = ubnds[0]/Sigmasqrt[0];
  c1 = lbnds[0]/Sigmasqrt[0];
  ncd1 = normcdf(d1,0,1);
  ncc1 = normcdf(c1,0,1);
  P = 0;
  for (i=0; i<nsamp; i++)
  {
    u = unirnd(0,1);
    y1 = norminv(u*ncd1+(1-u)*ncc1);
    d2 = (ubnds[1]-Sigmasqrt[1]*y1)/Sigmasqrt[3];
    c2 = (lbnds[1]-Sigmasqrt[1]*y1)/Sigmasqrt[3];
    P += (normcdf(d2,0,1)-normcdf(c2,0,1));
  }
  P *= (ncd1-ncc1)/((double) nsamp);
  
  return P;
}
  
double gamrnd(double alpha, double beta)
{
  // Uses Cheng's method to generate a Gamma rv with shape alpha and scale beta
  
  double u1, u2, v, a, b, c, x;
  
  a = 1/sqrt(2*alpha-1);
  b = alpha-log(4);
  c = alpha+1/a;
  
  u1 = unirnd(0,1);
  u2 = unirnd(0,1);
  v = a*log(u1/(1-u1));
  x = alpha*exp(v);
  
  while (b+c*v-x<log(u1*u1*u2))
  {
    u1 = unirnd(0,1);
    u2 = unirnd(0,1);
    v = a*log(u1/(1-u1));
    x = alpha*exp(v);
  }
  
  return x*beta;
  
}

int poissrnd(double lambda)
{
  double a, b, u;
  int x;
  
  a = exp(-lambda);
  b = a;
  u = unirnd(0,1);
  x = 0;
  while (u>a)
  {
    x++;
    b *= lambda/((double) x);
    a += b;
  }
  
  return x;
}

int poissrnd2(double lambda)
{
  double u, b, p, x;
  
  b = exp(-lambda);
  x = 0;
  p = 1;
  
  do
  {
    x++;
    u = unirnd(0,1);
    p *= u;
  }
  while (p>b);

  return (x-1);
  
}

double poisspdf(int x, double lambda)
{
  double p;
  
  p = pow(lambda,x)*exp(-lambda)/factorial(x);
  
  return p;
}

double vmrnd(double mu, double kappa)
{
  int accept, cnt=0;
  double tau, rho, r;
  double u1, u2, z, f, c, x;
  
  tau = 1+sqrt(1+4*kappa*kappa);
  rho = (tau-sqrt(2*tau))/(2*kappa);
  r = (1+rho*rho)/(2*rho);
  
  accept = 0;
  cnt = 0;
  x = 0;
  while (accept==0 & cnt<100)
  {
    cnt++;
    u1 = unirnd(0,1);
    u2 = unirnd(0,1);
    z = cos(PI*u1);
    f = (1+r*z)/(r+z);
    c = kappa*(r-f);
    if ((c*(2-c)-u2)>0)
    {
      x = sign(unirnd(0,1)-0.5)*acos(f);
      accept = 1;
    }
    else if ((log(c/u2)+1-c)>0)
    {
      x = sign(unirnd(0,1)-0.5)*acos(f);
      accept = 1;
    }
  }
 
  x += mu;
  x = mod(x+PI,2*PI)-PI;
    
  return x;
}

double unirnd(double a, double b)
{
  double x;
  
  x =(b-a)*rand()/((double) RAND_MAX +1)+a;

  return x;
}

double normrnd(double mu, double sigma)
{
    double x1, x2, y, w;

    y = 2.0;
    while (y>=1.0 || y==0)
    {
        x1 = unirnd(-1,1);
        x2 = unirnd(-1,1);
        y = x1*x1+x2*x2;
     }
     w = x1*sqrt(-2.0*log(y)/y);
     w = mu+sigma*w;
     return w;
}

double trnd(double v)
{
  double a, c, e, vp, u1, u2, x, xsqrd;
  int isok;

  c = -0.25*(v+1);
  a = 4*exp(c*log(v/(1+v)));
  e = 16/a;
  vp = sqrt(2*v)*exp(0.25*((v-1)*log(v-1)-(v+1)*log(v+1)));
  isok=0;
  while (!isok)
  {
    u1 = unirnd(0,1);
    u2 = unirnd(-vp,vp);
    x = u2/u1;
    xsqrd = x*x;
    if ((a*u1)<(5-xsqrd)) isok=1;
    else
    {
      if (xsqrd<(e/u1-3))
        if (u1<exp(c*log(1+xsqrd/v))) isok = 1;
    }
  }
  
  return x;
}



double genonerandn()
{
    double x1, x2, y, w;

    y = 2.0;
    while (y>=1.0)
    {
        x1 = unirnd(-1,1);
        x2 = unirnd(-1,1);
        y = x1*x1+x2*x2;
     }
     w = x1*sqrt(-2.0*log(y)/y);
     return w;
 }

double vmpdf(double x, double mu, double kappa)
{
  double p, sigmasq, ep;
  
  if (kappa<200)  p = exp(kappa*cos(x-mu))/(2*PI*bessi0(kappa));
  else
  {
    sigmasq = 1./kappa;
    ep = mod(x-mu+PI,2*PI)-PI;
    p = exp(-ep*ep/(2*sigmasq))/sqrt(2*PI*sigmasq);
  }
  
  return p;
}

double normpdf(double x, double mu, double sigma)
{
  double p;
  
  p = exp(-(x-mu)*(x-mu)/(2*sigma*sigma))*OneUponSqrt2Pi/sigma;
  return p;
}

double exppdf(double x, double lambda)
{
  double p;
  
  p = lambda*exp(-lambda*x);
 
  return p;
}

double norminv(double p)
{
	double q, r;
  /* Coefficients in rational approximations. */
  double a[] =
  {-3.969683028665376e+01,
  2.209460984245205e+02,
  -2.759285104469687e+02,
  1.383577518672690e+02,
  -3.066479806614716e+01, 
  2.506628277459239e+00};
  double b[] =
  {-5.447609879822406e+01,
	 1.615858368580409e+02,
	-1.556989798598866e+02,
	 6.680131188771972e+01,
	-1.328068155288572e+01};
  double c[] =
  {-7.784894002430293e-03,
	-3.223964580411365e-01,
	-2.400758277161838e+00,
	-2.549732539343734e+00,
	 4.374664141464968e+00,
	 2.938163982698783e+00};
  double d[] =
  {7.784695709041462e-03,
	3.224671290700398e-01,
	2.445134137142996e+00,
	3.754408661907416e+00};
	errno = 0;

	if (p < 0 || p > 1)
	{
		errno = EDOM;
		return 0.0;
	}
	else if (p == 0)
	{
		errno = ERANGE;
		return -HUGE_VAL /* minus "infinity" */;
	}
	else if (p == 1)
	{
		errno = ERANGE;
		return HUGE_VAL /* "infinity" */;
	}
	else if (p < LOW)
	{
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else if (p > HIGH)
	{
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else
	{
		/* Rational approximation for central region */
    		q = p - 0.5;
    		r = q*q;
		return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
}

double normcdf(double xx, double mu, double sigma)
{
	int sn;
	double R1, R2, y, y2, y3, y4, y5, y6, y7, x;
	double erf, erfc, z, z2, z3, z4;
	double phi;

  x = (xx-mu)/sigma;
	if (x < -UPPERLIMIT) return 0.0;
	if (x > UPPERLIMIT) return 1.0;

	y = x / SQRT2;
	if (y < 0) 
	{
		y = -y;
		sn = -1;
	}
	else  sn = 1;

	y2 = y * y;
	y4 = y2 * y2;
	y6 = y4 * y2;

	if (y < 0.46875) 
	{
		R1 = P10 + P11 * y2 + P12 * y4 + P13 * y6;
		R2 = Q10 + Q11 * y2 + Q12 * y4 + Q13 * y6;
		erf = y * R1 / R2;
		if (sn == 1) phi = 0.5 + 0.5*erf;
		else  phi = 0.5 - 0.5*erf;
	}
	else
	{
		if (y < 4.0)
		{
			y3 = y2 * y;
			y5 = y4 * y;
			y7 = y6 * y;
			R1 = P20 + P21 * y + P22 * y2 + P23 * y3 + P24 * y4 + P25 * y5 + P26 * y6 + P27 * y7;
			R2 = Q20 + Q21 * y + Q22 * y2 + Q23 * y3 + Q24 * y4 + Q25 * y5 + Q26 * y6 + Q27 * y7;
			erfc = exp(-y2) * R1 / R2;
			if (sn == 1)  phi = 1.0 - 0.5*erfc;
			else  phi = 0.5*erfc;
		}
		else
		{
			z = y4;
			z2 = z * z;
			z3 = z2 * z;
			z4 = z2 * z2;
			R1 = P30 + P31 * z + P32 * z2 + P33 * z3 + P34 * z4;
			R2 = Q30 + Q31 * z + Q32 * z2 + Q33 * z3 + Q34 * z4;
			erfc = (exp(-y2)/y) * (1.0 / SQRTPI + R1 / (R2 * y2));
			if (sn == 1) phi = 1.0 - 0.5*erfc;
			else  phi = 0.5*erfc;
		}
  }
	return phi;
}

void randperm(int m, int *p)
{
  int j;
  double *u;
  
  u = (double *) mxCalloc(m,sizeof(double));
  
  for (j=0; j<m; j++) u[j] = unirnd(0,1);
  insertionSort(u, m, p);
}

  