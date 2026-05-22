#define OneUponSqrt2Pi .39894228040143267794
#define twopi 6.283195307179587
#define LnSqrt2Pi -0.9189385332046727417803296 
#define SQRT2 1.414213562373095049
#define SQRTPI 1.772453850905516027

#define UPPERLIMIT 30.0 

#define P10 242.66795523053175
#define P11 21.979261618294152
#define P12 6.9963834886191355
#define P13 -.035609843701815385
#define Q10 215.05887586986120
#define Q11 91.164905404514901
#define Q12 15.082797630407787
#define Q13 1.0

#define P20 300.4592610201616005
#define P21 451.9189537118729422
#define P22 339.3208167343436870
#define P23 152.9892850469404039
#define P24 43.16222722205673530
#define P25 7.211758250883093659
#define P26 .5641955174789739711
#define P27 -.0000001368648573827167067
#define Q20 300.4592609569832933
#define Q21 790.9509253278980272
#define Q22 931.3540948506096211
#define Q23 638.9802644656311665
#define Q24 277.5854447439876434
#define Q25 77.00015293522947295
#define Q26 12.78272731962942351
#define Q27 1.0

#define P30 -.00299610707703542174
#define P31 -.0494730910623250734
#define P32 -.226956593539686930
#define P33 -.278661308609647788
#define P34 -.0223192459734184686
#define Q30 .0106209230528467918
#define Q31 .191308926107829841
#define Q32 1.05167510706793207
#define Q33 1.98733201817135256
#define Q34 1.0

#define LOW 0.02425
#define HIGH 0.97575

double normpdf(double x, double mu, double sigma);
double exppdf(double x, double lambda);
double norminv(double p);
double normcdf(double xx, double mu, double sigma);
double genonerandn();
double unirnd(double a, double b);
double trnd(double v);
double vmrnd(double mu, double kappa);
double vmpdf(double x, double mu, double kappa);
double normrnd(double mu, double sigma);
int poissrnd(double lambda);
int poissrnd2(double lambda);
double poisspdf(int x, double lambda);
double gamrnd(double alpha, double beta);
void randperm(int m, int *p);
double binormprob(double *lbnds, double *ubnds, double *Sigmasqrt, int nsamp);
double lnnormrectprob(double a, double b, int nreal);
double lnbinormprob(double *lbnds, double *ubnds, double *G, int nsamp);
double lntrectprob(double a, double b, double v, int nreal);
double lntrectprob2(double a, double b, double v);
