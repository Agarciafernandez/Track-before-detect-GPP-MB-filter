#include "mex.h"
#include "matrix.h"
#include "matrixutils.h"
#include <math.h>
#define SWAP(a,b) {temp=(a); (a)=(b); (b)=temp;}

double det(double *A, int n)
{
  double *L, *U, d;
  int i, k, j;
  
  L = (double *) mxCalloc(n*n,sizeof(double));
  U = (double *) mxCalloc(n*n,sizeof(double));
  for (i=0; i<n; i++) L[i+i*n] = 1;

  d = 1;
  for (k=0; k<n; k++)
  {
    for (j=k; j<n; j++)
    {
      U[k+j*n] = A[k+j*n];
      for (i=0; i<k; i++) U[k+j*n] -= L[k+i*n]*U[i+j*n];
    }
    d *= U[k+k*n];
    for (i=k+1; i<n; i++)
    {
      L[i+k*n] = A[i+k*n];
      for (j=0; j<k; j++) L[i+k*n] -= L[i+j*n]*U[j+k*n];
      L[i+k*n] /=U[k+k*n];
    }
  }
  
  return d;
}

          
double lngamma(double x)
{
  double z;
  
  x += 6;
  z = 1/(x*x);
  z = (((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z+0.083333333333333)/x;
  z = (x-0.5)*log(x)-x+0.918938533204673+z-log(x-1)-log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
  
  return z;
}

void starprod(double *A, int m, int n, double *B, int mp, int nq, double *C)
{
  int i, j, k, l, p, q;
  
  p = mp/m;
  q = (int) ((double) nq)/((double) q);

  for (k=0; k<p; k++)
  {
    for (l=0; l<q; l++)
    {
      C[k+l*p] = 0;
      for (i=0; i<m; i++)
      {
        for (j=0; j<n; j++) C[k+l*p] += A[i+j*m]*B[k+i*p+(l+j*q)*mp];
      }
    }
  }
  
}

int round2(double x)
{
  
  return floor(x+0.5);
  
}

double bessi0(double x)
{
  double ax,ans;
  double y, y2, y3, y4, e; 
  
  if ((ax=fabs(x)) < 3.75) // Polynomial fit.
  { 
    y=x/3.75;
    y*=y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  }
  else 
  {
    y=3.75/ax;
    e = exp(ax)/sqrt(ax);
    y2 = y*y;
    y3 = y2*y;
    y4 = y3*y;
    
    ans=e*(0.39894228+y*0.1328592e-1+y2*0.225319e-2-y3*0.157565e-2);
    ans += e*y4*(0.916281e-2-y*0.2057706e-1+y2*0.2635537e-1-y3*0.1647633e-1+y4*0.392377e-2);
  }
  
  return ans;
  
}

int max2int(int a, int b)
{
  if (b>a) return b;
  else return a;
}

int min2int(int a, int b)
{
  if (b>a) return a;
  else return b;
}

double max2(double a, double b)
{
  if (b>a) return b;
  else return a;
}

double min2(double a, double b)
{
  if (b>a) return a;
  else return b;
}

// double min(double *a, int n, int *pos)
// {
//   int i;
//   double mina;
//   
//   mina = a[0];
//   pos[0] = 0;
//   for (i=1; i<n; i++)
//   {
//     if (a[i]<mina)
//     {
//       mina = a[i];
//       pos[0] = i;
//     }
//   }
// 
//   return mina;
// }
// 
// double max(double *a, int n, int *pos)
// {
//   int i;
//   double maxa;
//   
//   maxa = a[0];
//   pos[0] = 0;
//   for (i=1; i<n; i++)
//   {
//     if (a[i]>maxa)
//     {
//       maxa = a[i];
//       pos[0] = i;
//     }
//   }
// 
//   return maxa;
// }

int maxint(int *a, int n)
{
  int maxa, i;
  
  maxa = a[0];
  for (i=1; i<n; i++)
  {
    if (a[i]>maxa) maxa = a[i];
  }
  
  return maxa;
}
  
int setdiff(int *a, int la, int *b, int lb, int *c)
{
  int c1, c2, isin, cnt;
  
  cnt = 0;
  for (c1=0; c1<la; c1++)
  {
    isin = 0;
    c2=0;
    while (c2<lb & ~isin) isin = (a[c1]==b[c2++]);
    if (!isin) c[cnt++] = a[c1];
  }
  
  return cnt;
}

int factorial (int n)
{
  int i, p;
  
  p = 1;
  for (i=2; i<=n; i++) p *= i;
  
  return p;
}

int binvecfind(int d, int J, int *b)
{
  int i, cnt;
  
  cnt = 0;
  for (i=0; i<J; i++)
  {
    if (modint(floor(d/pow(2,J-i-1)),2))
    {
      b[cnt] = i;
      cnt++;
    }
  }
  
  return cnt;
}


void binvec(int d, int J, int *b)
{
  int i;
  
  for (i=0; i<J; i++) b[i] = modint(floor(d/pow(2,J-i-1)),2);
}

void modkvec(int d, int J, int k, int *b)
{
  int i;
  
  for (i=0; i<J; i++) b[i] = modint(floor(d/pow(k,J-i-1)),k);
}

void modvec(int d, int J, int *k, int *b)
{
  int i, *cprod;
  
  cprod = (int *) mxCalloc(J,sizeof(int));
  
  if (J>0)
  {
    cprod[J-1] = 1;
    for (i=J-2; i>=0; i--) cprod[i] = cprod[i+1]*k[i+1];
  
    for (i=0; i<J; i++) b[i] = modint(floor(d/cprod[i]),k[i]);
  }

  mxFree(cprod);
}

int modint(int x, int y)
{
  return (x-floor(x/y)*y);
}

void kronint(int *X, int n1, int m1, int *Y, int n2, int m2, int *Z)
{
  int i, j, k, l;
  
  for (i=0; i<n1; i++)
  {
    for (j=0; j<m1; j++)
    {
      for (k=0; k<n2; k++)
      {
        for (l=0; l<m2; l++) Z[k+i*n2+(j*m2+l)*n1*n2] = X[i+j*n1]*Y[k+l*n2];
      }
    }
  }
}

void kron(double *X, int n1, int m1, double *Y, int n2, int m2, double *Z)
{
  int i, j, k, l;
  
  for (i=0; i<n1; i++)
  {
    for (j=0; j<m1; j++)
    {
      for (k=0; k<n2; k++)
      {
        for (l=0; l<m2; l++) Z[k+i*n2+(j*m2+l)*n1*n2] = X[i+j*n1]*Y[k+l*n2];
      }
    }
  }
}

int comb(int n, int k)
{
  int i, nCk;
  
  nCk = 1;
  if ((2*k)>n)
  {
    for (i=n; i>k; i--) nCk *= i;
    nCk /= factorial(n-k);
  }
  else
  {
    for (i=n; i>(n-k); i--) nCk *= i;
    nCk /= factorial(k);
  }
  
  
  return nCk;
}

void chol(double *X, int n, double *L)
/* Computes the Cholesky decomposition of the matrix X of size n.
The result is stored in L. Memory for L is assumed to be allocated.*/
{
  int i, j, k;
  double s;
  
  for (i=0; i<n; i++)
  {
    for (j=i; j<n; j++)
    {
      s = X[i*n+j];
      for (k=i-1; k>-1; k--)  s -= L[i*n+k]*L[j*n+k];
      if (i==j) L[j*n+i] = sqrt(s);
      else L[j*n+i] = s/L[i*n+i];
    }
  }
}

void cholsl(double *Ach, double *b, int n, double *x)
/* Solves Ax=b using Cholesky decomposition. 
     Ach is the upper triangular matrix satisfying Ach'*Ach=A*/
{
  int i, j;
  double *y;
  
  y = (double *) mxCalloc(n,sizeof(double));
  
  for (i=0; i<n; i++)
  {
    y[i] = b[i];
    for (j=0; j<i; j++) y[i] -= Ach[j+i*n]*y[j];
    y[i] = y[i]/Ach[i+i*n];
  }
  for (i=n-1; i>=0; i--)
  {
    x[i] = y[i];
    for (j=i+1; j<n; j++) x[i] -= Ach[i+j*n]*x[j];
    x[i] = x[i]/Ach[i+i*n];
  }
}

double inv(double *A, int n, double *Ainv)
{
    //Me da problemas si la matriz no es definida positiva->normal pues hay 
    //descomposicion de cholesqui
  double *b, *Ach, detA;
  int i,j;
  
  Ach = (double *) mxCalloc(n*n,sizeof(double));
  b = (double *) mxCalloc(n,sizeof(double));
  
  chol(A,n,Ach);
  detA = 1;
  for (i=0; i<n; i++)
  {
    detA *= Ach[i+i*n]*Ach[i+i*n];
    for(j=0; j<n; j++) b[j] = 0;
    b[i] = 1;
    cholsl(Ach,b,n,Ainv+i*n);
  }
  
  return detA;
}


void printmat(double *A,int rows, int cols)
/* Prints the elements of the matrix A which is of size rows x cols. */
{
    int c1,c2;
    
    for (c1=0; c1<rows; c1++)
    {
        for (c2=0; c2<cols; c2++) mexPrintf("%f ",A[c2*rows+c1]);
        mexPrintf("\n");
    }
}

void printintmat(int *A,int rows, int cols)
/* Prints the elements of the matrix A which is of size rows x cols. */
{
    int c1,c2;
    
    for (c1=0; c1<rows; c1++)
    {
        for (c2=0; c2<cols; c2++) mexPrintf("%d ",A[c2*rows+c1]);
        mexPrintf("\n");
    }
}

void matadd(double *A, double *B, int rows, int cols, double *C)
/* Matrix addition of A and B which are both of size rows x cols.
The result is stored in C. Memory for C is assumed to be allocated. */
{
    int c1,c2,ind;
    
    for (c1=0; c1<rows; c1++)
    {
        for (c2=0; c2<cols; c2++)
        {
            ind = c2*rows+c1;
            C[ind] = A[ind]+B[ind];
        }
    }
}

void matsub(double *A, double *B, int rows, int cols, double *C)
/* Matrix subtraction of B from A which are both of size rows x cols.
The result is stored in C. Memory for C is assumed to be allocated. */
{
    int c1,c2,ind;
    
    for (c1=0; c1<rows; c1++)
    {
        for (c2=0; c2<cols; c2++) 
        {
            ind = c2*rows+c1;
            C[ind] = A[ind]-B[ind];
         }
    }
}

void trans(double *A,int rows, int cols, double *B)
/* Matrix transpose of A which is of size rows x cols.
The result is stored in B. Memory for B is assumed to be allocated. */
{

    int c1,c2;
    
    for (c1=0; c1<rows; c1++)
    {
        for (c2=0; c2<cols; c2++) B[c1*cols+c2] = A[c2*rows+c1];
    }
}

void matmult(double *A,int na,int ma,double *B,int nb, int mb,double *C)
/* Matrix multiplication of A and B. A is na x ma and B is nb x mb. 
We must have ma=nb. There is no check for this.
The result is stored in C. Memory for C is assumed to be allocated. */
{
    int c1,c2,c3;
    
    for (c1=0; c1<na; c1++)
    {
        for (c2=0; c2<mb; c2++)
        {
            C[c2*na+c1] = 0;
            for (c3=0; c3<ma; c3++) C[c2*na+c1] = C[c2*na+c1]+A[c3*na+c1]*B[c2*nb+c3];
        }
    }
}

void compmatmult(double *Ar,double *Ai,int na,int ma,double *Br,double *Bi,int nb, int mb,double *Cr, double *Ci)
/* Matrix multiplication of A and B. A is na x ma and B is nb x mb. 
We must have ma=nb. There is no check for this.
The result is stored in C. Memory for C is assumed to be allocated. */
{
    int c1,c2;
    double *tmp1, *tmp2;
    
    tmp1 = (double *) mxCalloc(na*mb,sizeof(double));
    tmp2 = (double *) mxCalloc(na*mb,sizeof(double));
    
    matmult(Ar,na,ma,Br,nb,mb,Cr);
    matmult(Ai,na,ma,Bi,nb,mb,tmp1);
    matmult(Ar,na,ma,Bi,nb,mb,Ci);
    matmult(Ai,na,ma,Br,nb,mb,tmp2);

    for (c1=0; c1<na; c1++)
    {
      for (c2=0; c2<mb; c2++)
      {
        Cr[c1+c2*na] -= tmp1[c1+c2*na];
        Ci[c1+c2*na] += tmp2[c1+c2*na];
      }
    }
}

void invgj(double *A, int n, double *invA)
{
  int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  float big,dum,pivinv,temp;
  
  indxc = (int *) mxCalloc(n,sizeof(int));
  indxr = (int *) mxCalloc(n,sizeof(int));
  ipiv = (int *) mxCalloc(n,sizeof(int));
  
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) invA[j+i*n] = A[j+i*n];
    
  for (j=0; j<n; j++) ipiv[j]=0;
  for (i=0; i<n; i++)
  {  
    big=0.0;   
    for (j=0; j<n; j++)   
      if (ipiv[j] != 1)   
        for (k=0; k<n; k++)
        {
          if (ipiv[k] == 0)
          {
            if (fabs(invA[j+k*n]) >= big)
            {
              big=fabs(invA[j+k*n]);
              irow=j;
              icol=k;
            }
          }
        }
    ++(ipiv[icol]);
    if (irow != icol)
      for (l=0; l<n; l++) SWAP(invA[irow+l*n],invA[icol+l*n]);
    indxr[i]=irow;   
    indxc[i]=icol;   

    pivinv=1.0/invA[icol+icol*n];
    invA[icol+icol*n]=1.0;
    for (l=0; l<n; l++) invA[icol+l*n] *= pivinv;
    for (ll=0; ll<n; ll++)   
      if (ll != icol) 
      {    
        dum = invA[ll+icol*n];
        invA[ll+icol*n] = 0.0;
        for (l=0; l<n; l++) invA[ll+l*n] -= invA[icol+l*n]*dum;
      }
  }
  for (l=n-1; l>=0; l--) 
  {
    if (indxr[l] != indxc[l])
      for (k=0; k<n; k++)  SWAP(invA[k+indxr[l]*n],invA[k+indxc[l]*n]);
  }    
}

int find(double *zeta,int m,int *ind)
{
    int count,c;
    
    count = 0;
    for (c=0; c<m; c++)
    {
        if (zeta[c]>0.5) 
        {
            ind[count] = c;
            count++;
        }
    }
    return count;
}

int find_thresh(double *zeta, int m, double thresh, int *ind)
{
    int count,c;
    
    count = 0;
    for (c=0; c<m; c++)
    {
        if (zeta[c]>thresh) 
        {
            ind[count] = c;
            count++;
        }
    }
    return count;
}

int findint(int *zeta,int m,int *ind)
{
  int count,c;
    
  count = 0;
  for (c=0; c<m; c++)
  {
    if (zeta[c]!=0) 
    {
      ind[count] = c;
      count++;
    }
  }
  return count;
}

void insertionSort(double *numbers, int array_size, int *indices)
/* Sorts in ascending order */
{
  int i, j, intindex;
  double index;

  for (i=0; i<array_size; i++) indices[i] = i;
  
  for (i=1; i < array_size; i++)
  {
    index = numbers[i];
    intindex = i;
    j = i;
    while ((j > 0) && (numbers[j-1] > index))
    {
      numbers[j] = numbers[j-1];
      indices[j] = indices[j-1];
      j = j - 1;
    }
    numbers[j] = index;
    indices[j] = intindex;
  }
}

void insertionSort2(double *numbers, int array_size, int *indices)
/* Sorts in descending order */
{
  int i, j, intindex;
  double index;

  for (i=0; i<array_size; i++) indices[i] = i;
  
  for (i=1; i < array_size; i++)
  {
    index = numbers[i];
    intindex = i;
    j = i;
    while ((j > 0) && (numbers[j-1] < index))
    {
      numbers[j] = numbers[j-1];
      indices[j] = indices[j-1];
      j = j - 1;
    }
    numbers[j] = index;
    indices[j] = intindex;
  }
}

void sortint(int *numbers, int array_size, int *indices)
/* Sorts in ascending order */
{
  int i, j, intindex;
  double index;

  for (i=0; i<array_size; i++) indices[i] = i;
  
  for (i=1; i < array_size; i++)
  {
    index = numbers[i];
    intindex = i;
    j = i;
    while ((j > 0) && (numbers[j-1] > index))
    {
      numbers[j] = numbers[j-1];
      indices[j] = indices[j-1];
      j = j - 1;
    }
    numbers[j] = index;
    indices[j] = intindex;
  }
}

void matmultint(int *A,int na,int ma,int *B,int nb, int mb,int *C)
{
    int c1,c2,c3;
    
    for (c1=0; c1<na; c1++)
    {
        for (c2=0; c2<mb; c2++)
        {
            C[c2*na+c1] = 0;
            for (c3=0; c3<ma; c3++) C[c2*na+c1] = C[c2*na+c1]+A[c3*na+c1]*B[c2*nb+c3];
        }
    }
}

void eye(int n, double *I)
{
  int i, j;
  
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++) I[i+j*n] = 0;
    I[i*n+i] = 1;
  }
}

void ones(int n, double *o)
{
  int i;
  
  for (i=0; i<n; i++) o[i] = 1;
}

void repmat(int *Y, int nrows, int ncols, int n1, int m1, int *Z)
{
  int i, j, k, l;
  
  for (i=0; i<n1; i++)
  {
    for (j=0; j<m1; j++)  
    {
      for (k=0; k<nrows; k++)
      {
        for (l=0; l<ncols; l++) Z[i*nrows+k+(j*ncols+l)*n1*nrows] = Y[k+l*nrows];
      }
    }
  }
}

double mod(double x,double y)
{
    double z;
    
    if (y != 0)
        z = x-floor(x/y)*y;
    else
        z = x;
    
    return z;
}

int isinvec(int z, int *v, int m)
{
  int i, ind;
  
  i = 0;
  ind = 0;
  while (i<m && !ind)
  {
    if (z==v[i]) ind=1;
    i++;
  }
  
  return ind;
}

int isin(int *zeta, int **V, int *numvalid, int j)
{
    int i,ind;
    
    ind = 0;
    for (i=0; i<numvalid[j]; i++)
        ind = ind | (zeta[j]==V[j][i]);
    
    return ind;
}

void diff(int *a, int *b, int numels)
{
  int i;
  
  for (i=1; i<numels; i++) b[i-1] = a[i]-a[i-1];
}

int isunique(int *zeta,int j, int r)
{
    int i,ind;
    
    ind = 1;
    for (i=0; i<j; i++)
        ind = ind & (zeta[j]!=zeta[i]);
        
    for (i=j+1; i<r; i++)
        ind = ind & (zeta[j]!=zeta[i]);
        
    return ind;
}

void nchoosek(int *v, int d, int m, int *C)
{
  int i, j, l, numadd, *Ctmp, indmax, ind, nrows;

  if (m==0) return;
  
  if (m==1)
  {
    for (i=0; i<d; i++) C[i] = v[i];
  }
  else
  {
    nrows = comb(d,m);
    Ctmp = (int *) mxCalloc((m-1)*comb(d-1,m-1),sizeof(int));
    indmax = 0;
    for (i=0; i<d-m+1; i++)
    {
      nchoosek(v+i+1,d-i-1,m-1,Ctmp);
      numadd = comb(d-i-1,m-1);
      for (j=0; j<numadd; j++)
      {
        ind = indmax+j;
        C[ind] = v[i];
        for (l=0; l<m-1; l++) C[ind+(l+1)*nrows] = Ctmp[j+l*numadd];
      }
      indmax += numadd;
    }
  }
}

void perms(int *v, int n, int *P)
{
  int i, j, l, numadd, *Ptmp, *vtmp, indmax;
  
  if (n==0) return;
  
  if (n==1)
    P[0] = v[0];
  else
  {
    numadd = factorial(n-1);
    Ptmp = (int *) mxCalloc((n-1)*numadd,sizeof(int));
    vtmp = (int *) mxCalloc(n-1,sizeof(int));
    indmax = 0;
    for (i=0; i<n; i++)
    {
      for (j=0; j<i; j++) vtmp[j] = v[j];
      for (j=i+1; j<n; j++) vtmp[j-1] = v[j];
      perms(vtmp,n-1,Ptmp);
      for (j=0; j<numadd; j++)
      {
        P[j+indmax] = v[i];
        for (l=0; l<n-1; l++) P[j+indmax+(l+1)*factorial(n)] = Ptmp[j+l*numadd];
      }
      indmax += numadd;
    }
  }
}

// double atan2(double y, double x)
// {
//   double z;
//   
//   z = atan(y/x);
//   if (x<0)
//   {
//     if (y>0) z += PI;
//     else z -= PI;
//   }
//   
//   return z;
// }

int unionint(int *A, int na, int *B, int nb, int *C)
{
  int nc, i, j, *dummy, notin;
  
  for (i=0; i<na; i++) C[i] = A[i];
  nc = na;
  
  for (j=0; j<nb; j++)
  {
    notin = 1;
    i = 0;
    while(notin && (i<nc))
    {
      if (B[j]==C[i++]) notin=0;
    }
    if (notin)  C[nc++] = B[j];
  }
  
  dummy = (int *) mxCalloc(nc,sizeof(int));
  
  sortint(C,nc,dummy);
  
  return nc;
}

double sign(double x)
{
  double y;
  
  y = x/fabs(x);
  
  return y;
}

double sec(double x)
{
  double y;
  
  y = 1/cos(x);
  
  return y;
}

int intersect(int *x, int n, int *y, int m, int *z)
{
  int i, j, matched, cnt;
  
  cnt = 0;
  for (i=0; i<n; i++)
  {
    matched = 0;
    j = 0;
    while (j<m && !matched)
    {
      if (x[i]==y[j])
      {
        z[cnt] = x[i];
        cnt++;
        matched = 1;
      }
      else j++;
    }
  }
  
  return cnt;
}
