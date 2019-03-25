#include "functions.h"
#include <stdlib.h>
#include <stdio.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include <sys/stat.h>


/* argument of (x[1],x[2]) in [0,2\pai) or -1 (the origin) */
double arg( double *x )
{
  double eps=1.0e-10;
  double a;

  if ( x[1]*x[1]+x[2]*x[2] < eps )
    {
      return -1.0;
    }

  if ( x[1] >= fabs(x[2]) )
    {
      a = atan(x[2]/x[1]);
      return (a>=0.0)?a:(a+2*pi);
    }

  if ( x[2] >= fabs(x[1]) )
    {
      a = atan(-x[1]/x[2]);
      return a + 0.5*pi;
    }

  if ( x[1] <= - fabs(x[2]) )
    {
      a = atan(x[2]/x[1]);
      return a+pi;
    }

  if ( x[2] <= - fabs(x[1]) )
    {
      a = atan(-x[1]/x[2]);
      return a + 1.5*pi;
    }
  exit(0);
}


double *dvector(int i, int j)
{
  double *a,*q;
  int k;
  if ( (a=(double *)malloc( ((j-i+1)*sizeof(double))) ) == NULL )
    {
      printf("dvector: memory allocation is failed \n");
      exit(1);
    }
    q=a;
   for(k=i;k<=j;k++) q[k]=0.0;
  return(a-i);
}

void free_dvector(double *a, int i)
{
  free( (void *)(a + i) );
}

int *ivector(int i, int j)
{
  int *a,*q,k;

  if ( (a=(int *)malloc( ((j-i+1)*sizeof(int))) ) == NULL )
    {
      printf("ivector: memory allocation is failed \n");
      exit(1);
    }
    q=a; for(k=i;k<=j;k++) q[k]=0;
  return(a-i);
}

void free_ivector(int *a, int i)
{
  free( (void *)(a + i) );
}

double **dmatrix(int nr1, int nr2, int nl1, int nl2)
{
  int i, nrow, ncol,j;
  double **a;

  nrow = nr2 - nr1 + 1 ; /* number of row */
  ncol = nl2 - nl1 + 1 ; /* number of column */

  if ( ( a=(double **)malloc( nrow*sizeof(double *) ) ) == NULL )
  {
    printf("dmatrix: memory allocation is failed \n");
    exit(1);
  }
  a = a - nr1;

  for( i=nr1; i<=nr2; i++) a[i] = (double *)malloc(ncol*sizeof(double));
  for( i=nr1; i<=nr2; i++) a[i] = a[i]-nl1;

  for(i=nr1;i<=nr2;i++){
    for(j=nl1;j<=nl2;j++)a[i][j] = 0.0;
  }


  return(a);
}

void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2)
{
  int i;

  for ( i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
  free((void *)(a+nr1));
}

int **imatrix(int nr1, int nr2, int nl1, int nl2)
{
  int i, nrow, ncol,j;
  int **a;

  nrow = nr2 - nr1 + 1 ; /* number of row */
  ncol = nl2 - nl1 + 1 ; /* number of column */

  if ( ( a=(int **)malloc( nrow*sizeof(int *) ) ) == NULL )
    {
      printf("imatrix: memory allocation is failed \n");
      exit(1);
    }
  a = a - nr1;

  for( i=nr1; i<=nr2; i++) a[i] = (int *)malloc(ncol*sizeof(int));
  for( i=nr1; i<=nr2; i++) a[i] = a[i]-nl1;

  for(i=nr1;i<=nr2;i++){
    for(j=nl1;j<=nl2;j++)a[i][j] = 0;
  }


  return(a);
}

void free_imatrix(int **a, int nr1, int nr2, int nl1, int nl2)
{
  int i;

  for ( i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
  free((void *)(a+nr1));
}


void matrix_vector_product(double **a, double *b , double *c, int n)
{
  double wk;
  int i, j;

  for ( i = 1; i <= n; i++)
  {
    wk = 0.0;
    for ( j = 1; j <= n; j++ )
    {
      wk += a[i][j]*b[j];
    }
    c[i] = wk;
  }
}

double inner_product( int m, int n, double *a, double *b)
{
  int i;
  double s = 0.0;

  for( i = m; i <= n; i++) s += a[i]*b[i];

  return s ;
}

double vector_norm1( double *a, int m, int n )
{
  int i;
  double norm = 0.0;
  for ( i = m; i <= n; i++ )
  {
    norm += fabs(a[i]);
  }
  return norm;
}

void cg(double **a, double *b, double *x, int n, double EPS, int KMAX)
{
  double eps, *r, *p, *tmp, alpha, beta, work;
  int i, k=0;

  r = dvector(1,n); /* r[1...n] */
  p = dvector(1,n); /* p[1...n] */
  tmp = dvector(1,n); /* tmp[1...n] */

  matrix_vector_product( a, x, tmp, n );  /* tmp <- A x_0 */

  for( i = 1; i <= n; i++)
  {
    p[i] = b[i] - tmp[i] ; r[i] = p[i];
  }
  /* added by mk */
  eps = vector_norm1(r, 1, n);
  printf("initial residual=%e\n",eps);
  if ( eps < EPS )  goto OUTPUT;

  do
  {
    matrix_vector_product( a, p, tmp, n );  /* tmp <- A p_k */
    work = inner_product( 1, n, p, tmp); /* work <- (p,Ap_k) */
    alpha = inner_product( 1, n, p, r ) / work ;

    for( i = 1; i <= n; i++) x[i] = x[i] + alpha*p[i];
    for( i = 1; i <= n; i++) r[i] = r[i] - alpha*tmp[i];

    k++;
    eps = vector_norm1(r, 1, n);
    printf("%d-th residual =%e\n", k, eps);
    if ( eps < EPS )  goto OUTPUT;

    beta = - inner_product( 1, n, r, tmp) / work;
    for( i = 1; i <= n; i++) p[i] = r[i] + beta*p[i];
 }while( k < KMAX );

  OUTPUT:;

  free_dvector( r, 1 ); free_dvector( p, 1 ); free_dvector( tmp, 1 );
  if ( k == KMAX )
  {
    printf("CG did not converge within %d iterations\n", KMAX);
    exit(1);
  }
  else
  {
    printf("CG: %d iterations \n", k);
  }
}


void inv(double **a,double **inv_a,int n)
{
  double buf,tmp,c; //一時的なデータを蓄える
  int i,j,k; //カウンタ
//単位行列を作る
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++){
      inv_a[i][j]=(i==j)?1.0:0.0;
    }
  }

  //掃き出し法
  for(k=1;k<=n;k++){
    for(i=1;i<=n;i++){
      if(i != k){
        c = a[i][k]/a[k][k];
        for(j=1;j<=n;j++){
          a[i][j] -= c*a[k][j];
          inv_a[i][j] -= c*inv_a[k][j];
        }
      }
    }

  }

  for(i=1;i<=n;i++){
    c = 1.0/a[i][i];
    for(j=1;j<=n;j++){
      inv_a[i][j] = c*inv_a[i][j];
    }
  }

}

void inv4(double **a,double **inv_a,double *D)
{

  double det,buf;
  int n=4;  //配列の次数
  int i,j,k;


  det  = a[1][1]*a[2][2]*a[3][3]*a[4][4] + a[1][1]*a[2][3]*a[3][4]*a[4][2] + a[1][1]*a[2][4]*a[3][2]*a[4][3];
  det += a[1][2]*a[2][1]*a[3][4]*a[4][3] + a[1][2]*a[2][3]*a[3][1]*a[4][4] + a[1][2]*a[2][4]*a[3][3]*a[4][1];
  det += a[1][3]*a[2][1]*a[3][2]*a[4][4] + a[1][3]*a[2][2]*a[3][4]*a[4][1] + a[1][3]*a[2][4]*a[3][1]*a[4][2];
  det += a[1][4]*a[2][1]*a[3][3]*a[4][2] + a[1][4]*a[2][2]*a[3][1]*a[4][3] + a[1][4]*a[2][3]*a[3][2]*a[4][1];
  det -= a[1][1]*a[2][2]*a[3][4]*a[4][3] + a[1][1]*a[2][3]*a[3][2]*a[4][4] + a[1][1]*a[2][4]*a[3][3]*a[4][2];
  det -= a[1][2]*a[2][1]*a[3][3]*a[4][4] + a[1][2]*a[2][3]*a[3][4]*a[4][1] + a[1][2]*a[2][4]*a[3][1]*a[4][3];
  det -= a[1][3]*a[2][1]*a[3][4]*a[4][2] + a[1][3]*a[2][2]*a[3][1]*a[4][4] + a[1][3]*a[2][4]*a[3][2]*a[4][1];
  det -= a[1][4]*a[2][1]*a[3][2]*a[4][3] + a[1][4]*a[2][2]*a[3][3]*a[4][1] + a[1][4]*a[2][3]*a[3][1]*a[4][2];
  *D = det;

 inv_a[1][1] = a[2][2]*a[3][3]*a[4][4] + a[2][3]*a[3][4]*a[4][2] + a[2][4]*a[3][2]*a[4][3] - a[2][2]*a[3][4]*a[4][3] - a[2][3]*a[3][2]*a[4][4] - a[2][4]*a[3][3]*a[4][2];
 inv_a[1][2] = a[1][2]*a[3][4]*a[4][3] + a[1][3]*a[3][2]*a[4][4] + a[1][4]*a[3][3]*a[4][2] - a[1][2]*a[3][3]*a[4][4] - a[1][3]*a[3][4]*a[4][2] - a[1][4]*a[3][2]*a[4][3];
 inv_a[1][3] = a[1][2]*a[2][3]*a[4][4] + a[1][3]*a[2][4]*a[4][2] + a[1][4]*a[2][2]*a[4][3] - a[1][2]*a[2][4]*a[4][3] - a[1][3]*a[2][2]*a[4][4] - a[1][4]*a[2][3]*a[4][2];
 inv_a[1][4] = a[1][2]*a[2][4]*a[3][3] + a[1][3]*a[2][2]*a[3][4] + a[1][4]*a[2][3]*a[3][2] - a[1][2]*a[2][3]*a[3][4] - a[1][3]*a[2][4]*a[3][2] - a[1][4]*a[2][2]*a[3][3];
 inv_a[2][1] = a[2][1]*a[3][4]*a[4][3] + a[2][3]*a[3][1]*a[4][4] + a[2][4]*a[3][3]*a[4][1] - a[2][1]*a[3][3]*a[4][4] - a[2][3]*a[3][4]*a[4][1] - a[2][4]*a[3][1]*a[4][3];
 inv_a[2][2] = a[1][1]*a[3][3]*a[4][4] + a[1][3]*a[3][4]*a[4][1] + a[1][4]*a[3][1]*a[4][3] - a[1][1]*a[3][4]*a[4][3] - a[1][3]*a[3][1]*a[4][4] - a[1][4]*a[3][3]*a[4][1];
 inv_a[2][3] = a[1][1]*a[2][4]*a[4][3] + a[1][3]*a[2][1]*a[4][4] + a[1][4]*a[2][3]*a[4][1] - a[1][1]*a[2][3]*a[4][4] - a[1][3]*a[2][4]*a[4][1] - a[1][4]*a[2][1]*a[4][3];
 inv_a[2][4] = a[1][1]*a[2][3]*a[3][4] + a[1][3]*a[2][4]*a[3][1] + a[1][4]*a[2][1]*a[3][3] - a[1][1]*a[2][4]*a[3][3] - a[1][3]*a[2][1]*a[3][4] - a[1][4]*a[2][3]*a[3][1];
 inv_a[3][1] = a[2][1]*a[3][2]*a[4][4] + a[2][2]*a[3][4]*a[4][1] + a[2][4]*a[3][1]*a[4][2] - a[2][1]*a[3][4]*a[4][2] - a[2][2]*a[3][1]*a[4][4] - a[2][4]*a[3][2]*a[4][1];
 inv_a[3][2] = a[1][1]*a[3][4]*a[4][2] + a[1][2]*a[3][1]*a[4][4] + a[1][4]*a[3][2]*a[4][1] - a[1][1]*a[3][2]*a[4][4] - a[1][2]*a[3][4]*a[4][1] - a[1][4]*a[3][1]*a[4][2];
 inv_a[3][3] = a[1][1]*a[2][2]*a[4][4] + a[1][2]*a[2][4]*a[4][1] + a[1][4]*a[2][1]*a[4][2] - a[1][1]*a[2][4]*a[4][2] - a[1][2]*a[2][1]*a[4][4] - a[1][4]*a[2][2]*a[4][1];
 inv_a[3][4] = a[1][1]*a[2][4]*a[3][2] + a[1][2]*a[2][1]*a[3][4] + a[1][4]*a[2][2]*a[3][1] - a[1][1]*a[2][2]*a[3][4] - a[1][2]*a[2][4]*a[3][1] - a[1][4]*a[2][1]*a[3][2];
 inv_a[4][1] = a[2][1]*a[3][3]*a[4][2] + a[2][2]*a[3][1]*a[4][3] + a[2][3]*a[3][2]*a[4][1] - a[2][1]*a[3][2]*a[4][3] - a[2][2]*a[3][3]*a[4][1] - a[2][3]*a[3][1]*a[4][2];
 inv_a[4][2] = a[1][1]*a[3][2]*a[4][3] + a[1][2]*a[3][3]*a[4][1] + a[1][3]*a[3][1]*a[4][2] - a[1][1]*a[3][3]*a[4][2] - a[1][2]*a[3][1]*a[4][3] - a[1][3]*a[3][2]*a[4][1];
 inv_a[4][3] = a[1][1]*a[2][3]*a[4][2] + a[1][2]*a[2][1]*a[4][3] + a[1][3]*a[2][2]*a[4][1] - a[1][1]*a[2][2]*a[4][3] - a[1][2]*a[2][3]*a[4][1] - a[1][3]*a[2][1]*a[4][2];
 inv_a[4][4] = a[1][1]*a[2][2]*a[3][3] + a[1][2]*a[2][3]*a[3][1] + a[1][3]*a[2][1]*a[3][2] - a[1][1]*a[2][3]*a[3][2] - a[1][2]*a[2][1]*a[3][3] - a[1][3]*a[2][2]*a[3][1];

 for(i=1;i<=n;i++){
   for(j=1;j<=n;j++)inv_a[i][j] = inv_a[i][j]/det;
 }

 return ;

}

void inv3(double **a,double **inv_a,double *D){
  double det;
  double buf; //一時的なデータを蓄える
  int i,j,k; //カウンタ
  int n=3;  //配列の次数

  det  = a[1][1]*a[2][2]*a[3][3] + a[2][1]*a[3][2]*a[1][3] + a[3][1]*a[1][2]*a[2][3];
  det -= a[1][1]*a[3][2]*a[3][2] + a[2][1]*a[1][2]*a[3][3] + a[3][1]*a[2][2]*a[1][3];

  *D = det;

  inv_a[1][1] = a[2][2]*a[3][3] - a[2][3]*a[3][2]; inv_a[1][2] = a[1][3]*a[3][2] - a[1][2]*a[3][3]; inv_a[1][3] = a[1][2]*a[2][3] - a[1][3]*a[2][2];
  inv_a[2][1] = a[2][3]*a[3][1] - a[2][1]*a[3][3]; inv_a[2][2] = a[1][1]*a[3][3] - a[1][3]*a[3][1]; inv_a[2][3] = a[1][3]*a[2][1] - a[1][1]*a[2][3];
  inv_a[3][1] = a[2][1]*a[3][2] - a[2][2]*a[3][1]; inv_a[3][2] = a[1][2]*a[3][1] - a[1][1]*a[3][2]; inv_a[3][3] = a[1][1]*a[2][2] - a[1][2]*a[2][1];

  for(i=1;i<=3;i++){
    for(j=1;j<=3;j++)inv_a[i][j] = inv_a[i][j]/det;
  }

}

double det3(double **a)
{
  double x=0.0;

    x=a[1][1]*a[2][2]*a[3][3];
    x+=a[2][1]*a[3][2]*a[1][3];
    x+=a[3][1]*a[1][2]*a[2][3];
    x-=a[3][1]*a[2][2]*a[1][3];
    x-=a[2][1]*a[1][2]*a[3][3];
    x-=a[1][1]*a[3][2]*a[2][3];

  return x;
}

void start(char s[]){ printf("Start: %s ... ",s); fflush(stdout);}
void finish(void){ printf("Finish!\n"); fflush(stdout);}
void exit_func(char s[]){ printf("\n---\nExit: %s.\n---\n",s); fflush(stdout); exit(1);}

void init_imatrix(int **a,int i1,int i2,int j1,int j2){ int i,j; for(i=i1;i<=i2;i++){ for(j=j1;j<=j2;j++)a[i][j] = 0; } return;}
void init_dmatrix(double **a,int i1,int i2,int j1,int j2){int i,j; for(i=i1;i<=i2;i++){ for(j=j1;j<=j2;j++)a[i][j] = 0.0;} return;}

void   printdv(double *a,int i1,int i2){
  int i;
  for(i = i1;i<=i2;i++)printf("%f \n",a[i]);
}
void   printiv(int *a,int i1,int i2){
  int i;
  for(i = i1;i<=i2;i++)printf("%d \n",a[i]);
}
void   printdm(double **a,int i1,int i2,int j1,int j2){
  int i,j;
  for(i = i1;i<=i2;i++){
    for(j=j1;j<=j2;j++)printf("%f ",a[i][j]);
    puts("");
  }
}
void   printim(int **a,int i1,int i2,int j1,int j2){
  int i,j;
  for(i = i1;i<=i2;i++){
    for(j=j1;j<=j2;j++)printf("%d ",a[i][j]);
    puts("");
  }
}
