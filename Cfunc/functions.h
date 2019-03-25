#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define pi 3.14159265358979
#define Max(a,b) (((a)>(b))?(a):(b))
#define Min(a,b) (((a)<(b))?(a):(b))
#define Sgn(a) (((a)==0.0)?0.0:((a)/(fabs(a))))

double arg( double *x );
double *dvector(int i,int j);
int    *ivector(int i, int j);
double **dmatrix(int nr1, int nr2, int nl1, int nl2);
int    **imatrix(int nr1, int nr2, int nl1, int nl2);
void   free_dvector(double *a, int i);
void   free_ivector(int *a, int i);
void   free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
void   free_imatrix(int **a, int nr1, int nr2, int nl1, int nl2);
void   matrix_vector_product( double **a, double *b, double *c, int n );
double vector_norm1( double *a, int m, int n );
double inner_product( int m, int n, double *a, double *b);

void   cg( double **a, double *b, double *x, int n, double EPS, int KMAX );
double det3( double **a );
void   inv(double **a,double **inv_a, int n);
void   inv4(double **a,double **inv_a,double *D);
void   start(char s[]);
void   finish(void);
void   exit_func(char s[]);
void   printdv(double *a,int i1,int i2);
void   printiv(int *a,int i1,int i2);
void   printdm(double **a,int i1,int i2,int j1,int j2);
void   printim(int **a,int i1,int i2,int j1,int j2);

#endif
