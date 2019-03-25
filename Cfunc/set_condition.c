#include "set_condition.h"

void set_F(double *F,int np,int dim){
	int i;

	if(dim == 2){
	for(i=1;i<=np;i++){
		F[i]      =  0.0;
		F[i+np]   =  -1.0;
	}
}else	if(dim == 3){
	for(i=1;i<=np;i++){
	F[i]      =  0.0;
	F[i+np]   =  0.0;
	F[i+np+np]= -1.0;
	}
}else{
printf("Error : dim");
}
	return;
}

double g(double x,double y){
  return 0.0;
}

void make_gg(double *gg,int **dp,double **npxy){
	int i;
	for(i=1;i<=dp[0][0];i++){
		gg[i] = g(npxy[ dp[i][0] ][1],npxy[ dp[i][0] ][2]);
	}
}

void make_gg3D(double *gg1,double *gg2,double *gg3,int **dp,double **npxy){
	int i;

  for(i=1;i<=dp[0][0];i++){
    gg1[i] = 0.0;
    gg2[i] = 0.0;
    gg3[i] = 0.0;

    /*
    gg1[i] = g( npxy[ dp[i] ][1] , npxy[ dp[i] ][2] );
    gg2[i] = g( npxy[ dp[i] ][1] , npxy[ dp[i] ][2] );
    gg3[i] = g( npxy[ dp[i] ][1] , npxy[ dp[i] ][2] );
    */
  }

}

bc_t Inv(bc_t C){
	bc_t inv;
	int d = 2.0; // d=次元
	inv.la = -C.la/(2.0*C.mu*(d*C.la + 2*C.mu));
	inv.mu = 0.25/C.mu;
	return inv;
}

bc_t Prod(bc_t C1,bc_t C2){
	bc_t Prod;
	int d = 2.0; // d=次元
	Prod.la = d*C1.la*C2.la + 2.0*(C1.la*C2.mu + C2.la*C1.mu);
	Prod.mu = 2.0*C1.mu*C2.mu;
	return Prod;
}

bc_t Sum(bc_t C1,bc_t C2){
	bc_t sum;
	sum.la = C1.la + C2.la;
	sum.mu = C1.mu + C2.mu;
	return sum;
}

bc_t Diff(bc_t C1,bc_t C2){
	bc_t diff;
	diff.la = C1.la - C2.la;
	diff.mu = C1.mu - C2.mu;
	return diff;
}

bc_t Times(double al,bc_t C){
	bc_t times;
	times.la = al*C.la;
	times.mu = al*C.mu;
	return times;
}


int  ij2a(int i,int j,int *ia,int *ja){
	int k,z;
	z=ia[i+1]-ia[i];
		for(k=0;k<z;k++){
		if(ja[ia[i] + k] == j)   return ia[i]+k;
	}
		return -1;
}

void matrix_vector_product_crs(double *A,int *ia,int *ja,double *b , double *c, int n){
	double wk;
	int i, j;

	for ( i = 1; i <= n; i++)
	{
		wk = 0.0;
		for ( j = ia[i]; j < ia[i+1]; j++ )
		{
			wk += A[j] * b[ja[j]];
			//      wk += A[ ij2a(i,j,ia,ja) ]*b[j];
		}
		c[i] = wk;
	}
}

void cgcrs(double *A,int *ia,int *ja,double *b, double *x, int n, double EPS, int KMAX){
	double eps, *r, *p, *tmp, alpha, beta, work;
	int i, k=0;

	r = dvector(1,n); /* r[1...n] */
	p = dvector(1,n); /* p[1...n] */
	tmp = dvector(1,n); /* tmp[1...n] */
	for(i=1;i<=n;i++){
		r[i]=0.0;p[i]=0.0;tmp[i]=0.0;x[i]=0.0;
	}

	matrix_vector_product_crs( A, ia, ja, x, tmp, n );  /* tmp <- A x_0 */

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
		matrix_vector_product_crs( A, ia, ja, p, tmp, n );  /* tmp <- A p_k */
		work = inner_product( 1, n, p, tmp); /* work <- (p,Ap_k) */
		alpha = inner_product( 1, n, p, r ) / work ;

		for( i = 1; i <= n; i++) x[i] = x[i] + alpha*p[i];
		for( i = 1; i <= n; i++) r[i] = r[i] - alpha*tmp[i];

		k++;
		eps = vector_norm1(r, 1, n);
		if((k % 100) == 0) printf("%d-th residual =%e\n", k, eps);
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

void make_iaja(int **npnp,int *ia,int *ja,int np,int dim){
  int i,j=1,k,l=1,nnpnp;

	if(dim == 2){
  for(i=1;i<=np;i++){
    k=0;
    while(npnp[i][k+1] != 0) k++;
    for(j=1;j<=k;j++){
      ja[ ia[i] + j - 1 ] = npnp[i][j];
      ja[ ia[i] + j - 1 + k ] = npnp[i][j] + np;
    }
    ia[i+1] = ia[i] + 2*k;
  }

  for(i=np+1;i<=2*np;i++){
    k=0;
    while(npnp[i-np][k+1] != 0) k++;

    for(j=1;j<=k;j++){
      ja[ ia[i] + j - 1 ] = npnp[i-np][j];
      ja[ ia[i] + j - 1 + k ] = npnp[i-np][j] + np;
    }
    ia[i+1] = ia[i] + 2*k;
  	}
	}else if(dim == 3){
	for(i=1;i<=np;i++){
		k=0;
		while(npnp[i][k+1] != 0) k++;
		for(j=1;j<=k;j++){
			ja[ ia[i] + j - 1 ] = npnp[i][j];
			ja[ ia[i] + j - 1 + k ] = npnp[i][j] + np;
			ja[ ia[i] + j - 1 + 2*k]= npnp[i][j] + 2*np;
		}
		ia[i+1] = ia[i] + 3*k;
	}

	for(i=np+1;i<=2*np;i++){
		k=0;
		while(npnp[i-np][k+1] != 0) k++;

		for(j=1;j<=k;j++){
			ja[ ia[i] + j - 1 ] = npnp[i-np][j];
			ja[ ia[i] + j - 1 + k ] = npnp[i-np][j] + np;
			ja[ ia[i] + j - 1 + 2*k]= npnp[i-np][j] + 2*np;
		}
		ia[i+1] = ia[i] + 3*k;
	}

	for(i=2*np+1;i<=3*np;i++){
		k=0;
		while(npnp[i-2*np][k+1] != 0) k++;
		for(j=1;j<=k;j++){
			ja[ ia[i] + j - 1 ] = npnp[i-2*np][j];
			ja[ ia[i] + j - 1 + k ] = npnp[i-2*np][j] + np;
			ja[ ia[i] + j - 1 + 2*k]= npnp[i-2*np][j] + 2*np;
		}
		ia[i+1] = ia[i] + 3*k;
	}

}else{
	printf("Error : dim");
}
}

void make_npnp(int **elnp,int **npnp,int np,int ne,int *na,int dim){
	int i,j,k,l,x=0,tmp;

	if(dim == 2){
		for(i=1;i<=ne;i++){
			for(j=1;j<=3;j++){
				for(k=1;k<=3;k++){
					for(l=1;l<=30;l++){
						if(npnp[ elnp[i][j] ][l] == elnp[i][k]) break;
						if(npnp[ elnp[i][j] ][l] == 0){
							npnp[ elnp[i][j] ][l] = elnp[i][k]; x++;break;
						}
					}
				}
			}
		}
	} else if(dim == 3){
	for(i=1;i<=ne;i++){
		for(j=1;j<=4;j++){
			for(k=1;k<=4;k++){
				for(l=1;l<=30;l++){
					if(npnp[ elnp[i][j] ][l] == elnp[i][k]) break;

					if(npnp[ elnp[i][j] ][l] == 0){
						npnp[ elnp[i][j] ][l] = elnp[i][k]; x++;break;
					}
				}
			}
		}
	}
	}	else{
		printf("Error: dim");
		return;
	}

	//sort//
	for(i=1;i<=np;i++){
		for(j=1;j<=29;j++){
			for(k=30;k>=j+1;k--){
				if(npnp[i][k-1] > npnp[i][k] && npnp[i][k] != 0){
					tmp = npnp[i][k];
					npnp[i][k] = npnp[i][k-1];
					npnp[i][k-1] = tmp;
				}
			}
		}
	}
	*na=x;
	return;
}

void set_dbc(int **dp, int *ia, int *ja, double *A, double *b,double *g1,double *g2){
  int ndp=dp[0][0],j,k,i0,ii,ij;
  int np = b[0];

  for(ii=1;ii<=ndp;ii++)
  {
    i0=dp[ii][0];
    for(k=ia[i0];k<ia[i0+1];k++){
      j=ja[k];
      if(j==i0){
	      ij=ij2a(j,i0+np,ia,ja); b[j]-=A[ij]*g2[ii];
	      continue;
      }
      if(j==i0+np){
	     ij=ij2a(j,i0,ia,ja); b[j]-=A[ij]*g1[ii];
	      continue;
      }
      ij=ij2a(j,i0,ia,ja); b[j]-=A[ij]*g1[ii];
      ij=ij2a(j,i0+np,ia,ja); b[j]-=A[ij]*g2[ii];
    }

  }

 for(ii=1;ii<=ndp;ii++){
   i0=dp[ii][0]; b[i0]=g1[ii]; b[i0+np]=g2[ii+np];

    for(k=ia[i0];k<ia[i0+1];k++){
      j=ja[k];
      if(j==i0){
	      A[k]=1.0;
	      ij=ij2a(j,i0+np,ia,ja);A[ij]=0.0;
	      ij=ij2a(i0+np,j,ia,ja);A[ij]=0.0;
      }
      else if(j == i0+np){
	     ij=ij2a(j,j,ia,ja); A[ij]=1.0;
	     A[k]=0.0;
	     ij=ij2a(j,i0,ia,ja); A[ij]=0.0;
      }else{
	      A[k]=0.0;
	      ij=ij2a(j,i0,ia,ja); A[ij]=0.0;
	      ij=ij2a(j,i0+np,ia,ja);A[ij]=0.0;
	      ij=ij2a(i0+np,j,ia,ja);A[ij]=0.0;
      }
    }
  }
  return;
}
void set_dbc3D(int **dp, int *ia, int *ja, double *A, double *b,double *gg1,double *gg2,double *gg3){
  int ndp=dp[0][0],j,k,i0,ii,ij,np=b[0];

  for(ii=1;ii<=ndp;ii++)
  {
    i0=dp[ii][0];
    for(k=ia[i0];k<ia[i0+1];k++)
    {
      j=ja[k];
      if(j==i0){
	ij=ij2a(j,i0+np,ia,ja); b[j]-=A[ij]*gg2[ii];
	ij=ij2a(j,i0+np+np,ia,ja); b[j]-=A[ij]*gg3[ii];
	continue;
      }
      if(j==i0+np){
	ij=ij2a(j,i0,ia,ja); b[j]-=A[ij]*gg1[ii];
	ij=ij2a(j,i0+np+np,ia,ja); b[j]-=A[ij]*gg3[ii];
	continue;
      }
      if(j==i0+np+np){
	ij=ij2a(j,i0,ia,ja); b[j]-=A[ij]*gg1[ii];
	ij=ij2a(j,i0+np,ia,ja); b[j]-=A[ij]*gg2[ii];
      }
     ij=ij2a(j,i0,ia,ja); b[j]-=A[ij]*gg1[ii];
     ij=ij2a(j,i0+np,ia,ja); b[j]-=A[ij]*gg2[ii];
     ij=ij2a(j,i0+np+np,ia,ja); b[j]-=A[ij]*gg3[ii];
    }
 }

 for(ii=1;ii<=ndp;ii++){
   i0=dp[ii][0]; b[i0]=gg1[i0]; b[i0+np]=gg2[i0+np];b[i0+np+np]=gg3[i0+np+np];
   for(k=ia[i0];k<ia[i0+1];k++){
      j=ja[k];
      if(j==i0){
	A[k]=1.0;
	ij=ij2a(j,i0+np,ia,ja);    A[ij]=0.0;
	ij=ij2a(i0+np,j,ia,ja);    A[ij]=0.0;
	ij=ij2a(j,i0+np+np,ia,ja); A[ij]=0.0;
	ij=ij2a(i0+np+np,j,ia,ja); A[ij]=0.0;
      }
      else if(j == i0+np){
	ij=ij2a(j,j,ia,ja); A[ij]=1.0;
	A[k]=0.0;
	ij=ij2a(j,i0,ia,ja); A[ij]=0.0;
      }
      else if(j == i0+np+np){
	ij=ij2a(j,j,ia,ja); A[ij]=1.0;
	A[k]=0.0;
	ij=ij2a(j,i0,ia,ja); A[ij]=0.0;
      }else{
	A[k]=0.0;
	ij=ij2a(j,i0,ia,ja); A[ij]=0.0;

	ij=ij2a(j,i0+np,ia,ja);A[ij]=0.0;
	ij=ij2a(j,i0+np+np,ia,ja);A[ij]=0.0;

	ij=ij2a(i0+np,j,ia,ja);A[ij]=0.0;
	ij=ij2a(i0+np+np,j,ia,ja);A[ij]=0.0;
      }
  }
 }

 return;
}
