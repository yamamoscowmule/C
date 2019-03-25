#ifndef _SET_CONDITION_H_
#define _SET_CONDITION_H_

void set_F(double *F,int np,int dim); // 外力
double g(double x,double y); // ディリクレ境界条件
void make_gg(double *gg,int **dp,double **npxy);
void make_gg3D(double *gg1,double *gg2,double *gg3,int **dp,double **npxy);


typedef struct{ // 弾性テンソルをlambda,muの関数として表現　C(la,mu)
	double la;
	double mu;
}bc_t ;

bc_t Inv(bc_t C); // 逆テンソル　C^-1(la,mu)　の計算
bc_t Prod(bc_t C1,bc_t C2); //テンソル積の計算　C(la1,mu1) * C(la2,mu2)
bc_t Sum(bc_t C1,bc_t C2);
bc_t Diff(bc_t C1,bc_t C2);
bc_t Times(double al,bc_t C);

int  ij2a(int i,int j,int *ia,int *ja); // 座標をcrs表記に
void matrix_vector_product_crs(double *A,int *ia,int *ja,double *b , double *c, int n); // crs表記の行列の積
void cgcrs(double *A,int *ia,int *ja,double *b, double *x, int n, double EPS, int KMAX); //crsでのcg

void make_iaja(int **npnp,int *ia,int *ja,int np,int dim);
void make_npnp(int **elnp,int **npnp,int np,int ne,int *na,int dim); //接点接点対応表を作る

void set_dbc(int **dp, int *ia, int *ja, double *A, double *b,double *g1,double *g2); //ディリクレ2d
void set_dbc3D(int **dp, int *ia, int *ja, double *A, double *b,double *gg1,double *gg2,double *gg3); //ディリクレ３d

#endif
