# Cfunc

c言語で使う基本的なヘッダファイルおよび関数を定義しています。

## ヘッダファイル一覧

- stdio.h
- stdlib.h
- string.h
- math.h
- time.h
- sys/stat.h

## define

|||
|---|---|
|pi  |円周率　3.14159265358979  |
|Max(a,b)  |a,bの最大値  |
|Min(a,b) | a,bの最小値|
|sgn(a)|aの符号|

## function

- 配列生成関連

|||
|---|---|
|dvector(int i,int j)|番号iからjの実数配列を生成。double *aのように定義|
|ivector(int i,int j)|番号iからjの整数配列を生成。int *aのように定義|
|dmatrix(int i1,int i2,int j1,int j2)|i1からi2列, j1からj2列の2次元実数配列を生成。double **a|
|imatrix(int i1,int i2,int j1,int j2)|i1からi2列, j1からj2列の2次元整数配列を生成。int **a|
|free_dvector(double *a, int i)|メモリ解放|
|free_ivector(int *a, int i);|　　|
|free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);|　　|
|free_imatrix(int **a, int nr1, int nr2, int nl1, int nl2);|　　|

- 配列計算

|||
|---|---|
|matrix_vector_product( double **a, double *b, double *c, int n )|行列aとベクトルbの積をcに入れる。n次元|
|vector_norm1( double *a, int m, int n )|mからnまでのベクトルaのユークリッドノルムを返す|
|inner_product( int m, int n, double *a, double *b)|mからnまでのベクトルa,bの内積を返す |
|det3( double **a )|3次元正方行列aの行列式|
|inv(double \*\*a, double **inv_a, int n)|n次元正方行列の逆行列をinv_aに返す|
| cg( double \*\*a, double *b, double *x, int n, double EPS, int KMAX )|連立1次方程式ax=bの解をcg法で解く。終了条件EPS　ステップ数KMAX|
