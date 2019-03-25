#ifndef _SET_MESH_H_
#define _SET_MESH_H_


// find_npnenf
// FreeFem++で作成したメッシュファイル(Th.msh)から要素数、接点数、境界接点数を求める
void find_npnenf3D(char *fname,int *np,int *ne,int *nf);
void find_npnenf2D(char *fname,int *np,int *ne,int *nf);

// Output_mesh
// npxy.dat、elnp.dat、fanp.dat及び解sol.datを、フォルダ名「ステップ数」に格納
// paraviewで読み込めるファイルを作る目的
// npxy.datは動いた後の座標
// 3Dのみ

void output_mesh3D(double **npxy,int **elnp,int **fanp,double *c,int _i,int np,int ne,int nf);

// input_mesh
// FreeFemで作成したメッシュファイル(Th.msh)から行列npxy,elnp,fanpを作成
void input_mesh3D(char *fname,double **npxy,int **elnp,int **fanp,int np,int ne,int nf);
void input_mesh2D(char *fname,double **npxy,int **elnp,int **fanp,int np,int ne,int nf);

// output_sol
// 解をSOLフォルダに格納
// 3Dも作る予定
void output_sol2D(double **npxy,int **elnp,double *c,int np,int ne,char *outputname);

#endif
