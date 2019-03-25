#include "set_mesh.h"


void find_npnenf3D(char *fname,int *np,int *ne,int *nf){
	FILE *fp;
	int i,j,nnp,nne,nnf,c,dummy2;
	double **a;
	int **b;
	char dummy1[128];
	if((fp = fopen(fname,"r")) == NULL){
	  printf("file not found.\n");
		return;
	}else{
		fscanf(fp,"%s",dummy1); fscanf(fp,"%d",&dummy2);
		fscanf(fp,"%s",dummy1); fscanf(fp,"%d",&dummy2);
		fscanf(fp,"%s",dummy1);
	}
	fscanf(fp,"%d",&nnp);
	a=dmatrix(1,nnp,1,4);
	for(i=1;i<=nnp;i++) fscanf(fp,"%lf %lf %lf %lf",&a[i][1],&a[i][2],&a[i][3],&a[i][4]);
	free_dmatrix(a,1,nnp,1,4);
	fscanf(fp,"%s",dummy1); fscanf(fp,"%d",&nne); b = imatrix(1,nne,1,5);
	for(i=1;i<=nne;i++) fscanf(fp,"%d %d %d %d %d ",&b[i][1],&b[i][2],&b[i][3],&b[i][4],&b[i][5]);
	free_imatrix(b,1,nne,1,5);
	fscanf(fp,"%s",dummy1); fscanf(fp,"%d",&nnf); b=imatrix(1,nnf,1,4);
	for(i=1;i<=nnf;i++) fscanf(fp,"%d %d %d %d",&b[i][1],&b[i][2],&b[i][3],&b[i][4]);
	free_imatrix(b,1,nnf,1,4);
	*np=nnp;*ne=nne;*nf=nnf;
	return;
}

void find_npnenf2D(char *fname,int *np,int *ne,int *nf){
	FILE *fp;
	int _np,_ne,_nf;

	if((fp = fopen(fname,"r")) == NULL){
    printf("file not found.\n");
    return;
  }else{ fscanf(fp,"%d",&_np); fscanf(fp,"%d",&_ne); fscanf(fp,"%d",&_nf); }

	*np=_np;*ne=_ne,*nf=_nf;
	return;
}

void output_mesh3D(double **npxy,int **elnp,int **fanp,double *c,int _i,int np,int ne,int nf){
	char dname[1024],fname1[1024],fname2[1024],fname3[1024],fname4[1024];
	FILE *fin;
	int i,j,type;
	start("Output_mesh");
	mkdir("SOL/",0755);
	mkdir("SOL/3D/",0755);
	sprintf(dname,"SOL/3D/");
	mkdir(dname,0755);
	strcpy(fname1,dname); strcat(fname1,"npxy.dat");
	strcpy(fname2,dname); strcat(fname2,"elnp.dat");
	strcpy(fname3,dname); strcat(fname3,"fanp.dat");
	strcpy(fname4,dname); strcat(fname4,"sol.dat");

	if((fin=fopen(fname1,"w"))==NULL){
		printf("Can't open file: npxy.dat.\n");
		return;
	}
	fprintf(fin,"%d %d %d\n",np,3,1);
	for(i=1;i<=np;i++){
		for(j=1;j<=3;j++) fprintf(fin,"%lf ",npxy[i][j]);
		fprintf(fin,"\n");
	}
		fclose(fin);


	if((fin=fopen(fname2,"w"))==NULL){
		printf("Can't open file: elnp.dat.\n");
		return;
	}
	fprintf(fin,"%d %d %d\n",ne,4,0);
	for(i=1;i<=ne;i++){
		for(j=1;j<=4;j++) fprintf(fin,"%d ",elnp[i][j]);
		fprintf(fin,"\n");
	}
	fclose(fin);


	if((fin=fopen(fname3,"w"))==NULL){
		printf("Can't open file: fanp.dat.\n");
		return;
	}
	fprintf(fin,"%d %d %d\n",nf,4,0);
	for(i=1;i<=nf;i++){
		for(j=1;j<=4;j++) fprintf(fin,"%d ",fanp[i][j]);
		fprintf(fin,"\n");
	}
	fclose(fin);


	if((fin=fopen(fname4,"w"))==NULL){
		printf("Can't open file: npxy.dat.\n");
		return;
	}
	fprintf(fin,"%d %d %d\n",np,1,1);
	for(i=1;i<=np;i++)fprintf(fin,"c[i]\n");

	finish();
	return;
}

void input_mesh3D(char *fname,double **npxy,int **elnp,int **fanp,int np,int ne,int nf){
	FILE *fp;
	int i,j,nnp,nne,nnf,c,_int;
	char _char[128];
	start("Input_mesh");

	if((fp = fopen(fname,"r")) == NULL){
		printf("file not found.\n");
		return;
	}else{
		fscanf(fp,"%s",_char); fscanf(fp,"%d",&_int);
		fscanf(fp,"%s",_char); fscanf(fp,"%d",&_int);
		fscanf(fp,"%s",_char);
	}
	fscanf(fp,"%d",&_int);
	for(i=1;i<=np;i++) fscanf(fp,"%lf %lf %lf %lf",&npxy[i][1],&npxy[i][2],&npxy[i][3],&npxy[i][4]);
	fscanf(fp,"%s",_char); fscanf(fp,"%d",&_int);
	for(i=1;i<=ne;i++) fscanf(fp,"%d %d %d %d %d ",&elnp[i][1],&elnp[i][2],&elnp[i][3],&elnp[i][4],&elnp[i][5]);
	fscanf(fp,"%s",_char); fscanf(fp,"%d",&_int);
	for(i=1;i<=nf;i++) fscanf(fp,"%d %d %d %d",&fanp[i][1],&fanp[i][2],&fanp[i][3],&fanp[i][4]);
	fclose(fp);
	finish();
	return;
}

void input_mesh2D(char *fname,double **npxy,int **elnp,int **fanp,int np,int ne,int nf){
	FILE *fp;
	int i,dummy;
	start("input_mesh2D");
	if((fp = fopen(fname,"r")) == NULL){
		printf("file not found.\n");
		return;
	}else{
	fscanf(fp,"%d %d %d",&dummy,&dummy,&dummy);
	for(i=1;i<=np;i++)fscanf(fp,"%lf %lf %lf",&npxy[i][1],&npxy[i][2],&npxy[i][3]);
	for(i=1;i<=ne;i++)fscanf(fp,"%d %d %d %d",&elnp[i][1],&elnp[i][2],&elnp[i][3],&elnp[i][4]);
	for(i=1;i<=nf;i++)fscanf(fp,"%d %d %d ",&fanp[i][1],&fanp[i][2],&fanp[i][3]);
	}
	finish();
	return;
}

void output_sol2D(double **npxy,int **elnp,double *c,int np,int ne){
	FILE *fp;
	char fname[124],dname[124]= "SOL/2D/";
	int i,j,k,l;
	start("output_sol2D");
	mkdir("SOL/",0755);
	mkdir("SOL/2D/",0755);
	strcat(dname,"sol.dat");
	if((fp=fopen(dname,"w"))==NULL){
		printf("Can't open file: npxy.dat.\n");
		return;
	}
	for(i=1;i<=ne;i++){
		for(k=1;k<=3;k++){
			fprintf(fp,"%f %f %f\n",npxy[elnp[i][k]][1],npxy[elnp[i][k]][2],c[elnp[i][k]]);
		}
		fprintf(fp,"%f %f %f\n",npxy[elnp[i][1]][1],npxy[elnp[i][1]][2],c[elnp[i][1]]);
		fprintf(fp,"\n");
		fprintf(fp,"\n");
	}
	fclose(fp);
	finish();
	return;
}
