#ifndef DATA_CLEAN_H
#define DATA_CLEAN_H

#include "stdio.h"
#include "stdlib.h"
#include "DataClean.h"
#include "datatype.h"
#include "parameter.h"

double Em[3] = {1.0,0.0,0.0};
double w[4] = {-0.5625, 0.52083333333, 0.52083333333, 0.52083333333};
double c1[4] = {0.33333333333, 0.6, 0.2 ,0.2};
double c2[4] = {0.33333333333, 0.2, 0.6, 0.2};
double c3[4] = {0.33333333334, 0.2, 0.2, 0.6};
double kik0[3] = {0.0,0.0, -k0};
int maxedge = 0, maxnode = 0, maxpatch =0;
int *ipatch = NULL, *iedge = NULL;
double *xyznode = NULL  ;

cuDoubleComplex *ImpedanceMatrix = NULL, *Excitation = NULL; 


void ReadData_all_surface(void){

	int temp, i, j;
	//int *ipatchp = NULL;
	//1，读取all_surface文件
	FILE * fp = fopen("Input\\all_surface.txt", "r");

	fscanf(fp,"%d",&maxpatch);
	fscanf(fp,"%d",&temp);
	fscanf(fp,"%d",&maxpatch);
	fscanf(fp,"%d",&temp);


	ipatch = (int *)malloc(sizeof(int) * maxpatch * 3);
	//printf("ss = %d",maxpatch); *ipatch + i * 3 + j

	for ( i = 0; i < maxpatch; ++ i){
		fscanf(fp,"%d",&temp);
		for ( j = 0; j < 3; ++ j){
			fscanf(fp,"%d", ipatch + i * 3 + j );
		}
		fscanf(fp,"%d",&temp);
	}

	fclose(fp);

}


void ReadData_xyz_node(void){

	int temp, i, j;
	//int *ipatchp = NULL;
	//1，读取all_surface文件
	FILE * fp = fopen("Input\\total_Nodes.txt", "r");

	fscanf(fp,"%d", &maxnode);
	fscanf(fp,"%d", &temp);
	fscanf(fp,"%d", &maxnode);
	fscanf(fp,"%d", &temp);


	xyznode = (double *)malloc(sizeof(double) * maxnode * 3);
	//printf("ss = %d",maxpatch); *ipatch + i * 3 + j

	for ( i = 0; i < maxnode; ++ i){
		//fscanf(fp,"%d",&temp);
		for ( j = 0; j < 3; ++ j){
			fscanf(fp,"%lf", xyznode + i * 3 + j );
		}
		//fscanf(fp,"%d",&temp);
	}

	fclose(fp);

}



void ReadData_surface_edge(void){

	int temp, i, j;
	//int *ipatchp = NULL;
	//1，读取all_surface文件
	FILE * fp = fopen("Input\\surface_edge.txt", "r");

	fscanf(fp,"%d",&maxedge);
	fscanf(fp,"%d",&temp);
	fscanf(fp,"%d",&maxedge);
	fscanf(fp,"%d",&temp);


	iedge = (int *)malloc(sizeof(int) * maxedge * 6);
	//printf("ss = %d",maxpatch); *ipatch + i * 3 + j

	for ( i = 0; i < maxedge; ++ i){
		//fscanf(fp,"%d",&temp);
		for ( j = 0; j < 6; ++ j){
			fscanf(fp,"%d", iedge + i * 6 + j );
		}
		//fscanf(fp,"%d",&temp);
	}

	fclose(fp);
 	printf("%d %d\n",iedge[0],iedge[1]);
}
#endif