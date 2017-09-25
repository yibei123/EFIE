#include "stdio.h"
#include "cuda_runtime.h"
#include "datatype.h"
#include "parameter.h"
#include "DataClean.h"
#include "ImpedanceMatrix_Cal.h"
#include "CalcIntegral.h"
#include "mkl.h"
#include<string.h>
#include "LinearEquation.h"
 
int main(void)
{
	int n;
	cuDoubleComplex *RWG;
	double RCS[361];

	FILE * fp = fopen("Output\\current.txt", "w");
	FILE * fp1 = fopen("Output\\rcs.txt", "w");

	//L_L_mn_test();
	//Ele_Singular_Test();
 	//double I1P[3];
	//printf("Integral_1()=%lf\n",Integral_1_Test());
	//Integral_1P_Test(I1P);

	//1, Read Data
	ReadData_all_surface();
	ReadData_xyz_node();
	ReadData_surface_edge();

    //2, Compute ImpedanceMatrix
	printf("ss=%d\n",maxedge);
	ImpedanceMatrix_Cal();
	printComplex(ImpedanceMatrix[0]);

	//3, Compute Excit_Cal
	Excit_Cal();
	
    RWG = (cuDoubleComplex *)malloc(sizeof(cuDoubleComplex) * maxedge);

	DirectSolve(ImpedanceMatrix, Excitation, maxedge, RWG);

	

	GetRCS(RWG, RCS);

	for ( n = 0; n < maxedge; ++ n){
		fprintf(fp,"%lf,%lf\n", RWG[n].x,RWG[n].y);
	}
	fclose(fp);

	for ( n = 0; n < 361; ++ n){
		fprintf(fp1,"%lf\n", RCS[n]);
	}
	fclose(fp1);



	return 0;

}