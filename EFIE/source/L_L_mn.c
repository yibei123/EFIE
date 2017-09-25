#include "datatype.h"
#include "tinyfunc.h"
#include "parameter.h"
#include "CalcIntegral.h"
#include "stdio.h"

cuDoubleComplex L_L_mn(double *patch1, double *patch2){
	int m, n;
	double lm = distance(patch1 + 3), rm[4][3], rn[4][3], temp1[3], temp2[3];
	double ln = distance(patch2 + 3), Cent_1[3], Cent_2[3], edge[6], R;
	cuDoubleComplex Coeff_M = make_cuDoubleComplex(0.0, lm * ln * k11/ (16.0 * PI));
	cuDoubleComplex Coeff_E = make_cuDoubleComplex(0.0, - lm * ln / (4.0 * PI * k22));
	cuDoubleComplex Ele_1 = make_cuDoubleComplex(0.0, 0.0), Ele_2 = make_cuDoubleComplex(0.0, 0.0);



	for(m = 0; m < 4; ++ m){
		for(n = 0; n < 3; ++ n){
			rm[m][n] = c1[m] * patch1[n] + c2[m] * patch1[n + 3] + c3[m] * patch1[n + 6];
			rn[m][n] = c1[m] * patch2[n] + c2[m] * patch2[n + 3] + c3[m] * patch2[n + 6];
		}
	} 


	for (m = 0; m < 3; ++ m){
		Cent_1[m] = (patch1[m] + patch1[m + 3] + patch1[m + 6])/3.0;
		Cent_2[m] = (patch2[m] + patch2[m + 3] + patch2[m + 6])/3.0;
		Cent_1[m] -= Cent_2[m]; 
	}

	//printf("(dot_product(Cent_1, Cent_1) = %lf\n",dot_product(Cent_1, Cent_1));


	if ( (dot_product(Cent_1, Cent_1)) < 1e-6) {
	     //printf("sssss");
	     Ele_2 = Ele_Singular(patch1, patch2, rm, rn);
 /*        printf("sssss");*/
	     Ele_1 = Mag_Singular(patch1, patch2, rm, rn);

	}
	else{
		for(m = 0; m < 4; ++m){
			for(n = 0; n < 4; ++n){
			   	edge[0] = rm[m][0];
				edge[1] = rm[m][1];
				edge[2] = rm[m][2];
				edge[3] = rn[n][0];
				edge[4] = rn[n][1];
				edge[5] = rn[n][2];
				R = distance(edge);
				Point_Point_Vector(rm[m], patch1, temp1, 3);
				Point_Point_Vector(rn[n], patch2, temp2, 3);
				Ele_1 = cuCadd(Ele_1, cuCdmul(make_cuDoubleComplex(cos(k0*R),-sin(k0*R)),w[m]*w[n]*dot_product(temp1,temp2)/R ) );
                Ele_2 = cuCadd(Ele_2, cuCdmul(make_cuDoubleComplex(cos(k0*R),-sin(k0*R)),w[m]*w[n]/R));
			}
		}
	
	}

	return 	cuCadd(cuCmul(Ele_1, Coeff_M),cuCmul(Ele_2, Coeff_E)); 
	

}

void L_L_mn_test(void){
	int m, n;
    double patch1[9] = {0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0};//, rm[3] = {0.5,0.5,0.0};
	double patch2[9] = {1.0,0.0,0.0,2.0,0.0,0.0,1.0,1.0,0.0};
	double lm = distance(patch1 + 3), rm[4][3], rn[4][3], temp1[3], temp2[3];
	double ln = distance(patch2 + 3), Cent_1[3], Cent_2[3], edge[6], R;
	cuDoubleComplex result, Coeff_M = make_cuDoubleComplex(0.0, lm * ln * k11/ (16.0 * PI));
	cuDoubleComplex Coeff_E = make_cuDoubleComplex(0.0, - lm * ln / (4.0 * PI * k22));
	cuDoubleComplex Ele_1 = make_cuDoubleComplex(0.0, 0.0), Ele_2 = make_cuDoubleComplex(0.0, 0.0);



	for(m = 0; m < 4; ++ m){
		for(n = 0; n < 3; ++ n){
			rm[m][n] = c1[m] * patch1[n] + c2[m] * patch1[n + 3] + c3[m] * patch1[n + 6];
			rn[m][n] = c1[m] * patch2[n] + c2[m] * patch2[n + 3] + c3[m] * patch2[n + 6];
		}
	} 


	for (m = 0; m < 3; ++ m){
		Cent_1[m] = (patch1[m] + patch1[m + 3] + patch1[m + 6])/3.0;
		Cent_2[m] = (patch2[m] + patch2[m + 3] + patch2[m + 6])/3.0;
		Cent_1[m] -= Cent_2[m]; 
	}

	//printf("(dot_product(Cent_1, Cent_1) = %lf\n",dot_product(Cent_1, Cent_1));


	if ( (dot_product(Cent_1, Cent_1)) < 1e-6) {
	     //printf("sssss");
	     Ele_2 = Ele_Singular(patch1, patch2, rm, rn);
 /*        printf("sssss");*/
	     Ele_1 = Mag_Singular(patch1, patch2, rm, rn);

	}
	else{
		for(m = 0; m < 4; ++m){
			for(n = 0; n < 4; ++n){
			   	edge[0] = rm[m][0];
				edge[1] = rm[m][1];
				edge[2] = rm[m][2];
				edge[3] = rn[n][0];
				edge[4] = rn[n][1];
				edge[5] = rn[n][2];
				R = distance(edge);
				Point_Point_Vector(rm[m], patch1, temp1, 3);
				Point_Point_Vector(rn[n], patch2, temp2, 3);
				Ele_1 = cuCadd(Ele_1, cuCdmul(make_cuDoubleComplex(cos(k0*R),-sin(k0*R)),w[m]*w[n]*dot_product(temp1,temp2)/R ) );
                Ele_2 = cuCadd(Ele_2, cuCdmul(make_cuDoubleComplex(cos(k0*R),-sin(k0*R)),w[m]*w[n]/R));
			}
		}
	
	}

	result  = 	cuCadd(cuCmul(Ele_1, Coeff_M),cuCmul(Ele_2, Coeff_E)); 
	//result  = 	cuCadd(cuCmul(Ele_1, Coeff_M),cuCmul(Ele_2, Coeff_E)); 

}