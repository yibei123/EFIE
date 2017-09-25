#include "mkl.h"
#include "tinyfunc.h"
#include "parameter.h"
#include "ImpedanceMatrix_Cal.h"

void ScattCal_Patch(double patch[][3], double* vecUnit, cuDoubleComplex Single_RWG_Coeff, cuDoubleComplex *vecEstt);

void GetRCS(cuDoubleComplex *RWG, double *RCS){
	double Es, theta_RCS, theta_vec[361], phi_vec[361],vecUnit[3];
	double patch1[3][3], patch2[3][3];
	int i, j, vp, vn, vc1, vc2;
	cuDoubleComplex	 VecEs[3], vecEstt[3];

	for (i = 0; i < 361; ++i){
	    theta_vec[i]=theta_vec[0]+(i*1.0/180.0)*PI;
		phi_vec[i]=0.0;

	    vecUnit[0]=sin(theta_vec[i]) * cos(phi_vec[i]);
		vecUnit[1]=sin(theta_vec[i]) * sin(phi_vec[i]);
		vecUnit[2]=cos(theta_vec[i]);

		VecEs[0] = make_cuDoubleComplex(0.0, 0.0);
		VecEs[1] = make_cuDoubleComplex(0.0, 0.0);
		VecEs[2] = make_cuDoubleComplex(0.0, 0.0);

		for(j = 0; j < maxedge; ++j){
			vp = iedge[6 * j + 4] -1;
			vn = iedge[6 * j + 5] -1;
			vc1 = iedge[6 * j ]-1;
			vc2 = iedge[6 * j + 1]-1;
			patch1[0][0] = xyznode[vp * 3 ];
			patch1[0][1] = xyznode[vp * 3 + 1];
			patch1[0][2] = xyznode[vp * 3 + 2];
			patch1[1][0] = xyznode[vc1 * 3 ];
			patch1[1][1] = xyznode[vc1 * 3 + 1];
			patch1[1][2] = xyznode[vc1 * 3 + 2];
			patch1[2][0] = xyznode[vc2 * 3 ];
			patch1[2][1] = xyznode[vc2 * 3 + 1];
			patch1[2][2] = xyznode[vc2 * 3 + 2];
			patch2[0][0] = xyznode[vn * 3 ];
			patch2[0][1] = xyznode[vn * 3 + 1];
			patch2[0][2] = xyznode[vn * 3 + 2];
			patch2[2][0] = xyznode[vc1 * 3 ];
			patch2[2][1] = xyznode[vc1 * 3 + 1];
			patch2[2][2] = xyznode[vc1 * 3 + 2];
			patch2[1][0] = xyznode[vc2 * 3 ];
			patch2[1][1] = xyznode[vc2 * 3 + 1];
			patch2[1][2] = xyznode[vc2 * 3 + 2];

			ScattCal_Patch(patch1, vecUnit, RWG[j], vecEstt);
			VecEs[0] = cuCadd(VecEs[0], vecEstt[0]);
			VecEs[1] = cuCadd(VecEs[1], vecEstt[1]);
			VecEs[2] = cuCadd(VecEs[2], vecEstt[2]);

			ScattCal_Patch(patch2, vecUnit, RWG[j], vecEstt);
			VecEs[0] = cuCsub(VecEs[0], vecEstt[0]);
			VecEs[1] = cuCsub(VecEs[1], vecEstt[1]);
			VecEs[2] = cuCsub(VecEs[2], vecEstt[2]);
		} 
		Es = (VecEs[0].x)*(VecEs[0].x) + (VecEs[1].x)*(VecEs[1].x) + (VecEs[2].x)*(VecEs[2].x) +
			 (VecEs[0].y)*(VecEs[0].y) + (VecEs[1].y)*(VecEs[1].y) + (VecEs[2].y)*(VecEs[2].y);

	    Es=(k11*k11)*Es/(4.0*PI);
		RCS[i]=10.0*log10(Es); 
	}


}


void ScattCal_Patch(double patch[][3], double* vecUnit, cuDoubleComplex Single_RWG_Coeff, cuDoubleComplex *vecEstt){

	double ln = distance(patch[1]), rm[3], vec1[3], Temp_real, Temp_imag;
	int m, n;
	cuDoubleComplex	J1_Pos[3], vecJ_Posvec, Temp, CmplxJ;	vecEstt[0] = make_cuDoubleComplex(0.0, 0.0);
	vecEstt[1] = make_cuDoubleComplex(0.0, 0.0);
	vecEstt[2] = make_cuDoubleComplex(0.0, 0.0);
	CmplxJ = make_cuDoubleComplex(0.0, -1.0);

    for(m = 0; m < 4; ++ m){
		for (n = 0; n <3; ++ n){
		    rm[n] = c1[m] * patch[0][n] + c2[m] * patch[1][n] + c3[m] * patch[2][n];
		}

		Temp_real = w[m] * cos(k0 * dot_product(rm, vecUnit));
		Temp_imag = w[m] * sin(k0 * dot_product(rm, vecUnit));
		Temp = make_cuDoubleComplex(Temp_real,Temp_imag);

		Point_Point_Vector(rm, patch[0], vec1, 3);		
		MulScalar_Vector(vec1, vec1, ln/2.0, 3);
		J1_Pos[0] = cuCdmul(Single_RWG_Coeff, vec1[0]);
		J1_Pos[1] = cuCdmul(Single_RWG_Coeff, vec1[1]);
		J1_Pos[2] = cuCdmul(Single_RWG_Coeff, vec1[2]);

		vecJ_Posvec =  cuCadd(cuCadd(cuCdmul(J1_Pos[0], vecUnit[0]), cuCdmul(J1_Pos[1], vecUnit[1])) ,cuCdmul(J1_Pos[2], vecUnit[2]));
		J1_Pos[0] =cuCsub(J1_Pos[0], cuCdmul(vecJ_Posvec, vecUnit[0])) ;
		J1_Pos[1] =cuCsub(J1_Pos[1], cuCdmul(vecJ_Posvec, vecUnit[1])) ;
		J1_Pos[2] =cuCsub(J1_Pos[2], cuCdmul(vecJ_Posvec, vecUnit[2])) ;

		vecEstt[0] =  cuCadd(vecEstt[0], cuCmul(cuCmul(Temp, J1_Pos[0]),CmplxJ)) ;
		vecEstt[1] =  cuCadd(vecEstt[1], cuCmul(cuCmul(Temp, J1_Pos[1]),CmplxJ));
		vecEstt[2] =  cuCadd(vecEstt[2], cuCmul(cuCmul(Temp, J1_Pos[2]),CmplxJ)) ;


	}
	 

}