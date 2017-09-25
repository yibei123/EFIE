#include "datatype.h"
#include "tinyfunc.h"
#include "parameter.h"
#include "stdlib.h"

void Excit_Cal(void){
	int m, n, i, vp, vn, vc1, vc2;
	double patch1[3][3], patch2[3][3], rm[4][3], rn[4][3];
	double lm, temp[3], mid_temp_real, mid_temp_imag;
	Excitation = (cuDoubleComplex *)malloc(sizeof(cuDoubleComplex)*maxedge);

	for (i = 0; i < maxedge; ++i){
		vp = iedge[6 * i + 4] -1;
		vn = iedge[6 * i + 5] -1;
		vc1 = iedge[6 * i ]-1;
		vc2 = iedge[6 * i + 1]-1;
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

		Excitation[i] = make_cuDoubleComplex(0.0,0.0);

		lm = distance(patch1[1]);

		for(m = 0; m < 4; ++ m){
			for(n = 0; n < 3; ++ n){
				rm[m][n] = c1[m] * patch1[0][n] + c2[m] * patch1[1][n] + c3[m] * patch1[2][n];
				rn[m][n] = c1[m] * patch2[0][n] + c2[m] * patch2[1][n] + c3[m] * patch2[2][n];
			}
			Point_Point_Vector(rm[m], patch1[0], temp, 3);
			mid_temp_real = w[m]*lm*dot_product(temp, Em) * cos(dot_product(kik0, rm[m]))/2.0;
			mid_temp_imag = -w[m]*lm*dot_product(temp, Em) * sin(dot_product(kik0, rm[m]))/2.0;
			Excitation[i] = cuCadd(Excitation[i], make_cuDoubleComplex(mid_temp_real,mid_temp_imag));

			Point_Point_Vector(rn[m], patch2[0], temp, 3);
			mid_temp_real = -w[m]*lm*dot_product(temp, Em) * cos(dot_product(kik0, rn[m]))/2.0;
			mid_temp_imag = w[m]*lm*dot_product(temp, Em) * sin(dot_product(kik0, rn[m]))/2.0;
			Excitation[i] = cuCadd(Excitation[i], make_cuDoubleComplex(mid_temp_real,mid_temp_imag));
		} 
	}

}
	
