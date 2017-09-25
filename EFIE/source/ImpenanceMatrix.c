#include "parameter.h"
#include "datatype.h"
#include "stdlib.h"
#include "CalcIntegral.h"
#include "stdio.h"

void ImpedanceMatrix_Cal(void){
	double patch1[3][3], patch2[3][3], patch3[3][3], patch4[3][3];
	int vp, vn, vc1, vc2, m, n;
	ImpedanceMatrix = (cuDoubleComplex *)malloc(sizeof(cuDoubleComplex)*maxedge*maxedge);

	for (m = 0; m < maxedge; ++m){
		for (n = 0; n < maxedge; ++n){
			//printf("%d,%d\n", m,n);
			//printf("%f,%f,%f\n",xyznode[9],xyznode[10],xyznode[11]);
			vp = iedge[6 * m + 4] -1;
			vn = iedge[6 * m + 5] -1;
			vc1 = iedge[6 * m ]-1;
			vc2 = iedge[6 * m + 1]-1;
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


			vp = iedge[6 * n + 4]-1;
			vn = iedge[6 * n + 5]-1;
			vc1 = iedge[6 * n ]-1;
			vc2 = iedge[6 * n + 1]-1;
			patch3[0][0] = xyznode[vp * 3 ];
			patch3[0][1] = xyznode[vp * 3 + 1];
			patch3[0][2] = xyznode[vp * 3 + 2];
			patch3[1][0] = xyznode[vc1 * 3 ];
			patch3[1][1] = xyznode[vc1 * 3 + 1];
			patch3[1][2] = xyznode[vc1 * 3 + 2];
			patch3[2][0] = xyznode[vc2 * 3 ];
			patch3[2][1] = xyznode[vc2 * 3 + 1];
			patch3[2][2] = xyznode[vc2 * 3 + 2];
			patch4[0][0] = xyznode[vn * 3 ];
			patch4[0][1] = xyznode[vn * 3 + 1];
			patch4[0][2] = xyznode[vn * 3 + 2];
			patch4[2][0] = xyznode[vc1 * 3 ];
			patch4[2][1] = xyznode[vc1 * 3 + 1];
			patch4[2][2] = xyznode[vc1 * 3 + 2];
			patch4[1][0] = xyznode[vc2 * 3 ];
			patch4[1][1] = xyznode[vc2 * 3 + 1];
			patch4[1][2] = xyznode[vc2 * 3 + 2];
			//printf("ss = %d,%d,%d,%d\n",vc1,vc2,vp,vn);
		 //	printf("ss=%lf,%lf,%lf\n",patch1[0][0],patch1[0][1],patch1[0][2]);
			//printf("ss=%lf,%lf,%lf\n",patch2[1][0],patch2[1][1],patch2[1][2]);
			//printf("ss=%lf,%lf,%lf\n",patch3[0][0],patch3[0][1],patch3[0][2]);
			//printf("ss=%lf,%lf,%lf\n",patch4[1][0],patch4[1][1],patch4[1][2]);
			ImpedanceMatrix[m * maxedge + n] = cuCsub(cuCsub(L_L_mn(patch1[0], patch3[0]), L_L_mn(patch1[0], patch4[0]))  , cuCsub(L_L_mn(patch2[0], patch3[0]) , L_L_mn(patch2[0], patch4[0])));
		 //   printComplex(L_L_mn(patch1[0], patch3[0]));
			//printComplex(L_L_mn(patch1[0], patch4[0]));
			//printComplex(L_L_mn(patch2[0], patch3[0]));
			//printComplex(L_L_mn(patch2[0], patch4[0]));
			//if (m==1 && n==1) {
			//    printf("ss=%lf", cuCreal(L_L_mn(patch1[0], patch3[0])));
			//	printf("%f,%f,%f\n",patch1[0][0],patch1[0][1],patch1[0][2]);
			//}
		}
	}
}

