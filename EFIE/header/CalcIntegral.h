#ifndef CALCINTEGRAL_H
#define CALCINTEGRAL_H
#include "datatype.h"

cuDoubleComplex L_L_mn(double *patch1, double *patch2);
void L_L_mn_test(void);
cuDoubleComplex  Ele_Singular(double *patch1, double *patch2, double rm[][3], double rn[][3]);
cuDoubleComplex	 Mag_Singular(double *patch1, double *patch2, double rm[][3], double rn[][3]);
double	Integral_1_Test();
void Integral_1P_Test(double *p);
void Ele_Singular_Test(void);


void Excit_Cal(void);


#endif 