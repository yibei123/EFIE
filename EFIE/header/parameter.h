#ifndef PARAMETER_H
#define PARAMETER_H

#include "cuda_runtime.h"
#include "math.h"
#include "datatype.h"


#define	PI  3.141592653589793
#define MU1  4.0 * PI * 1.0e-7
#define EPS1  1.0 / (36.0 * PI) * 1.0e-9
#define RESISTANCE  sqrt(MU1 / EPS1);
#define FREQUENCY  3.0e8
#define LAMBDA  3.0e8 / FREQUENCY
#define OMEGA  2.0 * PI * FREQUENCY
#define k0  2.0 * PI / (LAMBDA*1.0)
#define k11  OMEGA * MU1 
#define k22  OMEGA * EPS1 

extern double Em[3];

   
extern double w[4];
extern double c1[4];
extern double c2[4];
extern double c3[4];
extern double kik0[3];

extern int maxedge, maxnode, maxpatch;
extern int *ipatch, *iedge;
extern double *xyznode;

typedef struct{
	double v1[3];
	double v2[3];
	int pt;
	int nt;
	double pv[3];
	double nv[3];
} RWG_Define;

extern cuDoubleComplex *ImpedanceMatrix, *Excitation; 





#endif