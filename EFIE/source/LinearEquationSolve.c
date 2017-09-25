#include "mkl.h"
#include "string.h"
#include "parameter.h"

void DirectSolve(cuDoubleComplex *impedmatrix, cuDoubleComplex *Excit, int dim, cuDoubleComplex *CurrentCoefficient){


	void *a1 = NULL, *b1 = NULL, *ipiv1 = NULL;
	int  m1, n ;

	a1 = (MKL_Complex16 *)malloc(sizeof(MKL_Complex16) * dim * dim);
	b1 = (MKL_Complex16 *)malloc(sizeof(MKL_Complex16) * dim);
	ipiv1 = (MKL_INT *)malloc(sizeof(MKL_INT) * dim);
	memcpy(a1, impedmatrix, sizeof(MKL_Complex16) * dim * dim);
	memcpy(b1, Excit, sizeof(MKL_Complex16) * dim);
	
	m1 = dim;
	LAPACKE_zgetrf (101 , m1 , m1 ,a1 , m1 , ipiv1 );
	LAPACKE_zgetrs (101 , 'N' , m1 , 1, a1 , m1 , ipiv1 ,b1 , 1);
	memcpy(CurrentCoefficient, b1, sizeof(cuDoubleComplex) * dim);
	
	free(a1);
	free(b1);
	free(ipiv1);

}