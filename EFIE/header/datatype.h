#ifndef DATATYPE_H
#define DATATYPE_H

#include "cuda_runtime.h"
#include "math.h"

//cite http://blog.csdn.net/crouse246/article/details/45967861
//Real Complex Definition


#ifndef CUFCMPLX
#define CUFCMPLX
typedef float2 cuFloatComplex; 
#endif

//1,real part of a complex number
extern float cuCrealf (cuFloatComplex x);   


//2,imaginary part of a complex number
extern float cuCimagf (cuFloatComplex x);   


//3,Creat a complex number
extern cuFloatComplex make_cuFloatComplex (float r, float i);  
 

//4,conjugate of a complex number
extern cuFloatComplex cuConjf (cuFloatComplex x);  
 


//5,Add operation
extern  cuFloatComplex cuCaddf (cuFloatComplex x, cuFloatComplex y);  



//6,Minus operation
extern  cuFloatComplex cuCsubf (cuFloatComplex x, cuFloatComplex y);  


//7,multiply operation
extern  cuFloatComplex cuCmulf (cuFloatComplex x, cuFloatComplex y);  


//8,Divide operation
extern  cuFloatComplex cuCdivf (cuFloatComplex x, cuFloatComplex y);  
  

//9,Abs operation
extern  float cuCabsf (cuFloatComplex x);  
 




// Double Complex Definition
typedef double2 cuDoubleComplex;  

//1,double real part
extern double cuCreal (cuDoubleComplex x);   


//2,double imaginary part
extern double cuCimag (cuDoubleComplex x);   
 

//3,create double complex
extern  cuDoubleComplex make_cuDoubleComplex (double r, double i);  


//4,conjugate of double complex
extern  cuDoubleComplex cuConj(cuDoubleComplex x);  
 

//5,Add opertion of  double complex
extern  cuDoubleComplex cuCadd (cuDoubleComplex x, cuDoubleComplex y);  


//6,Minus opertion of  double complex
extern  cuDoubleComplex cuCsub (cuDoubleComplex x, cuDoubleComplex y);  


//7,Multiply opertion of  double complex
extern  cuDoubleComplex cuCmul (cuDoubleComplex x, cuDoubleComplex y);  


//8,Divide opertion of  double complex
extern   cuDoubleComplex cuCdiv (cuDoubleComplex x, cuDoubleComplex y);  
  

//9,abs opertion of  double complex
extern   double cuCabs (cuDoubleComplex x);  

extern   cuDoubleComplex cuCdadd(cuDoubleComplex x, double y);
extern   cuDoubleComplex cuCdmul(cuDoubleComplex x, double y);
extern  void printComplex(cuDoubleComplex x);
#endif