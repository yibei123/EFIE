#include "cuda_runtime.h"
#include "math.h"
#include "stdio.h"

//cite http://blog.csdn.net/crouse246/article/details/45967861
//Real Complex Definition


typedef float2 cuFloatComplex;  

//1,real part of a complex number
__host__ __device__   float cuCrealf (cuFloatComplex x)   
{   
    return x.x;   
}  

//2,imaginary part of a complex number
__host__ __device__   float cuCimagf (cuFloatComplex x)   
{   
    return x.y;   
}  

//3,Creat a complex number
__host__ __device__   cuFloatComplex make_cuFloatComplex (float r, float i)  
{  
    cuFloatComplex res;  
    res.x = r;  
    res.y = i;  
    return res;  
}  

//4,conjugate of a complex number
__host__ __device__   cuFloatComplex cuConjf (cuFloatComplex x)  
{  
    return make_cuFloatComplex (cuCrealf(x), -cuCimagf(x));  
} 


//5,Add operation
__host__ __device__   cuFloatComplex cuCaddf (cuFloatComplex x, cuFloatComplex y)  
{  
    return make_cuFloatComplex (cuCrealf(x) + cuCrealf(y),   
                                cuCimagf(x) + cuCimagf(y));  
} 


//6,Minus operation
__host__ __device__   cuFloatComplex cuCsubf (cuFloatComplex x, cuFloatComplex y)  
{  
        return make_cuFloatComplex (cuCrealf(x) - cuCrealf(y),   
                                    cuCimagf(x) - cuCimagf(y));  
}  

//7,multiply operation
__host__ __device__   cuFloatComplex cuCmulf (cuFloatComplex x, cuFloatComplex y)  
{  
    cuFloatComplex prod;  
    prod = make_cuFloatComplex  ((cuCrealf(x) * cuCrealf(y)) - (cuCimagf(x) * cuCimagf(y)),  
                                 (cuCrealf(x) * cuCimagf(y)) + (cuCimagf(x) * cuCrealf(y)));  
    return prod;  
}  

//8,Divide operation
__host__ __device__   cuFloatComplex cuCdivf (cuFloatComplex x, cuFloatComplex y)  
{  
    cuFloatComplex quot;  
    float s = fabsf(cuCrealf(y)) + fabsf(cuCimagf(y));  
    float oos = 1.0f / s;  
    float ars = cuCrealf(x) * oos;  
    float ais = cuCimagf(x) * oos;  
    float brs = cuCrealf(y) * oos;  
    float bis = cuCimagf(y) * oos;  
    s = (brs * brs) + (bis * bis);  
    oos = 1.0f / s;  
    quot = make_cuFloatComplex (((ars * brs) + (ais * bis)) * oos, ((ais * brs) - (ars * bis)) * oos);  
    return quot;  
}  

//9,Abs operation
__host__ __device__   float cuCabsf (cuFloatComplex x)  
{  
    float a = cuCrealf(x);  
    float b = cuCimagf(x);  
    float v, w, t;  
    a = fabsf(a);  
    b = fabsf(b);  
    if (a > b)  
    {  
        v = a;  
        w = b;   
    }   
    else   
    {  
        v = b;  
        w = a;  
    }  
    t = w / v;//保证分母比分子大  
    t = 1.0f + t * t;  
    t = v * sqrtf(t);  
    if ((v == 0.0f) || (v > 3.402823466e38f) || (w > 3.402823466e38f))  
    {  
        t = v + w;  
    }  
    return t;  
}  




// Double Complex Definition
typedef double2 cuDoubleComplex;  

//1,double real part
__host__ __device__   double cuCreal (cuDoubleComplex x)   
{   
    return x.x;   
} 

//2,double imaginary part
__host__ __device__   double cuCimag (cuDoubleComplex x)   
{   
    return x.y;   
}  

//3,create double complex
__host__ __device__   cuDoubleComplex make_cuDoubleComplex (double r, double i)  
{  
    cuDoubleComplex res;  
    res.x = r;  
    res.y = i;  
    return res;  
} 

//4,conjugate of double complex
__host__ __device__   cuDoubleComplex cuConj(cuDoubleComplex x)  
{  
    return make_cuDoubleComplex (cuCreal(x), -cuCimag(x));  
}  

//5,Add opertion of  double complex
__host__ __device__   cuDoubleComplex cuCadd (cuDoubleComplex x, cuDoubleComplex y)  
{  
    return make_cuDoubleComplex (cuCreal(x) + cuCreal(y), cuCimag(x) + cuCimag(y));  
} 

//6,Minus opertion of  double complex
__host__ __device__   cuDoubleComplex cuCsub (cuDoubleComplex x, cuDoubleComplex y)  
{  
        return make_cuDoubleComplex (cuCreal(x) - cuCreal(y),   
                                    cuCimag(x) - cuCimag(y));  
} 

//7,Multiply opertion of  double complex
__host__ __device__   cuDoubleComplex cuCmul (cuDoubleComplex x, cuDoubleComplex y)  
{  
    cuDoubleComplex prod;  
    prod = make_cuDoubleComplex  ((cuCreal(x) * cuCreal(y)) - (cuCimag(x) * cuCimag(y)),  
                                 (cuCreal(x) * cuCimag(y)) + (cuCimag(x) * cuCreal(y)));  
    return prod;  
}  

//8,Divide opertion of  double complex
__host__ __device__   cuDoubleComplex cuCdiv (cuDoubleComplex x, cuDoubleComplex y)  
{  
    cuDoubleComplex quot;  
    double s = fabs(cuCreal(y)) + fabs(cuCimag(y));  
    double oos = 1.0 / s;  
    double ars = cuCreal(x) * oos;  
    double ais = cuCimag(x) * oos;  
    double brs = cuCreal(y) * oos;  
    double bis = cuCimag(y) * oos;  
    s = (brs * brs) + (bis * bis);  
    oos = 1.0 / s;  
    quot = make_cuDoubleComplex (((ars * brs) + (ais * bis)) * oos, ((ais * brs) - (ars * bis)) * oos);  
    return quot;  
}  

//9,abs opertion of  double complex
__host__ __device__   double cuCabs (cuDoubleComplex x)  
{  
    double a = cuCreal(x);  
    double b = cuCimag(x);  
    double v, w, t;  
    a = fabs(a);  
    b = fabs(b);  
    if (a > b)   
    {  
        v = a;  
        w = b;   
    } else   
    {  
        v = b;  
        w = a;  
    }  
    t = w / v;//保证分母比分子大  
    t = 1.0 + t * t;  
    t = v * sqrt(t);  
    if ((v == 0.0) || (v > 1.7e308) || (w > 1.7e308))  
    {  
        t = v + w;  
    }  
    return t;  
}

__host__ __device__  cuDoubleComplex cuCdadd(cuDoubleComplex x, double y){
	return make_cuDoubleComplex(cuCreal(x) + y, cuCimag(x));

}
__host__ __device__   cuDoubleComplex cuCdmul(cuDoubleComplex x, double y){
	return make_cuDoubleComplex(cuCreal(x) * y, cuCimag(x) * y);

}

void printComplex(cuDoubleComplex x){
	printf("real = %e, imag = %e\n",cuCreal(x),cuCimag(x));

}
