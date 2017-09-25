#ifndef TINY_FUNC_H
#define TINY_FUNC_H
  
double distance(double *edge);
double area(double *patch);

void cross_product(double *vecA, double *vecB, double *result);

double dot_product(double *vecA, double *vecB);

void Patch_Edge_Vector(double *patch, int m, int n, double* vector);
void Point_Point_Vector(double *vec1, double *vec2, double *vec, int m);
void NormalizeVector(double *vector, int m);
void MulScalar_Vector(double * vec, double * vec1, double num, int m);

#endif

