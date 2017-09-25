#include "math.h"
#include "tinyfunc.h"
#include "stdio.h"

double distance(double *edge){

	return sqrt((edge[0]-edge[3])*(edge[0]-edge[3])+(edge[1]-edge[4])*(edge[1]-edge[4])  +
		       (edge[2]-edge[5])*(edge[2]-edge[5]));

}

double area(double *patch){
    double a = distance(patch);
	double b = distance(patch + 3);
	double cp[6], c, p;

	cp[0] = patch[0];
	cp[1] = patch[1];
	cp[2] = patch[2];

	cp[3] = patch[6];
	cp[4] = patch[7];
	cp[5] = patch[8];
	c = distance(cp);



	p = (a+b+c)/2.0;


	return sqrt(p*(p-a)*(p-b)*(p-c));

}

void cross_product(double *vecA, double *vecB, double *result){
	result[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
	result[1] = vecA[2] * vecB[0] - vecA[0] * vecB[2];
	result[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
}

double dot_product(double *vecA, double *vecB){
	return vecA[0] * vecB[0] + vecA[1] * vecB[1] + vecA[2] * vecB[2];
}


void Patch_Edge_Vector(double *patch, int m, int n, double* vector){
     vector[0] = patch[m * 3] - patch[n * 3];
	 vector[1] = patch[m * 3 + 1] - patch[n * 3 + 1];
	 vector[2] = patch[m * 3 + 2] - patch[n * 3 + 2];
}

void Point_Point_Vector(double *vec1, double *vec2, double *vec, int m){
     vec[0] = vec1[0] - vec2[0];
	 vec[1] = vec1[1] - vec2[1];
	 vec[2] = vec1[2] - vec2[2];
}

void NormalizeVector(double *vector, int m){
	double sum = 0.0;
	int i;
	for(i = 0; i < m; ++i ){
		sum += vector[i] * vector[i];
	} 

	sum = sqrt(sum);
    for(i = 0; i < m; ++i ){
		vector[i] /= sum;
	} 
}

void MulScalar_Vector(double * vec, double * vec1, double num, int m){
	int i;
	for(i = 0; i < m; ++ i){
		vec1[i] = vec[i] * num;
	}
}