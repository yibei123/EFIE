#include "datatype.h"
#include "tinyfunc.h"
#include "parameter.h"
#include "stdio.h"


double	Integral_1(double *patch,double *rm);

cuDoubleComplex  Ele_Singular(double *patch1, double *patch2, double rm[][3], double rn[][3]){
	int m, n;
	double edge[6], R;
	cuDoubleComplex  result = make_cuDoubleComplex(0.0,0.0);


	for (m = 0; m < 4; ++ m){
		for (n = 0; n < 4; ++n){
			edge[0] = rm[m][0];
			edge[1] = rm[m][1];
			edge[2] = rm[m][2];
			edge[3] = rn[n][0];
			edge[4] = rn[n][1];
			edge[5] = rn[n][2];
			R = distance(edge);
			if (R < 1e-8){
				
				result = cuCsub(result, make_cuDoubleComplex(0.0, w[m] * w[n] * k0)); 
			}
			else {
				result = make_cuDoubleComplex(cuCreal(result)+ w[m] * w[n] * (cos(k0 * R)-1.0)/R, 
					      cuCimag(result) - w[m] * w[n] * sin(k0 * R)/R);
			
			}
		
		}
	    result = cuCdadd(result, Integral_1(patch2,rm[m])*w[m]/area(patch2));
	}


	//temp = cuCreal(result);
    //printf("%lf,%lf\n",cuCreal(result),cuCimag(result)); 


	return 	result;


}






double	Integral_1(double *patch,double *rm){
     double unit_Pi0[9], temp_1[3], temp_2[3], total_integral =  0.0;
	 double edge[6],  vecUnitN[3], temp[3];
	 double Ivec[9], UI[9], d, Pi0[3], Ri02[3], lip[3], lin[3], Rip[3], Rin[3];
	 int m;

	 Patch_Edge_Vector(patch, 1, 0, edge);
	 Patch_Edge_Vector(patch, 2, 1, edge + 3);
	 cross_product(edge, edge+3, vecUnitN);
	 NormalizeVector(vecUnitN, 3);

	 temp[0] = patch[0] - rm[0];
	 temp[1] = patch[1] - rm[1];
	 temp[2] = patch[2] - rm[2];

	 d = fabs(dot_product(temp, vecUnitN));
	 //printf("d = %lf\n",d);



	 Patch_Edge_Vector(patch, 1, 0, Ivec);
	 NormalizeVector(Ivec, 3);
	 cross_product(Ivec, vecUnitN, UI);
	 temp[0] = patch[0] - rm[0];
	 temp[1] = patch[1] - rm[1];
	 temp[2] = patch[2] - rm[2];
	 Pi0[0] = fabs(dot_product(temp, UI));
	 Ri02[0] = Pi0[0]*Pi0[0] + d * d;

	 Patch_Edge_Vector(patch, 0, 2, Ivec + 3);
	 NormalizeVector(Ivec + 3, 3);
	 cross_product(Ivec + 3, vecUnitN, UI + 3);
	 temp[0] = patch[6] - rm[0];
	 temp[1] = patch[7] - rm[1];
	 temp[2] = patch[8] - rm[2];
	 Pi0[1] = fabs(dot_product(temp, UI + 3));
	 Ri02[1] = Pi0[1]*Pi0[1] + d * d;

	 Patch_Edge_Vector(patch, 2, 1, Ivec + 6);
	 NormalizeVector(Ivec + 6, 3);
	 cross_product(Ivec + 6, vecUnitN, UI + 6);
	 temp[0] = patch[3] - rm[0];
	 temp[1] = patch[4] - rm[1];
	 temp[2] = patch[5] - rm[2];
	 Pi0[2] = fabs(dot_product(temp, UI + 6));
	 Ri02[2] = Pi0[2]*Pi0[2] + d * d;



	 Point_Point_Vector(patch, rm, temp, 3);
	 lip[1] = dot_product(Ivec + 3, temp);
	 lin[0] = dot_product(Ivec, temp);

	 Point_Point_Vector(patch + 3, rm, temp, 3);
	 lip[0] = dot_product(Ivec, temp);
	 lin[2] = dot_product(Ivec + 6, temp);

	 Point_Point_Vector(patch + 6, rm, temp, 3);
	 lip[2] = dot_product(Ivec + 6, temp);
	 lin[1] = dot_product(Ivec + 3, temp);
	
	 


	 Point_Point_Vector(patch, rm, temp, 3);
	 Rip[1] = sqrt(dot_product(temp, temp));
	 Rin[0]	= Rip[1];

	 Point_Point_Vector(patch + 3, rm, temp, 3);
	 Rip[0] = sqrt(dot_product(temp, temp));
	 Rin[2]	= Rip[0];

	 Point_Point_Vector(patch + 6, rm, temp, 3);
	 Rip[2] = sqrt(dot_product(temp, temp));
	 Rin[1]	= Rip[2];


	 //printf("Rip = %lf,%lf,%lf\n",Rip[0],Rip[1],Rip[2]);
	 //printf("Rin = %lf,%lf,%lf\n",Rin[0],Rin[1],Rin[2]);
	 //printf("lip = %lf,%lf,%lf\n",lip[0],lip[1],lip[2]);
	 //printf("lin = %lf,%lf,%lf\n",lin[0],lin[1],lin[2]);
	 //printf("Pi0 = %lf,%lf,%lf\n",Pi0[0],Pi0[1],Pi0[2]);
	 //printf("Ri02 = %lf,%lf,%lf\n",Ri02[0],Ri02[1],Ri02[2]);



	 for (m=0; m<3; m++){
		 if (fabs(Pi0[m]) > 1e-8){
			
			if (m == 2){
				
				Point_Point_Vector(patch + 3 * m, rm, temp, 3);
			}
			else{
				Point_Point_Vector(patch + 3 * (1-m), rm, temp, 3);
			}
		
			MulScalar_Vector(Ivec + 3 * m, temp_1, lip[m], 3);
			Point_Point_Vector(temp, temp_1, temp_2, 3);
			MulScalar_Vector(temp_2, unit_Pi0 + 3 * m, 1.0/Pi0[m], 3);		
		 }

	 }
	



	for (m = 0; m < 3; ++m){
		if (fabs(Pi0[m]) < 1e-8){
			total_integral += 0.0;
		}
		else{
		    //printf("dot_product(unit_Pi0 + 3 * m, UI + 3 * m) = %lf\n",dot_product(unit_Pi0 + 3 * m, UI + 3 * m));
			total_integral += dot_product(unit_Pi0 + 3 * m, UI + 3 * m) * (
				     Pi0[m] * log((Rip[m]+lip[m])/(Rin[m]+lin[m])) - 
					fabs(d)*(atan(Pi0[m]*lip[m]/(Ri02[m]+fabs(d)*Rip[m]))- atan(Pi0[m]*lin[m]/(Ri02[m]+fabs(d)*Rin[m]))) );

		}
	} 
   return total_integral;

}





void Integral_1P(double *patch,double *rm,double *I1P);
//==============================================================
cuDoubleComplex	 Mag_Singular(double *patch1, double *patch2, double rm[][3], double rn[][3]){

	int m, n;
	double edge[6], R, temp1[3], temp2[3];
	double I1P[3];
	cuDoubleComplex  result = make_cuDoubleComplex(0.0,0.0);


	for (m = 0; m < 4; ++ m){
			edge[0] = rm[m][0];
			edge[1] = rm[m][1];
			edge[2] = rm[m][2];
			Point_Point_Vector(rm[m], patch1, temp1, 3);
		for (n = 0; n < 4; ++n){
			edge[3] = rn[n][0];
			edge[4] = rn[n][1];
			edge[5] = rn[n][2];
			Point_Point_Vector(rn[n], patch2, temp2, 3);
			R = distance(edge);
			if (R < 1e-8){
				result = cuCsub(result, make_cuDoubleComplex(0.0, w[m] * w[n] * dot_product(temp1,temp2) * k0)); 
			}
			else {
				result = make_cuDoubleComplex(cuCreal(result)+ w[m] * w[n] * dot_product(temp1,temp2)* (cos(k0 * R)-1.0)/R, 
					      cuCimag(result) - w[m] * w[n] *dot_product(temp1,temp2) * sin(k0 * R)/R);
			
			}
		
		}
		
		Point_Point_Vector(rm[m], patch2, temp2, 3);
		MulScalar_Vector(temp2, temp2,Integral_1(patch2,rm[m]),3);
		Integral_1P(patch2,rm[m], I1P);

		temp2[0] += I1P[0];
		temp2[1] += I1P[1];
		temp2[2] += I1P[2];
	    result = cuCdadd(result, w[m] * dot_product(temp1, temp2) /area(patch2));
	}
		
	//printf("mag = %lf,%lf\n",cuCreal(result),cuCimag(result));
	return 	result;

}

void  Integral_1P(double *patch,double *rm,double *I1P){
     
	 double unit_Pi0[9], temp_1[3], temp_2[3];
	 double edge[6], vecUnitN[3], temp[3];
	 double Ivec[9], UI[9], d, Pi0[3], Ri02[3], lip[3], lin[3], Rip[3], Rin[3];
	 int m;

	Patch_Edge_Vector(patch, 1, 0, edge);
	 Patch_Edge_Vector(patch, 2, 1, edge + 3);
	 cross_product(edge, edge+3, vecUnitN);
	 NormalizeVector(vecUnitN, 3);

	 temp[0] = patch[0] - rm[0];
	 temp[1] = patch[1] - rm[1];
	 temp[2] = patch[2] - rm[2];

	 d = fabs(dot_product(temp, vecUnitN));
	 //printf("d = %lf\n",d);



	 Patch_Edge_Vector(patch, 1, 0, Ivec);
	 NormalizeVector(Ivec, 3);
	 cross_product(Ivec, vecUnitN, UI);
	 temp[0] = patch[0] - rm[0];
	 temp[1] = patch[1] - rm[1];
	 temp[2] = patch[2] - rm[2];
	 Pi0[0] = fabs(dot_product(temp, UI));
	 Ri02[0] = Pi0[0]*Pi0[0] + d * d;

	 Patch_Edge_Vector(patch, 0, 2, Ivec + 3);
	 NormalizeVector(Ivec + 3, 3);
	 cross_product(Ivec + 3, vecUnitN, UI + 3);
	 temp[0] = patch[6] - rm[0];
	 temp[1] = patch[7] - rm[1];
	 temp[2] = patch[8] - rm[2];
	 Pi0[1] = fabs(dot_product(temp, UI + 3));
	 Ri02[1] = Pi0[1]*Pi0[1] + d * d;

	 Patch_Edge_Vector(patch, 2, 1, Ivec + 6);
	 NormalizeVector(Ivec + 6, 3);
	 cross_product(Ivec + 6, vecUnitN, UI + 6);
	 temp[0] = patch[3] - rm[0];
	 temp[1] = patch[4] - rm[1];
	 temp[2] = patch[5] - rm[2];
	 Pi0[2] = fabs(dot_product(temp, UI + 6));
	 Ri02[2] = Pi0[2]*Pi0[2] + d * d;



	 Point_Point_Vector(patch, rm, temp, 3);
	 lip[1] = dot_product(Ivec + 3, temp);
	 lin[0] = dot_product(Ivec, temp);

	 Point_Point_Vector(patch + 3, rm, temp, 3);
	 lip[0] = dot_product(Ivec, temp);
	 lin[2] = dot_product(Ivec + 6, temp);

	 Point_Point_Vector(patch + 6, rm, temp, 3);
	 lip[2] = dot_product(Ivec + 6, temp);
	 lin[1] = dot_product(Ivec + 3, temp);
	
	 


	 Point_Point_Vector(patch, rm, temp, 3);
	 Rip[1] = sqrt(dot_product(temp, temp));
	 Rin[0]	= Rip[1];

	 Point_Point_Vector(patch + 3, rm, temp, 3);
	 Rip[0] = sqrt(dot_product(temp, temp));
	 Rin[2]	= Rip[0];

	 Point_Point_Vector(patch + 6, rm, temp, 3);
	 Rip[2] = sqrt(dot_product(temp, temp));
	 Rin[1]	= Rip[2];


	 //printf("Rip = %lf,%lf,%lf\n",Rip[0],Rip[1],Rip[2]);
	 //printf("Rin = %lf,%lf,%lf\n",Rin[0],Rin[1],Rin[2]);
	 //printf("lip = %lf,%lf,%lf\n",lip[0],lip[1],lip[2]);
	 //printf("lin = %lf,%lf,%lf\n",lin[0],lin[1],lin[2]);
	 //printf("Pi0 = %lf,%lf,%lf\n",Pi0[0],Pi0[1],Pi0[2]);
	 //printf("Ri02 = %lf,%lf,%lf\n",Ri02[0],Ri02[1],Ri02[2]);



	 for (m=0; m<3; m++){
		 if (fabs(Pi0[m]) > 1e-8){
			
			if (m == 2){
				//printf("%d\n",m);
				Point_Point_Vector(patch + 3 * m, rm, temp, 3);
			}
			else{
				Point_Point_Vector(patch + 3 * (1-m), rm, temp, 3);
			}
		
			MulScalar_Vector(Ivec + 3 * m, temp_1, lip[m], 3);
			Point_Point_Vector(temp, temp_1, temp_2, 3);
			MulScalar_Vector(temp_2, unit_Pi0 + 3 * m, 1.0/Pi0[m], 3);		
		 }

	 }

	I1P[0] = 0.0;
	I1P[1] = 0.0;
	I1P[2] = 0.0;

	for (m = 0; m < 3; ++m){

		temp[0] = Ri02[m]*log( (Rip[m]+lip[m])/(Rin[m]+lin[m]) )+lip[m]*Rip[m]-lin[m]*Rin[m]; 
		MulScalar_Vector(UI + 3 * m, temp_1, temp[0], 3);
	   	I1P[0] += temp_1[0];
		I1P[1] += temp_1[1];
		I1P[2] += temp_1[2];
	} 

	I1P[0] /= 2.0;
	I1P[1] /= 2.0;
	I1P[2] /= 2.0;


}






void Ele_Singular_Test(void){
	double patch1[9] = {0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0}, rm[4][3];
	double patch2[9] = {0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0}, rn[4][3];
	int m,n;
	cuDoubleComplex	 result ;
	for(m = 0; m < 4; ++ m){
		for(n = 0; n < 3; ++ n){
			rm[m][n] = c1[m] * patch1[n] + c2[m] * patch1[n + 3] + c3[m] * patch1[n + 6];
			rn[m][n] = c1[m] * patch2[n] + c2[m] * patch2[n + 3] + c3[m] * patch2[n + 6];
		}
	}



	result = Ele_Singular(patch1,patch2,rm,rn);

	result = Mag_Singular(patch1,patch2,rm,rn);


}







// Test for Integral 
double	Integral_1_Test(){
     double unit_Pi0[9], temp_1[3], temp_2[3], total_integral =  0.0;
	 double edge[6],  vecUnitN[3], temp[3];
	 double Ivec[9], UI[9], d, Pi0[3], Ri02[3], lip[3], lin[3], Rip[3], Rin[3];
	 int m;
	 double patch[9] = {0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0}, rm[3] = {0.5,0.5,0.0};



	 Patch_Edge_Vector(patch, 1, 0, edge);
	 Patch_Edge_Vector(patch, 2, 1, edge + 3);
	 cross_product(edge, edge+3, vecUnitN);
	 NormalizeVector(vecUnitN, 3);

	 temp[0] = patch[0] - rm[0];
	 temp[1] = patch[1] - rm[1];
	 temp[2] = patch[2] - rm[2];

	 d = fabs(dot_product(temp, vecUnitN));
	 //printf("d = %lf\n",d);



	 Patch_Edge_Vector(patch, 1, 0, Ivec);
	 NormalizeVector(Ivec, 3);
	 cross_product(Ivec, vecUnitN, UI);
	 temp[0] = patch[0] - rm[0];
	 temp[1] = patch[1] - rm[1];
	 temp[2] = patch[2] - rm[2];
	 Pi0[0] = fabs(dot_product(temp, UI));
	 Ri02[0] = Pi0[0]*Pi0[0] + d * d;

	 Patch_Edge_Vector(patch, 0, 2, Ivec + 3);
	 NormalizeVector(Ivec + 3, 3);
	 cross_product(Ivec + 3, vecUnitN, UI + 3);
	 temp[0] = patch[6] - rm[0];
	 temp[1] = patch[7] - rm[1];
	 temp[2] = patch[8] - rm[2];
	 Pi0[1] = fabs(dot_product(temp, UI + 3));
	 Ri02[1] = Pi0[1]*Pi0[1] + d * d;

	 Patch_Edge_Vector(patch, 2, 1, Ivec + 6);
	 NormalizeVector(Ivec + 6, 3);
	 cross_product(Ivec + 6, vecUnitN, UI + 6);
	 temp[0] = patch[3] - rm[0];
	 temp[1] = patch[4] - rm[1];
	 temp[2] = patch[5] - rm[2];
	 Pi0[2] = fabs(dot_product(temp, UI + 6));
	 Ri02[2] = Pi0[2]*Pi0[2] + d * d;



	 Point_Point_Vector(patch, rm, temp, 3);
	 lip[1] = dot_product(Ivec + 3, temp);
	 lin[0] = dot_product(Ivec, temp);

	 Point_Point_Vector(patch + 3, rm, temp, 3);
	 lip[0] = dot_product(Ivec, temp);
	 lin[2] = dot_product(Ivec + 6, temp);

	 Point_Point_Vector(patch + 6, rm, temp, 3);
	 lip[2] = dot_product(Ivec + 6, temp);
	 lin[1] = dot_product(Ivec + 3, temp);
	
	 


	 Point_Point_Vector(patch, rm, temp, 3);
	 Rip[1] = sqrt(dot_product(temp, temp));
	 Rin[0]	= Rip[1];

	 Point_Point_Vector(patch + 3, rm, temp, 3);
	 Rip[0] = sqrt(dot_product(temp, temp));
	 Rin[2]	= Rip[0];

	 Point_Point_Vector(patch + 6, rm, temp, 3);
	 Rip[2] = sqrt(dot_product(temp, temp));
	 Rin[1]	= Rip[2];


	 //printf("Rip = %lf,%lf,%lf\n",Rip[0],Rip[1],Rip[2]);
	 //printf("Rin = %lf,%lf,%lf\n",Rin[0],Rin[1],Rin[2]);
	 //printf("lip = %lf,%lf,%lf\n",lip[0],lip[1],lip[2]);
	 //printf("lin = %lf,%lf,%lf\n",lin[0],lin[1],lin[2]);
	 //printf("Pi0 = %lf,%lf,%lf\n",Pi0[0],Pi0[1],Pi0[2]);
	 //printf("Ri02 = %lf,%lf,%lf\n",Ri02[0],Ri02[1],Ri02[2]);



	 for (m=0; m<3; m++){
		 if (fabs(Pi0[m]) > 1e-8){
			
			if (m == 2){
				//printf("%d\n",m);
				Point_Point_Vector(patch + 3 * m, rm, temp, 3);
			}
			else{
				Point_Point_Vector(patch + 3 * (1-m), rm, temp, 3);
			}
		
			MulScalar_Vector(Ivec + 3 * m, temp_1, lip[m], 3);
			Point_Point_Vector(temp, temp_1, temp_2, 3);
			MulScalar_Vector(temp_2, unit_Pi0 + 3 * m, 1.0/Pi0[m], 3);		
		 }
		 temp_2[0] = Pi0[0];
		 temp_2[1] = Pi0[1];
		 temp_2[2] = Pi0[2];

	 }
	



	for (m = 0; m < 3; ++m){
		if (fabs(Pi0[m]) < 1e-4){
			total_integral += 0.0;
	     //   printf("Pi0 = %lf,%d\n",fabs(Pi0[m]),m);//,Pi0[1],Pi0[2]);
		    //printf("Pi0 = %lf,%lf,%lf\n",Pi0[0],Pi0[1],Pi0[2]);

		}
		else{
		    //printf("dot_product(unit_Pi0 + 3 * m, UI + 3 * m) = %lf\n",dot_product(unit_Pi0 + 3 * m, UI + 3 * m));
			total_integral += dot_product(unit_Pi0 + 3 * m, UI + 3 * m) * (
				     Pi0[m] * log((Rip[m]+lip[m])/(Rin[m]+lin[m])) - 
					fabs(d)*(atan(Pi0[m]*lip[m]/(Ri02[m]+fabs(d)*Rip[m]))- atan(Pi0[m]*lin[m]/(Ri02[m]+fabs(d)*Rin[m]))) );

		}
	} 
   return total_integral;

}



void  Integral_1P_Test(double *I1P){
     
	 double unit_Pi0[9], temp_1[3], temp_2[3];
	 double edge[6], vecUnitN[3], temp[3];
	 double Ivec[9], UI[9], d, Pi0[3], Ri02[3], lip[3], lin[3], Rip[3], Rin[3];
	 int m;
	 double patch[9] = {0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0}, rm[3] = {0.5,0.5,0.0};

	 Patch_Edge_Vector(patch, 1, 0, edge);
	 Patch_Edge_Vector(patch, 2, 1, edge + 3);
	 cross_product(edge, edge+3, vecUnitN);
	 NormalizeVector(vecUnitN, 3);

	 temp[0] = patch[0] - rm[0];
	 temp[1] = patch[1] - rm[1];
	 temp[2] = patch[2] - rm[2];

	 d = fabs(dot_product(temp, vecUnitN));
	 //printf("d = %lf\n",d);



	 Patch_Edge_Vector(patch, 1, 0, Ivec);
	 NormalizeVector(Ivec, 3);
	 cross_product(Ivec, vecUnitN, UI);
	 temp[0] = patch[0] - rm[0];
	 temp[1] = patch[1] - rm[1];
	 temp[2] = patch[2] - rm[2];
	 Pi0[0] = fabs(dot_product(temp, UI));
	 Ri02[0] = Pi0[0]*Pi0[0] + d * d;

	 Patch_Edge_Vector(patch, 0, 2, Ivec + 3);
	 NormalizeVector(Ivec + 3, 3);
	 cross_product(Ivec + 3, vecUnitN, UI + 3);
	 temp[0] = patch[6] - rm[0];
	 temp[1] = patch[7] - rm[1];
	 temp[2] = patch[8] - rm[2];
	 Pi0[1] = fabs(dot_product(temp, UI + 3));
	 Ri02[1] = Pi0[1]*Pi0[1] + d * d;

	 Patch_Edge_Vector(patch, 2, 1, Ivec + 6);
	 NormalizeVector(Ivec + 6, 3);
	 cross_product(Ivec + 6, vecUnitN, UI + 6);
	 temp[0] = patch[3] - rm[0];
	 temp[1] = patch[4] - rm[1];
	 temp[2] = patch[5] - rm[2];
	 Pi0[2] = fabs(dot_product(temp, UI + 6));
	 Ri02[2] = Pi0[2]*Pi0[2] + d * d;



	 Point_Point_Vector(patch, rm, temp, 3);
	 lip[1] = dot_product(Ivec + 3, temp);
	 lin[0] = dot_product(Ivec, temp);

	 Point_Point_Vector(patch + 3, rm, temp, 3);
	 lip[0] = dot_product(Ivec, temp);
	 lin[2] = dot_product(Ivec + 6, temp);

	 Point_Point_Vector(patch + 6, rm, temp, 3);
	 lip[2] = dot_product(Ivec + 6, temp);
	 lin[1] = dot_product(Ivec + 3, temp);
	
	 


	 Point_Point_Vector(patch, rm, temp, 3);
	 Rip[1] = sqrt(dot_product(temp, temp));
	 Rin[0]	= Rip[1];

	 Point_Point_Vector(patch + 3, rm, temp, 3);
	 Rip[0] = sqrt(dot_product(temp, temp));
	 Rin[2]	= Rip[0];

	 Point_Point_Vector(patch + 6, rm, temp, 3);
	 Rip[2] = sqrt(dot_product(temp, temp));
	 Rin[1]	= Rip[2];


	 //printf("Rip = %lf,%lf,%lf\n",Rip[0],Rip[1],Rip[2]);
	 //printf("Rin = %lf,%lf,%lf\n",Rin[0],Rin[1],Rin[2]);
	 //printf("lip = %lf,%lf,%lf\n",lip[0],lip[1],lip[2]);
	 //printf("lin = %lf,%lf,%lf\n",lin[0],lin[1],lin[2]);
	 //printf("Pi0 = %lf,%lf,%lf\n",Pi0[0],Pi0[1],Pi0[2]);
	 //printf("Ri02 = %lf,%lf,%lf\n",Ri02[0],Ri02[1],Ri02[2]);



	 for (m=0; m<3; m++){
		 if (fabs(Pi0[m]) > 1e-8){
			
			if (m == 2){
				//printf("%d\n",m);
				Point_Point_Vector(patch + 3 * m, rm, temp, 3);
			}
			else{
				Point_Point_Vector(patch + 3 * (1-m), rm, temp, 3);
			}
		
			MulScalar_Vector(Ivec + 3 * m, temp_1, lip[m], 3);
			Point_Point_Vector(temp, temp_1, temp_2, 3);
			MulScalar_Vector(temp_2, unit_Pi0 + 3 * m, 1.0/Pi0[m], 3);		
		 }

	 }

	I1P[0] = 0.0;
	I1P[1] = 0.0;
	I1P[2] = 0.0;

	for (m = 0; m < 3; ++m){

		temp[0] = Ri02[m]*log( (Rip[m]+lip[m])/(Rin[m]+lin[m]) )+lip[m]*Rip[m]-lin[m]*Rin[m]; 
		MulScalar_Vector(UI + 3 * m, temp_1, temp[0], 3);
	   	I1P[0] += temp_1[0];
		I1P[1] += temp_1[1];
		I1P[2] += temp_1[2];
	} 

	I1P[0] /= 2.0;
	I1P[1] /= 2.0;
	I1P[2] /= 2.0;
	temp_1[0] = I1P[0];
   	temp_1[1] = I1P[1];
   	temp_1[2] = I1P[2];


}



