#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void calcRotatedVector(double *init_vec, double *centre_vec, double *direc_vec,  double angle, double* final_vec)
	{

		double x = init_vec[0];
		double y = init_vec[1];
		double z = init_vec[2];

		double a = centre_vec[0];
		double b = centre_vec[1];
		double c = centre_vec[2];

		double u = direc_vec[0];
		double v = direc_vec[1];
		double w = direc_vec[2];

		double cos_angle = cos((2.*M_PI/360.) * angle);
		double sin_angle = sin((2.*M_PI/360.) * angle);

		final_vec[0] = (a*(v*v + w*w) - u*(b*v + c*w - u*x - v*y - w*z))*(1-cos_angle) + x*cos_angle + (- c*v + b*w - w*y + v*z)*sin_angle;
		final_vec[1] = (b*(u*u + w*w) - v*(a*u + c*w - u*x - v*y - w*z))*(1-cos_angle) + y*cos_angle + (c*u - a*w + w*x - u*z)*sin_angle;
		final_vec[2] = (c*(u*u + v*v) - w*(a*u + b*v - u*x - v*y - w*z))*(1-cos_angle) + z*cos_angle + (- b*u + a*v - v*x + u*y)*sin_angle;

	}

void makeVector(double a, double b, double c, double* vector)
	{
		vector[0] = a;
		vector[1] = b;
		vector[2] = c;
	}

int main () 
{

	double x = 10;
	double y = 0;
	double z = 0;

	double a = 0;
	double b = 0;
	double c = 0;

	double u = 0;
	double v = 0;
	double w = 1;

	double angle = 180;

	double* init_vec = (double*)malloc(sizeof(double)*3);
	double* centre_vec = (double*)malloc(sizeof(double)*3);
	double* direc_vec = (double*)malloc(sizeof(double)*3);
	double* final_vec = (double*)malloc(sizeof(double)*3);
	
	makeVector(x,y,z,init_vec);
	makeVector(a,b,c,centre_vec);
	makeVector(u,v,w,direc_vec);
	calcRotatedVector(init_vec, centre_vec, direc_vec, angle, final_vec);
	
	printf( "%f, %f, %f\n", final_vec[0], final_vec[1], final_vec[2]);
   return 0;
}