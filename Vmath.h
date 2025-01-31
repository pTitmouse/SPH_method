#pragma once
#include <iostream>
#include <vector>
#include<cmath>

#include "SPH.h"


class vec3
{
public:
	double x, y, z;
	
	vec3();
	vec3(double, double, double);
	vec3(double*);

	vec3 operator+ (vec3);
	vec3 operator- (vec3);
	vec3 operator* (double);
	vec3 operator/ (double);
	double operator* (vec3);
	vec3& operator= (const vec3&);
	vec3 operator= (double*); //!!!
	
	void set(double, double, double);
	void set(vec3);

	double abs();
	
	double proj(vec3); //Not checked!!!

};

double cos(vec3, vec3);
double abs(vec3);