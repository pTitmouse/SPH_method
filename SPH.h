#pragma once

#include "Vmath.h"


class Kernel
{
public:
	size_t D;
	size_t C = 0;
	
	Kernel(size_t);

	double W(vec3, vec3, double);
	vec3 gradW(vec3, vec3, double);

};
