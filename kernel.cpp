//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//#include <stdio.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "Vmath.h"
#include "Model.h"


int main()
{
	/*
	Kernel kernel(2);

	std::vector<Particle> particle(1);

	std::ofstream data;
	data.open("data_500.txt");


	double dt = 5 * 1e-3;
	double t = 0;
	while (t < 400)
	{
		std::cout << t << std::endl;
		//time(particle, dt);

		for (auto& p : particle)
		{
			data << p.pos.x << " " << p.pos.y << " " << p.pos.z << " ";
		}
		data << t << std::endl;

		t += dt;
	}

	data.close();
	*/
	const int N = 1000;
	double dt = 0.1;
	
	double t = 600;

	Kernel kernel(2);
	std::vector<Particle> particle(N);

	size_t i = 0;
	while (i < t / dt)
	{

		eiler_scheme(particle, kernel, dt, i);

	}


}



