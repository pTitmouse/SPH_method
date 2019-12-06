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

	const int N = 100;
	double dt = 0.01;
	
	double t = 10;

	Kernel kernel(2);
	std::vector<Particle> particle(N);
	particle[0].set_pos(0, 0, 0);
	particle[1].set_pos(1, 0, 0);
	double x = 0;
	double y = 0;

	for (size_t i = 0; i < particle.size(); ++i)
	{
		if (i % 9 == 1 && i != 0)
		{
			x += 1;
			y = 0;
		}
		particle[i].set_pos(x, y, 0);
		y += 1;
	}

	std::ofstream file("data_2");

	

	size_t i = 0;
	while (i < t / dt)
	{

		eiler_scheme(particle, kernel, dt, i);
		system("cls");
		std::cout << "Progress:" << i * dt * 100 / t << " % "<<std::endl
				  << "interior time = " << i*dt << " s";
		++i;

		for (auto& p : particle)
			file << p.pos.x << " " << p.pos.y << " " << p.pos.z << " ";
		file << std::endl;

	}
	

}



