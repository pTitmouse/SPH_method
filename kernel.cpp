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
	
	double t = 1;

	Kernel kernel(2);
	std::vector<Particle> particle(N);

	double x = 45;
	double y = 45;

	size_t c = 0;
	for (size_t i = 0; i < particle.size(); ++i)
	{

		if (i % 9 == 0 && i != 0)
		{	
			if (c % 2 == 0)
			{
				x = 45; 
			}
			else
				x = -0.5 + 45;
			y += 0.866;
			++c;
		}
		particle[i].set_pos(x, y, 0);
		x += 1;

		system("cls");
		std::cout << "Getenrating particles: " << i << "/" << N <<std::endl;
	}
	



	std::ofstream file("data_7");
	//std::ofstream energy("data_6_en");

	

	size_t i = 0;
	//double energy = 0;
	while (i < t / dt)
	{

		advanced_scheme(particle, kernel, dt, i);


		system("cls");
		std::cout << "Progress:" << i * dt * 100 / t << " % "<<std::endl
				  << "interior time = " << i*dt << " s";
		

		for (auto& p : particle)
			file << p.pos.x << " " << p.pos.y << " " << p.pos.z << " ";
		file << std::endl;

		++i;
	}

	file.close();

}



