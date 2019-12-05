#include "Vmath.h"
#include "Model.h"

Particle::Particle()
{
	set_prop(1, 1);
	pos.set(0, 0, 0);
	vel.set(0, 0, 0);
	ax.set(0, 0, 0);
}

double Particle::get_mass()
{
	return mass;
}

double Particle::get_rad()
{
	return rad;
}

void Particle::set_prop(double m, double R)
{
	mass = m;
	rad = R;
}

void Particle::set_pos(double X, double Y, double Z)
{
	pos.set(X, Y, Z);
}
void Particle::set_pos(vec3 other)
{
	pos.set(other);
}

void Particle::set_vel(double X, double Y, double Z)
{
	vel.set(X, Y, Z);
}
void Particle::set_vel(vec3 other)
{
	vel.set(other);
}

void Particle::set_ax(double X, double Y, double Z)
{
	ax.set(X, Y, Z);
}
void Particle::set_ax(vec3 other)
{
	ax.set(other);
}


vec3 ax(size_t, std::vector<Particle>&, Kernel&)
{
	return 0;
}

double dens(size_t n, std::vector<Particle>& particle, Kernel& kernel)
{
	double dens = 0;
	for (Particle& p : particle)
		dens += p.get_mass()* kernel.W(particle[n].pos, p.pos, particle[n].h);

	return dens;
}

void adapt_h(size_t n, double dt, std::vector<Particle>& particle, Kernel& kernel)
{
	particle[n].h = particle[n].h * (1 + divVel(n, particle, kernel)/kernel.D);
}

double divVel(size_t n, std::vector<Particle>& particle, Kernel& kernel)
{
	double divV = 0;
	for (Particle& p : particle)
		divV += (p.vel - particle[n].vel) * 
				kernel.gradW(particle[n].pos, p.pos, particle[n].h) * 
				p.get_mass();

	return divV / dens(n, particle, kernel);

}

vec3 rotVel(size_t n, double dt, std::vector<Particle>& particle, Kernel& kernel)
{
	vec3 rotV = 0;
	for (Particle& p : particle)
		rotV += (p.vel - particle[n].vel) /
		kernel.gradW(particle[n].pos, p.pos, particle[n].h) *
		p.get_mass();

	return rotV / dens(n, particle, kernel);
}



void eiler_scheme(std::vector<Particle>& particle, Kernel& kernel, double dt, size_t iteration)
{
	size_t n = 0;
	for (Particle& p : particle)
	{	
		p.ax = ax(n, particle, kernel);
		p.vel += p.ax * dt;
		p.pos += p.vel * dt;

		++n;
	}
}

void advanced_scheme(std::vector<Particle>& particle, Kernel& kernel, double dt, size_t iteration)
{
	size_t n = 0;
	for (Particle& p : particle)
	{
		p.ax = ax(n, particle, kernel);
		p.vel += p.ax * dt;
		p.pos += p.vel * dt;

		++n;
	}
}