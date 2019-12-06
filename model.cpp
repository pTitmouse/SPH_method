#include "Vmath.h"
#include "model.h"

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

vec3 ax(size_t n, std::vector<Particle>& particle, Kernel& kernel)
{
	vec3 ax (0, 0, 0);
	for (size_t i = 0; i < particle.size(); ++i)
	{
		if (i != n && neighbour(particle[n], particle[i]))
			ax += kernel.gradW(particle[n].pos, particle[i].pos, 1) *
			particle[i].get_mass() * U_lj(particle[n].pos, particle[i].pos) / dens(n, particle, kernel);
		else continue;
	}
	return ax * (-1) / particle[n].get_mass();
}

double U_lj(vec3 r1, vec3 r2)
{	
	double s = 1;
	double e = 1;
	double r = abs(r1 - r2);

	return 4 * e * (pow(s / r, 12) - pow(s / r, 6));
}

bool neighbour(Particle& p0, Particle& p1)
{
	if (abs(p1.pos.x - p0.pos.x) < 2 * p0.h
		&& abs(p1.pos.y - p0.pos.y) < 2 * p0.h
		&& abs(p1.pos.z - p0.pos.z) < 2 * p0.h)
		return 1;
	else
		return 0;
}