#include "Vmath.h"
#include "Model.h"

Particle::Particle()
{
	set_prop(EARTH_MASS, EARTH_RAD);
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


void time(std::vector<Particle>& particle, double dt)
{
	vec3 F(0, 0, 0);
	for (size_t i = 0; i < particle.size(); ++i)
	{
		//Computing total force
		/*
		for (size_t j = 0; j < particle.size(); ++j)
			if (i != j) F = F + Force(particle[i], particle[j]);
		*/
		particle[i].ax = ax(particle);
		particle[i].vel = particle[i].vel + particle[i].ax * dt;
		particle[i].pos = particle[i].pos + particle[i].vel * dt;

		F.set(0, 0, 0);
	}
}

vec3 ax(std::vector<Particle>& particle)
{
	return 0;
}

double dens(vec3 r, double h, std::vector<Particle>& particle, Kernel& kernel)
{
	double dens = 0;
	for (Particle& p : particle)
		dens += p.get_mass()* kernel.W(r, p.pos, h);

	return dens;
}

double h_value(size_t n, std::vector<Particle>& particle, Kernel& kernel)
{
	double h0 = 1;
	double dens0 = 1;

	return 1;

}