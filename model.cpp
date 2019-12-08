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
		if (neighbour(particle[n], p))
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
		if (neighbour(particle[n], p))
			divV += (p.vel - particle[n].vel) * 
					kernel.gradW(particle[n].pos, p.pos, particle[n].h) * 
					p.get_mass();

	return divV / dens(n, particle, kernel);

}

vec3 rotVel(size_t n, double dt, std::vector<Particle>& particle, Kernel& kernel)
{
	vec3 rotV = 0;
	for (Particle& p : particle)
		if (neighbour(particle[n], p))
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
	vec3 v_part;
	vec3 v_0;
	vec3 pos_0;

	for (Particle& p : particle)
	{
		p.ax = ax(n, particle, kernel);
		pos_0 = p.pos;
		v_0 = p.vel;

		v_part = v_0 + p.ax * dt * 0.5;
		p.pos = pos_0 + (v_0 + v_part) * dt * 0.25;

		p.ax = ax(n, particle, kernel);
		p.vel = v_0 + p.ax * dt;
		p.pos = pos_0 + (v_0 + v_part) * 0.5 *dt;

		++n;
	}
}

vec3 ax(size_t n, std::vector<Particle>& particle, Kernel& kernel)
{
	//100 remains magical number. Error occures while particle leaves the area of 100h x 100h
	Neighbour nei(100);
	nei.init(particle);

	vec3 ax (0, 0, 0);

	for (size_t i : nei(particle[n]))
	{
		if (i != n) // && neighbour(particle[n], particle[i]))
			ax += kernel.gradW(particle[n].pos, particle[i].pos, 1) *
			particle[i].get_mass() * U_lj(particle[n].pos, particle[i].pos) / dens(n, particle, kernel);
		else continue;
	}
	return ax * (-1) / particle[n].get_mass();
}

double U_lj(vec3 r1, vec3 r2)
{	
	double s = 0.891;
	double e = 300;
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


Neighbour::Neighbour(size_t size)
{
	web.resize(size);
	for (auto& num : web)
		num.resize(size);
}

void Neighbour::init(std::vector<Particle>& particle)
{
	for (auto& w1 : web)
		for (auto& w2 : w1)
			w2.resize(0);
	
	int x;
	int y;
	for (size_t i = 0; i < particle.size(); ++i)
	{
		 x = (int)(particle[i].pos.x / particle[i].h);
		 y = (int)(particle[i].pos.x / particle[i].h);
		web[x][y].push_back(i);
	}

}

std::vector<size_t> Neighbour::operator()(Particle& p)
{
	std::vector<size_t> res(0);
	size_t x = (size_t)(p.pos.x / p.h);
	size_t y = (size_t)(p.pos.y / p.h);
	for (int i = -1; i < 2; ++i)
		for (int j = -1; j < 2; ++j)
			res.insert(res.end(), web[x + i][y + j].begin(), web[x + i][y + j].end());
	
	return res;
}