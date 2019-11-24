#pragma once
#pragma once
#include "Vmath.h"

const double EARTH_MASS = 5.97e24;
const double EARTH_RAD = 6.4e6;
const double G = 6.7e-11;


class Particle
{
private:
	double mass, rad; //Mass and radius

public:

	Particle();

	vec3 pos; //Position radius-vector
	vec3 vel; //Velosity vector
	vec3 ax; //Axeleration vector

	double h; //Smooth radius
	
	double get_mass();
	double get_rad();

	void set_prop(double, double);

	void set_pos(double, double, double);
	void set_pos(vec3);

	void set_vel(double, double, double);
	void set_vel(vec3);

	void set_ax(double, double, double);
	void set_ax(vec3);
};

vec3 ax(std::vector<Particle>&); // dv/dt
double energy(std::vector<Particle>&); 
void time(std::vector<Particle>&, double);

double dens(vec3, double, std::vector<Particle>&, Kernel&); //density
double pres(vec3, double, std::vector<Particle>&, Kernel&); //pressure
double gamma(vec3, double, std::vector<Particle>&, Kernel&); //adiabatic exponent

double adapt_h(size_t, std::vector<Particle>&, Kernel&); //adaptive smooth radius