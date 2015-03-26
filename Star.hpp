#pragma once

#include "Density.hpp"
#include "Energy.hpp"
#include "Kappa.hpp"
#include "Luminosity.hpp"
#include "Mu.hpp"
#include "Mass.hpp"
#include "Tau.hpp"
#include "Pressure.hpp"
#include "Temperature.hpp"
#include "Opal.hpp"

#include <discpp.h>

using namespace std;


static Dislin dlin;

static char legendText[255] = { '/0' };
static int legendIndex = 1;

static void dislinInit(char* xName, char* yName, char* title) {
	dlin.disini();
	//dlin.nochek();
	dlin.pagera();
	dlin.name(xName, "X");
	dlin.name(yName, "Y");
	dlin.titlin(title, 1);

	legendIndex = 1;
	dlin.complx();
	dlin.texmod("ON");
};

class Star {
public:

	Density density;
	Energy energy;
	Kappa kappa;
	Luminosity luminosity;
	Mu mu;
	Mass mass;
	Tau tau;
	Pressure pressure;
	Temperature temperature;
	//Opal opal;

	long double R_star, T_star;
	long double T_c, rho_c;
	long double T_max, L_max, M_max, Rho_max,P_max;

	Star(long double T_c, long double rho_c, long double x, long double y, long double z/*, string opal_table*/) : /*opal(this, opal_table),*/ density(this, rho_c), energy(this), kappa(this), luminosity(this), mu(x,y,z), mass(this), tau(this), pressure(this), temperature(this, T_c), R_star(0), T_star(0), T_c(T_c), rho_c(rho_c) {
		//add plot points for values
		pushValues();
	}

	void solve();

	void graph();

	void push() {
		temperature.solver.push();
		density.solver.push();
		luminosity.solver.push();
		mass.solver.push();
	}

	void pop() {
		temperature.solver.pop();
		density.solver.pop();
		luminosity.solver.pop();
		mass.solver.pop();
	}

	void clear() {
		temperature.solver.clear();
		density.solver.clear();
		luminosity.solver.clear();
		mass.solver.clear();
	}

	void pushValues() {
		temperature.pushValues();
		density.pushValues();
		luminosity.pushValues();
		mass.pushValues();
		tau.pushValues();
		pressure.pushValues();
		kappa.pushValues();
	}

	long double frac_diff() {
		return (L_max - 4.0*pi*sigma_B*R_star*R_star*pow(T_star, 4)) / sqrt(4.0*pi*sigma_B*R_star*R_star*pow(T_star, 4));
	}

	void step();
	void iterate();
};